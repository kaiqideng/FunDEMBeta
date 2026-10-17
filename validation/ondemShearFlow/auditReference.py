#!/usr/bin/env python3
"""Audit the pinned Mercury reference snapshots without running a simulation.

This standard-library-only audit verifies complete periodic contact coverage and
tests a velocity-Verlet output-staggering interpretation. The inferred dashpot
coefficient is numerical evidence, not a substitute for the author's missing
simulation script. No model coefficients are fitted or changed by this script.
"""

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path


REFERENCE_SHA256 = {
    "data.0000": "f521cf8825134af59b32d86685cff60096679322ea3df4143bd570335affbb71",
    "data.0001": "ccb376ed5cce06f22734e2bd9e400f598d0bef70a3b764419ea40dbe3a4ee318",
    "data.0002": "da3c37c1f989b649a93b92b72884953466e69e74e97137f4a83efddaae1e240f",
    "fstat.0000": "91ed9184fb190f0dc858b7e77694e0bdf8e433cc004c1a9ff406f2f19cbdb5ad",
    "fstat.0001": "4d85e7244a589cfba64ec7ed12f67faa2015e13e70bc04ae42c4b97cf2fffe1e",
    "fstat.0002": "832d2e78e3d7b445228c6e9eeae125887145cce20165d23b522de1d88eb3f607",
}
PERIODS = (20.0, 10.3923048454133, 22.0)
WALL_COUNT = 480
DENSITY = 1.0
NORMAL_STIFFNESS = 200000.0
PUBLISHED_DAMPING = 50.0
OUTPUT_INTERVAL_STEPS = 50000


def dot(first, second):
    return sum(a * b for a, b in zip(first, second))


def norm(value):
    return math.sqrt(dot(value, value))


def read_verified(directory, suffix):
    path = directory / ("Mercury_out." + suffix)
    content = path.read_bytes()
    digest = hashlib.sha256(content).hexdigest()
    if digest != REFERENCE_SHA256[suffix]:
        raise ValueError("The pinned reference checksum does not match: " + str(path))
    return content.decode().splitlines(), digest


def geometric_contacts(particles):
    """Find every positive-overlap pair, excluding prescribed-wall/wall pairs."""
    diameter = 2.0 * particles[0][6]
    counts = tuple(max(1, math.floor(length / diameter)) for length in PERIODS)
    widths = tuple(length / count for length, count in zip(PERIODS, counts))
    bins = defaultdict(list)

    def cell(position):
        return tuple(min(counts[k] - 1, max(0, int(position[k] / widths[k]))) for k in range(3))

    for index, particle in enumerate(particles):
        bins[cell(particle)].append(index)
    found = {}
    for first, particle in enumerate(particles):
        base = cell(particle)
        for dz in (-1, 0, 1):
            z = base[2] + dz
            if z < 0 or z >= counts[2]:
                continue
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    neighbor = ((base[0] + dx) % counts[0], (base[1] + dy) % counts[1], z)
                    for second in bins[neighbor]:
                        if second <= first or second < WALL_COUNT:
                            continue
                        displacement = [particle[k] - particles[second][k] for k in range(3)]
                        for axis in (0, 1):
                            displacement[axis] -= round(displacement[axis] / PERIODS[axis]) * PERIODS[axis]
                        overlap = particle[6] + particles[second][6] - norm(displacement)
                        if overlap > 1.0e-12:
                            found[(first, second)] = overlap
    return found


def audit_frame(directory, frame):
    data_lines, data_digest = read_verified(directory, "data." + frame)
    force_lines, force_digest = read_verified(directory, "fstat." + frame)
    header = list(map(float, data_lines[0].split()))
    particles = [list(map(float, row.split())) for row in data_lines[1:] if row.strip()]
    contacts = [list(map(float, row.split())) for row in force_lines if row.strip() and not row.startswith("#")]
    if len(particles) != int(header[0]) or len(particles) != 3480 or any(len(row) != 14 for row in particles):
        raise ValueError("Invalid complete-particle frame " + frame)
    if any(len(row) != 16 for row in contacts):
        raise ValueError("Invalid contact record in frame " + frame)
    found = geometric_contacts(particles)
    expected_records = set()
    for first, second in found:
        expected_records.add((second, first))
        if first >= WALL_COUNT:
            expected_records.add((first, second))
    actual_records = [(int(row[1]), int(row[2])) for row in contacts]
    complete = set(actual_records) == expected_records and len(actual_records) == len(expected_records)
    if not complete:
        raise ValueError("Contact records do not cover all geometric overlaps in frame " + frame)

    # Reciprocal mobile/mobile records each supply the force on their own first
    # particle. A wall record supplies only the force on the mobile particle.
    forces = defaultdict(lambda: [0.0, 0.0, 0.0])
    for row in contacts:
        first = int(row[1])
        for axis in range(3):
            forces[first][axis] += row[8] * row[10 + axis] + row[9] * row[13 + axis]

    time = header[1]
    dt = time / (int(frame) * OUTPUT_INTERVAL_STEPS) if int(frame) else None
    diagnostics = []
    for row in contacts:
        first, second = int(row[1]), int(row[2])
        if first > second and second >= WALL_COUNT:
            continue  # Report each physical contact once, not reciprocal rows.
        normal = row[10:13]
        full_relative_velocity = [particles[first][k + 3] - particles[second][k + 3] for k in range(3)]
        normal_velocity = dot(full_relative_velocity, normal)
        correction = [0.0, 0.0, 0.0]
        for index, sign in ((first, 1), (second, -1)):
            if index >= WALL_COUNT:
                mass = DENSITY * 4.0 * math.pi * particles[index][6] ** 3 / 3.0
                for axis in range(3):
                    correction[axis] += sign * 0.5 * dt * forces[index][axis] / mass
        half_normal_velocity = normal_velocity - dot(correction, normal)
        coefficient = (NORMAL_STIFFNESS * row[6] - row[8]) / half_normal_velocity
        published_prediction = max(0.0, NORMAL_STIFFNESS * row[6] - PUBLISHED_DAMPING * half_normal_velocity)
        diagnostics.append({
            "particle_indices": [first, second],
            "overlap": row[6],
            "normal_force": row[8],
            "full_step_normal_velocity": normal_velocity,
            "reconstructed_half_step_normal_velocity": half_normal_velocity,
            "inferred_normal_dashpot": coefficient,
            "normal_force_using_published_dashpot": published_prediction,
            "published_prediction_relative_error": (published_prediction - row[8]) / row[8],
        })
    mass = DENSITY * 4.0 * math.pi * particles[0][6] ** 3 / 3.0
    mobile = particles[WALL_COUNT:]
    return {
        "frame": frame,
        "time": time,
        "particle_count": len(particles),
        "physical_contact_count": len(found),
        "contact_record_count": len(contacts),
        "complete_contact_coverage": complete,
        "data_sha256": data_digest,
        "fstat_sha256": force_digest,
        "inferred_sample_time_step": dt,
        "mobile_translational_kinetic_energy": 0.5 * mass * sum(dot(p[3:6], p[3:6]) for p in mobile),
        "mobile_maximum_velocity": max(norm(p[3:6]) for p in mobile),
        "contact_diagnostics": diagnostics,
    }


def collision_time_step(damping):
    effective_mass = math.pi / 12.0
    return math.pi / math.sqrt(NORMAL_STIFFNESS / effective_mass - (damping / (2.0 * effective_mass)) ** 2) / 50.0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, required=True, help="Folder with the six pinned Mercury_out.data/fstat files.")
    parser.add_argument("--output", type=Path, help="Write a new JSON report; existing files are not overwritten.")
    arguments = parser.parse_args()
    report = {
        "reference_commit": "56b7cfdfd774cf5ae057a1c5c90e5191a0deb97b",
        "assumptions": [
            "Unit material density and sphere diameter give mass pi/6.",
            "Particle rows 0 through 479 are prescribed wall grains; their force does not accelerate them.",
            "The reference output interval is 50000 steps; its timestamps determine the sample step size.",
            "The force snapshot was evaluated at half-step velocities and the exported velocities include the subsequent half kick.",
            "No missing external forces act on the mobile grains, and the reference snapshots include every positive-overlap contact.",
            "The inferred dashpot is conditional numerical evidence; the author's simulation script is not supplied.",
            "Tangential dashpot, history transport, and return mapping cannot be uniquely established by these scalar elongation samples.",
        ],
        "published_normal_dashpot": PUBLISHED_DAMPING,
        "collision_time_over_50_for_published_dashpot": collision_time_step(PUBLISHED_DAMPING),
        "collision_time_over_50_for_dashpot_15": collision_time_step(15.0),
        "frames": [audit_frame(arguments.reference_dir, frame) for frame in ("0000", "0001", "0002")],
    }
    text = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if arguments.output:
        with arguments.output.open("x", encoding="utf-8") as stream:
            stream.write(text)
    print(text, end="")


if __name__ == "__main__":
    main()
