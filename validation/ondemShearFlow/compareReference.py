#!/usr/bin/env python3
"""Compare time-matched particle snapshots without claiming steady-flow validation.

Profiles here use particle-center bins, not MercuryCG coarse graining, and are
instantaneous transient diagnostics. Particle rows keep the original ordering:
480 prescribed wall grains, then 3000 mobile grains. The number of force records
does not by itself indicate whether a snapshot contains every current contact.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import sys

import numpy as np


PARTICLE_COUNT = 3480
WALL_COUNT = 480
TIME_TOLERANCE = 1.0e-9
INITIAL_TOLERANCE = 1.0e-10
BOUNDS_TOLERANCE = 1.0e-8  # Mercury example headers round the Y extent.
MASS = math.pi / 6.0
INERTIA = math.pi / 60.0
PROFILE_EDGES = np.linspace(1.0, 21.0, 101)


@dataclass(frozen=True)
class Frame:
    """One complete Mercury-format particle snapshot, with stable row order."""

    path: Path
    time: float
    minimum: np.ndarray
    maximum: np.ndarray
    values: np.ndarray


def read_frame(path: Path) -> Frame:
    """Reject missing rows, malformed columns, and non-finite particle values."""
    with path.open(encoding="utf-8") as source:
        header = source.readline().split()
        if len(header) != 8:
            raise ValueError(f"{path}: expected an eight-value particle header")
        count_float = float(header[0])
        if not math.isfinite(count_float) or count_float != PARTICLE_COUNT:
            raise ValueError(f"{path}: expected {PARTICLE_COUNT} ordered particles")
        fields = np.asarray([float(value) for value in header[1:]], dtype=float)
        if not np.isfinite(fields).all() or np.any(fields[4:7] <= fields[1:4]):
            raise ValueError(f"{path}: non-finite time or invalid bounds")
        records = [line.split() for line in source if line.strip()]
    if len(records) != PARTICLE_COUNT or any(len(row) != 14 for row in records):
        raise ValueError(f"{path}: expected {PARTICLE_COUNT} complete fourteen-column rows; found {len(records)}")
    values = np.asarray(records, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"{path}: non-finite particle data")
    if not np.allclose(values[:, 6], 0.5, rtol=0.0, atol=INITIAL_TOLERANCE):
        raise ValueError(f"{path}: the normalized benchmark requires unit-diameter spheres")
    return Frame(path, float(fields[0]), fields[1:4], fields[4:7], values)


def minimum_image(displacement: np.ndarray, extent: np.ndarray) -> np.ndarray:
    """Apply periodic minimum images in X/Y only; Z is never wrapped."""
    result = displacement.copy()
    result[:, :2] -= extent[:2] * np.rint(result[:, :2] / extent[:2])
    return result


def initial_checks(reference: Frame, simulation: Frame) -> dict:
    """Establish initial row correspondence rather than matching nearest grains."""
    if abs(reference.time) > TIME_TOLERANCE or abs(simulation.time) > TIME_TOLERANCE:
        raise ValueError("Frame 0000 must represent time zero in both datasets")
    delta = minimum_image(simulation.values[:, :3] - reference.values[:, :3], reference.maximum - reference.minimum)
    result = {
        "position_max_absolute_error": float(np.abs(delta).max()),
        "velocity_max_absolute_error": float(np.abs(simulation.values[:, 3:6] - reference.values[:, 3:6]).max()),
        "radius_max_absolute_error": float(np.abs(simulation.values[:, 6] - reference.values[:, 6]).max()),
        "angular_velocity_max_absolute_error": float(np.abs(simulation.values[:, 10:13] - reference.values[:, 10:13]).max()),
        "tolerance": INITIAL_TOLERANCE,
        "stable_row_order": "Original rows retained; no nearest-particle remapping",
    }
    result["passed"] = all(result[key] <= INITIAL_TOLERANCE for key in (
        "position_max_absolute_error", "velocity_max_absolute_error", "radius_max_absolute_error", "angular_velocity_max_absolute_error"))
    for frame in (reference, simulation):
        if not np.allclose(frame.values[:240, 3:6], [-1.0, 0.0, 0.0], rtol=0.0, atol=INITIAL_TOLERANCE):
            raise ValueError(f"{frame.path}: lower wall velocity must be (-1, 0, 0)")
        if not np.allclose(frame.values[240:480, 3:6], [1.0, 0.0, 0.0], rtol=0.0, atol=INITIAL_TOLERANCE):
            raise ValueError(f"{frame.path}: upper wall velocity must be (1, 0, 0)")
    if not result["passed"]:
        raise ValueError(f"Initial row/kinematic agreement failed: {result}")
    return result


def state_metrics(frame: Frame, frame_index: str, dataset: str) -> dict:
    """Return mobile-grain diagnostics; prescribed walls never enter energies."""
    mobile = frame.values[WALL_COUNT:]
    velocity = mobile[:, 3:6]
    angular_velocity = mobile[:, 10:13]
    momentum = MASS * velocity.sum(axis=0)
    return {
        "frame": frame_index,
        "dataset": dataset,
        "time": frame.time,
        "mobile_count": len(mobile),
        "finite_mobile_count": int(np.isfinite(mobile).all(axis=1).sum()),
        "mean_velocity_x": float(velocity[:, 0].mean()),
        "translational_kinetic_energy": float(0.5 * MASS * np.square(velocity).sum()),
        "rotational_kinetic_energy": float(0.5 * INERTIA * np.square(angular_velocity).sum()),
        "minimum_center_z": float(mobile[:, 2].min()),
        "maximum_center_z": float(mobile[:, 2].max()),
        "momentum_x": float(momentum[0]),
        "momentum_y": float(momentum[1]),
        "momentum_z": float(momentum[2]),
    }


def difference_metrics(reference: Frame, simulation: Frame) -> dict:
    """Compare identical particle indices, not independently sorted positions."""
    displacement = minimum_image(simulation.values[WALL_COUNT:, :3] - reference.values[WALL_COUNT:, :3], reference.maximum - reference.minimum)
    velocity_difference = simulation.values[WALL_COUNT:, 3:6] - reference.values[WALL_COUNT:, 3:6]
    result = {
        "position_vector_rmse": float(np.sqrt(np.square(displacement).sum(axis=1).mean())),
        "position_max_vector_error": float(np.linalg.norm(displacement, axis=1).max()),
        "velocity_vector_rmse": float(np.sqrt(np.square(velocity_difference).sum(axis=1).mean())),
        "angular_velocity_vector_rmse": float(np.sqrt(np.square(simulation.values[WALL_COUNT:, 10:13] - reference.values[WALL_COUNT:, 10:13]).sum(axis=1).mean())),
    }
    for axis, label in enumerate("xyz"):
        result[f"position_rmse_{label}"] = float(np.sqrt(np.square(displacement[:, axis]).mean()))
        result[f"velocity_rmse_{label}"] = float(np.sqrt(np.square(velocity_difference[:, axis]).mean()))
    return result


def center_profiles(frame: Frame) -> dict:
    """Bin centers in 100 horizontal slices, explicitly without CG smoothing."""
    mobile = frame.values[WALL_COUNT:]
    indices = np.searchsorted(PROFILE_EDGES, mobile[:, 2], side="right") - 1
    indices[mobile[:, 2] == PROFILE_EDGES[-1]] = len(PROFILE_EDGES) - 2
    included = (indices >= 0) & (indices < 100)
    counts = np.bincount(indices[included], minlength=100)
    velocity_sum = np.zeros((100, 3))
    np.add.at(velocity_sum, indices[included], mobile[included, 3:6])
    means = np.divide(velocity_sum, counts[:, None], out=np.full((100, 3), np.nan), where=counts[:, None] > 0)
    squared_deviation = np.zeros(100)
    np.add.at(squared_deviation, indices[included], np.square(mobile[included, 3:6] - means[indices[included]]).sum(axis=1))
    temperature = np.divide(squared_deviation, 3.0 * counts, out=np.full(100, np.nan), where=counts > 0)
    bin_volume = float(np.prod((frame.maximum - frame.minimum)[:2]) * (PROFILE_EDGES[1] - PROFILE_EDGES[0]))
    return {
        "center_z": 0.5 * (PROFILE_EDGES[1:] + PROFILE_EDGES[:-1]),
        "count": counts,
        "velocity_x": means[:, 0],
        "center_binned_volume_fraction": counts * MASS / bin_volume,
        "granular_temperature": temperature,
        "centers_outside_profile_gap": int((~included).sum()),
    }


def profile_difference(reference: dict, simulation: dict) -> dict:
    """Measure field differences only where both bins contain particles."""
    occupied = (reference["count"] > 0) & (simulation["count"] > 0)
    result = {"jointly_occupied_bins": int(occupied.sum())}
    for quantity in ("velocity_x", "center_binned_volume_fraction", "granular_temperature"):
        selected = occupied if quantity != "center_binned_volume_fraction" else np.ones(100, dtype=bool)
        delta = simulation[quantity][selected] - reference[quantity][selected]
        result[f"{quantity}_rmse"] = float(np.sqrt(np.square(delta).mean())) if len(delta) else None
    return result


def write_csv(path: Path, rows: list[dict]) -> None:
    """Write a deterministic table with all observed field names."""
    names = list(dict.fromkeys(name for row in rows for name in row))
    with path.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)


def force_file_status(directory: Path, prefix: str) -> dict:
    """Count contact records without inferring incompleteness from sparse data."""
    counts = {}
    for path in sorted(directory.glob(f"{prefix}.fstat.*")):
        with path.open(encoding="utf-8") as source:
            counts[path.name] = sum(bool(line.strip()) and not line.lstrip().startswith("#") for line in source)
    return {"record_counts": counts, "mercury_cg_contact_stress_comparison_performed": False,
            "note": "Sparse contact records may be complete transient snapshots; record count alone does not establish missing data."}


def draw_comparison(path: Path, pairs: list[dict]) -> None:
    """Plot matched transient profiles, not a steady-state validation score."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    noninitial = [pair for pair in pairs if pair["frame"] != "0000"]
    if not noninitial:
        noninitial = pairs[:1]
    figure, axes = plt.subplots(len(noninitial), 3, figsize=(12.0, 4.0 * len(noninitial)), squeeze=False, constrained_layout=True)
    columns = (
        ("velocity_x", "Mean streamwise velocity"),
        ("center_binned_volume_fraction", "Center-binned volume fraction"),
        ("granular_temperature", "Granular temperature"),
    )
    for row, pair in enumerate(noninitial):
        for column, (quantity, label) in enumerate(columns):
            axis = axes[row, column]
            for dataset, color, style in (("reference", "black", "-"), ("simulation", "#0078D4", "--")):
                profile = pair[f"{dataset}_profile"]
                axis.plot(profile[quantity], profile["center_z"], color=color, linestyle=style, label="Reference" if dataset == "reference" else "FunDEM")
            axis.set(xlabel=label, ylabel="Center Z / diameter", ylim=(1.0, 21.0), title=f"t = {pair['time']:.9g}; frame {pair['frame']}")
            axis.grid(alpha=0.2)
            axis.legend(loc="best")
    figure.suptitle("Time-matched transient snapshots\n100 particle-center bins; not MercuryCG or steady-state validation")
    figure.savefig(path, dpi=180)
    plt.close(figure)


def compare(reference_directory: Path, simulation_directory: Path, output_directory: Path) -> dict:
    """Validate matching first, then write comparison artifacts to a new directory."""
    if output_directory.exists() and (not output_directory.is_dir() or any(output_directory.iterdir())):
        raise ValueError(f"Output must be absent or empty; refusing to overwrite {output_directory}")
    reference_paths = sorted(reference_directory.glob("Mercury_out.data.*"))
    if not reference_paths:
        raise ValueError(f"No Mercury_out.data snapshots in {reference_directory}")
    pairs = []
    partial = []
    initial_result = None
    metric_rows = []
    profile_rows = []
    summaries = []
    for reference_path in reference_paths:
        frame_index = reference_path.name.rsplit(".", 1)[-1]
        if not frame_index.isdigit():
            continue
        try:
            reference = read_frame(reference_path)
        except ValueError as error:
            if int(frame_index) == 0:
                raise
            partial.append({"frame": frame_index, "reason": str(error)})
            continue
        simulation_path = simulation_directory / f"ShearFlow.data.{frame_index}"
        if not simulation_path.is_file():
            raise ValueError(f"Missing index-matched FunDEM frame: {simulation_path}")
        simulation = read_frame(simulation_path)
        if abs(reference.time - simulation.time) > TIME_TOLERANCE:
            raise ValueError(f"Frame {frame_index} times differ ({reference.time:.16g} versus {simulation.time:.16g}); no nearest-time alignment is permitted")
        if not np.allclose(reference.minimum, simulation.minimum, rtol=0.0, atol=BOUNDS_TOLERANCE) or not np.allclose(reference.maximum, simulation.maximum, rtol=0.0, atol=BOUNDS_TOLERANCE):
            raise ValueError(f"Frame {frame_index}: domain bounds disagree")
        if int(frame_index) == 0:
            initial_result = initial_checks(reference, simulation)
        reference_profile = center_profiles(reference)
        simulation_profile = center_profiles(simulation)
        pair = {"frame": frame_index, "time": reference.time, "reference_profile": reference_profile, "simulation_profile": simulation_profile}
        pairs.append(pair)
        differences = difference_metrics(reference, simulation)
        for dataset, frame, profile in (("reference", reference, reference_profile), ("simulation", simulation, simulation_profile)):
            row = state_metrics(frame, frame_index, dataset)
            row.update(differences if dataset == "simulation" else {})
            row["centers_outside_profile_gap"] = profile["centers_outside_profile_gap"]
            metric_rows.append(row)
            for index in range(100):
                profile_rows.append({"frame": frame_index, "time": reference.time, "dataset": dataset, "bin": index, **{
                    key: float(profile[key][index]) for key in ("center_z", "count", "velocity_x", "center_binned_volume_fraction", "granular_temperature")}})
        summaries.append({"frame": frame_index, "reference_time": reference.time, "simulation_time": simulation.time,
                          "particle_differences": differences, "center_bin_profile_differences": profile_difference(reference_profile, simulation_profile)})
    if initial_result is None:
        raise ValueError("Complete frame 0000 is required to establish stable initial row correspondence")
    summary = {
        "comparison_scope": "Instantaneous time-matched transient particle snapshots; not completed ON-DEM steady-flow validation",
        "validation_pass_claimed": False,
        "partial_reference": bool(partial) or len(pairs) < 3,
        "skipped_reference_frames": partial,
        "matched_frame_count": len(pairs),
        "initial_checks": initial_result,
        "population": {"total": PARTICLE_COUNT, "excluded_prescribed_wall_rows": WALL_COUNT, "mobile": PARTICLE_COUNT - WALL_COUNT},
        "normalization": {"diameter": 1.0, "density": 1.0, "mobile_mass": MASS, "mobile_sphere_inertia": INERTIA},
        "matching": {"method": "Exact file-index and original row order; no nearest-particle or nearest-time remapping", "time_tolerance": TIME_TOLERANCE, "bounds_tolerance": BOUNDS_TOLERANCE},
        "profiles": {"method": "100 instantaneous particle-center bins in Z=[1,21]; horizontal averaging over periodic X/Y", "granular_temperature": "sum_i |v_i - mean_bin(v)|^2 / (3*N_bin)", "concentration": "N_bin * pi/6 / (Lx*Ly*bin_height); center-binned, not geometric slice intersection", "mercury_cg_equivalent": False},
        "reference_contacts": force_file_status(reference_directory, "Mercury_out"),
        "simulation_contacts": force_file_status(simulation_directory, "ShearFlow"),
        "frames": summaries,
        "limitations": [
            "MercuryCG contact pressure and shear-stress fields were not computed; sparse force-file counts alone are not evidence of incomplete data.",
            "Only transient particle snapshots are compared; no late-time averaging, time-step convergence, or MercuryCG validation was performed.",
            "Discrete chaotic trajectory divergence is diagnostic, not by itself a pass/fail criterion for bulk flow agreement.",
            "Empty profile bins have undefined mean velocity and temperature and are left as NaN in profiles.csv.",
        ],
    }
    output_directory.mkdir(parents=True, exist_ok=True)
    write_csv(output_directory / "metrics.csv", metric_rows)
    write_csv(output_directory / "profiles.csv", profile_rows)
    draw_comparison(output_directory / "comparison.png", pairs)
    with (output_directory / "summary.json").open("w", encoding="utf-8") as destination:
        json.dump(summary, destination, indent=2, allow_nan=False)
        destination.write("\n")
    return summary


def main() -> int:
    """Compare local files only; never download data or advance a simulation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True, help="Directory containing original Mercury_out.data/fstat snapshots")
    parser.add_argument("--simulation", type=Path, required=True, help="Directory containing time-matched ShearFlow.data/fstat snapshots")
    parser.add_argument("--output", type=Path, required=True, help="New or empty destination directory")
    arguments = parser.parse_args()
    try:
        result = compare(arguments.reference, arguments.simulation, arguments.output)
    except (OSError, ValueError) as error:
        print(f"Comparison refused: {error}", file=sys.stderr)
        return 1
    print(f"Compared {result['matched_frame_count']} exactly time-matched frames; initial correspondence passed.")
    print(f"Transient diagnostics written to {arguments.output}; no steady-flow validation pass is claimed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
