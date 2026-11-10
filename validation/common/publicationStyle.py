"""Publication figures from immutable historical checkpoint tables.

No solver is invoked. Sparse checkpoint symbols are never joined or fitted.
The generators compare every rendered byte in --check mode; the recorded
Matplotlib version and STIX font files identify the rendering environment.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np

RUNS = ("cpu", "gpu", "cpu-dt-half", "gpu-dt-half")
LABELS = (r"CPU, $\Delta t$", r"GPU, $\Delta t$", r"CPU, $\Delta t/2$", r"GPU, $\Delta t/2$")
COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7")
MARKERS = ("o", "x", "s", "+")
SIZES = (6.2, 5.5, 4.4, 4.2)
WIDTH_MM = 180.0
STYLE = {
    "font.family": "STIXGeneral", "font.size": 8.5,
    "mathtext.fontset": "stix", "axes.labelsize": 9,
    "axes.titlesize": 9, "axes.linewidth": .65,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#333333", "axes.labelcolor": "black",
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.major.width": .65, "ytick.major.width": .65,
    "grid.color": "#dddddd", "grid.linewidth": .4,
    "legend.fontsize": 8.5, "legend.frameon": False,
    "lines.linewidth": .8, "lines.markersize": 4.5,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "none", "svg.hashsalt": "fundem-publication-20261031",
    "savefig.facecolor": "white", "figure.facecolor": "white",
    "axes.formatter.useoffset": False, "axes.unicode_minus": True,
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def num(row, field):
    value = float(row[field])
    if not math.isfinite(value):
        raise ValueError(f"Nonfinite value: {field}")
    return value


def settings():
    plt.rcParams.update(STYLE)


def marker(axis, x, y, index, **kwargs):
    return axis.plot(x, y, linestyle="none", marker=MARKERS[index],
                     markersize=SIZES[index], markerfacecolor="none",
                     markeredgewidth=.85, color=COLORS[index], zorder=3 + index,
                     clip_on=False, **kwargs)[0]


def decorate(axis, letter):
    axis.text(0, 1.075, f"({letter})", transform=axis.transAxes,
              va="bottom", ha="left", fontsize=10)
    axis.grid(axis="y", zorder=0)
    axis.xaxis.set_major_locator(MaxNLocator(nbins=5))
    axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
    axis.ticklabel_format(style="plain", useOffset=False)


def canvas(height_mm=122.0, rows=2, cols=2, reference=True):
    fig, axes = plt.subplots(rows, cols, figsize=(WIDTH_MM / 25.4, height_mm / 25.4), squeeze=False)
    fig.subplots_adjust(left=.103, right=.983, bottom=.109, top=.875,
                        wspace=.34, hspace=.63)
    handles = [Line2D([], [], color=COLORS[i], marker=MARKERS[i], linestyle="none",
                      markerfacecolor="none", markersize=SIZES[i], markeredgewidth=.85)
               for i in range(4)]
    labels = list(LABELS)
    if reference:
        if reference == "analytical":
            handles.insert(0, Line2D([], [], color="black", linestyle="--", linewidth=.8))
            labels.insert(0, "Analytical reference")
        else:
            handles.insert(0, Line2D([], [], color="black", marker="D", linestyle="none",
                                      markerfacecolor="none", markersize=4, markeredgewidth=.75))
            labels.insert(0, "Reference")
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.54, .986),
               ncol=len(labels), columnspacing=1.1, handletextpad=.35,
               handlelength=1.0, borderaxespad=0)
    for i, ax in enumerate(axes.flat):
        decorate(ax, chr(ord("a") + i))
    return fig, axes


def bounds(fig, case):
    """Detect visible text outside the final canvas and overlapping labels."""
    fig.canvas.draw()
    renderer, frame = fig.canvas.get_renderer(), fig.bbox
    outside = set()
    for ax in fig.axes:
        for coordinate, limits in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
            low, high = sorted(limits)
            for tick in coordinate.get_major_ticks() + coordinate.get_minor_ticks():
                if not low <= tick.get_loc() <= high:
                    outside.update((tick.label1, tick.label2))
    texts = []
    for text in fig.findobj(matplotlib.text.Text):
        if text in outside or not text.get_visible() or not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        if box.x0 < frame.x0 - 1 or box.y0 < frame.y0 - 1 or box.x1 > frame.x1 + 1 or box.y1 > frame.y1 + 1:
            raise ValueError(f"Clipped text in {case}: {text.get_text()}")
        texts.append((text, box))
    for ax in fig.axes:
        tick_texts = [t for t in ax.get_xticklabels() + ax.get_yticklabels()
                      if t.get_visible() and t.get_text() and t not in outside]
        for i, text in enumerate(tick_texts):
            for other in tick_texts[i + 1:]:
                if text.get_window_extent(renderer).overlaps(other.get_window_extent(renderer)):
                    raise ValueError(f"Overlapping tick labels in {case}")
    occupied = [ax.get_tightbbox(renderer) for ax in fig.axes]
    for i, box in enumerate(occupied):
        if any(box.overlaps(other) for other in occupied[i + 1:]):
            raise ValueError(f"Overlapping panel labels or axes in {case}")
        if any(box.overlaps(entry.get_window_extent(renderer)) for entry in fig.legends):
            raise ValueError(f"Legend/panel overlap in {case}")


def render(fig, case):
    fig.tight_layout(pad=.85, w_pad=2.4, h_pad=2.6, rect=(0, 0, 1, .91))
    bounds(fig, case)
    result = {}
    for fmt in ("svg", "pdf", "png"):
        buffer = io.BytesIO()
        if fmt == "svg":
            metadata = {"Date": None, "Creator": "FunDEM publication checkpoint generator",
                        "Description": "Immutable historical checkpoints; no new solver execution"}
        elif fmt == "pdf":
            metadata = {"CreationDate": None, "ModDate": None,
                        "Creator": "FunDEM publication checkpoint generator",
                        "Producer": "Matplotlib", "Title": case,
                        "Subject": "Historical checkpoint comparison"}
        else:
            metadata = {"Software": "FunDEM publication checkpoint generator",
                        "Description": "Immutable historical checkpoints; 180 mm at 600 dpi"}
        fig.savefig(buffer, format=fmt, dpi=600, metadata=metadata)
        result[fmt] = buffer.getvalue()
    return result


def selected_body(data):
    names = list(dict.fromkeys(row.get("instance", "sphere") for row in data))
    return next((name for name in ("left", "sphere", "small", "lower") if name in names), names[0])


def component(data, quantity):
    candidates = [row for row in data if row.get("quantity") == quantity]
    return max("xyz", key=lambda axis: max((abs(num(r, "expected_" + axis)) for r in candidates), default=0))


def observable(data, case):
    """Select two actual physical fields, retaining every run's original time."""
    body = selected_body(data)
    chosen = [(i + 2, row) for i, row in enumerate(data) if row.get("instance", "sphere") == body]
    vector = "expected_x" in data[0]
    spin_case = any(x in case for x in ("dem04", "dem05", "dem06", "dem07", "dem13", "dem15"))
    if vector:
        selected = [row for _, row in chosen]
        axis = component(selected, "run_vel")
        spin_axis = component(selected, "run_omega")
        specs = [("run_vel", axis, "velocity"), ("run_omega", spin_axis, "angular velocity")] if spin_case else [
            ("run_pos", axis, "position"), ("run_vel", axis, "velocity")]
        result = []
        for quantity, comp, name in specs:
            points = []
            for index, row in chosen:
                if row["quantity"] == quantity:
                    points.append({"row": index, "run": row["run"], "time": num(row, "time_s"),
                                   "reference": num(row, "expected_" + comp), "actual": num(row, "actual_" + comp)})
            result.append({"name": name, "component": comp, "referenceField": "expected_" + comp,
                           "actualField": "actual_" + comp, "quantityFilter": quantity, "points": points})
        return body, result
    fields = data[0]
    if "reference_official_y_m" in fields:
        specs = [("reference_official_y_m", "fundem_native_z_m", "position", "z"),
                 ("reference_official_vy_m_s", "fundem_native_vz_m_s", "velocity", "z")]
    else:
        axis = "x" if "reference_vx_m_s" in fields else "z"
        spin = next((c for c in "yzx" if f"reference_omega_{c}_rad_s" in fields), None)
        if spin and spin_case:
            specs = [(f"reference_v{axis}_m_s", f"fundem_v{axis}_m_s", "velocity", axis),
                     (f"reference_omega_{spin}_rad_s", f"fundem_omega_{spin}_rad_s", "angular velocity", spin)]
            if "dem07" in case or "dem13" in case:
                specs = [(f"reference_omega_{spin}_rad_s", f"fundem_omega_{spin}_rad_s", "angular velocity", spin),
                         (f"reference_torque_{spin}_N_m", f"fundem_torque_{spin}_N_m", "torque", spin)]
        else:
            specs = [(f"reference_{axis}_m", f"fundem_{axis}_m", "position", axis),
                     (f"reference_v{axis}_m_s", f"fundem_v{axis}_m_s", "velocity", axis)]
    result = []
    for reference, actual, name, comp in specs:
        result.append({"name": name, "component": comp, "referenceField": reference, "actualField": actual,
                       "points": [{"row": index, "run": row["run"], "time": num(row, "time_s"),
                                   "reference": num(row, reference), "actual": num(row, actual)} for index, row in chosen]})
    return body, result


def scales(points, residual=False):
    values = [p["actual"] - p["reference"] for p in points] if residual else [p[k] for p in points for k in ("actual", "reference")]
    maximum = max(map(abs, values), default=0)
    exponent = math.floor(math.log10(maximum)) if maximum and (maximum < .01 or maximum >= 10000) else 0
    return 10.0 ** exponent, exponent


def labels(spec, exponent, residual=False):
    symbol, unit = {"position": (spec["component"], "m"),
                    "velocity": (r"v_" + spec["component"], r"m\,s^{-1}"),
                    "angular velocity": (r"\omega_" + spec["component"], r"rad\,s^{-1}"),
                    "torque": (r"M_" + spec["component"], r"N\,m")}[spec["name"]]
    prefix = r"\Delta " if residual else ""
    scale = (r"10^{" + str(exponent) + r"}\," ) if exponent else ""
    return "$" + prefix + symbol + r"\;(\mathrm{" + "}" + scale + r"\mathrm{" + unit + r"})$"


def physical(axis, spec, time_factor, time_unit, residual=False):
    points = spec["points"]
    factor, exponent = scales(points, residual)
    if residual:
        axis.axhline(0, color="black", linewidth=.6, zorder=1)
    else:
        reference = [p for p in points if p["run"] == "cpu"]
        axis.plot([p["time"] / time_factor for p in reference], [p["reference"] / factor for p in reference],
                  linestyle="none", marker="D", markersize=4, markerfacecolor="none",
                  markeredgewidth=.75, color="black", zorder=2, clip_on=False)
    for i, run in enumerate(RUNS):
        selected = [p for p in points if p["run"] == run]
        marker(axis, [p["time"] / time_factor for p in selected],
               [((p["actual"] - p["reference"]) if residual else p["actual"]) / factor for p in selected], i)
    axis.set_xlabel(r"Time, $t$ (" + time_unit + ")")
    axis.set_ylabel(labels(spec, exponent, residual))
    axis.margins(x=.13, y=.23)
    if residual and all(p["actual"] == p["reference"] for p in points):
        axis.set_ylim(-1, 1)
    return {"kind": "signed residual" if residual else "physical checkpoints",
            "quantity": spec["name"], "component": spec["component"],
            "referenceField": spec["referenceField"], "actualField": spec["actualField"],
            "quantityFilter": spec.get("quantityFilter"), "sourceCsvRows": [p["row"] for p in points],
            "timeUnit": time_unit, "timeDivisor_s": time_factor, "ordinateDivisorSI": factor,
            "derivation": "(actual-reference)/ordinateDivisorSI" if residual else "CSV value/ordinateDivisorSI",
            "pointCount": len(points), "allResidualsExactlyZero": all(p["actual"] == p["reference"] for p in points)}


def metrics(data, case):
    result = []
    for index, row in enumerate(data):
        actual, expected = num(row, "actual"), num(row, "expected")
        error, tolerance = num(row, "error"), num(row, "tolerance")
        if tolerance <= 0 or error < 0 or row["error_kind"] not in ("relative", "absolute"):
            raise ValueError(f"Invalid criterion in {case}")
        computed = abs(actual - expected)
        if row["error_kind"] == "relative":
            if expected == 0:
                raise ValueError(f"Zero relative target: {case}")
            computed /= abs(expected)
        if not math.isclose(computed, error, rel_tol=1e-8, abs_tol=2e-14):
            raise ValueError(f"Stored criterion mismatch: {case}")
        if (error < tolerance) != (row["passed"].lower() == "true"):
            raise ValueError(f"Stored pass flag mismatch: {case}")
        result.append({**row, "row": index + 2, "ratio": computed / tolerance})
    return result


def short_quantity(quantity):
    for old, new in (("Overlap at official check time (m)", r"$\delta$"),
                     ("Force at official check time (N)", r"$F_n$"),
                     ("Hertz peak energy (J)", r"$E$"),
                     ("Total x momentum", r"$p_x$"),
                     ("Rebound vx (m/s)", r"$v_x$"), ("Rebound vz (m/s)", r"$v_z$"),
                     ("Sliding vx (m/s)", r"$v_x$"), ("Sliding omega_y (rad/s)", r"$\omega_y$"),
                     ("Preserved omega_y (rad/s)", r"$\omega_y$"), ("Preserved omega_z (rad/s)", r"$\omega_z$"),
                     ("Freefall z (m)", r"$z$"), ("Freefall vz (m/s)", r"$v_z$")):
        quantity = quantity.replace(old, new)
    quantity = quantity.replace(" at official check time", "").replace(" (m/s)", "").replace(" (rad/s)", "")
    # Convert source-code identifiers only outside already formatted math.
    import re
    parts = re.split(r"(\$[^$]*\$)", quantity)
    def symbol(match):
        return r"$\omega_" + match[1] + "$" if match[1] else "$v_" + match[2] + "$"
    return "".join(part if part.startswith("$") else re.sub(r"\b(?:omega_([xyz])|v([xyz]))\b", symbol, part) for part in parts)


def acceptance(axis, data):
    quantities = list(dict.fromkeys(row["quantity"] for row in data))
    ratios = [r["ratio"] for r in data]
    positive = [r for r in ratios if r > 0]
    floor = min(positive, default=1e-6) / 4
    linear = min(positive, default=0) > .05 and max(ratios) > 1
    axis.axvline(1, color="black", linestyle="--", linewidth=.8, zorder=1)
    for i, run in enumerate(RUNS):
        selected = [r for r in data if r["run"] == run]
        marker(axis, [r["ratio"] if linear else max(r["ratio"], floor) for r in selected],
               [quantities.index(r["quantity"]) + (i - 1.5) * .14 for r in selected], i)
    if linear:
        axis.set_xlim(0, max(1, *ratios) * 1.23)
        axis.xaxis.set_major_locator(MaxNLocator(nbins=5))
    else:
        axis.set_xscale("log")
        axis.set_xlim(floor / 2, max(1, *ratios) * 2.3)
    axis.set_yticks(range(len(quantities)), [short_quantity(q) for q in quantities])
    axis.tick_params(axis="y", length=0, pad=5)
    axis.set_ylim(len(quantities) - .55, -.55)
    axis.set_xlabel(r"Error / tolerance, $\varepsilon/\tau$ (dimensionless)")
    axis.grid(False)
    axis.grid(axis="x", which="major", zorder=0)
    return {"kind": "acceptance ratios", "quantities": quantities,
            "rows": [{k: r[k] for k in ("row", "run", "quantity", "expected", "actual", "error_kind", "tolerance", "ratio", "source", "sourceFields") if k in r} for r in data],
            "scale": "linear" if linear else "logarithmic",
            "displayFloor": None if linear else floor, "exactZeroCount": ratios.count(0),
            "zeroConvention": "Linear axis; zeros retain their true ordinate." if linear else "Only exact-zero ratios are placed at the stated positive display floor; all positive ratios retain their true values.",
            "boundary": 1, "maximumErrorToToleranceRatio": max(ratios),
            "derivation": "abs(actual-expected)/tolerance; divide additionally by abs(expected) for a relative criterion"}


def historical(path):
    return {row["field"]: row["value"] for row in read(path)}


def hertz_figure(root, record):
    """Use only measured checkpoint fields and retained global scalar summaries."""
    case = record["case"]
    folder = root / "cases" / case
    source = folder / "checkpoints/official_checkpoint_comparison.csv"
    scalar_source = folder / "HISTORICAL-RESULTS.csv"
    data, values = read(source), historical(scalar_source)
    scalar_rows = {r["field"]: i + 2 for i, r in enumerate(read(scalar_source))}
    sphere_pair = "twosphere" in case
    selected = [(i + 2, r) for i, r in enumerate(data) if not sphere_pair or r["instance"] == "left"]
    # Pair-level contact and energy columns are repeated for each sphere. Keep
    # one copy, rather than averaging or counting the same pair twice.
    if sphere_pair:
        for _, r in selected:
            other = [q for q in data if q["run"] == r["run"] and q["official_step"] == r["official_step"] and q["instance"] == "right"]
            if len(other) != 1 or any(num(r, k) != num(other[0], k) for k in ("overlap_m", "contact_force_N", "total_energy_J")):
                raise ValueError(f"Inconsistent duplicated pair fields: {case}")
    initial = float(values["derived.initial_energy_J"])
    peak_overlap = float(values["derived.analytic_maximum_overlap_m"])
    peak_force = float(values["derived.analytic_peak_force_N"])
    check_step = int(float(values["official_analytic_step"]))
    for _, r in selected:
        if not math.isclose(num(r, "total_energy_J"), num(r, "kinetic_energy_J") + num(r, "hertz_elastic_energy_J"), rel_tol=2e-13, abs_tol=1e-15):
            raise ValueError(f"Inconsistent retained total energy: {case}")
    energy_deviation = [(num(r, "total_energy_J") - initial) / initial for _, r in selected]
    energy_exponent = math.floor(math.log10(max(map(abs, energy_deviation)))) if any(energy_deviation) else 0
    energy_scale = 10.0 ** energy_exponent
    force_scale = 1000.0 if peak_force >= 1000 else 1.0
    fig, axes = canvas(reference="analytical")
    panels = []
    definitions = [
        ("a", "Hertz compression checkpoints", "overlap_m", 1e-6, peak_overlap, r"$\delta\;(\mu\mathrm{m})$"),
        ("b", "Hertz normal-force checkpoints", "contact_force_N", force_scale, peak_force, r"$F_n\;(\mathrm{kN})$" if force_scale == 1000 else r"$F_n\;(\mathrm{N})$"),
        ("c", "relative total-energy drift", "total_energy_J", energy_scale, 0, r"$(E_{\mathrm{tot}}-E_0)/E_0\;(10^{" + str(energy_exponent) + r"})$")]
    for ax, (letter, kind, field, factor, reference, ylabel) in zip(list(axes.flat)[:3], definitions):
        ax.axhline(reference / factor if letter != "c" else 0, color="black", linestyle="--", linewidth=.8, zorder=1)
        for i, run in enumerate(RUNS):
            rows = [r for _, r in selected if r["run"] == run]
            ys = [((num(r, field) - initial) / initial if letter == "c" else num(r, field)) / factor for r in rows]
            marker(ax, [num(r, "time_s") / 1e-6 for r in rows], ys, i)
        ax.set_xlabel(r"Time, $t$ ($\mu$s)")
        ax.set_ylabel(ylabel)
        ax.margins(x=.13, y=.23)
        panels.append({"letter": letter, "kind": kind, "field": field,
                       "source": source.relative_to(root).as_posix(), "sourceCsvRows": [i for i, _ in selected],
                       "timeDivisor_s": 1e-6, "ordinateDivisor": factor,
                       "reference": reference, "referenceKind": "conserved total energy" if letter == "c" else "continuum analytical peak",
                       "referenceField": "derived.initial_energy_J" if letter == "c" else ("derived.analytic_maximum_overlap_m" if letter == "a" else "derived.analytic_peak_force_N"),
                       "referenceSource": scalar_source.relative_to(root).as_posix(),
                       "initialEnergy_J": initial,
                       "derivation": "(total_energy_J-E0)/E0/ordinateDivisor" if letter == "c" else "CSV field/ordinateDivisor",
                       "pointCount": len(selected), "numericalCheckpointLines": "none"})
    criteria = []
    for i, run in enumerate(RUNS):
        prefix = f"runs[{i}]."
        if values[prefix + "run"] != run:
            raise ValueError(f"Historical run order mismatch: {case}")
        at_check = [(index, r) for index, r in selected if r["run"] == run and int(r["official_step"]) == check_step]
        if len(at_check) != 1:
            raise ValueError(f"Missing official Hertz check: {case}/{run}")
        index, row = at_check[0]
        official_energy = .4 * num(row, "contact_force_N") * num(row, "overlap_m")
        if not math.isclose(official_energy, float(values[prefix + "official_check_elastic_energy_from_force_J"]), rel_tol=2e-12, abs_tol=1e-15):
            raise ValueError(f"Historical force-overlap criterion mismatch: {case}/{run}")
        definitions = [
            (r"$U_H(t_c)$", official_energy, initial, "relative", float(values["official_analytic_tolerance"]), source.relative_to(root).as_posix(), index, ["contact_force_N", "overlap_m", "derived.initial_energy_J"]),
            (r"$\delta_{\mathrm{max}}$", float(values[prefix + "maximum_overlap_m"]), peak_overlap, "relative", float(values["supplemental_peak_relative_tolerance"]), scalar_source.relative_to(root).as_posix(), scalar_rows[prefix + "maximum_overlap_m"], [prefix + "maximum_overlap_m", "derived.analytic_maximum_overlap_m"]),
            (r"$F_{\mathrm{max}}$", float(values[prefix + "peak_force_N"]), peak_force, "relative", float(values["supplemental_peak_relative_tolerance"]), scalar_source.relative_to(root).as_posix(), scalar_rows[prefix + "peak_force_N"], [prefix + "peak_force_N", "derived.analytic_peak_force_N"]),
            (r"$\max|\Delta E/E_0|$", float(values[prefix + "maximum_relative_total_energy_error"]), 0, "absolute", float(values["supplemental_energy_relative_tolerance"]), scalar_source.relative_to(root).as_posix(), scalar_rows[prefix + "maximum_relative_total_energy_error"], [prefix + "maximum_relative_total_energy_error"]) ]
        for quantity, actual, expected, error_kind, tolerance, filename, row_number, fields in definitions:
            error = abs(actual - expected) / (abs(expected) if error_kind == "relative" else 1)
            criteria.append({"run": run, "quantity": quantity, "row": row_number, "actual": actual,
                             "expected": expected, "error_kind": error_kind, "tolerance": tolerance,
                             "ratio": error / tolerance, "source": filename, "sourceFields": fields})
    panels.append({"letter": "d", **acceptance(axes[1, 1], criteria)})
    caption = (f"Retained Hertz-contact measurements for {case}. Panels (a) and (b) show contact overlap and normal-force magnitude at the three retained physical checkpoint times for each run. "
               r"Their black dashed lines denote the continuum analytical maxima $\delta_{\mathrm{max}}$ and $F_{\mathrm{max}}$, not LAMMPS saved-state trajectories. "
               r"Panel (c) shows signed total-energy drift $(E_{\mathrm{tot}}-E_0)/E_0$, where $E_{\mathrm{tot}}=K+U_H$ and the dashed zero line denotes conservation of the retained initial energy $E_0$. "
               f"Here $E_0={initial:.12g}\\,\\mathrm{{J}}$. "
               r"Blue open circles, orange crosses, green open squares and purple plus signs denote CPU at $\Delta t$, GPU at $\Delta t$, CPU at $\Delta t/2$ and GPU at $\Delta t/2$, respectively. "
               "All numerical checkpoint observations remain unconnected symbols at their original times; overlapping runs are neither averaged nor shifted. ")
    if sphere_pair:
        caption += "Pair-level contact and total-energy columns are duplicated for the left and right spheres in the retained table; one identical copy is plotted, without summing the pair twice. "
    caption += (r"Panel (d) shows error divided by its stated limit. The $U_H(t_c)$ row is the official sampled criterion $U_H(t_c)=(2/5)F_n(t_c)\delta(t_c)$ against $E_0$, at "
                f"$t_c={float(values['official_analytic_time_s']) / 1e-6:.12g}\\,\\mu\\mathrm{{s}}$. "
                r"The $\delta_{\mathrm{max}}$ and $F_{\mathrm{max}}$ rows are separately declared global-peak diagnostics; the final row is the retained whole-event maximum $|\Delta E/E_0|$. "
                "Each uses a 0.01 limit; only the first row is the named official energy assertion. Global scalar summaries are read from HISTORICAL-RESULTS.csv, not inferred from the three checkpoints or assigned guessed peak times. "
                f"The dashed vertical acceptance boundary is error/limit = 1. Exactly {panels[-1]['exactZeroCount']} zero ratios are displayed at the explicit logarithmic floor {panels[-1]['displayFloor']:.12g}; positive ratios are unchanged. "
                "Run markers have fixed vertical offsets within each categorical diagnostic row. These historical measurements preserve the recorded outcome and do not constitute a new simulation or continuous-time validation.")
    return fig, {"case": case, "historicalOutcome": record["historical_outcome"], "caption": caption,
                 "width_mm": WIDTH_MM, "height_mm": 122.0,
                 "selectedBody": "left/right contact pair (one duplicated scalar copy)" if sphere_pair else "sphere-wall contact",
                 "source": source.relative_to(root).as_posix(), "sourceSha256": sha(source), "rows": len(data),
                 "sources": [{"file": p.relative_to(root).as_posix(), "sha256": sha(p)} for p in (source, scalar_source)],
                 "panels": panels, "textBoundsAndLegendLayoutChecked": True}


def endpoint_criteria(data, values):
    result = []
    for run in RUNS:
        row = max((r for r in data if r["run"] == run), key=lambda r: num(r, "time_s"))
        for axis, field, key in (("vx", "fundem_vx_m_s", "derived.analytic_final_vx_m_s"),
                                 ("omega_y", "fundem_omega_y_rad_s", "derived.analytic_final_omega_y_rad_s")):
            expected, actual = float(values[key]), num(row, field)
            error = abs(actual - expected) / abs(expected)
            tolerance = float(values["official_analytic_tolerance"])
            result.append({"run": run, "quantity": r"$v_x$" if axis == "vx" else r"$\omega_y$",
                           "row": data.index(row) + 2, "expected": expected, "actual": actual,
                           "error_kind": "relative", "tolerance": tolerance, "ratio": error / tolerance})
    return result


def lammps_figure(root, record):
    case = record["case"]
    if case in ("dem02-hertz-twosphere", "dem02-hertzmaterial-twosphere", "dem03-hertzmaterial-wall"):
        return hertz_figure(root, record)
    folder = root / "cases" / case
    source = folder / "checkpoints/official_checkpoint_comparison.csv"
    if not source.exists():
        source = folder / "checkpoints/source_checkpoint_comparison.csv"
    data = read(source)
    if not data or set(r["run"] for r in data) != set(RUNS):
        raise ValueError(f"Missing historical runs: {case}")
    sources = [source]
    body, observables = observable(data, case)
    all_times = [p["time"] for spec in observables for p in spec["points"]]
    maximum_time = max(all_times)
    time_factor, time_unit = (1e-6, r"$\mu$s") if maximum_time <= .001 else ((.001, "ms") if maximum_time <= .1 else (1, "s"))
    fig, axes = canvas()
    panels = [physical(axes[0, i], spec, time_factor, time_unit) for i, spec in enumerate(observables)]
    panels.append(physical(axes[1, 0], observables[0], time_factor, time_unit, True))
    criterion_path = folder / "checkpoints/official_metrics_comparison.csv"
    criteria = None
    endpoint = case in ("dem05-mindlin-slip", "dem05-slip-mu02")
    if criterion_path.exists():
        sources.append(criterion_path)
        criteria = metrics(read(criterion_path), case)
    elif endpoint:
        values_path = folder / "HISTORICAL-RESULTS.csv"
        sources.append(values_path)
        criteria = endpoint_criteria(data, historical(values_path))
    if criteria:
        panels.append(acceptance(axes[1, 1], criteria))
    else:
        panels.append(physical(axes[1, 1], observables[1], time_factor, time_unit, True))
    for i, panel in enumerate(panels):
        panel["letter"] = chr(97 + i)
    caption = (f"Historical checkpoint comparison for {case}. Panels (a) and (b) show {observables[0]['name']} "
               f"({observables[0]['component']} component) and {observables[1]['name']} ({observables[1]['component']} component) "
               f"of the selected body '{body}'; (c) shows the signed FunDEM-minus-reference residual of (a). "
               "Black open diamonds denote pinned LAMMPS saved-state references; blue open circles, orange crosses, "
               "green open squares and purple plus signs denote CPU at Δt, GPU at Δt, CPU at Δt/2 and GPU at Δt/2, respectively. "
               "Every run uses its retained physical checkpoint times; no time shift, jitter, averaging, fitting or connecting trajectory is used. "
               "Coincident symbols are retained with different shapes and sizes. Axes state physical units and explicit powers of ten. ")
    if criteria:
        definition = "; ".join(f"{q}: " + ", ".join(sorted({f"{r['error_kind']} tolerance {float(r['tolerance']):.12g}" for r in criteria if r['quantity'] == q})) for q in dict.fromkeys(r["quantity"] for r in criteria))
        caption += ("Panel (d) shows recomputed error divided by its stated tolerance; the dashed vertical line is the acceptance boundary ε/τ = 1. "
                    "Absolute errors are |actual−expected|; relative errors divide this by |expected|. " + definition + ". ")
        floor = panels[-1]["displayFloor"]
        if floor is not None:
            caption += f"Exactly {panels[-1]['exactZeroCount']} zero ratios are displayed at the explicit logarithmic floor {floor:.12g}; every positive ratio is unmodified. "
        else:
            caption += "Panel (d) uses a linear error-ratio axis; zeros retain their true value. "
        caption += "The four run markers have deterministic vertical offsets within each categorical criterion row solely to distinguish coincident values. "
        if endpoint:
            caption += "The endpoint targets in (d) are the official analytical rolling limits retained in HISTORICAL-RESULTS.csv, not the LAMMPS saved-state ordinates in (a,b); the relative tolerance is 0.01. The historical FAIL remains visible above the acceptance boundary. "
    else:
        caption += "Panel (d) shows the signed residual of (b). Both residual panels use linear scales, with zero shown at its true ordinate. "
    if case == "dem10-stack":
        caption += "The source's vertical y coordinate maps to native z; the selected body's source and native ordinates retain this documented mapping. "
    if case == "dem01-hooke-2d-lj":
        caption = caption.replace("Axes state physical units", "This reduced-unit source retains its original numerical units; the inherited SI-style field suffixes do not assert a physical SI calibration. Axes state source coordinate units")
        # The source is genuinely in LJ reduced units: replace labels explicitly.
        for ax in axes.flat:
            ax.set_xlabel(r"Time, $t$ (reduced units)")
            ax.set_ylabel(ax.get_ylabel().replace(r"\mathrm{m\,s^{-1}}", r"\mathrm{reduced\ velocity}").replace(r"\mathrm{m}", r"\mathrm{reduced\ length}"))
        for panel in panels:
            if "timeUnit" in panel:
                panel["timeUnit"] = "reduced source time"
    caption += "These sparse historical checkpoints do not establish a new validation outcome for the current source or a continuous-time error bound."
    return fig, {"case": case, "historicalOutcome": record["historical_outcome"], "caption": caption,
                 "width_mm": WIDTH_MM, "height_mm": 122.0, "selectedBody": body,
                 "source": source.relative_to(root).as_posix(), "sourceSha256": sha(source),
                 "rows": len(data), "sources": [{"file": p.relative_to(root).as_posix(), "sha256": sha(p)} for p in sources],
                 "panels": panels, "textBoundsAndLegendLayoutChecked": True}


def pfc_figure(root, record):
    case = record["case"]
    folder = root / "cases" / case
    sweep = case.endswith("sweep")
    source = folder / "checkpoints" / ("sixteen_beta_comparison.csv" if sweep else "first_rebound_comparison.csv")
    data = read(source)
    expected_rows = 64 if sweep else 4
    if len(data) != expected_rows or set(r["run"] for r in data) != set(RUNS):
        raise ValueError(f"Incomplete PFC checkpoint table: {case}")
    fig, axes = canvas(height_mm=92.0, rows=1, cols=3, reference=False)
    fig.subplots_adjust(left=.084, right=.982, bottom=.20, top=.82, wspace=.48)
    beta_single = .2 if "beta02" in case else .6
    maximum_relative = 0.0
    for i, run in enumerate(RUNS):
        selected = sorted((r for r in data if r["run"] == run), key=lambda r: num(r, "beta") if sweep else 0)
        beta = [num(r, "beta") if sweep else beta_single for r in selected]
        target = [num(r, "analytic_restitution") if sweep else math.exp(-math.pi * b / math.sqrt(1 - b * b)) for r, b in zip(selected, beta)]
        measured = [num(r, "measured_restitution") for r in selected]
        errors = [num(r, "absolute_restitution_error") for r in selected]
        if any(not math.isclose(abs(m - a), e, rel_tol=1e-9, abs_tol=2e-14) for m, a, e in zip(measured, target, errors)):
            raise ValueError(f"PFC analytical error mismatch: {case}")
        x = beta if sweep else [i]
        marker(axes[0, 0], x, measured, i)
        marker(axes[0, 1], x, errors, i)
        relative = [100 * e / abs(a) if a != 0 else None for e, a in zip(errors, target)]
        marker(axes[0, 2], [xx for xx, rr in zip(x, relative) if rr is not None], [rr for rr in relative if rr is not None], i)
        maximum_relative = max(maximum_relative, *(r for r in relative if r is not None))
    if sweep:
        dense_beta = np.linspace(0, 1, 1001)
        analytic = np.zeros_like(dense_beta)
        analytic[:-1] = np.exp(-math.pi * dense_beta[:-1] / np.sqrt(1 - dense_beta[:-1] ** 2))
        axes[0, 0].plot(dense_beta, analytic, color="black", linewidth=.9, zorder=1)
    else:
        target = math.exp(-math.pi * beta_single / math.sqrt(1 - beta_single * beta_single))
        axes[0, 0].axhline(target, color="black", linewidth=.9, zorder=1)
    axes[0, 1].axhline(.01, color="black", linestyle="--", linewidth=.8)
    for ax in axes.flat:
        if sweep:
            ax.set_xlabel(r"Damping ratio, $\beta$ (dimensionless)")
            ax.set_xlim(-.04, 1.04)
            ax.set_xticks([0, .25, .5, .75, 1])
        else:
            ax.set_xticks(range(4), ["CPU", "GPU", "CPU", "GPU"])
            ax.set_xlim(-.45, 3.45)
            ax.set_xlabel(r"$\Delta t$                         $\Delta t/2$")
    axes[0, 0].set_ylabel(r"$e_h$ (dimensionless)")
    axes[0, 1].set_ylabel(r"$|e_h-e_0|$ (dimensionless)")
    axes[0, 2].set_ylabel(r"$100|e_h-e_0|/e_0$ (%)")
    axes[0, 0].set_ylim(-.035, 1.045) if sweep else axes[0, 0].margins(y=.3)
    axes[0, 1].set_ylim(-.00035, .0111)
    axes[0, 2].set_ylim(-.025 * max(1, maximum_relative), 1.17 * max(1, maximum_relative))
    caption = (f"Historical first-rebound restitution for {case}. Panel (a) compares the retained first-apex height estimator "
               r"$e_h=\sqrt{\max(z_{\mathrm{apex}}-R,0)/(z_{\mathrm{initial}}-R)}$ with the public PFC $M_d=0$ contact target "
               r"$e_0=\exp[-\pi\beta/\sqrt{1-\beta^2}]$ for $0\leq\beta<1$, with $e_0=0$ at $\beta=1$ (black solid analytical curve or horizontal line). "
               "Panel (b) shows absolute coefficient errors and the strict published absolute-error limit 0.01 (black dashed line). "
               "Panel (c) shows relative discrepancy as a separate diagnostic, not the acceptance criterion. "
               r"Blue open circles, orange crosses, green open squares and purple plus signs denote CPU at $\Delta t$, GPU at $\Delta t$, "
               r"CPU at $\Delta t/2$ and GPU at $\Delta t/2$, respectively; all numerical observations are unconnected symbols, "
               r"and coincident runs remain overlaid without averaging or jitter. $\Delta t=10^{-5}\,\mathrm{s}$; $\Delta t/2=5\times10^{-6}\,\mathrm{s}$. ")
    if sweep:
        caption += (r"Sixteen independent damping ratios are retained per run. Relative discrepancy is undefined at $e_0=0$ ($\beta=1$) "
                    f"and that point is omitted only from (c); the largest plotted relative discrepancy is {maximum_relative:.3f}%. "
                    r"At $\beta=0.9$ and $\beta=1$, the retained $e_h=0$ denotes non-rebound. The recorded extremum at $\beta=0.9$ is below the nominal radius and is clamped to zero by the height estimator; $\beta=1$ has no recorded apex, so its retained zero is a non-rebound convention rather than an observed apex-height measurement. "
                    "The black analytical curve is an evaluation of the stated closed-form equation, not interpolation of the numerical observations. ")
    else:
        caption += rf"Only four summary observations exist at $\beta={beta_single:g}$; the horizontal categories do not form a time history. "
    caption += ("All three ordinates are linear and exact zeros retain their true value. PFC was not executed locally: "
                "the reference is its published analytical relation, not a PFC numerical trajectory. The historical PASS records absolute-error acceptance only and does not imply small relative error at strong damping.")
    panels = [
        {"letter": "a", "kind": "restitution", "actualField": "measured_restitution", "reference": "exp(-pi*beta/sqrt(1-beta^2)), endpoint 0 at beta=1", "analyticalCurveSamples": 1001 if sweep else 1},
        {"letter": "b", "kind": "absolute error", "field": "absolute_restitution_error", "boundary": .01},
        {"letter": "c", "kind": "relative discrepancy", "derivation": "100*absolute_restitution_error/analytic_restitution", "undefinedAt": "analytic_restitution == 0", "maximumPercent": maximum_relative},
    ]
    for panel in panels:
        panel["sourceCsvRows"] = list(range(2, len(data) + 2))
        panel["scale"] = "linear"
    return fig, {"case": case, "historicalOutcome": record["historical_outcome"], "caption": caption,
                 "width_mm": WIDTH_MM, "height_mm": 92.0, "source": source.relative_to(root).as_posix(),
                 "sourceSha256": sha(source), "rows": len(data),
                 "sources": [{"file": source.relative_to(root).as_posix(), "sha256": sha(source)}],
                 "panels": panels, "textBoundsAndLegendLayoutChecked": True}


def generate(root, family, check=False):
    settings()
    root = Path(root).resolve()
    records = read(root / "HISTORICAL-INDEX.csv")
    expected_count = 35 if family == "LAMMPS" else 3
    if len(records) != expected_count:
        raise ValueError(f"Expected {expected_count} cases")
    font_path = Path(font_manager.findfont("STIXGeneral", fallback_to_default=False))
    manifest = {"schemaVersion": 2, "kind": "Publication presentation of immutable historical checkpoints; no new simulation",
                "presentationDate": "2026-10-31", "family": family,
                "style": {"width_mm": WIDTH_MM, "font": "STIXGeneral", "mathFont": "STIX",
                          "fontSize_pt": 8.5, "axisLabelSize_pt": 9, "tickSize_pt": 8,
                          "panelLetterSize_pt": 10, "lineWidth_pt": .8, "pngDpi": 600,
                          "pdfFontType": 42, "svgText": "editable text; installed STIX font needed for identical display",
                          "colors": list(COLORS), "markers": list(MARKERS),
                          "markerSizes_pt": list(SIZES), "checkpointLines": "none", "jitter": "none"},
                "renderer": {"matplotlib": matplotlib.__version__, "python": platform.python_version(),
                             "fontFile": font_path.name, "fontSha256": sha(font_path),
                             "deterministicMetadata": "No timestamps; fixed SVG hash salt and creator strings"},
                "cases": []}
    for record in records:
        fig, item = (lammps_figure if family == "LAMMPS" else pfc_figure)(root, record)
        before = {entry["file"]: entry["sha256"] for entry in item["sources"]}
        output = root / "cases" / record["case"] / "figures"
        if not check:
            output.mkdir(exist_ok=True)
        artifacts = []
        for fmt, content in render(fig, record["case"]).items():
            path = output / ("comparison." + fmt)
            if check:
                if path.read_bytes() != content:
                    raise ValueError(f"Non-reproducible figure bytes: {path}")
            else:
                path.write_bytes(content)
            artifacts.append({"file": path.relative_to(root).as_posix(), "sha256": hashlib.sha256(content).hexdigest(),
                              "format": fmt, "dpi": 600 if fmt == "png" else "vector"})
        plt.close(fig)
        for filename, digest in before.items():
            if sha(root / filename) != digest:
                raise ValueError(f"Immutable checkpoint table changed: {filename}")
        item["artifacts"] = artifacts
        manifest["cases"].append(item)
    encoded = json.dumps(manifest, indent=2, ensure_ascii=False) + "\n"
    destination = root / "FIGURES.json"
    if check:
        if destination.read_text(encoding="utf-8") != encoded:
            raise ValueError(f"Non-reproducible figure manifest: {destination}")
    else:
        destination.write_text(encoded, encoding="utf-8")
    print(f"FIGURES_{'CHECK_PASSED' if check else 'WRITTEN'}: {expected_count} {family} cases, {expected_count * 3} SVG/PDF/600-dpi PNG files; source CSV hashes unchanged")
