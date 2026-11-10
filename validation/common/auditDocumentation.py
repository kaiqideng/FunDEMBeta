#!/usr/bin/env python3
"""Check publication reports, immutable historical data, figures and current checksums.

This standard-library-only audit checks presentation and evidence integrity. It
does not execute a solver or confer a new physical-verification outcome.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import struct
from urllib.parse import unquote, urlsplit
import xml.etree.ElementTree as ET
import zlib

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = {"LAMMPSComparison": (35, 33, 2), "PFCComparison": (3, 3, 0)}
LEDGER = "EVIDENCE-SHA256SUMS.txt"
OLD_LEDGER = "HISTORICAL-EVIDENCE-SHA256SUMS.txt"
SECTIONS = (
    "Abstract",
    "1. Benchmark definition",
    "2. Numerical method and verification procedure",
    "3. Results and discussion",
    "4. Conclusions and limitations",
    "References",
    "Appendix A. Reproducibility and evidence",
)
FIGURE_FORMATS = {"svg", "pdf", "png"}
WIDTH_MM = 180.0
PNG_DPI = 600


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checksum_rows(path):
    records = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if not match:
            raise ValueError(f"Invalid checksum row: {path}")
        name = match[2]
        if Path(name).is_absolute() or ".." in Path(name).parts or name in records:
            raise ValueError(f"Unsafe or duplicated checksum path: {name}")
        records[name] = match[1]
    return records


def public_files(folder):
    return sorted(path for path in folder.rglob("*") if path.is_file()
                  and path.name != LEDGER and "__pycache__" not in path.parts and not path.is_symlink())


def write_ledgers():
    for family in FAMILIES:
        folder = ROOT / family
        lines = [f"{sha(path)}  {path.relative_to(folder).as_posix()}\n" for path in public_files(folder)]
        (folder / LEDGER).write_text("".join(lines), encoding="utf-8")


def finite_positive(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value) and value > 0


def local_manifest_file(folder, name):
    """Resolve only an existing regular file inside its manifest family."""
    if not isinstance(name, str) or not name or "\\" in name:
        raise ValueError(f"Invalid manifest path: {name!r}")
    path = Path(name)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"Unsafe manifest path: {name}")
    result = folder / path
    if not result.is_file() or result.is_symlink() or not result.resolve().is_relative_to(folder.resolve()):
        raise ValueError(f"Missing, linked or escaped manifest file: {name}")
    return result


def png_metadata(raw):
    """Read dimensions and physical resolution, checking PNG chunk structure/CRCs."""
    if not raw.startswith(b"\x89PNG\r\n\x1a\n"):
        raise ValueError("Invalid PNG signature")
    position = 8
    dimensions = resolution = None
    seen_idat = ended = False
    while position < len(raw):
        if position + 12 > len(raw):
            raise ValueError("Truncated PNG chunk")
        size = struct.unpack_from(">I", raw, position)[0]
        kind = raw[position + 4:position + 8]
        end = position + 12 + size
        if end > len(raw):
            raise ValueError("Truncated PNG chunk data")
        payload = raw[position + 8:position + 8 + size]
        expected_crc = struct.unpack_from(">I", raw, position + 8 + size)[0]
        if zlib.crc32(kind + payload) & 0xffffffff != expected_crc:
            raise ValueError(f"PNG CRC mismatch: {kind!r}")
        if kind == b"IHDR":
            if dimensions is not None or position != 8 or size != 13:
                raise ValueError("Invalid or repeated PNG IHDR")
            dimensions = struct.unpack(">II", payload[:8])
            if not all(dimensions):
                raise ValueError("Zero PNG dimensions")
        elif kind == b"pHYs":
            if resolution is not None or size != 9:
                raise ValueError("Invalid or repeated PNG pHYs")
            resolution = struct.unpack(">IIB", payload)
        elif kind == b"IDAT":
            seen_idat = True
        elif kind == b"IEND":
            if size or end != len(raw):
                raise ValueError("Invalid PNG IEND or trailing data")
            ended = True
            break
        position = end
    if dimensions is None or resolution is None or not seen_idat or not ended:
        raise ValueError("PNG missing dimensions, physical resolution, image data or end marker")
    x_ppm, y_ppm, unit = resolution
    if unit != 1 or not x_ppm or not y_ppm:
        raise ValueError("PNG pHYs must specify nonzero pixels per metre")
    return {"pixels": dimensions, "dpi": (x_ppm * 0.0254, y_ppm * 0.0254),
            "millimetres": (1000 * dimensions[0] / x_ppm, 1000 * dimensions[1] / y_ppm)}


def svg_millimetres(value):
    match = re.fullmatch(r"\s*([0-9]+(?:\.[0-9]*)?|\.[0-9]+)\s*(mm|cm|in|pt|px)?\s*", value or "")
    if not match:
        raise ValueError(f"Invalid SVG physical length: {value!r}")
    factor = {"mm": 1, "cm": 10, "in": 25.4, "pt": 25.4 / 72,
              "px": 25.4 / 96, None: 25.4 / 96}[match[2]]
    result = float(match[1]) * factor
    if not finite_positive(result):
        raise ValueError("Invalid SVG dimension")
    return result


def check_figure(path, width_mm, height_mm):
    """Validate standalone dimensions and static/physical metadata in each format."""
    suffix = path.suffix.lower()
    if suffix == ".svg":
        text = path.read_text(encoding="utf-8")
        root = ET.fromstring(text)
        if root.tag.rsplit("}", 1)[-1] != "svg":
            raise ValueError("Invalid SVG root")
        forbidden = {"script", "foreignObject", "animate", "animateMotion", "animateTransform", "set", "discard"}
        for element in root.iter():
            if element.tag.rsplit("}", 1)[-1] in forbidden:
                raise ValueError("Non-static SVG element")
            for attribute, value in element.attrib.items():
                name = attribute.rsplit("}", 1)[-1].lower()
                if name.startswith("on") or (name == "href" and not value.startswith("#")):
                    raise ValueError("Scripted or externally referenced SVG")
                if name == "style" and re.search(r"animation\s*:|@import|url\(\s*['\"]?(?:https?:|//|javascript:)", value, re.I):
                    raise ValueError("Dynamic or externally referenced SVG style")
            if element.tag.rsplit("}", 1)[-1] == "style" and re.search(r"animation\s*:|@import|url\(\s*['\"]?(?:https?:|//|javascript:)", element.text or "", re.I):
                raise ValueError("Dynamic or externally referenced SVG stylesheet")
        dimensions = (svg_millimetres(root.get("width")), svg_millimetres(root.get("height")))
        viewbox = [float(value) for value in (root.get("viewBox") or "").replace(",", " ").split()]
        if len(viewbox) != 4 or not all(math.isfinite(v) for v in viewbox) or viewbox[2] <= 0 or viewbox[3] <= 0:
            raise ValueError("Invalid SVG viewBox")
        if not math.isclose(dimensions[0] / dimensions[1], viewbox[2] / viewbox[3], rel_tol=1e-6):
            raise ValueError("SVG physical/viewBox aspect-ratio mismatch")
    elif suffix == ".png":
        metadata = png_metadata(path.read_bytes())
        if any(abs(dpi - PNG_DPI) > 0.02 for dpi in metadata["dpi"]):
            raise ValueError(f"PNG resolution is not {PNG_DPI} dpi: {metadata['dpi']}")
        if any(abs(pixels - mm * PNG_DPI / 25.4) > 1.01
               for pixels, mm in zip(metadata["pixels"], (width_mm, height_mm))):
            raise ValueError(f"PNG pixel dimensions disagree with publication size: {metadata['pixels']}")
        dimensions = metadata["millimetres"]
    elif suffix == ".pdf":
        raw = path.read_bytes()
        if not re.match(rb"%PDF-1\.[0-9]", raw) or b"%%EOF" not in raw[-1024:]:
            raise ValueError("Invalid PDF header or end marker")
        number = rb"([-+0-9.eE]+)"
        match = re.search(rb"/MediaBox\s*\[\s*" + rb"\s+".join([number] * 4) + rb"\s*\]", raw)
        if not match:
            raise ValueError("PDF physical page dimensions unavailable")
        x0, y0, x1, y1 = (float(value) for value in match.groups())
        dimensions = ((x1 - x0) * 25.4 / 72, (y1 - y0) * 25.4 / 72)
    else:
        raise ValueError(f"Unsupported figure format: {suffix}")
    # One raster pixel at 600 dpi is 0.0423 mm; metadata rounding needs <0.05 mm.
    tolerance = 0.05 if suffix == ".png" else 0.02
    if any(not math.isfinite(actual) or abs(actual - expected) > tolerance
           for actual, expected in zip(dimensions, (width_mm, height_mm))):
        raise ValueError(f"Figure size {dimensions} mm disagrees with manifest {(width_mm, height_mm)} mm")


def report_sections(text, label, failures):
    headings = re.findall(r"^## ([^\n]+)\s*$", text, flags=re.M)
    positions = []
    for heading in SECTIONS:
        if headings.count(heading) != 1:
            failures.append(f"Missing or duplicated publication section {heading}: {label}")
        else:
            positions.append(headings.index(heading))
    if positions != sorted(positions):
        failures.append(f"Publication sections out of order: {label}")


def figure_manifest(folder, cases, failures, counts):
    label = folder.name
    try:
        manifest = json.loads((folder / "FIGURES.json").read_text(encoding="utf-8"))
        expected_family = "LAMMPS" if label == "LAMMPSComparison" else "PFC"
        if manifest.get("schemaVersion") != 2 or manifest.get("family") != expected_family:
            raise ValueError("Expected schemaVersion 2 and matching family")
        style = manifest.get("style", {})
        if style.get("width_mm") != WIDTH_MM or style.get("pngDpi") != PNG_DPI:
            raise ValueError("Expected publication style width 180 mm and PNG resolution 600 dpi")
        entries = manifest.get("cases")
        if not isinstance(entries, list) or any(not isinstance(entry, dict) for entry in entries):
            raise ValueError("Invalid manifest cases array")
        by_case = {entry.get("case"): entry for entry in entries}
        if len(by_case) != len(entries) or set(by_case) != {case["case"] for case in cases}:
            raise ValueError("Manifest case inventory differs from historical index")
    except (OSError, ValueError, TypeError) as exc:
        failures.append(f"Invalid figure manifest: {label}: {exc}")
        return {}
    counts["figureManifests"] += 1
    return by_case


def audit_case_figure(folder, case, entry, summary, failures, counts):
    label = f"{folder.name}/{case['case']}"
    try:
        if entry.get("historicalOutcome") != case["historical_outcome"]:
            raise ValueError("Figure historical outcome differs from frozen index")
        width, height = entry.get("width_mm"), entry.get("height_mm")
        if width != WIDTH_MM or not finite_positive(height):
            raise ValueError("Figure requires width 180 mm and a positive finite height")
        caption = entry.get("caption")
        if not isinstance(caption, str) or not caption.strip():
            raise ValueError("Missing figure caption")
        captions = re.findall(r"^\*\*Figure 1\.\*\*\s*([^\n]+)", summary, flags=re.M)
        if captions != [caption]:
            raise ValueError("Report Figure 1 caption does not exactly match manifest")
        counts["matchedFigureCaptions"] += 1
        sources = entry.get("sources")
        if not isinstance(sources, list) or not sources:
            raise ValueError("Missing figure sources")
        source_names = set()
        primary_rows = None
        for source in sources:
            if not isinstance(source, dict):
                raise ValueError("Invalid figure source record")
            name = source.get("file")
            path = local_manifest_file(folder, name)
            # Some official analytical endpoint targets live in the immutable
            # historical table rather than the saved-state checkpoint table.
            allowed_source = (name.startswith(f"cases/{case['case']}/checkpoints/")
                              or name == f"cases/{case['case']}/HISTORICAL-RESULTS.csv")
            if name in source_names or not allowed_source or path.suffix != ".csv":
                raise ValueError(f"Invalid or duplicate checkpoint source: {name}")
            source_names.add(name)
            if source.get("sha256") != sha(path):
                raise ValueError(f"Figure source hash mismatch: {name}")
            counts["figureSourceHashes"] += 1
            if name == entry.get("source"):
                if source.get("sha256") != entry.get("sourceSha256"):
                    raise ValueError("Primary figure source hash inconsistent")
                with path.open(encoding="utf-8-sig", newline="") as stream:
                    primary_rows = sum(1 for _ in csv.DictReader(stream))
        if primary_rows is None or entry.get("rows") != primary_rows:
            raise ValueError("Primary source missing or CSV row count differs from manifest")
        panels = entry.get("panels")
        if not isinstance(panels, list) or not panels or any(not isinstance(panel, dict) for panel in panels):
            raise ValueError("Invalid figure panels")
        letters = [panel.get("letter") for panel in panels]
        if letters != [chr(ord("a") + i) for i in range(len(panels))]:
            raise ValueError("Figure panel letters must be ordered and unique")
        for panel in panels:
            if not isinstance(panel.get("kind"), str) or not panel["kind"]:
                raise ValueError("Figure panel kind missing")
            if "sourceCsvRows" in panel and (not isinstance(panel["sourceCsvRows"], list)
                    or not panel["sourceCsvRows"] or any(not isinstance(row, int) or isinstance(row, bool)
                    or not 2 <= row <= primary_rows + 1 for row in panel["sourceCsvRows"])):
                raise ValueError("Figure panel source row outside primary checkpoint table")
            counts["figurePanels"] += 1
        if entry.get("textBoundsAndLegendLayoutChecked") is not True:
            raise ValueError("Figure manifest lacks completed layout-check declaration")
        artifacts = entry.get("artifacts")
        if not isinstance(artifacts, list) or len(artifacts) != len(FIGURE_FORMATS):
            raise ValueError("Expected exactly three figure artifacts")
        if {artifact.get("format") for artifact in artifacts if isinstance(artifact, dict)} != FIGURE_FORMATS:
            raise ValueError("Figure artifacts must include SVG, PDF and PNG once each")
        for artifact in artifacts:
            suffix = artifact["format"]
            name = f"cases/{case['case']}/figures/comparison.{suffix}"
            if artifact.get("file") != name or artifact.get("dpi") != (PNG_DPI if suffix == "png" else "vector"):
                raise ValueError(f"Figure artifact path or resolution metadata invalid: {suffix}")
            path = local_manifest_file(folder, name)
            if artifact.get("sha256") != sha(path):
                raise ValueError(f"Figure artifact hash mismatch: {name}")
            check_figure(path, width, height)
            counts["figures"] += 1
            counts["figureArtifactHashes"] += 1
    except (OSError, ValueError, TypeError, KeyError, ET.ParseError, struct.error, OverflowError) as exc:
        failures.append(f"Invalid publication figure: {label}: {exc}")


def audit():
    failures = []
    historical_ledger_differences = []
    counts = {"markdownFiles": 0, "localLinks": 0, "currentChecksums": 0,
              "frozenCsvFiles": 0, "archivedReports": 0, "figures": 0, "cases": 0,
              "publicationReports": 0, "currentReportHashes": 0, "figureManifests": 0,
              "matchedFigureCaptions": 0, "figureSourceHashes": 0, "figureArtifactHashes": 0,
              "figurePanels": 0}
    for document in sorted(ROOT.rglob("*.md")):
        text = document.read_text(encoding="utf-8")
        counts["markdownFiles"] += 1
        if re.search(r"[\u3400-\u9fff]", text):
            failures.append(f"Non-English CJK text: {document.relative_to(ROOT)}")
        prose = re.sub(r"^```[^\n]*\n.*?^```[^\n]*$", "", text, flags=re.M | re.S)
        for address in re.findall(r"!?\[[^\]\n]*\]\(([^)\n]+)\)", prose):
            address = address.strip().strip("<>")
            uri = urlsplit(address)
            if uri.scheme or not uri.path:
                continue
            counts["localLinks"] += 1
            target = document.parent / unquote(uri.path)
            if not target.exists():
                failures.append(f"Missing link: {document.relative_to(ROOT)} -> {address}")
    for family, (expected_cases, passes, fails) in FAMILIES.items():
        folder = ROOT / family
        historical = checksum_rows(folder / OLD_LEDGER)
        for name, expected in historical.items():
            if not name.endswith(".csv"):
                continue
            path = folder / name
            if not path.is_file() or sha(path) != expected:
                failures.append(f"Changed historical CSV: {family}/{name}")
            counts["frozenCsvFiles"] += 1
        current = checksum_rows(folder / LEDGER)
        expected_names = {path.relative_to(folder).as_posix() for path in public_files(folder)}
        if set(current) != expected_names:
            failures.append(f"Current checksum inventory mismatch: {family}")
        for name, expected in current.items():
            path = folder / name
            if not path.is_file() or path.is_symlink() or sha(path) != expected:
                failures.append(f"Changed distributed file: {family}/{name}")
            counts["currentChecksums"] += 1
        with (folder / "HISTORICAL-INDEX.csv").open(encoding="utf-8-sig", newline="") as stream:
            cases = list(csv.DictReader(stream))
        outcomes = [case["historical_outcome"] for case in cases]
        if (len(cases), outcomes.count("PASS"), outcomes.count("FAIL")) != (expected_cases, passes, fails):
            failures.append(f"Historical case outcomes changed: {family}")
        manifest_entries = figure_manifest(folder, cases, failures, counts)
        for case in cases:
            label = f"{family}/{case['case']}"
            case_root = folder / "cases" / case["case"]
            archived = case_root / "HISTORICAL-REPORT.md"
            provenance = (case_root / "PROVENANCE.md").read_text(encoding="utf-8")
            archived_hash = re.search(r"Archived report SHA256 \(`?HISTORICAL-REPORT\.md`?\): `([0-9a-f]{64})`", provenance)
            expected = archived_hash[1] if archived_hash else ""
            if not archived.is_file() or sha(archived) != expected:
                failures.append(f"Archived report changed: {label}")
            ledger_hash = historical[f"cases/{case['case']}/REPORT.md"]
            if ledger_hash != expected:
                historical_ledger_differences.append({"case": case["case"], "family": family,
                    "oldLedgerReportSha256": ledger_hash, "archivedReportSha256": expected})
            counts["archivedReports"] += 1
            report = case_root / "REPORT.md"
            summary = report.read_text(encoding="utf-8")
            report_sections(summary, label, failures)
            counts["publicationReports"] += 1
            outcome_pattern = r"\*\*[^*\n]*\b" + re.escape(case["historical_outcome"]) + r"\b[^*\n]*\*\*"
            if not re.search(outcome_pattern, summary):
                failures.append(f"Historical outcome missing from report: {label}")
            current_hash = re.search(r"Current report SHA256(?: \(`?REPORT\.md`?\))?: `([0-9a-f]{64})`", provenance)
            if not current_hash or current_hash[1] != sha(report):
                failures.append(f"Current report hash missing or mismatched: {label}")
            counts["currentReportHashes"] += 1
            if case["case"] in manifest_entries:
                audit_case_figure(folder, case, manifest_entries[case["case"]], summary, failures, counts)
            counts["cases"] += 1
    ondem = ROOT / "ONDEMComparison/shearFlow/README.md"
    if ondem.is_file():
        report_sections(ondem.read_text(encoding="utf-8"), str(ondem.relative_to(ROOT)), failures)
        counts["publicationReports"] += 1
    return {"checks": counts, "failures": failures,
            "preexistingHistoricalLedgerReportDifferences": historical_ledger_differences,
            "scope": "Publication structure, preserved data, source hashes and standalone figure integrity only; no physical rerun"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-manifests", action="store_true", help="Explicitly refresh current presentation checksums; never change historical checksums")
    parser.add_argument("--report", type=Path, help="Optional external JSON audit report")
    args = parser.parse_args()
    if args.write_manifests:
        write_ledgers()
    result = audit()
    if args.report:
        args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    display = {key: value for key, value in result.items() if key != "preexistingHistoricalLedgerReportDifferences"}
    display["preexistingHistoricalLedgerReportDifferenceCount"] = len(result["preexistingHistoricalLedgerReportDifferences"])
    print(json.dumps(display, indent=2))
    raise SystemExit(1 if result["failures"] else 0)
