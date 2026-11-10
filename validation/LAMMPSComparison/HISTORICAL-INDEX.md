# LAMMPS historical comparisons

These are historical outcomes, not a new run of the current core.

[Family guide](README.md) explains the report structure and reproduction flow.
Each concise `REPORT.md` links its byte-preserved `HISTORICAL-REPORT.md`;
the numerical CSV records and original acceptance limits are unchanged.

| Field | Type | Value |
| --- | --- | --- |
| `reference_family` | string | `LAMMPS` |
| `historical_evidence_dates[0]` | string | `2026-10-28` |
| `historical_evidence_dates[1]` | string | `2026-10-29` |
| `migration_date` | string | `2026-10-30` |
| `case_count` | number | `35` |
| `historical_outcomes.PASS` | number | `33` |
| `historical_outcomes.FAIL` | number | `2` |

| Case | Historical outcome | C++ source | Report |
| --- | --- | --- | --- |
| `dem01-hertz-tsuji` | **PASS** | [Source](cases/dem01-hertz-tsuji/dem01-hertz-tsuji.cpp) | [Report](cases/dem01-hertz-tsuji/REPORT.md) |
| `dem01-hertz-twosphere` | **PASS** | [Source](cases/dem01-hertz-twosphere/dem01-hertz-twosphere.cpp) | [Report](cases/dem01-hertz-twosphere/REPORT.md) |
| `dem01-hooke-2d-lj` | **PASS** | [Source](cases/dem01-hooke-2d-lj/dem01-hooke-2d-lj.cpp) | [Report](cases/dem01-hooke-2d-lj/REPORT.md) |
| `dem01-hooke-massvel` | **PASS** | [Source](cases/dem01-hooke-massvel/dem01-hooke-massvel.cpp) | [Report](cases/dem01-hooke-massvel/REPORT.md) |
| `dem01-hooke-twosphere` | **PASS** | [Source](cases/dem01-hooke-twosphere/dem01-hooke-twosphere.cpp) | [Report](cases/dem01-hooke-twosphere/REPORT.md) |
| `dem01-hooke-twosphere-e02` | **PASS** | [Source](cases/dem01-hooke-twosphere-e02/dem01-hooke-twosphere-e02.cpp) | [Report](cases/dem01-hooke-twosphere-e02/REPORT.md) |
| `dem01-hooke-twosphere-e08` | **PASS** | [Source](cases/dem01-hooke-twosphere-e08/dem01-hooke-twosphere-e08.cpp) | [Report](cases/dem01-hooke-twosphere-e08/REPORT.md) |
| `dem01-hooke-veldamp` | **PASS** | [Source](cases/dem01-hooke-veldamp/dem01-hooke-veldamp.cpp) | [Report](cases/dem01-hooke-veldamp/REPORT.md) |
| `dem01-hookehistory-twosphere` | **PASS** | [Source](cases/dem01-hookehistory-twosphere/dem01-hookehistory-twosphere.cpp) | [Report](cases/dem01-hookehistory-twosphere/REPORT.md) |
| `dem01-legacy-hooke-twosphere` | **PASS** | [Source](cases/dem01-legacy-hooke-twosphere/dem01-legacy-hooke-twosphere.cpp) | [Report](cases/dem01-legacy-hooke-twosphere/REPORT.md) |
| `dem01-legacy-hookehistory-twosphere` | **PASS** | [Source](cases/dem01-legacy-hookehistory-twosphere/dem01-legacy-hookehistory-twosphere.cpp) | [Report](cases/dem01-legacy-hookehistory-twosphere/REPORT.md) |
| `dem02-hertz-peak` | **PASS** | [Source](cases/dem02-hertz-peak/dem02-hertz-peak.cpp) | [Report](cases/dem02-hertz-peak/REPORT.md) |
| `dem02-hertz-twosphere` | **PASS** | [Source](cases/dem02-hertz-twosphere/dem02-hertz-twosphere.cpp) | [Report](cases/dem02-hertz-twosphere/REPORT.md) |
| `dem02-hertzmaterial-twosphere` | **PASS** | [Source](cases/dem02-hertzmaterial-twosphere/dem02-hertzmaterial-twosphere.cpp) | [Report](cases/dem02-hertzmaterial-twosphere/REPORT.md) |
| `dem02-legacy-hertzhistory-twosphere` | **PASS** | [Source](cases/dem02-legacy-hertzhistory-twosphere/dem02-legacy-hertzhistory-twosphere.cpp) | [Report](cases/dem02-legacy-hertzhistory-twosphere/REPORT.md) |
| `dem03-hertz-peak-wall` | **PASS** | [Source](cases/dem03-hertz-peak-wall/dem03-hertz-peak-wall.cpp) | [Report](cases/dem03-hertz-peak-wall/REPORT.md) |
| `dem03-hertz-tsuji-wall` | **PASS** | [Source](cases/dem03-hertz-tsuji-wall/dem03-hertz-tsuji-wall.cpp) | [Report](cases/dem03-hertz-tsuji-wall/REPORT.md) |
| `dem03-hertz-wall` | **PASS** | [Source](cases/dem03-hertz-wall/dem03-hertz-wall.cpp) | [Report](cases/dem03-hertz-wall/REPORT.md) |
| `dem03-hertzmaterial-wall` | **PASS** | [Source](cases/dem03-hertzmaterial-wall/dem03-hertzmaterial-wall.cpp) | [Report](cases/dem03-hertzmaterial-wall/REPORT.md) |
| `dem03-legacy-hertzhistory-wall` | **PASS** | [Source](cases/dem03-legacy-hertzhistory-wall/dem03-legacy-hertzhistory-wall.cpp) | [Report](cases/dem03-legacy-hertzhistory-wall/REPORT.md) |
| `dem04-hertz-oblique` | **PASS** | [Source](cases/dem04-hertz-oblique/dem04-hertz-oblique.cpp) | [Report](cases/dem04-hertz-oblique/REPORT.md) |
| `dem04-hookehistory-oblique` | **PASS** | [Source](cases/dem04-hookehistory-oblique/dem04-hookehistory-oblique.cpp) | [Report](cases/dem04-hookehistory-oblique/REPORT.md) |
| `dem05-hookehistory-slip` | **PASS** | [Source](cases/dem05-hookehistory-slip/dem05-hookehistory-slip.cpp) | [Report](cases/dem05-hookehistory-slip/REPORT.md) |
| `dem05-mindlin-slip` | **FAIL** | [Source](cases/dem05-mindlin-slip/dem05-mindlin-slip.cpp) | [Report](cases/dem05-mindlin-slip/REPORT.md) |
| `dem05-slip-mu02` | **FAIL** | [Source](cases/dem05-slip-mu02/dem05-slip-mu02.cpp) | [Report](cases/dem05-slip-mu02/REPORT.md) |
| `dem05-sliptime` | **PASS** | [Source](cases/dem05-sliptime/dem05-sliptime.cpp) | [Report](cases/dem05-sliptime/REPORT.md) |
| `dem06-hookehistory-spin` | **PASS** | [Source](cases/dem06-hookehistory-spin/dem06-hookehistory-spin.cpp) | [Report](cases/dem06-hookehistory-spin/REPORT.md) |
| `dem06-spin-identical` | **PASS** | [Source](cases/dem06-spin-identical/dem06-spin-identical.cpp) | [Report](cases/dem06-spin-identical/REPORT.md) |
| `dem07-hooke-rolling` | **PASS** | [Source](cases/dem07-hooke-rolling/dem07-hooke-rolling.cpp) | [Report](cases/dem07-hooke-rolling/REPORT.md) |
| `dem09-viscous-linear` | **PASS** | [Source](cases/dem09-viscous-linear/dem09-viscous-linear.cpp) | [Report](cases/dem09-viscous-linear/REPORT.md) |
| `dem10-freefall` | **PASS** | [Source](cases/dem10-freefall/dem10-freefall.cpp) | [Report](cases/dem10-freefall/REPORT.md) |
| `dem10-stack` | **PASS** | [Source](cases/dem10-stack/dem10-stack.cpp) | [Report](cases/dem10-stack/REPORT.md) |
| `dem11-regionwall` | **PASS** | [Source](cases/dem11-regionwall/dem11-regionwall.cpp) | [Report](cases/dem11-regionwall/REPORT.md) |
| `dem13-twist-sds` | **PASS** | [Source](cases/dem13-twist-sds/dem13-twist-sds.cpp) | [Report](cases/dem13-twist-sds/REPORT.md) |
| `dem15-hooke-oblique` | **PASS** | [Source](cases/dem15-hooke-oblique/dem15-hooke-oblique.cpp) | [Report](cases/dem15-hooke-oblique/REPORT.md) |
