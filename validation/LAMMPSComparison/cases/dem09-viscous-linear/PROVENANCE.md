# Historical provenance: dem09-viscous-linear

Historical outcome: **PASS**. This is not a fresh acceptance result.

Official source: [LAMMPS reference](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem09-viscous-linear-3d-si.yaml).

Original numerical-record SHA256: `03b6ce64c6be669fd0cb682a706cf62707227d218fc2ed6c30d136f2660c876d`.
Original provenance-record SHA256: `3d8f01d55e0f147eaf6bc138067480124f5650e472929689f6a979e2ce8751f6`.
Original report SHA256 before link/image cleanup: `839ecdbcb7c211752da59cfbbbec06e91badbe1a92f7cb14c2bf9aacd1eeaf9a`.
Archived report SHA256 (HISTORICAL-REPORT.md): `d8ca118be619878363ff8a7565dd6543643d2ef49332489d0f32f1d93447b5cd`.

Current report SHA256: `66ac2698e2e4f386196de75c5a1c51015eb9ae9004cc11b0a734ef0bee4f20af`.

The table preserves the original provenance fields and numeric tokens. Historical format/hash statements refer to the archived original records; omitted artifact locations are not runtime inputs.

[Case summary](REPORT.md) · [Historical report](HISTORICAL-REPORT.md) · [Numerical record](HISTORICAL-RESULTS.md)

<details>
<summary>Preserved provenance fields</summary>

| Field | Type | Value |
| --- | --- | --- |
| `case` | string | `dem09-viscous-linear` |
| `reference_family` | string | `LAMMPS` |
| `historical_outcome` | string | `PASS` |
| `historical_evidence_dates[0]` | string | `2026-10-28` |
| `historical_evidence_dates[1]` | string | `2026-10-29` |
| `migration_date` | string | `2026-10-30` |
| `status_policy` | string | `Historical comparison, not a new run or a current-core acceptance result.` |
| `reference_execution_policy` | string | `Use the case report: official saved LAMMPS states and analytical/public PFC targets; no local PFC execution is claimed.` |
| `official_demo_id` | null | `null` |
| `official_source` | string | `https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem09-viscous-linear-3d-si.yaml` |
| `source_commit` | string | `29dd3fd7642a6fe9e8b4dd99120df4068b221250` |
| `official_git_blob` | string | `c1948eb56cac61f1c6159fc060d7d191c1ab897c` |
| `source_snapshot_sha256` | string | `49c4ead259f1f98b3a913e67665d862ff787041d1654ad0f0724ae20e3a5c419` |
| `official_criterion` | string | `0.5% relative terminal-speed assertion at 0.08 s: compare -vz with m*g/gamma (official tag unvetted)` |
| `original_report_md_sha256` | string | `ffaca3835b52475c58fa2931abf79bf16942e43abe18cfcc93000202619a0adb` |
| `original_analysis_sha256` | string | `399cbd584dad7bd79d60cddcdb0339ced3ea3f3a13c38175fedcdb189e4b1494` |
| `published_analysis_sha256` | string | `03b6ce64c6be669fd0cb682a706cf62707227d218fc2ed6c30d136f2660c876d` |
| `analysis_transformation` | string | `Numerical JSON retained; personal absolute artifact paths replaced with historical-archive labels. Full raw histories are not distributed.` |
| `compact_checkpoint_tables[0].file` | string | `checkpoints/official_checkpoint_comparison.csv` |
| `compact_checkpoint_tables[0].sha256` | string | `46c1e1364b922b4662705d88318e2f6bc86030c37076fa44f3bef7300e812d7d` |
| `compact_checkpoint_tables[0].bytes` | number | `1860` |
| `historical_runs[0].run` | string | `cpu` |
| `historical_runs[0].steps` | number | `800` |
| `historical_runs[0].stored_frames` | number | `801` |
| `historical_runs[0].status` | string | `completed` |
| `historical_runs[0].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[0].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[0].run_json_sha256` | string | `ebfbfca7e61eda929ba85847c9d55db41ccff9483035034448ece7825e2b255d` |
| `historical_runs[0].project_snapshot_sha256` | string | `5d197423fedea2fa16e506e7d8dc403e73891cf8c83b63a8eeaa179ef55b3ac7` |
| `historical_runs[0].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].trajectory_csv_sha256` | string | `1d2aef947ea25fba227a640aef22f385c2eb32ae98ce819ee4cd00207faf1873` |
| `historical_runs[0].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].contacts_csv_sha256` | string | `f390738cc9752ce6e39801ceadbbdca63598f6d92a0a88244b97100357dd9964` |
| `historical_runs[0].native_output_interval_s` | number | `0.0001` |
| `historical_runs[0].trajectory_csv_rows` | number | `801` |
| `historical_runs[0].contact_csv_rows` | number | `0` |
| `historical_runs[1].run` | string | `gpu` |
| `historical_runs[1].steps` | number | `800` |
| `historical_runs[1].stored_frames` | number | `801` |
| `historical_runs[1].status` | string | `completed` |
| `historical_runs[1].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[1].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[1].run_json_sha256` | string | `c8a826f1893efdeddffcc45e468e7d2e1e4ce8cb7634c0aa3538e546a61e0dcd` |
| `historical_runs[1].project_snapshot_sha256` | string | `e594384195c094cf8b002ca25e13dd037caa2d11467e49b14ee0e55dfa10bc1d` |
| `historical_runs[1].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].trajectory_csv_sha256` | string | `b7164f448fb2fd0a9671452729d82eb43f1c5b21962ef1fb0a1f5c00b2e9d1a1` |
| `historical_runs[1].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].contacts_csv_sha256` | string | `f390738cc9752ce6e39801ceadbbdca63598f6d92a0a88244b97100357dd9964` |
| `historical_runs[1].native_output_interval_s` | number | `0.0001` |
| `historical_runs[1].trajectory_csv_rows` | number | `801` |
| `historical_runs[1].contact_csv_rows` | number | `0` |
| `historical_runs[2].run` | string | `cpu-dt-half` |
| `historical_runs[2].steps` | number | `1600` |
| `historical_runs[2].stored_frames` | number | `1601` |
| `historical_runs[2].status` | string | `completed` |
| `historical_runs[2].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[2].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[2].run_json_sha256` | string | `96e45a0feaf7a35f3c0eb003fcf5e8d8ca46bf335f5879c006764e8928e30b1e` |
| `historical_runs[2].project_snapshot_sha256` | string | `df4f4131789bce6a2cf14603669e82e947e62d8ffbaf8f302ca3d8cc0d3bd22a` |
| `historical_runs[2].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].trajectory_csv_sha256` | string | `1ef4292647b51d1560b660d5fbf942cd99a7610b137fda26a8134ca587c8b3a9` |
| `historical_runs[2].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].contacts_csv_sha256` | string | `f390738cc9752ce6e39801ceadbbdca63598f6d92a0a88244b97100357dd9964` |
| `historical_runs[2].native_output_interval_s` | number | `5e-05` |
| `historical_runs[2].trajectory_csv_rows` | number | `1601` |
| `historical_runs[2].contact_csv_rows` | number | `0` |
| `historical_runs[3].run` | string | `gpu-dt-half` |
| `historical_runs[3].steps` | number | `1600` |
| `historical_runs[3].stored_frames` | number | `1601` |
| `historical_runs[3].status` | string | `completed` |
| `historical_runs[3].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[3].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[3].run_json_sha256` | string | `7e53c425877e076a167c0a6c2fc381c41c823841a33ce3f7b927c446d73f0392` |
| `historical_runs[3].project_snapshot_sha256` | string | `b85a90624b0a92af733d8ad6e850408d4447128c2130f53b330682784a9180a1` |
| `historical_runs[3].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].trajectory_csv_sha256` | string | `6056a3a88a484ab691aa2cea92ad2914c4c69ef6162c0d0da96929043a1ce74e` |
| `historical_runs[3].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].contacts_csv_sha256` | string | `f390738cc9752ce6e39801ceadbbdca63598f6d92a0a88244b97100357dd9964` |
| `historical_runs[3].native_output_interval_s` | number | `5e-05` |
| `historical_runs[3].trajectory_csv_rows` | number | `1601` |
| `historical_runs[3].contact_csv_rows` | number | `0` |
| `published_report_md_sha256` | string | `839ecdbcb7c211752da59cfbbbec06e91badbe1a92f7cb14c2bf9aacd1eeaf9a` |
| `report_transformation` | string | `Historical prose and numerical results preserved; add historical-scope banner, update relative evidence links, append native-rerun and compact-evidence links.` |

</details>
