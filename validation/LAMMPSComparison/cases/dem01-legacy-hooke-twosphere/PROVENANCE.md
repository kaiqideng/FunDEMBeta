# Historical provenance: dem01-legacy-hooke-twosphere

Historical outcome: **PASS**. This is not a fresh acceptance result.

Official source: [LAMMPS reference](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem01-legacy-hooke-twosphere-3d-si.yaml).

Original numerical-record SHA256: `d34e01ab8e11ce3b39fd7af791b90d33cd624a393cc5a91b58dc339cb050da9d`.
Original provenance-record SHA256: `5dc9ddba2131a0f16bd15f409f764a966ea7d9dad84f919efa44f1205f13e126`.
Original report SHA256 before link/image cleanup: `c7a1677cbca4a6bc0f46c2ffc378d96dd44e0e649e85a1c597523f1076c27ca9`.
Archived report SHA256 (HISTORICAL-REPORT.md): `b2f94f6fc45726887921baafa464e15235af7e0caacc5696cc9f2e7794ea74c2`.

Current report SHA256: `d215db568a767361c107d5b0f4e8e3a957ad98856810d308b8bf7e50ac6fe7f6`.

The table preserves the original provenance fields and numeric tokens. Historical format/hash statements refer to the archived original records; omitted artifact locations are not runtime inputs.

[Case summary](REPORT.md) · [Historical report](HISTORICAL-REPORT.md) · [Numerical record](HISTORICAL-RESULTS.md)

<details>
<summary>Preserved provenance fields</summary>

| Field | Type | Value |
| --- | --- | --- |
| `case` | string | `dem01-legacy-hooke-twosphere` |
| `reference_family` | string | `LAMMPS` |
| `historical_outcome` | string | `PASS` |
| `historical_evidence_dates[0]` | string | `2026-10-28` |
| `historical_evidence_dates[1]` | string | `2026-10-29` |
| `migration_date` | string | `2026-10-30` |
| `status_policy` | string | `Historical comparison, not a new run or a current-core acceptance result.` |
| `reference_execution_policy` | string | `Use the case report: official saved LAMMPS states and analytical/public PFC targets; no local PFC execution is claimed.` |
| `official_demo_id` | string | `dem01-legacy-hooke-twosphere` |
| `official_source` | string | `https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem01-legacy-hooke-twosphere-3d-si.yaml` |
| `source_commit` | string | `29dd3fd7642a6fe9e8b4dd99120df4068b221250` |
| `official_git_blob` | string | `fcb49be715ccc2b5efd6ed6c005c33d47256842c` |
| `source_snapshot_sha256` | string | `83615aaf6067ed57ae4d8a76e27ba6b6958537513e31173155fb3e78f1360188` |
| `official_criterion` | string | `Official collision_restitution selected-time analytic assertion at 8.385e-05 s; per-quantity error types/limits in analysis; source tags granular, legacy, unvetted` |
| `original_report_md_sha256` | string | `e3ed6b5504038a10d7c2707976270f88953a02ad72a96d1ddbae4db2623b81b0` |
| `original_analysis_sha256` | string | `f8c46e66677bc1a3c3bec9de8c53e47e8a7b280c5e2d92df9784106f1735236e` |
| `published_analysis_sha256` | string | `d34e01ab8e11ce3b39fd7af791b90d33cd624a393cc5a91b58dc339cb050da9d` |
| `analysis_transformation` | string | `Numerical JSON retained; personal absolute artifact paths replaced with historical-archive labels. Full raw histories are not distributed.` |
| `compact_checkpoint_tables[0].file` | string | `checkpoints/official_metrics_comparison.csv` |
| `compact_checkpoint_tables[0].sha256` | string | `9f879bf55b13a7ec2ee5bc9b173366447748c0fad4a743af64bf1c18761696df` |
| `compact_checkpoint_tables[0].bytes` | number | `1017` |
| `compact_checkpoint_tables[1].file` | string | `checkpoints/source_checkpoint_comparison.csv` |
| `compact_checkpoint_tables[1].sha256` | string | `5fa230de3b1470c568e1a39ca51b07aa0662fab9ad55f499d4137c16361d893a` |
| `compact_checkpoint_tables[1].bytes` | number | `8934` |
| `historical_runs[0].run` | string | `cpu` |
| `historical_runs[0].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[0].steps` | number | `2795` |
| `historical_runs[0].stored_frames` | number | `2796` |
| `historical_runs[0].status` | string | `completed` |
| `historical_runs[0].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[0].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[0].run_json_sha256` | string | `0407bb12618f665754f0c1c9de5618576c99cb6279c7a33273f0c8b7ddddb5f3` |
| `historical_runs[0].project_snapshot_sha256` | string | `b0483a8e4c883b987c76e31c07a0e39355485c29210f9427874abb1507086fdc` |
| `historical_runs[0].native_output_interval_s` | number | `3e-08` |
| `historical_runs[0].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].trajectory_csv_sha256` | string | `3a360e4d67de66e8395146f19fa0038cb6c7cb33bed30537d8602178c5722320` |
| `historical_runs[0].trajectory_csv_rows` | number | `5592` |
| `historical_runs[0].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].contacts_csv_sha256` | string | `0cb5c0350f17493d87b881b3e73e9cad9d8aa992f14f57899375d743e00f4b5a` |
| `historical_runs[0].contact_csv_rows` | number | `1008` |
| `historical_runs[1].run` | string | `gpu` |
| `historical_runs[1].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[1].steps` | number | `2795` |
| `historical_runs[1].stored_frames` | number | `2796` |
| `historical_runs[1].status` | string | `completed` |
| `historical_runs[1].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[1].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[1].run_json_sha256` | string | `543fc6297ff71854cad0ce22472a3ee56e2ac9ce89712d00bba064b3c1aec255` |
| `historical_runs[1].project_snapshot_sha256` | string | `ea835f6c92e9991c62ebf465b899a6f07eae13cb9ce5a646da191badad73f9d9` |
| `historical_runs[1].native_output_interval_s` | number | `3e-08` |
| `historical_runs[1].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].trajectory_csv_sha256` | string | `c07cd47706aada35876ad5a0ef3b6c60d955ceb40e3762d1fb9a97da2588f781` |
| `historical_runs[1].trajectory_csv_rows` | number | `5592` |
| `historical_runs[1].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].contacts_csv_sha256` | string | `02f5a3736e1afe3fd72ccc26e3b0cfaf0df6594163b5000e69b92c9f04ba0ae9` |
| `historical_runs[1].contact_csv_rows` | number | `1008` |
| `historical_runs[2].run` | string | `cpu-dt-half` |
| `historical_runs[2].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[2].steps` | number | `5590` |
| `historical_runs[2].stored_frames` | number | `5591` |
| `historical_runs[2].status` | string | `completed` |
| `historical_runs[2].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[2].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[2].run_json_sha256` | string | `3f530b979eaf8b1bc71506c335fd03ea1ea311c8c61e084228c61159c5711e3c` |
| `historical_runs[2].project_snapshot_sha256` | string | `70d0995331ffbb3cf9f493ce8d146c10a9fa4bb7c92379b652d5a5dd3a22590b` |
| `historical_runs[2].native_output_interval_s` | number | `1.5e-08` |
| `historical_runs[2].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].trajectory_csv_sha256` | string | `b57c4729362bf8247a8a30a493e8270ebdcadd6bba65991ae49a91589316a3c7` |
| `historical_runs[2].trajectory_csv_rows` | number | `11182` |
| `historical_runs[2].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].contacts_csv_sha256` | string | `d5a8c77f3c58c1d1470b9ba16941081530f7708006a8f48ca382df29ce74c616` |
| `historical_runs[2].contact_csv_rows` | number | `2016` |
| `historical_runs[3].run` | string | `gpu-dt-half` |
| `historical_runs[3].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[3].steps` | number | `5590` |
| `historical_runs[3].stored_frames` | number | `5591` |
| `historical_runs[3].status` | string | `completed` |
| `historical_runs[3].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[3].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[3].run_json_sha256` | string | `4829d3a8cac5e49ffa4c74993cbbf2395ca769adf91ea1ed6e19ad96306fcf2a` |
| `historical_runs[3].project_snapshot_sha256` | string | `d3b168afe46ac9e4b87d4968e0d1c960c2a81bfc478396a0f52745cb562238b8` |
| `historical_runs[3].native_output_interval_s` | number | `1.5e-08` |
| `historical_runs[3].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].trajectory_csv_sha256` | string | `13c31db2ff340fed3d4649fcf4163b8d3fa4c26bf067c0a64facb323a7ed3e89` |
| `historical_runs[3].trajectory_csv_rows` | number | `11182` |
| `historical_runs[3].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].contacts_csv_sha256` | string | `a8e37ecba9880e97bf9729e8e6bff7d456444df0409e7131c5589e6608967a5b` |
| `historical_runs[3].contact_csv_rows` | number | `2016` |
| `published_report_md_sha256` | string | `c7a1677cbca4a6bc0f46c2ffc378d96dd44e0e649e85a1c597523f1076c27ca9` |
| `report_transformation` | string | `Historical prose and numerical results preserved; add historical-scope banner, update relative evidence links, append native-rerun and compact-evidence links.` |

</details>
