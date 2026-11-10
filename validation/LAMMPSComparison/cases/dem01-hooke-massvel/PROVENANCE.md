# Historical provenance: dem01-hooke-massvel

Historical outcome: **PASS**. This is not a fresh acceptance result.

Official source: [LAMMPS reference](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem01-hooke-massvel-3d-si.yaml).

Original numerical-record SHA256: `d9a705a871e2ca9c91db17270888d29deab9c09485fba222a2dcb3bdb486ffa9`.
Original provenance-record SHA256: `980090e02968487e3ca61c5387f9ca8ec2983f748849136ef9e89b23ff6302ec`.
Original report SHA256 before link/image cleanup: `11ab313751de3d1f171e8105a47cb727d8d9ea91d66c552a26be13d9e6d9f1e1`.
Archived report SHA256 (HISTORICAL-REPORT.md): `d0a0b9ff736d3fa254aa105c4c6239707d8983f16698bb2aa0fbd2715cc8a9b9`.

Current report SHA256: `61f4333cc548408bbe6d2dc6671acf97b1c3b50d5148c2e2b8cdb1c80cbd37d8`.

The table preserves the original provenance fields and numeric tokens. Historical format/hash statements refer to the archived original records; omitted artifact locations are not runtime inputs.

[Case summary](REPORT.md) · [Historical report](HISTORICAL-REPORT.md) · [Numerical record](HISTORICAL-RESULTS.md)

<details>
<summary>Preserved provenance fields</summary>

| Field | Type | Value |
| --- | --- | --- |
| `case` | string | `dem01-hooke-massvel` |
| `reference_family` | string | `LAMMPS` |
| `historical_outcome` | string | `PASS` |
| `historical_evidence_dates[0]` | string | `2026-10-28` |
| `historical_evidence_dates[1]` | string | `2026-10-29` |
| `migration_date` | string | `2026-10-30` |
| `status_policy` | string | `Historical comparison, not a new run or a current-core acceptance result.` |
| `reference_execution_policy` | string | `Use the case report: official saved LAMMPS states and analytical/public PFC targets; no local PFC execution is claimed.` |
| `official_demo_id` | string | `dem01-hooke-massvel` |
| `official_source` | string | `https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem01-hooke-massvel-3d-si.yaml` |
| `source_commit` | string | `29dd3fd7642a6fe9e8b4dd99120df4068b221250` |
| `official_git_blob` | string | `9c634f2c0e3ee28173088739de4cfbba73bd76e4` |
| `source_snapshot_sha256` | string | `a14338f9fa8e9349f050f0494eafba8f6c6b350ab1a465465c5b1e3bb867f761` |
| `official_criterion` | string | `Official collision_restitution selected-time analytic assertion at 8.385e-05 s; per-quantity error types/limits in analysis; source tags granular, unvetted` |
| `original_report_md_sha256` | string | `06ee6be6833cb322cf0bc4e5294ee7abad2abb3e2b35a37f2db5669c5b0a0478` |
| `original_analysis_sha256` | string | `1ca41c0c267a603801152f36eb80903cf9c091506cadca22dfa026b53c32fe71` |
| `published_analysis_sha256` | string | `d9a705a871e2ca9c91db17270888d29deab9c09485fba222a2dcb3bdb486ffa9` |
| `analysis_transformation` | string | `Numerical JSON retained; personal absolute artifact paths replaced with historical-archive labels. Full raw histories are not distributed.` |
| `compact_checkpoint_tables[0].file` | string | `checkpoints/official_metrics_comparison.csv` |
| `compact_checkpoint_tables[0].sha256` | string | `5a072aaac7b03b37aae382451e1f8e7ccd8990bba058f94e9f14ef85dfbc90a0` |
| `compact_checkpoint_tables[0].bytes` | number | `1019` |
| `compact_checkpoint_tables[1].file` | string | `checkpoints/source_checkpoint_comparison.csv` |
| `compact_checkpoint_tables[1].sha256` | string | `6f90b325a1474bffd3d4956a97ab6ed782fff27267a0d62ae1f784ba5adecb8c` |
| `compact_checkpoint_tables[1].bytes` | number | `7005` |
| `historical_runs[0].run` | string | `cpu` |
| `historical_runs[0].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[0].steps` | number | `2795` |
| `historical_runs[0].stored_frames` | number | `2796` |
| `historical_runs[0].status` | string | `completed` |
| `historical_runs[0].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[0].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[0].run_json_sha256` | string | `1c691aaecf08527c60a2bd803f0ca637af43eaedb58100d4b9b5651546c64552` |
| `historical_runs[0].project_snapshot_sha256` | string | `d5f6d6b00730e476003248385b53e24dc22c73b3ead5e8fd089e5ee550ce29cf` |
| `historical_runs[0].native_output_interval_s` | number | `3e-08` |
| `historical_runs[0].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].trajectory_csv_sha256` | string | `3a8a9ef13a2d502768794233c51c5fd46b74938dec44d7e80a71d9ba2f9ca070` |
| `historical_runs[0].trajectory_csv_rows` | number | `5592` |
| `historical_runs[0].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].contacts_csv_sha256` | string | `b205dddd9a40816416e611f3077fa30be7dcf8ba06e0e4e7d4f163a125ae8d0c` |
| `historical_runs[0].contact_csv_rows` | number | `1007` |
| `historical_runs[1].run` | string | `gpu` |
| `historical_runs[1].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[1].steps` | number | `2795` |
| `historical_runs[1].stored_frames` | number | `2796` |
| `historical_runs[1].status` | string | `completed` |
| `historical_runs[1].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[1].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[1].run_json_sha256` | string | `c08fed40b923977d240891f483a1ca96213c0defadf081b4a0797eed2797d78b` |
| `historical_runs[1].project_snapshot_sha256` | string | `8e3973ccfc1b1947d3b4f3adad16693227744b1ecf31ced261c55c7016909c78` |
| `historical_runs[1].native_output_interval_s` | number | `3e-08` |
| `historical_runs[1].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].trajectory_csv_sha256` | string | `365cd541eb34e93c6470112bbecc690cec4f08046cd67f059ddcad7457084a36` |
| `historical_runs[1].trajectory_csv_rows` | number | `5592` |
| `historical_runs[1].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].contacts_csv_sha256` | string | `c175b8a1ba7fa8d2416f997cb018da337f2c67c403cd0383cd1eff1229e4feb5` |
| `historical_runs[1].contact_csv_rows` | number | `1007` |
| `historical_runs[2].run` | string | `cpu-dt-half` |
| `historical_runs[2].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[2].steps` | number | `5590` |
| `historical_runs[2].stored_frames` | number | `5591` |
| `historical_runs[2].status` | string | `completed` |
| `historical_runs[2].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[2].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[2].run_json_sha256` | string | `4306870dcee17a44138e4f7ce7b4ec9f67f40cecc3a4967e804668d709b74c6a` |
| `historical_runs[2].project_snapshot_sha256` | string | `0ae40736f909a4cea9dcfe3ad2c0c961821ab13bbff53fbfd5865bc8f95c36e3` |
| `historical_runs[2].native_output_interval_s` | number | `1.5e-08` |
| `historical_runs[2].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].trajectory_csv_sha256` | string | `b7fdd0be032332c418d3b6491aeeae9cd1cd2657e66c4cff03df17af511c5d67` |
| `historical_runs[2].trajectory_csv_rows` | number | `11182` |
| `historical_runs[2].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].contacts_csv_sha256` | string | `e5cf880f0565602037bfc86080d4f27193c703b8b843fe07a503a287881042fc` |
| `historical_runs[2].contact_csv_rows` | number | `2014` |
| `historical_runs[3].run` | string | `gpu-dt-half` |
| `historical_runs[3].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[3].steps` | number | `5590` |
| `historical_runs[3].stored_frames` | number | `5591` |
| `historical_runs[3].status` | string | `completed` |
| `historical_runs[3].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[3].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[3].run_json_sha256` | string | `c05d56c62ddaeb547dccc90cb906eeef1a71db162ebfec79b688e1475153a171` |
| `historical_runs[3].project_snapshot_sha256` | string | `82b79ac92752a59a8eeeafc636286516e564b4715a255efd76caed55b546029f` |
| `historical_runs[3].native_output_interval_s` | number | `1.5e-08` |
| `historical_runs[3].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].trajectory_csv_sha256` | string | `18aa8c0c3e54ee105954a39d17d9b8ec71032e05631c4425b5bbf5d0bbfbc5ea` |
| `historical_runs[3].trajectory_csv_rows` | number | `11182` |
| `historical_runs[3].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].contacts_csv_sha256` | string | `418c6efb60d0252ab137cc5c47fc3db66e926f5d4903ea7d58874f141780f417` |
| `historical_runs[3].contact_csv_rows` | number | `2014` |
| `published_report_md_sha256` | string | `11ab313751de3d1f171e8105a47cb727d8d9ea91d66c552a26be13d9e6d9f1e1` |
| `report_transformation` | string | `Historical prose and numerical results preserved; add historical-scope banner, update relative evidence links, append native-rerun and compact-evidence links.` |

</details>
