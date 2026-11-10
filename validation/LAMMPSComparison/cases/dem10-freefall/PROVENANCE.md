# Historical provenance: dem10-freefall

Historical outcome: **PASS**. This is not a fresh acceptance result.

Official source: [LAMMPS reference](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem10-freefall-3d-si.yaml).

Original numerical-record SHA256: `c16eacc9f11e6e9c048cd97606d35555eac2f3ff725d377e28d1e6f5342b9fe4`.
Original provenance-record SHA256: `1c36e99f11704465e373a9e18c86aad2eda55ef708b2aaacf056443d3ac0a78f`.
Original report SHA256 before link/image cleanup: `f7f97d9f2369df37ff32ee1f03cc7556302311b71ff1cab7f9f53e84176599bf`.
Archived report SHA256 (HISTORICAL-REPORT.md): `590eb1fc548b3831de9355891624d1cebeb651ffb5ddc140a2eef38a413251fa`.

Current report SHA256: `e0237c9a8e4f5246ea32450169aad10bf7639a4660746bb7a8c578ec0157b3cd`.

The table preserves the original provenance fields and numeric tokens. Historical format/hash statements refer to the archived original records; omitted artifact locations are not runtime inputs.

[Case summary](REPORT.md) · [Historical report](HISTORICAL-REPORT.md) · [Numerical record](HISTORICAL-RESULTS.md)

<details>
<summary>Preserved provenance fields</summary>

| Field | Type | Value |
| --- | --- | --- |
| `case` | string | `dem10-freefall` |
| `reference_family` | string | `LAMMPS` |
| `historical_outcome` | string | `PASS` |
| `historical_evidence_dates[0]` | string | `2026-10-28` |
| `historical_evidence_dates[1]` | string | `2026-10-29` |
| `migration_date` | string | `2026-10-30` |
| `status_policy` | string | `Historical comparison, not a new run or a current-core acceptance result.` |
| `reference_execution_policy` | string | `Use the case report: official saved LAMMPS states and analytical/public PFC targets; no local PFC execution is claimed.` |
| `official_demo_id` | string | `dem10-freefall` |
| `official_source` | string | `https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem10-freefall-3d-si.yaml` |
| `source_commit` | string | `29dd3fd7642a6fe9e8b4dd99120df4068b221250` |
| `official_git_blob` | string | `be3161683b34bb6cecb5bc8474e8a89bceebfa1b` |
| `source_snapshot_sha256` | string | `c8e1a18f434b724691bbbcd32d6636c264350b2318b83ed9953aa12a5aac8356` |
| `official_criterion` | string | `Official freefall selected-time analytic assertion at 0.08 s; per-quantity error types/limits in analysis; source tags granular, unvetted` |
| `original_report_md_sha256` | string | `3e32ba3db3dec0a32bc07c5983b3cb06adbf20c2cb76844871cf5d34232cee08` |
| `original_analysis_sha256` | string | `9ddb1ab588f5c7731a9ba8e29ce47193d31db27cc25c52c932029ddb7369ca83` |
| `published_analysis_sha256` | string | `c16eacc9f11e6e9c048cd97606d35555eac2f3ff725d377e28d1e6f5342b9fe4` |
| `analysis_transformation` | string | `Numerical JSON retained; personal absolute artifact paths replaced with historical-archive labels. Full raw histories are not distributed.` |
| `compact_checkpoint_tables[0].file` | string | `checkpoints/official_metrics_comparison.csv` |
| `compact_checkpoint_tables[0].sha256` | string | `4d0e7518c5652db9a592acbaf46f81ee33eaf448087ca8bf86ac583e4a6c26eb` |
| `compact_checkpoint_tables[0].bytes` | number | `1215` |
| `compact_checkpoint_tables[1].file` | string | `checkpoints/source_checkpoint_comparison.csv` |
| `compact_checkpoint_tables[1].sha256` | string | `fa2d1c8ac50ffec0af6cf276ecdd1f3f41f7c0222e73f598ed9eaed1f27120d2` |
| `compact_checkpoint_tables[1].bytes` | number | `3579` |
| `historical_runs[0].run` | string | `cpu` |
| `historical_runs[0].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[0].steps` | number | `8000` |
| `historical_runs[0].stored_frames` | number | `801` |
| `historical_runs[0].status` | string | `completed` |
| `historical_runs[0].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[0].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[0].run_json_sha256` | string | `1e8c8d4b28096ef4f39a7f963ede4d4c4d1d2993543f477645cd8bdd23bcb284` |
| `historical_runs[0].project_snapshot_sha256` | string | `bf05af8043974921b7f565fc2faa709decf7fdc522a7787bd3dcb13d1c1f2073` |
| `historical_runs[0].native_output_interval_s` | number | `0.0001` |
| `historical_runs[0].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].trajectory_csv_sha256` | string | `78c35f79bc992ef241041ecc4f5f7f71a02e50f4765acd4eb0392891db5b34a1` |
| `historical_runs[0].trajectory_csv_rows` | number | `8001` |
| `historical_runs[0].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].contacts_csv_sha256` | string | `4da48d470405d964a05f039b51c3db48d69a4f6e63f7e7ac1a1de6af953d74d1` |
| `historical_runs[0].contact_csv_rows` | number | `0` |
| `historical_runs[1].run` | string | `gpu` |
| `historical_runs[1].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[1].steps` | number | `8000` |
| `historical_runs[1].stored_frames` | number | `801` |
| `historical_runs[1].status` | string | `completed` |
| `historical_runs[1].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[1].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[1].run_json_sha256` | string | `8719362ba883d4ecf8b74abb91c305b4aafb3d9f19a5b23b9cbd4b56da0158a8` |
| `historical_runs[1].project_snapshot_sha256` | string | `986fc0d22c478751597c0e4c686b349d9260a5ad3caf89b12ed5ea40680d30ab` |
| `historical_runs[1].native_output_interval_s` | number | `0.0001` |
| `historical_runs[1].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].trajectory_csv_sha256` | string | `78c35f79bc992ef241041ecc4f5f7f71a02e50f4765acd4eb0392891db5b34a1` |
| `historical_runs[1].trajectory_csv_rows` | number | `8001` |
| `historical_runs[1].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].contacts_csv_sha256` | string | `4da48d470405d964a05f039b51c3db48d69a4f6e63f7e7ac1a1de6af953d74d1` |
| `historical_runs[1].contact_csv_rows` | number | `0` |
| `historical_runs[2].run` | string | `cpu-dt-half` |
| `historical_runs[2].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[2].steps` | number | `16000` |
| `historical_runs[2].stored_frames` | number | `801` |
| `historical_runs[2].status` | string | `completed` |
| `historical_runs[2].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[2].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[2].run_json_sha256` | string | `e5023cf595704c60a42822b0ac301ba109d084828fcf1f1cf65983edfc1cf0cc` |
| `historical_runs[2].project_snapshot_sha256` | string | `303d8ed7aee84a425be56aeb44fca7bb96d7e08171f1a23b037a82c6c1114a0f` |
| `historical_runs[2].native_output_interval_s` | number | `0.0001` |
| `historical_runs[2].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].trajectory_csv_sha256` | string | `d448a6a9d369b625e08c03482ff4c2b40461e541d9751b12b3d7962e3337ce52` |
| `historical_runs[2].trajectory_csv_rows` | number | `16001` |
| `historical_runs[2].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].contacts_csv_sha256` | string | `4da48d470405d964a05f039b51c3db48d69a4f6e63f7e7ac1a1de6af953d74d1` |
| `historical_runs[2].contact_csv_rows` | number | `0` |
| `historical_runs[3].run` | string | `gpu-dt-half` |
| `historical_runs[3].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[3].steps` | number | `16000` |
| `historical_runs[3].stored_frames` | number | `801` |
| `historical_runs[3].status` | string | `completed` |
| `historical_runs[3].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[3].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[3].run_json_sha256` | string | `23c83beebaf2b2cb21a017898788898d7c10f7009f7f0e3a760cfe8e7baabc55` |
| `historical_runs[3].project_snapshot_sha256` | string | `6fb1584d958e2790bcbfb737b1203d75505a3ee6ca2da8dde551079648342877` |
| `historical_runs[3].native_output_interval_s` | number | `0.0001` |
| `historical_runs[3].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].trajectory_csv_sha256` | string | `2776b31c7ff6524289164d18360fdd762ab1a0c44f81afea6d4f527c78405129` |
| `historical_runs[3].trajectory_csv_rows` | number | `16001` |
| `historical_runs[3].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].contacts_csv_sha256` | string | `4da48d470405d964a05f039b51c3db48d69a4f6e63f7e7ac1a1de6af953d74d1` |
| `historical_runs[3].contact_csv_rows` | number | `0` |
| `published_report_md_sha256` | string | `f7f97d9f2369df37ff32ee1f03cc7556302311b71ff1cab7f9f53e84176599bf` |
| `report_transformation` | string | `Historical prose and numerical results preserved; add historical-scope banner, update relative evidence links, append native-rerun and compact-evidence links.` |

</details>
