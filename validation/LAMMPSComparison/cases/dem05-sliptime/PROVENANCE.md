# Historical provenance: dem05-sliptime

Historical outcome: **PASS**. This is not a fresh acceptance result.

Official source: [LAMMPS reference](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem05-sliptime-3d-si.yaml).

Original numerical-record SHA256: `8329a1fda6b4d8fd4a707e123ec4f448d2f774a704034380c2d289a40b162666`.
Original provenance-record SHA256: `d5e7f4cc9e9a9a38731c7d14c7bd967a0ccc4456b4de7cbb4aa209ed06cb4e53`.
Original report SHA256 before link/image cleanup: `38358a0cafd74207bf4ac18b155f515cdf998a34b33dfe5d2e68058d1c46f115`.
Archived report SHA256 (HISTORICAL-REPORT.md): `c3b49e32e320bacb88bda6eff9e7161bee034d4710b5cc75aea6d2081fb0d3bd`.

Current report SHA256: `534b5f40490ed11d8287181988f8686cf5c5ff028494693719cef1e30144646b`.

The table preserves the original provenance fields and numeric tokens. Historical format/hash statements refer to the archived original records; omitted artifact locations are not runtime inputs.

[Case summary](REPORT.md) · [Historical report](HISTORICAL-REPORT.md) · [Numerical record](HISTORICAL-RESULTS.md)

<details>
<summary>Preserved provenance fields</summary>

| Field | Type | Value |
| --- | --- | --- |
| `case` | string | `dem05-sliptime` |
| `reference_family` | string | `LAMMPS` |
| `historical_outcome` | string | `PASS` |
| `historical_evidence_dates[0]` | string | `2026-10-28` |
| `historical_evidence_dates[1]` | string | `2026-10-29` |
| `migration_date` | string | `2026-10-30` |
| `status_policy` | string | `Historical comparison, not a new run or a current-core acceptance result.` |
| `reference_execution_policy` | string | `Use the case report: official saved LAMMPS states and analytical/public PFC targets; no local PFC execution is claimed.` |
| `official_demo_id` | string | `dem05-sliptime` |
| `official_source` | string | `https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem05-sliptime-3d-si.yaml` |
| `source_commit` | string | `29dd3fd7642a6fe9e8b4dd99120df4068b221250` |
| `official_git_blob` | string | `20f0a60a7ead5d669a420f48008aae1fd26fe928` |
| `source_snapshot_sha256` | string | `367bf8478d9b0c8d9d07ca41848d496fa63c8160412876038f3b66266bd5e97b` |
| `official_criterion` | string | `Official slip_transient selected-time analytic assertion at 0.4 s; per-quantity error types/limits in analysis; source tags granular, unvetted` |
| `original_report_md_sha256` | string | `0702122adcf0a20bea1f9da02426a6d778ed2e94f5d659232e2197b39db7fa0b` |
| `original_analysis_sha256` | string | `4114bc2b3821362a50a54bdf801964ec40b798d0208ed99e5df13bb434d38fc4` |
| `published_analysis_sha256` | string | `8329a1fda6b4d8fd4a707e123ec4f448d2f774a704034380c2d289a40b162666` |
| `analysis_transformation` | string | `Numerical JSON retained; personal absolute artifact paths replaced with historical-archive labels. Full raw histories are not distributed.` |
| `compact_checkpoint_tables[0].file` | string | `checkpoints/official_metrics_comparison.csv` |
| `compact_checkpoint_tables[0].sha256` | string | `c61071adafebd3d1e84dab064e9f5c9417fa98cdf2e71c7ac7b33b5ddef7de41` |
| `compact_checkpoint_tables[0].bytes` | number | `1219` |
| `compact_checkpoint_tables[1].file` | string | `checkpoints/source_checkpoint_comparison.csv` |
| `compact_checkpoint_tables[1].sha256` | string | `eba0bc51ff5c7a349380ea35139d9ccd65e9088bcaf59fad808ee0565da40a39` |
| `compact_checkpoint_tables[1].bytes` | number | `10775` |
| `historical_runs[0].run` | string | `cpu` |
| `historical_runs[0].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[0].steps` | number | `100000` |
| `historical_runs[0].stored_frames` | number | `1001` |
| `historical_runs[0].status` | string | `completed` |
| `historical_runs[0].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[0].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[0].run_json_sha256` | string | `ae582b25c5ef29fd7cd76b7cc64dc80297a0ee671fb96e54c8f805fc11d6893f` |
| `historical_runs[0].project_snapshot_sha256` | string | `e73c62baa05523f22337bf24d8e00f5d322972e6c3a46f2fc0147d209049f6d5` |
| `historical_runs[0].native_output_interval_s` | number | `0.001` |
| `historical_runs[0].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].trajectory_csv_sha256` | string | `1894265e29f2db3ac1f4d09c2394806d905d74b80c2b75486241e2da6bf043e7` |
| `historical_runs[0].trajectory_csv_rows` | number | `100001` |
| `historical_runs[0].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].contacts_csv_sha256` | string | `7fbfd0f95399e845267244c15f658f26e1123439ba027c7f0f5c25a83a460d20` |
| `historical_runs[0].contact_csv_rows` | number | `100000` |
| `historical_runs[1].run` | string | `gpu` |
| `historical_runs[1].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[1].steps` | number | `100000` |
| `historical_runs[1].stored_frames` | number | `1001` |
| `historical_runs[1].status` | string | `completed` |
| `historical_runs[1].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[1].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[1].run_json_sha256` | string | `3cf273f041f75cbcea85273fd736560fe34d1c6fd4a349bf3964988b71af3574` |
| `historical_runs[1].project_snapshot_sha256` | string | `95758ffbb9bceade557c96a14625f0614bae3f390143665828f53abefe0943ce` |
| `historical_runs[1].native_output_interval_s` | number | `0.001` |
| `historical_runs[1].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].trajectory_csv_sha256` | string | `53421136e5bc108cc1a11bf171cf17d5980b6ce553f1ce69e023d0a85acead0b` |
| `historical_runs[1].trajectory_csv_rows` | number | `100001` |
| `historical_runs[1].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].contacts_csv_sha256` | string | `0fbd108d8cf2904e0ec3d63076d23924601ed89991319483067bd82bf5e33e04` |
| `historical_runs[1].contact_csv_rows` | number | `100000` |
| `historical_runs[2].run` | string | `cpu-dt-half` |
| `historical_runs[2].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[2].steps` | number | `200000` |
| `historical_runs[2].stored_frames` | number | `1001` |
| `historical_runs[2].status` | string | `completed` |
| `historical_runs[2].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[2].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[2].run_json_sha256` | string | `d10bfee4f010c81bd542bcf80d81ad6d67c0453d7539ec6acba523725bc601cb` |
| `historical_runs[2].project_snapshot_sha256` | string | `83dbee403fcb33ba9aab69af36c162ce4e54ec0f2d5107fcef72cb70709bd430` |
| `historical_runs[2].native_output_interval_s` | number | `0.001` |
| `historical_runs[2].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].trajectory_csv_sha256` | string | `7730ecde2e11d139a4f1749b558c06ad7e622c5564f67661183ef577e88bc9b1` |
| `historical_runs[2].trajectory_csv_rows` | number | `200001` |
| `historical_runs[2].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].contacts_csv_sha256` | string | `9693df121a6b2ccf771fbef760108c0a941b1b08200401c11461753a0520983c` |
| `historical_runs[2].contact_csv_rows` | number | `200000` |
| `historical_runs[3].run` | string | `gpu-dt-half` |
| `historical_runs[3].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[3].steps` | number | `200000` |
| `historical_runs[3].stored_frames` | number | `1001` |
| `historical_runs[3].status` | string | `completed` |
| `historical_runs[3].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[3].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[3].run_json_sha256` | string | `018b32febed22e6cc13eecc8be0ad1c068e2f61efb1198dccf305e729db9e927` |
| `historical_runs[3].project_snapshot_sha256` | string | `7bcc2dbc5b5027bdae0ca0d25ebcbb3a2041e9def5391a84bd96c8acca1f1391` |
| `historical_runs[3].native_output_interval_s` | number | `0.001` |
| `historical_runs[3].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].trajectory_csv_sha256` | string | `bd5f78319e767d5a40b21376cd7e42086f56bcff1db610ab7b8060575d327359` |
| `historical_runs[3].trajectory_csv_rows` | number | `200001` |
| `historical_runs[3].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].contacts_csv_sha256` | string | `5259ee45edcfd43924b7093758c8c10e3692c62684c87cc5bcf7eebc75c351aa` |
| `historical_runs[3].contact_csv_rows` | number | `200000` |
| `published_report_md_sha256` | string | `38358a0cafd74207bf4ac18b155f515cdf998a34b33dfe5d2e68058d1c46f115` |
| `report_transformation` | string | `Historical prose and numerical results preserved; add historical-scope banner, update relative evidence links, append native-rerun and compact-evidence links.` |

</details>
