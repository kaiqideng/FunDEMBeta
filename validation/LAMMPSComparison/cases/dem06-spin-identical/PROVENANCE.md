# Historical provenance: dem06-spin-identical

Historical outcome: **PASS**. This is not a fresh acceptance result.

Official source: [LAMMPS reference](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem06-spin-identical-3d-si.yaml).

Original numerical-record SHA256: `2ac37a23a78fec3e0253066f2ba5d3f7ef8a760c095c63cdf3179ccad2a9c6d4`.
Original provenance-record SHA256: `5c96467fafb896c85f32d4f83bd01443d8b233b190248c66e0b8dc96d387dd8a`.
Original report SHA256 before link/image cleanup: `f64b54492e9c66cbb3a8a64d5cd44b1aa50a05853c8ab06f30f1fcdd591f4cff`.
Archived report SHA256 (HISTORICAL-REPORT.md): `5501e377f91a60398361cdce5e49c2c85ca433ddf749960700751a9bcf7b6ef9`.

Current report SHA256: `ddbca3aba87b757182e34a7609a9e9c70e6488ddc7b4591202ff79d3be5f4a0b`.

The table preserves the original provenance fields and numeric tokens. Historical format/hash statements refer to the archived original records; omitted artifact locations are not runtime inputs.

[Case summary](REPORT.md) · [Historical report](HISTORICAL-REPORT.md) · [Numerical record](HISTORICAL-RESULTS.md)

<details>
<summary>Preserved provenance fields</summary>

| Field | Type | Value |
| --- | --- | --- |
| `case` | string | `dem06-spin-identical` |
| `reference_family` | string | `LAMMPS` |
| `historical_outcome` | string | `PASS` |
| `historical_evidence_dates[0]` | string | `2026-10-28` |
| `historical_evidence_dates[1]` | string | `2026-10-29` |
| `migration_date` | string | `2026-10-30` |
| `status_policy` | string | `Historical comparison, not a new run or a current-core acceptance result.` |
| `reference_execution_policy` | string | `Use the case report: official saved LAMMPS states and analytical/public PFC targets; no local PFC execution is claimed.` |
| `official_demo_id` | string | `dem06-spin-identical` |
| `official_source` | string | `https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem06-spin-identical-3d-si.yaml` |
| `source_commit` | string | `29dd3fd7642a6fe9e8b4dd99120df4068b221250` |
| `official_git_blob` | string | `97f00c537ce8d677399d1fc390e6e5f70605c32f` |
| `source_snapshot_sha256` | string | `1ef6695681b3cf76e91737d9578839b4a767a23d7974b62e2c078d8b59f827a6` |
| `official_criterion` | string | `Official spin_no_friction selected-time analytic assertion at 0.00024 s; per-quantity error types/limits in analysis; source tags granular, unvetted` |
| `original_report_md_sha256` | string | `20d08667c0133d6009e3ccf62b85528beca6d45575e94f1de03d8bf902656018` |
| `original_analysis_sha256` | string | `f950cdce31f102c4f35f22406705236d29a1f9413294ef54e05b82116ac15830` |
| `published_analysis_sha256` | string | `2ac37a23a78fec3e0253066f2ba5d3f7ef8a760c095c63cdf3179ccad2a9c6d4` |
| `analysis_transformation` | string | `Numerical JSON retained; personal absolute artifact paths replaced with historical-archive labels. Full raw histories are not distributed.` |
| `compact_checkpoint_tables[0].file` | string | `checkpoints/official_metrics_comparison.csv` |
| `compact_checkpoint_tables[0].sha256` | string | `445d1d2074e43f75d3a99f9bf7d4a2f38d223a6781c2358ee1928ecf0aa57d13` |
| `compact_checkpoint_tables[0].bytes` | number | `1237` |
| `compact_checkpoint_tables[1].file` | string | `checkpoints/source_checkpoint_comparison.csv` |
| `compact_checkpoint_tables[1].sha256` | string | `4aba6593c6dda8bc92ea68bf26d595603af80d3a3ddf6ece758e2e541642d574` |
| `compact_checkpoint_tables[1].bytes` | number | `7321` |
| `historical_runs[0].run` | string | `cpu` |
| `historical_runs[0].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[0].steps` | number | `2400` |
| `historical_runs[0].stored_frames` | number | `2401` |
| `historical_runs[0].status` | string | `completed` |
| `historical_runs[0].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[0].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[0].run_json_sha256` | string | `6e2dd01dd454bcb45cbf3b27628f6fb01726cc5ca002da1b4435afc96783665c` |
| `historical_runs[0].project_snapshot_sha256` | string | `e0cbf7e364881ba15eb4ce6b00ab263d09ba73bdb2b8c4903d99fbe49b2308bd` |
| `historical_runs[0].native_output_interval_s` | number | `1e-07` |
| `historical_runs[0].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].trajectory_csv_sha256` | string | `d07b2213a081f42fd8229deffcbb95b3bb616c8a0c5e619c95619bd5970fdd0b` |
| `historical_runs[0].trajectory_csv_rows` | number | `4802` |
| `historical_runs[0].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].contacts_csv_sha256` | string | `f1f1051fbb2a86f21aae6190ebb540659c3205dd15bfdf3914c13b2342003565` |
| `historical_runs[0].contact_csv_rows` | number | `198` |
| `historical_runs[1].run` | string | `gpu` |
| `historical_runs[1].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[1].steps` | number | `2400` |
| `historical_runs[1].stored_frames` | number | `2401` |
| `historical_runs[1].status` | string | `completed` |
| `historical_runs[1].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[1].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[1].run_json_sha256` | string | `8b983bfdc78036ca23484b25c018e84f20881cd2409934db6376addc095ee7cd` |
| `historical_runs[1].project_snapshot_sha256` | string | `039377f170ccc9e71bc4ef241940d60a2fdf40475dae027b9b205cacba9a75d5` |
| `historical_runs[1].native_output_interval_s` | number | `1e-07` |
| `historical_runs[1].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].trajectory_csv_sha256` | string | `1cab5950c5df9f1aacd9f6e4901a3d6d937812942200c2777eb0730168ff43dc` |
| `historical_runs[1].trajectory_csv_rows` | number | `4802` |
| `historical_runs[1].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].contacts_csv_sha256` | string | `6a7f9fd89673981cba60b5179d807555acd22eddc1427d723f29e33c91385d97` |
| `historical_runs[1].contact_csv_rows` | number | `198` |
| `historical_runs[2].run` | string | `cpu-dt-half` |
| `historical_runs[2].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[2].steps` | number | `4800` |
| `historical_runs[2].stored_frames` | number | `4801` |
| `historical_runs[2].status` | string | `completed` |
| `historical_runs[2].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[2].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[2].run_json_sha256` | string | `26f03f62d43d92110a036e1e92d4ee544add5ccde40e2b7671e63e6b070e973f` |
| `historical_runs[2].project_snapshot_sha256` | string | `7ff2633ac353209206e4d6c67636423fdaea875c44510787b9d250771816d90d` |
| `historical_runs[2].native_output_interval_s` | number | `5e-08` |
| `historical_runs[2].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].trajectory_csv_sha256` | string | `32c26acbc35c24f4c81bea4b02ea7c0a1994da70b9055bae03a96ae4b04ff334` |
| `historical_runs[2].trajectory_csv_rows` | number | `9602` |
| `historical_runs[2].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].contacts_csv_sha256` | string | `f8d395a85b196056eefd9f49b344edb715b7861a7dc331cad35f5b1418686031` |
| `historical_runs[2].contact_csv_rows` | number | `396` |
| `historical_runs[3].run` | string | `gpu-dt-half` |
| `historical_runs[3].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[3].steps` | number | `4800` |
| `historical_runs[3].stored_frames` | number | `4801` |
| `historical_runs[3].status` | string | `completed` |
| `historical_runs[3].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[3].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[3].run_json_sha256` | string | `493d040d06287d78b6a0c346d600aa80778ab37338fa17bc9fba853d30374e47` |
| `historical_runs[3].project_snapshot_sha256` | string | `4fb8b5409f78c2af20124c0bff6c2be2460c63781863d34f09d7a342bfb82876` |
| `historical_runs[3].native_output_interval_s` | number | `5e-08` |
| `historical_runs[3].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].trajectory_csv_sha256` | string | `00a38232220855827a1593bc255b4a77c576d15908196f21dc6f58f4d2523158` |
| `historical_runs[3].trajectory_csv_rows` | number | `9602` |
| `historical_runs[3].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].contacts_csv_sha256` | string | `f1f75c6e9036d5458abe87d7f1cc7889d893c23ee3689bfa235849500c72af88` |
| `historical_runs[3].contact_csv_rows` | number | `396` |
| `published_report_md_sha256` | string | `f64b54492e9c66cbb3a8a64d5cd44b1aa50a05853c8ab06f30f1fcdd591f4cff` |
| `report_transformation` | string | `Historical prose and numerical results preserved; add historical-scope banner, update relative evidence links, append native-rerun and compact-evidence links.` |

</details>
