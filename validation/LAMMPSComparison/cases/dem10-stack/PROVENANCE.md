# Historical provenance: dem10-stack

Historical outcome: **PASS**. This is not a fresh acceptance result.

Official source: [LAMMPS reference](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem10-stack-3d-si.yaml).

Original numerical-record SHA256: `004cb89db352239b579a2466b5cea5bf647c39ad7ee273f5edb60808be1c998a`.
Original provenance-record SHA256: `d72b8cabe459297eb9adac379b0020b6929bc7cd4192e8de53dec3bde1d260f2`.
Original report SHA256 before link/image cleanup: `85f2a17cf104e22d02c8ffb949b2e2df7c46eb664376c22342e94eb707efc903`.
Archived report SHA256 (HISTORICAL-REPORT.md): `53830a799f0f87f6e3a6529b6417ff4bdb50f06aac5fbac2e6bcca511272a234`.

Current report SHA256: `50240eb6a3175d3c3bbbbfa8643f808886c2b946d2cc9316077dd29fdea3d1a2`.

The table preserves the original provenance fields and numeric tokens. Historical format/hash statements refer to the archived original records; omitted artifact locations are not runtime inputs.

[Case summary](REPORT.md) · [Historical report](HISTORICAL-REPORT.md) · [Numerical record](HISTORICAL-RESULTS.md)

<details>
<summary>Preserved provenance fields</summary>

| Field | Type | Value |
| --- | --- | --- |
| `case` | string | `dem10-stack` |
| `reference_family` | string | `LAMMPS` |
| `historical_outcome` | string | `PASS` |
| `historical_evidence_dates[0]` | string | `2026-10-28` |
| `historical_evidence_dates[1]` | string | `2026-10-29` |
| `migration_date` | string | `2026-10-30` |
| `status_policy` | string | `Historical comparison, not a new run or a current-core acceptance result.` |
| `reference_execution_policy` | string | `Use the case report: official saved LAMMPS states and analytical/public PFC targets; no local PFC execution is claimed.` |
| `official_demo_id` | string | `dem10-stack` |
| `official_source` | string | `https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem10-stack-3d-si.yaml` |
| `source_commit` | string | `29dd3fd7642a6fe9e8b4dd99120df4068b221250` |
| `official_git_blob` | string | `405583bd3e27294ab46506f6ba242485a756de33` |
| `source_snapshot_sha256` | string | `78c0e3ef4c79448cf1f6232369f30458d4a2e16c8136e989710ac493a5ae0a26` |
| `official_criterion` | string | `0.05% relative total mechanical-energy error at 20 ms (official tag unvetted)` |
| `original_report_md_sha256` | string | `ae454a6d2173f99ab1fbd8282cc909ea41a87694a7d98f8ee87836fdc4b8ac17` |
| `original_analysis_sha256` | string | `35103b0d003f770d7d57c3a2f5355188bcaef55570d7ab459732b7ed12db0515` |
| `published_analysis_sha256` | string | `004cb89db352239b579a2466b5cea5bf647c39ad7ee273f5edb60808be1c998a` |
| `analysis_transformation` | string | `Numerical JSON retained; personal absolute artifact paths replaced with historical-archive labels. Full raw histories are not distributed.` |
| `compact_checkpoint_tables[0].file` | string | `checkpoints/official_checkpoint_comparison.csv` |
| `compact_checkpoint_tables[0].sha256` | string | `ec6204d2dd159bda72a901936347417a6578474d5165b087b9c7994b048db45c` |
| `compact_checkpoint_tables[0].bytes` | number | `7286` |
| `historical_runs[0].run` | string | `cpu` |
| `historical_runs[0].steps` | number | `10000` |
| `historical_runs[0].stored_frames` | number | `10001` |
| `historical_runs[0].status` | string | `completed` |
| `historical_runs[0].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[0].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[0].run_json_sha256` | string | `1412eb3a746a27408f2764a2caed72845c41308a51860b9ba810c6728912eab0` |
| `historical_runs[0].project_snapshot_sha256` | string | `20bb56832c6d23df080e78b953ddf8f50248c2eb2700b1f475ff6dc83c722a03` |
| `historical_runs[0].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].trajectory_csv_sha256` | string | `5a870138201fe9e6888b3664f9e16fcf443e486f4c7667d960109925b324d209` |
| `historical_runs[0].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].contacts_csv_sha256` | string | `8b8163d9877a993bf649016b8bec56825bc7f8b97d40b1803baf15114c287e5b` |
| `historical_runs[0].native_output_interval_s` | number | `2e-06` |
| `historical_runs[0].trajectory_csv_rows` | number | `20002` |
| `historical_runs[0].contact_csv_rows` | number | `30000` |
| `historical_runs[1].run` | string | `gpu` |
| `historical_runs[1].steps` | number | `10000` |
| `historical_runs[1].stored_frames` | number | `10001` |
| `historical_runs[1].status` | string | `completed` |
| `historical_runs[1].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[1].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[1].run_json_sha256` | string | `4da18945e09405c3ffcdd39fcd06874c0248e21934b2733eba44eac051f9263b` |
| `historical_runs[1].project_snapshot_sha256` | string | `22aca48bc8ccece3ef2c35e57a4b56a6878ea3b6cf8546373b0b0a07a1c6b68d` |
| `historical_runs[1].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].trajectory_csv_sha256` | string | `4253c42caba1572630d4c6a0092a0f01d8217b0a3fd41511234616676aa5e086` |
| `historical_runs[1].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].contacts_csv_sha256` | string | `807e960915c270cb295ce482c7fbd81eade4ac0380968c932f059c38ab589433` |
| `historical_runs[1].native_output_interval_s` | number | `2e-06` |
| `historical_runs[1].trajectory_csv_rows` | number | `20002` |
| `historical_runs[1].contact_csv_rows` | number | `30000` |
| `historical_runs[2].run` | string | `cpu-dt-half` |
| `historical_runs[2].steps` | number | `20000` |
| `historical_runs[2].stored_frames` | number | `20001` |
| `historical_runs[2].status` | string | `completed` |
| `historical_runs[2].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[2].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[2].run_json_sha256` | string | `1798072a7b33a7f68e71bcddc899fb4f51e74081d6449659539ecc0921df41f5` |
| `historical_runs[2].project_snapshot_sha256` | string | `d1e5c972696e41ba2b62e68c6193db07b2320619095d82d27a7fc14e00dff51b` |
| `historical_runs[2].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].trajectory_csv_sha256` | string | `6adfce9d701a539e75341184342eba03eca6c879b18131b6a4c5a9be82b39341` |
| `historical_runs[2].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].contacts_csv_sha256` | string | `af1b43e9946fd127babc52d869bd5e9920963a7e642540a6372b1b7786f012a6` |
| `historical_runs[2].native_output_interval_s` | number | `1e-06` |
| `historical_runs[2].trajectory_csv_rows` | number | `40002` |
| `historical_runs[2].contact_csv_rows` | number | `60000` |
| `historical_runs[3].run` | string | `gpu-dt-half` |
| `historical_runs[3].steps` | number | `20000` |
| `historical_runs[3].stored_frames` | number | `20001` |
| `historical_runs[3].status` | string | `completed` |
| `historical_runs[3].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[3].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[3].run_json_sha256` | string | `6b0892c373dd2ba58d27684874810579ef6e3b7d3dc355a2ba5241da93762e24` |
| `historical_runs[3].project_snapshot_sha256` | string | `b5f2fb7a58ec69ce527622e1bc6d583dd5151416f2a8316cba4efcea23f95708` |
| `historical_runs[3].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].trajectory_csv_sha256` | string | `3d6e3026eddbd295cd2851e2645d900241cf91017f1e12b2d8a81efb3076bca5` |
| `historical_runs[3].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].contacts_csv_sha256` | string | `e1c0aa19eddf475812c2faa473e1f91be4dbb29b05a51a6ef1105ffa7c5aa975` |
| `historical_runs[3].native_output_interval_s` | number | `1e-06` |
| `historical_runs[3].trajectory_csv_rows` | number | `40002` |
| `historical_runs[3].contact_csv_rows` | number | `60000` |
| `published_report_md_sha256` | string | `85f2a17cf104e22d02c8ffb949b2e2df7c46eb664376c22342e94eb707efc903` |
| `report_transformation` | string | `Historical prose and numerical results preserved; add historical-scope banner, update relative evidence links, append native-rerun and compact-evidence links.` |

</details>
