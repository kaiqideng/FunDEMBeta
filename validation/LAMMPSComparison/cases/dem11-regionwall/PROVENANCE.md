# Historical provenance: dem11-regionwall

Historical outcome: **PASS**. This is not a fresh acceptance result.

Official source: [LAMMPS reference](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem11-regionwall-3d-si.yaml).

Original numerical-record SHA256: `6313b05cba8dedc46505d78d4325a5b85846e6752ed6c1f5a8b5deb7042280ce`.
Original provenance-record SHA256: `5df6e1eceeeb8e346f4cd323fddfca6d12c020f59a95f66508f3e2962cd01cca`.
Original report SHA256 before link/image cleanup: `cce699053d2e125069a7484481aafc75ba28458ece83885f162c4bef1e95e8ce`.
Archived report SHA256 (HISTORICAL-REPORT.md): `ab3faeb4f33dbfef0e627859d8c473a2fc5b1ad03515f242719f356a55be94ae`.

Current report SHA256: `5e32d5c8fbfe21602b7beae15053789f758fa1d1152981ac655d97d6d1966a5a`.

The table preserves the original provenance fields and numeric tokens. Historical format/hash statements refer to the archived original records; omitted artifact locations are not runtime inputs.

[Case summary](REPORT.md) · [Historical report](HISTORICAL-REPORT.md) · [Numerical record](HISTORICAL-RESULTS.md)

<details>
<summary>Preserved provenance fields</summary>

| Field | Type | Value |
| --- | --- | --- |
| `case` | string | `dem11-regionwall` |
| `reference_family` | string | `LAMMPS` |
| `historical_outcome` | string | `PASS` |
| `historical_evidence_dates[0]` | string | `2026-10-28` |
| `historical_evidence_dates[1]` | string | `2026-10-29` |
| `migration_date` | string | `2026-10-30` |
| `status_policy` | string | `Historical comparison, not a new run or a current-core acceptance result.` |
| `reference_execution_policy` | string | `Use the case report: official saved LAMMPS states and analytical/public PFC targets; no local PFC execution is claimed.` |
| `official_demo_id` | string | `dem11-regionwall` |
| `official_source` | string | `https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem11-regionwall-3d-si.yaml` |
| `source_commit` | string | `29dd3fd7642a6fe9e8b4dd99120df4068b221250` |
| `official_git_blob` | string | `e90c87ee07bed0d26ab25bb04320a43f73ef1b4c` |
| `source_snapshot_sha256` | string | `1dfd8202721fbad97cfc800734790168c1677f9a489f0c77b9022613a09f5d63` |
| `official_criterion` | string | `Official wall_restitution selected-time analytic assertion at 0.000205 s; per-quantity error types/limits in analysis; source tags granular, unvetted` |
| `original_report_md_sha256` | string | `1c838311daf503ac61f872fe40c2a4d09b69c627b2d0308139a54c9592b983e1` |
| `original_analysis_sha256` | string | `6c2ec769396bcaf5f06ed827b32270ffc25979b8560ef2eac6ff7a7999874370` |
| `published_analysis_sha256` | string | `6313b05cba8dedc46505d78d4325a5b85846e6752ed6c1f5a8b5deb7042280ce` |
| `analysis_transformation` | string | `Numerical JSON retained; personal absolute artifact paths replaced with historical-archive labels. Full raw histories are not distributed.` |
| `compact_checkpoint_tables[0].file` | string | `checkpoints/official_metrics_comparison.csv` |
| `compact_checkpoint_tables[0].sha256` | string | `1f1dba8f49fd90f4df51031752b462e7299a7f35447b4b37221742461ab7d276` |
| `compact_checkpoint_tables[0].bytes` | number | `613` |
| `compact_checkpoint_tables[1].file` | string | `checkpoints/source_checkpoint_comparison.csv` |
| `compact_checkpoint_tables[1].sha256` | string | `9d3fb8ea6e5b51106dcb07827e788de09f67da9bf1b04226967609c77c7ed684` |
| `compact_checkpoint_tables[1].bytes` | number | `3862` |
| `historical_runs[0].run` | string | `cpu` |
| `historical_runs[0].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[0].steps` | number | `20500` |
| `historical_runs[0].stored_frames` | number | `206` |
| `historical_runs[0].status` | string | `completed` |
| `historical_runs[0].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[0].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[0].run_json_sha256` | string | `c4a5eab31f5ec0825ccf1995a0472be85b98c545b1bbdb1443ab45ee83a9c46f` |
| `historical_runs[0].project_snapshot_sha256` | string | `c1ebce978f4c979b417e8643e5e4e13dcdb031445d6374052a3529edfc23c607` |
| `historical_runs[0].native_output_interval_s` | number | `1e-06` |
| `historical_runs[0].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].trajectory_csv_sha256` | string | `bd9321b69ba6bcae8926244fdcf58c2fa95c82ecc5a4c09068d13352405e7cb4` |
| `historical_runs[0].trajectory_csv_rows` | number | `20501` |
| `historical_runs[0].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].contacts_csv_sha256` | string | `aa710f299a4f706faa0554aedb3fb4dded3716765bc235fc55a92fce636c3f45` |
| `historical_runs[0].contact_csv_rows` | number | `2083` |
| `historical_runs[1].run` | string | `gpu` |
| `historical_runs[1].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[1].steps` | number | `20500` |
| `historical_runs[1].stored_frames` | number | `206` |
| `historical_runs[1].status` | string | `completed` |
| `historical_runs[1].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[1].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[1].run_json_sha256` | string | `dc813e40ebcf2106073763bf54ce1d41a625ed6a8502dbd51323b5365b97a56e` |
| `historical_runs[1].project_snapshot_sha256` | string | `abd850be4c18d359ecf6597db7d74430e7f3260511976b955d54fa947f4cd403` |
| `historical_runs[1].native_output_interval_s` | number | `1e-06` |
| `historical_runs[1].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].trajectory_csv_sha256` | string | `9f221e4ad809cb20a0473dab24fc8fac7607324a311232e36f74becf46b8f734` |
| `historical_runs[1].trajectory_csv_rows` | number | `20501` |
| `historical_runs[1].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].contacts_csv_sha256` | string | `e7870093c95119a2d3203c4e3e0797f3fdf91aa6f3b78b77269c331792e885b4` |
| `historical_runs[1].contact_csv_rows` | number | `2083` |
| `historical_runs[2].run` | string | `cpu-dt-half` |
| `historical_runs[2].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[2].steps` | number | `41000` |
| `historical_runs[2].stored_frames` | number | `206` |
| `historical_runs[2].status` | string | `completed` |
| `historical_runs[2].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[2].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[2].run_json_sha256` | string | `3e8ac8b960566bae8c964b49cf204cc39072b6644828f59802beaa2ec73fb07d` |
| `historical_runs[2].project_snapshot_sha256` | string | `160e0ed9083652df2c64fa4ed5882156fe03bc5e38d6fcfefabed4ecf6b57655` |
| `historical_runs[2].native_output_interval_s` | number | `1e-06` |
| `historical_runs[2].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].trajectory_csv_sha256` | string | `487e636d704dc22b43b9999a75f70271ecdab12bf01b86ed2f55380f6d9d051e` |
| `historical_runs[2].trajectory_csv_rows` | number | `41001` |
| `historical_runs[2].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].contacts_csv_sha256` | string | `8c5e9f2fe00a5930a6443102bda849d591ff028da071fed8402065bef26df40b` |
| `historical_runs[2].contact_csv_rows` | number | `4166` |
| `historical_runs[3].run` | string | `gpu-dt-half` |
| `historical_runs[3].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[3].steps` | number | `41000` |
| `historical_runs[3].stored_frames` | number | `206` |
| `historical_runs[3].status` | string | `completed` |
| `historical_runs[3].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[3].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[3].run_json_sha256` | string | `9319a7ac4b9e71ac05feb1140cb65a9896e84b3d2fa9a847570a34d506f85df8` |
| `historical_runs[3].project_snapshot_sha256` | string | `809b7deadc21f01da22e23c0e2918f06c31a6c947d4d3fa2b81bf1b441d2cc7e` |
| `historical_runs[3].native_output_interval_s` | number | `1e-06` |
| `historical_runs[3].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].trajectory_csv_sha256` | string | `c5c9a3d59b7adde2c9281189b98ff8700b31b68769f2c30e9e294aad127aeca5` |
| `historical_runs[3].trajectory_csv_rows` | number | `41001` |
| `historical_runs[3].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].contacts_csv_sha256` | string | `6251d9498a6b649d816cec34bf16be5b3238b929d36e80a858dea4f49ac4fe85` |
| `historical_runs[3].contact_csv_rows` | number | `4166` |
| `published_report_md_sha256` | string | `cce699053d2e125069a7484481aafc75ba28458ece83885f162c4bef1e95e8ce` |
| `report_transformation` | string | `Historical prose and numerical results preserved; add historical-scope banner, update relative evidence links, append native-rerun and compact-evidence links.` |

</details>
