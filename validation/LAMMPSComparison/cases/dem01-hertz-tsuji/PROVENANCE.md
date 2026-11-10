# Historical provenance: dem01-hertz-tsuji

Historical outcome: **PASS**. This is not a fresh acceptance result.

Official source: [LAMMPS reference](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem01-hertz-tsuji-3d-si.yaml).

Original numerical-record SHA256: `2d6b927d855883e8aa05a7aced941f509f8605f9775715cc69ea1a34dddd578c`.
Original provenance-record SHA256: `399a5d5b1578813cdaf651eebd7022e8f3e7a78e88ed53dd6c508ec01df49b70`.
Original report SHA256 before link/image cleanup: `8d52101ab48211725acd2586c82519d869b2786fc978dfa7a639dd84c0b47627`.
Archived report SHA256 (HISTORICAL-REPORT.md): `a8963a453bdda1a5a633385fea896b573f70764aa169731f0918068e358813a0`.

Current report SHA256: `47f9b74411eb9c8447bc1cca5f0d45eeb4fc288cd11d17839a20d607d717376e`.

The table preserves the original provenance fields and numeric tokens. Historical format/hash statements refer to the archived original records; omitted artifact locations are not runtime inputs.

[Case summary](REPORT.md) · [Historical report](HISTORICAL-REPORT.md) · [Numerical record](HISTORICAL-RESULTS.md)

<details>
<summary>Preserved provenance fields</summary>

| Field | Type | Value |
| --- | --- | --- |
| `case` | string | `dem01-hertz-tsuji` |
| `reference_family` | string | `LAMMPS` |
| `historical_outcome` | string | `PASS` |
| `historical_evidence_dates[0]` | string | `2026-10-28` |
| `historical_evidence_dates[1]` | string | `2026-10-29` |
| `migration_date` | string | `2026-10-30` |
| `status_policy` | string | `Historical comparison, not a new run or a current-core acceptance result.` |
| `reference_execution_policy` | string | `Use the case report: official saved LAMMPS states and analytical/public PFC targets; no local PFC execution is claimed.` |
| `official_demo_id` | string | `dem01-hertz-tsuji` |
| `official_source` | string | `https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem01-hertz-tsuji-3d-si.yaml` |
| `source_commit` | string | `29dd3fd7642a6fe9e8b4dd99120df4068b221250` |
| `official_git_blob` | string | `e908cea7546be4d6c0579269aef5a782b4104bd1` |
| `source_snapshot_sha256` | string | `08b0692502f99e62e2cb419faa76d730b513412a50444f4956720b5c410cc3fa` |
| `official_criterion` | string | `Official collision_restitution selected-time analytic assertion at 4.1e-05 s; per-quantity error types/limits in analysis; source tags granular, unvetted` |
| `original_report_md_sha256` | string | `0ca54ff66544990e0e5733529cabcf233ab1c801ad176c269646fd1336465de5` |
| `original_analysis_sha256` | string | `4e6953353fdf48d2776d589f564780b695a1ccf8921e1f946d0d03e8028c0b4a` |
| `published_analysis_sha256` | string | `2d6b927d855883e8aa05a7aced941f509f8605f9775715cc69ea1a34dddd578c` |
| `analysis_transformation` | string | `Numerical JSON retained; personal absolute artifact paths replaced with historical-archive labels. Full raw histories are not distributed.` |
| `compact_checkpoint_tables[0].file` | string | `checkpoints/official_metrics_comparison.csv` |
| `compact_checkpoint_tables[0].sha256` | string | `db0671aa8c290bd754222daf9e8320db06877ca61028bcd50ada8eff344166f9` |
| `compact_checkpoint_tables[0].bytes` | number | `1007` |
| `compact_checkpoint_tables[1].file` | string | `checkpoints/source_checkpoint_comparison.csv` |
| `compact_checkpoint_tables[1].sha256` | string | `c9bf282777d02f87fb958ecdb10f27208a9a64b03e44644363933a85428d9985` |
| `compact_checkpoint_tables[1].bytes` | number | `13778` |
| `historical_runs[0].run` | string | `cpu` |
| `historical_runs[0].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[0].steps` | number | `5100` |
| `historical_runs[0].stored_frames` | number | `5101` |
| `historical_runs[0].status` | string | `completed` |
| `historical_runs[0].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[0].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[0].run_json_sha256` | string | `b49b6014206de0547d46a925c422a8ab9d9f5ee22801def20660f5b83f853a01` |
| `historical_runs[0].project_snapshot_sha256` | string | `efcaf53b9e13b77fc72f08e5a60f18755d0194b00dcf595f48b053df23180c1d` |
| `historical_runs[0].native_output_interval_s` | number | `1e-08` |
| `historical_runs[0].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].trajectory_csv_sha256` | string | `07a5ea60cb8c15fc9fbc4e630955db8661ef91d08e7b2a004ecb2376302d4742` |
| `historical_runs[0].trajectory_csv_rows` | number | `10202` |
| `historical_runs[0].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].contacts_csv_sha256` | string | `19cdee5f1f2b2a15dc160e7b00f395bafe653f852d93e490f726a374fd477e04` |
| `historical_runs[0].contact_csv_rows` | number | `1251` |
| `historical_runs[1].run` | string | `gpu` |
| `historical_runs[1].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[1].steps` | number | `5100` |
| `historical_runs[1].stored_frames` | number | `5101` |
| `historical_runs[1].status` | string | `completed` |
| `historical_runs[1].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[1].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[1].run_json_sha256` | string | `59e573b3638d2a99e074e25ccccf8518681966f539f25d63cc3b5587d10e16c5` |
| `historical_runs[1].project_snapshot_sha256` | string | `5e77e398e010af82bbaa01f2c1f8c772bbd24a89750f5761347f079bd1326c01` |
| `historical_runs[1].native_output_interval_s` | number | `1e-08` |
| `historical_runs[1].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].trajectory_csv_sha256` | string | `98bed7edba245e90b29610635ef160f6e4eaca76babf2a56b5b5e3272c57b97a` |
| `historical_runs[1].trajectory_csv_rows` | number | `10202` |
| `historical_runs[1].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].contacts_csv_sha256` | string | `8fbd303bb314310cf79fef1c26da0631b4afcd951f3b8167a422f16432c341c3` |
| `historical_runs[1].contact_csv_rows` | number | `1251` |
| `historical_runs[2].run` | string | `cpu-dt-half` |
| `historical_runs[2].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[2].steps` | number | `10200` |
| `historical_runs[2].stored_frames` | number | `10201` |
| `historical_runs[2].status` | string | `completed` |
| `historical_runs[2].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[2].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[2].run_json_sha256` | string | `d5f435076a00a9af71f111583211cfb1619bc55033643d83e9acfcbd4aaa0a1b` |
| `historical_runs[2].project_snapshot_sha256` | string | `431ccbba3d5212de233c5431557152c6e81689f752973fcbfb712392c2266ee7` |
| `historical_runs[2].native_output_interval_s` | number | `5e-09` |
| `historical_runs[2].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].trajectory_csv_sha256` | string | `c0ac95e1d604173a4faf1cf5e92665f92aeb18cdcebbc4f469bc9541a793f88d` |
| `historical_runs[2].trajectory_csv_rows` | number | `20402` |
| `historical_runs[2].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].contacts_csv_sha256` | string | `cf9995669bd2977d6ec30e94244394efa91f9bb83ca2d4d9dfd5d5d536958da3` |
| `historical_runs[2].contact_csv_rows` | number | `2502` |
| `historical_runs[3].run` | string | `gpu-dt-half` |
| `historical_runs[3].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[3].steps` | number | `10200` |
| `historical_runs[3].stored_frames` | number | `10201` |
| `historical_runs[3].status` | string | `completed` |
| `historical_runs[3].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[3].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[3].run_json_sha256` | string | `8e438d244e550c1185e6feb0ca2f7c2abd4e5247646018e93a218a870249603a` |
| `historical_runs[3].project_snapshot_sha256` | string | `9d882c3549ddadf2f7c400ed2c74435a654853d7982131d8fbdb3f17dc50d2c9` |
| `historical_runs[3].native_output_interval_s` | number | `5e-09` |
| `historical_runs[3].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].trajectory_csv_sha256` | string | `abd3e3069580b2bd5747f1c946e13c3f7f1fe5d7f40cad9b75f0740c21f92166` |
| `historical_runs[3].trajectory_csv_rows` | number | `20402` |
| `historical_runs[3].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].contacts_csv_sha256` | string | `36fc60c248595cac818b5923c92182a319bbdfeefbf2e46e4b30b22b965535bf` |
| `historical_runs[3].contact_csv_rows` | number | `2502` |
| `published_report_md_sha256` | string | `8d52101ab48211725acd2586c82519d869b2786fc978dfa7a639dd84c0b47627` |
| `report_transformation` | string | `Historical prose and numerical results preserved; add historical-scope banner, update relative evidence links, append native-rerun and compact-evidence links.` |

</details>
