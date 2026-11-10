# Historical provenance: dem03-hertz-peak-wall

Historical outcome: **PASS**. This is not a fresh acceptance result.

Official source: [LAMMPS reference](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem03-hertz-peak-wall-3d-si.yaml).

Original numerical-record SHA256: `6c9973b9bbd67a2029791c9256c67dec6eb303de9c88b703bdf0f38eade813a7`.
Original provenance-record SHA256: `1bde88e838d5690c7856609b3cf979ac6f530bd0702a0f2d42701a2215163125`.
Original report SHA256 before link/image cleanup: `7d716c769367d41e395b3af66a32ad51597386d4e64cce67d5bffc828ac0c290`.
Archived report SHA256 (HISTORICAL-REPORT.md): `04f347633b180c6d5d83be905bfab5a416827767701f201e1f38384a44e7a226`.

Current report SHA256: `76a506c4aeda45538a7ae41b680c91e682ce2da74fe1be1068d5b3ec80af6d82`.

The table preserves the original provenance fields and numeric tokens. Historical format/hash statements refer to the archived original records; omitted artifact locations are not runtime inputs.

[Case summary](REPORT.md) · [Historical report](HISTORICAL-REPORT.md) · [Numerical record](HISTORICAL-RESULTS.md)

<details>
<summary>Preserved provenance fields</summary>

| Field | Type | Value |
| --- | --- | --- |
| `case` | string | `dem03-hertz-peak-wall` |
| `reference_family` | string | `LAMMPS` |
| `historical_outcome` | string | `PASS` |
| `historical_evidence_dates[0]` | string | `2026-10-28` |
| `historical_evidence_dates[1]` | string | `2026-10-29` |
| `migration_date` | string | `2026-10-30` |
| `status_policy` | string | `Historical comparison, not a new run or a current-core acceptance result.` |
| `reference_execution_policy` | string | `Use the case report: official saved LAMMPS states and analytical/public PFC targets; no local PFC execution is claimed.` |
| `official_demo_id` | string | `dem03-hertz-peak-wall` |
| `official_source` | string | `https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem03-hertz-peak-wall-3d-si.yaml` |
| `source_commit` | string | `29dd3fd7642a6fe9e8b4dd99120df4068b221250` |
| `official_git_blob` | string | `531d86bb3e2993d869041734a6dce426e3737758` |
| `source_snapshot_sha256` | string | `7b2405044d4a79f791eb5b79eb4564ee1e4d9b5ee6be6ba09271ce49d184b26b` |
| `official_criterion` | string | `Official hertz_peak selected-time analytic assertion at 2.138e-05 s; per-quantity error types/limits in analysis; source tags granular, unvetted` |
| `original_report_md_sha256` | string | `734ba16a33c4fc96a430c6df062f262fd273236eaf6c3b851126d61199338a2d` |
| `original_analysis_sha256` | string | `9a430f8d56092829ae5377ffeb43f26570cdd4fe780451343e749106e55d72ae` |
| `published_analysis_sha256` | string | `6c9973b9bbd67a2029791c9256c67dec6eb303de9c88b703bdf0f38eade813a7` |
| `analysis_transformation` | string | `Numerical JSON retained; personal absolute artifact paths replaced with historical-archive labels. Full raw histories are not distributed.` |
| `compact_checkpoint_tables[0].file` | string | `checkpoints/official_metrics_comparison.csv` |
| `compact_checkpoint_tables[0].sha256` | string | `8c1745e5ded49f8d742749b0a1cd354a5880997944246b46805c776063c5c8d7` |
| `compact_checkpoint_tables[0].bytes` | number | `1363` |
| `compact_checkpoint_tables[1].file` | string | `checkpoints/source_checkpoint_comparison.csv` |
| `compact_checkpoint_tables[1].sha256` | string | `90ac29f2fc4d1a30d35cb6ff4d11bdd3e22886769abfb9d12dec82881a1b75e4` |
| `compact_checkpoint_tables[1].bytes` | number | `5570` |
| `historical_runs[0].run` | string | `cpu` |
| `historical_runs[0].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[0].steps` | number | `3038` |
| `historical_runs[0].stored_frames` | number | `3039` |
| `historical_runs[0].status` | string | `completed` |
| `historical_runs[0].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[0].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[0].run_json_sha256` | string | `15590f35219c6eb5f88803aa7fdb9902bc8d1f08affbfb57a55caeca68773c56` |
| `historical_runs[0].project_snapshot_sha256` | string | `02b90cc028ca3ef06e8567f4c0f6a7175623de701ae81c446d1d7dca5d8f898b` |
| `historical_runs[0].native_output_interval_s` | number | `1e-08` |
| `historical_runs[0].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].trajectory_csv_sha256` | string | `c6237df4a3ce3ee366a5fe01fb80b87757e1feb2f5d117085e85bd821c602f84` |
| `historical_runs[0].trajectory_csv_rows` | number | `3039` |
| `historical_runs[0].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].contacts_csv_sha256` | string | `3295f7702203f7bf7c3c68eb027cf28af2c26c318750b39e87a15418ab221194` |
| `historical_runs[0].contact_csv_rows` | number | `1711` |
| `historical_runs[1].run` | string | `gpu` |
| `historical_runs[1].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[1].steps` | number | `3038` |
| `historical_runs[1].stored_frames` | number | `3039` |
| `historical_runs[1].status` | string | `completed` |
| `historical_runs[1].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[1].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[1].run_json_sha256` | string | `bd1f30c131356ad6f76011f76a9ba459f43d23eb2c463c35db8ab46545e1ef5d` |
| `historical_runs[1].project_snapshot_sha256` | string | `1cd41105b94f7e370fbce91734fedb43f6b06b17cb7009bcb1a8a72335de27de` |
| `historical_runs[1].native_output_interval_s` | number | `1e-08` |
| `historical_runs[1].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].trajectory_csv_sha256` | string | `dcb5671f591990522ce925df3e5f25d5d8b0b8817a18806a3e92b33fadb09e46` |
| `historical_runs[1].trajectory_csv_rows` | number | `3039` |
| `historical_runs[1].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].contacts_csv_sha256` | string | `3295f7702203f7bf7c3c68eb027cf28af2c26c318750b39e87a15418ab221194` |
| `historical_runs[1].contact_csv_rows` | number | `1711` |
| `historical_runs[2].run` | string | `cpu-dt-half` |
| `historical_runs[2].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[2].steps` | number | `6076` |
| `historical_runs[2].stored_frames` | number | `6077` |
| `historical_runs[2].status` | string | `completed` |
| `historical_runs[2].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[2].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[2].run_json_sha256` | string | `3efa5a56fb61f20e661ee260ceb2f8306a9ff1627fa4502e09759cb0cbc23264` |
| `historical_runs[2].project_snapshot_sha256` | string | `984701d5c3880e74da7c3a2b21f1345bf8e31de2311f31aa29a078ef29077cb5` |
| `historical_runs[2].native_output_interval_s` | number | `5e-09` |
| `historical_runs[2].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].trajectory_csv_sha256` | string | `57dfb666d947d76e8849c206d6240387700f5ab48f6d0cdd3a8512506fa55c7e` |
| `historical_runs[2].trajectory_csv_rows` | number | `6077` |
| `historical_runs[2].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].contacts_csv_sha256` | string | `81acc947c1d49b7758546cd8a469c88c14e95c2e88d4ee37d36bc45b787e3712` |
| `historical_runs[2].contact_csv_rows` | number | `3423` |
| `historical_runs[3].run` | string | `gpu-dt-half` |
| `historical_runs[3].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[3].steps` | number | `6076` |
| `historical_runs[3].stored_frames` | number | `6077` |
| `historical_runs[3].status` | string | `completed` |
| `historical_runs[3].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[3].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[3].run_json_sha256` | string | `e2c00bd24b23fa840abd6b2d622ff6680c3e062a70325546a620e134ce4579bb` |
| `historical_runs[3].project_snapshot_sha256` | string | `01ba5890ecec994149a7656a20093a28930dea34f742197adcb79503b63d941b` |
| `historical_runs[3].native_output_interval_s` | number | `5e-09` |
| `historical_runs[3].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].trajectory_csv_sha256` | string | `99fe8fc7910862863cd60e6c9863df7ef62f9669fe89e7265f4e72ef200205c3` |
| `historical_runs[3].trajectory_csv_rows` | number | `6077` |
| `historical_runs[3].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].contacts_csv_sha256` | string | `81acc947c1d49b7758546cd8a469c88c14e95c2e88d4ee37d36bc45b787e3712` |
| `historical_runs[3].contact_csv_rows` | number | `3423` |
| `published_report_md_sha256` | string | `7d716c769367d41e395b3af66a32ad51597386d4e64cce67d5bffc828ac0c290` |
| `report_transformation` | string | `Historical prose and numerical results preserved; add historical-scope banner, update relative evidence links, append native-rerun and compact-evidence links.` |

</details>
