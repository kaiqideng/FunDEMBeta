# Historical provenance: dem01-hooke-veldamp

Historical outcome: **PASS**. This is not a fresh acceptance result.

Official source: [LAMMPS reference](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem01-hooke-veldamp-3d-si.yaml).

Original numerical-record SHA256: `0f28fc0d08c5d03696a4aa0c2fa57a99496779d71daefe620065b1a22f823395`.
Original provenance-record SHA256: `d3cab7106383d8dc169b12f500e1f43cb819ae11d724dd9a326a12955ec0482b`.
Original report SHA256 before link/image cleanup: `b69305837e0c33d3b554689bef290e40140babf810dfce23b16fbae871d90beb`.
Archived report SHA256 (HISTORICAL-REPORT.md): `e5296a8f05b31815665a369537f7c9d7a0df1fb6f20cddd55b3c01f2cbe2f0c7`.

Current report SHA256: `553eb407d606a6f8d60f33b7f96afc74eeed0f9d033dc62e33cdc637d060ff6e`.

The table preserves the original provenance fields and numeric tokens. Historical format/hash statements refer to the archived original records; omitted artifact locations are not runtime inputs.

[Case summary](REPORT.md) · [Historical report](HISTORICAL-REPORT.md) · [Numerical record](HISTORICAL-RESULTS.md)

<details>
<summary>Preserved provenance fields</summary>

| Field | Type | Value |
| --- | --- | --- |
| `case` | string | `dem01-hooke-veldamp` |
| `reference_family` | string | `LAMMPS` |
| `historical_outcome` | string | `PASS` |
| `historical_evidence_dates[0]` | string | `2026-10-28` |
| `historical_evidence_dates[1]` | string | `2026-10-29` |
| `migration_date` | string | `2026-10-30` |
| `status_policy` | string | `Historical comparison, not a new run or a current-core acceptance result.` |
| `reference_execution_policy` | string | `Use the case report: official saved LAMMPS states and analytical/public PFC targets; no local PFC execution is claimed.` |
| `official_demo_id` | string | `dem01-hooke-veldamp` |
| `official_source` | string | `https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem01-hooke-veldamp-3d-si.yaml` |
| `source_commit` | string | `29dd3fd7642a6fe9e8b4dd99120df4068b221250` |
| `official_git_blob` | string | `18dcda1ac6edf4c667c9bcbd5283be1f0455d2e2` |
| `source_snapshot_sha256` | string | `689481bef71bc7f1966f9b47afe3f8398604aaeeaeaa63f98df891f2cd06a444` |
| `official_criterion` | string | `Official collision_restitution selected-time analytic assertion at 8.385e-05 s; per-quantity error types/limits in analysis; source tags granular, unvetted` |
| `original_report_md_sha256` | string | `35283dad8701ddfa2e7b4df557d527dc804d5d209b19826a5519c47b583a5b9e` |
| `original_analysis_sha256` | string | `4e754312ae6d04ccae404cf54bf4c7a731e2517185d9eaf55cd736c2d6430e22` |
| `published_analysis_sha256` | string | `0f28fc0d08c5d03696a4aa0c2fa57a99496779d71daefe620065b1a22f823395` |
| `analysis_transformation` | string | `Numerical JSON retained; personal absolute artifact paths replaced with historical-archive labels. Full raw histories are not distributed.` |
| `compact_checkpoint_tables[0].file` | string | `checkpoints/official_metrics_comparison.csv` |
| `compact_checkpoint_tables[0].sha256` | string | `3cde993cc9ec48995d970ccc305c21290c9b279e963d54f14829779e356c579b` |
| `compact_checkpoint_tables[0].bytes` | number | `1015` |
| `compact_checkpoint_tables[1].file` | string | `checkpoints/source_checkpoint_comparison.csv` |
| `compact_checkpoint_tables[1].sha256` | string | `b31e31117a0db924e857ea840f5c7d7484c5fbb585ea0e2f01b787b4fcbef10e` |
| `compact_checkpoint_tables[1].bytes` | number | `7033` |
| `historical_runs[0].run` | string | `cpu` |
| `historical_runs[0].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[0].steps` | number | `2795` |
| `historical_runs[0].stored_frames` | number | `2796` |
| `historical_runs[0].status` | string | `completed` |
| `historical_runs[0].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[0].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[0].run_json_sha256` | string | `b294060f7226dc763ca98fc3049c5758925774f5583e95664bc1cc42770f767c` |
| `historical_runs[0].project_snapshot_sha256` | string | `030c0dc206208fe4252d52409733a06b746d2c52aeda16314452c2cf2e8cb767` |
| `historical_runs[0].native_output_interval_s` | number | `3e-08` |
| `historical_runs[0].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].trajectory_csv_sha256` | string | `6f553e936ea837202e7d96f955d58dff6a63327e254a6cddbd56007fe5ace5c5` |
| `historical_runs[0].trajectory_csv_rows` | number | `5592` |
| `historical_runs[0].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].contacts_csv_sha256` | string | `c63394b7606abff1fbca95b936b0dc5cf37f2ff98f823fd5370cf9b52decc598` |
| `historical_runs[0].contact_csv_rows` | number | `1008` |
| `historical_runs[1].run` | string | `gpu` |
| `historical_runs[1].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[1].steps` | number | `2795` |
| `historical_runs[1].stored_frames` | number | `2796` |
| `historical_runs[1].status` | string | `completed` |
| `historical_runs[1].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[1].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[1].run_json_sha256` | string | `51918716ac9639e790531abd83ff68a84e500c4313b9c5202dac24af2db10623` |
| `historical_runs[1].project_snapshot_sha256` | string | `edacde67b649423479850bb4b0834c4a69ad9be2dcb863c46f5889f9a20f314a` |
| `historical_runs[1].native_output_interval_s` | number | `3e-08` |
| `historical_runs[1].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].trajectory_csv_sha256` | string | `7c439d4e374f860ac4187b968a28ab0ab3cf60d7a2f0f491c78b0bddda16c1b0` |
| `historical_runs[1].trajectory_csv_rows` | number | `5592` |
| `historical_runs[1].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].contacts_csv_sha256` | string | `f069328493b2b5a05715ce9cefb84450522e0570b5d0c534b4a483b2ef61b596` |
| `historical_runs[1].contact_csv_rows` | number | `1008` |
| `historical_runs[2].run` | string | `cpu-dt-half` |
| `historical_runs[2].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[2].steps` | number | `5590` |
| `historical_runs[2].stored_frames` | number | `5591` |
| `historical_runs[2].status` | string | `completed` |
| `historical_runs[2].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[2].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[2].run_json_sha256` | string | `3a8d0dff975d042518885b4d5f3ef2cfc4899d9ab7ba2ba9b52c551fe470d51a` |
| `historical_runs[2].project_snapshot_sha256` | string | `a17d835d370ff2d37eac3a8c9a8005aa42aa6e39c5857c376507c0f4520fa62e` |
| `historical_runs[2].native_output_interval_s` | number | `1.5e-08` |
| `historical_runs[2].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].trajectory_csv_sha256` | string | `b8829807686b828c50d2fbf1657cdf1d55e204bfd0f330a63c6ba7bd756ec6fa` |
| `historical_runs[2].trajectory_csv_rows` | number | `11182` |
| `historical_runs[2].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].contacts_csv_sha256` | string | `6148b5b9d93707008390525dc695d01902ab2f5772cba5618cc57611150d2db7` |
| `historical_runs[2].contact_csv_rows` | number | `2016` |
| `historical_runs[3].run` | string | `gpu-dt-half` |
| `historical_runs[3].role` | string | `main_cpu_gpu_timestep` |
| `historical_runs[3].steps` | number | `5590` |
| `historical_runs[3].stored_frames` | number | `5591` |
| `historical_runs[3].status` | string | `completed` |
| `historical_runs[3].passes_official_analytic_tolerance` | boolean | `true` |
| `historical_runs[3].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[3].run_json_sha256` | string | `dc931614c9f2ec9614bc67095e0b70adfea50954eed936da8bfcd515e58bcce5` |
| `historical_runs[3].project_snapshot_sha256` | string | `f3ae765adaf6ac04ef73e033388e7e9c03ff389ce0d2256218c415df51efe5da` |
| `historical_runs[3].native_output_interval_s` | number | `1.5e-08` |
| `historical_runs[3].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].trajectory_csv_sha256` | string | `2cbb18d4fdd5ed634231e6f415139a49d48d52b14512807fde37e8079496140a` |
| `historical_runs[3].trajectory_csv_rows` | number | `11182` |
| `historical_runs[3].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].contacts_csv_sha256` | string | `393fc35b72a0a2fa1e88d81fb8f699e76cc2badd105fd020509ce2f35d52f11f` |
| `historical_runs[3].contact_csv_rows` | number | `2016` |
| `published_report_md_sha256` | string | `b69305837e0c33d3b554689bef290e40140babf810dfce23b16fbae871d90beb` |
| `report_transformation` | string | `Historical prose and numerical results preserved; add historical-scope banner, update relative evidence links, append native-rerun and compact-evidence links.` |

</details>
