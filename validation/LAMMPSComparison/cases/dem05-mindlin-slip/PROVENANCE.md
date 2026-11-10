# Historical provenance: dem05-mindlin-slip

Historical outcome: **FAIL**. This is not a fresh acceptance result.

Official source: [LAMMPS reference](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem05-mindlin-slip-3d-si.yaml).

Original numerical-record SHA256: `82bfc36ac50ed3259546e6aa2364917d33c72d6ca577152f96a99458b29fa2b6`.
Original provenance-record SHA256: `b142c61a03c4c2d2ae25456192218d04afcdbea900825780ecaec7d276d07336`.
Original report SHA256 before link/image cleanup: `18c0093b1b3a3435a8982e2510aad09736cb1e78cc692f5e732a9a1466a74f97`.
Archived report SHA256 (HISTORICAL-REPORT.md): `27cb8a4ef2506abd2d49e170b2e2a2701b51f6ad0b163a4cf56432166f46334c`.

Current report SHA256: `89984a260a5499cd90a6b68da6bc5442b7cfb4869296fc149e61372c223f99b7`.

The table preserves the original provenance fields and numeric tokens. Historical format/hash statements refer to the archived original records; omitted artifact locations are not runtime inputs.

[Case summary](REPORT.md) · [Historical report](HISTORICAL-REPORT.md) · [Numerical record](HISTORICAL-RESULTS.md)

<details>
<summary>Preserved provenance fields</summary>

| Field | Type | Value |
| --- | --- | --- |
| `case` | string | `dem05-mindlin-slip` |
| `reference_family` | string | `LAMMPS` |
| `historical_outcome` | string | `FAIL` |
| `historical_evidence_dates[0]` | string | `2026-10-28` |
| `historical_evidence_dates[1]` | string | `2026-10-29` |
| `migration_date` | string | `2026-10-30` |
| `status_policy` | string | `Historical comparison, not a new run or a current-core acceptance result.` |
| `reference_execution_policy` | string | `Use the case report: official saved LAMMPS states and analytical/public PFC targets; no local PFC execution is claimed.` |
| `official_demo_id` | null | `null` |
| `official_source` | string | `https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem05-mindlin-slip-3d-si.yaml` |
| `source_commit` | string | `f9720ce2d176d615495d1d6c31349b1f71be1c73` |
| `official_git_blob` | string | `864e4de3731c20e3ee0941e946439adfd40bfbbe` |
| `source_snapshot_sha256` | null | `null` |
| `official_criterion` | string | `1% final vx and omega_y` |
| `original_report_md_sha256` | string | `4324459397b6f2094ba6bee6935846999ec37a68aa1eb467bc46afe074b3896e` |
| `original_analysis_sha256` | string | `5e9986d9a1d5945734b309b95b529bdab7d41a7a2f4be731b1ec54a5824191f3` |
| `published_analysis_sha256` | string | `82bfc36ac50ed3259546e6aa2364917d33c72d6ca577152f96a99458b29fa2b6` |
| `analysis_transformation` | string | `Numerical JSON retained; personal absolute artifact paths replaced with historical-archive labels. Full raw histories are not distributed.` |
| `compact_checkpoint_tables[0].file` | string | `checkpoints/official_checkpoint_comparison.csv` |
| `compact_checkpoint_tables[0].sha256` | string | `3f8ecbdc78fea97470860c52e6b125f2d36753a56328d5a698fb27cc8488b771` |
| `compact_checkpoint_tables[0].bytes` | number | `5032` |
| `historical_runs[0].run` | string | `cpu` |
| `historical_runs[0].steps` | number | `75000` |
| `historical_runs[0].stored_frames` | number | `751` |
| `historical_runs[0].status` | string | `completed` |
| `historical_runs[0].passes_official_analytic_tolerance` | boolean | `false` |
| `historical_runs[0].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[0].run_json_sha256` | string | `468fd47c296db68b99653998ac1448af4dc6334927e00a52ab51864d69689fe4` |
| `historical_runs[0].project_snapshot_sha256` | string | `73e13a31c4246f65cdda1ceca15233ad5578be9b2fcc5ae3840762bbf6bb567c` |
| `historical_runs[0].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].trajectory_csv_sha256` | string | `de5a3d7721567d4ac4716f3b8cdedd6d8f0b259941de9d1c3e2c3502a85d551b` |
| `historical_runs[0].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[0].contacts_csv_sha256` | string | `888fd6ba6ef06b99b1b325eb705c3ef3ee0f699f6f0ee8f96202dc352e2e79c8` |
| `historical_runs[1].run` | string | `gpu` |
| `historical_runs[1].steps` | number | `75000` |
| `historical_runs[1].stored_frames` | number | `751` |
| `historical_runs[1].status` | string | `completed` |
| `historical_runs[1].passes_official_analytic_tolerance` | boolean | `false` |
| `historical_runs[1].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[1].run_json_sha256` | string | `7a60136d072f6a098f41fc2b0968b7c48b5617627fc572061e63cfebfa1a13f0` |
| `historical_runs[1].project_snapshot_sha256` | string | `c57ed46cd1df3a3750ff720c52a0b3081fa473d110ae8379c1ef1f9f12c1517a` |
| `historical_runs[1].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].trajectory_csv_sha256` | string | `117543a9c034f21710b4c7cb19cb9eaa8db7161f359868cc886f092b7cb3d2ae` |
| `historical_runs[1].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[1].contacts_csv_sha256` | string | `bc9ba9b1e2615fbaa938bd288256c7d941567d8464973f60411b6426ae4eea3a` |
| `historical_runs[2].run` | string | `cpu-dt-half` |
| `historical_runs[2].steps` | number | `150000` |
| `historical_runs[2].stored_frames` | number | `751` |
| `historical_runs[2].status` | string | `completed` |
| `historical_runs[2].passes_official_analytic_tolerance` | boolean | `false` |
| `historical_runs[2].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[2].run_json_sha256` | string | `018379c8743c7b96891ca71184eceb998daf40edcc30b9412c823a1754f83507` |
| `historical_runs[2].project_snapshot_sha256` | string | `a680b1284a1e0a396f93211dbc760f783f59ed26d6271cc8ee0ce30010542655` |
| `historical_runs[2].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].trajectory_csv_sha256` | string | `65a81dadaf2050d23da4e9f21b3c0548f4fda2c4b7cd1946abced3f796907344` |
| `historical_runs[2].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[2].contacts_csv_sha256` | string | `883a98fcf67c0b426aa2270db874e9ea49393e5342c7fc5ad124d6c8b77ff4a4` |
| `historical_runs[3].run` | string | `gpu-dt-half` |
| `historical_runs[3].steps` | number | `150000` |
| `historical_runs[3].stored_frames` | number | `751` |
| `historical_runs[3].status` | string | `completed` |
| `historical_runs[3].passes_official_analytic_tolerance` | boolean | `false` |
| `historical_runs[3].result_directory` | string | `Omitted historical artifact location.` |
| `historical_runs[3].run_json_sha256` | string | `06a5033e78305ff4d0836262f1234655d532d816d653a6712ad202b6f6408deb` |
| `historical_runs[3].project_snapshot_sha256` | string | `58127c38114b8e3d3df7b61bc310a30df9d684b0cb1a9d775a4fea276c180ee6` |
| `historical_runs[3].trajectory_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].trajectory_csv_sha256` | string | `c243c3a86b7d32fafa2b907e6d5d6e3b9eb9e31960c4f629c8713d9bd1ca9729` |
| `historical_runs[3].contacts_csv` | string | `Omitted historical artifact location.` |
| `historical_runs[3].contacts_csv_sha256` | string | `54177c245e41f7dd431c1f505057c4e933c2d1928161bad94dda05f09810ca45` |
| `published_report_md_sha256` | string | `18c0093b1b3a3435a8982e2510aad09736cb1e78cc692f5e732a9a1466a74f97` |
| `report_transformation` | string | `Historical prose and numerical results preserved; add historical-scope banner, update relative evidence links, append native-rerun and compact-evidence links.` |

</details>
