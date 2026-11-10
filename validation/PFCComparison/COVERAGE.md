# Historical PFC coverage

This inventory covers the fifteen PFC 7.0 verification-index entries, checked
against the PFC 9.0 index on 29 October 2026. It is a record of that campaign's
configured source/application path, **not a current-core feature matrix**.
PFC was not executed locally. The three retained comparisons cover one official
example's Md=0 row; they do not imply fifteen reproduced examples.

| Official example | Campaign scope or exclusion |
| --- | --- |
| [Measure Logic](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/measure_logic/PorosityMeasurement.html) | Exact clipped sphere/circle-region porosity operator unavailable; nominal packing porosity is different. |
| [FCC Array Strength](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/ArrayStrength/ArrayStrength.html) | Finite-particle spin/component constraints, loads and PFC local damping unavailable on the configured path. |
| [Tip-Loaded Cantilever](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/Cantilever/CantileverBeam.html) | Nonzero reference stages require finite-body tip forces/moments; full parallel-bond equivalence unproven. |
| [Linear Restitution](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/linear_contact_model/cmlinear_restitution.html) | Sixteen Md=0 ratios in CPU/GPU at two steps: 64 historical PASS checks. Md=1 no-tension row excluded. |
| [Hertz Complex Loading](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/hertz_contact_model/cmhertz_single.html) | Prescribed finite-sphere path not representable; using two infinite-mass spheres changes the model and detection. |
| [Wave Propagation](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/WavePropagation/WaveProp.html) | Finite-mass half-cosine force driver unavailable; prescribed motion changes boundary inertia/impedance. |
| [Rolling: Flat Surface](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/rrlinear_contact_model/flat_surface/cmrrlinear_flatsurface.html) | Active PFC2D disk inertia I/(mR²)=1/2 differs from sphere 2/5; dp_mode=3 also differs. |
| [Rolling: Repose Angle](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/rrlinear_contact_model/repose_angle/cmrrlinear_reposeangle.html) | rrlinear dp_mode=3, elastic-force friction cap, slip-cut damping and calm/local-damping initialization differ. |
| [Adhesive Rolling](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/arrlinear_contact_model/repose_angle/cmarrlinear_reposeangle.html) | Positive-gap attraction and modified contact creation unavailable. |
| [Burger Stress Relaxation](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/burger_contact_model/stress_relaxation/cmburger_stressrelaxation.html) | Kelvin/Maxwell internal rheology unavailable; a Hooke spring/dashpot is not equivalent. |
| [Sliding Wedge](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/RigidBlocks/SlidingWedge/rigid_blocks_sliding_wedge.html) | Rigid-block contact, preload, local damping and restore/bisection protocol unavailable. |
| [Spring Network](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/springnetwork/springnetwork.html) | RBSN/soft-bond fictitious-stress correction, weakening and healing differ from B1–B4/BK bonds. |
| [Permanent Magnets](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/squareMagnet/squareMagnet.html) | Long-range magnetic dipole forces/moments unavailable. |
| [Cable Anchor Pull-Test](https://docs.itascacg.com/pfc700/pfc/rblocksel/test3d/cable/pullout/pullout_pfc.html) | Cable/grout structural elements and rigid-block coupling unavailable. |
| [Hybrid Bolt Shear](https://docs.itascacg.com/pfc700/pfc/rblocksel/test3d/hybrid/shear_pull/shear_pull.html) | Hybrid cable/dowel, grout, plastic yield and rupture operators unavailable. |

The normal comparison preserves the source law and inputs. For the sweep,
`e_ball=e_contact/(2−e_contact)` exactly inverts the harmonic mixing with a floor
restitution of 1; it is not empirical fitting. The critical β=1 branch uses
e_contact=0. Disk rotation is active in the excluded 2D example; changing density
cannot turn sphere inertia into disk inertia. PFC local damping depends on force
magnitude and velocity sign, unlike velocity-proportional viscous damping.

See the [sweep report](cases/pfc-linear-restitution-md0-sweep/REPORT.md),
[bond feasibility](diagnostics/audits/pfc-bond-feasibility-review.md),
[prescribed Hertz feasibility](diagnostics/audits/pfc-hertz-normal-feasibility.md)
and [original full coverage notes](HISTORICAL-COVERAGE.md). Original coverage
JSON, vendor source snapshots and raw archives remain separately retained.
