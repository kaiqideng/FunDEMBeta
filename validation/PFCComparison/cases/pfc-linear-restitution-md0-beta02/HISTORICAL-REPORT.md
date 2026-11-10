# PFC linear-contact restitution: verification report

> Historical evidence from 28–29 October 2026, migrated on 30 October 2026. The outcome and numerical results below are preserved, not a new run of the current core. Historical raw-frame/archive retention statements describe the original campaign, not files distributed here.

**Reference:** PFC3D 7.0 official linear-contact restitution example  
**Selected subcase:** original ball 8, dashpot mode Md=0, normal damping ratio beta_n=0.2  
**Date:** 29 October 2026  
**Published comparison target: PASS in the four retained runs.**

## Official source, scope and unchanged parameters

The [official example](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/linear_contact_model/cmlinear_restitution.html) contains the complete `restitution.p3dat` input. This verification selects **one independent original ball** with x=0.75 m, y=0.05 m, initial center z=1 m and radius 0.025 m from the Md=0 row. The other balls are independent parameter-sweep subcases. It does **not** claim the full 32-ball sweep, Md=1 no-tension clipping, PFC solver execution, or a new numbered LAMMPS demo.

The versioned PFC7 documentation text, embedded source and mapping audit are retained with SHA-256 identity. There is no invented Git commit or blob. Official-page web-connector text snapshot SHA-256: **3dc3744ccb9cc04fb3ec1db4a2d44f27f83b50b8e187020b03c54ab80f00e7e5**. [Source/contact audit](../../diagnostics/audits/pfc-linear-beta02-contact-audit.md).

| Quantity | Official input / native mapping |
|---|---|
| Particle | R=0.025 m; density=2650 kg/m^3; mass=0.173442094417 kg |
| Initial state | Center=(0.75,0.05,1) m; zero velocity and spin |
| Wall / gravity | Fixed z=0 plane; gravity=(0,0,-10) m/s^2 throughout |
| Normal contact | Linear kn=50000 N/m; beta_n=0.2; Md=0 full signed dashpot |
| Other channels | Zero tangential stiffness, friction/damping, rolling/twist resistance, local/global damping |
| Automatic damping | e_input=0.526620599330; derived cn=37.2496544324 N s/m; manual overrides off |
| Time | Official fixed dt=1e-5 s, 100000 steps to 1 s; half-step 5e-6 s, 200000 steps |
| Retention | Every numerical step in particle/contact CSV; native binary output every 1 ms, 1001 frames/run |

The initial center height is 1 m; the **free-fall clearance is 1-R=0.975 m**, not 1 m. Both sphere and wall material restitution are set from the published damping ratio. `linearStiffnessEnabled=true` preserves kn; `dampingEnabled=false` disables manual overrides and uses automatic contact damping. All stored manual dashpots are zero. Sphere-wall effective mass is the **full sphere mass**, not half mass. Tangential channels remain inactive.

## Analytic coefficient, published target and measured rebound

The official damped-oscillator relation is:

    e0 = exp[-beta_n pi / sqrt(1-beta_n^2)] = 0.526620599330
    beta_n = -ln(e0)/sqrt(pi^2+ln(e0)^2)
    cn = 2 beta_n sqrt(m kn)
    Fn = kn delta - cn vn_assembly

No coefficient was fitted and no solver source was changed. Md=0 retains the signed viscous unloading tail while geometrical contact exists. A negative normal force in that tail is not an adhesion law. **Md=1**, which clips the dashpot to prevent tensile total force, is excluded because it is not this native contact mode.

The official FISH height-apex method evaluates the rebound after first impact:

    e_height = sqrt[(z_first_rebound_apex - R)/(z0 - R)]
    coefficient_error = abs(e_height - e0)
    published comparison target: coefficient_error < 0.01

The page reports a maximum difference below 1%, but its `make_tables` expression actually computes **max(abs(e_num-e_ana))** over the Md=0 row. Thus this report uses an **absolute coefficient error below 0.01**, not relative 1%. The example does not supply an automated assertion; the target is the published comparison target applied to this one subcase. Reference e0 neglects gravity during contact, while the actual simulation retains gravity; a finite contact-time offset can remain under dt refinement.

| Run | Apex-derived e | Absolute coefficient error | First rebound apex z (m) | Published target |
|---|---:|---:|---:|---|
| CPU dt | 0.5240174395 | 0.0026032 | 0.292729420 | PASS |
| GPU dt | 0.5240174395 | 0.0026032 | 0.292729420 | PASS |
| CPU dt/2 | 0.5237242200 | 0.0028964 | 0.292429882 | PASS |
| GPU dt/2 | 0.5237242200 | 0.0028964 | 0.292429882 | PASS |

The initial high center and the falling branch are excluded from the first-rebound apex measurement. Only the first rebound is used; later bounces within the 1 s run do not replace it. Native output frames every 1 ms are retained for viewing; the quantitative apex analysis uses every-step CSV rather than that coarser binary cadence.

An independent gravity-inclusive continuous contact solution gives supplemental e_height=0.524032804605, first apex z=0.292745120795 m at t=0.678989309193 s. It solves the forced damped contact oscillator followed by ballistic rebound, without fitting native results. The official unforced coefficient remains the published comparator. Halving dt also increases the error to this gravity-inclusive solution; this two-step check does not establish monotonic convergence.

## Native force, damping and impulse checks

Stored velocities follow the final half-kick. The force-law audit accounts for gravity when recovering the force-assembly velocity:

    vn_assembly = vz_saved - (dt/2)(Fn/m - 10)
    cn_inferred = (kn delta - Fn)/vn_assembly

Near-zero denominators are excluded from coefficient division, while direct force reconstruction covers valid contact states. The coefficient is inferred from native states; it is not a separate native output channel. Linear spring energy is counted once, and gravitational potential energy is included in mechanical-energy diagnostics. The dashpot intentionally dissipates energy, so energy conservation is not the acceptance criterion.

| Run | Inferred cn (N s/m) | Min. signed Fn (N) | Max. force residual (N) | Max. impulse residual (kg m/s) |
|---|---:|---:|---:|---:|
| CPU dt | 37.249654432 | -85.27055 | 0 | 6e-13 |
| GPU dt | 37.249654432 | -85.27055 | 5.7e-14 | 6e-13 |
| CPU dt/2 | 37.249654432 | -86.10769 | 0 | 2.6e-12 |
| GPU dt/2 | 37.249654432 | -86.10769 | 5.7e-14 | 2.6e-12 |

CPU/GPU whole-history maximum position / velocity / signed-force differences are 3.29e-15 m / 1.34e-12 m/s / 1.69e-10 N.
Halving dt changes apex-derived restitution by -0.0002932196. The absolute differences from the published analytic coefficient are 0.00260316 / 0.002896379. The half-step error increases; this comparison does not demonstrate monotonic convergence. Both time steps use the same physical parameters and 1 s endpoint.
The fixed wall supplies external impulse, while gravity acts throughout. The maximum residual in m(v-v0) = integral(F_contact + mg)dt is 2.62e-12 kg m/s. Sphere momentum alone is not conserved.

Numerical criteria and measurements are retained below and in the companion CSV tables.

Analysis JSON (retained locally; not included) | Independent physics audit (retained locally; not included)

## Complete retained results and verified scope

All four new runs use real **FunDEM SimulationSession** linked to **FunDEMSoftware 0.4.3 CPU/CUDA native libraries**. Installed GUI executable byte identity is not established. Each run retains 1001 configured native output frames, all numerical particle/contact CSV records, input projects and backend logs: **4004 native frames total**.

- **CPU dt**: loadable project (retained locally; not included), complete native result (retained locally; not included); 100,000 numerical steps / 1,001 initial-inclusive native frames.
- **GPU dt**: loadable project (retained locally; not included), complete native result (retained locally; not included); 100,000 numerical steps / 1,001 initial-inclusive native frames.
- **CPU dt/2**: loadable project (retained locally; not included), complete native result (retained locally; not included); 200,000 numerical steps / 1,001 initial-inclusive native frames.
- **GPU dt/2**: loadable project (retained locally; not included), complete native result (retained locally; not included); 200,000 numerical steps / 1,001 initial-inclusive native frames.

Run manifest (retained locally; not included) | Complete native ZIP (retained locally; not included) | Frame/hash retention audit (retained locally; not included)

**Verified scope:** one published full-normal-dashpot linear sphere-wall subcase, beta_n-to-restitution-to-automatic-damping conversion, signed unloading, gravity-aware wall impulse and CPU/CUDA/time-step comparison. Full PFC parameter sweep, no-tension Md=1, active friction/rotation, and many-particle dynamics require separate cases. The existing numbered LAMMPS series and all its earlier pass/fail results remain separate and unchanged.

## Repository evidence and native rerun input

- [Portable native-core input](pfc-linear-restitution-md0-beta02.cpp). The category README documents the runner and model-mapping limits.
- [Compact historical numerical summary](HISTORICAL-RESULTS.md).
- [Source pins, original outcome and retained-artifact hashes](PROVENANCE.md).
- [Compact checkpoint table: first_rebound_comparison.csv](checkpoints/first_rebound_comparison.csv).

A new runner result must be compared independently with these targets; opening or completing the input is not an automatic PASS. Full native frames, per-step histories, execution logs, vendor inputs, proprietary programs and PDF reports are not included.
