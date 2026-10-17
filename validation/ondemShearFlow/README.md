# ON-DEM granular shear flow between bumpy walls

This is a **validation-only CPU adaptation**, not a claim of completed ON-DEM
validation or a project that the current Workbench can run. It reuses FunDEM
spheres, their stable container indices, velocity-Verlet integration helpers,
solver clock, and binary VTU writer. Periodic search and the benchmark pair law
are confined to the external-force component in this directory. Default `SphereDEM`, GPU kernels, and existing
Workbench project semantics are unchanged.

## External contact-force component

`periodicContactForce.h/.cpp` implements the benchmark contact model separately
from the time integrator and material definitions. Its `apply(host, historyDt)`
adds normal elastic/damping force, tangential elastic/damping force, and the
associated contact moments to existing host particle loads. It does not clear
unrelated forces, move particles, or advance the solver clock. The explicit
`kn`, `kt`, `cn`, `ct`, and friction coefficient belong to this component; no
restitution-to-dashpot or material-stiffness combination is performed.

The case adapter overrides the existing
`addSphereExternalForceAndTorque(host)` hook and invokes the component after
clearing loads once. Initialization uses zero history time; each DEM advance
uses the full current history step once. The existing periodic velocity-Verlet
sequence is retained. It does not call the native nonperiodic contact search
and then run another periodic search: only the external component evaluates
contacts. This avoids both duplicate force assembly and duplicate searches.

The sphere material now supplies **density and mass/inertia only**, with zero
native contact stiffness/friction. Setting native material stiffness to zero
does not disable the external component. Periodic wrapping, stable pair
history, fixed-wall roles, source IDs, and the scientific contact output remain
unchanged. This is an in-process C++ case component, not a DLL module or a
new Workbench model selector. It implements the same declared benchmark law;
the refactor does not resolve the source parameter inconsistency by changing
the force coefficients.

The external-force extraction passed all three ON-DEM checks (periodic
mechanics, Mercury IO, and standalone additive force/hook behavior). A separate
1000-step original-particle run produced byte-identical particle and contact
files at steps 0, 500, and 1000 before and after extraction. This verifies the
refactor, not agreement with the author's reference results.

## Sources and frozen input

- [Benchmark definition](https://on-dem.atlassian.net/wiki/spaces/Index/pages/244744204/Granular+shear+flows+between+bumpy+walls+case).
- [Author's repository at commit 56b7cfdfd774cf5ae057a1c5c90e5191a0deb97b](https://github.com/dalilavescovi/GranularShearFlows/tree/56b7cfdfd774cf5ae057a1c5c90e5191a0deb97b).
- [Original initial coordinates](https://raw.githubusercontent.com/dalilavescovi/GranularShearFlows/56b7cfdfd774cf5ae057a1c5c90e5191a0deb97b/InitialGeometricalConfiguration/GranularShearFlow.data).
- [Author's initial particle output](https://raw.githubusercontent.com/dalilavescovi/GranularShearFlows/56b7cfdfd774cf5ae057a1c5c90e5191a0deb97b/outputs_example/Mercury_out.data.0000),
  used to confirm the directions of the wall velocities.

The repository does not supply a license or the actual simulation script. The
dataset is therefore fetched from its author rather than redistributed in the
FunDEM source. Fetching is explicit, pinned to that commit, TLS-verified, and
SHA256-checked:

```text
947fec877c67e9ba7e8d1d46b66f47f66df6d713e0b33014556f967ce7ca099d
```

The initial file's one-based source IDs are retained for traceability; runtime
indices and Mercury contact indices are zero-based. No particle sorting,
regeneration, or settling step changes the initial coordinates.

## Published setup

All quantities use the reference normalization `rho_p = d = V = 1`; they are not
calibrated real-world material values.

| Quantity | Value |
| --- | --- |
| Mobile spheres | 3000 |
| Prescribed spheres | 240 bottom + 240 top |
| Diameter; material density | 1; 1 |
| Domain | `[0,20] x [0,10.3923048454133] x [0,22]` |
| Periodic directions | X and Y only |
| Bottom / top center Z | 0.5 / 21.5 |
| Clear wall gap | 20 |
| Bottom / top velocity | `(-1,0,0)` / `(1,0,0)` |
| Gravity | Zero |
| Pair `kn`, `kt` | `2e5`, `5.7143e4` |
| Dashpot inputs `cn`, `ct` | `50`, `14.2857` |
| Sliding friction coefficient | 0.3 |
| Time step | `7.2320e-5` |
| Full benchmark | `200000000` steps |
| Output interval | 50000 steps |

Wall grains retain the same radius/material description, but their velocities
are prescribed and their accelerations are disabled. Wall-wall contacts do not
contribute. Every mobile contact uses the tabulated **pair** coefficients;
native series-combination rules are deliberately not applied a second time.

## Build and check

Keep the terminal in the parent `workplace/` containing `FunDEMBeta/`:

```bash
cmake -S FunDEMBeta -B build -DFUNDEM_ENABLE_CUDA=OFF
cmake --build build --target ondemShearFlow checkOndemShearFlow --parallel
cmake -DDESTINATION=build/validation/ondemShearFlow_input/GranularShearFlow.data -P FunDEMBeta/validation/ondemShearFlow/fetchInitialConfiguration.cmake
build/validation/ondemShearFlow --check
```

This CPU-only configure command is for the independent benchmark; it does not
require changing a separate GPU Workbench build. With a multi-configuration
generator, add `--config Release` to the build and use
`build/validation/Release/ondemShearFlow.exe --input` with the downloaded file.

`checkOndemShearFlow` runs small numerical, external-force assembly, and file-format regression cases.
`--check` loads all original particles and initializes the periodic force state
without advancing time. Neither operation is a steady-flow validation.

A short opt-in run, with optional binary VTU for ParaView:

```bash
build/validation/ondemShearFlow --steps 100 --output-interval 50 --vtu --output build/ondem-short-run
```

The full reference duration must be requested explicitly:

```bash
build/validation/ondemShearFlow --steps 200000000 --output-interval 50000 --output build/ondem-reference-run
```

No steps are run by default. Existing nonempty output directories are refused,
not cleared. `--no-output` disables frame writing for timing measurements;
`--vtu` only adds visualization files, not another force calculation.

For an explicitly requested reference-snapshot diagnostic, the runner also
accepts `--normal-damping`, `--sliding-damping`, and `--sliding-stiffness`.
These override only that run and leave the published defaults unchanged. The
startup log prints the coefficients actually passed to both the force component
and energy accounting. The inferred diagnostic set is not an author-confirmed
replacement for the published parameter table:

```bash
OMP_NUM_THREADS=2 build/validation/ondemShearFlow --steps 100000 --output-interval 50000 --time-step 7.192545032e-5 --normal-damping 15 --sliding-damping 4.285714285714286 --sliding-stiffness 57142.857142857145 --output build/ondem-snapshot-coefficients-run
```

Here `ct=30/7` and `kt=400000/7`; the decimal arguments retain their double
precision. Initial coordinates, particle order, normal stiffness, friction,
wall motion, domain and gravity are unchanged. Keep diagnostic outputs separate
from the earlier page-coefficient run.

The requested diagnostic completed 100000 steps on 2026-10-17. Mobile total
kinetic-energy differences at the two reference times were -0.323035% and
+0.487018%, down from -51.077875% and -64.970898% with the page coefficients.
Individual trajectories and instantaneous contact pairs still differed; this
is not exact snapshot reproduction or a steady-state validation pass. The
diagnostic run reused the external component without changing its force law.

## Periodic contacts and history

The linked-cell search is rebuilt every step. Neighbor cells wrap in X/Y, while
Z is nonperiodic. Each unordered particle pair is evaluated once; minimum-image
separation is used for the normal, both contact lever arms, relative contact
velocity, and torques. Tangential histories use stable original pair indices,
not grid positions or periodic-image indices. Wrapping through a side therefore
does not create a new physical contact. Open contacts lose their old history.

The implementation uses direct-force dashpots:

```text
n  = minimumImage(x_i - x_j) / distance
vn = dot(v_contact_i - v_contact_j, n)
Fn = max(0, kn * overlap - cn * vn)
Ft_trial = -kt * tangentialSpring - ct * tangentialVelocity
|Ft| <= mu * Fn
```

The tangential spring is transported into the current tangent plane, integrated,
and corrected after capping the **total** spring-plus-dashpot force. This explicit
Cundall-style return differs from the default FunDEM law, which limits damping
separately. It is documented here rather than silently changing existing cases.
Rolling and torsional resistance are absent.

## Reference ambiguities and validation status

The page does not provide a force-equation convention for its damping inputs.
Here `cn` and `ct` multiply relative velocities directly, rather than multiplying
them by effective mass. For unit diameter and density, the mobile grain mass is
`pi/6`, not one. With `m_eff = pi/12`, `kn = 2e5`, and `cn = 50`, the analytical
damped collision time divided by 50 is approximately `7.2319727693e-5`, consistent
with the page's rounded time step. This is a documented inference, not confirmation
from the benchmark author's simulation script.

The supplied example frame 0001 has time `3.596272516`; divided by the stated
50000-step interval, it implies `7.192545032e-5`, not the page's time step. This
runner follows the parameter table, not that inconsistent sample timestamp.
Its Coulomb return, history rotation, and noncohesive normal-force truncation
must also be checked against the authors' actual contact implementation before
claiming exact cross-code equivalence.

Checked on 2026-10-17 with GCC 13.3, Release, CPU: both added CTest regressions
passed, including a 120-particle linked-cell/brute-force comparison. The frozen
3480-particle input passed `--check` and a 100-step smoke run to time `0.007232`.
Its initial output matched all source coordinates and radii, and all three
particle-output frames were finite. These are implementation checks only.

The full steady-state run, time-step convergence, and continuum-field comparisons
remain required. Compilation, an initial-state check, or a short smoke run is
not evidence that those comparisons have passed.

## Executed transient comparison

On 2026-10-17, a separate 100000-step CPU run used the published pair
coefficients unchanged (`kn=2e5`, `kt=5.7143e4`, `cn=50`, `ct=14.2857`,
`mu=0.3`). Its time step was explicitly overridden to `7.192545032e-5` so
frames 50000 and 100000 matched the author's example timestamps. This is
an example-time comparison, not a run with the page's default time step.
`--time-step` changes only the integration step; the default remains
`7.2320e-5`. The original positions, radii and translational velocities match
the reference exactly; its negligible initial angular velocities differ by at
most `3.970502315e-12`.

The instantaneous **mobile-particle translational plus rotational kinetic
energy**, excluding all prescribed wall spheres, did not agree:

| Normalized time | Author's example | FunDEM | Relative difference |
| --- | --- | --- | --- |
| 3.596272516 | 5.20357556 | 2.54569974 | -51.08% |
| 7.192545032 | 8.47848212 | 2.96993614 | -64.97% |

The run completed in 291.01 wall-clock seconds on the development host with
`OMP_NUM_THREADS=2`. The full 200-million-step benchmark was **not** run.
Constant-cost extrapolation would be about 6.7 days on this setup; changing
contact density and system load can change the actual cost.

To reproduce the comparison, first obtain the six pinned reference files
described below. Use new, empty output directories:

```bash
OMP_NUM_THREADS=2 build/validation/ondemShearFlow --steps 100000 --output-interval 50000 --time-step 7.192545032e-5 --output build/ondem-example-time-run
python3 -m pip install numpy==2.2.6 matplotlib==3.10.3
python3 FunDEMBeta/validation/ondemShearFlow/compareReference.py --reference build/ondem-reference --simulation build/ondem-example-time-run --output build/ondem-comparison
```

Install the analysis packages in an isolated Python environment when required
by the host's package policy. The Python tools are optional analysis utilities,
not dependencies of the C++ solver. The comparison writes `metrics.csv`,
`profiles.csv`, `summary.json`, and `comparison.png`; frame indices, particle
ordering and timestamps are checked before comparison. Profiles are simple
instantaneous particle-center bins, explicitly **not** MercuryCG fields or
late-time averaged benchmark results. These observed differences are not a
validation pass. The reference force audit below provides strong evidence that
the supplied examples do not use the page's stated normal damping coefficient;
the model has not been tuned to erase that discrepancy.

## Auditing the author's reference snapshots

The pinned repository supplies three complete 3480-particle snapshots, at times
`0`, `3.596272516`, and `7.192545032`. Their corresponding contact files are
also complete for those snapshots: an independent periodic geometric search
finds respectively zero, one, and three physical contacts, exactly accounting
for the zero, two, and five exported contact records. Small contact files must
not be mistaken for truncated examples. These early transient snapshots are
not the late-time averaged continuum profiles required by the benchmark.

After placing the six upstream `Mercury_out.data.0000` through `.0002` and
`Mercury_out.fstat.0000` through `.0002` files in a reference directory, run:

```bash
python3 FunDEMBeta/validation/ondemShearFlow/auditReference.py --reference-dir build/ondem-reference --output build/ondem-reference-audit.json
```

This standard-library-only script verifies every file against its pinned
SHA256, reconstructs all periodic contacts, and records the numerical audit and
its assumptions in JSON. An existing output report is never overwritten.

There is a second source inconsistency beyond the sample timestamps. Assuming
the author's output uses velocity Verlet, the contact force is evaluated before
the final velocity half kick. Recover the force-evaluation velocity by
`v_half = v_output - dt * F / (2*m)` for mobile grains, leaving prescribed wall
velocities unchanged. With all exported forces included, the four unique
sample contacts independently give
`cn = (kn*overlap - Fn) / vn_half = 15.000000` to within `6.5e-8`, including
the mobile-wall contact. The sample step inferred from its timestamps,
`7.192545032e-5`, also agrees with the analytical collision time divided by 50
for `cn=15`, rather than the page's `cn=50`.

This is strong conditional evidence that the example files used a different
normal dashpot from the published parameter table, not confirmation from the
missing author simulation script. The tangential dashpot and history update
cannot be uniquely identified from these scalar elongation samples. The
FunDEM benchmark continues to use the published `cn=50`, `ct=14.2857`, and
`dt=7.2320e-5`; the audit does not tune the model to match the examples.

## Output and comparison

Frames are paired `ShearFlow.data.0000` and `ShearFlow.fstat.0000` files. The
particle file includes mobile and prescribed wall spheres, preserving their
array indices. Each mobile-mobile contact is written twice, reversing indices
and vector directions; a mobile-wall contact is written once with the mobile
particle first. Fixed-fixed records are absent. No dummy force records are added.
The normal points from the second particle to the first; the tangential direction
is the direction of the force acting on the first particle. Tangential elongation
is exported as the negative spring norm, matching the sign in the reference
output examples (the page itself describes this field as a length).

These discrete files are intended for the reference MercuryCG workflow. The
sample repository files are format examples, not complete experimental profiles.
For the benchmark comparison, generate vertical profiles of streamwise velocity,
concentration, granular temperature, pressure, and shear stress. Average over the
periodic horizontal X/Y directions and the prescribed late-time output window,
using the reference coarse-graining definition and 100 horizontal slices. The
page's phrase about averaging X/Z is inconsistent with a vertical profile; X/Y
is the explicit interpretation used here. No unverified MercuryCG command line
or substitute coarse-graining kernel is imposed by this runner.

Optional VTU particle frames include radius, velocity, angular velocity, force,
original source ID, and an infinite-mass flag. In ParaView, use the radius with
a sphere glyph. These are visualization outputs, not the benchmark's continuum
profiles. The default software project serializer has no periodic-boundary
contract, so no misleading runnable `.fundem.json` is generated.
