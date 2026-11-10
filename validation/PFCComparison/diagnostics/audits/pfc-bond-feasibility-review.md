# PFC bond verification: reviewed against the implemented FunDEM bond path

> Historical diagnostic/coverage evidence from the 28–29 October 2026 campaign. Capability statements concern that configured source snapshot and application path, not an exhaustive current-core feature list. Migration changes links and artifact labels, not the recorded outcome.

Date: 29 October 2026. This review corrects the initial overly broad description of FunDEM as unbonded. FunDEM **does support** sphere-sphere, sphere-LS and LS-LS virtual bonds through loadable `eventBonds`, on CPU and CUDA. Current source snapshots and hashes (retained locally; not included).

The configured bond supplies four physical stiffnesses and a positive equivalent length. `bond::setStiffness` maps them to B1=kn, B2=ks*L^2, B4=kt and B3=kb-B2/4-B4/2. `calculateBondForceAndTorque` evaluates endpoint-position/orientation forces, couples and four elastic energies. Irreversible mixed-mode Benzeggagh-Kenane energy damage is available; zero cross-sectional area disables damage evolution. These are real implemented capabilities, not unsupported merely because a PFC example is bonded.

## Wave propagation: axial elasticity is mappable; original pulse loading is unavailable

The [official PFC wave example](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/WavePropagation/WaveProp.html) has fifty collinear spheres in 3D, R=0.5, material density1500 and first-ball density750. Only axial motion is excited. Thus transverse/spin constraints are inactive invariants and do not by themselves prevent a faithful mapping.

Its active elastic `linearcbond` branch can reduce to a FunDEM axial B1 bond with the separate overlap-contact stiffness set to zero, avoiding double-counted compression. Other bond stiffnesses can be zero; damage can be disabled, matching the source's elastic-only solve. Half-mass first-ball assignment is directly supported. The pair stiffness derived from the two source particle stiffnesses in series is 7,853,981.63397 N/m. This axial mapping does not establish general PFC contact-bond shear or fracture equivalence.

The actual blocker is the prescribed external force on that **finite-mass** first ball. In SI units it is

    Z = pi*R^2*rho_bar*C = 78539.8163397 kg/s
    F(t,v1) = 7853.98163397*(1-cos(8*pi*t)) - Z*v1, 0<=t<=0.25 s
    F(t,v1) = -Z*v1, t>0.25 s

The configured particle-damping module can reproduce the quiet term exactly with rate Z/m_first=200 1/s. The force registry has no finite-particle additive half-cosine loading operator for the first term. Prescribing first-ball motion requires Wall/infinite mass and changes its acoustic impedance and reflected-wave solution. Auxiliary particles, force-cancellation mechanisms or injected poses would change the official configuration; none were introduced. No full wave case is counted as executed.

## Cantilever: bonds exist; the official tip forces and moments are missing

The [official cantilever](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/Cantilever/CantileverBeam.html) uses a contact spring and a parallel bond in parallel. Four bond stiffnesses may reproduce an axial or linearized elastic branch, and a fully fixed root sphere can be represented by an infinite-mass body. A blanket exclusion based on absent bonds is incorrect.

Every nonzero official subcase applies force and/or moment to the tip sphere. The current configured force registry provides viscous damping, sphere hydrodynamics and periodic boundaries, without an arbitrary tip-force or tip-moment operator. The source also uses PFC local damping during static convergence. Changing force-controlled loading to a fitted tip displacement would be a different protocol and is not done. The high source strengths keep fracture inactive, so differences between PFC stress-strength breakage and FunDEM BK damage are **not** used as the decisive blocker for this elastic test. General finite-rotation equivalence of the incremental PFC parallel bond and the endpoint-frame virtual bond is not asserted.

## Spring-network capabilities: specific operators differ

The [RBSN verification](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/springnetwork/springnetwork.html) exercises more than generic bonded elasticity. Its [model](https://docs.itascacg.com/pfc700/common/contactmodel/rbsn/doc/manual/cmspringnetwork.html) constructs interface-area/moment stiffness and adds body-stress-dependent fictitious forces for a target Poisson ratio; the verification activates a nonzero ratio. Other stages use displacement-based tensile/cohesive tables, frictional weakening and healing.

The configured B1..B4/BK bond has no body-stress Poisson correction or those table/healing state operators. They cannot be replaced by relabeling the four stiffnesses or BK fracture energies. This identifies a model-specific gap while acknowledging the implemented FunDEM bonds.

The inventory reasons in [PFC-COVERAGE.md](../../COVERAGE.md) and coverage JSON (retained locally; not included) now reflect these reviewed distinctions. Unsupported entries are neither simulation failures nor completed validations.
