# Research only: official PFC Hertz prescribed normal-loading stage

> Historical diagnostic/coverage evidence from the 28–29 October 2026 campaign. Capability statements concern that configured source snapshot and application path, not an exhaustive current-core feature list. Migration changes links and artifact labels, not the recorded outcome.

The retrieved PFC7 page and its three embedded datafiles are retained as research provenance. No FunDEM validation case, execution result, pass/fail conclusion, or report was created for this candidate.

The exact official first stage has two R=0.5 spheres at y=0 and y=1, density2500, G=1e9, nu=0.3, mu=0.5. All degrees of freedom are fixed, but the upper sphere has prescribed vy=-0.01 for0.1 time unit, dt0.001. The resulting normal overlap0.001 gives Hertz force30116.9300968417; the official page prints3.0117e4. The conservative normal law alone maps through E=2G(1+nu)=2.6e9 without fitting.

Native feasibility review on2026-10-29 found that the available FunDEM SimulationSession cannot express both official constrained trajectories while computing their contact:

- FunDEMSoftware/src/model/projectModel.cpp requires a prescribed-motion target to use a Wall material.
- FunDEMSoftware/src/simulation/native/nativeEventSolver.cpp rejects prescribed particles with nonzero inverse mass.
- FunDEMBeta/src/execution/cuda/contactDetectionKernel.cu, detectSphereContact, skips a sphere pair if both inverse masses are zero.

Therefore simply assigning both spheres prescribed motion removes their contact from native detection. Leaving one body dynamically free, injecting compensating forces, manually overwriting poses, or editing the solver would change the official case or validation method. Those substitutions are not used.

The complete official loading path also contains rolling, axial twisting, shear-history scaling and Md3 damping, so it must not be claimed as equivalent from normal-law agreement alone. This candidate is excluded from executed-case counts and manifests.
