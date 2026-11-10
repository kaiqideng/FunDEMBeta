# PFC Hertz prescribed normal stage: excluded candidate

> Historical diagnostic/coverage evidence from the 28–29 October 2026 campaign. Capability statements concern that configured source snapshot and application path, not an exhaustive current-core feature list. Migration changes links and artifact labels, not the recorded outcome.

Research sources and the detailed feasibility note are retained in [reference/pfc/hertz-normal-source](../reference/pfc/hertz-normal-source/FEASIBILITY.md). No native case was executed or counted.

FunDEMSoftware `projectModel.cpp` requires prescribed-motion targets to use Wall material; `nativeEventSolver.cpp` requires their inverse mass to be zero. FunDEMBeta `contactDetectionKernel.cu` skips sphere contacts when both inverse masses are zero. The official two constrained/prescribed sphere trajectories therefore cannot produce their contact through the current configured SimulationSession path. No solver edit, force cancellation or pose injection is used to bypass this restriction.

The next executed official subcase instead uses the original freely falling PFC linear-contact ball12, Md0 and beta_n0.6.
