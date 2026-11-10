# ON-DEM comparisons

| Case | Implementation | Comparison status |
| --- | --- | --- |
| [Bumpy-wall shear flow](shearFlow/README.md) | Isolated periodic CPU adapter; target `ondemShearFlow` | Transient discrepancies retained; full benchmark not run |

This category contains external benchmark adaptations, not analytical unit
checks. The source path is `ONDEMComparison/shearFlow`; existing executable
names (`ondemShearFlow`, `checkOndemShearFlow`) are unchanged.

The case README separates published parameters, diagnostic coefficients,
recorded checks and unresolved differences. Full calculation requires an
explicitly obtained initial configuration and an explicit long-run command.
Building the source or running the three small checks neither downloads
third-party data nor starts the 200-million-step benchmark. Original particle
coordinates are not bundled because no redistribution license is supplied.
