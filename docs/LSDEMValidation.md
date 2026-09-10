# LS-DEM validation benchmark

## Selected case

The recommended first validation case for FunDEMBeta is the
**two-sphere central-compression benchmark** from Feldfogel, Karapiperis,
Andrade, and Kammer:

- [Open-access journal article](https://doi.org/10.1002/nme.7400)
- [Open reference dataset](https://doi.org/10.3929/ethz-b-000640908)
- [Preprint](https://arxiv.org/abs/2209.15431)

This is a better first code-verification case than a full triaxial specimen.
It isolates LS contact quadrature, has one prescribed displacement and one
measured reaction, and directly tests whether surface refinement changes the
physical response. The paper then extends the same test to a 27-sphere
assembly and interlocking blocks.

The foundational experimental-scale LS-DEM validation remains the XRCT-based
triaxial study by Kawamoto, Ando, Viggiani, and Andrade:

- [Level set discrete element method for three-dimensional computations with
  triaxial case study](https://doi.org/10.1016/j.jmps.2016.02.021)

That study is a valuable second-stage validation, but reproducing it requires
the segmented particle geometries, specimen fabric, boundary conditions, and
calibrated material parameters. It is therefore not a compact regression case.

## Published configuration

| Parameter | Value |
| --- | ---: |
| Shape | Two identical LS spheres |
| Sphere radius, \(R\) | \(5\ \mathrm{mm}\) |
| Initial state | Just touching on the center line |
| Grain 0 | Fixed |
| Grain 1 | Prescribed center-line displacement toward Grain 0 |
| Penetration, \(\delta\) | \(0\) to \(0.1\ \mathrm{mm}\) |
| Distributed normal stiffness, \(k_n^*\) | \(1\ \mathrm{GPa/mm}\) |
| Response | Center-line reaction \(P(\delta)\) |
| Primary study | Five progressively refined surface meshes |

The stiffness is traction per penetration, not the ordinary DEM spring
stiffness. In FunDEM SI units,

\[
1\ \mathrm{GPa/mm}=10^{12}\ \mathrm{N/m^3}.
\]

The adapted formulation evaluated by the paper is

\[
\mathbf F_n=\sum_a k_n^*\,d_a\,\mathbf n_a\,A_a,
\]

where \(A_a\) is the tributary surface area of node \(a\). This is the same
quantity stored by FunDEM in `LSSurfaceNode::area_` and applied to LS nodal
contact stiffness. The benchmark therefore targets the current implementation
directly.

## FunDEM realization

Use `levelset::Sphere{0.005}`, zero gravity, zero tangential stiffness, and
zero friction. Build the same analytic sphere at several independent
resolutions:

| Series | Surface subdivision levels | Surface-node counts |
| --- | --- | --- |
| Surface convergence | 2, 3, 4, 5 | 162, 642, 2,562, 10,242 |
| Volume-grid convergence | Hold level 5 | 20, 40, 80 nodes per diameter |

The node counts follow the icosphere relation used by
`levelset::Sphere::buildSurfaceNode`. They do not need to equal the paper's
mesh counts; convergence is tested by the limiting response, not by identical
triangulations.

For each resolution:

1. Fix Grain 0 and give Grain 1 the deformable LS material with
   `normalStiffnessPerUnitArea = 1.0e12`.
2. Move Grain 1 from \(\delta=0\) to \(0.1\ \mathrm{mm}\) in small,
   uniform displacement increments.
3. At every increment, clear stale force and contact history, perform contact
   detection and one elastic force evaluation, and record the center-line
   reaction. Do not advance Grain 1 dynamically.
4. Write `penetration`, `reactionForce`, surface-node count, volume-grid
   spacing, and backend to a DAT or CSV file.

The moving particle must remain finite-mass while its pose is prescribed:
FunDEM intentionally gives two infinite-mass bodies zero effective stiffness.
A validation fixture should therefore impose Grain 1's position kinematically
rather than mark both particles as fixed.

For small penetration, one-sided surface quadrature has the leading-order
continuum check

\[
P\simeq \frac{\pi}{2}R k_n^*\delta^2.
\]

At \(R=5\ \mathrm{mm}\) and \(\delta=0.1\ \mathrm{mm}\), this is about
\(78.5\ \mathrm N\). The paper curve and its raw dataset remain the primary
reference because the exact discrete response also contains finite-curvature
and grid-interpolation effects.

## Acceptance checks

The future executable should fail clearly unless all of these hold:

- the two finest surface-discretization curves differ by less than a declared
  tolerance, initially 2% over \(0.02\leq\delta\leq0.1\ \mathrm{mm}\);
- refining only the level-set volume grid also approaches a stable curve;
- reaction is proportional to \(k_n^*\) for the paper's
  \(0.25, 0.50, 0.75, 1.00\ \mathrm{GPa/mm}\) stiffness series;
- CPU and GPU curves agree within floating-point reduction tolerance;
- the finest FunDEM curve is compared point-by-point with the paper's open
  reference data, reporting RMSE, peak error, and the maximum relative error
  away from zero load.

The 2% threshold is a proposed FunDEM regression criterion, not a value claimed
by the paper. It should be finalized after the first resolution study.
