# plate2D_timeVaryingBC_explicitBoundary

**Solver:** `electroPotentialMultiRegionFoam`  
**Case:** 2D two-region plate with time-varying Dirichlet boundary and explicit interface coupling 
**Author:** Rention Pasolari  
**Toolkit:** SoPLASMA  
**License:** GPLv3  
**Date:** November 2025 

> **Note:** This tutorial requires **Gmsh** to generate the multi-region mesh.  
> If Gmsh is not available, a custom `blockMeshDict` may be created instead to
> construct the two-region geometry manually.
---

## Description

This tutorial solves the electrostatic Poisson equation in the **left fluid region** and
the Laplace equation in the **right dielectric region** of a 2D plate. The two regions
are coupled **explicitly**: each region solves its own equation independently, and the
interface potential is passed staggeredly from one region to the other. An internal
sub-iteration loop is performed at every time step to ensure that both regions reach a
consistent interface potential before proceeding.

A linearly rising voltage (ramp) is applied on the left boundary of the fluid region,
while the right boundary of the dielectric region is grounded. All other boundaries use
`zeroGradient`.

Each time step represents a fully converged steady electrostatic solution, enabling the
simulation of the system response under a continuously increasing applied voltage.

---

## Problem Setup

- **Geometry:** 2D two-region plate  
  - Left: fluid region (solves Poisson equation)
  - Right: dielectric region (solves Laplace equation)

- **Materials:**
  - Fluid region: εᵣ defined in `constant/leftFluid`
  - Dielectric region: εᵣ defined in `constant/rightDielectric`

- **Coupling:**
  - Explicit interface coupling
  - Interface potential passed staggeredly between regions
  - Internal sub-iterations each time step until interface convergence  
    (controlled via `system/<region>/fvSolution` files)

- **Boundary Conditions:**
  - Left boundary of fluid region: linearly rising Dirichlet voltage (ramp)
  - Right boundary of dielectric region: fixedValue 0 V (ground)
  - All other boundaries: `zeroGradient`

- **Charge density:** zero (no space charge in fluid region)

- **Surface charge:** zero (no surface charge at the interface)

- **Time stepping:** each time step corresponds to a fully converged steady electrostatic state

---

## Run Instructions

This case can be executed:

### Serial
```bash
./Allrun-serial
```

### Parallel
```bash
./Allrun-prallel
```
