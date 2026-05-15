# plate2D_timeVaryingBC_implicitBoundary

**Solver:** `electroPotentialMultiRegionFoam`  
**Case:** 2D two-region plate with time-varying Dirichlet boundary and implicit interface coupling  
**Author:** Rention Pasolari  
**Toolkit:** SoPLASMA  
**License:** GPLv3  
**Date:** November 2025 

> **Note:** This tutorial requires **Gmsh** to generate the multi-region mesh.  
> If Gmsh is not available, a custom `blockMeshDict` may be created instead to
> construct the two-region geometry manually.
---

## Description

This tutorial solves the electrostatic Poisson equation in the **left fluid region**
and the Laplace equation in the **right dielectric region** of a 2D plate. Unlike the
explicit case, here the coupling between the two regions is **implicit (monolithic)**:
both domains are assembled into a single global matrix and solved simultaneously.

The implicit coupling is activated by setting `useImplicit true` in the coupled patch
boundary condition for the electric potential (`ePotential`). This ensures fully
consistent interface continuity without requiring sub-iterations between regions.

A linearly rising voltage (ramp) is applied on the left boundary of the fluid region,
while the right boundary of the dielectric region is grounded. All other boundaries use
`zeroGradient`.

Each time step represents a fully converged steady electrostatic solution, enabling the
simulation of the system response under a continuously increasing applied voltage.

---

## Problem Setup

- **Geometry:** 2D two-region plate  
  - Left: fluid region (Poisson equation)
  - Right: dielectric region (Laplace equation)

- **Materials:**
  - Fluid region: εᵣ defined in `constant/leftFluid`
  - Dielectric region: εᵣ defined in `constant/rightDielectric`

- **Coupling:**
  - Implicit (monolithic) interface coupling
  - Single global matrix including both regions
  - Enabled by `useImplicit true` in the boundary condition
  - No explicit sub-iterations required

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
> **Parallel execution note:**  
> For implicit (monolithic) multi-region coupling, the interface must be decomposed
> consistently so that matching portions of the interface lie on the same processor.  
> This is enforced using the `constraint` entry in `system/decomposeParDict`.

### Solver Options

The case can run using:

- OpenFOAM linear solvers (default)
- PETSc linear solvers (if PETSc support is compiled)

If PETSc is installed and you want to use it:

1. Use the PETSc fvSolution file:  
   `cp system/fvSolution-petsc system/fvSolution`

2. Make sure `controlDict` includes the PETSc library:  
   `libs ("libpetscFoam.so");`
