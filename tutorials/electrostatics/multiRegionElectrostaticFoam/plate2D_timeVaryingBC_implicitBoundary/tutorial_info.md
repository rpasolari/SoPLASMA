# plate2D_timeVaryingBC_implicitBoundary

**Solver:** `multiRegionElectrostaticFoam`  
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

This tutorial solves the electrostatic Poisson equation in the **left gas region** and
the Laplace equation in the **right dielectric region** of a 2D plate. The two regions
are coupled **implicitly**.
A linearly rising voltage (ramp) is applied on the left boundary of the gas region,
while the right boundary of the dielectric region is grounded. All other boundaries use
`zeroGradient`.

Each time step represents a fully converged steady electrostatic solution, enabling the
simulation of the system response under a continuously increasing applied voltage.

---

## Problem Setup

- **Geometry:** 2D two-region plate  
  - Left: gas region (solves Poisson equation)
  - Right: dielectric region (solves Laplace equation)

- **Materials:**
  - Gas region: εᵣ defined in `constant/leftGas`
  - Dielectric region: εᵣ defined in `constant/rightDielectric`

- **Coupling:**
  - Implicit interface coupling

- **Boundary Conditions:**
  - Left boundary of gas region: linearly rising Dirichlet voltage (ramp)
  - Right boundary of dielectric region: fixedValue 0 V (ground)
  - All other boundaries: `zeroGradient`

- **Charge density:** zero (no space charge in gas region)

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

---

### Solver Options

The case can run using:

- OpenFOAM linear solvers (default)
- PETSc linear solvers (if PETSc support is compiled)

If PETSc is installed and you want to use it:

1. Use the PETSc fvSolution file (can be set in `Allrun-*` scripts)  
   `cp system/fvSolution-petsc system/fvSolution`

2. Make sure `controlDict` includes the PETSc library:  
   `libs ("libpetscFoam.so");`

