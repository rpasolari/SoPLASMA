# plate2D_timeVaryingBC

**Solver:** `singleRegionElectrostaticFoam`  
**Case:** 2D plate with time-varying Dirichlet boundary condition    
**Author:** Rention Pasolari  
**Toolkit:** SoPLASMA  
**License:** GPLv3  
**Date:** November 2025 

---

## Description

This tutorial solves the electrostatic Poisson equation in a 2D plate. 
A linearly rising voltage (ramp) is applied on the top boundary, while the other boundaries
are grounded.  

Each time step represents a fully converged steady electrostatic solution, allowing simulation of 
an increasing applied voltage.

---

## Problem Setup

- Geometry: 2D plate
- Material: vacuum (ε₀, relative permittivity εᵣ = 1)
- Potential:
  - Top: ramp Dirichlet (0 → Vmax over 10 s)
  - Other walls: fixedValue 0 V (ground)
- Charge density: 0 (no space charge)

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
