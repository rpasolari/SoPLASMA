# positiveStreamer_AMR

**Solver:** `soPlasmaFoam`  
**Case:** 2D axisymmetric positive streamer in dry air with adaptive mesh refinement  
**Author:** Rention Pasolari  
**Toolkit:** SoPLASMA  
**License:** GPLv3  
**Date:** June 2026

---

## Description

This tutorial simulates a positive streamer discharge in dry air using adaptive mesh
refinement (AMR) via the blastAMR library. The streamer propagates from anode to
cathode in a 2D axisymmetric domain, driven by a background electric field enhanced
by an initial Gaussian ion seed.

The AMR dynamically refines the mesh at the streamer head and channel edges based on
the Townsend ionization criterion, and coarsens it in the quiescent bulk regions,
keeping the total cell count low while maintaining high resolution where needed.

---

## Problem Setup

- **Geometry:** 2D axisymmetric domain, 1.25 × 1.25 cm²

- **Species:** electrons and positive ions (ions assumed immobile)

- **Background field:** −15 kV/cm (below breakdown threshold)

- **Initial condition:** Gaussian ion seed at z₀ = 1 cm

- **Boundary Conditions:**
  - Top electrode: fixed voltage 18.75 kV
  - Bottom electrode: ground (0 V)
  - Axis and outer boundary: `zeroGradient`

- **AMR settings:**
  - Refinement criterion: α·Δx > 0.5
  - Unrefinement criterion: α·Δx < 0.1
  - Maximum refinement levels: 3 (base 24 µm → minimum 3 µm)

- **Time stepping:** fixed 1 ps time step

---

## Run Instructions

### Serial
```bash
./Allrun-serial
```

### Parallel
```bash
./Allrun-parallel
```
