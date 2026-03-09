# Drift-Diffusion Approximation

The **Drift-Diffusion Approximation** is a simplified model used to describe the transport of charged particles (like electrons and ions) within a medium, such as a plasma. It is primarily used when the system is in a "high-collisional" regime, where the particle's momentum relaxes much faster than the time scales of interest.

---

## 1. Theory

The core assumption of this approximation is the **neglect of inertia**. In this regime, we assume that particles instantly reach their terminal velocity because the forces acting on them (from electric fields or concentration gradients) are perfectly balanced by frequent collisions with the background medium.

### A. The Flux Equation
The total particle flux $\mathbf{\Gamma}$ is expressed as the sum of two distinct physical mechanisms: **Drift** and **Diffusion**.

$$\mathbf{\Gamma} = \mathbf{\Gamma}_{drift} + \mathbf{\Gamma}_{diff}$$

#### a. Drift Component
Drift describes the motion of charged particles in response to an external force, typically an electric field $\mathbf{E}$.
* **Formula:** $\mathbf{\Gamma}_{drift} = \pm n \mu \mathbf{E}$
* **$n$:** Particle density.
* **$\mu$:** Mobility of the particle.
* **Sign:** Positive for cations, negative for electrons and anions.

#### b. Diffusion Component
Diffusion describes the motion of particles from regions of high concentration to regions of low concentration due to random thermal motion (Brownian motion).
* **Formula:** $\mathbf{\Gamma}_{diff} = -D \nabla n$
* **$D$:** Diffusion coefficient.
* **$\nabla n$:** The gradient of the particle density.

### B. Continuity Equation
To find the evolution of the particle density over time, the drift-diffusion flux is plugged into the continuity equation:

$$\frac{\partial n}{\partial t} + \nabla \cdot \mathbf{\Gamma} = S - R$$

Where **$S$** represents source terms (generation) and **$R$** represents sink terms (recombination or loss).

---

## 2. Numerical Implementation & Controls

The drift-diffusion can be selected in `constant\plasmaSpecies`, for each participating particle, or by setting the default choices:

``` cpp
// constant/<regionName>/plasmaSpecies
species (e pIon);

speciesProperties
{

  defaultProperties
  {
    driftDiffusionFluxScheme  driftDiffusion;
    advectionMode             implicit;
  }

  e
  {
    transportModel            driftDiffusion;

    charge                    -1;
    mass                      9.1093837e-31;

    mobilityModel             constant;
    mobilityCoeffs
    {
      mu                      1e-2;
    }

    diffusivityModel          constant;
    diffusivityCoeffs
    {
      D                       1e-6;
    }
  }

  pIon
  {
    transportModel            driftDiffusion;

    charge                    +1;
    mass                      1.67e-26;

    mobilityModel             constant;
    mobilityCoeffs
    {
      mu                      1e-4;
    }

    diffusivityModel          constant;
    diffusivityCoeffs
    {
      D                       1e-6;
    }
  }
}

```

---
