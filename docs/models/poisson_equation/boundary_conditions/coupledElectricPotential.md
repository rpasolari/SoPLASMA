# Coupled Electric Potential (`coupledElectricPotential`)

The `coupledElectricPotential` is a specialized multi-region boundary condition designed to solve the Poisson equation across interfaces with differing dielectric permittivities ($\epsilon_r$). It is primarily used for plasma–dielectric or dielectric–dielectric interfaces, where it enforces potential continuity and accounts for the physical effects of surface charge accumulation ($\sigma_s$).

Based on OpenFOAM's `compressible::turbulentTemperatureRadCoupledMixed` logic, this BC ensures physical consistency of the electrostatic field across multi-region meshes.

---

## 1. Theory

At any interface between two regions (Region 1 and Region 2), the solver enforces two fundamental physical constraints derived from Maxwell's equations:

### A. Potential Continuity
The electric potential must be continuous across the boundary to prevent non-physical infinite electric fields:
$$\phi_1 = \phi_2$$

### B. Electric Flux Jump (Gauss's Law)
The displacement field ($\mathbf{D} = \epsilon \mathbf{E}$) must account for any surface charge density ($\sigma_s$) accumulated at the interface (e.g., electron or ion accumulation on a dielectric wall):
$$(\mathbf{D}_1 - \mathbf{D}_2) \cdot \mathbf{n}_{1 \to 2} = \sigma_s$$



Substituting the potential ($\mathbf{E} = -\nabla \phi$) and absolute permittivity ($\epsilon = \epsilon_0 \epsilon_r$), the jump condition becomes:
$$-\epsilon_1 \nabla \phi_1 \cdot \mathbf{n}_1 - \epsilon_2 \nabla \phi_2 \cdot \mathbf{n}_2 = \sigma_s$$

Where:
* $\epsilon_0$: Vacuum permittivity ($\approx 8.854 \times 10^{-12}$ F/m).
* $\epsilon_r$: Relative permittivity (dielectric constant) defined in each region.
* $\sigma_s$: Surface charge density ($C/m^2$).

---

## 2. Numerical Formulation

OpenFOAM implements this coupling using a `mixed` boundary condition. The patch value $\phi_f$ is calculated as a weighted average between a Dirichlet-like term ($\text{refValue}$) and a Neumann-like term ($\text{refGrad}$):

$$\phi_f = f \cdot \text{refValue} + (1 - f) \cdot \left( \phi_c + \frac{\text{refGrad}}{\Delta} \right)$$

Where:
* $f$ is the `valueFraction`.
* $\phi_c$ is the internal cell-center value.
* $\Delta$ is the inverse of the distance between cell-center and patch-center (`patch().deltaCoeffs()`).

### Derivation: From Physics to OpenFOAM Coefficients

To map the flux jump condition into the finite volume framework, we discretize the normal gradients. The flux balance at the interface is discretized as:
$$-\epsilon_{loc} \Delta_{loc} (\phi_f - \phi_{loc}) - \epsilon_{nbr} \Delta_{nbr} (\phi_f - \phi_{nbr}) = \sigma_s$$

#### a. Solving for the Patch Potential ($\phi_f$)
Rearranging the terms to isolate $\phi_f$:
$$\phi_f (\epsilon_{loc} \Delta_{loc} + \epsilon_{nbr} \Delta_{nbr}) = \epsilon_{loc} \Delta_{loc} \phi_{loc} + \epsilon_{nbr} \Delta_{nbr} \phi_{nbr} + \sigma_s$$

Which gives:
$$\phi_f = \frac{\epsilon_{loc} \Delta_{loc} \phi_{loc} + \epsilon_{nbr} \Delta_{nbr} \phi_{nbr} + \sigma_s}{\epsilon_{loc} \Delta_{loc} + \epsilon_{nbr} \Delta_{nbr}}$$

#### b. Identifying `mixed` Coefficients

To find the specific expressions for $f$ (`valueFraction`) and $g$ (`refGrad`), we compare our derived physics equation (1) with the standard OpenFOAM `mixed` BC equation (2).

**Equation (1) - Discretized Physics:**
$$\phi_f = \frac{\epsilon_{loc} \Delta_{loc} \phi_{loc} + \epsilon_{nbr} \Delta_{nbr} \phi_{nbr} + \sigma_s}{\epsilon_{loc} \Delta_{loc} + \epsilon_{nbr} \Delta_{nbr}}$$

**Equation (2) - OpenFOAM `mixed` BC:**
$$\phi_f = f \cdot \phi_{nbr} + (1 - f) \cdot \left( \phi_{loc} + \frac{g}{\Delta_{loc}} \right)$$

#### c. Solving for `valueFraction` ($f$)
We rewrite Equation (1) by splitting the numerator to isolate the $\phi_{nbr}$ term:
$$\phi_f = \left( \frac{\epsilon_{nbr} \Delta_{nbr}}{\epsilon_{loc} \Delta_{loc} + \epsilon_{nbr} \Delta_{nbr}} \right) \phi_{nbr} + \frac{\epsilon_{loc} \Delta_{loc} \phi_{loc} + \sigma_s}{\epsilon_{loc} \Delta_{loc} + \epsilon_{nbr} \Delta_{nbr}}$$

By directly comparing the coefficient of $\phi_{nbr}$ with Equation (2), we identify:
$$f = \frac{\epsilon_{nbr} \Delta_{nbr}}{\epsilon_{loc} \Delta_{loc} + \epsilon_{nbr} \Delta_{nbr}}$$

#### d. Solving for `refGrad` ($g$)
Now we look at the second term of Equation (2). Since $1 - f$ must be:
$$1 - f = 1 - \frac{\epsilon_{nbr} \Delta_{nbr}}{\epsilon_{loc} \Delta_{loc} + \epsilon_{nbr} \Delta_{nbr}} = \frac{\epsilon_{loc} \Delta_{loc}}{\epsilon_{loc} \Delta_{loc} + \epsilon_{nbr} \Delta_{nbr}}$$

We can substitute this back into the second half of Equation (1):
$$\frac{\epsilon_{loc} \Delta_{loc} \phi_{loc} + \sigma_s}{\epsilon_{loc} \Delta_{loc} + \epsilon_{nbr} \Delta_{nbr}} = (1 - f) \left( \phi_{loc} + \frac{\sigma_s}{\epsilon_{loc} \Delta_{loc}} \right)$$

Comparing this result to the $(1 - f) \left( \phi_{loc} + \frac{g}{\Delta_{loc}} \right)$ term in Equation (2), it is clear that:
$$\frac{g}{\Delta_{loc}} = \frac{\sigma_s}{\epsilon_{loc} \Delta_{loc}}$$

Resulting in:
$$g = \text{refGrad} = \frac{\sigma_s}{\epsilon_{loc}}$$

---

## 3. Usage

Define the boundary condition in the `0/ePotential` (or your potential field) file for the interface patch:

```cpp
// 0/<regionName>/ePotential
<patchName>
{
    type            coupledElectricPotential;

    // --- Coupling Settings ---
    phiNbr          ePotential;    // Name of potential field in neighbor region
    useImplicit     true;          // Enable monolithic matrix coupling
    
    // --- Surface Charge ---
    // References a volScalarField in the 0/ folder
    surfCharge      surfCharge;    // Local side charge density [C/m^2]
    surfChargeNbr   none;          // Neighbor side charge density

    // --- Monitoring & Debugging ---
    verbose         true;          // Print stats to console
    logInterval     10;            // Write to file every N steps

    // --- Standard Fields ---
    value           uniform 0;     // Initial/restart value
}
```

## 4. Requirements

To ensure the physical and numerical consistency of the electric potential coupling, the following requirements must be met:

### A. Patch Type and Topology
The boundary condition relies on mapping fields between separate regions. Therefore, the interface patches must be of a geometric type that supports topology mapping or topological connection:

* **`mappedWall`**
* **`cyclic` / `AMI`**

### B. Property Dictionary (`electricProperties`)
Each region involved in the coupling (e.g., the plasma domain and the dielectric solid domain) must have its own `electricProperties` dictionary located in its `constant/` directory. 

The BC looks for the `dielectricConstant` keyword to determine the relative permittivity ($\epsilon_r$) of that specific material:

```cpp
// Example: constant/<regionName>/electricProperties for a dielectric region
dielectricConstant 4.0; // Relative permittivity (dimensionless)
```

### C. Surface Charge Management

When the `surfCharge` and/or `surfChargeNbr` keywords are specified, the solver expects the corresponding fields to be present in the `0/` directory of their respective regions. These are standard `volScalarField` files and must be initialized with the correct dimensions ($[0, -2, 0, 0, 0, 1, 0]$ for $C/m^2$).

#### Important Note: Avoiding Surface Charge Double-Counting

The `coupledElectricPotential` boundary condition is designed to be highly flexible, allowing both adjacent regions to contribute to the interface's net surface charge. This requires careful configuration to avoid numerical double-counting.

* **One-Sided Assignment (Recommended):** In most plasma–dielectric simulations, surface charge accumulation is calculated only within the plasma region. 
    * **Plasma Region Setup:** Set `surfCharge` to the name of your charge field (e.g., `surfCharge`) and `surfChargeNbr` to `none`.
    * **Dielectric Region Setup:** Mirror this by setting the local `surfCharge` to `none` and pointing `surfChargeNbr` to the field in the plasma region.

* **Net Charge Logic:** The b.c. internally sums the contributions. If you specify a non-zero field for both `surfCharge` and `surfChargeNbr` on the same interface, the solver will calculate the total density as:
  $$\sigma_{total} = \sigma_{loc} + \sigma_{nbr}$$
  If both fields contain the same physical data, this will lead to a non-physical doubling of the potential jump.


### D. Implicit - Monolithic Coupling

Most solvers developed in the current toolkit allow for an **implicit coupling** of different regions. In a standard segregated approach, each region's Poisson or Laplace equation is solved separately, and the solver iterates to converge the interface value which is computationally expensive.

By setting the `useImplicit` keyword to `true` in both regions, the toolkit enables **monolithic coupling**. This constructs a single, large global matrix containing the equations for both regions, solving the entire system in one step.

#### Important Notes on Monolithic Coupling

a. **Parallel Decomposition Constraints:**
   For implicit coupling to work in parallel, the coupled patch-pairs must reside on the same processor (this applies to the neighboring cells across the interface). Additionally, they must have a matching number of faces. To ensure this during the decomposition process, you must use a constraint in your `decomposeParDict`:

```cpp
// system/decomposeParDict
constraints
    {
    faces
    {
        type    preserveFaceZones;
        zones   (interfaceZone);
        enabled true;
    }
}
```

b. **PETSc Solver Configuration:**
    If the `PETSc` library is used to solve the electric potential, you must explicitly specify `implicitRegionCoupling true` in the `fvSolution` dictionary:

``` cpp
// system/fvSolution
ePotential
{
    solver                  petsc;
    implicitRegionCoupling  true;

    petsc
    {
        options
        {
            ksp_type  cg;
            pc_type   hypre;
            mat_type  aij;
        }
    }
    tolerance   1e-8;
    relTol      0.0;
}
```

c. **Solver Compatibility:**
    If using native OpenFOAM solvers, the `GAMG` (Geometric-Algebraic Multi-Grid) solver will not work with monolithic coupling. This is because the geometric grid agglomeration logic cannot handle the discontinuous connectivity of the coupled region matrices.

d. **Segregated Mode (useImplicit False):**
    If `useImplicit` is set to `false`, a non-monolithic (segregated) coupling is used. In this case, you must specify the convergence controls for the interface iterations in `fvSolution`:

```cpp
// system/fvSolution
ePotentialControls
{
    nonCoupledResidualControl
    {
        tolerance       1e-12;
        maxIter         1000;
    }
}
```

---

