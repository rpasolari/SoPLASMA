# How to use AMR in the foamPlasmaToolkit

In **foamPlasmaToolkit**, the **blastAMR** library has been incorporated for effective Adaptive Mesh Refinement. The `blastAMR` library comes directly with the toolkit in the `ThirdParty` directory. 

> ⚠️ **IMPORTANT:**
> You **must** use the version provided in the `ThirdParty` directory, as specific modifications have been applied to ensure compatibility with the toolkit.

`blastAMR` supports refinement in 2D, 3D, and 2D-axisymmetric meshes. Unlike the built-in AMR tool in OpenFOAM, it is not limited to hexahedral cells and can handle complex geometries.

Currently, `blastAMR` is incorporated only into the **plasmaDielectricFoam** solver and is not available for solvers that solve electrostatics only.

---

## The Basic Files
Implementing AMR with `blastAMR` requires editing several files. In multi-region cases, these must be handled for each specific region:

* **dynamicMeshDict**: Located in `constant/$region/dynamicMeshDict`.
* **plasmaSimulationControls**: A global file in `system/plasmaSimulationControls`.
* **decomposeParDict**: Located in each region’s system directory (must be copied to all processor directories).
* **refineField**: Located in the `0/` folder of each region participating in the AMR process.

---

## The dynamicMeshDict File
This file controls the core AMR options. A typical configuration for `foamPlasmaToolkit` using a coded error estimator based on electron density ($n_e$) is shown below:

```cpp
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      dynamicMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dynamicFvMesh   adaptiveFvMesh;
refiner         polyRefiner;

errorEstimator  coded;
name            electronDensity;

code
#{
    // Initialize error field to 0
    error_ = 0.0;

    if (mesh_.foundObject<volScalarField>("n_e"))
    {
        const volScalarField& ne = mesh_.lookupObject<volScalarField>("n_e");
        
        const scalar threshold = 1e16;

        forAll(error_, i)
        {
            // Binary refinement: If above threshold, set to 1.0
            if (ne[i] > threshold)
            {
                error_[i] = 1.0;
            }
            else
            {
                error_[i] = 0.0;
            }
        }
    }
#};

balance             yes;

lowerRefineLevel    0.5; 
upperRefineLevel    1.0;
unrefineLevel       0.1;

refineInterval      20;   
unrefineInterval    20;
balanceInterval     20;

maxRefinement       4;   
maxCellsLimit       1000000;

nBufferLayers       1; 

correctFluxes
( 
    (phi none)
);
```
For more information, refer to the [blastAMR tutorials](https://github.com/STFS-TUDa/blastAMR).

---

## The plasmaSimulationControls
This file controls global AMR behavior, specifically for multi-region coordination.

```cpp
AMR
{
    enable      true;
    AMRControls
    {
        restrictAMRToPlasma false;
        conformalAMRSync    true;
    }
}
```

---

## The refineField
The `refineField` is a custom field used to synchronize refinement across regions. It passes refinement information from the **master region** (the gas region) to communication regions (e.g., a dielectric). This ensures interface patches remain conformal.

**Example for a neighboring (dielectric) region:**
```cpp
dimensions      [ 0 0 0 0 0 0 0 ];

internalField   uniform 0;

boundaryField
{
    ".*"
    {
        type            calculated;
        value           uniform 0;
    }
    axis
    {
        type            empty;
    }
    wedge_0
    {
        type            wedge;
    }
    wedge_1
    {
        type            wedge;
    }
    dielectricRegion_to_gasRegion
    {
        type            mappedFixedInternalValue;
        fieldName       refineField;
        sampleRegion    gasRegion;
        samplePatch     gasRegion_to_dielectricRegion;
        sampleMode      nearestPatchFace;
        offset          ( 0 0 0 );
        offsetMode      uniform;
        value           uniform 0;
    }
}
```

---

## The decomposeParDict
To ensure successful multi-region AMR, the cells of two neighboring regions sharing a face **must** stay on the same processor. This is achieved using the `preserveFaceZones` and `preserveMappedInterface` constraints.

```cpp
numberOfSubdomains  4;

method          scotch;

constraints
{
    faces
    {
        type preserveFaceZones;
        zones
        (
            interfaceZone // Should be defined by topoSet or meshing tool
        );
        enabled true;
    }

    syncInterface
    {
        type    preserveMappedInterface;
        patches (gasRegion_to_dielectricRegion dielectricRegion_to_gasRegion);
    }
}
```

---

## Simulation Case Scenarios

### 1. Single Region Case
This is the simplest setup.
* No constraints in `decomposeParDict` are needed.
* No `refineField` is required.
* Enable AMR in `plasmaSimulationControls`:
```cpp
AMR
{
    enable      true;
}
```

### 2. Multi-Region Case (Explicit Coupling)
Used when regions are not implicitly coupled (Poisson solved per region, communicating via boundary conditions).
* Strict conformality is not mandatory as BCs handle mapping via interpolation.
* `refineField` is recommended but not mandatory.
* Regions can refine independently.
* No specific `decomposeParDict` constraints are required.

> ⚠️ **WARNING:**
> If you wish to restrict AMR **only** to the gas region, you must delete the `dynamicMeshDict` from the other regions, otherwise OpenFOAM will still attempt to apply refinement.

### 3. Multi-Region Case (Implicit Coupling)
Implicit coupling requires two strict conditions:
1.  Neighboring cells sharing a face **MUST** lie on the same processor.
2.  The interface meshes **MUST** be conformal.

**For Condition 1:** Use the `decomposeParDict` configuration shown previously. `preserveFaceZones` is applied during splitting; `preserveMappedInterface` is used during `blastAMR` balancing. If balancing is disabled, the latter may be omitted (though keeping it is recommended).

**For Condition 2:** You can ensure conformality in two ways:

#### Option A: Synchronized Refinement
Use `conformalAMRSync true` in `plasmaSimulationControls` and create `refineField` files for all regions.

**Example `refineField` for gasRegion:**
```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2412                                  |
|   \\  /    A nd           | Website:  [www.openfoam.com](https://www.openfoam.com)                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       volScalarField;
    location    "0/gasRegion";
    object      refineField;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [ 0 0 0 0 0 0 0 ];
internalField   uniform 0;

boundaryField
{
    ".*"
    {
        type            calculated;
        value           uniform 0;
    }
    axis   { type empty; }
    wedge_0 { type wedge; }
    wedge_1 { type wedge; }
    gasRegion_to_dielectricRegion
    {
        type            zeroGradient;
    }
}
```

**Example `refineField` for dielectricRegion:**
```cpp
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       volScalarField;
    location    "0/dielectricRegion";
    object      refineField;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [ 0 0 0 0 0 0 0 ];
internalField   uniform 0;

boundaryField
{
    ".*"
    {
        type            calculated;
        value           uniform 0;
    }
    axis   { type empty; }
    wedge_0 { type wedge; }
    wedge_1 { type wedge; }
    dielectricRegion_to_gasRegion
    {
        type            mappedFixedInternalValue;
        fieldName       refineField;
        sampleRegion    gasRegion;
        samplePatch     gasRegion_to_dielectricRegion;
        sampleMode      nearestPatchFace;
        offset          ( 0 0 0 );
        offsetMode      uniform;
        value           uniform 0;
    }
}
```

#### Option B: Patch Protection
If your mesh is not hexahedral (e.g., tetrahedral), refinement patterns are not always predictable and can break conformality. In this case, set `conformalAMRSync false` and protect the interface patches in the `dynamicMeshDict` of **both** regions.

**For gas/plasma region:**
```cpp
protectedPatches (gasRegion_to_dielectricRegion);
```
**For dielectric region:**
```cpp
protectedPatches (dielectricRegion_to_gasRegion);
```
This freezes the interface while allowing the rest of the domain to refine. **Processor adjacency constraints still apply.**

---

## The protectedPatches Feature
`protectedPatches` is also useful for complex multi-region setups. For example, in a three-region case (Gas -> Dielectric A -> Dielectric B), if you sync Gas and Dielectric A, you might want to protect the interface between Dielectric A and B to prevent unwanted refinement propagation.

---

## Extra Refinement in a Dielectric Region
You can add extra refinement conditions to a dielectric region beyond those inherited from the gas region. 
* **Explicit Coupling:** Straightforward independence.
* **Implicit Coupling:** If patches are not protected, we recommend using **coded refinement** for the dielectric. Ensure the `refineField` condition is always satisfied and overwrites any local criteria to maintain interface integrity.

---

## References

1.  **blastAMR**
    * [https://github.com/STFS-TUDa/blastAMR](https://github.com/STFS-TUDa/blastAMR)
2.  **blastFoam**
    * [https://github.com/synthetik-technologies/blastfoam](https://github.com/synthetik-technologies/blastfoam)
