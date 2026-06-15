# blastAMR Installation Guide

`SoPLASMA` optionally supports **blastAMR**, which provides adaptive mesh refinement
on hexahedral and polyhedral meshes in 1D, 2D, 2D-axisymmetric, and 3D configurations,
with dynamic load-balancing for parallel simulations.

blastAMR is **not required** for basic use. Use it only if you need dynamic adaptive
mesh refinement capabilities.

## References

This documentation and functionality rely on the following project:

- **blastAMR** - Adaptive mesh refinement library for OpenFOAM  
  https://github.com/STFS-TUDa/blastAMR

---

## blastAMR Installation / Setup

> **Note:** *`blastAMR` is automatically compiled during the build of `SoPLASMA`.
> When you install the toolkit, it will build `blastAMR` automatically.*

> *If `SoPLASMA` was installed **before** blastAMR, you need to compile it manually
> by following the steps below.*

### 1. Go to blastAMR ThirdParty directory (SoPLASMA)

```bash
cd $SoPLASMA_THIRD_PARTY/blastAMR
```

### 2. Build blastAMR

```bash
./Allwmake -prefix=openfoam
```

### 3. Ensure that blastAMR is installed

```bash
foamHasLibrary -verbose amrDynamicFvMesh   # should print: Can load "amrDynamicFvMesh"
```

---

## Environment Setup

### Add blastAMR to your environment

Open your `~/.bashrc` in an editor:

```bash
gedit ~/.bashrc
```

Add the following lines:

```bash
export AMRLB_PROJECT=$SoPLASMA_THIRD_PARTY/blastAMR
export FOAM_CODE_TEMPLATES=$AMRLB_PROJECT/etc/codeTemplates/dynamicCode
```

> ⚠️ **IMPORTANT:** The `FOAM_CODE_TEMPLATES` export is only required if you use
> **coded** (run-time compiled) AMR indicators or error estimators. It is not needed
> for standard built-in AMR functionality.

Reload:

```bash
source ~/.bashrc
```

---

## Use blastAMR

### 1. Add blastAMR libraries to `controlDict`

Add the following to the `libs` entry of your `system/controlDict`:

```cpp
libs
(
    "libamrDynamicFvMesh.so"
    "libamrDynamicMesh.so"
    "libamrLoadPolicies.so"
    "libamrCPULoad.so"
    "libamrFunctionObjects.so"
    "libamrIndicators.so"
    "libamrLagrangian.so"
);
```

> **Note:** Not all libraries are required for every case. At minimum,
> `libamrDynamicFvMesh.so` and `libamrDynamicMesh.so` are needed for basic AMR.
> The remaining libraries enable load balancing, function objects, field-based
> indicators, and Lagrangian particle support respectively.

### 2. Set dynamic mesh in `constant/dynamicMeshDict`

The mesh type must be set to one of the blastAMR refiners. For example, for a
2D axisymmetric wedge mesh:

```cpp
dynamicFvMesh   adaptiveFvMesh;
```
