# Dependencies

This document lists all dependencies required or optionally used by the `foamPlasmaToolkit`, along with their purpose, importance, and the versions tested on the developer’s system.

---

## Tested Configurations

The toolkit has been tested on the following configurations:

| Ubuntu Version | OpenFOAM Version (OpenCFD Ltd.) |
|----------------|---------------------------------|
| 24.04 LTS      | v2412                           |
| 22.04 LTS      | v2412                           |

Additional Linux distributions and OpenFOAM (OpenCFD Ltd.) versions may work, but are not officially validated at this time.

## 1. Core Dependencies

These are required for compiling and running the main toolkit.

| Dependency              | Importance | Purpose                                        | Tested Version                |
|-------------------------|------------|------------------------------------------------|-------------------------------|
| **GCC / C++ compiler**  | Required   | Compiling OpenFOAM-dependent code              | 13.3.0                        |
| **CMake**               | Required   | Build system for optional utilities / examples | 3.28.3                        |
| **Make / Build tools**  | Required   | Used by wmake and OpenFOAM build system        | Ubuntu default (build-essential) |
| **MPI (OpenMPI)**       | Required   | Parallel execution support (OpenFOAM runs)     | Ubuntu default package  |

---

## 2. Optional Dependencies

These extend functionality but are *not* required for standard use.

| Dependency | Importance | Purpose | Tested Version |
|-----------|------------|---------|----------------|
| **PETSc** | Optional   | Advanced linear solvers and HPC features | 3.24.0 |
| **Python 3** | Optional | Utilities, scripts, or post-processing | 3.13.9 |
| **CUDA** | Optional*   | GPU acceleration for PETSc (GPU-enabled builds only) | 12.8 |

\* **CUDA is only required if PETSc is compiled with GPU support.** If PETSc is built in CPU-only mode, CUDA is not needed.

---

## 3. Tutorial-Specific Tools

Some tutorials require additional software not needed by the core toolkit.

| Dependency | Purpose | Notes |
|-----------|---------|-------|
| **ParaView** | Visualization of example cases | Any modern release should work |
| **gmsh** | Mesh generation for certain tutorials | Only needed for gmsh-based examples |

---

## 4. Python Dependencies

Some utilities and post-processing scripts in the `foamPlasmaToolkit` rely on a small set of Python packages. These are optional and only needed if you plan to use the Python-based tools.

### Python Packages

| Package        | Purpose                        | Tested Version |
|----------------|--------------------------------|----------------|
| **numpy**      | Numerical operations            | 2.3.4          |
| **scipy**      | Scientific functions           | 1.16.3         |
| **matplotlib** | Plotting and visualization     | 3.10.7         |
| **pandas**     | Data handling for analysis     | 2.3.3          |
| **gmsh**       | Mesh generation for specific tutorials | 4.15.0 |

These packages are listed in the file: [`docs/python_dependencies.txt`](./python_dependencies.txt)

### Installing the Python Dependencies

To install all Python requirements at once, run:

#### Using pip
```bash
python3 -m pip install -r docs/python_dependencies.txt
```

#### Using conda
```bash
conda install --file docs/python_dependencies.txt
```

---

## 5. Notes

- Optional dependencies should only be installed if you need the corresponding features.  
- PETSc and petsc4Foam require correct linking with your OpenFOAM installation; see [`docs/petsc4foam.md`](./petsc4foam.md) for details.
