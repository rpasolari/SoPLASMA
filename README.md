# SoPLASMA

**SoPLASMA** is a collection of solvers, models, boundary conditions, and utilities developed for plasma applications in **OpenFOAM**. It aims to extend the capabilities of OpenFOAM to simulate plasma-related phenomena, including plasma-fluid interactions, plasma chemistry, and  coupled electromagnetic effects.


##	Description

This toolkit provides an evolving framework for implementing plasma physics features within
OpenFOAM. It is designed for research and development purposes, focusing on modularity and compatibility with the OpenFOAM v2412 (OpenCFD Ltd.).  

##	Status

The project is currently in its **initial development phase**. Functionality, structure, and documentation are expected to change as new components are added.

##	License

This project is distributed under the terms of the **GNU General Public License v3.0 (GPLv3)**.  
For full license details, see the [LICENSE](LICENSE) file.

© 2026 Rention Pasolari

This software is not part of OpenFOAM, but it is developed using the OpenFOAM framework and linked against OpenFOAM libraries (v2412).

## Requirements

### Core Dependency
- **OpenFOAM (OpenCFD Ltd. release)**  
  *(Developed and primarily tested with version **v2412**.)*

This toolkit relies on the OpenCFD Ltd. version of OpenFOAM for compilation and compatibility.

### Additional Dependencies

The toolkit includes a number of additional components that may be required depending on what you intend to use:

- **Crucial dependencies** needed for building the core functionality.
- **Optional dependencies** that enable extended features. 
- **Tutorial-specific tools** required only for certain example cases.  
- **External libraries**, such as **petsc4foam** *(see [`docs/getting_started/petsc4foam.md`](docs/getting_started/petsc4foam.md))*

For a complete and recommended overview of all dependencies, see **[`docs/getting_started/dependencies.md`](docs/getting_started/dependencies.md)**.

## Tested Platforms

The toolkit has been tested on the following configurations:

| Ubuntu Version | OpenFOAM Version (OpenCFD Ltd.) |
|----------------|---------------------------------|
| 24.04 LTS      | v2412                           |
| 22.04 LTS      | v2412                           |

Additional Linux distributions and OpenFOAM (OpenCFD Ltd.) versions may work, but are not officially validated at this time.

---

## Installation / Setup

## 1. Install Required Dependencies

Before installing OpenFOAM and the SoPLASMA, make sure all required dependencies are installed. Refer to **[`docs/getting_started/dependencies.md`](docs/getting_started/dependencies.md)** and install all packages marked as *required*.


## 2. Install OpenFOAM v2412 (required)

This toolkit is designed **for OpenFOAM v2412 (OpenCFD Ltd.)**. It has not been tested in other  versions.

If you already have OpenFOAM installed, skip to **Step 2 (Download and build SoPLASMA)**.


### A. Install OpenFOAM prerequisites

```bash
sudo apt-get update
sudo apt-get install build-essential autoconf autotools-dev cmake gawk gnuplot
sudo apt-get install flex libfl-dev libreadline-dev zlib1g-dev openmpi-bin libopenmpi-dev mpi-default-bin mpi-default-dev
sudo apt-get install libgmp-dev libmpfr-dev libmpc-dev
sudo apt-get install libfftw3-dev libscotch-dev libptscotch-dev libboost-system-dev libboost-thread-dev libcgal-dev
```

### B. Get the source code

Choose the directory where you want OpenFOAM v2412 to be installed.

Common system-wide installation paths are:

- `/usr/lib/openfoam/openfoam2412/`
- `/opt/openfoam2412/`

However, both locations require **sudo** privileges. This can be inconvenient when compiling ThirdParty packages or running tests, as some steps may require elevated permissions.

To avoid this, you can install OpenFOAM in your **home directory**, which allows you to compile and work without using `sudo`.

#### Example: Install OpenFOAM v2412 in your home directory

```bash
cd ~
mkdir OpenFOAM
cd OpenFOAM
```

You can download the official OpenFOAM v2412 source code directly from **openfoam.com**:

```bash
wget -O -  http://dl.openfoam.com/source/v2412/OpenFOAM-v2412.tgz | tar xvz
wget -O -  http://dl.openfoam.com/source/v2412/ThirdParty-v2412.tgz | tar xvz
```

If `wget` fails to establish a connection, you can download the archives manually from the official mirrors and extract them into your chosen OpenFOAM directory.

Manual download links:

- [OpenFOAM-v2412.tgz (source code)](https://sourceforge.net/projects/openfoam/files/v2412/OpenFOAM-v2412.tgz/download)  
- [ThirdParty-v2412.tgz (third-party libraries)](https://sourceforge.net/projects/openfoam/files/v2412/ThirdParty-v2412.tgz/download)


After downloading, extract both archives inside your OpenFOAM installation folder:

```bash
tar -xvzf OpenFOAM-v2412.tgz
tar -xvzf ThirdParty-v2412.tgz
```

### C. Set the environment variables

Open the `.bashrc` file in your home directory with your preferred editor (note the leading dot):

```bash
nano ~/.bashrc
# or
gedit ~/.bashrc
```

Add the following line at the end of the file:

```bash
source $HOME/OpenFOAM/OpenFOAM-v2412/etc/bashrc
```

Save the file, then reload your shell configuration:

```bash
source ~/.bashrc
```

### D. Install ThirdPart and OpenFOAM

First, install the ThirdParty packages:

```bash
cd ~/OpenFOAM/ThirdParty-v2412
./Allwmake -j
```

After the ThirdParty build completes, install OpenFOAM itself:

```bash
cd ~/OpenFOAM/OpenFOAM-v2412
./Allwmake -j
```

When the compilation finishes, verify that OpenFOAM is correctly installed:

```bash
foamVersion   # should print: OpenFOAM-v2412
```

## 3. Install SoPLASMA

### A. Get the source code

Choose the directory where you want to install **SoPLASMA**, then download the toolkit by cloning the GitHub repository.

#### Example: Install SoPLASMA in your home directory

```bash
cd ~
git clone https://github.com/rpasolari/SoPLASMA.git
cd SoPLASMA
```

Activate the toolkit environment (so OpenFOAM can detect its extensions):
```bash
source etc/bashrc
```

⚠️ **IMPORTANT:** To avoid having to run this manually every time, add it to your `~/.bashrc` as shown below.

Open the `.bashrc` file in your home directory using any text editor (note the leading dot):

```bash
gedit ~/.bashrc
```

Then add the following line at the end of the file (adjust the path to match where you installed the toolkit):

```bash
source $HOME/SoPLASMA/etc/bashrc
```

Save the file, then reload your shell configuration:

```bash
source ~/.bashrc
```

To verify that the environment was successfully loaded, open a **new terminal** and run:

```bash
echo $SoPLASMA   # should print the toolkit folder
```

### B. Build SoPLASMA
From inside the SoPLASMA directory, compile everything:

```bash
./Allwmake
```
If the compilation finishes without errors, the installation is complete.

⚠️ **IMPORTANT:** Some ThirdParty packages included with the SoPLASMA (but not required as prerequisites) may fail to build (for example, **petsc4foam**). To install these optional components, refer to the corresponding instructions in the `docs/` directory.


## Contributors

**Rention Pasolari** <r.pasolari@gmail.com>  
Software architecture, implementation, testing, validation, documentation, maintenance

**Konstantinos Kourtzanidis** <kkourtza@gmail.com>  
Conceptual design, physical modeling, debugging, testing, validation, scientific supervision
