# CUDA Toolkit Installation Guide

The `foamPlasmaToolkit` can optionally be used on systems with **NVIDIA GPUs** to accelerate linear
algebra operations when PETSc is compiled with CUDA support.

CUDA is **not required** for normal CPU-only builds.

## References

This guide follows the official NVIDIA installation documentation for CUDA and driver setup:

- **NVIDIA CUDA Toolkit Documentation**  
  <https://docs.nvidia.com/cuda/>

- **NVIDIA CUDA Downloads**  
  <https://developer.nvidia.com/cuda-toolkit-archive>

- **NVIDIA Driver Downloads**  
  <https://www.nvidia.com/Download/index.aspx>

## Tested Configuration

CUDA support has been tested with:

- **CUDA Toolkit 12.8**
- **NVIDIA Driver 573.24**
- **Ubuntu 24.04 LTS**

Other versions *may* work, but are not officially verified.

## 1. Install NVIDIA GPU Drivers

Make sure NVIDIA drivers are installed (version **573.24** or newer):

```bash
nvidia-smi
```

If the command prints driver info (e.g., version 573.24), you are ready to proceed. 
Example expected output:

```bash
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 570.152      Driver Version: 573.24       CUDA Version: 12.8     |
+-----------------------------------------------------------------------------+
```

If nvidia-smi returns command not found or errors, install the recommended driver:

```bash
sudo ubuntu-drivers autoinstall
sudo reboot
```

After reboot, verify again:

```bash
nvidia-smi
```

If you prefer manual installation or need a specific driver version, see the official NVIDIA page:

➡ https://www.nvidia.com/Download/index.aspx

## 2. Install CUDA Toolkit

⚠️ **IMPORTANT:** CUDA version 12.8 has been tested. Other versions may work but are not guaranteed.

### Option A — Ubuntu Native

If you are using **Ubuntu 24.04 LTS** and want to install **CUDA 12.8**, follow the commands below. Otherwise, visit the **[CUDA Toolkit Archive](https://developer.nvidia.com/cuda-toolkit-archive)** and select the appropriate version for your system.


```bash
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/cuda-ubuntu2404.pin
sudo mv cuda-ubuntu2404.pin /etc/apt/preferences.d/cuda-repository-pin-600
wget https://developer.download.nvidia.com/compute/cuda/12.8.0/local_installers/cuda-repo-ubuntu2404-12-8-local_12.8.0-570.86.10-1_amd64.deb
sudo dpkg -i cuda-repo-ubuntu2404-12-8-local_12.8.0-570.86.10-1_amd64.deb
sudo cp /var/cuda-repo-ubuntu2404-12-8-local/cuda-*-keyring.gpg /usr/share/keyrings/
sudo apt-get update
sudo apt-get -y install cuda-toolkit-12-8
```

### Option B — WSL (If installing inside WSL2)

If you are using **Ubuntu 24.04 LTS** and want to install **CUDA 12.8**, follow the commands below. Otherwise, visit the **[CUDA Toolkit Archive](https://developer.nvidia.com/cuda-toolkit-archive)** and select the appropriate version for your system.

```bash
wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-wsl-ubuntu.pin
sudo mv cuda-wsl-ubuntu.pin /etc/apt/preferences.d/cuda-repository-pin-600
wget https://developer.download.nvidia.com/compute/cuda/12.8.0/local_installers/cuda-repo-wsl-ubuntu-12-8-local_12.8.0-1_amd64.deb
sudo dpkg -i cuda-repo-wsl-ubuntu-12-8-local_12.8.0-1_amd64.deb
sudo cp /var/cuda-repo-wsl-ubuntu-12-8-local/cuda-*-keyring.gpg /usr/share/keyrings/
sudo apt-get update
sudo apt-get -y install cuda-toolkit-12-8
```

## 3. Add CUDA to Environment

Open the `.bashrc` file in your home directory using a text editor (note the leading dot):

```bash
gedit ~/.bashrc
```

Then add the following lines:

```bash
export PATH=/usr/local/cuda-12.8/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-12.8/lib64:$LD_LIBRARY_PATH
```

⚠️ **IMPORTANT:** Adjust the paths to match the CUDA version you installed and its installation location.

Reload:

```bash
source ~/.bashrc
```

Verify CUDA installation:

```bash
nvcc --version
```

Expected result includes `release 12.8`.


