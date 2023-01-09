# Nanoarchitectural ANalysis and Organization (NANO)
## Compiling
### Prerequisits
Project build is managed using the [cmake](https://cmake.org/) utility, and requires [ITK](https://itk.org/) as a dependency.
Make sure both these pre-requisits are satisfied before you start to build this repository.
All the instructions below are for a *nix based operating systems. 
We now list step-by-step process to copy, compile and then run this programme.
### Clone this repository
As a first step clone this repo to you local machine as follows.
```bash
$ git clone git@github.com:murailab/NANO.git
```
### Prepare for Build
As explained above compilation of this project is handled using `cmake`. To start the build create a build directory inside the 
downloaded repository
```bash
$ cd NANO
$ mkdir cmake-build
```
Point the `ITK_DIR` variable in `CMakeLists.txt` to appropriate location based on your ITK install. If you installed ITK 5.3
to its default location, then you can set `ITK_DIR` to something as follows.
```bash
$ sed -i '8s/.*/set(ITK_DIR \/usr\/local\/lib\/ITK\/lib\/cmake\/ITK-5.3)/' CMakeLists.txt
```
Make sure to adjust the path above based on your ITK install location.
### Generate Makefiles
Generate the required makefiles using cmake. 
```bash
$ cd cmake-build
$ cmake ..
```
### Compile
Finally, compile the project using the `make` utility. 
```bash
$ make -j 4
```
Depending on the number of cores available to you, instead of `4` you may use a higher number for faster compilation.

Once you complete the build process, an executable called `nano` should be generated in the `cmake-build` folder.
This is a commandline utility and the main entry point for running various computational pipelines
used in [[1]](#1). The [examples](#Examples) section below demonstrates how to perform some of the computations.
## Usage
Tools For Astrocyte Analysis:
```bash
$ nano -input <input_file_path> <-tool, -module n>  -mode 0
```

### Options::
	 -input
		 path to input file
	 -output
		 path to output file
	 -tool
		  If present  will run tool modes instead of submodules
	 -module
		 -module 0 : Computation related to distances  
		 -module 0 -mode 0 : compute inter-psd distances along astrocytic volume
		 -module 0 -mode 1 : compute distance from mitochondria to the astrocytic surface and beyond
	 -tool -mode n
		 Stand alone scripts/tools modes. Following tools are available:
		 -tool -mode 0 : Resample input image to isotropic grid (default)

## Examples

### Resample astrocyte volume to an isotropic grid
```bash
nano -tool -mode 0 -input <input_file_path>
```

### Distance from mitochondria to astrocyte surface and outside
```bash
nano -module 0 -mode 0 -input <input_folder_path>
```
This step assumes that `input_folder_path` contains 3D image stacks `astrocyte.tif`
containing the segmented astrocyte and `mito.tif` containing all the mitochondria 
inside the astrocyte.

### Table of various distance metrics for PSD regions
```bash
nano -module 0 -mode 1 <input_folder_path>
```
Make sure you have already generated various distance function before running this step.
(See step [above](#Distance-from-mitochondria-to-astrocyte-surface-and-outside))

## Citation
If you use any parts of this code, pleas consider citing the associated paper [[1]](#1)
## References
[1] Salmon, Christopher K., et al.
"Organizing Principles of Astrocytic Nanoarchitecture in the Mouse Cerebral Cortex."
Current Biology (2023).
