# GEOSldas Fixture

This document explains how to build, set up, and run the GEOS land modeling and data assimilation system (`GEOSldas`).

## How to Build GEOSldas

### Step 1: Load the Build Modules  

Load the `GEOSenv` module provided by the GMAO Software Infrastructure team.  It contains the latest `git`, `CMake`, and `manage_externals` modules and must be loaded in any interactive window that is used to check out and build the model.

```
module use -a (path)
module load GEOSenv
```

where `(path)` depends on the computer and operating system: 

| System        | Path                                              |
| ------------- |---------------------------------------------------|
| NCCS SLES11   | `/discover/swdev/gmao_SIteam/modulefiles-SLES11`  |
| NCCS SLES12   | `/discover/swdev/gmao_SIteam/modulefiles-SLES12`  |
| NAS           | `/nobackup/gmao_SIteam/modulefiles`               |
| GMAO desktops | `/ford1/share/gmao_SIteam/modulefiles`            |


For NCCS, you can add the following to your `.cshrc`:
```
if ( ! -f /etc/os-release ) then
   module use -a /discover/swdev/gmao_SIteam/modulefiles-SLES11
else
   module use -a /discover/swdev/gmao_SIteam/modulefiles-SLES12
endif
module load GEOSenv
```


### Step 2: Obtain the Model

For development work, clone the _entire_ repository and use the `develop` branch as your starting point (equivalent to the `UNSTABLE` tag in the old CVS repository):
```
git clone -b develop git@github.com:GEOS-ESM/GEOSldas.git
```
For science runs, you can also obtain a specific tag or branch _only_ (as opposed to the _entire_ repository), e.g.: 
```
git clone -b v17.9.0-beta.3 --single-branch git@github.com:GEOS-ESM/GEOSldas.git
```


### Step 3: Build the Model

To build the model in a single step, do the following:
```
cd ./GEOSldas
parallel_build.csh
``` 
from a head node. Doing so will checkout all the external repositories of the model and build it. When done, the resulting model build will be found in `build/` and the installation will be found in `install/`, with setup scripts like `ldas_setup` in `install/bin`. 

To obtain a build that is suitable for debugging, you can run `parallel_build.csh -debug`, which will build in `build-Debug/` and install in `install-Debug/`.

See below for how to build the model in multiple steps.

---

## How to Set Up and Run GEOSldas

a) Set up the job as follows:

```
cd (build_path)/GEOSldas/install/bin
source g5_modules
./ldas_setup setup [-v] [--runmodel]  (exp_path)  ("exe"_input_filename)  ("bat"_input_filename)
```  

where

| Parameter              | Description                                              |
| -----------------------|----------------------------------------------------------|
| `build_path`           | path to build directory                                  |
| `exp_path`             | path of desired experiment directory                     |
| `"exe"_input_filename` | filename (with path) of "experiment" inputs              |
| `"bat"_input_filename` | filename (with path) of "batch" (job scheduler) inputs   |

The three arguments for `ldas_setup` are positional and must be ordered as indicated above.

The latter two files contain essential information about the experiment setup. 
Sample files can be generated as follows:
```        
ldas_setup sample --exeinp > YOUR_exeinp.txt
ldas_setup sample --batinp > YOUR_exeinp.txt
```

Edit these sample files following the examples and comments within the sample files.  

The ldas_setup script creates a run directory and other directories at:
`[exp_path]/[exp_name]`

Configuration input files will be created at:
`[exp_path]/[exp_name]/run`

For more options and documentation, use any of the following:
```
ldas_setup        -h
ldas_setup sample -h
ldas_setup setup  -h
```

b) Configure the experiment output by editing the ```./run/HISTORY.rc``` file as needed.

c) Run the job:
```
cd [exp_path]/[exp_name]/run/
sbatch lenkf.j
```

For more information, see the README files and ppt tutorial in `./src/Applications/LDAS_App/doc/`.

-----------------------------------------------------------------------------------

## Additional Information

### How to Build the Model in Multiple Steps

The steps detailed below are essentially those performed by `parallel_build.csh` in Step 3 above. Either method should yield identical builds.

##### Checkout externals
```
cd GEOSldas
checkout_externals
```

##### Load Compiler, MPI Stack, and Baselibs
On tcsh:
```
source @env/g5_modules
```
or on bash:
```
source @env/g5_modules.sh
```

##### Create Build Directory
We currently do not allow in-source builds of GEOSldas. So we must make a directory:
```
mkdir build
```
The advantages of this is that you can build both a Debug and Release version with the same clone if desired.

##### Run CMake
CMake generates the Makefiles needed to build the model.
```
cd build
cmake .. -DBASEDIR=$BASEDIR/Linux -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_INSTALL_PREFIX=../install
```
This will install to a directory parallel to your `build` directory. If you prefer to install elsewhere change the path in:
```
-DCMAKE_INSTALL_PREFIX=<path>
```
and CMake will install there.

##### Build and Install with Make
```
make -j6 install
```
If you are using SLES12 at NCCS, you **should** run `make -j6 install` on an interactive _compute_ node.  

