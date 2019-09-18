# GEOS LDAS Fixture

## How to build GEOS LDAS

### Preliminary Steps

#### Load Build Modules

In your `.bashrc` or `.tcshrc` or other rc file add a line:

##### NCCS (SLES11)

```
module use -a /discover/swdev/gmao_SIteam/modulefiles-SLES11
```

##### NAS
```
module use -a /nobackup/gmao_SIteam/modulefiles
```

##### GMAO Desktops
On the GMAO desktops, the SI Team modulefiles should automatically be
part of running `module avail` but if not, they are in:

```
module use -a /ford1/share/gmao_SIteam/modulefiles
```

Also do this in any interactive window you have. This allows you to get module files needed to correctly checkout and build the model.

Now load the `GEOSenv` module:
```
module load GEOSenv
```
which obtains the latest `git`, `CMake`, and `manage_externals` modules.

#### Obtain the Model

```
git clone git@github.com:GEOS-ESM/GEOSldas.git
```

---

### Single Step Building of the Model

If all you wish is to build the model, you can run `parallel_build.csh` from a head node. Doing so will checkout all the external repositories of the model and build it. When done, the resulting model build will be found in `build/` and the installation will be found in `install/` with setup scripts like `ldas_setup` in `install/bin`.

#### Develop Version of LDAS

The user will notice two files in the main directory: `Externals.cfg` and `Develop.cfg`. The difference between these two is that `Externals.cfg` always refers to stable tested released subrepositories. The `Develop.cfg` points to the `develop` branches of the subrespositories, or more recent tags not ready for prime time. This is equivalent in the CVS days of the difference between a stable tag and the development `-UNSTABLE` tag. In order to build the `Develop.cfg` version of the model with `parallel_build.csh` do:

```
git checkout develop
parallel_build.csh -develop
```
Note you have to first checkout the develop branch, and then also use Develop.cfg (This may be simplified in the future)

#### Debug Version of LDAS

To obtain a debug version, you can run `parallel_build.csh -debug` which will build with debugging flags. This will build in `build-Debug/` and install into `install-Debug/`.

---

### Multiple Steps for Building the Model

The steps detailed below are essentially those that `parallel_build.csh` performs for you. Either method should yield identical builds.

##### Checkout externals
```
cd GEOSldas
checkout_externals
```

#### Build the Model

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

---

## Setup up a run
```
cd ../(some_architecture)/bin
source g5_modules
./ldas_setup setup [-v] [--runmodel]  exp_path  "exe"_input_filename  "bat"_input_filename
```  
where

>exp_path             = path of desired experiment directory

>"exe"_input_filename = filename (with path) of "experiment" inputs

>"bat"_input_filename = filename (with path) of "batch" inputs

must be ordered as above (positional arguments).

The latter two files contain essential information about the experiment setup. 
Sample files can be generated as follows:
```        
ldas_setup sample --exeinp > YOUR_exeinp.txt
ldas_setup sample --batinp > YOUR_exeinp.txt
```

Edit these sample files (see comments within sample files).  See README files
and ppt tutorial (in ./src/Applications/LDAS_App/doc/) for more information.

The ldas_setup script creates a run directory and other directories at:
```[exp_path]/[exp_name]```

Configuration input files will be created at:

```[exp_path]/[exp_name]/run```

For more options and documentation run any of the following:
```
ldas_setup        -h
ldas_setup sample -h
ldas_setup setup  -h
```

Configure experiment output by editing the ```HISTORY.rc``` file.

---

## Run a job:

	cd [exp_path]/[exp_name]/run/

	sbatch lenkf.j

See ppt tutorial (in ./src/Applications/LDAS_App/doc/) for more information about how to run GEOSldas.

