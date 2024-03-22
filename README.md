# GEOSldas Fixture

This document explains how to build, set up, and run the GEOS land modeling and data assimilation system (`GEOSldas`) on the most common systems used by GMAO.  Additional steps are needed on other systems.

## How to Build GEOSldas

### Step 1: Load the Build Modules

Load the `GEOSenv` module provided by the GMAO Software Infrastructure team.  It contains the latest `git`, `CMake`, and `mepo` modules and must be loaded in any interactive window that is used to check out and build the model.

```
module use -a (path)
module load GEOSenv
```

where `(path)` depends on the computing system; at NCCS, `(path)` also depends on the operating system (SLES12 on Skylake and Cascade Lake nodes; SLES15 on Milan nodes, as of Jan. 2024):

| System        | Path                                              |
| ------------- |---------------------------------------------------|
| NCCS Discover | `/discover/swdev/gmao_SIteam/modulefiles-SLES12`  |
|               | `/discover/swdev/gmao_SIteam/modulefiles-SLES15`  |
| NAS           | `/nobackup/gmao_SIteam/modulefiles`               |
| GMAO desktops | `/ford1/share/gmao_SIteam/modulefiles`            |

Step 1 can be coded into the user's shell configuration file (e.g., `.bashrc` or `.cshrc`). See the [GEOSgcm Wiki](https://github.com/GEOS-ESM/GEOSgcm/wiki/) for sample shell configuration files.

### Step 2: Obtain the Model

For development work, clone the _entire_ repository and use the `develop` branch as your starting point:
```
git clone -b develop git@github.com:GEOS-ESM/GEOSldas.git
```
For science runs, you can also obtain a specific tag or branch _only_ (as opposed to the _entire_ repository), e.g.:
```
git clone -b v17.9.1 --single-branch git@github.com:GEOS-ESM/GEOSldas.git
```


### Step 3: Build the Model

To build the model in a single step, do the following from a head node:
```
cd ./GEOSldas
parallel_build.csh
```
This checks out all the external repositories of the model (albeit only on the first run, [see subsection on mepo below](#mepo)!) and then builds and installs the model. 

At **NCCS**, the default is to build GEOSldas on SLES12 (Skylake or Cascade Lake nodes); to build GEOSldas on SLES15 (Milan nodes), use `parallel_build.csh -mil`.

The resulting model build is found in `build[-SLESxx]/`, and the installation is found in `install[-SLESxx]/`, with setup scripts like `ldas_setup` in `install[-SLESxx]/bin`.

To obtain a build that is suitable for debugging, use `parallel_build.csh -debug`, which builds in `build-Debug[-SLESxx]/` and installs in `install-Debug[-SLESxx]/`.  There is also an option for aggressive  optimization.  For details, see the [GEOSldas Wiki](https://github.com/GEOS-ESM/GEOSldas/wiki).

Instructions for building the model in multiple steps are provided below.

---

## How to Set Up (Configure) and Run GEOSldas


a) At **NCCS**, GEOSldas must be built, configured, and run on the same operating system. To run GEOSldas on Milan nodes (SLES15), start with `ssh discover-mil`.

b) Set up the job as follows:

```
cd (build_path)/GEOSldas/install[-SLESxx]/bin
source g5_modules                        [for bash or zsh: source g5_modules.[z]sh]
./ldas_setup setup [-v]  (exp_path)  ("exe"_input_filename)  ("bat"_input_filename)
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
ldas_setup sample --batinp > YOUR_batinp.txt
```

Edit these sample files following the examples and comments within the sample files.

The ldas_setup script creates a run directory and other directories at:
`[exp_path]/[exp_name]`

Configuration input files are created at:
`[exp_path]/[exp_name]/run`

For more options and documentation, use any of the following:
```
ldas_setup        -h
ldas_setup sample -h
ldas_setup setup  -h
```

c) Configure the experiment output by editing the ```./run/HISTORY.rc``` file as needed.

d) Run the job:
```
cd [exp_path]/[exp_name]/run/
sbatch lenkf.j
```

At **NCCS**, the appropriate SLURM directive `#SBATCH --constraint=[xxx]` is automatically added into `lenkf.j` depending on the operating system.

For more information, see the files in `./doc/`. Moreover, descriptions of the configuration (resource) parameters are included in the sample "exeinp" and "batinp" files that can be generated using `ldas_setup`.



-----------------------------------------------------------------------------------

## Additional Information

### How to Build the Model in Multiple Steps

The steps detailed below are essentially those performed by `parallel_build.csh` in Step 3 above. Either method should yield identical builds.

#### mepo

The GEOSldas is comprised of a set of sub-repositories. These are
managed by a tool called [mepo](https://github.com/GEOS-ESM/mepo). To
clone all the sub-repositories, you can run `mepo clone` inside the fixture:
```
cd GEOSldas
mepo init
mepo clone
```
External sub-repositories are stored in directories pre-faced with `@`. After `parallel_build.csh` has run once and created `./@env/`, `parallel_build.csh` skips `mepo clone` in subsequent runs. This means that the sub-repositories in your sandbox could get out of sync with the GEOSldas repository while you are working on your sandbox, which may result in a difficult-to-understand build error when `parallel_build.csh` is used. If this happens, try a fresh clone, or use [mepo commands](https://github.com/GEOS-ESM/mepo/wiki) to update the sub-repositories manually.

#### Load Compiler, MPI Stack, and Baselibs
On tcsh:
```
source @env/g5_modules
```
or on bash:
```
source @env/g5_modules.sh
```

#### Create Build Directory
We currently do not allow in-source builds of GEOSldas. So we must make a directory:
```
mkdir build
```

#### Run CMake
CMake generates the Makefiles needed to build the model.
```
cd build
cmake .. -DBASEDIR=$BASEDIR/Linux -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_INSTALL_PREFIX=../install
```
This installs into a directory parallel to your `build` directory. If you prefer to install elsewhere change the path in:
```
-DCMAKE_INSTALL_PREFIX=<path>
```
and CMake will install there.

#### Build and Install with Make
```
make -j6 install
```
If you are at NCCS, you **should** run `make -j6 install` on an interactive _compute_ node.


## Contributing

Please check out our [contributing guidelines](CONTRIBUTING.md).

## License

All files are currently licensed under the [Apache-2.0 license (`LICENSE`)](LICENSE).

Previously, the code was licensed under the [NASA Open Source Agreement, Version 1.3 (`LICENSE-NOSA`)](LICENSE-NOSA).
