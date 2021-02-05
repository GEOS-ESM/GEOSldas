# GEOSldas Configuration Files

This document describes the files involved in the specification and pre-processing of the GEOSldas configuration parameters. 

Note that the values of the configuration parameters that are used during the execution of `GEOSldas.x` are written into the GEOSldas log file (`./output/*/rc_out/Y*/M*/*.ldas_log.*.txt`).

---
### `"exeinp"` and `"batinp"` files

Inputs to `ldas_setup` that contain user-defined configuration information.  

Use `ldas_setup` to create sample files that contain descriptions of  the configuration parameters.  

For details, see [README.md](https://github.com/GEOS-ESM/GEOSldas/blob/main/README.md).


---
### `LDASsa_SPECIAL_inputs_*.nml`
 
Optional Fortran namelist (nml) files that contain additional user-defined configuration information for ensemble perturbations and data assimilation.  

The path to these nml files is specified in the `"exeinp"` file. 

During `ldas_setup`, **default** configuration files (`LDASsa_DEFAULT_inputs_*.nml`) are created in the experiment `./run` directory.  These files contain a complete set of the available configuration parameters with descriptions.  The default configuration is a single-member land model simulation without perturbations and without data assimilation. The "DEFAULT" nml files must be present in the experiment `./run` directory and should not be edited.  

To run GEOSldas with ensemble perturbations and data assimilation, users must create "SPECIAL" nml files (`LDASsa_SPECIAL_inputs_*.nml`) that contain the desired settings of the parameters.  Only the nml parameters that are different from those in the "DEFAULT" files need to be included in the "SPECIAL" nml files.

Parameters are grouped into three separate nml files:
* `ensprop` : Perturbations applied during the land model ensemble propagation step. 
* `ensupd`  : Assimilated observations and parameters of the ensemble-based analysis.
* `catbias` : Dynamic bias estimation (defunct).


---
### `CAP.rc` and `cap_restart`

Created by `ldas_setup`.

`CAP.rc` contains the experiment start/end times and the number/length of job segments. See documentation in `"exeinp"` file.

`cap_restart` contains the start time of the next job segment.  Note that `GEOSldas.x` reads the start time from this file, not from `CAP.rc` file.

To extend a simulation past the originally specified end date used during `ldas_setup`, users can change the end date in `CAP.rc` and resubmit the `lenkf.j` job script.  In this particular case, it is not necessary to run `ldas_setup` again.

If an experiment did not complete successfully (e.g., because of a system downtime or because a job exceeded the wall-time limit), users must ensure that the start time in `cap_restart` matches the time stamps of the linked restart files in `../input/restart/` and then resubmit `lenkf.j`.


---
### `LDAS.rc`

Contains parameters needed to run the main executable (`GEOSldas.x`), including information on land model configuration, version, grid, time steps, processor layout.

Created in the experiment `./run` directory during `ldas_setup`.  Merges information from the `"exeinp"` and `"batinp"` files and the resource parameter template files for GEOSldas, i.e., [GEOSldas_LDAS.rc](https://github.com/GEOS-ESM/GEOSldas/blob/main/src/Applications/LDAS_App/GEOSldas_LDAS.rc) and  [GEOS_SurfaceGridComp.rc](https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/main/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Shared/GEOS_SurfaceGridComp.rc) (from the linked GEOSgcm_GridComp repository).

**Users are generally discouraged from editing** `LDAS.rc`.  For example, editing `MET_PATH` does not change the associated directory link.  Instead, users should run `ldas_setup` to create `LDAS.rc`.

`LDAS.rc` can be useful as a compact overview of the experiment configuration. 

---
### `lenkf.j`

GEOSldas job script created by `ldas_setup`.  The beginning of `lenkf.j` contains configuration parameters for the resource manager (e.g., SBATCH).

**Users are generally discouraged from editing** `lenkf.j`, but the number of processors (`ntasks`) can be edited safely without going through `ldas_setup`.


