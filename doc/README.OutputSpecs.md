# GEOSldas Output Specifications

### Overview

Most diagnostic output from GEOSldas is generated through the **HISTORY** functionality provided by [MAPL](https://github.com/GEOS-ESM/MAPL).

Restart/checkpoint files are also written through MAPL (although not through HISTORY).  The key configuration parameters that control checkpoint file output are described in the sample "exeinp" files that can be generated using `ldas_setup`.  For details, see [README.md](https://github.com/GEOS-ESM/GEOSldas/blob/main/README.md).

Notable exceptions from the MAPL-generated output include:
* log files [ASCII], 
* observation-space "ObsFcstAna" data assimilation diagnostics [BINARY], and
* SMAP L4_SM-specific "aup" data assimilation diagnostics files (available _**only**_ for simulations in EASEv2_M09 tile space) [BINARY].

Output of the latter two sets of files can be turned on/off in the `[NML_INPUT_PATH]/LDASsa_SPECIAL_inputs_ensupd.nml` configuration file, and Matlab readers are available in the 
[src/Applications/LDAS_App/util](https://github.com/GEOS-ESM/GEOSldas/tree/main/src/Applications/LDAS_App/util) directory.


### MAPL HISTORY output

As part of `ldas_setup`, a sample `HISTORY.rc` configuration file is created in the experiment's `./run` directory.  Users specify the desired output by editing `HISTORY.rc`.  

`HISTORY.rc` defines a number of output file "Collections", each of which contains one or more output variables.  Output can be gridded ("2d") or in the native tile space ("1d").  

All variables contained in a given Collection are written:
* on the same ("2d") grid (if gridded), 
* at the same frequency, and 
* with either time-average ("tavg") or instantaneous ("inst") sampling mode.  

In the following example, the `tavg3_2d_lnd_Nx` Collection contains time-average ("tavg"), 3-hourly ("3"), gridded ("2d") data, and the inst1_1d_lfs_Nt Collection contains snapshot/instantaneous ("inst"), 1-hourly, tile-space ("1d") data. 
```
COLLECTIONS:
            'tavg3_2d_lnd_Nx'
            'inst1_1d_lfs_Nt'
           ::
```

To be available for output through MAPL HISTORY, a variable must be defined as an `ExportSpec` in a `GEOS_*GridComp.F90` file. 
The list of variables ("fields") in each Collection consists of three columns: 
- (column 1) variable name in `GEOS_[GCNAME]GridComp.F90` file,
- (column 2) GridComp name [GCNAME], and 
- (column 3) user-specified variable name that appears in nc4 output (optional).

For example, to write the "WCSF" and "WCRZ" variables (surface and root-zone soil moisture) from the Catchment model ("CATCH") GridComp under the names "SFMC" and "RZMC", respectively, the definition of the `tavg3_2d_lnd_Nx` Collection in `HISTORY.rc` must include the following:
```
 tavg3_2d_lnd_Nx.mode:        'time-averaged',
 tavg3_2d_lnd_Nx.frequency:   030000,
 tavg3_2d_lnd_Nx.fields:      'WCSF'       , 'CATCH'  , 'SFMC'         ,
                              'WCRZ'       , 'CATCH'  , 'RZMC'         ,
                              ::
```

The same variable can be written in more than one Collection (i.e., at different temporal and/or spatial resolutions), and there is no limit to the number of Collections that are defined for and written by a simulation.  

Gridded ("2d") output can be on a grid other than the "native" grid that is associated with the tile space used in the simulation, as long as the output grid is defined in `HISTORY.rc`.  If the native is EASE grid, we dont have the "2d" option at this moment. If the native is Cubed-Sphere grid, an output gird can be defined for example as following:

```
GRID_LABELS: PC720x361-DC
   ::
PC720x361-DC.GRID_TYPE: LatLon
PC720x361-DC.IM_WORLD: 720
PC720x361-DC.JM_WORLD: 361
PC720x361-DC.POLE: PC
PC720x361-DC.DATELINE: DC
PC720x361-DC.LM: 1
```

and in `HISTORY.rc`, a line "VERSION: 1" should be added at the beginning. Finally, the definition of a collection can reference the output grid as:

```
tavg3_2d_lnd_Nx.grid_label:  PC720x361-DC,
```

MAPL HISTORY can generally write gridded ("2d") output in binary or netcdf-4 (nc4) format _**except**_ for GEOSldas simulations in EASE-grid tile space.  Output in tile space ("1d") must be written in binary format.

**Special considerations for GEOSldas**

1. Gridded ("2d") output _**cannot**_ be written for simulations in any EASE-grid tile space.  But since each EASE-grid cell contains at most one tile, it is straightforward to convert tile-space ("1d") output into gridded output (2d arrays) on the fly using the `i_indg` and `j_indg` data in the `./rc_out/*tilecoord*` file.  See Matlab readers and scripts in `./src/Application/LDAS_App/util`.

2. When running in EASE-grid tile space, `GEOSldas.x` always writes tile-space ("1d") output in binary format using MAPL HISTORY.  If nc4 output is requested in `HISTORY.rc`, the binary output is then automatically converted into nc4 format in post-processing as part of the `lenkf.j` job script using the `tile_bin2nc4.F90` utility. 

3. GEOSldas can bundle sub-daily nc4 output into daily nc4 files and write monthly-average output through the `POSTPROC_HIST` configuration option. 

The output from the MERRA-2 and GEOS-FP products is also written with MAPL HISTORY.  The corresponding "File Specification" documents, which are published as GMAO ["Office Notes"](https://gmao.gsfc.nasa.gov/pubs/), contain further documentation of many aspects of MAPL HISTORY.





