# src/Components Fortran Dev Memory

This file is durable project memory for future coding sessions in this repo.
It focuses on `src/Components` with bias toward Fortran development.

## Snapshot

- Captured: `2026-06-14 16:39:07 MDT`
- Fixture commit: `363a42b`
- `@GEOSldas_GridComp` commit: `1a39324`
- `@GEOSgcm_GridComp` commit: `ff111ad8`

## Quick Facts

- Top-level component trees:
  - `src/Components/@GEOSldas_GridComp`
  - `src/Components/@GEOSgcm_GridComp`
- Non-`.git` files under `src/Components`: `441`
- Fortran source files (`.F90/.f90/.F/.f/.F95/.F03/.F08`): `225`
- `CMakeLists.txt` files: `38`

## Topology

`src/Components/CMakeLists.txt` includes:

- `esma_add_subdirectory (GEOSldas_GridComp)`
- `esma_add_subdirectory (GEOSgcm_GridComp)`

### `@GEOSldas_GridComp` structure

- `GEOS_LdasGridComp.F90` is the top composite LDAS grid component.
- Subcomponents wired in `src/Components/@GEOSldas_GridComp/CMakeLists.txt`:
  - `GEOSmetforce_GridComp`
  - `GEOSlandpert_GridComp`
  - `GEOSens_GridComp`
  - `GEOSlandassim_GridComp`
  - `LDAS_Shared` (as subdir)
- App/executables in `GEOSldas_App`:
  - `GEOSldas.x`
  - `preprocess_ldas.x`
  - `mwrtm_bin2nc4.x`
  - `ascat_mask_maker.x` (under `util/inputs/ASCAT_sm_mask`)

Subdir file counts (`@GEOSldas_GridComp`, excluding `.git`):

- `GEOSldas_App`: `90` files, `5` Fortran
- `GEOSlandassim_GridComp`: `18` files, `16` Fortran
- `LDAS_Shared`: `16` files, `15` Fortran
- `GEOSlandpert_GridComp`: `9` files, `8` Fortran
- `GEOSmetforce_GridComp`: `5` files, `4` Fortran
- `GEOSens_GridComp`: `2` files, `1` Fortran

### `@GEOSgcm_GridComp` structure in this checkout

Important: this checkout is surface/land-focused and partial relative to full GEOSgcm.

- Present:
  - `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/...`
- Missing (directory absent):
  - `GEOSdataatm_GridComp`
  - `GEOSmkiau_GridComp`
  - `GEOSogcm_GridComp`
  - `GEOSwgcm_GridComp`
- Also missing in-tree aggregate source files:
  - `GEOS_GcmGridComp.F90`
  - `GEOS_AgcmGridComp.F90`
  - `GEOS_AgcmSimpleGridComp.F90`
  - `GEOS_PhysicsGridComp.F90`
  - `GEOS_SurfaceGridComp.F90`

Surface subtree counts:

- `GEOSsurface_GridComp`: `161` Fortran files (recursive)
- `GEOSland_GridComp`: `102` Fortran files (recursive)
- `Utils/Raster/makebcs`: `18` Fortran files
- `Utils/mk_restarts`: `21` Fortran files (includes `obsolete/`)

## CMake Target Map

### Core libraries

LDAS side:

- `GEOSldas_GridComp` from `GEOS_LdasGridComp.F90`
- `GEOS_LdasShared` from `LDAS_Shared/*.F90` + legacy `.f` sources
- `GEOSmetforce_GridComp`
- `GEOSlandpert_GridComp`
- `GEOSens_GridComp`
- `GEOSlandassim_GridComp`
- `GEOSexportcatchincr_GridComp`

GCM/surface side (present subtree):

- `GEOSland_GridComp`
- `GEOSlandice_GridComp`
- `GEOScatch_GridComp`
- `GEOScatchCN_GridComp`
- `GEOS_CatchCNShared`
- `GEOS_CatchCNCLM40_GridComp` + `CLM40`
- `GEOS_CatchCNCLM45_GridComp` + `CLM45`
- `GEOSroute_GridComp`
- `GEOSigni_GridComp`
- `GEOSvegdyn_GridComp`
- `GEOS_LandShared`
- `GEOS_SurfaceShared`
- `makebcs`
- raster routing helper library (from `routing_constant.f90`)
- mk_restarts helper library

### Executables

LDAS:

- `GEOSldas.x`
- `preprocess_ldas.x`
- `mwrtm_bin2nc4.x`
- `ascat_mask_maker.x`

Surface utilities:

- `makebcs` utilities:
  - `CombineRasters.x`
  - `mkCatchParam.x`
  - `mkCubeFVRaster.x`
  - `mkLandRaster.x`
  - `mkLatLonRaster.x`
  - `mkMITAquaRaster.x`
  - `mkMOMAquaRaster.x`
  - `FillMomGrid.x`
  - `mk_runofftbl.x`
  - `mkEASETilesParam.x`
  - `TileFile_ASCII_to_nc4.x`
- Raster preproc:
  - `loss_during_day.x`
  - `loss_surf_5cm_gensoil.x`
  - Routing `exe_srcs` converted to `.x`:
    - `get_finalID_msk.x`
    - `get_landocean_Greenland_real.x`
    - `get_outlets_land_allcat.x`
    - `get_sinkxy_land.x`
    - `get_outlets_catchindex.x`
    - `get_outlets_land.x`
    - `Pfaf_to_2d_30s_land.x`
  - Topography utils:
    - `generate_scrip_cube_topo.x`
    - `convert_bin_to_netcdf_topo.x`
    - `convert_to_gmao_output_topo.x`
- `mk_restarts` `exe_srcs` converted to `.x`:
  - `Scale_Catch.x`
  - `Scale_CatchCN.x`
  - `cv_SaltRestart.x`
  - `SaltIntSplitter.x`
  - `SaltImpConverter.x`
  - `mk_CICERestart.x`
  - `mk_CatchCNRestarts.x`
  - `mk_CatchRestarts.x`
  - `mk_LakeLandiceSaltRestarts.x`
  - `mk_RouteRestarts.x`
  - `mk_GEOSldasRestarts.x`
  - `mk_catchANDcnRestarts.x`

## LDAS Module Inventory (from `module` declarations)

- `GEOS_LdasGridCompMod`
- `GEOS_MetforceGridCompMod`
- `GEOS_LandPertGridCompMod`
- `GEOS_EnsGridCompMod`
- `GEOS_LandAssimGridCompMod`
- `GEOS_ExportCatchIncrGridCompMod`
- `LDAS_ConvertMod`
- `LDAS_DriverTypes`
- `LDAS_ExceptionsMod`
- `LDAS_ForceMod`
- `LDAS_HashTable`
- `LDAS_InterpMod`
- `LDAS_PertRoutinesMod`
- `LDAS_PertTypes`
- `LDAS_TileCoordRoutines`
- `LDAS_TileCoordType`
- `LDAS_ensdrv_Globals`
- `LDAS_ensdrv_functions`
- `LDAS_ensdrv_mpi`
- `RepairForcingMod`
- `enkf_types`
- `catch_types`
- `enkf_general`
- `adapt_types`
- `mwRTM_types`
- `mwRTM_routines`
- `clsm_ensupd_glob_param`
- `clsm_ensupd_read_obs`
- `clsm_ensupd_upd_routines`
- `clsm_ensupd_enkf_update`
- `clsm_ensdrv_drv_routines`
- `clsm_ensdrv_out_routines`
- `clsm_bias_routines`
- `clsm_adapt_routines`
- `catch_bias_types`
- `land_pert_routines`
- `force_and_cat_progn_pert_types`
- `my_matrix_functions`
- `io_hdf5`
- `preprocess_ldas_routines`
- `nr_fft`
- `nr_jacobi`
- `nr_ran2_gasdev`
- `random_fieldsMod`
- `StringRandom_fieldsMapMod`

## Coupling Notes

`GEOS_LdasGridComp.F90` directly `use`s:

- `GEOS_MetforceGridCompMod`
- `GEOS_LandGridCompMod`
- `GEOS_LandPertGridCompMod`
- `GEOS_EnsGridCompMod`
- `GEOS_LandAssimGridCompMod`
- `GEOS_LandiceGridCompMod`
- `LDAS_*` shared modules
- `catch_constants`, `StieglitzSnow`, `SurfParams`

CMake-level LDAS dependencies include:

- `GEOSland_GridComp`
- `GEOSlandice_GridComp`
- `makebcs`
- `MAPL`
- `HDF5`, `NetCDF` (via `GEOSlandassim_GridComp`)

This means many LDAS edits have coupling risk with surface shared modules and `makebcs` utilities.

## Practical Development Notes

- These component trees are nested git repos (`.git` exists under both `@GEOSldas_GridComp` and `@GEOSgcm_GridComp`).
- Local env script for Fortran build tasks:
  - `scripts/dev_env.sh`
  - Exports `FC`, `LDFLAGS`, `CPPFLAGS`, and `GEOSLDAS_MPI_VERSION_STRING`.
- VS Code tasks (see `.vscode/tasks.json`) include:
  - release configure/build
  - debug build
  - format current file
  - syntax-check current file
  - smoke target build
- If configuring with local OpenMPI and CMake tests fail due empty MPI version strings, pass:
  - `-DMPI_Fortran_LIBRARY_VERSION_STRING=Unknown`

Companion generated inventories in `doc/`:

- `doc/src_components_file_inventory.txt` (all non-`.git` files under `src/Components`)
- `doc/src_components_fortran_inventory.txt` (all Fortran files under `src/Components`)
- `doc/src_components_cmake_targets_inventory.txt` (CMake target/subdirectory statements)
- Refresh all of the above plus snapshot/count lines with:
  - `scripts/refresh_components_memory.sh`

## Appendix: Full CMakeLists Index (src/Components)

```text
src/Components/@GEOSgcm_GridComp/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSland_GridComp/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSland_GridComp/GEOScatchCN_GridComp/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSland_GridComp/GEOScatchCN_GridComp/GEOScatchCNCLM40_GridComp/CLM40/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSland_GridComp/GEOScatchCN_GridComp/GEOScatchCNCLM40_GridComp/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSland_GridComp/GEOScatchCN_GridComp/GEOScatchCNCLM45_GridComp/CLM45/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSland_GridComp/GEOScatchCN_GridComp/GEOScatchCNCLM45_GridComp/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSland_GridComp/GEOScatchCN_GridComp/Shared/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSland_GridComp/GEOScatch_GridComp/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSland_GridComp/GEOSigni_GridComp/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSland_GridComp/GEOSroute_GridComp/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSland_GridComp/GEOSvegdyn_GridComp/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSland_GridComp/Shared/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSlandice_GridComp/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Shared/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Utils/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Utils/Raster/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Utils/Raster/makebcs/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Utils/Raster/preproc/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Utils/Raster/preproc/routing/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Utils/Raster/preproc/soil/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Utils/Raster/preproc/topography/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Utils/Raster/preproc/topography/utils_topo/CMakeLists.txt
src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Utils/mk_restarts/CMakeLists.txt
src/Components/@GEOSldas_GridComp/CMakeLists.txt
src/Components/@GEOSldas_GridComp/GEOSens_GridComp/CMakeLists.txt
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/CMakeLists.txt
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/GEOSexportcatchincr_GridComp/CMakeLists.txt
src/Components/@GEOSldas_GridComp/GEOSlandpert_GridComp/CMakeLists.txt
src/Components/@GEOSldas_GridComp/GEOSldas_App/CMakeLists.txt
src/Components/@GEOSldas_GridComp/GEOSldas_App/util/inputs/ASCAT_sm_mask/CMakeLists.txt
src/Components/@GEOSldas_GridComp/GEOSmetforce_GridComp/CMakeLists.txt
src/Components/@GEOSldas_GridComp/LDAS_Shared/CMakeLists.txt
src/Components/CMakeLists.txt
```

## Appendix: Full LDAS Fortran File Index

```text
src/Components/@GEOSldas_GridComp/GEOS_LdasGridComp.F90
src/Components/@GEOSldas_GridComp/GEOSens_GridComp/GEOS_EnsGridComp.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/GEOS_LandAssimGridComp.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/GEOSexportcatchincr_GridComp/GEOS_ExportCatchIncrGridComp.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/adapt_types.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/catch_bias_types.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_adapt_routines.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_bias_routines.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensdrv_drv_routines.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensdrv_out_routines.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_enkf_update.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_glob_param.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_read_obs.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_upd_routines.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/enkf_general.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/io_hdf5.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/mwRTM_routines.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/mwRTM_types.F90
src/Components/@GEOSldas_GridComp/GEOSlandpert_GridComp/GEOS_LandPertGridComp.F90
src/Components/@GEOSldas_GridComp/GEOSlandpert_GridComp/LDAS_PertRoutines.F90
src/Components/@GEOSldas_GridComp/GEOSlandpert_GridComp/force_and_cat_progn_pert_types.F90
src/Components/@GEOSldas_GridComp/GEOSlandpert_GridComp/land_pert.F90
src/Components/@GEOSldas_GridComp/GEOSlandpert_GridComp/nr_fft.F90
src/Components/@GEOSldas_GridComp/GEOSlandpert_GridComp/nr_jacobi.F90
src/Components/@GEOSldas_GridComp/GEOSlandpert_GridComp/nr_ran2_gasdev.F90
src/Components/@GEOSldas_GridComp/GEOSlandpert_GridComp/random_fields.F90
src/Components/@GEOSldas_GridComp/GEOSldas_App/GEOSldas.F90
src/Components/@GEOSldas_GridComp/GEOSldas_App/preprocess_ldas.F90
src/Components/@GEOSldas_GridComp/GEOSldas_App/preprocess_ldas_routines.F90
src/Components/@GEOSldas_GridComp/GEOSldas_App/util/inputs/ASCAT_sm_mask/ascat_mask_maker.F90
src/Components/@GEOSldas_GridComp/GEOSldas_App/util/inputs/mwRTM_params/mwrtm_bin2nc4.F90
src/Components/@GEOSldas_GridComp/GEOSmetforce_GridComp/GEOS_MetforceGridComp.F90
src/Components/@GEOSldas_GridComp/GEOSmetforce_GridComp/LDAS_Forcing.F90
src/Components/@GEOSldas_GridComp/GEOSmetforce_GridComp/LDAS_HashTable.F90
src/Components/@GEOSldas_GridComp/GEOSmetforce_GridComp/LDAS_Interp.F90
src/Components/@GEOSldas_GridComp/LDAS_Shared/LDAS_Convert.F90
src/Components/@GEOSldas_GridComp/LDAS_Shared/LDAS_DriverTypes.F90
src/Components/@GEOSldas_GridComp/LDAS_Shared/LDAS_Exceptions.F90
src/Components/@GEOSldas_GridComp/LDAS_Shared/LDAS_PertTypes.F90
src/Components/@GEOSldas_GridComp/LDAS_Shared/LDAS_RepairForcing.F90
src/Components/@GEOSldas_GridComp/LDAS_Shared/LDAS_TileCoordRoutines.F90
src/Components/@GEOSldas_GridComp/LDAS_Shared/LDAS_TileCoordType.F90
src/Components/@GEOSldas_GridComp/LDAS_Shared/LDAS_ensdrv_Globals.F90
src/Components/@GEOSldas_GridComp/LDAS_Shared/LDAS_ensdrv_functions.F90
src/Components/@GEOSldas_GridComp/LDAS_Shared/LDAS_ensdrv_mpi.F90
src/Components/@GEOSldas_GridComp/LDAS_Shared/catch_types.F90
src/Components/@GEOSldas_GridComp/LDAS_Shared/enkf_types.F90
src/Components/@GEOSldas_GridComp/LDAS_Shared/my_lu_decomp.f
src/Components/@GEOSldas_GridComp/LDAS_Shared/my_matrix_functions.F90
src/Components/@GEOSldas_GridComp/LDAS_Shared/nr_indexx.f
```

## Appendix: GEOSgcm Fortran Density (present subtree)

```text
102 src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSland_GridComp
 55 src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Utils
  3 src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Shared
  1 src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSlandice_GridComp
```

## Fast Commands To Rebuild This Map

```bash
# Counts
find src/Components -path '*/.git' -prune -o -path '*/.git/*' -prune -o -type f -print | wc -l
find src/Components -path '*/.git' -prune -o -path '*/.git/*' -prune -o -type f -print | rg -i '\.(F90|f90|F|f|F95|f95|F03|f03|F08|f08)$' | wc -l
find src/Components -path '*/.git' -prune -o -path '*/.git/*' -prune -o -name CMakeLists.txt -print | wc -l

# CMake target statements
find src/Components -name CMakeLists.txt -print0 | xargs -0 awk '/esma_add_library[[:space:]]*\(|ecbuild_add_executable[[:space:]]*\(|add_executable[[:space:]]*\(|add_library[[:space:]]*\(|esma_add_subdirectory[[:space:]]*\(|esma_add_subdirectories[[:space:]]*\(|add_subdirectory[[:space:]]*\(/ {print FILENAME ":" FNR ":" $0}'

# Fortran modules
rg -n -i '^[[:space:]]*module[[:space:]]+[A-Za-z][A-Za-z0-9_]*' src/Components -g '*.[Ff]*' | rg -vi 'end[[:space:]]+module|module[[:space:]]+procedure'
```
