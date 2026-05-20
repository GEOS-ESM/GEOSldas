# CYGNSS Preprocessed Operator Working Notes

This is the GEOSldas-side working document for integrating the preprocessed
CYGNSS coefficient operator. It links back to the upstream handoff note and
tracks decisions, implementation steps, validation notes, and unresolved
questions as the work evolves.

## Source Handoff

Primary handoff note:

```text
/Users/amfox/Desktop/CYGNSS_operator/docs/geosldas_cygnss_preprocessed_operator_handoff.md
```

Python simulator reference:

```text
/Users/amfox/Desktop/CYGNSS_operator/scripts/run_geosldas_cygnss_preprocessed_simulator.py
```

Current test coefficient product:

```text
/Users/amfox/Desktop/CYGNSS_operator/artifacts/out_images/cygnss_coeff_preprocessor_m36_selected50/cygnss_coefficients_m36_20191028_selected50.nc4
```

Current simulator output:

```text
/Users/amfox/Desktop/CYGNSS_operator/artifacts/out_images/cygnss_preprocessed_geosldas_simulator_m36_20191028_selected50/cygnss_preprocessed_sim_obs_summary.csv
```

## Current Build Context

Working GEOSldas checkout:

```text
/Users/amfox/Desktop/GEOSldas_develop/GEOSldas
```

Current build/install tree used during this work:

```text
build-develop-20260520
install-develop-20260520
```

The Mac build notes for this workspace are maintained in:

```text
doc/README.BuildNotes.Mac.md
```

## Accepted Initial Design

Use the low-invasive path from the handoff:

- evaluate CYGNSS `H(x)` inside `get_obs_pred()` for both forecast and analysis
- keep the normal GEOSldas observation shell in `obs_type`
- keep ragged CYGNSS coefficient support in a separate module-level cache
- avoid changing `obs_type`
- avoid changing the MPI derived-type plumbing in `LDAS_ensdrv_mpi.F90`
- use one preprocessed CYGNSS observation per owner tile for the first proof of concept
- use `owner_tilenum = sp_nearest_tile_index0 + 1`
- use dB for both `Observations_l%obs` and `Obs_pred_l`
- use `FOV = 100 km` initially so existing halo communication is likely to include support tiles
- set `assim = .false.` until forecast/analysis diagnostics validate cleanly

The distinct preprocessed operator should not be conflated with existing
CYGNSS Level 3 soil-moisture species such as `CYGNSS_SM_6hr` and
`CYGNSS_SM_daily`.

## Existing GEOSldas Hooks

Forecast observation equivalents are computed before the EnKF update in:

```text
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_enkf_update.F90
```

Relevant call:

```fortran
call get_obs_pred(.true., ...)
```

Analysis observation equivalents are computed after the EnKF update in the same
file.

Relevant call:

```fortran
call get_obs_pred(.false., ...)
```

The main implementation target is:

```text
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_upd_routines.F90
```

Key places in that file:

- `read_ens_upd_inputs()` fills `obs_param%fcstvarname` and `obs_param%fcstunits`
- `get_obs_pred()` decides which model diagnostics are needed
- `get_obs_pred()` maps model tiles to each observation and fills `Obs_pred_l`

## Likely File Touches

Expected first implementation files:

```text
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/CMakeLists.txt
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_read_obs.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_upd_routines.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/mwRTM_routines.F90
src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/cygnss_preprocessed_obs.F90
src/Components/@GEOSldas_GridComp/GEOSldas_App/LDASsa_DEFAULT_inputs_ensupd.nml
```

Files to avoid unless the design changes:

```text
src/Components/@GEOSldas_GridComp/LDAS_Shared/enkf_types.F90
src/Components/@GEOSldas_GridComp/LDAS_Shared/LDAS_ensdrv_mpi.F90
```

## Implementation Checklist

- [x] Add a `CYGNSS_L1_DDM3X5_CROP_SCALAR` namelist species distinct from existing CYGNSS L3 soil moisture.
- [x] Allow `obs_param_nml(:)%varname = 'cygl1scal'` in `read_ens_upd_inputs()`.
- [x] Fill `fcstvarname = 'cygl1scal'` and `fcstunits = 'dB'`.
- [ ] Add a `cygnss_preprocessed_obs.F90` module-level cache for the coefficient product.
- [x] Read NetCDF product type `cygnss_tile_coefficient_preprocessor_netcdf` in the observation reader.
- [x] Validate schema version `0.3` or fail clearly in the observation reader.
- [x] Abort if the coefficient product grid tag inferred from filename/metadata does not match `tile_grid_d%gridtype`.
- [ ] Convert zero-based support `tile_index0` to one-based full-domain tile numbers in the coefficient cache.
- [x] Convert zero-based `sp_nearest_tile_index0` to one-based owner tile numbers.
- [x] Enforce one observation per owner tile for the first implementation.
- [x] Drop duplicate owner-tile observations deterministically by smallest `sp_nearest_tile_distance_km`.
- [x] Log read/kept/dropped observation counts.
- [ ] Add a public LR reflectivity helper in `mwRTM_routines.F90`.
- [ ] Match the Python simulator first: SFMC clipped to MWRTM porosity, Mironov dielectric, LR Fresnel reflectivity, no roughness.
- [x] Request `sfmc_l` and `sfmc_lH` for `varname='cygl1scal'` in `get_obs_pred()`.
- [ ] Make support-tile lookup exact within the local-plus-halo tile arrays.
- [ ] Evaluate `Hx_linear = sum_t C_t * R_t`.
- [ ] Convert `Hx_db = 10 * log10(Hx_linear)` and store in `Obs_pred_l`.
- [ ] Abort loudly or set nodata with a clear diagnostic if support tiles are missing from the halo during development.
- [ ] Build and install.
- [ ] Compare GEOSldas forecast/analysis `H(x)` against the Python simulator on the selected-50 product.

## Support Tile Lookup Note

The current `get_obs_pred()` path usually finds tiles by ellipse/FOV search and
stores local-plus-halo tile positions in `ind_tmp(:)`. The preprocessed CYGNSS
operator is different: its support tiles are explicit.

The likely clean approach is to maintain a full-domain tile-number vector for
the local-plus-halo arrays, for example:

```fortran
integer, allocatable :: tilenum_lH(:)
```

The first `N_catl` entries should map local tiles through `l2f` or equivalent
full-domain numbering; halo entries should be populated when remote
`tile_coord_lH` entries are appended. The CYGNSS evaluator can then map each
support tile number to a local-plus-halo index before reading `sfmc_lH`,
`mwRTM_param`, and other tile fields.

This is probably the first implementation risk to burn down.

## Reflectivity Helper Note

`mwRTM_routines.F90` already contains the relevant Mironov dielectric and
Fresnel reflectivity logic inside `mwRTM_get_Tb()`.

The new helper should return LR reflectivity directly, without running the full
tau-omega brightness-temperature model. Initial validation should prioritize
parity with the Python simulator over existing roughness or vegetation
conventions.

Conceptual API:

```fortran
call mwRTM_get_lr_reflectivity( &
     N_tile, freq, inc_angle, mwRTM_param, sfmc, reflectivity )
```

## Initial Namelist Sketch

Use a new species number when implementing. Values below are placeholders to
capture the initial design.

```fortran
obs_param_nml(?)%descr          = 'CYGNSS_L1_DDM3X5_CROP_SCALAR'
obs_param_nml(?)%orbit          = 3
obs_param_nml(?)%pol            = 0
obs_param_nml(?)%N_ang          = 0
obs_param_nml(?)%freq           = 1.57542e9
obs_param_nml(?)%FOV            = 100.
obs_param_nml(?)%FOV_units      = 'km'
obs_param_nml(?)%assim          = .false.
obs_param_nml(?)%scale          = .false.
obs_param_nml(?)%getinnov       = .false.
obs_param_nml(?)%RTM_ID         = 0
obs_param_nml(?)%nodata         = -9999.
obs_param_nml(?)%varname        = 'cygl1scal'
obs_param_nml(?)%units          = 'dB'
obs_param_nml(?)%fcstvarname    = 'NULL'
obs_param_nml(?)%fcstunits      = 'NULL'
obs_param_nml(?)%path           = '/path/to/preprocessed/netcdf/directory'
obs_param_nml(?)%name           = 'cygnss_coefficients_m36_20191028_selected50.nc4'
obs_param_nml(?)%errstd         = 3.
obs_param_nml(?)%std_normal_max = 2.5
obs_param_nml(?)%zeromean       = .true.
obs_param_nml(?)%coarsen_pert   = .true.
obs_param_nml(?)%xcorr          = 0.25
obs_param_nml(?)%ycorr          = 0.25
obs_param_nml(?)%adapt          = 0
```

## Open Questions

- Should first validation use exactly the Python simulator's SFMC-to-porosity convention, or the existing GEOSldas `catch2mwRTM_vars()` convention?
- Should missing support tiles abort in all development builds, or only when a debug flag is enabled?
- Where should the CYGNSS preprocessed cache lifecycle live so it is read once per relevant file/time and remains synchronized with `Observations_l` after sorting and compaction?
- How should production runs handle multiple CYGNSS observations per owner tile once we move past the selected-50 proof of concept?
- Which `update_type` path should include CYGNSS once `assim=.true.` is enabled?

## Work Log

### 2026-05-20

- Read the CYGNSS handoff note and Python simulator reference.
- Confirmed `get_obs_pred()` is the right forecast/analysis hook.
- Confirmed `read_ens_upd_inputs()` and `get_obs_pred()` both need a new `cygl1scal` varname case.
- Confirmed existing CYGNSS L3 soil-moisture support is separate and should not be reused as the preprocessed operator.
- Created this working note to track the GEOSldas-side implementation.
- Added disabled default namelist species 56, `CYGNSS_L1_DDM3X5_CROP_SCALAR`, with `varname='cygl1scal'`.
- Bumped `N_obs_species_nml` from 55 to 56.
- Added a CYGNSS L1 scalar observation reader that validates product type/schema, reads `observed_y_db`, assigns one observation per owner tile, filters to the assimilation window, and logs duplicate/drop counts.
- Added a model-grid guard so an M36 coefficient product cannot be read into a non-M36 tile grid.
- Added a `get_obs_pred()` `cygl1scal` resource request for `sfmc_l`/`sfmc_lH`; the actual coefficient-weighted forward operator is still the next implementation step.
- Verified `GEOSlandassim_GridComp` builds in `build-develop-20260520`.
