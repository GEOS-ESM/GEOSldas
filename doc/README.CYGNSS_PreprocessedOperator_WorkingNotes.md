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
- [x] Add a `cygnss_preprocessed_obs.F90` module-level cache for the coefficient product.
- [x] Read NetCDF product type `cygnss_tile_coefficient_preprocessor_netcdf` in the observation reader.
- [x] Validate schema version `0.3` or fail clearly in the observation reader.
- [x] Abort if the coefficient product grid tag inferred from filename/metadata does not match `tile_grid_d%gridtype`.
- [x] Convert zero-based support `tile_index0` to one-based full-domain tile numbers during support lookup.
- [x] Convert zero-based `sp_nearest_tile_index0` to one-based owner tile numbers.
- [x] Enforce one observation per owner tile for the first implementation.
- [x] Drop duplicate owner-tile observations deterministically by smallest `sp_nearest_tile_distance_km`.
- [x] Log read/kept/dropped observation counts.
- [x] Resolve daily CYGNSS L1 scalar products from `path/Yyyyy/Mmm/name-with-yyyymmdd`.
- [x] Read every daily product touched by the current assimilation window.
- [x] Load the matching daily coefficient-product set in the forward-operator cache.
- [x] Add a public LR reflectivity helper in `mwRTM_routines.F90`.
- [x] Match the Python simulator first: SFMC clipped to MWRTM porosity, Mironov dielectric, LR Fresnel reflectivity, no roughness.
- [x] Request `sfmc_l` and `sfmc_lH` for `varname='cygl1scal'` in `get_obs_pred()`.
- [x] Request MWRTM `clay` and `poros` in the local-plus-halo tile arrays for `varname='cygl1scal'`.
- [x] Make support-tile lookup exact within the local-plus-halo tile arrays.
- [x] Evaluate `Hx_linear = sum_t C_t * R_t`.
- [x] Convert `Hx_db = 10 * log10(Hx_linear)` and store in `Obs_pred_l`.
- [x] Abort loudly or set nodata with a clear diagnostic if support tiles are missing from the halo during development.
- [x] Include `cygl1scal` in update type 12 and 13 soil-moisture observation selection.
- [x] Document that CYGNSS-only update type 12/13 uses the non-Tb soil-moisture state vector.
- [x] Build `GEOSlandassim_GridComp`.
- [ ] Install.
- [ ] Compare GEOSldas forecast/analysis `H(x)` against the Python simulator on the selected-50 product.

## Support Tile Lookup Note

The current `get_obs_pred()` path usually finds tiles by ellipse/FOV search and
stores local-plus-halo tile positions in `ind_tmp(:)`. The preprocessed CYGNSS
operator is different: its support tiles are explicit.

The first implementation uses `tile_coord_lH(:)%f_num` as the full-domain tile
number and scans it for each coefficient support tile. This is simple and exact
for the selected-50 proof of concept. If this becomes expensive in larger runs,
replace the scan with a cached tile-number-to-`lH` index map.

`mwRTM_param` itself remains local-only. For CYGNSS, `get_obs_pred()` bundles
the small static fields needed at support tiles, currently `clay` and `poros`,
through the existing local-plus-halo communication helper alongside `sfmc_lH`.

## Daily Product Layout

CYGNSS L1 scalar products now follow a CYGNSS-L3-like daily directory layout.
The namelist path is the collection root, not the exact day directory:

```fortran
obs_param_nml(56)%path = '/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1'
obs_param_nml(56)%name = 'cygnss_l1_ddm3x5_crop_scalar_m36_yyyymmdd_cyg02.nc4'
```

For a 2019-11-01 run this resolves to:

```text
/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1/Y2019/M11/cygnss_l1_ddm3x5_crop_scalar_m36_20191101_cyg02.nc4
```

Only the explicit `yyyymmdd` token is replaced. Do not use a generic `dd`
token here because the filename itself contains `ddm3x5`.

The reader scans every daily file touched by the EnKF assimilation window:

```text
(date_time - dtstep_assim/2, date_time + dtstep_assim/2]
```

This means a window centered near midnight can read both adjacent daily files.
Each observation still carries its own timestamp, and the existing half-open
time-window check decides whether it belongs in the current cycle. Duplicate
owner tiles across files are handled the same way as duplicates within one file:
keep the observation with the smallest `sp_nearest_tile_distance_km`.

The forward-operator coefficient cache uses the same date-window file list as
the reader. When multiple daily products are loaded, their ragged support arrays
are concatenated and `tile_start` is offset into the combined support array.
This keeps `Observations_l` and `H(x)` synchronized for windows crossing
midnight.

## Assimilation Update Type

CYGNSS L1 scalar observations use `varname='cygl1scal'`. They have their own
observation operator, but the increments should use the soil-moisture update
state vector.

The implementation includes `cygl1scal` in the soil-moisture observation
selection for both update type 12 and update type 13:

- `update_type=12`: joint 3d soil moisture/Tskin/ght(1) analysis plus 1d snow
  analysis. This is the right path for the current all-sensors setup because it
  already includes the snow-analysis branch.
- `update_type=13`: 3d soil moisture/Tskin/ght(1) analysis without the snow
  branch.

For CYGNSS-only soil-moisture updates, no Tb observations are present, so the
state vector is:

- mineral tiles: `srfexc`, `rzexc`
- PEATCLSM tiles: `srfexc`, `rzexc`, `catdef`

Temperature states `tc1`, `tc2`, `tc4`, and `ght(1)` are added only when Tb
observations are also selected in the same local update.

## Runtime Behavior Learned

CYGNSS L1 scalar observations inherit the generic model-based satellite SFMC
QC through `qc_model_based_for_sat_sfmc()`. This is expected behavior, not a
reader loss, when the number of CYGNSS observations in `ldas_ObsFcstAna` is
smaller than the number of observations read from the preprocessed product.

The QC can set `Obs_pred` to nodata when the model state is unsuitable for a
soil-moisture-like update, including precipitation, snow, frozen/cold surface
conditions, or missing model/static inputs needed by the forward operator.
For the 2019-11-01 selected CYGNSS test, Discover diagnostics showed:

```text
09z: total=4  valid=1  nodata=3
12z: total=15 valid=6  nodata=9
15z: total=10 valid=7  nodata=3
```

The 12z `15 -> 6` behavior was explained by model-QC nodata in the forecast
operator, mostly SFMC/QC rejection plus one clay/porosity support issue. If this
comes up again, first check the log for read/kept counts, then check
`ldas_ObsFcstAna`; do not assume observations vanished in MPI output or tile
mapping.

The 2019-11-01/2019-11-02 two-day innovation-only test validated daily file
selection. Species 56 was configured with `ASSIM=F`, `GETINNOV=T`, and the
daily name template:

```text
cygnss_l1_ddm3x5_crop_scalar_m36_yyyymmdd_cyg02.nc4
```

At the cycle crossing midnight, the reader loaded both adjacent daily files:

```text
Reading .../CYGNSS_L1/Y2019/M11/cygnss_l1_ddm3x5_crop_scalar_m36_20191101_cyg02.nc4
Reading .../CYGNSS_L1/Y2019/M11/cygnss_l1_ddm3x5_crop_scalar_m36_20191102_cyg02.nc4
CYGNSS preprocessed obs daily files read: 2
CYGNSS preprocessed obs read: 100
```

CYGNSS rows with valid forecast and analysis `H(x)` appeared in `ObsFcstAna` at:

```text
20191101_0900z   2
20191101_1200z  10
20191101_1500z  11
20191102_0900z   6
20191102_1200z  10
20191102_1500z   9
```

This run still had `obsparam_assim=0` for all species and
`GEOSldas_update_type=12`, so it validates innovation diagnostics and daily file
handling, not CYGNSS assimilation increments.

The recurring Discover test workflow is:

```bash
cd /discover/nobackup/projects/land_da/GEOSldas_develop/GEOSldas/src/Components/@GEOSldas_GridComp
git pull --ff-only

cd /discover/nobackup/projects/land_da/GEOSldas_develop/GEOSldas/build
make -j6 install

cd ../install/bin
rm -rf /discover/nobackup/projects/land_da/all_sensors/LS_DAv8_M36_as_20191101_cygl1

./ldas_setup setup \
  /discover/nobackup/projects/land_da/all_sensors/ \
  /discover/nobackup/projects/land_da/all_sensors/LS_DAv8_M36_all_sensors_20191101_cygl1.txt \
  /discover/nobackup/projects/land_da/all_sensors/bat_inp_debug.txt

cd /discover/nobackup/projects/land_da/all_sensors/LS_DAv8_M36_as_20191101_cygl1/run
sbatch lenkf.j
```

## Reflectivity Helper Note

`mwRTM_routines.F90` already contains the relevant Mironov dielectric and
Fresnel reflectivity logic inside `mwRTM_get_Tb()`.

The new helper should return LR reflectivity directly, without running the full
tau-omega brightness-temperature model. Initial validation should prioritize
parity with the Python simulator over existing roughness or vegetation
conventions.

Implemented API:

```fortran
call mwRTM_get_lr_reflectivity( &
     freq, inc_angle, clay, poros, sfmc, reflectivity )
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

## Next Steps Toward Assimilation

To run CYGNSS in innovation-only mode, set the special namelist species 56 path
and file name, keep `assim = .false.`, and set `getinnov = .true.` with
`out_ObsFcstAna = 2` or `3`.

To actually assimilate CYGNSS, the first namelist switch is:

```fortran
obs_param_nml(56)%assim    = .true.
obs_param_nml(56)%getinnov = .true.
```

For the all-sensors experiment, use `update_type=12` so the existing snow branch
is preserved while CYGNSS enters the soil-moisture branch. Use `update_type=13`
only for a soil-moisture/Tb-only experiment without the 1d snow analysis.

Before using assimilation for science, confirm:

- `obs_param_nml(56)%errstd` is defensible. The initial value of `3.0 dB` is a
  placeholder.
- `obs_param_nml(56)%xcorr` and `%ycorr` are at least the effective FOV radius
  checks expected by `read_ens_upd_inputs()`. For `FOV = 100 km`, the old
  `0.25 deg` setting is too small.
- `bias_Npar`, `scale`, and `zeromean` are a deliberate choice. The first
  assimilation test should probably remain unscaled and unbiased unless a
  calibration strategy is ready.
- The test run still writes `ObsFcstAna` so forecast and analysis `H(x)` can be
  checked before trusting increments.
- Multi-file daily selection has been tested on Discover before moving from
  innovation-only diagnostics to assimilation.

## Work Log

### 2026-05-21

- Validated the selected-50 2019-11-01 innovation-only run well enough to move on.
- Decided to use a daily CYGNSS-L3-style product layout:
  `CYGNSS_L1/Yyyyy/Mmm/cygnss_l1_ddm3x5_crop_scalar_m36_yyyymmdd_cyg02.nc4`.
- Added date-window file selection for CYGNSS L1 scalar reads so cycles near
  midnight can use both adjacent daily files.
- Updated the coefficient-product cache to load the same date-window file set
  used by the reader.
- Verified `GEOSlandassim_GridComp` builds in `build-develop-20260520` after
  the daily-file changes.
- Optimized the coefficient-product cache so warm-cache calls do not reopen
  NetCDF files.
- Added `cygl1scal` to update type 13, then to update type 12 for the
  all-sensors setup.
- Clarified the update type 12 and 13 comments: CYGNSS-only updates use the
  soil-moisture state vector (`srfexc`, `rzexc`, plus `catdef` on PEATCLSM
  tiles), while temperature states are added only when Tb observations are
  selected.
- Validated a two-day innovation-only run using 50 observations per day. The
  midnight cycle read both daily files and `ObsFcstAna` contained valid CYGNSS
  forecast/analysis equivalents.

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
- Added a `get_obs_pred()` `cygl1scal` resource request for `sfmc_l`/`sfmc_lH`.
- Verified `GEOSlandassim_GridComp` builds in `build-develop-20260520`.
- Added `cygnss_preprocessed_obs.F90` to cache the coefficient product and evaluate exact support-tile sums.
- Added an MWRTM LR reflectivity helper using clipped SFMC, Mironov dielectric, and the Python operator LR Fresnel formula.
- Extended the observation halo communication helper to carry MWRTM `clay` and `poros` for CYGNSS support tiles.
- Hooked `cygl1scal` observations into `get_obs_pred()` so `Obs_pred_l` receives coefficient-weighted `H(x)` in dB.
- Re-verified `GEOSlandassim_GridComp` builds in `build-develop-20260520`; install and Discover-side runtime validation remain open.
