# H SAF ASCAT SSM Reader (H121/H139) — Implementation Notes

## Overview

New reader for the H SAF ASCAT Soil Surface Moisture (SSM) products on the 12.5 km
Fibonacci grid (H121 CDR v8 and H139 ICDR). These replace / supplement the existing
EUMET BUFR reader (`read_obs_sm_ASCAT_EUMET`), which reads NRT BUFR data on a coarser
25 km sampling.

| Product | Type | Period | Satellites |
|---------|------|--------|-----------|
| H121    | CDR v8 (reprocessed) | 2007-01-01 – 2024-12-31 | Metop-A, -B, -C |
| H139    | ICDR (extension)     | 2025-01-01 – ongoing    | Metop-B, -C |

Both products share the same netCDF schema and filename structure. The reader handles
both transparently via the same flist mechanism.

---

## H121/H139 file format

**Filename pattern:**

```
W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP{A|B|C}-12.5km-H{121|139}_C_LIIB_{creation14}_{start14}_{end14}____.nc
```

- `creation14`: 14-char creation timestamp (`YYYYMMDDHHMMSS`) — **unpredictable value**, cannot be used for filename construction
- `start14` / `end14`: sensing window, fixed offset from end of filename string
- Files are per-satellite, per-hour (24 files/day/satellite for H121)
- ~233,000 observations per file (globally, all land)

**Key netCDF variables:**

| Variable | Type | Description |
|----------|------|-------------|
| `soil_moisture` | `int16` (`_FillValue=-32767`) | SSM × 10000 (0–10000 → 0–100%) |
| `latitude` | `int32` | degrees × 1e6 |
| `longitude` | `int32` | degrees × 1e6 |
| `time` | `float64` | days since 1970-01-01 00:00:00 |
| `surface_flag` | `uint8` | bit 0: open water |
| `processing_flag` | `uint8` | bits 0–1: model/backscatter unusable |
| `snow_cover_probability` | `uint8` | % (0–100, 255=fill) |
| `frozen_soil_probability` | `uint8` | % (0–100, 255=fill) |
| `wetland_fraction` | `uint8` | % (0–100, 255=fill) |
| `topographic_complexity` | `uint8` | % (0–100, 255=fill) |
| `subsurface_scattering_probability` | `uint8` | % (0–100, 255=fill) |

`uint8` variables are read into Fortran `integer(1)`. Fill value 255 maps to -1 (signed),
which is rejected by all `>= threshold` QC checks — no special fill handling needed.

---

## GEOSldas integration

### New obs species

Three species added (indices 56, 57, 58 in `LDASsa_DEFAULT_inputs_ensupd.nml`):

| Index | Name | Satellite | Period |
|-------|------|-----------|--------|
| 56 | `ASCAT_HSAF_META_SM` | Metop-A | 2007-01-01 – 2021-11-15 |
| 57 | `ASCAT_HSAF_METB_SM` | Metop-B | 2013-06-01 – present |
| 58 | `ASCAT_HSAF_METC_SM` | Metop-C | 2019-04-01 – present |

All use `varname = 'sfds'` (degree of saturation), `units = '%'`, `FOV = 12.5 km`.
The reader dispatches on species name to select which satellite's files to open and
to enforce the satellite operating-period date filter.

### Files changed

| File | Change |
|------|--------|
| `clsm_ensupd_glob_param.F90` | `N_obs_species_nml` bumped 55 → 58 |
| `clsm_ensupd_read_obs.F90` | New subroutine `read_obs_sm_ASCAT_HSAF`; dispatch `case` added |
| `LDASsa_DEFAULT_inputs_ensupd.nml` | Species 56–58 added |

**No changes needed to:**
- `get_obs_pred` — already handles `case ('sfmc','sfds')`
- `read_ens_upd_inputs` — already maps `varname='sfds'` → `fcstvarname='sfmc'`
- Scaling: uses existing `scale_obs_sfmc_zscore` (same as EUMET ASCAT)

### N_obs_species_nml cascade

Changing this parameter causes a full recompile of everything that `use`s
`clsm_ensupd_glob_param`, which is a large fraction of the ensemble update routines.
Expect a long first build after this change.

---

## Reader design

**Subroutine:** `read_obs_sm_ASCAT_HSAF` in `clsm_ensupd_read_obs.F90`

Signature is identical to `read_obs_sm_ASCAT_EUMET` — same arguments in/out.

### File discovery

Uses the standard GEOSldas flist mechanism:
- Flist path: `{flistpath}/Y{YYYY}/M{MM}/D{DD}/{flistname}`
- Each line in the flist is a full filename (without directory path)
- The reader prepends `{path}/Y{YYYY}/M{MM}/` to each flist entry

Sensing-start time is extracted from the filename at a fixed offset from the end:
```fortran
str_start_time = tmpfname(n_fn-35:n_fn-22)   ! 14-char YYYYMMDDHHMMSS
```
This works because the creation timestamp (14 chars) and end timestamp (14 chars)
plus separators and suffix have a known length from the end, even though the creation
timestamp value itself is unpredictable.

### Time conversion

`time` variable is days since Unix epoch (1970-01-01):
```fortran
real*8, parameter :: unix_to_J2000_s = 10957.5d0 * 86400.0d0 + 64.184d0
obs_j2000 = time_raw(ii) * 86400.0d0 - unix_to_J2000_s
```

### QC logic

Applied per observation in order:

1. `ssm < 0 or ssm > 10000` → reject (invalid/fill)
2. `obs_j2000` outside assimilation window → reject
3. `surface_flag` bit 0 set → reject (open water)
4. `processing_flag` bits 0–1 set → reject (model or backscatter not usable)
5. `wetland_fraction >= 10%` → reject
8. `topographic_complexity >= 10%` → reject
9. `subsurface_scattering_probability >= 10%` → reject

Snow and frozen soil are **not** screened here — `qc_model_based_for_sat_sfmc`
handles those using model state (SWE, surface temperature).

Thresholds as parameters:
```fortran
real, parameter :: thr_wetland=10., thr_topo=10., thr_subsfc=10.
```

### Memory approach

Processes files one at a time (allocate/free per file), unlike the EUMET reader
which loads all obs into a large `tmp_data(max_obs, 14)` array. H121/H139 has
~233k obs/file × up to 30 files/day → per-file processing avoids large allocations.

Tile super-obs are accumulated across files into `ASCAT_sm`, `ASCAT_sm_std`,
`ASCAT_lon`, `ASCAT_lat`, `ASCAT_time` arrays (dimensioned by `N_catd`).

---

## Discover data paths

Paths set in `LDASsa_DEFAULT_inputs_ensupd.nml`:

| Species | path | flistpath | flistname |
|---------|------|-----------|-----------|
| META | `/discover/nobackup/amfox/ASCAT_HSAF/H121/metop_a/` | `/discover/nobackup/amfox/ASCAT_fname_lists/H121/metop_a/` | `H121_METOPA.txt` |
| METB | `/discover/nobackup/amfox/ASCAT_HSAF/H121_H139/metop_b/` | `/discover/nobackup/amfox/ASCAT_fname_lists/H121_H139/metop_b/` | `H121_H139_METOPB.txt` |
| METC | `/discover/nobackup/amfox/ASCAT_HSAF/H121_H139/metop_c/` | `/discover/nobackup/amfox/ASCAT_fname_lists/H121_H139/metop_c/` | `H121_H139_METOPC.txt` |

Data source: H121/H139 downloaded from H SAF FTP. H29 (NRT, 3-minute granules,
same schema) is **not** targeted by this reader.

---

## Pending / TODO

- [ ] Flist generation script for Discover: populate `Y{YYYY}/M{MM}/D{DD}/` trees
      from the staged netCDF files
- [ ] Stage H121 data at Discover paths above (download via H SAF FTP)
- [ ] Stage H139 data for METB/C (2025-onward)
- [ ] End-to-end test run with `assim = .false., getinnov = .true.` to check O-F stats
- [ ] Tune QC thresholds (especially `thr_subsfc`) based on initial O-F diagnostics
- [ ] Consider whether `thr_wetland = 10%` is too aggressive over high-wetland tiles

---

## Branch

`feature/amfox/ascat-hsaf-v8` (branched from `develop`)
