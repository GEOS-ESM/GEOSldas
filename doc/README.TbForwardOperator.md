
# GEOSldas Tb Forward Operator — Reference for Claude

## What it does

The Tb forward operator maps Catchment land surface model state to L-band microwave
brightness temperatures (h- and v-polarization), used as the observation operator
H(x) in the ensemble Kalman filter (EnKF) data assimilation of SMOS/SMAP observations.

The model is an **L-band tau-omega radiative transfer model** (RTM). Four configurations
exist, selected via `RTM_ID`:

| RTM_ID | Dielectric model | Atm correction   | Pol mixing Q     | Target sensor |
|--------|-----------------|-----------------|-----------------|---------------|
| 1      | Wang-Schmugge   | Pellarin        | Q = 0           | SMOS          |
| 2      | Wang-Schmugge   | none            | Q = 0           | SMAP          |
| 3      | Mironov         | Pellarin        | Q = 0.1771*h_mc | SMOS          |
| 4      | Mironov         | none            | Q = 0.1771*h_mc | SMAP L4_SM v8 |

---

## Key source files

| File | Role |
|------|------|
| `src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/mwRTM_routines.F90` | All RTM physics |
| `src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/mwRTM_types.F90`   | `mwRTM_param_type` definition |
| `src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_upd_routines.F90` | Calling context (EnKF loop) |

Public interface from `mwRTM_routines`:
- `mwRTM_get_Tb`        — main forward operator
- `catch2mwRTM_vars`    — converts Catchment model state to RTM inputs

---

## Parameters: `mwRTM_param_type`

Defined in `mwRTM_types.F90`. Read from a netCDF file at runtime (via
`get_mwrtm_param()` in `GEOS_LandAssimGridComp.F90`).

| Field | Description | Units |
|-------|-------------|-------|
| `vegcls`, `soilcls` | Land cover and soil class (may differ from Catchment params) | integer |
| `sand`, `clay` | Soil texture fractions | [0–1] |
| `poros` | Soil porosity (may differ from Catchment `poros`) | m3/m3 |
| `wang_wt`, `wang_wp` | Wang model transition SM and wilting point | m3/m3 |
| `rgh_hmin`, `rgh_hmax` | Min/max roughness parameter h | dim-less |
| `rgh_wmin`, `rgh_wmax` | SM thresholds for h transitions | m3/m3 |
| `rgh_Nrh`, `rgh_Nrv` | Incidence angle exponents for h/v pol roughness | dim-less |
| `rgh_polmix` | Polarization mixing parameter | dim-less |
| `omega` | Single scattering albedo | dim-less |
| `bh`, `bv` | Vegetation b parameter (tau = b*VWC); h and v pol | dim-less |
| `lewt` | VWC = lewt * LAI | kg/m2 per unit LAI |
| `vegopacity` | Alternative: time-varying opacity from file (replaces bh/bv/lewt) | dim-less |

Either `(bh, bv, lewt)` **or** `vegopacity` must be non-nodata for a tile; not both.

---

## Inputs to `mwRTM_get_Tb`

After `catch2mwRTM_vars` is called:

| Variable | Source | Notes |
|----------|--------|-------|
| `soilmoist` [m3/m3] | Catchment `sfmc` scaled by `poros_mwRTM/poros_catch` | RTM uses its own porosity |
| `soiltemp` [K] | Catchment layer-1 soil temp `tp1` (converted from °C) | Used for both soil and canopy temperature |
| `LAI` | From forcing/model | Only needed when using `(bh, bv, lewt)` |
| `SWE` [kg/m2] | Catchment snow state | Used for masking |
| `Tair` [K] | Near-surface air temperature | Only needed for Pellarin atm correction |
| `elev` [m] | Tile elevation above sea level | Only needed for Pellarin atm correction |
| `freq` [Hz] | Observation frequency | Set per obs type |
| `inc_angle` [deg] | Incidence angle | Set per obs type |
| `RTM_ID` | Integer 1–4 | Selects physics configuration |

---

## Physics walkthrough (`mwRTM_get_Tb`)

For each tile, Tb is only computed when ALL of:
- `SWE < 1e-4 kg/m2` (snow-free)
- `soiltemp > 273.35 K` (non-frozen, MAPL_TICE + 0.2)
- mwRTM parameters are not nodata

Otherwise `Tb_h = Tb_v = nodata_generic`.

### Step 1 — Soil dielectric constant

**RTM_ID 1, 2** → `DIELWANG` (Wang-Schmugge 1980):
- Calls `DIEL_WAT` for free water dielectric constant (Debye relaxation, Dobson/Ulaby)
- Mixes ice, air, rock, and water dielectric constants
- Piecewise linear in SM with transition at `wang_wt`
- Adds conductivity loss at low frequency

**RTM_ID 3, 4** → `MIRONOV` (Mironov et al. IEEE TGRS 2009):
- Distinguishes bound soil water (BSW) and free soil water (FSW) phases
- Clay-fraction-based regression coefficients for relaxation parameters
- Piecewise linear in SM with transition at `mvt` (clay-dependent)

Both return a complex dielectric constant `c_er`.

### Step 2 — Fresnel reflectivity (smooth surface)

```
roh = |( cos(inc) - sqrt(c_er - sin²(inc)) ) / ( cos(inc) + sqrt(c_er - sin²(inc)) )|²
rov = |( c_er*cos(inc) - sqrt(c_er - sin²(inc)) ) / ( c_er*cos(inc) + sqrt(c_er - sin²(inc)) )|²
```

### Step 3 — Roughness correction (h-Q model, Choudhury et al. 1979)

h parameter varies linearly with soil moisture between `rgh_hmin` and `rgh_hmax`:
- `sm <= rgh_wmin` → `h = rgh_hmax`
- `sm >= rgh_wmax` → `h = rgh_hmin`
- linear interpolation in between

Polarization mixing Q (only at L-band, i.e. `freq < 2 GHz`):
- RTM_ID 1, 2: Q = 0
- RTM_ID 3, 4: Q = 0.1771 * h_mc

Rough surface reflectivity:
- RTM_ID 1, 2: `rsh = ((1-Q)*roh + Q*rov) * exp(-h * cos^Nrh(inc))`
- RTM_ID 3, 4: `rsh = ((1-Q)*roh + Q*rov) * exp(-h * cos²(inc))`
(same for v-pol with rov and Nrv)

### Step 4 — Vegetation attenuation (tau-omega)

**Using static params (bh, bv, lewt)**:
```
VWC = lewt * LAI
exptauh = exp(-bh * VWC / cos(inc))
exptauv = exp(-bv * VWC / cos(inc))
```

**Using vegopacity from file**:
```
exptauh = exptauv = exp(-vegopacity)
```

Canopy temperature: `Tc = soiltemp`

Vegetation emission factor:
```
Ah = Tc * (1 - omega) * (1 - exptauh)
Av = Tc * (1 - omega) * (1 - exptauv)
```

### Step 5 — Top-of-vegetation Tb (Crow et al. 2005, Eq. 1)

```
Tb_h = soiltemp * (1 - rsh) * exptauh  +  Ah * (1 + rsh * exptauh)
Tb_v = soiltemp * (1 - rsv) * exptauv  +  Av * (1 + rsv * exptauv)
```

### Step 6 — Atmospheric correction (RTM_ID 1, 3 only; Pellarin model)

```
tau_atm = exp(-3.926 - 0.2211*Z[km] - 0.00369*Tair[K])
GOSSAT  = exp(-tau_atm / cos(inc))
TAEQ    = exp(4.927 + 0.002195 * Tair)
Tb_ad   = TAEQ*(1 - GOSSAT) + 2.7*GOSSAT    ! downwelling
Tb_au   = TAEQ*(1 - GOSSAT)                  ! upwelling

! include downwelling reflected off vegetation into top-of-vegetation Tb
Tb_h += Tb_ad * rsh * exptauh²
Tb_v += Tb_ad * rsv * exptauv²

! propagate through atmosphere to top-of-atmosphere
Tb_h = Tb_h * exp(-tau_atm/cos(inc)) + Tb_au
Tb_v = Tb_v * exp(-tau_atm/cos(inc)) + Tb_au
```

RTM_ID 2 and 4 skip this step (SMAP geometry makes atmospheric correction negligible).

---

## How it is called in the EnKF loop

In `clsm_ensupd_upd_routines.F90`, inside the per-ensemble-member loop:

```fortran
! 1. Convert Catchment state to RTM inputs
call catch2mwRTM_vars( N_catl, cat_param%vegcls, cat_param%poros, &
     mwRTM_param%poros, sfmc_l(:,n_e), tsurf_l(:,n_e), tp_l(1,:), &
     smoist, stemp_l(:,n_e) )

! 2. Loop over unique (freq, inc_angle, RTM_ID) combinations
do j = 1, N_TbuniqFreqAngRTMid
   freq      = Tb_freq_ang_RTMid(j,1)
   inc_angle = Tb_freq_ang_RTMid(j,2)
   RTM_id    = nint(Tb_freq_ang_RTMid(j,3))

   call mwRTM_get_Tb( N_catl, freq, inc_angle, mwRTM_param, &
        tile_coord_l%elev, lai, smoist, stemp_l(:,n_e),     &
        SWE, met_force%Tair, RTM_ID, Tb_h_vec, Tb_v_vec )

   Tb_h_l(:,j,n_e) = Tb_h_vec
   Tb_v_l(:,j,n_e) = Tb_v_vec
end do
```

Output arrays `Tb_h_l` and `Tb_v_l` are shaped `(N_tile, N_TbFreqAngRTMid, N_ens)`.
These become the predicted observations `Hx` in the EnKF update.

---

## Key references

- De Lannoy et al. 2013, J. Hydromet., doi:10.1175/JHM-D-12-092.1 — full description of RTM_ID=1
- Crow et al. 2005 — tau-omega equation (Eq. 1)
- Wang and Schmugge 1980, IEEE TGRS — soil dielectric model (RTM_ID 1,2)
- Mironov et al. 2009, IEEE TGRS, doi:10.1109/TGRS.2008.2011631 — soil dielectric model (RTM_ID 3,4)
- Pellarin (in CMEMv3.0) — atmospheric correction
