# CYGNSS Preprocessed Operator Slide Draft

## Slide 1: What We Implemented

We added a new CYGNSS L1 scalar observation capability to GEOSldas.

- Added a new observation species:
  - `CYGNSS_L1_DDM3X5_CROP_SCALAR`
- Added native and forecast-equivalent variable support for:
  - `cygl1scal`, in `dB`
- Built a dedicated reader for preprocessed CYGNSS L1 daily NetCDF products.
- Reader supports daily archive layout:
  - `CYGNSS_L1/Yyyyy/Mmm/cygnss_l1_ddm3x5_crop_scalar_m36_yyyymmdd_cyg02.nc4`
- Reader handles assimilation windows that cross midnight by reading multiple daily files.
- Reader filters observations to the current EnKF assimilation window.
- Reader rejects bad owner-tile entries.
- Reader keeps one observation per owner tile, choosing the closest specular point.

We implemented the CYGNSS forward operator inside GEOSldas.

- For each CYGNSS observation, GEOSldas identifies the owner tile and local support-tile halo.
- The operator gathers model soil moisture and mwRTM-related static inputs for the support tiles.
- It uses preprocessed coefficients from the CYGNSS operator pipeline.
- It reuses the mwRTM low-resolution reflectivity calculation for each support tile, using L-band frequency, specular-point incidence angle, clay, porosity, and model surface soil moisture.
- It combines nearby support-tile contributions into one predicted CYGNSS scalar value.
- It calculates the model-equivalent CYGNSS reflectivity-like scalar as `H(x) = sum_t C_t * R_t(x)`, then converts the linear reflectivity prediction to `dB` with `10*log10`.
- It applies model-state QC for missing soil moisture, missing clay/porosity, invalid reflectivity calls, and nonphysical predicted values.

We connected CYGNSS to the existing EnKF machinery.

- Added `cygl1scal` handling in the GEOSldas observation-prediction pathway.
- Threaded assimilation date/time and timestep information into the prediction call.
- Added efficient product caching so daily coefficient files are not repeatedly reopened for each observation.
- Wrote CYGNSS forecast and analysis equivalents to `ObsFcstAna`.
- Added CYGNSS L1 scalar to `update_type = 12`.
- CYGNSS-only updates affect soil-moisture state variables, not temperature states.
- Preserved existing `update_type = 12` behavior for snow and other sensors.

## Slide 2: What We Validated So Far

We ran a 2-day Discover smoke test using example daily CYGNSS files for 2019-11-01 and 2019-11-02.

- The GEOSldas run completed successfully:
  - `GEOSldas Run Status: 0`
- The run used `update_type = 12`.
- CYGNSS L1 scalar appeared in `ObsFcstAna` as an assimilated species.
- The daily file logic worked:
  - one daily file was read for normal windows
  - two daily files were read for the midnight-crossing window
- CYGNSS observations were read, filtered by time window, and mapped into tile space.
- CYGNSS observations were assimilated at expected 3-hour windows.

Assimilated CYGNSS observation counts in `ObsFcstAna`:

- 2019-11-01 09z: 2
- 2019-11-01 12z: 10
- 2019-11-01 15z: 11
- 2019-11-02 09z: 6
- 2019-11-02 12z: 10
- 2019-11-02 15z: 9

We confirmed that assimilation produced land increments.

- Nonzero increments appeared in:
  - `SRFEXC_INCR`
  - `RZEXC_INCR`
  - `CATDEF_INCR`
- Soil temperature and heat-content increments remained zero for CYGNSS-only updates, as intended.
- Obs-space residuals generally moved in the right direction from forecast to analysis.
- This is a successful technical smoke test.
- It is not yet a full science validation.

## Slide 3: What Still Needs To Be Done

Science tuning and sensitivity testing:

- Tune CYGNSS observation error standard deviation.
- Tune localization length scale.
- Revisit `xcorr`, `ycorr`, `xcompact`, and `ycompact` settings.
- Decide whether bias correction or scaling is needed.
- Evaluate whether current dB residual magnitudes are reasonable.
- Test whether increments are physically plausible in space and time.

Known caveat:

- The current test used a temporary localization choice:
  - `xcorr = ycorr = 1.25`
  - `xcompact = ycompact = 2.5`
- This allowed the run to proceed, but it has science implications.
- We need to revisit this before making science claims.

Product and workflow expansion:

- Test more than one CYGNSS spacecraft/day pair.
- Test larger daily archives with more observations.
- Decide final daily product naming and archive conventions.
- Decide how to handle multiple spacecraft and/or multiple files per day.
- Strengthen product metadata requirements.
- Continue documenting Discover setup and run workflow.

Validation next steps:

- Compare CYGNSS assimilation runs against no-CYGNSS baseline runs.
- Map increments relative to CYGNSS observation locations and support halos.
- Examine impacts on:
  - surface soil moisture
  - root-zone soil moisture
  - profile soil moisture
  - ensemble spread
  - observation-minus-forecast statistics
- Run longer experiments after the technical behavior is stable.
- Decide whether the operator is ready for broader CYGNSS experiment design.
