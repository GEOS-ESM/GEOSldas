
README.metforcing_and_bcs
=====================

Description:
---------------
Information about surface meteorological forcing and boundary conditions for GEOSldas.

Author:	
---------------
- reichle (5 Oct 2017) - first version; updated regularly
- jperket (10 Dec 2019) - converted to markdown


Surface Meteorological Forcing
================================================================================

The forcing time step is controlled with the configurable resource variable __FORCE_DTSTEP__ and __must match the frequency of the input forcing files__. 
Specify FORCE_DTSTEP=3600 [seconds] for most GEOS products, incl. MERRA, MERRA-2, FP, FP-IT/RP-IT, and GEOS-IT, which are 1-hourly datasets.

The spatial (horizontal) interpolation method for the met forcing is controlled
by `MET_HINTERP` (see optional parameters in `exeinp` input file to ldas_setup).

`MET_PATH` and `MET_TAG` must be consistent:

 - `MET_PATH` is the full path to the forcing data set.

 - `MET_TAG` is an identifier for the forcing data set.

 - Available non-GEOS forcing data sets are pre-defined in subroutine get_forcing().

 - For GEOS forcing datasets, see special `MET_TAG`
    parsing conventions in subroutines parse_MERRA_met_tag(), parse_MERRA2_met_tag()
    and parse_G5DAS_met_tag().

 - For details on corrected precipitation data see also (https://gmao.gsfc.nasa.gov/pubs/): 
   - Reichle and Liu (2014), 
     Observation-Corrected Precipitation Estimates in GEOS-5, 
     NASA Technical Report Series on Global Modeling and Data Assimilation, 
     NASA/TM-2014-104606, Vol. 35, National Aeronautics and Space Administration, 
     Goddard Space Flight Center, Greenbelt, Maryland, USA, 18pp.
   - Reichle, R. H., and Q. Liu (2021), 
     Observation-Corrected Precipitation for the SMAP Level 4 Soil Moisture (Version 6) Product and the GEOS R21C Reanalysis, 
     NASA Technical Report Series on Global Modeling and Data Assimilation, 
     NASA/TM-2021-104606, Vol. 59, National Aeronautics and Space Administration, 
     Goddard Space Flight Center, Greenbelt, Maryland, USA, 28pp.

COMMONLY USED values for `MET_PATH`:
------------------------------------

#### Legacy datasets
```
  MET_PATH : [XXX]/l_data/ECMWF/GRID/CORRECTED/netcdf/
  MET_PATH : [XXX]/l_data/GLDAS/netcdf/
  MET_PATH : [XXX]/l_data/GSWP-2/Baseline_Forcings/
  MET_PATH : [XXX]/l_data/RedArk/RedArk_subbasin_forcing/red_ark_forc/
```

#### ERA5 (LDAS-Monde via NASA LIS group)
```
  MET_PATH : /discover/nobackup/projects/lis/MET_FORCING/ERA5/
```

### GEOS datasets

#### MERRA  forcing (including precip-corrected MERRA forcing)
```
  MET_PATH : /discover/nobackup/projects/gmao/merra/iau/merra_land/MERRA_land_forcing/ 
```

#### MERRA2 forcing (including precip-corrected MERRA2 forcing)
```
  MET_PATH : /discover/nobackup/projects/gmao/merra/iau/merra_land/MERRA2_land_forcing/ 
```

Note: Used for SMAP Nature Run v05, v7.2, v8.1, v8.3, v9.1, v10.0 (before 1/1/2015).

#### GEOS FP forcing with "seamless" file names, for use with MET_TAG=GEOS.fp.asm[__prec*] (__PREFERRED__)
```
  MET_PATH : /discover/nobackup/projects/gmao/smap/SMAP_L4/GEOS/FP/
```

Note: Used for SMAP Nature Run v8.1, v8.3, v9.1, v10.0 (after 1/1/2015); SMAP L4_SM Version 5, 6, 7.

#### GEOS forcing with experiment-specific file names, incl. FP (__DEPRECATED__), FP-IT/RP-IT, and precip-corrected GEOS forcing
```
  MET_PATH : /discover/nobackup/projects/gmao/merra/iau/merra_land/GEOS5_land_forcing/
```

Note: Used for SMAP Nature Run v03, v04, v04.1, v05, v7.2 (after 1/1/2015); SMAP L4_SM Version 4.

#### GEOS-IT forcing (including precip-corrected GEOS-IT forcing)
```
  MET_PATH : /discover/nobackup/projects/gmao/merra/iau/merra_land/GEOSIT_land_forcing/
```

#### Forcing from post-processed output of the GEOS S2S system (FCST, AODAS)

```
  MET_PATH : [check with GMAO S2S group]
```



COMMONLY USED values for `MET_TAG`:
------------------------------------

#### Legacy datasets
```
  MET_TAG  : Berg_netcdf
  MET_TAG  : GLDAS_2x2_5_netcdf
  MET_TAG  : GSWP2_1x1_netcdf
  MET_TAG  : RedArk_ASCII
```

#### ERA5 (LDAS-Monde via NASA LIS group)
```
  MET_TAG  : ERA5_LIS
```

### GEOS datasets

#### MERRA
```
  MET_TAG  : d5_merra_cross__GEOSdas-2_1_4
```

#### MERRA-Land
```
  MET_TAG  : d5_merra_cross__GEOSdas-2_1_4__precCPCU
```

#### MERRA-2

 - All MERRA-2 options use native MERRA-2 "lfo" files for all surface met forcing fields *except* precipitation.

 - Option 1a: 
   Use corrected precip seen by the land w/in the MERRA-2 system (i.e., native MERRA-2 "lfo" files).
   Precip is as corrected within the MERRA-2 system.  Closest to land-only MERRA-2 replay. 
   
   ```
   MET_TAG  : M2COR_cross                                       ! all streams
   MET_TAG  : M2COR_100                                         ! stream 1 only
   MET_TAG  : M2COR_200                                         ! stream 2 only
   MET_TAG  : M2COR_300                                         ! stream 3 only
   MET_TAG  : M2COR_400                                         ! stream 4 only
   ```

 - Option 1b: 
   Use corrected precip forcing constructed in post-processing using MERRA-2 as background. 
   Background precip is typically from MERRA-2 "int" data, but corrected precip is stored
   in files that look like MERRA-2 "lfo" files.

   For example, select corrected precip version "CPCUGPCP22clim_MERRA2_BMTXS" as follows:
   ```
   MET_TAG  : M2COR_cross__precCPCUGPCP22clim_MERRA2_BMTXS      ! all streams
   MET_TAG  : M2COR_100__precCPCUGPCP22clim_MERRA2_BMTXS        ! stream 1 only
   MET_TAG  : M2COR_200__precCPCUGPCP22clim_MERRA2_BMTXS        ! stream 2 only
   MET_TAG  : M2COR_300__precCPCUGPCP22clim_MERRA2_BMTXS        ! stream 3 only
   MET_TAG  : M2COR_400__precCPCUGPCP22clim_MERRA2_BMTXS        ! stream 4 only
   ```
   
   This particular version uses as background the MERRA-2 model ("int") precip, rescaled to 
   the GPCPv2.2 climatology (indicated by "S" in "BMTXS").  Outside of Africa and the high 
   latitudes (tapered between 42.5 and 62.5 deg lat), this background precipitation is corrected 
   with CPCU data (also rescaled to the GPCPv2.2 climatology).  See Reichle and Liu (2014)
   GMAO Tech Memo #35 for more information on methods "BMTX".

 - Option 2:
   Use uncorrected precip generated by the AGCM w/in the MERRA-2 system (i.e., native 
   MERRA-2 "int" files).
   ```
   MET_TAG  : M2INT_cross                                       ! all streams
   MET_TAG  : M2INT_100                                         ! stream 1 only
   MET_TAG  : M2INT_200                                         ! stream 2 only
   MET_TAG  : M2INT_300                                         ! stream 3 only
   MET_TAG  : M2INT_400                                         ! stream 4 only
   ```
    
#### RP-IT/FP-IT (d591)
  ```
  MET_TAG  : cross_d591_RPFPIT                                  ! all streams
  MET_TAG  : d591_rpit1_jan00                                   ! stream 1 only
  MET_TAG  : d591_rpit2_jun06                                   ! stream 2 only
  MET_TAG  : d591_rpit3_jan11                                   ! stream 3 only
  MET_TAG  : d591_fpit                                          ! stream 4 only
  ```

#### RP-IT/FP-IT (d5124)
  ```
  MET_TAG  : cross_d5124_RPFPIT                                 ! all streams   - uses "late-look" through present
  MET_TAG  : d5124_rpit1_jan00                                  ! stream 1 only
  MET_TAG  : d5124_rpit2_jun04                                  ! stream 2 only
  MET_TAG  : d5124_rpit3_jan12                                  ! stream 3 only - updated through present
  ```

#### GEOS FP
  ```
  MET_TAG  : GEOS.fp.asm   ! PREFERRED:  "seamless" FP files (published/generic file names, ~same result as cross_FP)

  MET_TAG  : cross_FP      ! DEPRECATED: stitch FP experiment names across years
  MET_TAG  : e5110_fp      ! starting 11 Jun 2013
  MET_TAG  : e5130_fp      ! starting 20 Aug 2014
  MET_TAG  : e5131_fp      ! starting  1 May 2015
  MET_TAG  : f516_fp       ! starting 24 Jan 2017
  MET_TAG  : f517_fp       ! starting  1 Nov 2017
  MET_TAG  : f521_fp       ! starting 11 Jul 2018
  MET_TAG  : f522_fp       ! starting 13 Mar 2019
  MET_TAG  : f525_fp       ! starting 30 Jan 2020
  MET_TAG  : f525_p5_fp    ! starting  7 Apr 2020
  ...
  ```

#### GEOS-IT
  ```
  MET_TAG  : cross_d5294_GEOSIT
  ```

#### Forcing from post-processed output of the GEOS S2S system 

 - Forcing derived through post-processing of daily average output from the GEOS S2S system,
   including S2S hindcasts/forecasts ("FCST") and the "AODAS" used for S2S initialization. 

   S2S output is from the geosgcm_vis2d and geosgcm_surf Collections for FCST and from the 
   geosgcm_rad and geosgcm_surf Collections for AODAS (see GMAO Office Note No. 16).

   For FCST, post-processing includes a monthly bias correction to the MERRA-2 climatology.

   Daily data are disaggregated to 6-hourly (FCST) or 1-hourly (AODAS) using the MERRA-2 
   climatological diurnal cycle.

   For FCST, MET_TAG must specify S2S ensemble member ('ensX'; currently: 'ens1', 'ens2', 
   'ens3', or 'ens4') and month/day of forecast initialization ('MMMDD'; e.g., 'jan01'), 
   separated by double underscores.
   
   As of 14 Jun 2021:
   - Preparation of S2S forcing data ignores the 3-hour offset between S2S daily averages 
   (21z-21z) and the MERRA-2 daily averages (0z-0z) used for the temporal disaggregration.
   - The processing of the S2S output incorrectly partitioned total precipitation into snowfall 
   and convective precipitation.  Therefore, rainfall and snowfall are determined in the 
   S2S forcing reader from total precipitation and air temperature.  Convective rainfall is 
   set to 0.  (As of now, only total rainfall is used by Catchment.)

   ```
   MET_TAG  : GEOSs2sFCST__[ensX]__[MMMDD]
   MET_TAG  : GEOSs2sAODAS
   ```


#### SMAP L4_SM

#### Pre-beta SMAP L4_SM
```
  MET_TAG  : cross_FP__precCPCUG5FPv2
```

#### SMAP Nature Run v03
```
  MET_TAG  : cross_RPFPIT__precCPCUG5RPFPITv1                ! before 1/1/2014
  MET_TAG  : cross_FP__precCPCUG5FPv1                        ! after  1/1/2014
```

#### SMAP Nature Run v04
```
  MET_TAG  : cross_d591_RPFPIT__precCPCUG5RPFPITv2           ! before 1/1/2014
  MET_TAG  : cross_FP__precCPCUG5FPv2                        ! after  1/1/2014
```

#### SMAP Nature Run v04.1
```
  MET_TAG  : cross_d5124_RPFPIT__precCPCUG5RPFPITv2.1        ! before 1/1/2015
  MET_TAG  : cross_FP__precCPCUG5FPv2                        ! after  1/1/2015
```

#### SMAP Nature Run v05, v7.2, v8.1, v8.3;  SMAP L4_SM Version 4
```
  MET_TAG  : M2COR_cross__precCPCUGPCP22clim_MERRA2_BMTXS    ! before 1/1/2015
  MET_TAG  : cross_FP__precCPCUG5FPv3                        ! after  1/1/2015
```

#### SMAP L4_SM Version 5
```
  MET_TAG  : GEOS.fp.asm__precCPCULLKG5FPv3                  ! (precip corr with late-look CPCU)
  MET_TAG  : GEOS.fp.asm__precCPCUFLKG5FPv3                  ! (precip corr with first-look CPCU)
```

#### SMAP Nature Run v9.1; SMAP L4_SM Version 6
```
  MET_TAG  : M2COR_cross__precSMAPv6_CGIM2                   ! before 01/01/2015
  MET_TAG  : GEOS.fp.asm__precCPCU_IMGFinal_IMGFclim_G5FP    ! before 06/30/2021 (precip corr with IMERG-Final and late-look  CPCU) 
  MET_TAG  : GEOS.fp.asm__precCPCULLK_IMERGLateV06b_fp_v1    ! before 10/27/2021 (precip corr with IMERG-Late  and late-look  CPCU)
  MET_TAG  : GEOS.fp.asm__precCPCUFLK_IMERGLateV06b_fp_v1    ! before 05/09/2022 (precip corr with IMERG-Late  and first-look CPCU)
  MET_TAG  : GEOS.fp.asm__precCPCUFLK_IMERGLateV06c_fp_v1    !                   (precip corr with IMERG-Late  and first-look CPCU)
```

#### SMAP Nature Run v10.0*; SMAP L4_SM Version 7
```
  MET_TAG  : M2COR_cross__precSMAPv6_CGIM2                   ! before 01/01/2015 
  MET_TAG  : GEOS.fp.asm__precCPCU_IMGFinal_IMGFclim_G5FP    ! before 09/30/2021 (precip corr with IMERG-Final and late-look  CPCU) 
  MET_TAG  : GEOS.fp.asm__precCPCUFLK_IMERGLateV06b_fp_v1    ! before 05/09/2022 (precip corr with IMERG-Late  and first-look CPCU)
  MET_TAG  : GEOS.fp.asm__precCPCUFLK_IMERGLateV06c_fp_v1    ! before 07/02/2023 (precip corr with IMERG-Late  and first-look CPCU)
  MET_TAG  : GEOS.fp.asm__precCPCUFLK_IMERGLateV06d_fp_v1    ! before 11/07/2023 (precip corr with IMERG-Late  and first-look CPCU)
  MET_TAG  : GEOS.fp.asm__precCPCUFLK_IMERGLateV06e_fp_v1    !                   (precip corr with IMERG-Late  and first-look CPCU)

```
*Transitions dates are for L4_SM Version 7.  Transition dates for NRv10.0 may differ somewhat.


Boundary Conditions
================================================================================

Boundary conditions (bcs) are tile-space model parameters that are provided in a 
  set of files located in `BCS_PATH` for a given `BCS_RESOLUTION`. 

For "land" tiles, the discretization (tile-space) is constructed in one of two
  different ways:

  1. By intersecting watersheds with a regular grid (typically that used in
      the atmospheric model), e.g., `DC0576xPC0361`, where "DC" indicates that
      the grid cells are centered on the date line and "PC" indicates that
      the grid cells are centered on the poles.  Another example is the 
	  0.5-degree ("c180") cube-sphere grid used by the atmospheric model in the 
	  MERRA-2 reanalysis: `CF0180x6C_DE1440xPE0720`.

  2. Directly on a regular grid, e.g., `EASEv2_M09`.


**IMPORTANT**: Beginning with GEOSldas release v18.0.0, bcs must be provided in a revised
  directory layout and naming convention.  On NCCS/Discover, use:
```
  BCS_PATH : /discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/
```



COMMONLY USED boundary conditions (bcs):
----------------------------------------

#### MERRA-Land (and MERRA)
```
  BCS_PATH : /discover/nobackup/projects/gmao/ssd/land/l_data/geos5/bcs/SiB2_V2_bad_lon_onDL/DC/
```

#### Ganymed-4_0 (GM4): MERRA2
```
  BCS_PATH : /discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Ganymed-4_0/Ganymed-4_0_MERRA-2/
```
  Note: The same land bcs (for cube-shere only) are also in /discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus/Icarus_MERRA-2/


#### SMAP Nature Run v03
```
  BCS_PATH : /discover/nobackup/projects/gmao/ssd/land/l_data/geos5/bcs/CLSM_params/mkCatchParam_v15/
```

#### SMAP Nature Run v04, v04.1
```
  BCS_PATH : /discover/nobackup/projects/gmao/ssd/land/l_data/geos5/bcs/CLSM_params/mkCatchParam_SMAP_L4SM_v001/
```

#### SMAP Nature Run v05
```
  BCS_PATH : /discover/nobackup/projects/gmao/ssd/land/l_data/geos5/bcs/CLSM_params/mkCatchParam_SMAP_L4SM_v002/
```

#### Icarus-NL ("New Land"): SMAP Nature Run v7.2
```
  BCS_PATH : /discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NL/Icarus-NL_[XXXX]/
```

Notes:
- _DO NOT USE_ unless replicating previous experiments. There is "missing" data in green*.data, nirdf*.dat, and visdf*.dat files.
- This path remains in place to permit recreating experiments that have used this path.
- The sub-directory "Icarus-NL_MERRA-2/" contains the "new land" bcs.  The string "MERRA-2" in this sub-directory name refers to ocean bcs that are not relevant for GEOSldas.

#### Icarus-NLv2: SMAP L4_SM Version 4
```
  BCS_PATH : /discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_[XXXX]/
```

Notes: 
- Icarus-NLv2 is a update to Icarus-NL bcs. A patch has been applied to files green*.data, nirdf*.dat, and visdf*.dat. 
	
#### Icarus-NLv3 (NL3): SMAP Nature Run v8.1, GEOS-FP 5.25, 5.27, 5.29, GEOS-IT
```
  BCS_PATH : /discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv3/Icarus-NLv3_[XXXX]/
```

Notes: 
- Soil parameters for a small fraction (< 0.05%) of tiles changed to correct "Mali" bug.
- Vegdyn.data now netcdf4; reverts to using Dorman/Sellers veg heights (abandons JPL/Simard et al. 2011 Lidar data). 
- Some underlying ASCII data files are now grouped into a netcdf4 file. I.e., data in ar.new, bf.dat, ts.dat, etc are now in:
    - clsm/catch_params.nc4
    - clsm/catchcn_params.nc4  (additional parameters for CatchCN)

- Generated with cvs tag Jason-3_0_LANDBCS

#### Icarus-NLv4 (NL4): SMAP Nature Run v8.3, v9.1, SMAP L4_SM Version 5, 6
```
  BCS_PATH : /discover/nobackup/projects/gmao/smap/SMAP_L4/L4_SM/bcs/CLSM_params/Icarus-NLv4_EASE/
```

Notes: 
- Icarus-NLv4 is identical to Icarus-NLv3 except that NLv4 reinstates veg heights from JPL/Simard et al. 2011 Lidar data.
- Generated with GEOSldas tag v17.9.0-beta.7 under SLES11 O/S.


#### Icarus-NLv5 (NL5): SMAP Nature Run v10.0, SMAP L4_SM Version 7
```
  BCS_PATH : /discover/nobackup/projects/gmao/smap/SMAP_L4/L4_SM/bcs/CLSM_params/Icarus-NLv5_EASE/
```

Notes: 
- Icarus-NLv5 is identical to Icarus-NLv4 except for parameters associated with PEATCLSM.


#### v11
```
  BCS_PATH : /discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/v11/
```

Notes: 
- v11 is identical to Icarus-NLv5 (within roundoff) except for new MODIS-based snow albedo v2 parameters.










