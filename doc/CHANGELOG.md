
GEOSldas History / Changelog
===============================

Description:
------------
History of GEOS LDAS ("LDASsa"; "GEOSldas") development in Git and CVS

Author:
--------
- reichle (20 Jul 2010), first version; updated regularly
- jperket (10 Dec 2019), converted to Markdown

Summary and Objective:
----------------------
The development of the off-line (land-only) GEOS LDAS land modeling and assimilation 
system started with the "land EnKF driver" for the GEOS Catchment model in the early 
2000s.  

From the mid 2000s, development was formalized as the "LDASsa" project under CVS version 
control.

Between 2016 and 2019, LDASsa was rewritten for full ESMF compliance under the "GEOSldas" 
project, with the ultimate goal of bringing the land model updates from the off-line 
system into the GEOS AGCM and integrating the GEOS LDAS into the atmospheric DAS.  During 
this period, LDASsa and GEOSldas development continued in parallel.

In 2019, GEOS LDAS version control transferred from CVS to Git.

This README file contains the history of stable GEOSldas versions ("tags") in Git, followed by older, CVS LDASsa and GEOSldas versions and change logs.


Overview of Git Releases:
============================

[v17.9.0-beta.5](https://github.com/GEOS-ESM/GEOSldas/releases/tag/v17.9.0-beta.5) - 2020-05-11
------------------------------
- Pre-release meant for use under SLES12 at NCCS.  Still works for SLES11.

- New/Updated Science Functionality:

  - Forecast error covariance inflation with scalar (globally constant) factor.

- New/Updated Infrastructure:

  - Support for GEOS FP forcing with generic ("seamless") file names.
  - Resource parameter changes:
    - Renamed NUM_ENSEMBLE to NUM_LDAS_ENSEMBLE in "exeinp" file to be consistent with LDAS.rc.
    - Renamed MONTHLY_OUTPUT to POSTPROC_HIST.
  - Updated utilities to MAPL v2.1.3, ESMA_env v2.1.3+intel19.1.0.
  
- Bug Fixes and Other Minor Changes:

  - Added basic protections for concatenation of sub-daily into daily nc4 files and for generation of monthly-mean nc4 files.
  - Write ObsFcstAna and smapL4SMaup files into ./scratch, then move to ana/ens_avg/year/month dir in postprocessing.
  - Some cleanup of obsolete LDASsa code.

------------------------------
[v17.9.0-beta.4-SLES12](https://github.com/GEOS-ESM/GEOSldas/releases/tag/v17.9.0-beta.4-SLES12) - 2020-04-24
------------------------------
- Pre-release meant for use under SLES12 at NCCS, otherwise identical to v17.9.0-beta.4-SLES11.
- Works under SLES12 using the Intel-19 compiler.
- Also works under SLES11 using the Intel-18 compiler but is not zero-diff across compilers/operating systems.

------------------------------
[v17.9.0-beta.4-SLES11](https://github.com/GEOS-ESM/GEOSldas/releases/tag/v17.9.0-beta.4-SLES11) - 2020-04-23
------------------------------
- Pre-release meant for use under SLES11 at NCCS. Under SLES12, use v17.9.0-beta.4-SLES12 (or newer).
- Uses the Intel-18 compiler and also appears to work under SLES12. However, LDASsa with Intel-18 under SLES12 was found to create bad Fortran sequential binary files out of a subroutine that is very similar in LDASsa and GEOSldas.
- Zero-diff vs. v17.9.0-beta.3 for Catchment only (except SMAP L1C Tb fore-minus-aft check).
- Not zero-diff for CatchCN (via v1.8.3 of GEOS_GCMGridComp).

- New/Updated Science Functionality:

  - Resurrected SMAP L1C Tb fore-minus-aft check.

- New/Updated Infrastructure:

  - Updated utilities to MAPL v2.1.1, ESMA_env v2.1.1., ESMA_cmake v3.0.1.
  - New GEOS_SurfaceGridComp.rc file (via v1.8.3 of GEOS_GCMGridComp).
  - Parallel post-processing.
  - Cross-stream support for FP f525_p5 forcing.
  - ~sbatch~ submission for pre-processing of restarts to comply with SLES12 requirements.
  - Subdaily-to-daily concatenation processes before month is complete.
  - Temporary solution to create directories for ObsFcstAna files to enable extending an existing GEOSldas run without going through setup.

- Bug Fixes and Other Minor Changes:

  - Updated README.md.
  - ~obspertrseed~ restart file name when restarting from existing run.
  - Subdaily-to-daily nc4 concatenation (indent error).
  - Fixes for GNU compiler in debug mode.
  - Fixed ~landpert~ checkpoint output when on cube-sphere tiles.

------------------------------
[v17.9.0-beta.3](https://github.com/GEOS-ESM/GEOSldas/releases/tag/v17.9.0-beta.3) - 2020-03-18
------------------------------
- Additional RESTART options, incl. from re-tiling MERRA-2, FP, or other restarts on different tile space or with different boundary conditions
- Bug fixes

------------------------------
[v17.9.0-beta.2](https://github.com/GEOS-ESM/GEOSldas/releases/tag/v17.9.0-beta.2) - 2020-02-26
------------------------------
- New/Updated Science Functionality:

  - Assimilation when running on cube-sphere tiles.
  - Read forcing from cube-sphere grid when running on matching cube-sphere tiles.
  - Output of Catchment analysis increments via HISTORY.  
  - Added FP-5.25 upgrade (30 Jan 2020) to "cross-stream" forcing option.
  - Functionality to create regional (non-global) nc4 vegdyn restart file.
  - Configuration option to add extra variables into catch restart files (as needed by GCM).
  - Allows processing of (assimilation) observations for innovations output *without* perturbations turned on.

- New/Updated Infrastructure:

  - Support for SLES 12 in addition to SLES11 (ESMA_env v2.0.2).
  - Updated to MAPL v2.0.
  - Removed dycore and FMS
  - Conforms to GNU compiler (gcc-9.1).
  - Post-processing compression (gzip) of landpert restart files (except final time).
  - Added LDAS_app/mk_GEOSldasRestarts.F90 (adapted from GCM GridComp's mk_LDASsaRestarts.F90 in preparation for re-tiling changes).
  - Fixed output log file name and location.

- Bug Fixes and Other Minor Changes: 

  - Bug fix in select-update_type 9 (abs(deltaT)>0.)
  - Bug fix for local mwRTM and time dimension restart.
  - Replaced copy ("cp") with link ("ln") for catparam and mwrtm diagnostic output files.

------------------------------
[v17.9.0-beta.1](https://github.com/GEOS-ESM/GEOSldas/releases/tag/v17.9.0-beta.1) - 2020-01-17
------------------------------
- Commented out call to check_catch_progn in apply_progn_pert.   

------------------------------
[v17.9.0-beta.0] - 2019-12-20
------------------------------
- First tag for SMAP L4_SM Version 5, Catchment model consistent with f525land_fpp
  - Changed default Z0_FORMULATION to 4 (incl. addition of simple tree SAI)
    (to be used with Icarus-NLv3 as default BCs; reverting to look-up veg heights)
  - Fix PAR perturbations bug affecting assim
  - Fix GEOS forcing stream boundaries bug 
  - Changes to time stepping for sun angle

------------------------------
[v17.8.0] - 2019-12-10
------------------------------
- Closest match to LDASsa CVS tag LDASsa_m3-16_6_p2, the LDASsa tag used for generating the Version 4 L4_SM product. 
- v17.8.0 is a debugged version of GEOSldas_m4-17_8.
			    
------------------------------

			    

Overview of Previous CVS tags:
=======================================

Tags are ordered by *version*, beginning with LDASsa tags and followed by GEOSldas tags.  

Between 4 Oct 2017 and 7 Mar 2019, LDASsa and GEOSldas development overlapped and 
new tags were created for both software systems.  During this period, there is some 
overlap between LDASsa and GEOSldas tags in terms of science versions, with LDASsa
tags primarily intended for the SMAP L4_SM ops system. 

```text

=====================================================================================
CVS tag					Created      Application/Comments
=====================================================================================
reichle-LDASsa_m2-10                       Oct 2012  Soil Parameter Testbed
-------------------------------------------------------------------------------------
reichle-LDASsa_m2-11                       May 2013  SMAP Delivery 4
reichle-LDASsa_m2-SMAP_L4_SM_D00400
reichle-LDASsa_m2-SMAP_L4_SM_D00400_p1     Jun 2013
-------------------------------------------------------------------------------------
reichle-LDASsa_m2-12                       Aug 2013  MERRA-2 restarts 
                                                      do *NOT* use!
                                                      (bug in veg/alb time interpolation)
-------------------------------------------------------------------------------------
reichle-LDASsa_m2-13                    12 Dec 2013  SMOS assimilation	
reichle-LDASsa_m2-13_p1                 13 Jan 2013
reichle-LDASsa_m2-13_p2                 24 Jan 2013
-------------------------------------------------------------------------------------
reichle-LDASsa_m2-SMAP_Nature_v03_spin  10 Feb 2014  Spin-up of SMAP Nature Run v03
                                                      zero-diff vs. tag "SMAP_Nature_v03" 
                                                      except for SMOS obs handling etc
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m2-SMAP_L4_SM_D00500     19 Feb 2014  SMAP Delivery 5 
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m2-SMAP_Nature_v03        4 Mar 2014  GMAO SMAP Nature Run v03 
                                                      do *NOT* use for data assimilation!
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m2-14                     6 Mar 2014    
reichle-LDASsa_m2-SMAP_L4_SM_D00500_p1               SMAP Delivery 5, p1 (same as m2-14)  
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m2-14_1                   2 Apr 2014  
reichle-LDASsa_m2-SMAP_L4_SM_D00500_p2               SMAP Delivery 5, p2 (same as m2-14_1)  
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m2-14_1_p1               23 May 2014  (fixed MPI load-balanced analysis 
                                                      for 1d update types)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m2-SMAP_L4_SM_ORT        15 Jul 2014  SMAP Operational Readiness Test
                                                     (same as "m2-14_2" but with bad
                                                      oscillations fix)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m2-14_2                   7 Aug 2014  surface energy balance oscillations fix
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m2-14_2_p1               17 Oct 2014  fixes for IO-related hangs and 
                                                      obs error covariance
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m2-SMAP_L4_SM_ORT2       13 Nov 2014  SMAP Operational Readiness Test Phase 2
                                                      (same as "ORT" but with temp incr fix)
-------------------------------------------------------------------------------------
reichle-LDASsa_m3-15_0                  25 Nov 2014  GEOSsurface_GridComp as in Ganymed-4_0
                                                      but with oscillations fix;
                                                     major clean-up;
                                                     scientifically equiv. to "m2-14_2_p1" 
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-SMAP_Nature_v04_spin   5 Dec 2014  Spin-up of SMAP Nature Run v04
                                                     (Ganymed-4_1 roughness length)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-SMAP_Nature_v04_beta   9 Dec 2014  SMAP Nature Run v04 (except Tb FOV)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-15_1                   2 Apr 2015  revised Tb FOV
reichle-LDASsa_m3-SMAP_Nature_v04                  
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-15_1_p1               23 Apr 2015  minor bug fixes
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-15_1_p2                6 May 2015  plumbing for ascending L2_SM_AP
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-15_2                   1 Jul 2015  new compiler, minor fixes 
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-15_2_p1                8 Oct 2015  forcing updates (added MERRA-2)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-15_2_p2                1 Jun 2016  L1C_TB check for bad half-orbits
                                                     (patch for SMAP L4_SM ops)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-15_2_p2a              19 Jan 2017  FP-5.16 cross-stream functionality
                                                     (patch for SMAP L4_SM ops)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-15_2_p3                8 Jun 2016  Revised ObsFcstAna and L4_SM aup output
                                                     (patch for SMAP L4_SM testing)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-15_2_p4               28 Dec 2016  Added d5124 RPIT/FPIT functionality
                                                     (patch targeted for SMAP L4_SM ops)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-15_2_p5                9 Feb 2017  Added GEOS-5.16 FP functionality
                                                     (patch targeted for NRv4.1 and SMAP L4_SM ops)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-15_2_p6               17 Apr 2017  Fixed bug with Tb_fcst=0 output
                                                     (patch targeted for NRv4.1 and SMAP L4_SM ops)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-15_2_p7                8 Jun 2017  Added run-time input of obs file name lists
                                                     (for non L4-ops SMAP assim experiment)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-15_2_p8               30 Oct 2017  Added GEOS-5.17 FP functionality
                                                     (patch targeted for NRv4.1 and SMAP L4_SM ops)
-------------------------------------------------------------------------------------
reichle-LDASsa_m3-SMAP_Nature_v05_spinA 20 Nov 2015  spin-up of SMAP Nature Run v05 
                                                     (same as "reichle-LDASsa_m3-SMAP_Nature_v05"
                                                      EXCEPT not using dampen_tc_oscillations() 
                                                      for veg type 1, plus zero-diff changes)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-SMAP_Nature_v05        7 Dec 2015  SMAP Nature Run v05 (same model as 
                                                      "m3-16_0" except MERRA-2 SWGDN from "lfo"
                                                      and mwRTM soiltemp change)
-------------------------------------------------------------------------------------
reichle-LDASsa_m3-16_0                  23 Dec 2015  Revised Catchment/snow model constants
                                                      (WEMIN, CSOIL)
                                                     New roughness length formulation,
                                                      (JPL) veg height from file
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-16_0_p1                3 Feb 2016  Minor fixes.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-16_0_p2               25 Mar 2016  Minor fixes & speed-up.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-16_5                  31 Oct 2017  Revised Catchment model 
                                                      Reduced upward flow from srfexc to rzexc
                                                      SRTM-based tile space
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-SMAP_Nature_v07       31 Oct 2017  Identical to "reichle-LDASsa_m3-16_5"
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-16_6                  23 Jan 2018  Revised Catchment model
                                                      Increased upward flow from srfexc to rzexc
                                                     L1C_TB reader fix (use avg fore/aft timestamp)
                                                     Added update_type=10, L1C_TB_E reader
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-SMAP_Nature_v07_2      2 Mar 2018  Identical to "reichle-LDASsa_m3-16_6"
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-16_6_p1                9 Jul 2018  Added GEOS-5.21 FP functionality
                                                     (patch targeted for NRv7.2 and SMAP L4_SM ops)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-16_6_p2                7 Mar 2019  Added GEOS-5.22 FP functionality
                                                     (patch targeted for NRv7.2 and SMAP L4_SM ops)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-16_6_p3               30 Jan 2020  Added GEOS-5.25 FP functionality
                                                     (patch targeted for NRv7.2 and SMAP L4_SM ops)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-16_6_p4                3 Apr 2020  Added GEOS-5.25_p5 FP functionality
                                                     (patch targeted for NRv7.2 and SMAP L4_SM ops)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reichle-LDASsa_m3-16_6_p4_SLES12        14 Apr 2020  SLES12 version of *_p4 tag -- NOT zero-diff!!
                                                     (patch targeted for NRv7.2 and SMAP L4_SM ops)

-------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------
GEOSldas_m4-17_0                         4 Oct 2017  ESMF-compliant driver; new model version
                                                     - no assimilation
                                                     - does not incl LDASsa changes post 1Jun2016
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
GEOSldas_m4-17_6                        15 Mar 2018  Default Catchment model configuration matches 
                                                      science of reichle-LDASsa_m3-16_6
					             - incl some land assimilation components, but 
                                                       not functional
						     - improved load balancing on cube-sphere grid 
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
GEOSldas_m4-17_7                        16 Aug 2018  Added land assimilation and CatchmentCN 
                                                     - land assimilation not fully functional
						     - configurable LAND_PARAMS input replaces 
                                                       LAND_UPD compiler flag 
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
GEOSldas_m4-17_8                        13 Mar 2019  Land assimilation functional for SMAP L4_SM
                                                      but not yet for cube-sphere tile space.
                                                     CatchmentCN updates towards S2S v3:
						     - prescribed LAI/SAI parameters
						     - atmospheric CO2 runtime options
						     - scaled albedo and FPAR options
						     New batchjobs.j for submission of jobs in a 
                                                      way that reduces queue wait time (optional)
                                                     New irrigation elements (under development)
                                                     
===================================================================================== 



Change log and historic README files (in reverse chronological order):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Change log is in reverse chronological order beginning with GEOSldas and 
followed by LDASsa and precursor change log files.
 
Note that between 4 Oct 2017 and 7 Mar 2019, LDASsa and GEOSldas development 
overlapped.

***********************************************************************************
*                                                                                 *
*                           GEOSldas CVS change log                               *
*                                                                                 *
***********************************************************************************


===================================================================
===  Approximate stable point of GEOSldas_m4-17_8 tag           ===
===================================================================

===================================================================
===                                                             ===
===               02Feb2019                                     ===
===                                                             ===
===================================================================
- new function optionally applies increments to catchment
   during coupled runs, controlled by "LDAS_IAU" (not in rc files)
   
===================================================================
===                                                             ===
===               19Nov2018                                     ===
===                                                             ===
===================================================================
- Multiple modifications to GEOS_CatchCNGridComp from branch b_pLAI2:
  - added irrigation (for future,experimental)
  - added PRESCRIBE_DVG
  - replaced DO_CO2SC with ATM_CO2, giving mulptiple CO2 options:
    0: use fixed value
    1: use CT tracker monthly mean diurnal cycle
    2: use scaled CT tracker monthly mean diurnal cycle to match
       global EEA average
    3: use spatially fixed, interannually varying CMIP values
	(AGCM only)
    4: use AGCM values (AGCM only)
           
===================================================================
===                                                             ===
===               31Oct2018                                     ===
===                                                             ===
===================================================================
- introduced prescribed LAI and SAI for CatchCN

===================================================================
===                                                             ===
===               31Oct2018                                     ===
===                                                             ===
===================================================================
- introduced batchrun.j for submitting lenkf.j jobs with slurm 
  dependency (to reduce queue wait time)

===================================================================
===                                                             ===
===               01Oct2018                                     ===
===                                                             ===
===================================================================
- introduced irrigation model (experimental) activated by 'RUN_IRRIG' 
  - added module to process / write irrigation restarts
  - new file: GEOSsurface_GridComp/Shared/Raster/src/irrg_model.F90

===================================================================
===                                                             ===
===               25Sep2018                                     ===
===                                                             ===
===================================================================
- patched bug by removing moisture/energy corrections under snow
- added debugging tool for catchCN:
   GEOScatchCN_GridComp/dbg_cnlsm_offline.F90
- added modified albedo scaling scheme for catchCN:
   GEOScatchCN_GridComp/compute_FPAR_CDF_M09.F90
- added math routines in catchCN to derive FPAR and Albedo scale
  parameters in order to match MODIS FPAR and VISDF and NIRDF:
   GEOScatchCN_GridComp/math_routines.F90

===================================================================
===                                                             ===
===               14Sep2018                                     ===
===                                                             ===
===================================================================
- number of sub-daily files now checked before concatenated and
   monthly means created


===================================================================
===  Approximate stable point of GEOSldas_m4-17_7 tag           ===
===================================================================

===================================================================
===                                                             ===
===               07Sep2018                                     ===
===                                                             ===
===================================================================
- change MONTHLY_ONLY to MONTHLY_OUTPUT in GEOSldas_LDAS.rc,
  ldas_setup, lenkf.j
- change CATCHCN_DT to DTCN in GEOSldas_LDAS.rc
- updates to tutorial

===================================================================
===                                                             ===
===               16Aug2018                                     ===
===                                                             ===
===================================================================
- add debugger choice command argument to lenkf.j
- add optional CLSM debugging, activated by WRITE_LAND_BUDGET ifdef

===================================================================
===                                                             ===
===               06Aug2018                                     ===
===                                                             ===
===================================================================
- fix ldas_setup to be run from Linux/bin

===================================================================
===                                                             ===
===               02Aug2018                                     ===
===                                                             ===
===================================================================
- rename RC files from LDASsa to GEOSldas. Now files are:
   Applications/LDAS_App/GEOSldas_CAP.rc
   Applications/LDAS_App/GEOSldas_ExtData.rc
   Applications/LDAS_App/GEOSldas_HIST.rc
   Applications/LDAS_App/GEOSldas_LDAS.rc
- add configurable LAND_PARAMS rc options to quickly switch physics 
  parameter combinations, replacing compile-time LAND_UPD ifdefs 
- preprocess_ldas.F90 makes sure each PE has >=2 lats or lons
- regrid tool now uses Icarus-NLv2, which fixes phenology bug
- added description of Icarus-NLv2 bcs to README.GEOSldas_metforcing_and_bcs,
  discourage use of Icarus-NL bcs

===================================================================
===                                                             ===
===               25Jul2018                                     ===
===                                                             ===
===================================================================
- add matlab scripts from reichle-LDASsa_m3_16_5:
   read_ObsFcstAna.m
   read_catparam.m
   read_obslog.m
   read_smapL4SMaup.m
   read_smapL4SMlmc.m
- update to create cube-sphere ocean tiles in boundary coundtions

===================================================================
===                                                             ===
===               18Jul2018                                     ===
===                                                             ===
===================================================================
- add SMAP support to tile_bin2nc4.F90, and some name corrections:
   NLAND  -> tile
   TIME   -> time
   SNOMAS -> SNOWMASS

===================================================================
===                                                             ===
===               13Jul2018                                     ===
===                                                             ===
===================================================================
- add ability to read GEOS FP v5.21 forcing data

===================================================================
===                                                             ===
===               13Jul2018                                     ===
===                                                             ===
===================================================================
- bypass ensemble GridComp for catchCN in GEOS_LdasGridComp.F90
  (catchCN does not currently work with assim)

===================================================================
===                                                             ===
===               11Jul2018                                     ===
===                                                             ===
===================================================================
- add two default namelist files for assim

===================================================================
===                                                             ===
===               07Jul2018                                     ===
===                                                             ===
===================================================================
- add two default namelist files for assim
- add profiling timer to GEOS_LandAssimGridComp.F90

===================================================================
===                                                             ===
===               26Jun2018                                     ===
===                                                             ===
===================================================================
- changed Z0_FORMULATION to 3
- wrtm_bin2nc4.F90 checks for no data

===================================================================
===                                                             ===
===               20Jun2018                                     ===
===                                                             ===
===================================================================
- further updates to conform LDASsa in GEOS_LdasGridComp.F90,
  GEOS_LandAssimGridComp.F90
- add LDAS_Forcing print statements to log
  
===================================================================
===                                                             ===
===               14Jun2018                                     ===
===                                                             ===
===================================================================
- add LDASsa output subroutines for assim, clsm_ensdrv_out_routines.F90
- add land assim support in ldas_setup
- add cat param and mwRTM param in lenkf.j, preprocess_ldas.F90
- add subroutine to read/write pertubation restarts in LDAS_PertRoutines.F90

===================================================================
===                                                             ===
===               19Apr2018                                     ===
===                                                             ===
===================================================================
- add support of AppGridCreate to FVdycoreCubed_GridComp/c2l_CFIO_offline.F90 
  (with JMS.rc)

===================================================================
===                                                             ===
===               10Apr2018                                     ===
===                                                             ===
===================================================================
- add support for different ensemble member restart file
  and land assim in ldas_setup
- add support for JMS.rc (cube-sphere layout)

===================================================================
===                                                             ===
===               05Apr2018                                     ===
===                                                             ===
===================================================================
- Assimilation updates:
  - updating assim results and pert_seed in GEOS_LdasGridComp.F90
  - add update_assim to GEOS_LandAssimGridComp.F90
  - Return LDASsa get_pert subroutine
  
===================================================================
===                                                             ===
===               11Mar2018                                     ===
===                                                             ===
===================================================================
- change mwrtm param restart file name in ldas_setup

===================================================================
===  Approximate stable point of GEOSldas_m4-17_6 tag           ===
===================================================================

===================================================================
===                                                             ===
===               03Mar2018                                     ===
===                                                             ===
===================================================================
- fix bug in LDAS_PertRoutines.F90: N_force => N_progn

===================================================================
===                                                             ===
===               20Dec2017                                     ===
===                                                             ===
===================================================================
- update matlab utilties to read binaries from GEOSldas

===================================================================
===                                                             ===
===               15Dec2017                                     ===
===                                                             ===
===================================================================
- remove debugging print statements in met forcing

===================================================================
===                                                             ===
===               13Dec2017                                     ===
===                                                             ===
===================================================================
- add fix for CF grid ntasks to lenkf.j
- move LDAS_PertTypes.F90 from landpert_gridcomp to GEOSldas_GridComp/Shared/

===================================================================
===                                                             ===
===               11Dec2017                                     ===
===                                                             ===
===================================================================
- add land assimilation files:
   GEOSldas_GridComp/GEOSlandassim_GridComp/adapt_types.F90
   GEOSldas_GridComp/GEOSlandassim_GridComp/catch_bias_types.F90
   GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_adapt_routines.F90
   GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_bias_routines.F90
   GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensdrv_drv_routines.F90
   GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_enkf_update.F90
   GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_glob_param.F90
   GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_read_obs.F90
   GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_upd_routines.F90
   GEOSldas_GridComp/GEOSlandassim_GridComp/enkf_general.F90
   GEOSldas_GridComp/GEOSlandassim_GridComp/io_hdf5.F90
   GEOSldas_GridComp/GEOSlandassim_GridComp/mwRTM_routines.F90
   GEOSldas_GridComp/GEOSlandassim_GridComp/mwRTM_types.F90
   GEOSldas_GridComp/Shared/enkf_types.F90
   GEOSldas_GridComp/Shared/my_matrix_functions.F90

===================================================================
===                                                             ===
===               03Nov2017                                     ===
===                                                             ===
===================================================================
- incorporate change from tag reichle-LDASsa_m3-16_5 
  - added functionality for d5124 RPIT/FPIT forcing
  - added FP-5.16 transition (24 Jan 2017) to "cross-stream" dates
  - added FP-5.17 transition (1 Nov 2017) to "cross-stream" dates
- add description of MET_HINTERP
- add GEOSmwRTM_GridComp
- ldas_setup now looks for names of the restart files from LDASsa
- some cleanup and move around files

===================================================================
===                                                             ===
===                4Oct2017                                     ===
===                                                             ===
===================================================================
- Inaugural GEOSldas stable tag (GEOSldas_m4-17_0):
 - new CVS module "m4"
 - updated land model version ("17"; a.k.a. "v24_c05")
 - updated to Intel 17.0.4.196 and OpenMPI 2.1.1
 - major revisions in ldas_setup
   .preprocess_ldas : optimize the grid distribution and recreate tile, BCs
   .process_rst     : make restart 
   .process_hist    : preprocess history rc file
   .tile_bin2nc4    : post process output, convert binary to nc4
 - added GEOScatchCN_GridComp
 - tutorial w/ summary of model changes in ./src/Applications/LDAS_App/doc
 - output via MAPL History
 - functionality to create land boundary condtions (make_bcs)
 - functionality to create restart files (mk_restarts)
 - assimilation components do NOT yet work!!!
 - GEOSldas development history:
    - the following tags are nearly zero-diff:
        GEOSldas: GEOSldas_m4-16_0               (never released) 
	LDASsa:   reichle-LDASsa_m3-16_UNSTABLE  (approx. reichle-LDASsa_m3-16_5)
    - GEOSldas_m4-17_0 does not incl LDASsa changes after 1 Jun 2016!!
        (see LDASsa change log below)


***********************************************************************************
*                                                                                 *
*                             LDASsa change log                                   *
*                                                                                 *
***********************************************************************************

===================================================================
===                                                             ===
===                7Mar2019                                     ===
===                                                             ===
===================================================================

- added FP-5.22 transition (13 Mar 2019) to "cross-stream" dates 

===================================================================
===                                                             ===
===                9Jul2018                                     ===
===                                                             ===
===================================================================

- added FP-5.21 transition (11 Jul 2018) to "cross-stream" dates 

===================================================================
===                                                             ===
===               23Jan2018                                     ===
===                                                             ===
===================================================================

- slightly increased upward flow from rzexc to srfexc again by  
    changing "alpha" from 0.01 to 0.04  (MODEL CHANGE!)
- SMAP L1C_TB reader: 
  - removed stats check for L1C fore-minus-aft Tb differences!
  - use avg fore/aft timestamp so that fore and aft Tbs for same 
      location are never used in different assimilation windows

===================================================================
===                                                             ===
===               28Dec2017                                     ===
===                                                             ===
===================================================================

- added functionality to read and thin SMAP L1C TB E (Enhanced) obs
- added update_type = 10: 3d soil moisture excl catdef/Tskin/ght(1); TB obs


===================================================================
===                                                             ===
===               31Oct2017                                     ===
===                                                             ===
===================================================================

- added functionality for SRTM-based tiles (revised subroutine reorder_tiles())
- reduced upward flow from rzexc to srfexc (MODEL CHANGE!)
- bug fixes

===================================================================
===                                                             ===
===               30Oct2017                                     ===
===                                                             ===
===================================================================

- added FP-5.17 transition (1 Nov 2017) to "cross-stream" dates 


===================================================================
===                                                             ===
===                8Jun2017                                     ===
===                                                             ===
===================================================================

- added "%flistpath" and "%flistname" to "obs_param_type"
  for run-time input of obs file name lists


===================================================================
===                                                             ===
===               17Apr2017                                     ===
===                                                             ===
===================================================================

- changed protection against having obs that made it through reader but for
   which no tiles are found in get_obs_pred(), which happened when ellipse
   used in tile search straddled date line (fixes Tb_fcst=0 output)


===================================================================
===                                                             ===
===               19Jan2017                                     ===
===                                                             ===
===================================================================

- added FP-5.16 transition (24 Jan 2017) to "cross-stream" dates 


===================================================================
===                                                             ===
===               28Dec2016                                     ===
===                                                             ===
===================================================================

- added functionality for d5124 RPIT/FPIT forcing


===================================================================
===                                                             ===
===               08Jun2016                                     ===
===                                                             ===
===================================================================

- added SMAP & SMOS Tb obs and forecasts to ObsFcstAna and SMAP L4_SM aup output
   for obs that are not assimilated because they lack scaling params


===================================================================
===                                                             ===
===               01Jun2016                                     ===
===                                                             ===
===================================================================

- added stats check for fore-minus-aft Tb differences to SMAP L1C_TB reader


===================================================================
===                                                             ===
===               24Mar2016                                     ===
===                                                             ===
===================================================================

- new and more efficient work-around to make GEOS-5 forcing consistent 
   with LDASsa convention (speed-up by ~10% depending on configuration)
- update to Intel 15.0.3 and Open MPI 1.10.0
- bug fixes for reading perturbation std-dev from file
- added vectorized version "nr_ran2_2d()" 
- replaced Fortran "stop" with call to "ldas_abort()" in "enkf_general()"


===================================================================
===                                                             ===
===                3Feb2016                                     ===
===                                                             ===
===================================================================

- do *not* allow obs outside tile_grid (applies to all obs readers
   except outdated readers 'RedArkOSSE', 'VivianaOK', and 'tskin_ceop3n4')
- in ldsetup script removed requirement of N_ens=1 for "model" runs
- revised Aug 2014 cross-over date for FP forcing to Aug 20 at 0z 
- edited calculation of LAI-weighted greenness to avoid division by zero
- when mapping SMOS or SMAP EASEv2 gridded obs to tiles, apply requirement that
  obs and tile fall within same M36 grid cell only if tile space is EASEv2_M09
  or EASEv2_M36, but drop requirement for all other tile spaces.

===================================================================
===                                                             ===
===               23Dec2015                                     ===
===                                                             ===
===================================================================

- revised soil temperature that is used in mwRTM


===================================================================
===                                                             ===
===                7Dec2015                                     ===
===                                                             ===
===================================================================

- in SMAP "lmc" output set dztsurf=0.05 m everywhere (b/c of revised CSOIL_2)
- bug fix in date/time variables for stream transitions of MERRA-2 forcing
- use SWGDN from the MERRA-2 "rad" Collection instead of the "lfo" Collection


===================================================================
===                                                             ===
===               24Nov2015                                     ===
===                                                             ===
===================================================================

- use engineering fix for all veg types (dampen oscillations in off-line mode)
- LDASsa_DEFAULT_inputs_driver.nml: expanded lists of "commonly used" parameter values
- LDASsa_DEFAULT_inputs_ensupd.nml: updated SMOS paths to latest v620 obs
- added checks on parsing of MERRA, MERRA2, and G5DAS met_tag
- added protection against having obs that made it through reader but for
   which no tiles are found in get_obs_pred() (which may happen with certain 
   FOV and tile space settings)
- argument list change for catchment() to match AGCM revisions for debugging


===================================================================
===                                                             ===
===               20Nov2015                                     ===
===                                                             ===
===================================================================

- changed Stieglitz snow model parameter WEMIN back to 13 kg/m2 (MERRA value)
- changed Catchment model parameter CSOIL_2 back to 70000 J/K (pre-MERRA value)
- Catchment model "zbar" bug fix
- read vegetation height from boundary condition file
- added "veghght" to LDASsa "catparam" output file and to SMAP L4_SM "lmc" output
- revised turbulent roughness length (z0) formulation, new default formulation
- refined log messages for echo_pert_param()
- SMAP obs readers: added clean exit if h5 file to be read does not exist
- updated DAS utilities to "GEOSadas-5_13_1"

===================================================================
===                                                             ===
===                8Oct2015                                     ===
===                                                             ===
===================================================================

- minor MPI bug fix (avoid MPI crash under certain conditions)

===================================================================
===                                                             ===
===               22Sep2015                                     ===
===                                                             ===
===================================================================

- added MERRA2 forcing options
- revised FP "cross-stream" dates 
- bug fix in Princeton forcing routine

===================================================================
===                                                             ===
===                1Jul2015                                     ===
===                                                             ===
===================================================================

- new compiler (Intel 15.0.2.164; OpenMPI 1.8.4)
  (one of the few compiler/MPI combinations that worked ok with LDASsa derived types)
- bug fix in subroutine get_tile_num_from_latlon()  (needed for "LatLon" grid types)
- clarified units of LDASsa obs time stamp as J2000 seconds with "TT12" epoch
- added "monthsperjob" option to "ldsetup" script
- updated "cross-stream" dates for FP (added transition to GEOS-5.13.1 on 1 May 2015)
- updated comments on met forcing data in LDASsa_DEFAULT_inputs_driver.nml
- updated default values for SMAP Tb obs err std in LDASsa_DEFAULT_inputs_ensupd.nml


===================================================================
===                                                             ===
===                6May2015                                     ===
===                                                             ===
===================================================================

- added plumbing for ascending SMAP L2_SM_AP files


===================================================================
===                                                             ===
===               23Apr2015                                     ===
===                                                             ===
===================================================================

- Minor bug fix in subroutine read_obs_SMAP_halforbit_Tb().
- Minor bug fix in ldsetup.


===================================================================
===                                                             ===
===                2Apr2015                                     ===
===                                                             ===
===================================================================

- revised approach to field-of-view (FOV):
  - FOV in units of [km]:  get_obs_pred() uses Gaussian kernel out to 2*FOV
  - FOV in units of [deg]: constant averaging kernel out to 1*FOV
  - revisions also impact: (1) algorithm that determines which tile administers 
     an obs, (2) xhalo and yhalo used in get_obs_pred, (3) selection of obs
     that impact 1d FT analysis, (4) checks of xcorr vs. FOV, and (5) checks of 
     compact support lengths
- new "RTM_ID" field in obs_param_type to facilitate specific RTM configurations
- EASE grid only: screen Tb obs predictions for excessive water fraction
- constrain t2m, q2m between values at surface and in lowest atmospheric model layer
- new out_collection_ID=11 for SMOS pre-processing 
- added nine new soil parameters from most recent version of "soil_param.dat" 
   to cat_param_type; added comments with units of fields in cat_param_type
- added nml input "alb_from_SWnet" (switch for inferring albedo from SWnet
   instead from look-up table and MODIS albedo scaling factors),
   changed default to "alb_from_SWnet=.false."
- removed "ignore_SWnet_for_snow" which was meant to address a sibalb bug in the
   MERRA snow albedo but never worked as intended in the LDASsa MPI parallel 
   implementation because of a bug in LDASsa
- moved function shift_forcing_date() from clsm_ensdrv_functions.F90
   to clsm_ensdrv_force_routines.F90
- matlab reader for LDASsa perturbations restart files (read_pert_ldas_rst.m)
- by default only build MPI executables ("LDASsa_mpi.x" and "LDASsa_assim_mpi.x)


===================================================================
===                                                             ===
===                9Dec2014                                     ===
===                                                             ===
===================================================================

- clean-up of subroutine pmonth() (round-off differences)


===================================================================
===                                                             ===
===                5Dec2014                                     ===
===                                                             ===
===================================================================

- Ganymed-4_1 revisions of roughness length
- updated "cross-stream" transition dates for RP-IT, FP-IT, and FP


===================================================================
===                                                             ===
===               25Nov2014                                     ===
===                                                             ===
===================================================================

- change in CVS module to "LDASsa_m3"
- use GEOSsurface_GridComp as in Ganymed-4_0 but with oscillations fix and
    after clean-up of catchment.F90, catch_constants.f90, and StieglitzSnow.F90
- updated utilities to GEOSadas-5_13_0
  - IMPORTANT: values of some MAPL_Constants changed!
- clean-up of GEOSlana_GridComp/*.F90 files
  - split clsm_ensdrv_drv_routines.F90 into four files 
     ["drv", "init", "vegalb", "out_*]
  - make sure all modules include "private" statement
  - follow all "use" statements with "ONLY"
  - removed unused files (esat_qsat.F90, nr_sort.f, ldlnbcs)
  - removed unused variables
  - removed "domain_ID" (obsolete with MPI, originally meant to allow splitting
      global integrations into six entirely separate runs, one for each continent)
  - removed "out_wetness" from LDASsa nml inputs and hard-wired to "false"
      (applies only to collection_ID=1,2,3,10)
  - merged leap_year.F90 into date_time_util.F90
- completely revised ldsetup (python)
- if requested in general, write "incr" and/or "ObsFcstAna" files whenever it 
    was time for assimilation, even if there were no observations
- for SMAP "aup" Collection, changed units of soil moisture output
   from wetness [dimensionless] to volumetric soil moisture [m3/m3].
- explicitly select species that are appropriate for each update_type
- update type for stand-alone FT analysis (not yet tested)
- new obs species: freeze-thaw-fraction obs from SMAP L2_SM_AP granules
- renamed force_pert_type fields for consistency w/ met_force_type
- IMPORTANT: re-interpreted progn_pert as perturbation flux forcing!


===================================================================
===                                                             ===
===               15Oct2014                                     ===
===                                                             ===
===================================================================

- bug fix: fixed off-diagonal terms of obs error covariance for spatially variable 
           obs error variance in conjunction with spatial obs error correlations
- bug fix to address IO-related hangs:
   files opened by master process now only closed by master process
- avoid Fortran "system()" call in get_io_filename() to create directories;
   directories now pre-created by setup script


===================================================================
===                                                             ===
===               15Jul2014                                     ===
===                                                             ===
===================================================================

- disallow model time step greater than 7.5 minutes (model_dtstep>450s)
- clarified log messages about reading tavg lfo file for G5DAS file specs 
   when initializing forcing data
- cleaned up (non-critical) log messages re. EnKF analysis
- clarified log messages about number of observations
- check nml input path and file names in ldsetup (i.e., lenkf.pl) script
   to detect common job setup errors at the time of submission
- revised error handling (see "error codes" and subroutines
   ldas_abort(), ldas_warn() in ldas_exceptions.F90 


===================================================================
===                                                             ===
===               13Jun2014                                     ===
===                                                             ===
===================================================================

- engineering fix for surface energy balance oscillations in off-line mode
  - revised time interpolation for forcing "states"
  - new accounting term for energy balance
  - revised derivatives for turbulent flux terms


===================================================================
===                                                             ===
===               27May2014                                     ===
===                                                             ===
===================================================================

- optionally read perturbations std-dev from netcdf-4 file
- revised out_collection_ID=6 (SMAP "gph" Collection)     
- changed wilting point output in SMAP "lmc" Collection from wetness units 
   ("clsm_wpwet", [dimensionless]) to volumetric units ("clsm_wp", [m3 m-3])
- updated h5 element names in SMAP Tb reader again to work with *revised*
    GLOSIM3 simulated data


===================================================================
===                                                             ===
===               23May2014                                     ===
===                                                             ===
===================================================================

- bug fix: revised xcompact and ycompact parameters so that
           load-balanced analysis works with 1d update types


===================================================================
===                                                             ===
===                5May2014                                     ===
===                                                             ===
===================================================================

- fixed typo in matlab reader for SMAP L4_SM "lmc" file
- updated h5 element names in SMAP Tb reader to work with GLOSIM3 simulated data
  NOTE: This was only a temporary fix. GLOSIM3 files changed again ~22 May 2014.
        The tag "reichle-LDASsa_m2-SMAP_L4_SM_D00500_p3" created for Release 5
        on 7 May 2014 was DELETED again on 23 May 2014!


===================================================================
===                                                             ===
===                2Apr2014                                     ===
===                                                             ===
===================================================================

- revised parallelization of EnKF analysis for improved load balancing
- rewritten cat_enkf_increments() for update_type=2 (3d soil moisture)
   and update_type=7 (3d Tskin)
- removed "update_region" used with old update_type=2 (3d soil moisture)
- added t2m, q2m to catch_diagF type; new out_collection_ID=10
- revised obslog output: one file per run, nml variable to switch on/off
- apply all model-based QC only before EnKF update
- change default model-based QC of Tskin obs to also include frozen soil
- prior to assimilation, sort Observations_l by tilenum and then species 
   to fix lay-out dependency for MPI parallel execution
- added compilation flags "-g -traceback" by default
- added native SLURM directive for sponsor ID into ldsetup
- revisions to subroutine  read_LaRC_Tskin_nc4(), incl QC for solar zenith angle


===================================================================
===                                                             ===
===                6Mar2014                                     ===
===                                                             ===
===================================================================

- fixed obs_type derived type after addition of real*8 "time" field
  so that corresponding MPI STRUCT works with -align compiler flag

===================================================================
===                                                             ===
===                4Mar2014                                     ===
===                                                             ===
===================================================================

- added tpsn(1) to inst output for out_collection_ID=8

===================================================================
===                                                             ===
===               19Feb2014                                     ===
===                                                             ===
===================================================================

- added *.qa output for land-model-constants ("lmc") collection
- fixed tp units used for Tb output in collections 7,8 (SMAP Nature Run v03)

===================================================================
===                                                             ===
===               14Feb2014                                     ===
===                                                             ===
===================================================================

- added "time" field to Observations structure (obs_type)
- SMOS reader: make sure center-of-mass of tile that administers
   the obs is within EASEv2 M36 obs grid cell (discard obs otherwise)
- SMAP aup output: 
  - write only those obs that are assimilated
  - added output of "tb_[h/v]_obs_time_sec" and "tb_[h/v]_orbit_flag"
  - added super-obs computation for obs from more than one half-orbit
- new matlab function "write_smapL4SMqa.m" to generate *.qa files 
   from SMAP L4_SM "gph" or "aup" granules
- changed default Catchment "model_dtstep" to 7.5 min (450 s)
- disabled check in scale_obs_Tb_zscore() to avoid having to store redundant
   info in scaling files
- write "ldas_domdecomp.txt" output file to "rc_out/Yyyyy/Mmm/" 
   subdirectory instead of "rc_out/" (i.e., provide restart time tag)


===================================================================
===                                                             ===
===               10Feb2014                                     ===
===                                                             ===
===================================================================

- changed details of Wind time interpolation for G5DAS forcing 
- changed default surface turbulence scheme back to "Louis"
- updated transition dates for G5DAS RP/FP-IT forcing
- revised out_collection_IDs 7,8,9 for SMAP Nature Run v03 and
   mwRTM parameter calibration


===================================================================
===                                                             ===
===               24Jan2014                                     ===
===                                                             ===
===================================================================

- bilinear interpolation option for GEOS surface met forcing
- first version of the SMAP L1C_TB/L2_SM_AP Tb obs reader
- enable use of boundary condition files for higher-resolution *ocean*
  ("MERRA-2" 1/4 deg, "Ostia" 1/8 deg) in addition to "Reynolds" (1 deg)
- moved assemble_obs_cov() into cat_enkf_increments() for update_type=8
   (3d Tb analysis) to avoid running out of memory when assimilating
   Tbs from SMAP L2_SM_AP


===================================================================
===                                                             ===
===                8Jan2014                                     ===
===                                                             ===
===================================================================

- changed from ASCII to binary files for *tilecoord* and *tilegrid*
   domain files
- read latest-generation EASE-grid *.til file ok even if it is 
   *not* within a sub-dir named "SiB2_V2" 
- look for corrected precip data in "Yyyyy/Mmm/" and "Yyyyy/" dirs
- added output collection 8 (for mwRTM param calibration) with 
   different numbers of output variables for "tavg" and "inst" files
- added "cross-stream" capability for GEOS-5 RP-IT/FP-IT and FP forcing
- added "obslog" output files (SMOS only, not yet implemented for all readers!)
- limited the number of warnings in repair_forcing() to N_tile_warn_max
- cleaned up write statements into log file (avoid line breaks)


===================================================================
===                                                             ===
===               16Dec2013                                     ===
===                                                             ===
===================================================================

- revised max cloud fraction for QC of LaRC geostat Tskin retrievals
   from 5 to 20 percent in LaRC reader
- use FOV as the maximum distance allowed between obs lat/lon and tile
   com_lat/com_lon when searching for a tile to which the obs will be
   assigned
- added "fitted" SMOS Tb observations (angular fit); changed file names
  of preprocessed input files
- updated default paths to preprocessed SMOS soil moisture and Tb data
- control Intel 13 buffered I/O via "-assume buffered_io" compiler flag
  (as opposed to FORT_BUFFERED environment variable),
  and capture log output into "ldas_log" file via stdout
- added "progn_pert_type" to control prognostics perturbations
  (no longer use "cat_progn_type")
- added option to "coarsen" grid on which model, forcing, and observations
  perturbations are computed
- revised subroutine generate_white_field() to avoid if statement 
   within double loop
- added out_collection_ID=7 for SMAP Nature Run v03
  (write *different* output fields for "inst" and "tavg" output!)
- added stop if microwave radiative transfer (mwRTM) parameters are needed but missing
- bug fix in Tb observations predictions (get_obs_pred())
- allow output of SMAP L4_SM "aup" collection with "SMOS_fit_Tb*" obs species 
   (but *not* simultaneously with "SMAP_L*_Tb*" species!)
- bug fix: corrected MPI type used to MPI_BCAST integers "N_out_fields"
   and "out_collection_ID"
- bug fix: time interpolation of vegetation (GRN, LAI) and albedo


===================================================================
===                                                             ===
===               31Oct2013                                     ===
===                                                             ===
===================================================================

- split "cat_diagn" structure into "cat_diagS" and "cat_diagF"
- Changed how Catchment diagnostics ("cat_diagn") are recomputed
  after perturbations, EnKF updates, and bias corrections.
  The change also primarily how ensemble average diagnostics are output
  in "open loop" integrations (with perturbations, but without EnKF
  updates). 
- bug fix for instantaneous output of snow temperatures (tpsn)

===================================================================
===                                                             ===
===               26Sep2013                                     ===
===                                                             ===
===================================================================

- removed OpenMP (except for GNUmakefiles in GEOScatch_GridComp),
  renamed "omp_driver_routines.F90" to "ens_driver_routines.F90"
- switched to Intel 13 compiler on Discover
- minor optimization
- MKL library revisions (use MKL instead of NR ludcmp if available)
- MPI Barrier clean-up
- revised update_type=4
- obs bias restructuring and cleanup
- provide increments in LIAU format
- minor updates in "ldsetup" script
- added deallocate statements before exiting to avoid memory leak
  check failures


===================================================================
===                                                             ===
===               26Aug2013                                     ===
===                                                             ===
===================================================================

- major revisions for vegetation (LAI, GRN) and albedo scaling parameters
  for compatibility with "MAPL_readforcing" format use in Ganymed and
  flexible time stepping
- no longer read "tile_id" from *.til file; instead, assign IDs to
  tiles in the order in which they are read from the *.til file
  (ie, the order in which they are used in the GCM)
- added *.til file names from Ganymed-4 b.c. directories
- MERRA-2 "bug" fix for Catchment model "rzave" calculation
  (a.k.a. "c19" code change in soil parameter test bed)
- added possibility of using lat/lon of obs from reader
   (rather than using lat/lon from model tile_coord)
- in SMOS obs reader, temporarily shift lat/lon of obs for computation 
  of nearest tile to avoid ambiguous assignment of M09 model tile within 
  M36 obs grid cell


===================================================================
===                                                             ===
===               19Jul2013                                     ===
===                                                             ===
===================================================================

- added subroutine parse_G5DAS_met_tag() for using precip corrections
  with G5DAS forcing (as opposed to MERRA forcing)


===================================================================
===                                                             ===
===               21Jun2013                                     ===
===                                                             ===
===================================================================

- bug fix for SMOS 36 km EASE grid observations handling when
  model is on 9 km EASE tile space
- fixes for avoiding lay-out dependency of 3d analysis
- revised FFT to use Intel Math Kernel Library (MKL) if available
- overlapping rectangle bug fix in get_halo_tiles() and related fixes
- L1bas bug fix in reorder_tiles()
- added output of system date/time with each error message


===================================================================
===                                                             ===
===               29May2013                                     ===
===                                                             ===
===================================================================

- added 3d update type for brightness temperature assimilation
  (soil moisture, Tskin, and ght(1) increments)
- fixed get_obs_pert() for 3d updates and parallel perturbations
- added pert time steps (progn_pert_dtstep, force_pert_dtstep) to ens prop
  nml inputs (incl addition of basic checks for consistency of pert time steps


===================================================================
===                                                             ===
===               24May2013                                     ===
===                                                             ===
===================================================================

- MPI parallelization of perturbations
- switch MPI library to MVAPICH2 
  (link s/w development version of baselibs as defined in g5_modules)
- minor updates to Applications/LDAS_App/GNUmakefile


===================================================================
===                                                             ===
===                3May2013                                     ===
===                                                             ===
===================================================================

- added six SMAP Tb species (using SMOS reader for now)
- added output of SMAP L4_SM land-model-constants ("lmc") and 
  analysis update ("aup") collections in binary, tile-space format
- fixed default settings for GOES-E


===================================================================
===                                                             ===
===               11Apr2013                                     ===
===                                                             ===
===================================================================

- revised treatment of LDASsa output collections
  - added out_collection_ID=6 for SMAP L4_SM gph collection
- removed modis_alb_param_type fields "sc_albvr" and "sc_albnr"
  (these fields were not used and input parameters for them are
   no longer produced)
- added EASEv2 (version 2) grid  (cylindrical only)

===================================================================
===                                                             ===
===               27Mar2013                                     ===
===                                                             ===
===================================================================

- revised and updated model bias estimation
- updated utilities to CVS tag GEOSadas-5_9_1_p8
  - abandoned OpenMP libs for MAPL_Base and GEOS_Shared and therefore
    disabled OpenMP in Applications/LDAS_App/GNUmakefile
  - GEOScatch_GridComp was NOT updated at this time
    
===================================================================
===                                                             ===
===               15Feb2013                                     ===
===                                                             ===
===================================================================

- bug fix in get_obs_pert()
- added 3d Tskin update routines from Clara (update_type=7)
- bug fix in 3d soil moisture update (update_type=2)
- added recomputation of soil temperature to recompute_diagnostic()
  [important for Tb assimilation]
- model-based QC for Tsurf: added hard-coded option to eliminate frozen 
  conditions in QC based on tsurf and top layer soil temperature (tp1)
- bug fix: adjusted dimension in call to get_pert() in get_obs_pert()
  [important if a mix of assimilated and non-assimilated species is used]
- minor changes to improve behavior with "-check bounds"
- registered GEOS-5 1/4-by-5/16 deg DC grid in "read_til_file()"

===================================================================
===                                                             ===
===                7Feb2013                                     ===
===                                                             ===
===================================================================

- edited subroutine scale_obs_sfmc_cdf() to enable use of stats files 
   that do not perfectly match current domain (as in tile_coord)
- removed mwRTM parameter look-up table option (too complicated with 
    "new" (200+) soil classes)
- added warning (but no "stop"!) if all mwRTM parameters are no-data values
- added tile_id check when reading mwRTM params from file
- removed interception water from mwRTM 
- scaled surface soil moisture from Catchment model prior to input into mwRTM
- added subroutine for zscore scaling of Tb obs to model climatology
- fixed case N_rows=0 upon input in unique_rows()


===================================================================
===                                                             ===
===                6Jul2012                                     ===
===                                                             ===
===================================================================

- revised get_obs_pred() using halos for better memory management w/ MPI
- revised domain decomposition and addressed date line issue;
  added output of "*domdecomp.txt" file


===================================================================
===                                                             ===
===                3Apr2012                                     ===
===                                                             ===
===================================================================

- moved Catchment diagnostic routines to GEOScatch_GridComp/catchment.F90 and
  deleted catch_diagn_routines.F90
- moved subroutine check_catch_progn() to GEOScatch_GridComp/catch_iau.F90 and
  replaced with wrapper
- added utilities to convert from "cat_progn" type to regular arrays
- added "centered_update" option for assimilation (EnKF) updates


===================================================================
===                                                             ===
===               29Feb2012                                     ===
===                                                             ===
===================================================================

- always write restart file at end of time loop
- write appropriate "status" message to file prior to completion and if
   exiting through subroutine "stop_it()"
- clean exit from MPI if exiting through subroutine "stop_it()"
- revised output routines
  - switched from 49/58 MERRA-Land variables to 50/59 MERRA-Land variables
    (added TSURF as separate output field)
  - time-averaged snow temperatures weighted by snow cover fraction
  - set sub-tile temperatures to no-data-values when corresponding fraction
    is zero
  - range-check for FRWLT output
- revised "DAS" definitions in subroutine that reads GEOS-5 forcing
   (now called "get_GEOS()")
- parse "met_tag" to decide whether to use "MERRA_defs"
    (rather than check for presence of "diag_sfc" file)


===================================================================
===                                                             ===
===               28Dec2011                                     ===
===                                                             ===
===================================================================

- updated GEOScatch_GridComp to "Fortuna-2_5_p3" GCM tag
- updated GMAO_Shared to GEOSadas-5_7_2_p3
- revised treatment of water/energy balance and lai/greenness
  in preparation for MERRA-Land file specifications
  (new structures "bal_diagn_type" and "veg_param_type")
- added "MERRA-Land" file specs (major revision of output subroutines)
- added bit-shaving option for better gzip compression of output files
- adding "PARdrct" and "PARdffs" to met_force for MERRA-Land file specs
- removed "out_Pplus" (was read from "ensupd" nml file)
- simplified output of "totalb" (eliminated field from "cat_diagn" structure)
- accidental discovery and fix of bug in "ens_inst" output
- read_grid_elev(): corrected format of GEOS-5 gridded binary elevation file
- read_catchment_def(): fixed check for mismatch between tile_coord_file and
        catchment_def_file


===================================================================
===                                                             ===
===                1Dec2011                                     ===
===                                                             ===
===================================================================

- added QC for Tb based on model *soil* temp (motivated by RFI)
- skip computation of increments unless "assim" flag is true for at least
  one element of obs_param
- added matlab reader for "obsparam" file


===================================================================
===                                                             ===
===               23Nov2011                                     ===
===                                                             ===
===================================================================

- deleted "LDASsa_DAS_inputs_*.nml" files from current development tag 
  (reichle-LDASsa_m2-10_MPI_UNSTABLE_v2)
- removed "./LDAS_*" paths from LDASsa_DEFAULT_inputs_driver.nml
- read mwRTM parameters from file if mwRTM_param_path is not empty, use
  look-up table otherwise
- changed tsurf_threshold and precip_threshold b/c QC now done for 
  individual ensemble members rather than the ensemble mean


===================================================================
===                                                             ===
===               22Nov2011                                     ===
===                                                             ===
===================================================================

- renamed the scale_obs_*() subroutines in clsm_ensupd_read_obs.F90
  that handle the cdf-matching and std-normal-deviate scaling.
  IMPORTANT: The obs and the model data must be in the SAME UNITS *before* 
  invoking the matlab functions "get_model_and_obs_stats.m" and 
  "get_cdf_match.m" that generate the scaling files.
- added the field "units" to obs_param_type

===================================================================
===                                                             ===
===                4Oct2011                                     ===
===                                                             ===
===================================================================

- major overhaul for assimilation with MPI and downscaling
  - revised enkf_types (obs_param_type, obs_type)
      new fields, renamed fields
  - implemented obs "FOV" and "downscaling"
      obs no longer need match tiles
      revised subroutines get_obs_pred() and get_obs_pert()
      moved model-based QC from subroutine read_obs() to get_obs_pred()
  - output of observation-space data into "ObsFcstAna" file, new matlab reader 
        (replaces "innov" ("OminusF") and "OminusA" files)
  - new output file *ldas_obsparam.txt
  - deleted variables Obs_pred_minus and Obs_pred_plus (not needed)
  - new file src/Components/GEOSlana_GridComp/GNUmakefile_parallel
  - enabled compilation of GEOSlana_GridComp files with
     and without compile-time "LDAS_MPI" flag
- continued major changes for SMOS angles
  - revised handling of ens upd namelist inputs
      define one Tb pol/orbit type for all angles
      species id from ens upd nml now changes in subsequent processing
- bug fix re. call to calc_tp in get_obs_pred()
- minor bug fix for reading elevation from tile coord file
- minor bug fixes in mwRTM_routines
- added reader for SMOS soil moisture obs (as part of SMOS Tb reader)
- accomodate empty files in SMOS reader

===================================================================
===                                                             ===
===                1Jun2011                                     ===
===                                                             ===
===================================================================

- added low-frequency microwave radiative transfer model (mwRTM)
- added elevation to tile_coord structure (needed for mwRTM)
- added "ldas_mwRTMparam" output file and reader 
- output file format changed for *_tilecoord.txt
- added SMOS multi-angle Tb obs
- added possibility of no-data-values in Obs_pred and innov, OminusA
- change calc_tp to return tp in CELSIUS
- fixed no-data-value checks


===================================================================
===                                                             ===
===               13May2011                                     ===
===                                                             ===
===================================================================

- include EASE grid tools (map from i,j index to lat/lon and vice versa)
- revised mapping of observations to model grid
- updated to Fortuna-2_4_p2 (incl RDC bug fix and interface change of
  subroutine helfsurface())
- bug fix related to renaming of "UVA" to "LPRM"
- include subroutine calc_tp()
- updated utilities to GEOSadas-5_7_1


===================================================================
===                                                             ===
===               3Mar2011                                      ===
===                                                             ===
===================================================================

- updated to revised file format for corrected precip netcdf input files
  (now lon-by-lat, previously lat-by-lon)
- adjusted definition of MERRA "main"stream (when using MERRA forcing)


===================================================================
===                                                             ===
===              10Feb2011                                      ===
===                                                             ===
===================================================================

- added T2m, Q2m diagnostic for Louis and for Helfand-Monin-Obukhov schemes
- moved subroutine pmonth() into module propagate_cat (in file process_cat.F90)


===================================================================
===                                                             ===
===              12Jan2011                                      ===
===                                                             ===
===================================================================

- fixed bug in io_catch_internal_rstrt() -- previously tc4 was mistakenly
  initialized from tpsn1, and qc4 was mistakenly initialized from q over
  snow
- revised AMSR-E/NSIDC "Surface Type" QC flag
- added obs pert checks for AMSR-E/LPRM soil moisture retrievals
- replaced "UVA" with "LPRM"
- added prelim capability to read nc4 forcing files as generated by Fortuna GCM
- revised initialization of MERRA/DAS forcing (for improved restarting of
   off-line replay runs)


===================================================================
===                                                             ===
===              14Dec2010                                      ===
===                                                             ===
===================================================================

- renamed *default* namelist inputs files
- revamped treatment of input paths (resolution now appended in LDASsa)
- search for input parameter files within select sub-directories of input path
- further re-organization of ldsetup
  - added "source g5_modules" to ldsetup
- added "noopenmp" and "openmp" compilation for MAPL_Base and GEOS_Shared
  (required for OpenMP changes of 2 Dec 2010)


===================================================================
===                                                             ===
===              02Dec2010                                      ===
===                                                             ===
===================================================================

- updated GMAO_Shared to GEOSadas-5_6_2
- updated Catchment model (and surface turbulent boundary layer) 
   to Fortuna-2_3_UNSTABLE
- enabled compilation of sequential, OpenMP, and MPI versions of LDASsa
- fix snow albedo bug in MERRA SWNET forcing for MERRA-Land integrations
- added sfc_turb_scheme (choose Louis or Helfand Monin-Obukhov)
- added wesn perturbation also to htsn and sndz to make snow 
   variables consistent
- reorganization of catch_constants module and constants in subroutine 
   catchment()
- revisions to MERRA and Fortuna replay capability
- clean-up of "ldsetup" (and therefore of "lenkf.pl" run script)


===================================================================
===                                                             ===
===              29Oct2010                                      ===
===                                                             ===
===================================================================

- added soilcls30 and soilcls100
- optimized restart-to-exp-domain mapping in initialize_model()
- renamed N_gndtmp to N_gt
- adapted to new interface of subroutines catchment(), calc_soil_moist(),
   and partition()
  [new input arg "dzsf", for catchment() new optional outputs]



===================================================================
===                                                             ===
===              README.update_14May2010                        ===
===                                                             ===
===================================================================


Notes on updating LDASsa to "GEOSadas-5_5_2"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

reichle, 14 May 2010

-------------------------------------------
NOT updated (kept at "reichle-ldas-m2-08"):
-------------------------------------------

  src/GMAO_Shared/GEOS_Shared/*  

    "GEOS_Shared" as tagged in "reichle-ldas-m2-08" is equivalent 
    to MERRA/GEOSdas-2_1_4 GCM tag (but was done differently in MERRA GCM tag)

    do NOT update for compatibility with LDASsa and MERRA

  src/GMAO_Shared/MAPL_Base/MAPL_Constants.F90

    do NOT update for compatibility with MERRA

  src/Components/GEOScatch_GridComp  

    do NOT update for compatibility with MERRA

----------------------------------------------------
Additional changes from GEOSadas-5_5_2:
----------------------------------------------------

For hdf4, nc3, and nc4 compatibility, use "Baselibs-3_1_9" at 

  /discover/nobackup/projects/gmao/share/dao_ops/Baselibs/v3.1.9_build1

After checkout, change basedir in in src/g5_modules to:

    set basedir = /discover/nobackup/projects/gmao/share/dao_ops/Baselibs/v3.1.9_build1


    
        


===================================================================
===                                                             ===
===              README.update_4Dec2009                         ===
===                                                             ===
===================================================================



Notes on updating LDASsa to "GEOSadas-5_4_0_p4_UNSTABLE"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

reichle, 4 Dec 2009

-------------------------------------------
NOT updated (kept at "reichle-ldas-m2-08"):
-------------------------------------------

  src/GMAO_Shared/GEOS_Shared/*  

    "GEOS_Shared" as tagged in "reichle-ldas-m2-08" is equivalent 
    to MERRA/GEOSdas-2_1_4 GCM tag (but was done differently in MERRA GCM tag)

    do NOT update for compatibility with LDASsa and MERRA

  src/GMAO_Shared/MAPL_Base/MAPL_Constants.F90

    do NOT update for compatibility with MERRA

  src/Components/GEOScatch_GridComp  

    do NOT update for compatibility with MERRA

----------------------------------------------------
Additional changes from GEOSadas-5_4_0_p4_UNSTABLE":
----------------------------------------------------

For hdf4, nc3, and nc4 compatibility, use "Baselibs-3_1_7" at 

  /discover/nobackup/dnadeau/basedir/Baselibs-3_1_7/x86_64-unknown-linux-gnu/ifort/

This requires the following changes for LDASsa (need to be done manually!!!):

 In src/Config/ESMA_base.mk change
 FROM:
    ...
    DIR_NETCDF = $(BASEDIR)/$(ARCH)
    INC_NETCDF = $(DIR_NETCDF)/include/netcdf
    LIB_NETCDF = $(BASELIB)/libnetcdf.a $(LIB_HDF5)
    ...
 TO:
    DIR_NETCDF = $(BASEDIR)/$(ARCH)
    INC_NETCDF = $(shell $(DIR_NETCDF)/bin/nc-config --includedir)/netcdf
    LIB_NETCDF = $(shell $(DIR_NETCDF)/bin/nc-config --libs)

 In src/g5_modules change basedir to:
    set basedir = /discover/nobackup/dnadeau/basedir/Baselibs-3_1_7/x86_64-unknown-linux-gnu/ifort






===================================================================
===                                                             ===
===              README_CLSM-update-Sep-2007.txt                ===
===                                                             ===
===================================================================



List of changes to update CLSM in LDASsa to latest GEOS5 version

Sep 10, 2007 - reichle


Update to latest GEOS5 version of CLSM is based on Sarith's changes
to an early version of LDASsa in 
SARITHPATH=/discover/home/sarith/Model_Comparison/Catch/src_geos5_v1.9/
(based on halem:~/reichle/CLSM_ensdrv/V3/)

Catchment model files (catchment.F90, catch_constants.f90, sibalb_coeff.f90)
are now read directly from GEOScatch_GridComp

1.) 
CVS tag for LDASsa just before updating CLSM to latest GEOS5 version:
tag=reichle-ldas-before-CLSM-update

2.) 
Added new type "modis_alb_type" to driver_types.F90:
- based on SARITHPATH/modis_alb_types.f90

3.) 
Added read statements for MODIS alb param to subroutine read_land_parameters()
- based on SARITHPATH/clsm_ensdrv_modis_routines.f90
- updated call statement for read_land_parameters() in clsm_ensdrv_main.f90
  (along with appropriate variable declarations)

4.)
Added "modis_alb_param_path" to subroutine read_driver_inputs() 

5.)
Added "modis_alb_param_path" to namelist file clsm_ensdrv_default_inputs.nml

6.) 
Added time interpolation of modis_alb_param to interpolate_to_timestep()
(includes changes to call statement)

7.) 
Deleted albed_ntp from clsm_ensdrv_main.F90, omp_driver_routines.F90,
process_cat.F90

8.) 
Moved computation of albedo from subroutine interpolate_to_timestep() to 
subroutine propagate_cat()
(includes changes to call statement)

9.) 
In subroutine read_land_parameters() changed re-setting of dzrz back to 
earlier version for consistency with Sarith's Catchment parameters

10.)
Replaced catchment.f with catchment.F90 from GEOScatch_GridComp
NOTE: repair "quick fix" that is still missing in GEOScatch_GridComp version
(see "reichle" entry dated 30 July 2007)
NOTE: make calc_soil_moist() and get_tf_nd() public
Add "use catchment_model" to catch_diagn_routines.F90, process_cat.F90, 
clsm_ensupd_upd_routines.F90, and clsm_ensdrv_drv_routines.F90.

11.)
Added wtot, wchange, etot, echange, hsnacc to catch_diagn_type
(also for operators)
Added same to call statement for subroutine catchment()
Added tile_coord%tile_id to call statement for subroutine catchment()

12.)
Remove sibalb.f (subroutine sibalb() is now in catchment.F90)
Add sibalb_coeff.f90

13.) 
Replace pmonth_stage4.f with Sarith/LIS pmonth_stage4.f90
NOTE: Now the roughness length is actually provided as an output
from subroutine pmonth() and used in turb calc's
(previously, the roughness length from the lookup table
that was read in read_land_parameters() was used)
Delete vgd_ntp, zol_ntp, vgd_lookup, zol_lookup

14.)
Added catch_constants.f90
Use values from module catch_constants throughout the driver
via "use catch_constants, ONLY: ..." statements
(including clsm_ensdrv_drv_routines.F90, catch_types.F90, 
clsm_ensdrv_pert_routines.F90, process_cat.F90, esat_qsat.F90, 
catch_diagn_routines.F90, clsm_ensdrv_glob_param.F90)

15.) 
Moved subroutine turb() from turb.f into module propagate_cat()
in process_cat.F90

16.) 
Changes to CVS module "LDASsa" so that Catchment model files are
read directly from GEOScatch_GridComp.   Note that a special
GNUmakefile is used in GEOScatch_GridComp that does NOT compile
GEOScatch_GridComp.F90 - consequently, ESMF and other GMAO shared
infrastructure is NOT needed to compile LDASsa.
***********************************************************************
***  VERY IMPORTANT!!! ************************************************
*                                                                     *
* When updating GEOScatch_GridComp from CVS, must make corresponding  *
* modification to GEOScatch_GridComp/GNUmakefile                      *
*                                                                     *
***********************************************************************

17.)
Removed all files in GMAO_Shared/GMAO_cfio because cfio is not used 
in LDASsa.




===================================================================
===                                                             ===
===              release_notes.txt                              ===
===                                                             ===
===================================================================



V1.1: file format change : The (exp_run)_tile_coord.dat file now contains 
                           minlon, maxlon, minlat, maxlat for each tile.
                           Note the change in the format of the file.
                           - reichle, 19 Aug 2005


V2:    bug fix: restart for smaller domain from larger domain
                fixed bug in subroutine initialize_model()
       bug fix: added "if" statement in subroutine
                get_tile_num_from_latlon() in tile_coord.f90
                (takes care of some special cases)
       GSWP2 grid: added "pole on center" check for grids in 
                   subroutine read_tile_coord()
       Assimilation: included all files/subroutines for assimilation,
                     make sure to "#define ASSIM" in clsm_ensdrv_main.f90
                     and to set "ASSIM = 1" in Makefile
                     Tskin innovations should work, Tskin assimilation not
                     ready
           
       - reichle, 27 Sep 2005


V3:  includes Tskin scaling, assimilation, and (model) bias estimation
     changed file names of default namelist files
  
     - reichle, 8 Nov 2005

V3.1: - fixed minor bug in tile_coord.f90 (done 2005/11/17)
      - added minimal "inst" output capability
      - add RUNTIMELIMIT to bsub options in run_job_global.sh 
        for "background" queue

     - reichle, 6 Dec 2005




===================================================================
===                                                             ===
===              version_info.txt                               ===
===                                                             ===
===================================================================

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tag:  ldas-1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
branch:  stassi-ldas-1 (20051201)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tag:  stassi-ldas-2 (20051205)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


SUMMARY
=======
In this version, the baseline test is run in a separate "work"
directory rather than in the checkout directory.  Also, namelist and
data inputs are linked into the work directory rather than searching
for them in other directories.  A setup script, ldsetup, creates the
run script, lenkf.pl,

In this version, lenkf.pl, is hard-coded with the baseline inputs and
is not modified by the inputs entered through the ldsetup script.


ADDITIONS
=========
The following namelist files were copied from
~reichle/CLSM_ensdrv/V3/etc and added to the CVS repository:

 * ens_upd_gswp2_scale_assim.nml
 * ens_prop_fpert_gswp.nml
 * driver_inputs_gswp2_1by1.nml

(NOTE: These are NON-DEFAULT namelist files specified in the run
script.)

These files were needed for the initial test run.  It looks like they
were in the Attic, which would mean that they were previously in the
repository.

The following files were added to this version:

 * ldlnbcs : script for linking input to local directory

 * ldsetup : script for setting up the work directory and writing the
   run script, ldlnbcs script, and namelist files there for testing.


MODIFICATIONS
=============
The following files were modified in order to make the initial
baseline run:

 * driver_inputs_gswp2_1b1.nml
   ---------------------------
    The file was modified to get data from the
    /share/todling/fvInput/ldas/DE/ directory instead of from the
    /land/l_data/geos5/bcs/SiB2_V2/DE/ directory.  The /land/...
    data were being modified and could not provide stable data for
    creating baseline results.

   (NOTE: New files are more correct.  Copy to /share/todling/ and 
   re-create baseline.)

 * run_job_global.sh
   -----------------

   The following modifications were made in order to be able to run
   the baseline case from the test directory:
   - modified WORKDIR location to point to /u1/stassi/ldas/stassi-ld1/src
   - modified EXPPATH location to point to /u1/stassi/scratch
   - modified the SPONSOR (charge number)
   - added ETCDIR variable definition for convenience
   - create ETCDIR if necessary and copy namelist files there

   NOTE: Anyone wishing to duplicate the baseline test will need to
         repeat the first two modifications above, except pointing to
         their own directories.  The SPONSOR value may also need
         changing.

In the following namelist and source code files, the input directory
location was changed to local ('./') as much as possible:

 * clsm_ensdrv_default_inputs.nml
 * driver_inputs_gswp2_1by1.nml
 * clsm_bias_routines.f90
 * clsm_ensdrv_drv_routines.f90
 * clsm_ensdrv_pert_routines.f90
 * clsm_ensupd_upd_routines.f90

(NOTE: LAI_GRNyy directories all have "green.dat" files for year "yy".
Check with gcm, where the data come from.)

(NOTE: Time invariant files, e.g. ar.new and bf.dat, in GCM these
parameters are put into restart files?)

(NOTE: Look at how the restart files are being handled.)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tag:  stassi-ldas-3 (20051206)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUMMARY
=======
This version is basically the same as version stassi-ldas-2, except in
this version, the values input by ldsetup get put into the lenkf.pl
script.

Note: In order to run the test case, the following values need to be
modified in the lenkf.pl script created by ldsetup using default
responses:

my $ENDMONTH     = 6;    ->    my $ENDMONTH     = 9;
my $ENDHOUR      = 6;    ->    my $ENDHOUR      = 0;


MODIFICATIONS
=============
The following data path modification was made to reflect a
restructuring of the input directories:

 /land/l_data/ISCCP/GSWP2_grid_V1/DX1/tskin_d32/
 ->  /land/l_data/ISCCP/GSWP2_1by1_V1/

in the following two namelist files:

 * clsm_ensupd_default_inputs.nml
 * ens_upd_gswp2_scale_assim.nml


Other modifications were made:

 * ldsetup
  -------
  - $expid variable added
  - changed sub set_defaults() to sub default(), which acts more like
    a function, returning individual default values, rather than
    setting all default values globally.
  - replace hard-coded WORKDIR value with $fvwork
  - replace START and END times with variables
  - replace hard-coded path of EXECUTABLE with variable
  - replace hard-coded value of $EXPRUN with $expid
  - add sub twoDigit() and use this sub to simplify the creation of
    $STARTDATESTRING

(NOTE: switch [optional flag?] in ldsetup, either specify end date/time
or duration.)

(NOTE: need to add DOMAIN, perhaps passed down from controlling
script.)

(NOTE: set up to not prompt if all the input parameters are supplied.)

(NOTE: for testing on halem, can use background queue instead of
gmao_long. No charge for jobs sent to background.  Only one node )


 * version_info.txt
   ----------------
   Updated with the latest version information.


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tag:  stassi-ldas-4/5 (20051206)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUMMARY
=======



MODIFICATIONS
=============

 * Makefile : Removed references to ../etc directory; also removed
      trailing blanks.

 * Makefile.conf.OSF1 : Changed compile flag: "-g" -> "-O4"

 * ldsetup
   - Added $SUBMITDRIVERINPUTS and $SUBMITENSPROPINPUTS to runscript
   - Set all $SUBMIT.. variables to 0 (zero)

 * clsm_ensdrv_default_pert_inputs.nml : replaced with ens_prop_fpert_gswp.nml

 * ldlnbcs
   - replaced value of LDASBCS with a variable that gets set during
     ldsetup

 * ldsetup

   - prompt for LDAS input directory location
   - prompt for location of initial restart files
   - sed the ldlnbcs file to set the value for LDASBCS
   - set value of restart directory in the lenkf.x file
   - add $SUBMITDRIVERINPUTS and $SUBENSPROPINPUTS to lenkf.x script


REMOVED
=======
The following non-default namelist files were removed from the repository

 * driver_inputs_gswp2_1by1.nml
 * ens_prop_fpert_gswp.nml
 * ens_upd_gswp2_scale_assim.nml

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tag:  stassi-ldas-6 (20060125)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


SUMMARY
=======
Modifications were made to run the LDAS code on palm.


ADDITIONS
=========
More input data were added to the /share/todling/fvInput/ldas/ directory:
 * GSWP-2/ 
 * GSWP2_1by1_V1/
 * stats/


MODIFICATIONS
=============

 * clsm_ensdrv_default_inputs.nml: met_path value changed for new
      location of the GSWP-2/Baseline_Forcings/ data directory.

 * clsm_ensupd_default_inputs.nml:
   - obs_param(3)%path value changed for new location of GSWP2_1by1_V1
     data
   - obs_param(3)%scalepath value changed for new location of stats
     data

 * GNUmakefile
   - Added targets for "scripts" and "namelists"
   - Modified the "perl" target

 * ldlnbcs
   - Changed "#!/usr/bin/csh -f" -> "#!/bin/csh -f"
   - Changed variable name: LDASBCS -> LDASINP
   - Added entries for the input data directories: GSWP-2,
     GSWP2_1by1_V1, and stats

 * ldsetup
   - Added machine dependencies for halem vs. palm, specifically the
     fvhome and the rstrtdir directory locations, and the #BSUB and
     #PBS directives in the job file.
   - Replaced hard-coded "sed" command with "@DASSED"
 
 * Makefile.conf.Linux
   - Used many of the same options in the Makefile.conf.OSF1 file






===================================================================
===                                                             ===
===              README                                         ===
===                                                             ===
===================================================================


README for land EnKF driver

- reichle, 8 Nov 2005

Location of "frozen" version prior to inclusion in CVS:
halem:/u1/reichle/CLSM_ensdrv/V3/

All paths below are on halem.

See also "src/release_notes.txt".


Source code:
---------------------------
 For complete list of source files see "Makefile".
 
 "catchment.f" should match same from GEOS5.

 Parallelization: 
 OpenMP directives in 
 - clsm_ensdrv_main.f90
 - omp_driver_routines.f90
 - clsm_ensupd_upd_routines.f90
 Need benchmark for "1/2 deg tile space", MPI parallelization
 might be required.


Build:
---------------------------
 Makefile
 configure.sh -> Makefile.conf

 Variables "ASSIM" and "BIAS" in Makefile must be set to 1,
 see also corresponding cpp directives in clsm_ensdrv_main.f90.


"Global" fortran parameter files:
------------------------------------------------------  
 clsm_ensdrv_glob_param.f90       
 clsm_ensupd_glob_param.f90


Run script:
---------------------------
 run_job_global.sh

 
Hierarchy of command-line and namelist inputs:
------------------------------------------------------ 
 1. Command line arguments are used.
 2. If not available, namelist inputs from "special" namelist files (path
    and file name specified as command line argument) are used.
 3. If still not available, inputs from default namelist files are used.

 Note that "final" namelist files are written to the output directory for each run.


Namelist files:
---------------------------
 Default namelist files: 

 reading *default* driver inputs from ./../etc//clsm_ensdrv_default_inputs.nml
 reading *default* bias inputs from ./../etc//clsm_bias_default_inputs.nml
 reading *default* EnKF inputs from ./../etc//clsm_ensupd_default_inputs.nml
 reading *default* ens prop inputs from ./../etc//clsm_ensdrv_default_pert_inputs.nml

 For each of these, a "special" namelist path/file can be specified at the command
 line.  See "Hierarchy" above and "etc/" directory for sample files.

 "driver" inputs:        start time, end time, model parameter file paths, ...
 "ens prop/pert" inputs: std, correlation scales of perturbations to meteorological forcing 
                         data and prognostic variables
 "bias" inputs:          bias parameters
 "ens upd" inputs:       assimilation time step, anything related to observations


Restart files:
---------------------------
 run1.ensXXXX.rstrt.bkg.19860601_0000  (Catchment model prognostic variables, first time
                                        around only ensemble member 0 is required)

 run1.ensXXXX.pert.rstrt.19860601_0000 (ensemble perturbation restart inputs,
                                        not needed first time around)

 run1.ens_avg.rstrt.bias.19860601_0000 (bias restart inputs,
                                        not needed first time around)


Domain/grid/tile-space input files:
------------------------------------------------------
 Paths are specified in clsm_driver_inputs.nml

 Example files below are for special 1 deg-by-1 deg "GSWP" grid.  For 1/2 deg MERRA
 replace ".../DE/GSWP2_1by1/..." with ".../DC/FV_576x361/..."

 tile_coord file (same as GEOS5 (even though the land group maintain their own 
 directories in /land/l_data)!

 /land/l_data/geos5/bcs/SiB2_V2/DE/GSWP2_1by1//FV_360x180_DE_360x180_DE.til

 Catchment definitions file (additional inputs not available from above tile_coord file
 - this file is NOT part of GEOS5 but should be...):	

 /land/l_data/geos5/bcs/SiB2_V2/DE/GSWP2_1by1//catchment.def

 "Blacklist/whitelist" files are used for definition of smaller domains (not global)
 and to exclude tiles/catchments for which forcing data are not available.
 They should not be necessary for coupled integrations - GEOS5-DAS should provide
 complete global coverage of forcing data.

 Note that "final" "tile_coord.dat" and "domain.dat" files are written to the output 
 directory of each run.


Catchment model parameter input files (soil, vegetation, etc):
---------------------------------------------------------------------------------
 Paths are specified in clsm_driver_inputs.nml - all files are same as GEOS5 (even though 
 the land group maintain their own directories in /land/l_data/...)

 /land/l_data/geos5/bcs/SiB2_V2/DE/GSWP2_1by1/VEGETATION-GSWP2/mosaic_veg_typs_fracs
     
 /land/l_data/geos5/bcs/SiB2_V2/DE/GSWP2_1by1/VEGETATION-GSWP2/LAI_GRN_CLIMATOLOGY/green.dat 
 /land/l_data/geos5/bcs/SiB2_V2/DE/GSWP2_1by1/VEGETATION-GSWP2/LAI_GRN_CLIMATOLOGY/lai.dat
 
 /land/l_data/geos5/bcs/SiB2_V2/DE/GSWP2_1by1/VEGETATION-GSWP2/FIXED/zol.dat
 /land/l_data/geos5/bcs/SiB2_V2/DE/GSWP2_1by1/VEGETATION-GSWP2/FIXED/vgd.dat
 
 /land/l_data/geos5/bcs/SiB2_V2/DE/GSWP2_1by1//soil_param.dat
 /land/l_data/geos5/bcs/SiB2_V2/DE/GSWP2_1by1//tau_param.dat
 /land/l_data/geos5/bcs/SiB2_V2/DE/GSWP2_1by1//ar.new
 /land/l_data/geos5/bcs/SiB2_V2/DE/GSWP2_1by1//bf.dat
 /land/l_data/geos5/bcs/SiB2_V2/DE/GSWP2_1by1//ts.dat

 Note that "final" binary "cat_param.dat" file is written to the output 
 directory of each run.


Meteorological forcing input files:
-------------------------------------
  Path to surface meteorological forcing files ("met_path") is specified via namelist.
  Need to point to GEOS5-DAS "inst2d" and "tavg2d" output files.

  /land/l_data/GSWP-2/Baseline_Forcings/


Observations input files:
---------------------------
  Path to observations and related scaling files is specified via namelist.

  (File names for observations are not echoed to screen at this time.)

  Example scaling file:

  /land/reichle/NSIPP/catch/output/Tskin_sarith/run1tmp/GLOBAL_GSWP-2/stats/run1tmp.mean_std.1986-1995.Jun_03z.bin
  
  ```
  
  [Unreleased]: https://github.com/GEOS-ESM/GEOSldas/tree/develop
  [v17.8.0]: https://github.com/GEOS-ESM/GEOSldas/releases/tag/v17.8.0
  [v17.9.0-beta.0]:  https://github.com/GEOS-ESM/GEOSldas/releases/tag/v17.9.0-beta.0
  
