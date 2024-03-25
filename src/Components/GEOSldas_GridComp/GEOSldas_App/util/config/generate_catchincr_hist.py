#!/usr/bin/env python
#
# This code works with the python version loaded by g5_modules associated with GEOSldas v17.13.1:
#    python/GEOSpyD/Min4.11.0_py3.9_AND_Min4.8.3_py2.7
#
# This script generates a sample HISTORY.rc file for GEOSldas to write Catchment
#    model analysis increments in ensemble space, as needed in the weakly-coupled
#    Hybrid-4DEnVar land-atm DAS (LADAS).

import os
import glob
import subprocess as sp

# ------------------------------------------------------------------
#
# specify number of ensemble members here:

nens = 32

# ------------------------------------------------------------------
#
# some definitions of text elements in HISTORY.rc file

heads = """
#
#  Sample GEOSldas HISTORY.rc file for LADAS (atm ensemble)
#
#  - This sample HISTORY.rc is for the GEOSldas instance that is weakly coupled with the
#       atmospheric ensemble component of the Hybrid-4DEnVar ADAS (ADASens).
#
#  - The sample file was generated with the utility script
#       "GEOSldas/src/Components/GEOSldas_GridComp/GEOSldas_App/util/config/generate_catchincr_hist.py".
#
#  - The sample file triggers output of the GEOSldas "catch_progn_incr" collection in
#       ensemble space, which is needed by ADASens.
#
#  - The IDs of the ensemble members and their total number in GEOSldas must match
#       those of ADASens.
#
#  - The "catch_progn_incr" output is in tile space, which must be the same for
#       GEOSldas and ADASens.
#
#
##################################################################################

EXPID:  MyGEOSldasAtmEns

COLLECTIONS:
"""

label = """
::
"""

hist_template = """
'catch_progn_incr'
::
descr:       'Tile-space,3-Hourly,Instantaneous,Single-Level,Assimilation, Land Prognostics Increments',
template:    '%y4%m2%d2_%h2%n2z.bin',
mode:        'instantaneous',
frequency:   030000,
ref_time:    013000,
fields:     'TCFSAT_INCR'    , 'CATCHINCR_e'      ,
            'TCFTRN_INCR'    , 'CATCHINCR_e'      ,
            'TCFWLT_INCR'    , 'CATCHINCR_e'      ,
            'QCFSAT_INCR'    , 'CATCHINCR_e'      ,
            'QCFTRN_INCR'    , 'CATCHINCR_e'      ,
            'QCFWLT_INCR'    , 'CATCHINCR_e'      ,
            'CAPAC_INCR'     , 'CATCHINCR_e'      ,
            'CATDEF_INCR'    , 'CATCHINCR_e'      ,
            'RZEXC_INCR'     , 'CATCHINCR_e'      ,
            'SRFEXC_INCR'    , 'CATCHINCR_e'      ,
            'GHTCNT1_INCR'   , 'CATCHINCR_e'      ,
            'GHTCNT2_INCR'   , 'CATCHINCR_e'      ,
            'GHTCNT3_INCR'   , 'CATCHINCR_e'      ,
            'GHTCNT4_INCR'   , 'CATCHINCR_e'      ,
            'GHTCNT5_INCR'   , 'CATCHINCR_e'      ,
            'GHTCNT6_INCR'   , 'CATCHINCR_e'      ,
            'WESNN1_INCR'    , 'CATCHINCR_e'      ,
            'WESNN2_INCR'    , 'CATCHINCR_e'      ,
            'WESNN3_INCR'    , 'CATCHINCR_e'      ,
            'HTSNNN1_INCR'   , 'CATCHINCR_e'      ,
            'HTSNNN2_INCR'   , 'CATCHINCR_e'      ,
            'HTSNNN3_INCR'   , 'CATCHINCR_e'      ,
            'SNDZN1_INCR'    , 'CATCHINCR_e'      ,
            'SNDZN2_INCR'    , 'CATCHINCR_e'      ,
            'SNDZN3_INCR'    , 'CATCHINCR_e'      ,
"""

# ------------------------------------------------------------------
#
# write file "HISTORY.rc" with nens "catch_progn_incrXXXX" collections,
#   one for each ensemble member

with open('HISTORY.rc', 'w') as f:
    f.write(heads)
    collection, body = hist_template.split('::\n')
    collection = collection.strip('\n').strip("'")
    for i in range(nens):
       i = i +1
       sfx = '%04d'%(i)
       ids = "'"+collection+sfx+"'"
       f.write(ids+'\n')
    f.write(label)
    lines = body.split('\n')
    for i in range(nens):
       i = i+1
       sfx = '%04d'%(i)
       collect= collection+sfx+"."
       for line in lines:
          newline = line
          if ":" in line :
             newline = collect+line
          if "CATCHINCR_e" in newline:
             sfx = '%04d'%(i)
             frep = 'CATCHINCR_e'+sfx
             newline = newline.replace('CATCHINCR_e',frep)
          f.write(newline+'\n')
       f.write('::\n')

# ====================== EOF =====================================
