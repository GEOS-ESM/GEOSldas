#!/usr/bin/env python
#
# module load python/GEOSpyD/Ana2019.03_py3.7
#
import os
import glob
import subprocess as sp

heads = """
VERSION: 1
EXPID:  GEOSldas_expid
EXPDSC: GEOSldas_output
EXPSRC: GEOSldas

COLLECTIONS:
"""

label = """
::
GRID_LABELS: PC720x361-DC
::
PC720x361-DC.GRID_TYPE: LatLon
PC720x361-DC.IM_WORLD: 720
PC720x361-DC.JM_WORLD: 361
PC720x361-DC.POLE: PC
PC720x361-DC.DATELINE: DC
PC720x361-DC.LM: 1
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

nens = 32
with open('HISTORY.rc', 'w') as f:
    f.write(heads)
    collection, body = hist_template.split('::\n')
    collection = collection.strip('\n').strip("'")  
    for i in range(nens):
       i = i+1
       ids = collection+f'{i:04}'
       f.write(f"'{ids}'\n")
    f.write(label)
    lines = body.split('\n')
    for i in range(nens):
       i = i+1
       collect= collection+f'{i:04}.'
       for line in lines:
          newline = line
          if ":" in line :
             newline = collect+line
          if "CATCHINCR_e" in newline:
             newline = newline.replace('CATCHINCR_e',f'CATCHINCR_e{i:04}')
          f.write(newline+'\n')
       f.write('::\n') 
