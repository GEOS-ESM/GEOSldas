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
"""

hist_template = """
'catch_progn_incr'
::
descr:       'Tile-space,3-Hourly,Instantaneous,Single-Level,Assimilation, Land Prognostics Increments',
template:    '%y4%m2%d2_%h2%n2z.bin',
mode:        'instantaneous',
frequency:   030000,
ref_time:    000000,
fields:     'TCFSAT_INCR'    , 'CATCHINCR'      ,
            'TCFTRN_INCR'    , 'CATCHINCR'      ,
            'TCFWLT_INCR'    , 'CATCHINCR'      ,
            'QCFSAT_INCR'    , 'CATCHINCR'      ,
            'QCFTRN_INCR'    , 'CATCHINCR'      ,
            'QCFWLT_INCR'    , 'CATCHINCR'      ,
            'CAPAC_INCR'     , 'CATCHINCR'      ,
            'CATDEF_INCR'    , 'CATCHINCR'      ,
            'RZEXC_INCR'     , 'CATCHINCR'      ,
            'SRFEXC_INCR'    , 'CATCHINCR'      ,
            'GHTCNT1_INCR'   , 'CATCHINCR'      ,
            'GHTCNT2_INCR'   , 'CATCHINCR'      ,
            'GHTCNT3_INCR'   , 'CATCHINCR'      ,
            'GHTCNT4_INCR'   , 'CATCHINCR'      ,
            'GHTCNT5_INCR'   , 'CATCHINCR'      ,
            'GHTCNT6_INCR'   , 'CATCHINCR'      ,
            'WESNN1_INCR'    , 'CATCHINCR'      ,
            'WESNN2_INCR'    , 'CATCHINCR'      ,
            'WESNN3_INCR'    , 'CATCHINCR'      ,
            'HTSNNN1_INCR'   , 'CATCHINCR'      ,
            'HTSNNN2_INCR'   , 'CATCHINCR'      ,
            'HTSNNN3_INCR'   , 'CATCHINCR'      ,
            'SNDZN1_INCR'    , 'CATCHINCR'      ,
            'SNDZN2_INCR'    , 'CATCHINCR'      ,
            'SNDZN3_INCR'    , 'CATCHINCR'      ,
"""

nens = 32
with open('HISTORY.rc', 'w') as f:
    f.write(heads)
    collection, body = hist_template.split('::\n')
    collection = collection.strip('\n').strip("'")  
    for i in range(nens):
       ids = collection+f'{i:04}'
       f.write(f"'{ids}'\n")
    f.write(label)
    lines = body.split('\n')
    for i in range(nens):
       collect= collection+f'{i:04}.'
       for line in lines:
          newline = line
          if ":" in line :
             newline = collect+line
          if "CATCHINCR" in newline:
             newline = newline.replace('CATCHINCR',f'CATCHINCR{i:04}')
          f.write(newline+'\n')
       f.write('::\n') 
