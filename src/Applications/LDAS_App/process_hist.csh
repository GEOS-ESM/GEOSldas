#!/bin/csh -f

## I am changed the CUBE/EASE logic
## if CUBE we produce 2D
## anything else, SMAP and other offline grids we produce tile space

setenv    LSM_CHOICE $1
setenv    AEROSOL_DEPOSITION $2
setenv    GRID $3
setenv    GRIDNAME $4
setenv    HISTRC $5
setenv    RUN_IRRIG $6
setenv    ASSIM  $7
setenv    NENS $8

echo $GRIDNAME

if($ASSIM == 1) then
   sed -i 's|\#ASSIM|''|g' $HISTRC
   sed -i '/^\#EASE/d' $HISTRC
   sed -i '/^\#CUBE/d' $HISTRC
else
   sed -i '/^\#ASSIM/d' $HISTRC
endif

if($GRID == CUBE) then
   sed -i '/^\#EASE/d' $HISTRC
   sed -i 's|\#CUBE|''|g' $HISTRC
   sed -i 's|GRIDNAME|'"$GRIDNAME"'|g' $HISTRC
else
   sed -i '/^\#CUBE/d' $HISTRC
   sed -i 's|\#EASE|''|g' $HISTRC
   sed -i 's|GRIDNAME|'"$GRIDNAME"'|g' $HISTRC
endif

if($LSM_CHOICE == 1) then
   set GridComp = CATCH
   sed -i '/^>>>HIST_CATCHCN<<</d' $HISTRC
   sed -i '/^>>>HIST_CATCHCNCLM45<<</d' $HISTRC
endif

if($LSM_CHOICE == 2) then
   set GridComp = CATCHCN
   sed -i '/^>>>HIST_CATCHCNCLM45<<</d' $HISTRC
   sed -i 's/>>>HIST_CATCHCN<<</''/g' $HISTRC
endif

if($LSM_CHOICE == 3) then
   set GridComp = CATCHCN
   sed -i 's/>>>HIST_CATCHCN<<</''/g' $HISTRC
   sed -i 's/>>>HIST_CATCHCNCLM45<<</''/g' $HISTRC
endif

if($NENS > 1) then
   set GridComp = ENSAVG
   sed -i 's|VEGDYN|'VEGDYN_e0000'|g' $HISTRC
   sed -i 's|TP1|'TSOIL1TILE'|g' $HISTRC
   sed -i 's|TP2|'TSOIL2TILE'|g' $HISTRC
   sed -i 's|TP3|'TSOIL3TILE'|g' $HISTRC
   sed -i 's|TP4|'TSOIL4TILE'|g' $HISTRC
   sed -i 's|TP5|'TSOIL5TILE'|g' $HISTRC
   sed -i 's|TP6|'TSOIL6TILE'|g' $HISTRC
#   sed -i 's|DATAATM|'DATAATM0000'|g' $HISTRC
endif

sed -i 's|GridComp|'$GridComp'|g' $HISTRC

if($AEROSOL_DEPOSITION == 0) then
   sed -i '/^>>>HIST_AEROSOL<<</d' $HISTRC
else
   sed -i 's/>>>HIST_AEROSOL<<</''/g' $HISTRC
endif

if($RUN_IRRIG == 0) then
   sed -i '/^>>>HIST_IRRIG<<</d' $HISTRC
else
   sed -i 's/>>>HIST_IRRIG<<</''/g' $HISTRC
endif
