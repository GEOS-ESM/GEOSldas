esma_set_this ()

esma_add_subdirectories (@GEOSgcm_GridComp @FVdycoreCubed_GridComp)

set (alldirs
  GEOSmetforce_GridComp
  GEOSlandpert_GridComp
  GEOSens_GridComp
  GEOSlandassim_GridComp
  )

esma_add_library(${this}
  SRCS GEOS_LdasGridComp.F90
  SUBCOMPONENTS ${alldirs}
  SUBDIRS Shared
  DEPENDENCIES GEOSland_GridComp FVdycoreCubed_GridComp MAPL_Base
  INCLUDES ${INC_ESMF})
