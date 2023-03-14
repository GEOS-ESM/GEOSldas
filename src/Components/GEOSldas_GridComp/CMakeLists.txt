esma_set_this ()

esma_add_subdirectory(GEOSgcm_GridComp )

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
  DEPENDENCIES GEOSland_GridComp makebcs MAPL
  INCLUDES ${INC_ESMF})
