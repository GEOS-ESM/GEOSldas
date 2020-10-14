esma_set_this ()

esma_add_library(${this}
  SRCS GEOS_EnsGridComp.F90
  DEPENDENCIES GEOSland_GridComp GEOS_LdasShared MAPL
  INCLUDES ${INC_ESMF})
