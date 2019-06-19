esma_set_this ()

set (SRCS
  LDAS_HashTable.F90 LDAS_Forcing.F90 LDAS_Interp.F90 GEOS_MetforceGridComp.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES GEOSsurface_GridComp GEOSland_GridComp GEOS_LdasShared MAPL_Base
  INCLUDES ${INC_ESMF})
