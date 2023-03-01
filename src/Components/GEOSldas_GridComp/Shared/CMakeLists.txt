set (this GEOS_LdasShared)

set (SRCS
  enkf_types.F90 catch_types.F90 LDAS_ensdrv_Globals.F90 LDAS_DateTimeMod.F90 LDAS_DriverTypes.F90
  LDAS_Convert.F90 LDAS_Exceptions.F90 LDAS_TileCoordType.F90 LDAS_PertTypes.F90
  LDAS_ensdrv_functions.F90 my_matrix_functions.F90
  LDAS_TileCoordRoutines.F90
  LDAS_RepairForcing.F90
  LDAS_ensdrv_mpi.F90
  )

list (APPEND SRCS
  nr_indexx.f my_lu_decomp.f
  )

esma_add_library(${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL  GEOS_Shared GEOS_LandShared makebcs
  INCLUDES ${INC_ESMF} ${INC_NETCDF})
