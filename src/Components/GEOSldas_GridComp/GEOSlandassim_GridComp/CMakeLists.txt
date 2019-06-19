esma_set_this ()

set (SRCS
  io_hdf5.F90 enkf_general.F90 adapt_types.F90 mwRTM_types.F90 clsm_ensupd_glob_param.F90
  mwRTM_routines.F90 clsm_ensupd_upd_routines.F90 clsm_ensdrv_drv_routines.F90
  clsm_ensupd_read_obs.F90 catch_bias_types.F90 clsm_bias_routines.F90 clsm_adapt_routines.F90
  clsm_ensupd_enkf_update.F90 clsm_ensdrv_out_routines.F90 GEOS_LandAssimGridComp.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES GEOS_LdasShared GEOSens_GridComp GEOSlandpert_GridComp GEOSland_GridComp MAPL_Base
  INCLUDES ${INC_ESMF} ${INC_HDF5})

target_compile_definitions (${this} PRIVATE LDAS_MPI)
