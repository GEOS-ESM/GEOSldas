esma_set_this ()

set (alldirs
  GEOSexportcatchincr_GridComp 
  )

set (SRCS
  io_hdf5.F90 enkf_general.F90 adapt_types.F90 mwRTM_types.F90 clsm_ensupd_glob_param.F90
  mwRTM_routines.F90 clsm_ensupd_upd_routines.F90 clsm_ensdrv_drv_routines.F90
  clsm_ensupd_read_obs.F90 catch_bias_types.F90 clsm_bias_routines.F90 clsm_adapt_routines.F90
  clsm_ensupd_enkf_update.F90 clsm_ensdrv_out_routines.F90 GEOS_LandAssimGridComp.F90
  )

find_package(HDF5 REQUIRED COMPONENTS Fortran)

esma_add_library (${this}
  SRCS ${SRCS}
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES GEOS_LdasShared GEOSens_GridComp GEOSlandpert_GridComp GEOSland_GridComp makebcs MAPL NCEP_bufr_r4i4 GMAO_gfio_r4 hdf5hl_fortran hdf5_fortran ${NETCDF_LIBRARIES}
  INCLUDES ${INC_ESMF} ${INC_HDF5})

target_compile_definitions (${this} PRIVATE LDAS_MPI)
