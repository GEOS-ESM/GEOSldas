esma_set_this ()

set (SRCS
  nr_ran2_gasdev.F90 nr_jacobi.F90 nr_fft.F90
  random_fields.F90 land_pert.F90 force_and_cat_progn_pert_types.F90 LDAS_PertRoutines.F90
  GEOS_LandPertGridComp.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES GEOS_LdasShared GEOSens_GridComp GEOSland_GridComp MAPL ${MKL_LIBRARIES}
  INCLUDES ${INC_ESMF} ${MKL_INCLUDE_DIRS})

if (NOT CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
   target_compile_definitions(${this} PRIVATE MKL_AVAILABLE)
endif ()
