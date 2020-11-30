esma_set_this ()

set (SRCS
  GEOS_OutputIncrGridComp.F90 
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL
  )

target_compile_definitions (${this} PRIVATE LDAS_MPI)
