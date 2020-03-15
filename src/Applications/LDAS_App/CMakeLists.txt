ecbuild_add_executable (
   TARGET GEOSldas.x 
   SOURCES GEOSldas.F90
   LIBS GEOSldas_GridComp MAPL)

set(executables 
   preprocess_ldas 
   tile_bin2nc4 
   mwrtm_bin2nc4
   mk_GEOSldasRestarts)

foreach (prog ${executables})
   ecbuild_add_executable (
      TARGET ${prog}.x 
      SOURCES ${prog}.F90
      LIBS GEOSldas_GridComp mk_restarts)
endforeach ()

install(
   PROGRAMS ldas_setup process_hist.csh process_rst.csh write_nightlytest_execfiles GEOSldas_nightlytest
   DESTINATION bin
   )

file(GLOB rc_files GEOSldas_*rc)
file(GLOB nml_files LDASsa_DEFAULT*nml)
install(
   FILES ${rc_files} ${nml_files} lenkf.j.template
   DESTINATION etc
   )
