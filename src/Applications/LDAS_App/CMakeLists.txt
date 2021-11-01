
ecbuild_add_executable (
   TARGET GEOSldas.x 
   SOURCES GEOSldas.F90
   LIBS GEOSldas_GridComp MAPL)

ecbuild_add_executable (
   TARGET preprocess_ldas.x 
   SOURCES preprocess_ldas.F90 preprocess_ldas_routines.F90
   LIBS GEOSldas_GridComp MAPL)

ecbuild_add_executable (
   TARGET tile_bin2nc4.x 
   SOURCES tile_bin2nc4.F90
   LIBS MAPL)

ecbuild_add_executable (
   TARGET mwrtm_bin2nc4.x
   SOURCES mwrtm_bin2nc4.F90
   LIBS GEOSlandassim_GridComp)
   
install(
   PROGRAMS ldas_setup process_hist.csh process_rst.csh util/average_ensemble_forcing.py
   DESTINATION bin
   )

file(GLOB rc_files GEOSldas_*rc)
file(GLOB nml_files LDASsa_DEFAULT*nml)
install(
   FILES ${rc_files} ${nml_files} lenkf.j.template
   DESTINATION etc
   )
