
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
   SOURCES util/inputs/mwRTM_params/mwrtm_bin2nc4.F90
   LIBS GEOSlandassim_GridComp)
   
set (scripts
   process_hist.csh 
   remap_config_ldas.py 
   ens_forcing/average_ensemble_forcing.py 
   ens_forcing/ensemble_forc.py 
   ens_forcing/regrid_forc.csh 
   ens_forcing/enpert_forc.csh
   util/config/rewind_GEOSldas.csh
   )

install (
   PROGRAMS ${scripts}
   DESTINATION bin
   )

set(file ldas_setup)
configure_file(${file} ${file} @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${file} DESTINATION bin)

file(GLOB rc_files  GEOSldas_*rc)
file(GLOB nml_files LDASsa_DEFAULT*nml)

install(
   FILES ${rc_files} ${nml_files} lenkf.j.template
   DESTINATION etc
   )

esma_add_subdirectories(util/inputs/ASCAT_sm_mask)   
