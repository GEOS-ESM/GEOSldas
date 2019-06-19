add_executable (preprocess_ldas.x preprocess_ldas.F90)
target_link_libraries (preprocess_ldas.x GEOSldas_GridComp)

add_executable (tile_bin2nc4.x tile_bin2nc4.F90)
target_link_libraries (tile_bin2nc4.x GEOSldas_GridComp)

add_executable (mwrtm_bin2nc4.x mwrtm_bin2nc4.F90)
target_link_libraries (mwrtm_bin2nc4.x GEOSldas_GridComp)

