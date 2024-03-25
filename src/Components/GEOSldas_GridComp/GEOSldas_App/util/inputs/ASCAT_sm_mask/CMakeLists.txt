# build without installation

add_executable(ascat_mask_maker.x ascat_mask_maker.F90)
target_link_libraries(ascat_mask_maker.x MAPL)
