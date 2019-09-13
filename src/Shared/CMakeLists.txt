add_subdirectory (@MAPL)
add_subdirectory (@GMAO_Shared)

# Special case - GFDL_fms is built twice with two
add_subdirectory (@GFDL_fms GFDL_fms_r4)
add_subdirectory (@GFDL_fms GFDL_fms_r8)
