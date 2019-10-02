add_subdirectory (@MAPL)
add_subdirectory (@GMAO_Shared)

# Special case - FMS is built twice with two
add_subdirectory (@FMS fms_r4)
add_subdirectory (@FMS fms_r8)
