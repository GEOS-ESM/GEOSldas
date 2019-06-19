cmake_minimum_required (VERSION 3.13)
cmake_policy (SET CMP0048 NEW)  # VERSION in project()
cmake_policy (SET CMP0053 NEW)
cmake_policy (SET CMP0054 NEW)

project (
  ESMA
  VERSION 5.3
  LANGUAGES Fortran CXX C)  # Note - CXX is required for ESMF

set(geos_configuration "gcm" CACHE STRING "GEOS Configuration to build." FORCE)
set_property(CACHE geos_configuration PROPERTY STRINGS gcm adas)


enable_testing()

if ("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are disabled. Please
           issue cmake command in separate build directory.")
endif ("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

set (DOING_GEOS5 YES)

# Should find a better place for this - used in Chem component
set (ACG_FLAGS -v)

set (F2PYEXT .so)
set (F2PY_SUFFIX .so)

set (FV_PRECISION R4)

# Set the default build type to release with debug info
if (NOT CMAKE_BUILD_TYPE)
  message (STATUS "Setting build type to 'Debug' as none was specified.")
  set (CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/@cmake")
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/@cmake/check_compiler_support")

include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()

include ("${CMAKE_Fortran_COMPILER_ID}")


# include various macros that are used in the build process
include (check_fortran_support)

include (esma_check_if_debug)
include (esma_set_this)
include (esma_add_subdirectories)
include (esma_add_library)
include (esma_generate_automatic_code)
include (esma_create_stub_component)
include (esma_fortran_generator_list)

include (UseProTeX)
set (protex_flags -g -b -f)

set (LATEX_COMPILER pdflatex)
include (UseLatex)


if (APPLE)
  include(osx_extras)
endif ()


# OpenMP support
find_package (OpenMP)

# Threading support
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

# Position independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Generic DFLAGS
add_definitions(-Dsys${CMAKE_SYSTEM_NAME} -DESMA64)


# MPI Support - only invoked from Fortran sources in GEOS-5.
# But some BASEDIR packages use MPI from C/C++.
find_package (MPI REQUIRED)

if (APPLE)
  if (DEFINED ENV{MKLROOT})
    set (MKL_Fortran)
    find_package (MKL REQUIRED)
  else ()
    if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
      #USE FRAMEWORK
      message(STATUS "Found macOS and gfortran, using framework Accelerate")
      link_libraries("-framework Accelerate")
    endif ()
  endif ()
else ()
  find_package (MKL REQUIRED)
endif ()

add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})

# Unit testing
set (PFUNIT OFF CACHE BOOL "Activate pfunit based tests")
if (PFUNIT)
   add_custom_target(tests COMMAND ${CMAKE_CTEST_COMMAND})
endif ()

# Baselibs ...
include (FindBaselibs)


# Recursively build subdirectories
# List subdirectories
add_subdirectory (src)


add_custom_target (doc
  COMMAND
  )

# Need to copy various files into build directory.  Source tree cannot know where the build tree is ...
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g5_modules ${CMAKE_CURRENT_BINARY_DIR}/g5_modules COPYONLY)
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g5_modules ${CMAKE_CURRENT_BINARY_DIR}/bin/g5_modules COPYONLY)
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g5_modules.sh ${CMAKE_CURRENT_BINARY_DIR}/g5_modules.sh COPYONLY)
