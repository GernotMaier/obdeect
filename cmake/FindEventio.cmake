# Locate an installed EventIO C library with the CORSIKA IACT photon API, or
# build the adjacent EventIO source checkout when no binary is installed.
if(NOT Eventio_SOURCE_ROOT AND NOT Eventio_ROOT AND NOT Eventio_LIBRARY
   AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../hessio/eventio_basic/eventio.c")
  set(Eventio_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../hessio")
endif()

if(Eventio_SOURCE_ROOT AND EXISTS "${Eventio_SOURCE_ROOT}/corsika/io_simtel.c")
  enable_language(C)
  add_library(obdeect_eventio_c STATIC
    "${Eventio_SOURCE_ROOT}/eventio_basic/eventio.c"
    "${Eventio_SOURCE_ROOT}/corsika/io_simtel.c"
    "${Eventio_SOURCE_ROOT}/corsika/mc_atmprof.c"
    "${Eventio_SOURCE_ROOT}/common/warning.c"
    "${Eventio_SOURCE_ROOT}/common/fileopen.c"
    "${Eventio_SOURCE_ROOT}/common/straux.c")
  target_include_directories(obdeect_eventio_c PUBLIC
    "${Eventio_SOURCE_ROOT}/common"
    "${Eventio_SOURCE_ROOT}/eventio_basic"
    "${Eventio_SOURCE_ROOT}/eventio_addon"
    "${Eventio_SOURCE_ROOT}/corsika")
  set_target_properties(obdeect_eventio_c PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED YES)
  add_library(Eventio::Eventio ALIAS obdeect_eventio_c)
  set(Eventio_FOUND TRUE)
  return()
endif()

find_path(Eventio_CORSIKA_INCLUDE_DIR mc_tel.h
  HINTS ${Eventio_ROOT} PATH_SUFFIXES include corsika)
find_path(Eventio_BASIC_INCLUDE_DIR io_basic.h
  HINTS ${Eventio_ROOT} PATH_SUFFIXES include eventio_basic)
find_path(Eventio_COMMON_INCLUDE_DIR initial.h
  HINTS ${Eventio_ROOT} PATH_SUFFIXES include common)
find_library(Eventio_LIBRARY NAMES hessio eventio
  HINTS ${Eventio_ROOT} PATH_SUFFIXES lib lib64)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eventio REQUIRED_VARS Eventio_LIBRARY
  Eventio_CORSIKA_INCLUDE_DIR Eventio_BASIC_INCLUDE_DIR Eventio_COMMON_INCLUDE_DIR)
if(Eventio_FOUND AND NOT TARGET Eventio::Eventio)
  add_library(Eventio::Eventio UNKNOWN IMPORTED)
  set_target_properties(Eventio::Eventio PROPERTIES
    IMPORTED_LOCATION "${Eventio_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES
      "${Eventio_COMMON_INCLUDE_DIR};${Eventio_BASIC_INCLUDE_DIR};${Eventio_CORSIKA_INCLUDE_DIR}")
endif()
