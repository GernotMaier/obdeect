# Optional external library; support both installed headers and source layouts.
find_path(Hessio_CORSIKA_INCLUDE_DIR mc_tel.h
  HINTS ${Hessio_ROOT} PATH_SUFFIXES include corsika)
find_path(Hessio_EVENTIO_INCLUDE_DIR io_basic.h
  HINTS ${Hessio_ROOT} PATH_SUFFIXES include eventio_basic)
find_path(Hessio_COMMON_INCLUDE_DIR initial.h
  HINTS ${Hessio_ROOT} PATH_SUFFIXES include common)
find_library(Hessio_LIBRARY NAMES hessio
  HINTS ${Hessio_ROOT} PATH_SUFFIXES lib lib64)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Hessio REQUIRED_VARS Hessio_LIBRARY
  Hessio_CORSIKA_INCLUDE_DIR Hessio_EVENTIO_INCLUDE_DIR Hessio_COMMON_INCLUDE_DIR)
if(Hessio_FOUND AND NOT TARGET Hessio::Hessio)
  add_library(Hessio::Hessio UNKNOWN IMPORTED)
  set_target_properties(Hessio::Hessio PROPERTIES
    IMPORTED_LOCATION "${Hessio_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES
      "${Hessio_COMMON_INCLUDE_DIR};${Hessio_EVENTIO_INCLUDE_DIR};${Hessio_CORSIKA_INCLUDE_DIR}")
endif()
