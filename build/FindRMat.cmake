# FindRmat.cmake
#
# Finds the Rmat library
#
# This will define the following variables
#
#    RMat_FOUND
#    RMat_INCLUDE_DIRS
#    RMat_LIBRARIES
#
# and the following imported targets
#
#     RMat::RMat  
#
# Author: Jonas Refsgaard (jonas.refsgaard@kuleuven.be)

find_package(PkgConfig)
pkg_check_modules(PC_RMat QUIET RMat)

find_path(RMat_INCLUDE_DIR
  NAMES rmat.h
  PATHS ${PC_RMat_INCLUDE_DIRS}
  PATH_SUFFIXES rmat
)
find_library(RMat_LIBRARY
  NAMES RMat
  PATHS ${PC_RMat_LIBRARY_DIRS}
)

set(RMat_VERSION ${PC_RMat_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(RMat
  FOUND_VAR RMat_FOUND
  REQUIRED_VARS
    RMat_LIBRARY
    RMat_INCLUDE_DIR
  VERSION_VAR RMat_VERSION
)

#The old-school way
if(RMat_FOUND)
  set(RMat_LIBRARIES ${RMat_LIBRARY})
  set(RMat_INCLUDE_DIRS ${RMat_INCLUDE_DIR})
  set(RMat_DEFINITIONS ${PC_RMat_CFLAGS_OTHER})
endif()

#The modern way
if(RMat_FOUND AND NOT TARGET RMat::RMat)
  add_library(RMat::RMat UNKNOWN IMPORTED)
  set_target_properties(RMat::RMat PROPERTIES
    IMPORTED_LOCATION "${RMat_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_RMat_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${RMat_INCLUDE_DIR}"
  )
endif()

mark_as_advanced(
  RMat_INCLUDE_DIR
  RMat_LIBRARY
)

