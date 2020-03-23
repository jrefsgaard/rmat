# FindCPC.cmake
#
# Finds the CPC library
#
# This will define the following variables
#
#    CPC_FOUND
#    CPC_INCLUDE_DIRS
#    CPC_LIBRARIES
#
# and the following imported targets
#
#     CPC::CPC  
#
# Author: Jonas Refsgaard (jonas.refsgaard@kuleuven.be)

find_package(PkgConfig)
pkg_check_modules(PC_CPC QUIET CPC)

find_path(CPC_INCLUDE_DIR
  NAMES cpc/cwfcomp.h
  PATHS ${PC_CPC_INCLUDE_DIRS}
)
find_library(CPC_LIBRARY
  NAMES CPC
  PATHS ${PC_CPC_LIBRARY_DIRS}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CPC
  FOUND_VAR CPC_FOUND
  REQUIRED_VARS
    CPC_LIBRARY
    CPC_INCLUDE_DIR
)

#The old-school way
if(CPC_FOUND)
  set(CPC_LIBRARIES ${CPC_LIBRARY})
  set(CPC_INCLUDE_DIRS ${CPC_INCLUDE_DIR})
  set(CPC_DEFINITIONS ${PC_CPC_CFLAGS_OTHER})
endif()

#The modern way
if(CPC_FOUND AND NOT TARGET CPC::CPC)
  add_library(CPC::CPC UNKNOWN IMPORTED)
  set_target_properties(CPC::CPC PROPERTIES
    IMPORTED_LOCATION "${CPC_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_CPC_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${CPC_INCLUDE_DIR}"
  )
endif()

mark_as_advanced(
  CPC_INCLUDE_DIR
  CPC_LIBRARY
)
