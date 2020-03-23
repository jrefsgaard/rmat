# FindLOGFT.cmake
#
# Finds the LOGFT library
#
# This will define the following variables
#
#    LOGFT_FOUND
#    LOGFT_INCLUDE_DIRS
#    LOGFT_LIBRARIES
#
# and the following imported targets
#
#     LOGFT::LOGFT  
#
# Author: Jonas Refsgaard (jonas.refsgaard@kuleuven.be)

find_package(PkgConfig)
pkg_check_modules(PC_LOGFT QUIET LOGFT)

find_path(LOGFT_INCLUDE_DIR
  NAMES logft/logft.h
  PATHS ${PC_LOGFT_INCLUDE_DIRS}
)
find_library(LOGFT_LIBRARY
  NAMES logft
  PATHS ${PC_LOGFT_LIBRARY_DIRS}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LOGFT
  FOUND_VAR LOGFT_FOUND
  REQUIRED_VARS
    LOGFT_LIBRARY
    LOGFT_INCLUDE_DIR
)

#The old-school way
if(LOGFT_FOUND)
  set(LOGFT_LIBRARIES ${LOGFT_LIBRARY})
  set(LOGFT_INCLUDE_DIRS ${LOGFT_INCLUDE_DIR})
  set(LOGFT_DEFINITIONS ${PC_LOGFT_CFLAGS_OTHER})
endif()

#The modern way
if(LOGFT_FOUND AND NOT TARGET LOGFT::LOGFT)
  add_library(LOGFT::LOGFT UNKNOWN IMPORTED)
  set_target_properties(LOGFT::LOGFT PROPERTIES
    IMPORTED_LOCATION "${LOGFT_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_LOGFT_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${LOGFT_INCLUDE_DIR}"
  )
endif()

mark_as_advanced(
  LOGFT_INCLUDE_DIR
  LOGFT_LIBRARY
)

