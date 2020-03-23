# Tries to determine the location of AUSAlib.
# It will try to look in the systems include/lib paths
# Furthermore you can custimize it with the following environment variables
#
#  AUSALIB_INCLUDE_DIRS  Include directory
#  AUSALIB_LIBRARIES     Directory containing the build library

find_package(PkgConfig)
pkg_check_modules(PC_Rmat QUIET Rmat)

find_path(AUSALIB_INCLUDE_DIR
  NAMES ausa/AUSA.h
  PATHS ${PC_AUSALIB_INCLUDE_DIRS}
)
find_library(AUSALIB_LIBRARY
  NAMES AUSA
  PATHS ${PC_AUSALIB_LIBRARY_DIRS}
)

set(AUSALIB_VERSION ${PC_Rmat_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AUSALIB
  FOUND_VAR AUSALIB_FOUND
  REQUIRED_VARS
    AUSALIB_LIBRARY
    AUSALIB_INCLUDE_DIR
  VERSION_VAR AUSALIB_VERSION
)

#The old-school way
if(AUSALIB_FOUND)
  set(AUSALIB_LIBRARIES ${AUSALIB_LIBRARY})
  set(AUSALIB_INCLUDE_DIRS ${AUSALIB_INCLUDE_DIR})
  set(AUSALIB_DEFINITIONS ${PC_AUSALIB_CFLAGS_OTHER})
endif()

#The modern way
if(AUSALIB_FOUND AND NOT TARGET AUSALIB::AUSALIB)
  add_library(AUSALIB::AUSALIB UNKNOWN IMPORTED)
  set_target_properties(AUSALIB::AUSALIB PROPERTIES
    IMPORTED_LOCATION "${AUSALIB_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_AUSALIB_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${AUSALIB_INCLUDE_DIR}"
  )
endif()

mark_as_advanced(
  AUSALIB_INCLUDE_DIR
  AUSALIB_LIBRARY
)
