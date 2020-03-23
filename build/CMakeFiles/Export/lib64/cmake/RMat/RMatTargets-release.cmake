#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "RMat::RMat" for configuration "Release"
set_property(TARGET RMat::RMat APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(RMat::RMat PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "ROOT::TreePlayer"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libRMat.so"
  IMPORTED_SONAME_RELEASE "libRMat.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS RMat::RMat )
list(APPEND _IMPORT_CHECK_FILES_FOR_RMat::RMat "${_IMPORT_PREFIX}/lib64/libRMat.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
