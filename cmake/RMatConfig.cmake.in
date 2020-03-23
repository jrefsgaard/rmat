get_filename_component(RMat_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(CMakeFindDependencyMacro)
list(APPEND CMAKE_MODULE_PATH ${RMat_CMAKE_DIR})

#find_dependency (Boost 1.64 REQUIRED)
find_dependency (GSL REQUIRED)
if(NOT TARGET ROOT::Libraries)
  #This is dirty, but hell breaks loose if FindROOT is run two times.
find_dependency (ROOT REQUIRED)
endif()
find_dependency (Armadillo REQUIRED)
find_dependency (OpenMP)
find_dependency (AUSALIB REQUIRED)
find_dependency (LOGFT REQUIRED)
find_dependency (CPC REQUIRED)

if(NOT TARGET RMat::RMat)
    include("${RMat_CMAKE_DIR}/RMatTargets.cmake")
endif()
