cmake_minimum_required(VERSION 3.8)
project(libRMat VERSION 1.0 LANGUAGES CXX)

# Tell cmake where to look for special cmake modules.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#Specify default build type (This is important for optimisation)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")
endif()

#Make sure the library finds its dependencies, even after installation.
################################
#Make sure to have the full runtime-path to dependencies that may be installed in non-standard folders.
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")
################################

# Find files to compile
file(GLOB_RECURSE HEADERS "include/*")
file(GLOB_RECURSE SOURCE    "src/*.cpp")
add_library(RMat SHARED ${HEADERS} ${SOURCE})

target_include_directories(RMat
  PUBLIC 
    $<INSTALL_INTERFACE:include>    
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

#We use modern features.
target_compile_features(RMat PUBLIC cxx_std_11)

#Probably we only need Boost math-module, but there is a bug in FindBOOST.cmake
#find_package (Boost 1.64 REQUIRED)

find_package (GSL REQUIRED)

#New Find-module providing the imported target ROOT::Libraries.
#  ROOT::Libraries has the most common targets:
#    Core RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics 
#    MathCore Thread MultiProc
find_package (ROOT REQUIRED)

#I modified the FindArmadillo.cmake to provide an imported target Armadillo::Armadillo
find_package (Armadillo REQUIRED)

find_package (OpenMP)

#Modified Find-module to provide imported target AUSALIB::AUSALIB
find_package (AUSALIB REQUIRED)

#Modified Find-module to provide imported target LOGFT::LOGFT
find_package (LOGFT REQUIRED)

#Modified Find-module to provide imported target CPC::CPC
find_package (CPC REQUIRED)

target_link_libraries(RMat
  PRIVATE
#    Boost::boost
    GSL::gsl
    OpenMP::OpenMP_CXX
    AUSALIB::AUSALIB
    LOGFT::LOGFT
    CPC::CPC
    ROOT::TreePlayer
  PUBLIC
    ROOT::Libraries
    Armadillo::Armadillo   
)

#Install the library and headers
include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/RMat)

install(TARGETS RMat
    EXPORT rmat-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#Export the targets to a script
install(EXPORT rmat-targets
  FILE
    RMatTargets.cmake
  NAMESPACE
    RMat::
  DESTINATION
    ${INSTALL_CONFIGDIR}
)

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/RMatConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/RMatConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

#Install the config, and custom find modules
install(FILES
    ${CMAKE_CURRENT_LIST_DIR}/cmake/FindRMat.cmake
    ${CMAKE_CURRENT_LIST_DIR}/cmake/FindArmadillo.cmake
    ${CMAKE_CURRENT_LIST_DIR}/cmake/FindROOT.cmake
    ${CMAKE_CURRENT_LIST_DIR}/cmake/FindAUSALIB.cmake
    ${CMAKE_CURRENT_LIST_DIR}/cmake/FindLOGFT.cmake
    ${CMAKE_CURRENT_LIST_DIR}/cmake/FindCPC.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/RMatConfig.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

##############################################
## Exporting from the build tree
configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindRMat.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/FindRMat.cmake
COPYONLY)

export(EXPORT rmat-targets FILE ${CMAKE_CURRENT_BINARY_DIR}/RMatTargets.cmake NAMESPACE RMat::)
