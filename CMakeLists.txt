# Required CMake Version
###############
cmake_minimum_required(VERSION 2.8)

# Project setup
###############
PROJECT(ViennaLS)

set(VERSION_MAJOR 2)
set(VERSION_MINOR 3)
set(VERSION_PATCH 2)
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

ADD_DEFINITIONS(-DVIENNATS_VERSION="${VERSION}")

# Default build type
##############
IF (NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release")
ENDIF()

OPTION(VIENNATS_STATIC_BUILD "" OFF)

# Default compile options
#################
SET(CMAKE_CXX_FLAGS "-Wall -pedantic")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 ")
SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

# Check flags for C++11 (C++0x)
##############
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()


# Specify general include directory
################
INCLUDE_DIRECTORIES(".")
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/viennats/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/viennats/LSlib/)

# Specify that the output executable is placed in the root build directory
################
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Find prerequisites
####################

# Look for additional CMake modules in the cmake subdirectory ..
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

###############################################
#  compiler dependent settings for OpenMP
###############################################
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    #MESSAGE(STATUS ${CMAKE_CXX_FLAGS})
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Werror -ferror-limit=2")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # disable-new-dtags is necessary for ld > 2.27
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Werror -fmax-errors=2")
    FIND_PACKAGE(OpenMP REQUIRED)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp  ")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fmax-errors=2")
endif()

# disable-new-dtags sets RPATH which searches for libs recursively, instead of RUNPATH which does not
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--disable-new-dtags")

# Boost ..
SET(BOOST_MIN_VERSION 1.46.0)
FIND_PACKAGE(Boost ${BOOST_MIN_VERSION} REQUIRED COMPONENTS filesystem system)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

# HDF5
IF(NOT VIENNATS_STATIC_BUILD)
FIND_PACKAGE(HDF5 COMPONENTS CXX)
IF(HDF5_FOUND)
 IF(HDF5_CXX_FOUND)
   add_definitions(-DUSE_HDF5)
   add_definitions(${HDF5_DEFINITIONS})
   INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
   SET(LIBRARIES ${LIBRARIES} ${HDF5_CXX_LIBRARIES})
 ELSE(HDF5_CXX_FOUND)
   message(STATUS "Missing HDF5 C++ library: ViennaTS will be built without TDR file support")
 ENDIF(HDF5_CXX_FOUND)
ELSE(HDF5_FOUND)
 message(STATUS "Missing HDF5 library: ViennaTS will be built without TDR file support")
ENDIF(HDF5_FOUND)
ENDIF(NOT VIENNATS_STATIC_BUILD)

# VTK
# VTK_DIR needs to point to the VTK installation and to VTK/lib/cmake/vtk-{vtk-version}
IF(VIENNATS_STATIC_BUILD)
  set(VTK_DIR $ENV{VTK_STATIC_DIR})
ELSE(VIENNATS_STATIC_BUILD)
  set(VTK_DIR $ENV{VTK_DIR})
ENDIF(VIENNATS_STATIC_BUILD)

find_package(VTK REQUIRED)

IF(VTK_FOUND)
  message(STATUS "VTK found (include: ${VTK_INCLUDE_DIRS}, lib: ${VTK_LIBRARY_DIRS})")
  include(${VTK_USE_FILE})
  #only link necessary VTK libraries
  set(VTK_LIBRARIES vtksys;vtkFiltersCore;vtkFiltersGeneral;vtkFiltersGeometry;vtkIOXML)
  message(STATUS "Including following VTK Libraries: ${VTK_LIBRARIES}")
  IF(VIENNATS_STATIC_BUILD)
    message(STATUS "BUILDING STATIC BINARY")
  ENDIF(VIENNATS_STATIC_BUILD)

ELSE(VTK_FOUND)
  message(STATUS "VTK not found.")
ENDIF(VTK_FOUND)


# Subdirectories
################
#add_subdirectory(src)

# Add visibility of headers
# Necessary for Qt-Creator usage.
################
file( GLOB_RECURSE MAKE_HEADERS_VISIBLE_SRC *.hpp *.h)
add_custom_target( MAKE_HEADERS_VISIBLE SOURCES ${MAKE_HEADERS_VISIBLE_SRC} )


# add executable and link
ADD_EXECUTABLE(viennats-${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH} src/viennats.cpp)

IF(VIENNATS_STATIC_BUILD)
 TARGET_LINK_LIBRARIES(viennats-${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH} ${LIBRARIES} ${Boost_LIBRARIES} ${VTK_LIBRARIES} "-static")
ELSE()
 TARGET_LINK_LIBRARIES(viennats-${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH} ${LIBRARIES} ${Boost_LIBRARIES} ${VTK_LIBRARIES})
ENDIF(VIENNATS_STATIC_BUILD)
