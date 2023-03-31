cmake_minimum_required(VERSION 3.16)

project(dyablo-analysis LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenMP)
find_package(HDF5 REQUIRED)
find_package(VTK REQUIRED NO_MODULE)
if (NOT VTK_FOUND)
  message("Error : ${VTK_NOT_FOUND_MESSAGE}")
  return()
endif()
#include(${VTK_USE_FILE})
message (STATUS "VTK_VERSION: ${VTK_VERSION}")


include_directories(${HDF5_INCLUDE_DIRS})
include_directories(includes)
include_directories(external/pugixml/src)
include_directories(external/pybind11/include)
include_directories(external/stb)

# Defining optional stuff
set(CXX_FLAGS 
  "-DPUGIXML_HEADER_ONLY" 
  "-DPUGIXML_NO_XPATH"
  "-fPIC")

list(JOIN CXX_FLAGS " " CMAKE_CXX_FLAGS)

set(SRCS srcs/XdmfReader.cpp
         srcs/Snapshot.cpp
         srcs/Utils.cpp)

add_library(dyablo-analysis ${SRCS})
target_link_libraries(dyablo-analysis OpenMP::OpenMP_CXX ${VTK_LIBRARIES})

set(LIBS  dyablo-analysis 
          ${HDF5_LIBRARIES}
          pybind11::module
          pybind11::lto
          pybind11::windows_extras
          OpenMP::OpenMP_CXX)

# Testing stuff
add_executable(dyablo-extract dyablo-extract.cpp)
target_link_libraries(dyablo-extract 
      PRIVATE 
      dyablo-analysis 
      ${HDF5_LIBRARIES}
      OpenMP::OpenMP_CXX)

include_directories(external/HighFive/include)
add_executable(dyablo-project dyablo-project.cpp)
target_link_libraries(dyablo-project 
      PRIVATE 
      dyablo-analysis 
      ${HDF5_LIBRARIES}
      OpenMP::OpenMP_CXX)

##### Python module
add_subdirectory(external/pybind11)
add_library(pydy MODULE srcs/pydy.cpp)
target_link_libraries(pydy PRIVATE ${LIBS})

pybind11_extension(pydy)
pybind11_strip(pydy)

set_target_properties(pydy PROPERTIES CXX_VISIBILITY_PRESET "hidden"
                           CUDA_VISIBILITY_PRESET "hidden")

