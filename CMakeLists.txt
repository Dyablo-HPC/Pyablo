cmake_minimum_required(VERSION 3.16)

project(dyablo-analysis LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(HDF5 REQUIRED)
message("${HDF5_INCLUDE_DIRS}")
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(includes)
include_directories(external/pugixml/src)

# Defining optional stuff
set(CXX_FLAGS 
  "-DPUGIXML_HEADER_ONLY" 
  "-DPUGIXML_NO_XPATH")

list(JOIN CXX_FLAGS " " CMAKE_CXX_FLAGS)

set(SRCS srcs/XdmfReader.cpp
         srcs/Snapshot.cpp)

add_library(dyablo-analysis ${SRCS})

# Testing stuff
add_executable(test main.cpp)
target_link_libraries(test dyablo-analysis ${HDF5_LIBRARIES})