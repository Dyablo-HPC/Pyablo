cmake_minimum_required(VERSION 3.16)

project(pyablo LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)

set(PYBIND11_FINDPYTHON ON CACHE BOOL "Use FindPython with pybind11")

find_package(OpenMP)
find_package(HDF5 REQUIRED)
set( Pyablo_HDF5_TARGETS hdf5::hdf5 )
if( ${HDF5_IS_PARALLEL} )
  find_package(MPI REQUIRED)
  list( APPEND Pyablo_HDF5_TARGETS MPI::MPI_CXX )
endif()
message(  ${Pyablo_HDF5_TARGETS} )

find_package(pybind11 REQUIRED)

include_directories(includes)
include_directories(external/pugixml/src)
include_directories(external/stb)

# Defining optional stuff
set(CXX_FLAGS 
  "-DPUGIXML_HEADER_ONLY" 
  "-DPUGIXML_NO_XPATH"
  "-fPIC")

list(JOIN CXX_FLAGS " " CMAKE_CXX_FLAGS)

set(SRC src/XdmfReader.cpp
         src/Snapshot.cpp
         src/Utils.cpp)

add_library(pyablo_lib ${SRC})
target_link_libraries(pyablo_lib PUBLIC OpenMP::OpenMP_CXX ${Pyablo_HDF5_TARGETS})

set(LIBS  pyablo_lib 
          ${Pyablo_HDF5_TARGETS}
          pybind11::module
          pybind11::lto
          pybind11::windows_extras)

##### Python module
pybind11_add_module(_pyablo src/pyablo.cpp)
target_link_libraries(_pyablo PUBLIC ${LIBS})

install(TARGETS _pyablo DESTINATION pyablo) # Here _pyablo is the name of the binary, while pyablo is the folder under which to place it
pybind11_strip(_pyablo)

set_target_properties(_pyablo PROPERTIES CXX_VISIBILITY_PRESET "hidden"
                              CUDA_VISIBILITY_PRESET "hidden")
