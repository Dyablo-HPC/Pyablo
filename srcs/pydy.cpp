#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "XdmfReader.h"
#include "Snapshot.h"

namespace py = pybind11;
using namespace dyablo;

/**
 * Binding stuff to the pydy module
 **/
PYBIND11_MODULE(pydy, m) {
  /**
   *  XdmfReader
   **/
  py::class_<XdmfReader>(m, "XdmfReader")
    .def(py::init())
    .def("readSnapshot", &XdmfReader::readSnapshot);
  
  /**
   * Snapshot
   **/
  py::class_<dyablo::Snapshot>(m, "Snapshot")
    .def(py::init())
    .def("close",   &Snapshot::close)
    .def("setName", &Snapshot::setName)
    .def("setNDim", &Snapshot::setNDim)
    .def("print",   &Snapshot::print)

    .def("getCellFromPosition",   &Snapshot::getCellFromPosition)
    .def("getCellsFromPositions", &Snapshot::getCellsFromPositions)
    .def("getCellBoundingBox",    &Snapshot::getCellBoundingBox)
    .def("getCellCenter",         &Snapshot::getCellCenter)
    .def("getNCells",             &Snapshot::getNCells)
    .def("getUniqueCells",        &Snapshot::getUniqueCells)

    .def("probeDensity",       static_cast<float (Snapshot::*)(Vec)>(&Snapshot::probeDensity))
    .def("probeTotalEnergy",   static_cast<float (Snapshot::*)(Vec)>(&Snapshot::probeTotalEnergy))
    .def("probeMomentum",      static_cast<Vec   (Snapshot::*)(Vec)>(&Snapshot::probeMomentum))
    .def("probeVelocity",      static_cast<Vec   (Snapshot::*)(Vec)>(&Snapshot::probeVelocity))
    .def("probeRank",          static_cast<int   (Snapshot::*)(Vec)>(&Snapshot::probeRank))
    .def("probeLevel",         static_cast<int   (Snapshot::*)(Vec)>(&Snapshot::probeLevel))
    .def("probeOctant",        static_cast<int   (Snapshot::*)(Vec)>(&Snapshot::probeOctant))

    .def("probeDensities",     static_cast<std::vector<float> (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeDensity))
    .def("probeTotalEnergies", static_cast<std::vector<float> (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeTotalEnergy))
    .def("probeMomenta",       static_cast<std::vector<Vec>   (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeMomentum))
    .def("probeVelocities",    static_cast<std::vector<Vec>   (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeVelocity))
    .def("probeRanks",         static_cast<std::vector<int>   (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeRank))
    .def("probeLevels",        static_cast<std::vector<int>   (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeLevel))
    .def("probeOctants",       static_cast<std::vector<int>   (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeOctant))
   ;

  m.doc() = "dyablo-analysis python bindings"; // optional module docstring
}
