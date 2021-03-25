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
    .def("getNCells",             &Snapshot::getNCells)
    .def("getUniqueCells",        &Snapshot::getUniqueCells)
    .def("hasAttribute",          &Snapshot::hasAttribute)
    .def("getDomainBoundingBox",  &Snapshot::getDomainBoundingBox)

    .def("getCellCenter",   static_cast<Vec (Snapshot::*)(uint)>(&Snapshot::getCellCenter))
    .def("getCellsCenters", static_cast<std::vector<Vec> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getCellCenter))
    .def("getCellSize",     static_cast<Vec (Snapshot::*)(uint)>(&Snapshot::getCellSize))
    .def("getCellsSizes",   static_cast<std::vector<Vec> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getCellSize))
    .def("getCellVolume",   static_cast<float (Snapshot::*)(uint)>(&Snapshot::getCellVolume))
    .def("getCellsVolumes", static_cast<std::vector<float> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getCellVolume))

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
  
    .def("getDensity",  static_cast<std::vector<float> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getDensity))
    .def("getEnergy",   static_cast<std::vector<float> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getTotalEnergy))
    .def("getMomentum", static_cast<std::vector<Vec>   (Snapshot::*)(std::vector<uint>)>(&Snapshot::getMomentum))
    .def("getVelocity", static_cast<std::vector<Vec>   (Snapshot::*)(std::vector<uint>)>(&Snapshot::getVelocity))
    .def("getLevel",    static_cast<std::vector<int>   (Snapshot::*)(std::vector<uint>)>(&Snapshot::getLevel))
    .def("getRank",     static_cast<std::vector<int>   (Snapshot::*)(std::vector<uint>)>(&Snapshot::getRank))
    .def("getOctant",   static_cast<std::vector<int>   (Snapshot::*)(std::vector<uint>)>(&Snapshot::getOctant))
   
    .def("getRefinementCriterion", static_cast<float (Snapshot::*)(Vec)>(&Snapshot::getRefinementCriterion))
    .def("getRefinementCriteria",  static_cast<std::vector<float> (Snapshot::*)(std::vector<Vec>)>(&Snapshot::getRefinementCriterion))
  
    .def("getTotalMass",   static_cast<float (Snapshot::*)()>(&Snapshot::getTotalMass))
    .def("getTotalEnergy", static_cast<float (Snapshot::*)()>(&Snapshot::getTotalEnergy))
    .def("getTime",        static_cast<float (Snapshot::*)()>(&Snapshot::getTime))
   ;

  m.doc() = "dyablo-analysis python bindings"; // optional module docstring
}
