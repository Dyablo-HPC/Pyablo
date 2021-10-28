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

    .def("getCellsCenter", static_cast<std::vector<Vec> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getCellCenter))
    .def("getCellsSize",   static_cast<std::vector<Vec> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getCellSize))
    .def("getCellsVolume", static_cast<std::vector<float> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getCellVolume))

    .def("probeDensity",    static_cast<std::vector<float> (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeDensity))
    .def("probePressure",   static_cast<std::vector<float> (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probePressure))
    .def("probeEnergy",     static_cast<std::vector<float> (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeTotalEnergy))
    .def("probeMach",       static_cast<std::vector<float> (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeMach))
    .def("probeMomentum",   static_cast<std::vector<Vec>   (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeMomentum))
    .def("probeVelocity",   static_cast<std::vector<Vec>   (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeVelocity))
    .def("probeRank",       static_cast<std::vector<int>   (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeRank))
    .def("probeLevel",      static_cast<std::vector<int>   (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeLevel))
    .def("probeOctant",     static_cast<std::vector<int>   (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeOctant))
  
    .def("getDensity",  static_cast<std::vector<float> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getDensity))
    .def("getPressure", static_cast<std::vector<float> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getPressure))
    .def("getEnergy",   static_cast<std::vector<float> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getTotalEnergy))
    .def("getMach",     static_cast<std::vector<float> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getMach))
    .def("getMomentum", static_cast<std::vector<Vec>   (Snapshot::*)(std::vector<uint>)>(&Snapshot::getMomentum))
    .def("getVelocity", static_cast<std::vector<Vec>   (Snapshot::*)(std::vector<uint>)>(&Snapshot::getVelocity))
    .def("getLevel",    static_cast<std::vector<int>   (Snapshot::*)(std::vector<uint>)>(&Snapshot::getLevel))
    .def("getRank",     static_cast<std::vector<int>   (Snapshot::*)(std::vector<uint>)>(&Snapshot::getRank))
    .def("getOctant",   static_cast<std::vector<int>   (Snapshot::*)(std::vector<uint>)>(&Snapshot::getOctant))
   
    .def("getRefinementCriterion",  static_cast<std::vector<float> (Snapshot::*)(std::vector<Vec>)>(&Snapshot::getRefinementCriterion))
  
    .def("getTotalMass",                 static_cast<float (Snapshot::*)()>(&Snapshot::getTotalMass))
    .def("getTotalEnergy",               static_cast<float (Snapshot::*)()>(&Snapshot::getTotalEnergy))
    .def("getTotalKineticEnergy",        static_cast<float (Snapshot::*)()>(&Snapshot::getTotalKineticEnergy))
    .def("getTotalInternalEnergy",       static_cast<float (Snapshot::*)(double)>(&Snapshot::getTotalInternalEnergy))
    .def("getMaxMach",                   static_cast<float (Snapshot::*)()>(&Snapshot::getMaxMach))
    .def("getAverageMach",               static_cast<float (Snapshot::*)()>(&Snapshot::getAverageMach))
    .def("getTime",                      static_cast<float (Snapshot::*)()>(&Snapshot::getTime))

    .def("readAllFloat",     static_cast<std::vector<float> (Snapshot::*)(std::string)>(&Snapshot::readAllFloat))
    .def("mortonSort2d",     static_cast<std::vector<float> (Snapshot::*)(std::vector<float>, uint, uint, uint)>(&Snapshot::mortonSort2d))
    .def("mortonSort3d",     static_cast<std::vector<float> (Snapshot::*)(std::vector<float>, uint, uint, uint, uint)>(&Snapshot::mortonSort3d))
    .def("getSortingMask2d", static_cast<std::vector<uint64_t> (Snapshot::*)(uint, uint, uint)>(&Snapshot::getSortingMask2d))
    .def("getSortingMask3d", static_cast<std::vector<uint64_t> (Snapshot::*)(uint, uint, uint, uint)>(&Snapshot::getSortingMask3d))
   ;

  m.doc() = "dyablo-analysis python bindings"; // optional module docstring
}
