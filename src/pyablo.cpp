#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "XdmfReader.h"
#include "Snapshot.h"

namespace py = pybind11;
using namespace dyablo;

/**
 * Binding stuff to the pydy module
 **/
PYBIND11_MODULE(_pyablo, m) {
  /**
   *  XdmfReader
   **/
  py::class_<dyablo::XdmfReader>(m, "XdmfReader")
    .def(py::init())
    .def("readSnapshot", &XdmfReader::readSnapshot)
    .def("readTimeSeries", &XdmfReader::readTimeSeries);
  
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
    .def("getCellsVolume", static_cast<std::vector<double> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getCellVolume))

    .def("probeQuantity",   static_cast<std::vector<double> (Snapshot::*)(std::vector<Vec>, std::string)>(&Snapshot::probeQuantity))
    .def("probeDensity",    static_cast<std::vector<double> (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeDensity))
    .def("probePressure",   static_cast<std::vector<double> (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probePressure))
    .def("probeEnergy",     static_cast<std::vector<double> (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeTotalEnergy))
    .def("probeMach",       static_cast<std::vector<double> (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeMach))
    .def("probeMomentum",   static_cast<std::vector<Vec>   (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeMomentum))
    .def("probeVelocity",   static_cast<std::vector<Vec>   (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeVelocity))
    .def("probeRank",       static_cast<std::vector<int>   (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeRank))
    .def("probeLevel",      static_cast<std::vector<int>   (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeLevel))
    .def("probeOctant",     static_cast<std::vector<int>   (Snapshot::*)(std::vector<Vec>)>(&Snapshot::probeOctant))
  
    .def("getQuantity", static_cast<std::vector<double> (Snapshot::*)(std::vector<uint>, std::string)>(&Snapshot::getQuantity))
    .def("getDensity",  static_cast<std::vector<double> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getDensity))
    .def("getPressure", static_cast<std::vector<double> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getPressure))
    .def("getEnergy",   static_cast<std::vector<double> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getTotalEnergy))
    .def("getMach",     static_cast<std::vector<double> (Snapshot::*)(std::vector<uint>)>(&Snapshot::getMach))
    .def("getMomentum", static_cast<std::vector<Vec>   (Snapshot::*)(std::vector<uint>)>(&Snapshot::getMomentum))
    .def("getVelocity", static_cast<std::vector<Vec>   (Snapshot::*)(std::vector<uint>)>(&Snapshot::getVelocity))
    .def("getLevel",    static_cast<std::vector<int>   (Snapshot::*)(std::vector<uint>)>(&Snapshot::getLevel))
    .def("getRank",     static_cast<std::vector<int>   (Snapshot::*)(std::vector<uint>)>(&Snapshot::getRank))
    .def("getOctant",   static_cast<std::vector<int>   (Snapshot::*)(std::vector<uint>)>(&Snapshot::getOctant))
   
    .def("getRefinementCriterion",  static_cast<std::vector<double> (Snapshot::*)(std::vector<Vec>)>(&Snapshot::getRefinementCriterion))
  
    .def("getTotalMass",                 static_cast<double (Snapshot::*)()>(&Snapshot::getTotalMass))
    .def("getTotalEnergy",               static_cast<double (Snapshot::*)()>(&Snapshot::getTotalEnergy))
    .def("getTotalKineticEnergy",        static_cast<double (Snapshot::*)()>(&Snapshot::getTotalKineticEnergy))
    .def("getTotalInternalEnergy",       static_cast<double (Snapshot::*)(double)>(&Snapshot::getTotalInternalEnergy))
    .def("getMaxMach",                   static_cast<double (Snapshot::*)()>(&Snapshot::getMaxMach))
    .def("getAverageMach",               static_cast<double (Snapshot::*)()>(&Snapshot::getAverageMach))
    .def("getTime",                      static_cast<double (Snapshot::*)()>(&Snapshot::getTime))

    .def("readAllFloat",     static_cast<std::vector<double> (Snapshot::*)(std::string)>(&Snapshot::readAllFloat))
    .def("mortonSort2d",     static_cast<std::vector<double> (Snapshot::*)(std::vector<double>, uint, uint, uint)>(&Snapshot::mortonSort2d))
    .def("mortonSort3d",     static_cast<std::vector<double> (Snapshot::*)(std::vector<double>, uint, uint, uint, uint)>(&Snapshot::mortonSort3d))
    .def("getSortingMask2d", static_cast<std::vector<uint64_t> (Snapshot::*)(uint, uint, uint, uint, uint)>(&Snapshot::getSortingMask2d))
    .def("getSortingMask2d", static_cast<std::vector<uint64_t> (Snapshot::*)(uint, uint, uint)>(&Snapshot::getSortingMask2d))
    .def("getSortingMask3d", static_cast<std::vector<uint64_t> (Snapshot::*)(uint, uint, uint, uint)>(&Snapshot::getSortingMask3d))
    .def("getSortingMask3d", static_cast<std::vector<uint64_t> (Snapshot::*)(uint, uint, uint, uint, uint, uint, uint)>(&Snapshot::getSortingMask3d))

    .def("fillLine",       &Snapshot::fillLine)
    .def("fillLineUnique", &Snapshot::fillLineUnique)
   ;

   /**
    * Line structure
    **/
   py::class_<dyablo::Line>(m, "Line")
     .def(py::init())
     .def_readwrite("Nl",     &Line::Nl)
     .def_readwrite("start",  &Line::start)
     .def_readwrite("end",    &Line::end)
     .def_readonly("pos",     &Line::pos)
     .def_readonly("rho",     &Line::rho)
     .def_readonly("prs",     &Line::prs)
     .def_readonly("E",       &Line::E)
     .def_readonly("vel",     &Line::vel)
     .def_readonly("cellIds", &Line::cellIds);

  m.doc() = "pyablo python bindings"; // optional module docstring
}
