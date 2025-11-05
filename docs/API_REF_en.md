# Pyablo Module Documentation

This document provides a detailed reference for the `pyablo` module, which exposes C++ classes and functions to Python using Pybind11.

---

## Table of Contents
- [XdmfReader Class](#xdmfreader-class)
- [Snapshot Class](#snapshot-class)
- [Line Class](#line-class)

---

## XdmfReader Class

The `XdmfReader` class is designed to read XDMF (eXtensible Data Model and Format) files, which are commonly used for scientific data storage.

### Methods

| Method Signature | Description |
|-----------------|-------------|
| `XdmfReader()` | Constructor for the `XdmfReader` class. |
| `void readSnapshot()` | Reads a snapshot from an XDMF file. |
| `void readTimeSeries()` | Reads a time series from an XDMF file. |

---

## Snapshot Class

The `Snapshot` class represents a snapshot of simulation data, providing methods to access and manipulate cell data, physical quantities, and other simulation attributes.

### Methods

| Method Signature | Description |
|-----------------|-------------|
| `Snapshot()` | Constructor for the `Snapshot` class. |
| `void close()` | Closes the snapshot. |
| `void setName()` | Sets the name of the snapshot. |
| `void setNDim()` | Sets the number of dimensions. |
| `void print()` | Prints information about the snapshot. |
| `void getCellFromPosition()` | Gets a cell from a specified position. |
| `void getCellsFromPositions()` | Gets cells from specified positions. |
| `void getCellBoundingBox()` | Gets the bounding box of a cell. |
| `uint getNCells()` | Gets the total number of cells. |
| `void getUniqueCells()` | Gets unique cells. |
| `bool hasAttribute()` | Checks if an attribute exists. |
| `void getDomainBoundingBox()` | Gets the bounding box of the domain. |
| `std::vector<Vec> getCellsCenter(std::vector<uint>)` | Gets the centers of specified cells. |
| `std::vector<Vec> getCellsSize(std::vector<uint>)` | Gets the sizes of specified cells. |
| `std::vector<double> getCellsVolume(std::vector<uint>)` | Gets the volumes of specified cells. |
| `std::vector<double> probeQuantity(std::vector<Vec>, std::string)` | Probes a physical quantity from positions. |
| `std::vector<double> probeDensity(std::vector<Vec>)` | Probes density from positions. |
| `std::vector<double> probePressure(std::vector<Vec>)` | Probes pressure from positions. |
| `std::vector<double> probeEnergy(std::vector<Vec>)` | Probes total energy from positions. |
| `std::vector<double> probeMach(std::vector<Vec>)` | Probes Mach number from positions. |
| `std::vector<Vec> probeMomentum(std::vector<Vec>)` | Probes momentum from positions. |
| `std::vector<Vec> probeVelocity(std::vector<Vec>)` | Probes velocity from positions. |
| `std::vector<int> probeRank(std::vector<Vec>)` | Probes rank from positions. |
| `std::vector<int> probeLevel(std::vector<Vec>)` | Probes level from positions. |
| `std::vector<int> probeOctant(std::vector<Vec>)` | Probes octant from positions. |
| `std::vector<double> getQuantity(std::vector<uint>, std::string)` | Gets a physical quantity for specified cells. |
| `std::vector<double> getDensity(std::vector<uint>)` | Gets density for specified cells. |
| `std::vector<double> getPressure(std::vector<uint>)` | Gets pressure for specified cells. |
| `std::vector<double> getEnergy(std::vector<uint>)` | Gets total energy for specified cells. |
| `std::vector<double> getMach(std::vector<uint>)` | Gets Mach number for specified cells. |
| `std::vector<Vec> getMomentum(std::vector<uint>)` | Gets momentum for specified cells. |
| `std::vector<Vec> getVelocity(std::vector<uint>)` | Gets velocity for specified cells. |
| `std::vector<int> getLevel(std::vector<uint>)` | Gets level for specified cells. |
| `std::vector<int> getRank(std::vector<uint>)` | Gets rank for specified cells. |
| `std::vector<int> getOctant(std::vector<uint>)` | Gets octant for specified cells. |
| `std::vector<double> getRefinementCriterion(std::vector<Vec>)` | Gets refinement criterion for positions. |
| `double getTotalMass()` | Gets the total mass. |
| `double getTotalEnergy()` | Gets the total energy. |
| `double getTotalKineticEnergy()` | Gets the total kinetic energy. |
| `double getTotalInternalEnergy(double)` | Gets the total internal energy. |
| `double getMaxMach()` | Gets the maximum Mach number. |
| `double getAverageMach()` | Gets the average Mach number. |
| `double getTime()` | Gets the time. |
| `std::vector<double> readAllFloat(std::string)` | Reads all floats from a file. |
| `std::vector<double> mortonSort2d(std::vector<double>, uint, uint, uint)` | Sorts using Morton curve in 2D. |
| `std::vector<double> mortonSort3d(std::vector<double>, uint, uint, uint, uint)` | Sorts using Morton curve in 3D. |
| `std::vector<uint64_t> getSortingMask2d(uint, uint, uint, uint, uint)` | Gets a sorting mask in 2D. |
| `std::vector<uint64_t> getSortingMask2d(uint, uint, uint)` | Gets a sorting mask in 2D. |
| `std::vector<uint64_t> getSortingMask3d(uint, uint, uint, uint)` | Gets a sorting mask in 3D. |
| `std::vector<uint64_t> getSortingMask3d(uint, uint, uint, uint, uint, uint, uint)` | Gets a sorting mask in 3D. |
| `void fillLine()` | Fills a line with data. |
| `void fillLineUnique()` | Fills a line with unique data. |

---

## Line Class

The `Line` class represents a line of data, typically used for plotting or analysis.

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `Nl` | `uint` | Number of points in the line. |
| `start` | `Vec` | Start position. |
| `end` | `Vec` | End position. |
| `pos` | `Vec` | Positions of points. |
| `rho` | `double` | Density. |
| `prs` | `double` | Pressure. |
| `E` | `double` | Energy. |
| `vel` | `Vec` | Velocity. |
| `cellIds` | `std::vector<uint>` | Cell identifiers. |

---

