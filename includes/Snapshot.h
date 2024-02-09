#pragma once

#include <bits/stdc++.h>

#include "Utils.h"

namespace dyablo {

constexpr size_t CoordSize = 3; // Coordinates are stored as 3-float even in 2D

/**
 * Class storing info about dyablo snapshots
 * 
 * A general note : The strategy here is to never store arrays directly in
 * ram except for coordinates and connectivity for rapid lookup. 
 * Future optimization might include streaming or partial buffering
 * but we do this for the moment to avoid any memory explosion for big runs
 * 
 * @todo buffer storing for probing ? Maybe do a region extraction as in yt ?
 * @todo I'm pretty sure probeLocation can be written in a nice templated without
 *       having to explicitely define the template when calling it ...
 * @todo Improve getRefinementCriterion by enabling the user to provide a lambda for 
 *       the calculation
 * @todo Finish getBlock when iOct writing has been corrected
 * @todo switch static allocation of the selecters for probe/get to a heap allocation
 *       to avoid segfaults when having large datasets 
 * @todo How to detach EOS stuff from hardcoding ? -> Maybe the whole thing
 *       should be directly based on the dyablo code to reuse modules as possible
 **/
class Snapshot {
 private:
  /**
   * Hdf5 related stuff
   **/
  std::string name;                     //!< Name of the set
  std::map<std::string, hid_t> handles; //!< Map of all the opened file handles
  std::vector<hid_t> data_handles;      //!< List of all the opened dataset handles

  hid_t connectivity;                          //!< Handle to connectivity dataset
  hid_t coordinates;                           //!< Handle to coordinates dataset
  std::map<std::string, Attribute> attributes; //!< Map of attributes

  IntArray  index_buffer;  //!< Connectivity info (HEAVY !)
  RealArray vertex_buffer; //!< Coordinates info  (HEAVY !)

  int nDim;       //!< Number of dimensions of the dataset   
  int nElems;     //!< Number of indices per cell
  int nCells;     //!< Number of cells stored in the file
  int nVertices;  //!< Number of vertices stored in the file

  float time; //!< Current time of the snapshot

  static std::map<std::string, hid_t> type_corresp; //!< Mapping between type names and hid equivalents

  bool show_progress_bars;
 public:
  Snapshot() = default;
  ~Snapshot();

  void close();

  void print();

  /** Snapshot reading/construction from Hdf5 **/
  void setName(std::string name);
  void setTime(float time);
  void setNDim(int nDim);
  void addH5Handle(std::string handle, std::string filename);
  void setConnectivity(std::string handle, std::string xpath, int nCells);
  void setCoordinates(std::string handle, std::string xpath, int nVertices);
  void addAttribute(std::string handle, std::string xpath, std::string name, std::string type, std::string center);

  /** Cell info and access **/
  int getCellFromPosition(Vec v);
  BoundingBox getCellBoundingBox(uint iCell);
  Vec getCellCenter(uint iCell);
  Vec getCellSize(uint iCell);
  float getCellVolume(uint iCell);
  float getTime();

  /** Vector access 
   * @note: Please use these for large query as most of them are made in parallel
   **/
  UIntArray  getCellsFromPositions(VecArray v);
  VecArray   getCellCenter(UIntArray iCells);
  VecArray   getCellSize(UIntArray iCells);
  RealArray  getCellVolume(UIntArray iCells);
  VecArray   getUniqueCells(VecArray pos);

  /** Domain info **/
  int getNCells();
  bool hasAttribute(std::string attribute);
  BoundingBox getDomainBoundingBox();
  
  /** Generic attribute probing methods **/
  template<typename T>
  T probeLocation(Vec pos, std::string attribute);
  template<typename T>
  std::vector<T> probeLocation(VecArray pos, std::string attribute);
  template<typename T>
  std::vector<T> probeCells(UIntArray iCells, std::string attribute);

  /** High-level probing methods **/
  float probeQuantity(Vec pos, std::string attribute);
  float probeDensity(Vec pos);
  float probePressure(Vec pos);
  float probeTotalEnergy(Vec pos);
  float probeMach(Vec pos);
  Vec   probeMomentum(Vec pos);
  Vec   probeVelocity(Vec pos);
  int   probeLevel(Vec pos);
  int   probeRank(Vec pos);
  int   probeOctant(Vec pos);

  // Integrated quantities
  float getTotalMass();
  float getTotalEnergy();
  float getTotalInternalEnergy(double gamma);
  float getTotalKineticEnergy();
  float getMaxMach();
  float getAverageMach();

  float     getRefinementCriterion(Vec pos);
  RealArray getRefinementCriterion(VecArray pos);  
  
  // Vector functions
  RealArray probeQuantity(VecArray pos, std::string attribute);
  RealArray probeDensity(VecArray pos);
  RealArray probePressure(VecArray pos);
  RealArray probeTotalEnergy(VecArray pos);
  RealArray probeMach(VecArray pos);
  VecArray  probeMomentum(VecArray pos);
  VecArray  probeVelocity(VecArray pos);
  IntArray  probeLevel(VecArray pos);
  IntArray  probeRank(VecArray pos);
  IntArray  probeOctant(VecArray pos);

  // By cell
  RealArray getDensity(UIntArray iCells);
  RealArray getPressure(UIntArray iCells);
  RealArray getTotalEnergy(UIntArray iCells);
  RealArray getMach(UIntArray iCells);
  VecArray  getMomentum(UIntArray iCells);
  VecArray  getVelocity(UIntArray iCells);
  IntArray  getLevel(UIntArray iCells);
  IntArray  getRank(UIntArray iCells);
  IntArray  getOctant(UIntArray iCells);

  IntArray  getBlock(uint iOct);

  RealArray   readAllFloat(std::string attribute);
  RealArray   mortonSort2d(RealArray vec, uint iLevel, uint bx, uint by);
  RealArray   mortonSort3d(RealArray vec, uint iLevel, uint bx, uint by, uint bz);
  UInt64Array getSortingMask2d(uint iLevel, uint bx, uint by);
  UInt64Array getSortingMask3d(uint iLevel, uint bx, uint by, uint bz);
  UInt64Array getSortingMask2d(uint iLevel, uint bx, uint by, uint coarse_res_x, uint coarse_res_y);
  UInt64Array getSortingMask3d(uint iLevel, uint bx, uint by, uint bz, uint coarse_res_x, uint coarse_res_y, uint coarse_res_z);

  void fillLine(Line &line);
  void fillSlice(Slice &slice);

  // Static variable for vectorized reading
  static int vec_size;
};

}