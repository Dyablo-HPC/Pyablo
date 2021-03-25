#pragma once

#include <bits/stdc++.h>
#include "hdf5.h"

#include "Utils.h"

namespace dyablo {

struct Attribute {
  std::string name;
  std::string type;
  std::string center;
  hid_t handle;
};

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
 * @todo Time series
 * @todo Improve getRefinementCriterion by enabling the user to provide a lambda for 
 *       the calculation
 * @todo Finish getBlock when iOct writing has been corrected
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

  std::vector<int>   index_buffer;  //!< Connectivity info (HEAVY !)
  std::vector<float> vertex_buffer; //!< Coordinates info  (HEAVY !)

  int nDim;       //!< Number of dimensions of the dataset   
  int nElems;     //!< Number of indices per cell
  int nCells;     //!< Number of cells stored in the file
  int nVertices;  //!< Number of vertices stored in the file

  float time; //!< Current time of the snapshot

  static std::map<std::string, hid_t> type_corresp; //!< Mapping between type names and hid equivalents

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
   * @note: Please use thes for large query as most of them are made in parallel
   **/
  std::vector<int>   getCellsFromPositions(std::vector<Vec> v);
  std::vector<Vec>   getCellCenter(std::vector<uint> iCells);
  std::vector<Vec>   getCellSize(std::vector<uint> iCells);
  std::vector<float> getCellVolume(std::vector<uint> iCells);
  std::vector<Vec>   getUniqueCells(std::vector<Vec> pos);

  /** Domain info **/
  int getNCells();
  bool hasAttribute(std::string attribute);
  BoundingBox getDomainBoundingBox();
  
  /** Generic attribute probing methods **/
  template<typename T>
  T probeLocation(Vec pos, std::string attribute);
  template<typename T>
  std::vector<T> probeLocation(std::vector<Vec> pos, std::string attribute);
  template<typename T>
  std::vector<T> probeCells(std::vector<uint> iCells, std::string attribute);

  /** High-level probing methods **/
  float probeDensity(Vec pos);
  float probeTotalEnergy(Vec pos);
  Vec   probeMomentum(Vec pos);
  Vec   probeVelocity(Vec pos);
  int   probeLevel(Vec pos);
  int   probeRank(Vec pos);
  int   probeOctant(Vec pos);

  // Integrated quantities
  float getTotalMass();
  float getTotalEnergy();
  float getTotalKineticEnergy();

  float getRefinementCriterion(Vec pos);
  std::vector<float> getRefinementCriterion(std::vector<Vec> pos);  
  
  // Vector functions
  std::vector<float> probeDensity(std::vector<Vec> pos);
  std::vector<float> probeTotalEnergy(std::vector<Vec> pos);
  std::vector<Vec>   probeMomentum(std::vector<Vec> pos);
  std::vector<Vec>   probeVelocity(std::vector<Vec> pos);
  std::vector<int>   probeLevel(std::vector<Vec> pos);
  std::vector<int>   probeRank(std::vector<Vec> pos);
  std::vector<int>   probeOctant(std::vector<Vec> pos);

  // By cell
  std::vector<float> getDensity(std::vector<uint> iCells);
  std::vector<float> getTotalEnergy(std::vector<uint> iCells);
  std::vector<Vec>   getMomentum(std::vector<uint> iCells);
  std::vector<Vec>   getVelocity(std::vector<uint> iCells);
  std::vector<int>   getLevel(std::vector<uint> iCells);
  std::vector<int>   getRank(std::vector<uint> iCells);
  std::vector<int>   getOctant(std::vector<uint> iCells);

  std::vector<int>   getBlock(uint iOct);

  // Static variable for vectorized reading
  static int vec_size;
};

}