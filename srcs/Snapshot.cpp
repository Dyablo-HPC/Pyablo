#include "Snapshot.h"
#include <omp.h>

namespace dyablo {

int Snapshot::vec_size = 1000000;

/**
 * Static mapping between type names in xdmf file and H5T native types
 **/
std::map<std::string, hid_t> Snapshot::type_corresp = {
  {"Int", H5T_NATIVE_INT},
  {"Float", H5T_NATIVE_FLOAT}
};

/**
 * Destructor
 * Closes all HDF5 handles opened during the lifetime of the Snapshot
 **/
Snapshot::~Snapshot() {
}

/**
 * Closes all the handles in the file
 * This is not done in the destructor to ensure compatibility
 * with pybind11 !
 **/
void Snapshot::close() {
  for (auto dh: data_handles)
    H5Dclose(dh);

  for (auto [k, h]: handles)
    H5Fclose(h);
}

/**
 * Pretty prints a brief summary of the snapshot
 **/
void Snapshot::print() {
  std::cout << "Snapshot : " << name << std::endl;
  std::cout << " . Number of dimensions : " << nDim << std::endl;
  std::cout << " . Grid has " << nVertices << " vertices and " << nCells << " cells" << std::endl;
  std::cout << " . Attribute list : " << std::endl;
  for (auto [name, att]: attributes)
    std::cout << "   o " << name << " (" << att.type << ")" << std::endl;
  std::cout << " . Analysing with " << omp_get_max_threads() << " threads" << std::endl;
}

/**
 * Sets the associated name of the snapshot
 * @param name the name to give to the snapshot
 **/
void Snapshot::setName(std::string name) {
  this->name = name;
}

/**
 * Sets the associated time of the snapshot
 * @param time the time of the current snapshot
 **/
void Snapshot::setTime(float time) {
  this->time = time;
}

/**
 * Defines the number of dimensions of the snapshot.
 * This also sets the size of an element connectivity in memory
 * @param ndim the number of dimensions. Should be 2 or 3.
 **/
void Snapshot::setNDim(int nDim) {
  this->nDim = nDim;
  nElems = (nDim == 2 ? 4 : 8);
}

/**
 * Adds a handle to an HDF5 file if not already opened
 * @param handle the string reference to the actual file
 * @param filename the complete path to the hdf5 file
 **/
void Snapshot::addH5Handle(std::string handle, std::string filename) {
  if (handles.count(handle) == 0) {
    hid_t hid = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (hid < 0) {
      std::cerr << "ERROR : Could not open file " << filename.c_str() << std::endl;
      std::exit(1);
    }
    handles[handle] = hid;
  }
}

/**
 * Defines the connectivity of the snapshot using an hdf5 reference.
 * The connectivity is stored in index_buffer
 * @param handle the handle from which the dataset will be read
 * @param xpath the path to the dataset in hdf5 space
 * @param nCells the number of cells in the dataset
 **/
void Snapshot::setConnectivity(std::string handle, std::string xpath, int nCells) {
  connectivity = H5Dopen2(handles[handle], xpath.c_str(), H5P_DEFAULT);
  if (connectivity < 0) {
    std::cerr << "ERROR : Could not access connectivity info at " << handle << ":" << xpath << std::endl;
    std::exit(1);
  }
  this->nCells = nCells;
  data_handles.push_back(connectivity);

  // We read all the indices in memory
  uint nElem = (nDim == 2 ? 4 : 8);
  index_buffer.resize(nCells*nElem);
  herr_t status = H5Dread(connectivity, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, index_buffer.data());
  if (status < 0) {
    std::cerr << "ERROR while reading coordinates !" << std::endl;
    std::exit(1);
  }
}

/**
 * Defines the coordinates of all vertices in the mesh using an hdf5 reference
 * The coordinates are stored in vertex_buffer
 * @param handle the handle from which the dataset will be read
 * @param xpath the path to the dataset in hdf5 space
 * @param nVertices the number of vertices in the dataset
 **/
void Snapshot::setCoordinates(std::string handle, std::string xpath, int nVertices) {
  coordinates = H5Dopen2(handles[handle], xpath.c_str(), H5P_DEFAULT);
  if (coordinates < 0) {
    std::cerr << "ERROR : Could not access coordinates info at " << handle << "/" << xpath << std::endl;
    std::exit(1);
  }
  this->nVertices = nVertices;
  data_handles.push_back(coordinates);

  // We read all the coords in memory
  vertex_buffer.resize(nVertices*CoordSize);
  herr_t status = H5Dread(coordinates, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vertex_buffer.data());
  if (status < 0) {
    std::cerr << "ERROR while reading coordinates !" << std::endl;
    std::exit(1);
  }
}

/**
 * Adds a reference to an attribute in an HDF5 file
 * @param handle the handle from which the dataset will be read
 * @param xpath the path to the dataset in hdf5 space
 * @param name the name of the attribute and how it will be referenced later on
 * @param type the type of data stored in the file
 * @param center where is located the data in the cell
 **/
void Snapshot::addAttribute(std::string handle, std::string xpath, std::string name, std::string type, std::string center) {
  hid_t att_handle = H5Dopen2(handles[handle], xpath.c_str(), H5P_DEFAULT);
  if (att_handle < 0) {
    std::cout << "ERROR : Could not access attribute info at " << handle << "/" << xpath << std::endl;
    std::exit(1);
  }
  
  Attribute att {name, type, center, att_handle};
  attributes[name] = att;
  data_handles.push_back(att_handle);
}

/**
 * Finds which element corresponds to a specific location
 * @param pos the position of the element to probe
 * @return an integer corresponding to the cell id holding position pos
 **/
int Snapshot::getCellFromPosition(Vec pos) {
  for (uint iCell=0; iCell < nCells; ++iCell) {
    auto [min, max] = getCellBoundingBox(iCell);

    if (min[0] <= pos[0] && max[0] >= pos[0] 
      && min[1] <= pos[1] && max[1] >= pos[1]
      && (nDim==2 || (min[2] <= pos[2] && max[2] >= pos[2])))
      return iCell;
  }

  return -1;
}

/**
 * Finds which elements correspond to specific locations
 * This function can be especially long on large data sets !
 * @param pos a vector of positions to identify
 * @return a vector of integers corresponding to the cell ids holding the positions
 **/
std::vector<int> Snapshot::getCellsFromPositions(std::vector<Vec> pos) {
  std::vector<int> res(pos.size());

  int nPos = pos.size();
  #pragma omp parallel for shared(res), schedule(dynamic)
  for (uint iCell=0; iCell < nCells; ++iCell) {
    auto [min, max] = getCellBoundingBox(iCell);

    for (int iPos=0; iPos < nPos; ++iPos) {
      Vec &p = pos[iPos];
      if (min[0] <= p[0] && max[0] >= p[0] 
        && min[1] <= p[1] && max[1] >= p[1]
        && (nDim==2 || (min[2] <= p[2] && max[2] >= p[2])))
      res[iPos] = iCell;
    }
  }  
  return res;
}

/**
 * Returns the bounding box corresponding to a cell
 * @param iCell the id of the cell to probe
 * @return a pair of Vec storing the minimum 
 *         and maximum coordinates of the bounding box
 **/
BoundingBox Snapshot::getCellBoundingBox(uint iCell) {
  Vec min, max;

  int c0 = index_buffer[iCell*nElems];
  for (int i=0; i < nDim; ++i) {
    min[i] = vertex_buffer[c0*CoordSize + i];
    max[i] = min[i]; 
  }

  for (uint i=1; i < nElems; ++i) {
    int ci = index_buffer[iCell*nElems+i];
    float* coords = &vertex_buffer[ci*CoordSize];
    for (int i=0; i < nDim; ++i) { 
      min[i] = std::min(min[i], coords[i]);
      max[i] = std::max(max[i], coords[i]);
    }
  }

  return std::make_pair(min, max);
}

/**
 * Returns the center of the given cell
 * @param iCell the index of the cell to probe
 * @return a vector giving the center position of the cell iCell
 **/
Vec Snapshot::getCellCenter(uint iCell) {
  Vec out{0.0};
  for (int i=0; i < nElems; ++i) {
    int ci = index_buffer[iCell*nElems+i];
    float *coords = &vertex_buffer[ci*CoordSize];
    out[0] += coords[0];
    out[1] += coords[1];
    if (nDim == 3)
      out[2] += coords[2];
  }

  out[0] /= nElems;
  out[1] /= nElems;
  if (nDim == 3)
    out[2] /= nElems;
  return out;
}

/**
 * Returns the centers of given cells
 * @param iCells the indices of the cells to probe
 * @return a vector of positions corresponding to the center of the cells
 **/
std::vector<Vec> Snapshot::getCellCenter(std::vector<uint> iCells) {
  uint nPos = iCells.size();
  std::vector<Vec> out(nPos);
  
  #pragma omp parallel for schedule(dynamic)
  for (int i=0; i < nPos; ++i) {
    uint iCell = iCells[i];
    for (int j=0; j < nElems; ++j) {
      int ci = index_buffer[iCell*nElems+j];
      float *coords = &vertex_buffer[ci*CoordSize];
      out[i][0] += coords[0];
      out[i][1] += coords[1];
      if (nDim == 3)
        out[i][2] += coords[2];
    }

    out[i][0] /= nElems;
    out[i][1] /= nElems;
    if (nDim == 3)
      out[i][2] /= nElems;
  }

  return out;
}

/**
 * Returns the size of a cell
 * @param iCell the index of the cell to probe
 * @return a vector indicating the size of the cell along each axis
 **/
Vec Snapshot::getCellSize(uint iCell) {
  BoundingBox bb = getCellBoundingBox(iCell);
  Vec out;
  for (int i=0; i < nDim; ++i)
    out[i] = bb.second[i] - bb.first[i];
  return out;
}

/**
 * Returns the size of a list of cells
 * @param iCells a vector of cells to probe
 * @return a vector of Vec indicating the size of each cell probed
 **/
std::vector<Vec> Snapshot::getCellSize(std::vector<uint> iCells) {
  uint nSizes = iCells.size();
  std::vector<Vec> out(nSizes);
  
  #pragma omp parallel for schedule(dynamic)
  for (int i=0; i < nSizes; ++i) {
    uint iCell = iCells[i];
    out[i] = getCellSize(iCell);
  }

  return out;
}

/**
 * Returns the volume/surface of a cell
 * @param iCell the index of the cell to probe
 * @return a float indicating the surface in 2D or the volume in 3D of the cell
 **/
float Snapshot::getCellVolume(uint iCell) {
  BoundingBox bb = getCellBoundingBox(iCell);
  float out = bb.second[0] - bb.first[0];
  for (int i=1; i < nDim; ++i)
    out *= bb.second[i] - bb.first[i];
  return out;
}

/**
 * Returns the volume/surface of a cell
 * @param iCells a vector of cells to probe
 * @return a vector of floats indicating the volume/surface of each cell probed
 **/
std::vector<float> Snapshot::getCellVolume(std::vector<uint> iCells) {
  uint nSizes = iCells.size();
  std::vector<float> out(nSizes);
  
  #pragma omp parallel for schedule(dynamic)
  for (int i=0; i < nSizes; ++i) {
    uint iCell = iCells[i];
    out[i] = getCellVolume(iCell);
  }

  return out;
}

/**
 * Gets the associated time of the snapshot
 * @return the time of the current snapshot
 **/
float Snapshot::getTime() {
  return time;
}

/**
 * Returns the number of cells in the domain
 * @return the number of cells in the domain
 **/
int Snapshot::getNCells() {
  return nCells;
}

/**
 * Returns whether or not the Snapshot contains an attribute
 * @param attribute the name of the attribute to test
 * @returns a boolean indicating if attribute is stored in the Snapshot
 **/
bool Snapshot::hasAttribute(std::string attribute) {
  return attributes.count(attribute) > 0;
}

/**
 * Returns the bounding box of the domain
 * @note This is a very naive implementation where we consider only a cartesian-box
 *       domain ! This will not work using other geometries
 **/
BoundingBox Snapshot::getDomainBoundingBox() {
  BoundingBox out{getCellCenter(0), getCellCenter(nCells-1)};
  Vec s0 = getCellSize(0);
  Vec s1 = getCellSize(nCells-1);
  for (int i=0; i < nDim; ++i) {
    out.first[i]  -= s0[i]*0.5;
    out.second[i] += s1[i]*0.5;
  }
  return out;
}

/**
 * Probes a location for an attribute
 * @param T the type of result
 * @param pos the position to probe
 * @param attribute the name of the attribute to probe
 * @return a value of type T corresponding to the attribute at location pos
 * 
 * @todo Implement non-centered probing
 * @todo Add smoothing/interpolation
 * @todo Check that attribute actually exists
 * @todo Check for h5 errors while reading
 **/
template<typename T>
T Snapshot::probeLocation(Vec pos, std::string attribute) {
  if (attributes.count(attribute) == 0) {
    std::cerr << "ERROR : Attribute " << attribute << " is not stored in file !" << std::endl;
    return 0;
  }

  // Retrieving cell index
  int iCell = getCellFromPosition(pos);

  Attribute &att = attributes[attribute];
  hid_t data_type = type_corresp[att.type];

  // Selecting the element in the dataset
  herr_t status;
  hid_t space = H5Dget_space(att.handle);
  hsize_t select[1][1] = {(hsize_t)iCell};
  status = H5Sselect_elements(space, H5S_SELECT_SET, 1, (const hsize_t *)&select);

  hsize_t d=1;
  hid_t memspace = H5Screate_simple(1, &d, NULL);

  T value;

  if (data_type == H5T_NATIVE_INT) {
    herr_t status = H5Dread(att.handle, data_type, memspace, space, H5P_DEFAULT, &value);
  }
  else if (data_type == H5T_NATIVE_FLOAT) {
    herr_t status = H5Dread(att.handle, data_type, memspace, space, H5P_DEFAULT, &value);
  }

  return value;
}

/**
 * Probes multiple cells for an attribute
 * @param T the type of result
 * @param iCells a vector ids of the cells to probe
 * @param attribute the name of the attribute to probe
 * @return a vector of type T corresponding to the attribute at location pos
 * 
 * @todo Implement non-centered probing
 * @todo Add smoothing/interpolation
 * @todo Check that attribute actually exists
 * @todo Check for h5 errors while reading
 **/
template<typename T>
std::vector<T> Snapshot::probeCells(std::vector<uint> iCells, std::string attribute) {
  std::vector<T> out;
  if (attributes.count(attribute) == 0) {
    std::cerr << "ERROR : Attribute " << attribute << " is not stored in file !" << std::endl;
    return out;
  }
  
    // Retrieving cell indices and adding them to selection array
  int nCells = iCells.size();

  hsize_t* select = new hsize_t[nCells];
  for (int i=0; i < nCells; ++i)
    select[i] = (hsize_t)iCells[i];


  Attribute &att = attributes[attribute];
  hid_t data_type = type_corresp[att.type];

  // Selecting the element in the dataset
  herr_t status;
  hid_t space = H5Dget_space(att.handle);
  status = H5Sselect_elements(space, H5S_SELECT_SET, nCells, select);

  hsize_t d=nCells;
  hid_t memspace = H5Screate_simple(1, &d, NULL);

  out.resize(nCells);

  if (data_type == H5T_NATIVE_INT) {
    herr_t status = H5Dread(att.handle, data_type, memspace, space, H5P_DEFAULT, out.data());
  }
  else if (data_type == H5T_NATIVE_FLOAT) {
    herr_t status = H5Dread(att.handle, data_type, memspace, space, H5P_DEFAULT, out.data());
  }

  delete [] select;
  H5Sclose(memspace);
  H5Sclose(space);

  return out;
}

/**
 * Probes multiple locations for an attribute
 * @param T the type of result
 * @param pos a vector of positions to probe
 * @param attribute the name of the attribute to probe
 * @return a vector of type T corresponding to the attribute at location pos
 * 
 * @todo Implement non-centered probing
 * @todo Add smoothing/interpolation
 * @todo Check that attribute actually exists
 * @todo Check for h5 errors while reading
 **/
template<typename T>
std::vector<T> Snapshot::probeLocation(std::vector<Vec> pos, std::string attribute) {
  std::vector<T> out;
  if (attributes.count(attribute) == 0) {
    std::cerr << "ERROR : Attribute " << attribute << " is not stored in file !" << std::endl;
    return out;
  }
  
    // Retrieving cell indices and adding them to selection array
  int nPos = pos.size();
  std::vector<int> iCells = getCellsFromPositions(pos);

  hsize_t* select = new hsize_t[nPos];
  for (int i=0; i < nPos; ++i)
    select[i] = (hsize_t)iCells[i];

  Attribute &att = attributes[attribute];
  hid_t data_type = type_corresp[att.type];

  // Selecting the element in the dataset
  herr_t status;
  hid_t space = H5Dget_space(att.handle);
  status = H5Sselect_elements(space, H5S_SELECT_SET, nPos, select);

  hsize_t d=nPos;
  hid_t memspace = H5Screate_simple(1, &d, NULL);

  out.resize(nPos);

  if (data_type == H5T_NATIVE_INT) {
    herr_t status = H5Dread(att.handle, data_type, memspace, space, H5P_DEFAULT, out.data());
  }
  else if (data_type == H5T_NATIVE_FLOAT) {
    herr_t status = H5Dread(att.handle, data_type, memspace, space, H5P_DEFAULT, out.data());
  }

  delete [] select;
  H5Sclose(memspace);
  H5Sclose(space);

  return out;
}

/** 
 * Probes a location for density
 * @param pos the position to probe
 * @return the density at position pos
 **/
float Snapshot::probeDensity(Vec pos) {
  return probeLocation<float>(pos, "rho");
}

/** 
 * Probes a location for total energy
 * @param pos the position to probe
 * @return the total energy at position pos
 **/
float Snapshot::probeTotalEnergy(Vec pos) {
  return probeLocation<float>(pos, "e_tot");
}

/** 
 * Probes a location for the pressure
 * @param pos the position to probe
 * @return the total energy at position pos
 **/
float Snapshot::probePressure(Vec pos) {
  return probeLocation<float>(pos, "P");
}

/** 
 * Probes a location for the Mach number
 * @param pos the position to probe
 * @return the Mach number of the flow at position pos
 **/
float Snapshot::probeMach(Vec pos) {
  return probeLocation<float>(pos, "Mach");
}

/** 
 * Probes a location for momentum
 * @param pos the position to probe
 * @return the momentum at position pos
 **/
Vec Snapshot::probeMomentum(Vec pos) {
  Vec res;
  res[0] = probeLocation<float>(pos, "rho_vx");
  res[1] = probeLocation<float>(pos, "rho_vy");
  if (nDim == 3)
    res[2] = probeLocation<float>(pos, "rho_vz");

  return res;
}

/** 
 * Probes a location for velocity
 * @param pos the position to probe
 * @return the velocity at position pos
 **/
Vec Snapshot::probeVelocity(Vec pos) {
  Vec res = probeMomentum(pos);
  float rho = probeLocation<float>(pos, "rho");
  for (int i=0; i < nDim; ++i)
    res[i] /= rho;
  return res;
}

/** 
 * Probes a location for the refinement level
 * @param pos the position to probe
 * @return the refinement level at position pos
 **/
int Snapshot::probeLevel(Vec pos) {
  return probeLocation<int>(pos, "level");
}

/** 
 * Probes a location for the mpi-rank
 * @param pos the position to probe
 * @return the mpi-rank at position pos
 **/
int Snapshot::probeRank(Vec pos) {
  return probeLocation<int>(pos, "rank");
}

/**
 * Probes a location for the octant index
 * @param pos a position to probe
 * @return the octant index at position pos
 **/
int Snapshot::probeOctant(Vec pos) {
  return probeLocation<int>(pos, "iOct");
}

/** 
 * Probes multiple locations for density
 * @param pos a vector of positions to probe
 * @return the density at positions pos
 **/
std::vector<float> Snapshot::probeDensity(std::vector<Vec> pos) {
  return probeLocation<float>(pos, "rho");
}

/** 
 * Probes multiple locations for pressure
 * @param pos a vector of positions to probe
 * @return the pressure at positions pos
 **/
std::vector<float> Snapshot::probePressure(std::vector<Vec> pos) {
  return probeLocation<float>(pos, "P");
}

/** 
 * Probes multiple locations for total energy
 * @param pos a vector of positions to probe
 * @return the total energy at positions pos
 **/
std::vector<float> Snapshot::probeTotalEnergy(std::vector<Vec> pos) {
  return probeLocation<float>(pos, "e_tot");
}

/**
 * Probes multiple locations for the Mach number
 * @param pos a vector of positions to probe
 * @return the Mach number of the flow at the position
 **/
std::vector<float> Snapshot::probeMach(std::vector<Vec> pos) {
  return probeLocation<float>(pos, "Mach");
}

/** 
 * Probes multiple locations for momentum
 * @param pos a vector of positions to probe
 * @return the momentum at positions pos
 **/
std::vector<Vec> Snapshot::probeMomentum(std::vector<Vec> pos) {
  std::vector<float> res[3];
  res[0] = probeLocation<float>(pos, "rho_vx");
  res[1] = probeLocation<float>(pos, "rho_vy");
  if (nDim == 3)
    res[2] = probeLocation<float>(pos, "rho_vz");

  // Ugly af ...
  std::vector<Vec> out(pos.size());
  for (int i=0; i < pos.size(); ++i)
    for (int j=0; j < nDim; ++j)
      out[i][j] = res[j][i];

  return out;
}

/** 
 * Probes multiple locations for velocity
 * @param pos a vector of positions to probe
 * @return the velocity at positions pos
 **/
std::vector<Vec> Snapshot::probeVelocity(std::vector<Vec> pos) {
  std::vector<Vec>   res = probeMomentum(pos);
  std::vector<float> rho = probeLocation<float>(pos, "rho");
  for (int i=0; i < pos.size(); ++i) {
    for (int j=0; j < nDim; ++j)
      res[i][j] /= rho[i];
  }
  return res;
}

/** 
 * Probes multiple locations for the refinement level
 * @param pos a vector of positions to probe
 * @return the refinement level at positions pos
 **/
std::vector<int> Snapshot::probeLevel(std::vector<Vec> pos) {
  return probeLocation<int>(pos, "level");
}

/** 
 * Probes multiple locations for the mpi-rank
 * @param pos a vector of positions to probe
 * @return the mpi-rank at positions pos
 **/
std::vector<int> Snapshot::probeRank(std::vector<Vec> pos) {
  return probeLocation<int>(pos, "rank");
}

/**
 * Probes multiple locations for the octant index
 * @param pos a vector of positions to probe
 * @return the octant indices at positions pos
 **/
std::vector<int> Snapshot::probeOctant(std::vector<Vec> pos) {
  return probeLocation<int>(pos, "iOct");
}

/**
 * Returns the positions corresponding to the center of the cells 
 * traversed by the points along pos
 * @param pos a vector of positions
 * @return the center of the unique cells along the trajectory pos
 **/
std::vector<Vec> Snapshot::getUniqueCells(std::vector<Vec> pos) {
  std::vector<int> iCells = getCellsFromPositions(pos);
  std::vector<Vec> positions;

  // We put in the first element
  int last = iCells[0];
  positions.push_back(getCellCenter(last));
  for (int i=1; i < iCells.size(); ++i) {
    int iCell = iCells[i];
    if (iCell != last) {
      last = iCell;
      positions.push_back(getCellCenter(iCell));
    }
  }

  return positions;
}

/**
 * Returns the ids of all the cells corresponding to oct iOct.
 * @param iOct the id of the octant to retrieve
 * @return a vector of indices corresponding to the cell ids in octant iOct
 **/
std::vector<int> Snapshot::getBlock(uint iOct) {
  std::vector<int> out;
  if (attributes.count("iOct") == 0) {
    std::cerr << "ERROR : Cannot retrieve block as octant information is not available" << std::endl;
    std::cerr << "        in this run. To get octant information, make sure that " << std::endl;
    std::cerr << "        output/writeVariables has iOct defined in your configuration file" << std::endl;
    return out;
  }
  
  return out;
}

/**
 * Returns the value of the refinement criterion at position pos
 * @param pos the position where to probe the refinement criterion
 * @return a floating point value corresponding to the error for refinement at pos.
 *         
 * @note Result will be 0.0f if pos is at the edge of the domain
 * 
 * @todo abstract to any variable type
 * @todo abstract to any refinement error calculation
 **/
float Snapshot::getRefinementCriterion(Vec pos) {
  // Retrieving cells, sizes and domain bounding box
  BoundingBox bb = getDomainBoundingBox();
  uint iCell = getCellFromPosition(pos);
  Vec  h     = getCellSize(iCell);

  // Probing current location
  float rho = probeDensity(pos);
  float en  = probeTotalEnergy(pos);

  // Spatial offsets
  Vec off_x {h[0]*0.75f, 0.0f, 0.0f};
  Vec off_y {0.0f, h[1]*0.75f, 0.0f};

  Vec pxm = pos - off_x;
  Vec pxp = pos + off_x;
  Vec pym = pos - off_y;
  Vec pyp = pos + off_y;

  // We check that the positions are in the domain
  if ( !inBoundingBox(bb, pxm, 1) || !inBoundingBox(bb, pxp, 1)
    || !inBoundingBox(bb, pym, 2) || !inBoundingBox(bb, pyp, 2))
    return 0.0f;

  // And in 3D
  Vec pzp, pzm;
  if (nDim == 3) {
    Vec off_z {0.0f, 0.0, h[2]*0.75f};
    pzp = pos + off_z;
    pzm = pos - off_z;

    if (!inBoundingBox(bb, pzm, 3) || !inBoundingBox(bb, pzp, 3))
      return 0.0f;
  }
  
  // Probing density values on adjacent cells in 2D
  float rho_xp = probeDensity(pxp);
  float rho_xm = probeDensity(pxm);
  float rho_yp = probeDensity(pym);
  float rho_ym = probeDensity(pyp);

  // And energy
  float en_xp = probeTotalEnergy(pxp);
  float en_xm = probeTotalEnergy(pxm);
  float en_yp = probeTotalEnergy(pyp);
  float en_ym = probeTotalEnergy(pym);
  
  // Calculating error
  auto err_calc = [&] (float ui, float uim, float uip, float eps=0.01) {
    float grad_L = fabs(ui - uim);
    float grad_R = fabs(uip - ui);
    float cd = fabs(2.0*ui) + fabs(uip) + fabs(uim);
    float d2 = fabs(uip + uim - 2.0f*ui);
    float Ei = d2 / (grad_L + grad_R + eps * cd);
    return Ei;
  };

  float err = std::max({err_calc(rho, rho_xm, rho_xp),
                        err_calc(rho, rho_ym, rho_yp),
                        err_calc(en, en_xm, en_xp),
                        err_calc(en, en_ym, en_yp)});
  if (nDim == 3) {
    float rho_zp = probeDensity(pzp);
    float rho_zm = probeDensity(pzm);
    float en_zp = probeTotalEnergy(pzp);
    float en_zm = probeTotalEnergy(pzm);

    err = std::max({err,
                    err_calc(rho, rho_zm, rho_zp),
                    err_calc(en, en_zm, en_zp)});
  }

  return err;
}

/**
 * Returns the integrated total mass over the domain
 * @return the total mass in the domain
 **/
float Snapshot::getTotalMass() {
  // We require the density and the volume of each cell in the domain
  // which means 8 bytes per cell (1 float for volume, 1 for density)
  // We read everything by vectors to avoid filling the memory 
  float total_mass = 0.0f;
  std::vector<float> densities;
  std::vector<float> cell_volumes;

  int base_id = 0;
  int end_id;
  while (base_id < nCells) {
    end_id = std::min(base_id + vec_size, nCells);

    // Filling the id array
    std::vector<uint> cid;
    for (int i=base_id; i < end_id; ++i)
      cid.push_back(i);

    densities = probeCells<float>(cid, "rho");
    cell_volumes = getCellVolume(cid);

    uint nV = end_id - base_id;
    for (int i=0; i<nV; ++i)
      total_mass += densities[i] * cell_volumes[i];

    base_id += vec_size;
  }

  return total_mass;
}

/**
 * Returns the integrated total energy over the domain
 * @return the total mass in the domain
 **/
float Snapshot::getTotalEnergy() {
  // We require the density and the volume of each cell in the domain
  // which means 8 bytes per cell (1 float for volume, 1 for density)
  // We read everything by vectors to avoid filling the memory 
  float total_energy = 0.0f;
  std::vector<float> energies;
  std::vector<float> density;
  std::vector<float> cell_volumes;

  int base_id = 0;
  int end_id;
  while (base_id < nCells) {
    end_id = std::min(base_id + vec_size, nCells);
    
    // Filling the id array
    std::vector<uint> cid;
    for (int i=base_id; i < end_id; ++i)
      cid.push_back(i);

    energies = getTotalEnergy(cid);
    cell_volumes = getCellVolume(cid);

    uint nV = end_id - base_id;
    for (int i=0; i<nV; ++i)
      total_energy += energies[i] * cell_volumes[i];

    base_id += vec_size;
  }

  return total_energy;  
}

/**
 * Returns the integrated total internal energy over the domain
 * @return the total internal energy density in the domain
 **/
float Snapshot::getTotalInternalEnergy(double gamma) {
  // We require the density and the volume of each cell in the domain
  // which means 8 bytes per cell (1 float for volume, 1 for density)
  // We read everything by vectors to avoid filling the memory 
  float total_energy = 0.0f;
  std::vector<float> pressures;
  std::vector<float> cell_volumes;

  BoundingBox bb = getDomainBoundingBox();
  float tot_vol = 1.0;
  for (int i=0; i < nDim; ++i)
    tot_vol *= bb.second[i] - bb.first[i];

  int base_id = 0;
  int end_id;
  while (base_id < nCells) {
    end_id = std::min(base_id + vec_size, nCells);
    
    // Filling the id array
    std::vector<uint> cid;
    for (int i=base_id; i < end_id; ++i)
      cid.push_back(i);

    pressures = getPressure(cid);
    cell_volumes = getCellVolume(cid);

    uint nV = end_id - base_id;
    for (int i=0; i<nV; ++i)
      total_energy += cell_volumes[i] * pressures[i] / (gamma - 1.0);

    base_id += vec_size;
  }

  return total_energy;  
}

/**
 * Returns the integrated kinetic energy over the domain
 * @return the total kinetic energy in the domain
 **/
float Snapshot::getTotalKineticEnergy() {
  float total_Ek = 0.0;
  std::vector<float> densities;
  std::vector<Vec> momenta;
  std::vector<float> cell_volumes;

  int base_id = 0;
  int end_id;
  while (base_id < nCells) {
    end_id = std::min(base_id + vec_size, nCells);

    // Filling the id array
    std::vector<uint> cid;
    for (int i=base_id; i < end_id; ++i)
      cid.push_back(i);

    densities = getDensity(cid);
    momenta   = getMomentum(cid);
    cell_volumes = getCellVolume(cid);
    
    uint nV = end_id - base_id;
    auto norm2 = [](Vec v) {
      return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
    };

    for (int i=0; i<nV; ++i)
      total_Ek += 0.5 * cell_volumes[i] * norm2(momenta[i]) / densities[i];

    base_id += vec_size;
  }

  return total_Ek;
}

/**
 * Returns the maximum Mach number of the domain
 * @return the maximum Mach number of the flowin the domain
 **/
float Snapshot::getMaxMach() {
  float max_Mach = 0.0;
  std::vector<float> Mach;

  int base_id = 0;
  int end_id;
  while (base_id < nCells) {
    end_id = std::min(base_id + vec_size, nCells);

    // Filling the id array
    std::vector<uint> cid;
    for (int i=base_id; i < end_id; ++i)
      cid.push_back(i);

    Mach = probeCells<float>(cid, "Mach");
    max_Mach = std::max(max_Mach, *std::max_element(Mach.begin(), Mach.end()));
    
    base_id += vec_size;
  }

  return max_Mach;
}

/**
 * Returns the average Mach number of the domain
 * @return the average Mach number of the flowin the domain
 **/
float Snapshot::getAverageMach() {
  float avg_Mach = 0.0;
  std::vector<float> Mach;

  int base_id = 0;
  int end_id;
  while (base_id < nCells) {
    end_id = std::min(base_id + vec_size, nCells);

    // Filling the id array
    std::vector<uint> cid;
    for (int i=base_id; i < end_id; ++i)
      cid.push_back(i);

    Mach = probeCells<float>(cid, "Mach");
    avg_Mach += std::accumulate(Mach.begin(), Mach.end(), 0.0);
    
    base_id += vec_size;
  }

  return avg_Mach / nCells;
}

/**
 * Returns the value of the refinement criterion at a set of positions
 * @param pos the position vector where to probe the refinement criterion
 * @return a vector of floating point values corresponding to the error for
 *         refinement at each position. 
 * 
 * @note Result will be 0.0f for each position at the edge of the domain
 * 
 * @todo abstract to any variable type
 * @todo abstract to any refinement error calculation
 **/
std::vector<float> Snapshot::getRefinementCriterion(std::vector<Vec> pos) {
  uint nPos = pos.size(); 
  std::vector<float> out(nPos);

  #pragma omp parallel for schedule(dynamic)
  for (int i=0; i < nPos; ++i)
    out[i] = getRefinementCriterion(pos[i]);
  
  return out;
}

/**
 * Extracts density from a list of cells
 * @param iCells the ids of the cells to extract
 * @return a vector of densities for each cell
 **/
std::vector<float> Snapshot::getDensity(std::vector<uint> iCells) {
  return probeCells<float>(iCells, "rho");
}

/**
 * Extracts the pressure from a list of cells
 * @param iCells the ids of the cells to extract
 * @return a vector of pressures for each cell
 **/
std::vector<float> Snapshot::getPressure(std::vector<uint> iCells) {
  return probeCells<float>(iCells, "P");
}

/**
 * Extracts the total energy from a list of cells
 * @param iCells the ids of the cells to extract
 * @return a vector of energies for each cell
 **/
std::vector<float> Snapshot::getTotalEnergy(std::vector<uint> iCells) {
  return probeCells<float>(iCells, "e_tot");
}

/**
 * Extracts the Mach number from a list of cells
 * @param iCells the ids of the cells to extract
 * @return a vector of Mach number corresponding to the flow each cell
 **/
std::vector<float> Snapshot::getMach(std::vector<uint> iCells) {
  return probeCells<float>(iCells, "Mach");
}

/**
 * Extracts the momentum from a list of cells
 * @param iCells the ids of the cells to extract
 * @return a vector of momenta for each cell
 **/
std::vector<Vec> Snapshot::getMomentum(std::vector<uint> iCells) {
  std::vector<float> res[3];
  res[0] = probeCells<float>(iCells, "rho_vx");
  res[1] = probeCells<float>(iCells, "rho_vy");
  if (nDim == 3)
    res[2] = probeCells<float>(iCells, "rho_vz");

  // Ewwwww ...
  std::vector<Vec> out(iCells.size());
  for (uint i=0; i < iCells.size(); ++i)
    for (int j=0; j < nDim; ++j)
      out[i][j] = res[j][i];

  return out;
}

/**
 * Extracts the velocities from a list of cells
 * @param iCells the ids of the cells to extract
 * @return a vector of velocities for each cell
 **/
std::vector<Vec> Snapshot::getVelocity(std::vector<uint> iCells) {
  std::vector<Vec> momentum = getMomentum(iCells);
  std::vector<float> density = getDensity(iCells);
  std::vector<Vec> out(iCells.size());

  for (uint i=0; i < iCells.size(); ++i)
    for (uint j=0; j < nDim; ++j)
      out[i][j] = momentum[i][j] / density[i];

  return out;
}

/**
 * Extracts the AMR level from a list of cells
 * @param iCells the ids of the cells to extract
 * @return a vector of levels for each cell
 **/
std::vector<int> Snapshot::getLevel(std::vector<uint> iCells) {
  return probeCells<int>(iCells, "level");
}

/**
 * Extracts the MPI rank from a list of cells
 * @param iCells the ids of the cells to extract
 * @return a vector of rank for each cell
 **/
std::vector<int> Snapshot::getRank(std::vector<uint> iCells) {
  return probeCells<int>(iCells, "rank");
}

/**
 * Extracts the mesh octant from a list of cells
 * @param iCells the ids of the cells to extract
 * @return a vector of octant ids for each cell
 **/
std::vector<int> Snapshot::getOctant(std::vector<uint> iCells) {
  return probeCells<int>(iCells, "iOct");
}

/**
 * Reads all the data corresponding to a field
 * @param attribute the name of the field to read
 * @param sort indicates if the array should be sorted by dimensions
 * @return an array of float correponding to the linearized array of the field
 **/
std::vector<float> Snapshot::readAllFloat(std::string attribute) {
  std::vector<float> res;

  if (attributes.count(attribute) == 0) {
    std::cerr << "ERROR : Attribute " << attribute << " is not stored in file !" << std::endl;
    return res;
  }
  Attribute &att = attributes[attribute];
  hid_t data_type = type_corresp[att.type];
  if (data_type != H5T_NATIVE_FLOAT) {
    std::cerr << "ERROR : Datatype of " << attribute << " is not float !" << std::endl;
    return res;
  }

  herr_t status;
  hid_t space = H5Dget_space(att.handle);

  hsize_t *select = new hsize_t[nCells];
  for (int i=0; i < nCells; ++i) {
    select[i] = i;
  }

  status = H5Sselect_elements(space, H5S_SELECT_SET, nCells, select);
  delete [] select;

  hsize_t d=nCells;
  float* values = new float[nCells];
  hid_t memspace = H5Screate_simple(1, &d, NULL); 

  status = H5Dread(att.handle, data_type, memspace, space, H5P_DEFAULT, values);
  res.resize(nCells);

  std::copy(values, values+nCells, res.begin());

  delete [] values;
  return res;
}

/**
 * Sorts an array in 2D extracted from readAllXXXX along dimensions
 * using morton index.
 * @param vec the vector to sort
 * @param iLevel the level at which we are (only works for regular grids)
 * @param bx number of blocks along x
 * @param by number of blocks along y
 * @result the sorted vector
 **/
std::vector<float> Snapshot::mortonSort2d(std::vector<float> vec,
                          uint iLevel, 
                          uint bx, uint by) {
  uint nOctPerDim = (1 << iLevel);
  uint nOct = nOctPerDim*nOctPerDim;
  uint chunkSize = bx*by;
  uint Nx = bx * nOctPerDim;
  uint Ny = by * nOctPerDim;

  std::vector<float> result;

  if (chunkSize * nOct != nCells) {
    std::cerr << "ERROR : Number of cells is not corresponding to level/block size for Morton sort !" << std::endl;
    std::cerr << " . nCells = " << nCells << std::endl;
    std::cerr << " . cell count = " << chunkSize * nOct << "(chunk=" << chunkSize << "; nOct=" << nOct << ")" << std::endl;
    return result;
  }

  // Morton encoding
  auto splitBy2 = [](uint32_t z) {
    uint64_t x = z & 0xffffffff;
    x = (x | x << 16) & 0xffff0000ffff;
    x = (x | x << 8) & 0xff00ff00ff00ff;
    x = (x | x << 4) & 0xf0f0f0f0f0f0f0f;
    x = (x | x << 2) & 0x3333333333333333;
    x = (x | x << 1) & 0x5555555555555555;

    return x;
  };

  for (uint iy=0; iy < Ny; ++iy) {
    uint64_t iOy = iy / by;
    uint iBy = iy % by;
    for (uint ix=0; ix < Nx; ++ix) {
      uint64_t iOx = ix / bx;
      uint iBx = ix % bx;

      uint64_t key = splitBy2(iOx) | splitBy2(iOy) << 1;
      key *= chunkSize;
      key += iBy * bx + iBx;
      result.push_back(vec[key]);
    }
  }
  return result;
}

std::vector<float> Snapshot::mortonSort3d(std::vector<float> vec,
                                          uint iLevel, 
                                          uint bx, uint by, uint bz) {
  uint nOctPerDim = (1 << iLevel);
  uint nOct = nOctPerDim*nOctPerDim*nOctPerDim;
  uint chunkSize = bx*by*bz;
  uint Nx = bx * nOctPerDim;
  uint Ny = by * nOctPerDim;
  uint Nz = bz * nOctPerDim;

  std::vector<float> result;

  if (chunkSize * nOct != nCells) {
    std::cerr << "ERROR : Number of cells is not corresponding to level/block size for Morton sort !" << std::endl;
    std::cerr << " . nCells = " << nCells << std::endl;
    std::cerr << " . cell count = " << chunkSize * nOct << std::endl;
    return result;
  }

  // Morton encoding
  auto splitBy3 = [](uint32_t z) {
    uint64_t x = z & 0xffffffff;
    x = (x | x << 32) & 0x1f00000000ffff; 
    x = (x | x << 16) & 0x1f0000ff0000ff; 
    x = (x | x << 8) & 0x100f00f00f00f00f;
    x = (x | x << 4) & 0x10c30c30c30c30c3; 
    x = (x | x << 2) & 0x1249249249249249;

    return x;
  };

  uint bxby = bx*by;
  for (uint iz=0; iz < Nz; ++iz) {
    uint64_t iOz = iz / by;
    uint iBz = iz % by;
    for (uint iy=0; iy < Ny; ++iy) {
      uint64_t iOy = iy / by;
      uint iBy = iy % by;
      for (uint ix=0; ix < Nx; ++ix) {
        uint64_t iOx = ix / bx;
        uint iBx = ix % bx;

        uint64_t key = splitBy3(iOx) | splitBy3(iOy) << 1 | splitBy3(iOz) << 2;
        key *= chunkSize;
        key += iBz * bxby + iBy * bx + iBx;
        result.push_back(vec[key]);
      }
    }
  }
  return result;
}

std::vector<uint64_t> Snapshot::getSortingMask2d(uint iLevel, uint bx, uint by) {
  std::vector<uint64_t> result;

  uint nOctPerDim = (1 << iLevel);
  uint nOct = nOctPerDim*nOctPerDim;
  uint chunkSize = bx*by;
  uint Nx = bx * nOctPerDim;
  uint Ny = by * nOctPerDim;

  if (chunkSize * nOct != nCells) {
    std::cerr << "ERROR : Number of cells is not corresponding to level/block size for Morton sort !" << std::endl;
    std::cerr << " . nCells = " << nCells << std::endl;
    std::cerr << " . cell count = " << chunkSize * nOct << std::endl;
    return result;
  }
  
  // Morton encoding
  auto splitBy2 = [](uint32_t z) {
    uint64_t x = z & 0xffffffff;
    x = (x | x << 16) & 0xffff0000ffff;
    x = (x | x << 8) & 0xff00ff00ff00ff;
    x = (x | x << 4) & 0xf0f0f0f0f0f0f0f;
    x = (x | x << 2) & 0x3333333333333333;
    x = (x | x << 1) & 0x5555555555555555;

    return x;
  }; 

  for (uint iy=0; iy < Ny; ++iy) {
    uint64_t iOy = iy / by;
    uint iBy = iy % by;
    for (uint ix=0; ix < Nx; ++ix) {
      uint64_t iOx = ix / bx;
      uint iBx = ix % bx;

      uint64_t key = splitBy2(iOx) | splitBy2(iOy) << 1;
      key *= chunkSize;
      key += iBy * bx + iBx;
      result.push_back(key);
    }
  }
  return result;
}

std::vector<uint64_t> Snapshot::getSortingMask3d(uint iLevel, uint bx, uint by, uint bz) {
  std::vector<uint64_t> result;

  uint nOctPerDim = (1 << iLevel);
  uint nOct = nOctPerDim*nOctPerDim*nOctPerDim;
  uint chunkSize = bx*by*bz;
  uint Nx = bx * nOctPerDim;
  uint Ny = by * nOctPerDim;
  uint Nz = bz * nOctPerDim;
  
  // Morton encoding
  auto splitBy3 = [](uint32_t z) {
    uint64_t x = z & 0xffffffff;
    x = (x | x << 32) & 0x1f00000000ffff; 
    x = (x | x << 16) & 0x1f0000ff0000ff; 
    x = (x | x << 8) & 0x100f00f00f00f00f;
    x = (x | x << 4) & 0x10c30c30c30c30c3; 
    x = (x | x << 2) & 0x1249249249249249;

    return x;
  }; 

  uint bxby = bx*by;
  for (uint iz=0; iz < Nz; ++iz) {
    uint64_t iOz = iz / by;
    uint iBz = iz % by;
    for (uint iy=0; iy < Ny; ++iy) {
      uint64_t iOy = iy / by;
      uint iBy = iy % by;
      for (uint ix=0; ix < Nx; ++ix) {
        uint64_t iOx = ix / bx;
        uint iBx = ix % bx;

        uint64_t key = splitBy3(iOx) | splitBy3(iOy) << 1 | splitBy3(iOz) << 2;
        key *= chunkSize;
        key += iBz * bxby + iBy * bx + iBx;
        result.push_back(key);
      }
    }
  }
  return result;
}


}