#include "Snapshot.h"
#include <omp.h>

namespace dyablo {

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
  if (!handles.contains(handle)) {
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

  if (data_type == H5T_NATIVE_INT) {
    int value;
    herr_t status = H5Dread(att.handle, data_type, memspace, space, H5P_DEFAULT, &value);
    return value;
  }
  else if (data_type == H5T_NATIVE_FLOAT) {
    float value;
    herr_t status = H5Dread(att.handle, data_type, memspace, space, H5P_DEFAULT, &value);
    return value;
  }

  return 0;
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
  // Retrieving cell indices and adding them to selection array
  std::cout << "Probing multiple locations [this might be long] ..." << std::endl;
  int nPos = pos.size();
  std::vector<int> iCells = getCellsFromPositions(pos);

  hsize_t select[nPos][1];
  for (int i=0; i < nPos; ++i)
    select[i][0] = (hsize_t)iCells[i];

  Attribute &att = attributes[attribute];
  hid_t data_type = type_corresp[att.type];

  // Selecting the element in the dataset
  herr_t status;
  hid_t space = H5Dget_space(att.handle);
  status = H5Sselect_elements(space, H5S_SELECT_SET, nPos, (const hsize_t *)&select);

  hsize_t d=nPos;
  hid_t memspace = H5Screate_simple(1, &d, NULL);

  std::vector<T> out;
  out.resize(nPos);

  if (data_type == H5T_NATIVE_INT) {
    herr_t status = H5Dread(att.handle, data_type, memspace, space, H5P_DEFAULT, out.data());
  }
  else if (data_type == H5T_NATIVE_FLOAT) {
    herr_t status = H5Dread(att.handle, data_type, memspace, space, H5P_DEFAULT, out.data());
  }

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
 * Probes multiple locations for density
 * @param pos a vector of positions to probe
 * @return the density at positions pos
 **/
std::vector<float> Snapshot::probeDensity(std::vector<Vec> pos) {
  return probeLocation<float>(pos, "rho");
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
  std::vector<Vec> out;
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

}