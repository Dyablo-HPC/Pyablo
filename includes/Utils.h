#pragma once

#include <vector>
#include <array>
#include <string>
#include <map>
#include <stdexcept>

#include "hdf5.h"

namespace dyablo {

using Vec = std::array<float, 3>;
using BoundingBox = std::pair<Vec, Vec>;

using RealArray   = std::vector<float>;
using VecArray    = std::vector<Vec>;
using IntArray    = std::vector<int>;
using UIntArray   = std::vector<uint>;
using BoolArray   = std::vector<bool>;
using UInt64Array = std::vector<uint64_t>;
using RealTable   = std::vector<RealArray>;
using VecTable    = std::vector<VecArray>;

struct Attribute {
  std::string name;
  std::string type;
  std::string center;
  hid_t handle;
};

enum Direction : uint8_t {
  IX = 0,
  IY = 1,
  IZ = 2
};

struct VarState {
  float rho;
  float vx;
  float vy;
  float vz;
  float prs; 
};

struct Line {
  int Nl;
  Vec start, end;
  VecArray pos;
  RealArray rho, prs, E;
  VecArray vel;
  UIntArray cellIds;
};

enum class SliceDir {
  XY,
  XZ, 
  YZ
};

struct Slice {
  int Nx, Ny;
  SliceDir dir;
  float origin;

  VecArray pos;
  RealArray rho, prs, E;
  RealArray vx, vy, vz;
};

template<typename Tout, typename Tin>
Tout reshapeArray(const Tin &array, int Nx, int Ny);

/** Basic algebraic operations on vectors **/
Vec operator+(const Vec &v1, const Vec &v2);
Vec operator-(const Vec &v1, const Vec &v2);
Vec& operator+=(Vec &v1, const Vec &v2);
Vec& operator-=(Vec &v1, const Vec &v2);
Vec operator*(const Vec &v, float q);
Vec& operator*=(Vec &v, float q);
Vec operator/(const Vec &v, float q);
Vec& operator/=(Vec &v, float q);

/** Bounding box helpers **/
bool inBoundingBox(BoundingBox bb, Vec pos, int nDim);

class HDF5_store{
  std::map<std::string, hid_t> file_handles;
  std::map<std::string, hid_t> data_handles;
public:
  ~HDF5_store()
  {
    for(const auto& p : file_handles)
      H5Fclose( p.second );
    file_handles.clear();
    for(const auto& p : data_handles)
      H5Dclose( p.second );
    data_handles.clear();
  }

  hid_t open_file( const std::string& handle, const std::string& filename )
  {
    if( file_handles.count(handle) == 0 )
    {
      hid_t hid = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      if (hid < 0) {
        throw std::runtime_error( std::string("ERROR : Could not open file ") + filename);
      }
      file_handles[handle] = hid;
    }
    return file_handles.at(handle);
  }

  hid_t get_file(const std::string& handle)
  {
    return file_handles.at(handle);
  }

  hid_t open_data( const std::string& file_handle, const std::string& xpath )
  {
    std::string data_handle = file_handle + ":" + xpath;

    if( data_handles.count(data_handle) == 0 )
    {
      hid_t data_id = H5Dopen2(get_file(file_handle), xpath.c_str(), H5P_DEFAULT);
      if (data_id < 0) {
        throw std::runtime_error( std::string("ERROR : Could not open data ") + data_handle);
      }
      data_handles[data_handle] = data_id;
    }

    return data_handles.at(data_handle);
  }



};
}