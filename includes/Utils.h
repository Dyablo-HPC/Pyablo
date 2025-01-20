#pragma once

#include <vector>
#include <array>
#include <string>

#include "hdf5.h"

namespace dyablo {

#if 1
  #define H5T_NATIVE_REAL_T H5T_NATIVE_FLOAT
  using real_t = float;
#else
  #define H5T_NATIVE_REAL_T H5T_NATIVE_DOUBLE
  using real_t = double;
#endif

#define USE_CELL_CENTROID

using Vec = std::array<real_t, 3>;
using BoundingBox = std::pair<Vec, Vec>;

using RealArray   = std::vector<real_t>;
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
  real_t rho;
  real_t vx;
  real_t vy;
  real_t vz;
  real_t prs; 
};

struct Line {
  int Nl;
  Vec start, end;
  VecArray pos;
  RealArray rho, prs, E;
  VecArray vel;
};

enum class SliceDir {
  XY,
  XZ, 
  YZ
};

struct Slice {
  int Nx, Ny;
  SliceDir dir;
  real_t origin;

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
Vec operator*(const Vec &v, real_t q);
Vec& operator*=(Vec &v, real_t q);
Vec operator/(const Vec &v, real_t q);
Vec& operator/=(Vec &v, real_t q);

/** Bounding box helpers **/
bool inBoundingBox(BoundingBox bb, Vec pos, int nDim);
}