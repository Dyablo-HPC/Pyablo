#pragma once

#include <vector>
#include <array>
#include <string>

#include "hdf5.h"

namespace dyablo {

using Vec = std::array<double, 3>;
using BoundingBox = std::pair<Vec, Vec>;

using RealArray   = std::vector<double>;
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
  double rho;
  double vx;
  double vy;
  double vz;
  double prs; 
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
  double origin;

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
Vec operator*(const Vec &v, double q);
Vec& operator*=(Vec &v, double q);
Vec operator/(const Vec &v, double q);
Vec& operator/=(Vec &v, double q);

/** Bounding box helpers **/
bool inBoundingBox(BoundingBox bb, Vec pos, int nDim);
}