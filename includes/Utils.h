#pragma once

#include <vector>
#include <array>
#include <string>

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
}