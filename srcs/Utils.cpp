#include "Utils.h"

namespace dyablo {
  
/** Basic algebraic operations on vectors **/
Vec operator+(const Vec &v1, const Vec &v2) {
  return Vec{v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]};
};

Vec operator-(const Vec &v1, const Vec &v2) {
  return Vec{v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]};
}

Vec& operator+=(Vec &v1, const Vec &v2) {
  v1[0] += v2[0];
  v1[1] += v2[1];
  v1[2] += v2[2];
  return v1;
}

Vec& operator-=(Vec &v1, const Vec &v2) {
  v1[0] -= v2[0];
  v1[1] -= v2[1];
  v1[2] -= v2[2];
  return v1; 
}

Vec operator*(const Vec &v, float q) {
  return Vec{v[0]*q, v[1]*q, v[2]*q};
}

Vec& operator*=(Vec &v, float q) {
  v[0] *= q;
  v[1] *= q;
  v[2] *= q;
  return v;
}


/**
 * Returns if the given position is inside the bounding box
 * @param bb the bounding box
 * @param pos a position to check
 * @param nDim number of dimensions to test
 **/
bool inBoundingBox(BoundingBox bb, Vec pos, int nDim) {
  for (int i=0; i < nDim; ++i) {
    if (pos[i] < bb.first[i] || pos[i] > bb.second[i])
      return false;
  }
  return true;
}
}