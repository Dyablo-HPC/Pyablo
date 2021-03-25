#pragma once

#include <bits/stdc++.h>

namespace dyablo {

using Vec = std::array<float, 3>;
using BoundingBox = std::pair<Vec, Vec>;

/** Basic algebraic operations on vectors **/
Vec operator+(const Vec &v1, const Vec &v2);
Vec operator-(const Vec &v1, const Vec &v2);
Vec& operator+=(Vec &v1, const Vec &v2);
Vec& operator-=(Vec &v1, const Vec &v2);
Vec operator*(const Vec &v, float q);
Vec& operator*=(Vec &v, float q);

/** Bounding box helpers **/
bool inBoundingBox(BoundingBox bb, Vec pos, int nDim);
}