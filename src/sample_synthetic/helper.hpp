// MIT License
// 
// Copyright (c) 2019 Ivan Eichhardt
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
// Please contact the author of this library if you have any questions.
// Author: Ivan Eichhardt (ivan.eichhardt@sztaki.mta.hu)
#pragma once

#include <higher_order/types.hpp>

using namespace Eigen;
using namespace higher_order;

inline Vec2d project(const Mat3d& K, const Vec3d& X) {
  Vec2d hnorm = X.hnormalized();
  return K.topLeftCorner<2, 2>() * hnorm + K.topRightCorner<2, 1>();
}

inline Mat23d project_gradient(const Mat3d& K, const Vec3d& X) {
  Vec2d hnorm = X.hnormalized();
  Mat23d hnorm_gradient = (Mat23d() <<
    1 / X.z(), 0, -hnorm.x() / X.z(),
    0, 1 / X.z(), -hnorm.y() / X.z()).finished();
  return K.topLeftCorner<2, 2>() * hnorm_gradient;
}


inline Mat3d LookAt_direction(const Vector3d& direction, const Vector3d& up = Vector3d::UnitY()) {
  Mat3d R;
  R.row(2) = direction.normalized(); // z
  R.row(0) = up.cross(R.row(2)).normalized(); // x
  R.row(1) = R.row(2).cross(R.row(0)); // y
  return R;
}

inline Mat3d LookAt_target(const Vec3d& eye, const Vec3d& target, const Vec3d& up = Vec3d::UnitY()) {
  return LookAt_direction(Vec3d(target - eye), up);
}

// <openMVG/multiview/essential.cpp>
inline void RelativeCameraMotion(const Mat3d& R1,
  const Vec3d& t1,
  const Mat3d& R2,
  const Vec3d& t2,
  Mat3d* R,
  Vec3d* t) {
  *R = R2 * R1.transpose();
  *t = t2 - (*R) * t1;
}

// <openMVG/numeric/numeric.cpp>
inline Mat3d CrossProductMatrix(const Vec3d& x) {
  Mat3d X;
  X << 0, -x(2), x(1),
    x(2), 0, -x(0),
    -x(1), x(0), 0;
  return X;
}

// <openMVG/multiview/essential.cpp>
inline void EssentialFromRt(const Mat3d& R1,
  const Vec3d& t1,
  const Mat3d& R2,
  const Vec3d& t2,
  Mat3d* E) {
  Mat3d R;
  Vec3d t;
  RelativeCameraMotion(R1, t1, R2, t2, &R, &t);
  const Mat3d Tx = CrossProductMatrix(t);
  *E = Tx * R;
}

template <typename T>
inline Matrix<T, 3, 2> Nullspace(const Matrix<T, 3, 1> & normal) {
  return normal.jacobiSvd(Eigen::ComputeFullU).matrixU().block<3, 2>(0, 1);
}