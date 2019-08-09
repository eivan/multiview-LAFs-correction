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