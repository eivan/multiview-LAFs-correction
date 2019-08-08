#pragma once

#include <higher_order/types.hpp>

using namespace Eigen;
using namespace higher_order;

using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;
using Mat2 = Eigen::Matrix2d;
using Mat3 = Eigen::Matrix3d;

inline Vec2 project(const Mat3& K, const Vec3& X) {
  Vec2 hnorm = X.hnormalized();
  return K.topLeftCorner<2, 2>() * hnorm + K.topRightCorner<2, 1>();
}

inline Mat23 project_gradient(const Mat3& K, const Vec3& X) {
  Vec2 hnorm = X.hnormalized();
  Mat23 hnorm_gradient = (Mat23() <<
    1 / X.z(), 0, -hnorm.x() / X.z(),
    0, 1 / X.z(), -hnorm.y() / X.z()).finished();
  return K.topLeftCorner<2, 2>() * hnorm_gradient;
}


inline Mat3 LookAt_direction(const Vector3d& direction, const Vector3d& up = Vector3d::UnitY()) {
  Mat3 R;
  R.row(2) = direction.normalized(); // z
  R.row(0) = up.cross(R.row(2)).normalized(); // x
  R.row(1) = R.row(2).cross(R.row(0)); // y
  return R;
}

inline Mat3 LookAt_target(const Vec3& eye, const Vec3& target, const Vec3& up = Vec3::UnitY()) {
  return LookAt_direction(Vec3(target - eye), up);
}

// <openMVG/multiview/essential.cpp>
inline void RelativeCameraMotion(const Mat3& R1,
  const Vec3& t1,
  const Mat3& R2,
  const Vec3& t2,
  Mat3* R,
  Vec3* t) {
  *R = R2 * R1.transpose();
  *t = t2 - (*R) * t1;
}

// <openMVG/numeric/numeric.cpp>
inline Mat3 CrossProductMatrix(const Vec3& x) {
  Mat3 X;
  X << 0, -x(2), x(1),
    x(2), 0, -x(0),
    -x(1), x(0), 0;
  return X;
}

// <openMVG/multiview/essential.cpp>
inline void EssentialFromRt(const Mat3& R1,
  const Vec3& t1,
  const Mat3& R2,
  const Vec3& t2,
  Mat3* E) {
  Mat3 R;
  Vec3 t;
  RelativeCameraMotion(R1, t1, R2, t2, &R, &t);
  const Mat3 Tx = CrossProductMatrix(t);
  *E = Tx * R;
}

template <typename T>
inline Matrix<T, 3, 2> Nullspace(const Matrix<T, 3, 1> & normal) {
  return normal.jacobiSvd(Eigen::ComputeFullU).matrixU().block<3, 2>(0, 1);
}