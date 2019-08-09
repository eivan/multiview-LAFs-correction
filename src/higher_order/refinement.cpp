#include "refinement.hpp"

#include <vector>
#include <unordered_map>

void higher_order::stereo::CreateFirstOrderConstraint(
  const Mat3d& F,
  const Vec2d& p1, const Vec2d& p2,
  Vec2d& a, Vec2d& b)
{
  // See Eq. (2) and (3) of [1]
  a = (F * p1.homogeneous()).head<2>();
  b = (F.transpose() * p2.homogeneous()).head<2>();
}

void higher_order::stereo::CreateFirstOrderConstraint(
  const Mat3d& E,
  const Vec3d& q1, const Vec3d& q2,
  const Mat32d& nabla_q1, const Mat32d& nabla_q2,
  Vec2d& a, Vec2d& b)
{
  // See Eq. (2) and (3) of [1]
  a = nabla_q2.transpose() * E * q1;
  b = nabla_q1.transpose() * E.transpose() * q2;
}

void higher_order::stereo::RefineAffineCorrespondence(
  Mat2d& A,
  const Vec2d& a, const Vec2d& b)
{
  A -= (1 / a.dot(a)) * a * (a.transpose() * A + b.transpose());
}