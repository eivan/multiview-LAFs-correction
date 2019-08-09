#pragma once

#include <vector>

#include "types.hpp"

// References:
//
// [1] I.Eichhardt, D.Barath, Optimal Multi-view Correction of Local Affine Frames.
//     In Proc. British Machine Vision Conf., 2019.
//
// [2] D.Barath, L.Hajder, and J.Matas. Accurate closed-form estimation of
//     local affine transformations consistent with the epipolar geometry.
//     In Proc. British Machine Vision Conf., 2016.

namespace higher_order::stereo {

  // Creates first order constraint on a pair of local affine frames.
  // Note: This method works for perspective views only.
  // Input: the fundamental matrix (F) and image points (p1, p2).
  // Output: vector-coefficients (a, b) of the constraint.
  // See Eq. (2) and (3) of [1]
  void CreateFirstOrderConstraint(
    const Mat3d& F,
    const Vec2d& p1, const Vec2d& p2,
    Vec2d& a, Vec2d& b);

  // Creates first order constraint on a pair of local affine frames.
  // Note: This method works for central views.
  // Input: the essential matrix (F) and bearing vectors (q1, q2) with corresponding gradients.
  // Output: vector-coefficients (a, b) of the constraint.
  // See Eq. (2) and (3) of [1]
  void CreateFirstOrderConstraint(
    const Mat3d& E,
    const Vec3d& q1, const Vec3d& q2,
    const Mat32d& nabla_q1, const Mat32d& nabla_q2,
    Vec2d& a, Vec2d& b);

  // A compact and efficient solution to [2]
  void RefineAffineCorrespondence(
    Mat2d& A,
    const Vec2d& a, const Vec2d& b);

} // higher_order::stereo

namespace higher_order::multiview {

  // L2-Optimal Multi-view Correction of Local Affine Frames
  void RefineLocalAffineFrames(
    const Eigen::MatrixXd& epipolar_constraints,
    MatX2d& local_affine_frames);

  // L2-Optimal Multi-view Correction of Local Affine Frames
  void RefineLocalAffineFrames(
    const EpipolarConstraints& epipolar_constraints,
    ObservedLAFs& local_affine_frames);

} // higher_order::multiview

