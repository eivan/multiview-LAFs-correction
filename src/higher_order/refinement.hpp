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

  /*!
   * Creates first order constraint on a pair of local affine frames.
   * Note: This method works for perspective views only.
   * See Eq. (2) and (3) of I. Eichhardt, D. Barath, Optimal Multi-view Correction of Local Affine Frames. In Proc. British Machine Vision Conf., 2019.
   * \param F       the fundamental matrix
   * \param p1      image point 1
   * \param p2      image point 2
   * \param[out] a  first vector coefficient of the differential constraint
   * \param[out] b  second vector coefficient of the differential constraint
   */
  void CreateFirstOrderConstraint(
    const Mat3d& F,
    const Vec2d& p1, const Vec2d& p2,
    Vec2d& a, Vec2d& b);

  /*!
   * Creates first order constraint on a pair of local affine frames.
   * Note: This method works for central views.
   * See Eq. (2) and (3) of I. Eichhardt, D. Barath, Optimal Multi-view Correction of Local Affine Frames. In Proc. British Machine Vision Conf., 2019.
   * \param F       the fundamental matrix
   * \param p1      bearing vector on image 1
   * \param p2      bearing vector on image 2
   * \param[out] a  first vector coefficient of the differential constraint
   * \param[out] b  second vector coefficient of the differential constraint
   */
  void CreateFirstOrderConstraint(
    const Mat3d& E,
    const Vec3d& q1, const Vec3d& q2,
    const Mat32d& nabla_q1, const Mat32d& nabla_q2,
    Vec2d& a, Vec2d& b);

  /*!
   * Refines the linear transformation part of an Affine Correspondence.
   * Note: This is an efficient/improved solution to D.Barath, L.Hajder, and J.Matas. Accurate closed-form estimation of local affine transformations consistent with the epipolar geometry. In Proc. British Machine Vision Conf., 2016.
   * \param[in,out] A   the linear transformation part of an Affine Correspondences
   * \param[in] a       first vector coefficient of the differential constraint
   * \param[in] b       second vector coefficient of the differential constraint
   */
  void RefineAffineCorrespondence(
    Mat2d& A,
    const Vec2d& a, const Vec2d& b);

} // higher_order::stereo

namespace higher_order::multiview {

  /*!
   * Performs L2-Optimal Multi-view Correction of Local Affine Frames.
   * See: I. Eichhardt, D. Barath, Optimal Multi-view Correction of Local Affine Frames. In Proc. British Machine Vision Conf., 2019.
   * \param[in] epipolar_constraints the matrix containing the epopolar constraints
   * \param[in,out] local_affine_frames the two-column matrix containing the linear transformation part of Local Affine Frames
   */
  inline void RefineLocalAffineFrames(
    const Eigen::MatrixXd& epipolar_constraints,
    MatX2d& local_affine_frames)
  {
    assert(local_affine_frames.rows() >= 4);
    assert(epipolar_constraints.rows() > 0);
    assert(epipolar_constraints.cols() == local_affine_frames.rows());

    auto svd = epipolar_constraints.transpose().bdcSvd(Eigen::ComputeThinU);
    // Omega is m×n, assume rank is r, then Omega = U*Sigma*V^T
    //   U is m×r 
    //   Sigma is r×r
    //   V is n×r

    const auto r = std::min({
      svd.rank(),
      epipolar_constraints.cols() - 3,
      epipolar_constraints.rows() });

    local_affine_frames -= svd.matrixU().leftCols(r) *
      (svd.matrixU().leftCols(r).transpose() * local_affine_frames);

    // If rank can be guaranteed to be epipolar_constraints.cols() - 3, then:
    /* local_affine_frames -=
         epipolar_constraints.transpose() *
         epipolar_constraints.transpose().colPivHouseholderQr().solve(local_affine_frames);*/
  }

  /*!
   * Performs L2-Optimal Multi-view Correction of Local Affine Frames.
   * See: I. Eichhardt, D. Barath, Optimal Multi-view Correction of Local Affine Frames. In Proc. British Machine Vision Conf., 2019.
   * \param[in] epipolar_constraints the container listing the vector coefficients of the differential constraints
   * \param[in,out] local_affine_frames the container of the linear transformation part of Local Affine Frames 
   */
  inline void RefineLocalAffineFrames(
    const EpipolarConstraints& epipolar_constraints,
    ObservedLAFs& local_affine_frames)
  {
    assert(local_affine_frames.size() >= 2);
    assert(epipolar_constraints.size() > 0);

    const size_t N_V = local_affine_frames.size();
    const size_t N_C = epipolar_constraints.size();

    std::unordered_map<ViewID, size_t> inverse_view_indices;
    for (size_t i = 0; i < N_V; ++i) {
      const auto view = local_affine_frames[i].viewID;
      inverse_view_indices[view] = i;
    }

    // Epipolar constraints    
    Eigen::MatrixXd C;
    C.resize(N_C, 2 * N_V);
    C.setZero();

    for (size_t i = 0; i < N_C; ++i) {
      const auto& epi = epipolar_constraints[i];
      const auto& [view_i, view_j] = epi.view_pair;

      C.block<1, 2>(i, inverse_view_indices.at(view_i) * 2) = epi.b.transpose();  // "b" for view i
      C.block<1, 2>(i, inverse_view_indices.at(view_j) * 2) = epi.a.transpose();  // "a" for view j
    }

    // Fill up right hand side with the linear transformation (as in LAFs)
    MatX2d lafs;
    lafs.resize(2 * N_V, Eigen::NoChange);
    for (size_t i = 0; i < local_affine_frames.size(); ++i) {
      lafs.block<2, 2>(i * 2, 0) = local_affine_frames[i].M;
    }

    RefineLocalAffineFrames(C, lafs);

    // Copy solution back to input
    for (size_t i = 0; i < local_affine_frames.size(); ++i) {
      local_affine_frames[i].M = lafs.block<2, 2>(i * 2, 0);
    }
  }

} // higher_order::multiview

