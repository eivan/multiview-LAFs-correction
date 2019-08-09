#include "refinement.hpp"

#include <vector>
#include <unordered_map>

void higher_order::stereo::CreateFirstOrderConstraint(
  const Eigen::Matrix3d& F,
  const Eigen::Vector2d& p1, const Eigen::Vector2d& p2,
  Eigen::Vector2d& a, Eigen::Vector2d& b)
{
  // See Eq. (2) and (3) of [1]
  a = (F * p1.homogeneous()).head<2>();
  b = (F.transpose() * p2.homogeneous()).head<2>();
}

void higher_order::stereo::CreateFirstOrderConstraint(
  const Eigen::Matrix3d& E,
  const Eigen::Vector3d& q1, const Eigen::Vector3d& q2,
  const Mat32& nabla_q1, const Mat32& nabla_q2,
  Eigen::Vector2d& a, Eigen::Vector2d& b)
{
  // See Eq. (2) and (3) of [1]
  a = nabla_q2.transpose() * E * q1;
  b = nabla_q1.transpose() * E.transpose() * q2;
}

void higher_order::stereo::RefineAffineCorrespondence(
  Eigen::Matrix2d& A,
  const Eigen::Vector2d& a, const Eigen::Vector2d& b)
{
  A -= (1 / a.dot(a)) * a * (a.transpose() * A + b.transpose());
}

void higher_order::multiview::RefineLocalAffineFrames(
  const Eigen::MatrixXd& epipolar_constraints, MatX2& local_affine_frames)
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

void higher_order::multiview::RefineLocalAffineFrames(
  const EpipolarConstraints& epipolar_constraints,
  std::vector<std::pair<ViewID, Eigen::Matrix2d>>& local_affine_frames) 
{
  assert(local_affine_frames.size() >= 2);
  assert(epipolar_constraints.size() > 0);

  const size_t N_V = local_affine_frames.size();
  const size_t N_C = epipolar_constraints.size();

  std::unordered_map<ViewID, size_t> inverse_view_indices;
  for (size_t i = 0; i < local_affine_frames.size(); ++i) {
    const auto& view = local_affine_frames[i].first;
    inverse_view_indices[view] = i;
  }

  // Epipolar constraints    
  Eigen::MatrixXd C;
  C.resize(N_C, 2 * N_V);
  C.setZero();

  for (size_t i = 0; i < epipolar_constraints.size(); ++i) {
    const auto& epi = epipolar_constraints[i];
    const auto& [view_i, view_j] = epi.view_pair;

    C.block<1, 2>(i, inverse_view_indices.at(view_i) * 2) = epi.b.transpose();  // "b" for view i
    C.block<1, 2>(i, inverse_view_indices.at(view_j) * 2) = epi.a.transpose();  // "a" for view j
  }

  // Fill up right hand side with the linear transformation (as in LAFs)
  MatX2 lafs;
  lafs.resize(2 * N_V, Eigen::NoChange);
  for (size_t i = 0; i < local_affine_frames.size(); ++i) {
    lafs.block<2, 2>(i * 2, 0) = local_affine_frames[i].second;
  }

  RefineLocalAffineFrames(C, lafs);

  // Copy solution back to input
  for (size_t i = 0; i < local_affine_frames.size(); ++i) {
    local_affine_frames[i].second = lafs.block<2, 2>(i * 2, 0);
  }
}
