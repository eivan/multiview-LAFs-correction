#include "structure_estimation.hpp"

Eigen::Vector3d higher_order::stereo::EstimateNormal_FNE(
  const Mat23& gradient1, const Mat23& gradient2,
  const Eigen::Matrix2d& affine) {

  const  Eigen::Vector3d p = affine(1, 1) * gradient1.row(1).cross(gradient2.row(0)) -
    affine(0, 0) * gradient2.row(1).cross(gradient1.row(0));
  const  Eigen::Vector3d q = affine(1, 0) * gradient2.row(0).cross(gradient1.row(0)) -
    affine(0, 1) * gradient1.row(1).cross(gradient2.row(1));
  return (q.cross(p)).normalized();
}
