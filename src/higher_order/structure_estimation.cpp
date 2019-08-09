#include "structure_estimation.hpp"

higher_order::Vec3d higher_order::stereo::EstimateNormal_FNE(
  const Mat23d& gradient1, const Mat23d& gradient2,
  const Mat2d& affine) {

  const  Vec3d p = affine(1, 1) * gradient1.row(1).cross(gradient2.row(0)) -
    affine(0, 0) * gradient2.row(1).cross(gradient1.row(0));
  const  Vec3d q = affine(1, 0) * gradient2.row(0).cross(gradient1.row(0)) -
    affine(0, 1) * gradient1.row(1).cross(gradient2.row(1));
  return (q.cross(p)).normalized();
}
