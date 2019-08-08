#pragma once

#include "types.hpp"

namespace higher_order::stereo {

  // TODO: reference
  Eigen::Vector3d EstimateNormal_FNE(
    const Mat23& gradient1, const Mat23& gradient2,
    const Eigen::Matrix2d& affine);

} // higher_order::stereo