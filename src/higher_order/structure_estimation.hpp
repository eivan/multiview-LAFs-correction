#pragma once

#include "types.hpp"

namespace higher_order::stereo {

  // TODO: reference
  Vec3d EstimateNormal_FNE(
    const Mat23d& gradient1, const Mat23d& gradient2,
    const Mat2d& affine);

} // higher_order::stereo