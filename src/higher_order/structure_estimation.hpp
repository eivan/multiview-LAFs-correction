#pragma once

#include "types.hpp"

namespace higher_order::stereo {

/*!
 * Performs simple surface normal estimation using Fast Normal Estimation (FNE).
 * See: D. Barath, J. Molnar and L. Hajder, Optimal surface normal from affine transformation, Int. Conf. on Computer Vision Theory and Applications. (2015)
 * \param gradient1   the gradient of the world-to-camera (3D->2D) projection of the first view
 * \param gradient2   the gradient of the world-to-camera (3D->2D) projection of the second view
 * \param affine      the linear transformation part of an Affine Correspondences
 * \return            estimated 3D surface normal
 */
  Vec3d EstimateNormal_FNE(
    const Mat23d& gradient1, const Mat23d& gradient2,
    const Mat2d& affine);

} // higher_order::stereo