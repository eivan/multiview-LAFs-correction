// MIT License
// 
// Copyright (c) 2019 Ivan Eichhardt
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
// Please contact the author of this library if you have any questions.
// Author: Ivan Eichhardt (ivan.eichhardt@sztaki.mta.hu)
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