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
