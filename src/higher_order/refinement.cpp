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
#include "refinement.hpp"

#include <vector>
#include <unordered_map>

void higher_order::stereo::CreateFirstOrderConstraint(
  const Mat3d& F,
  const Vec2d& p1, const Vec2d& p2,
  Vec2d& a, Vec2d& b)
{
  // See Eq. (2) and (3) of [1]
  a = (F * p1.homogeneous()).head<2>();
  b = (F.transpose() * p2.homogeneous()).head<2>();
}

void higher_order::stereo::CreateFirstOrderConstraint(
  const Mat3d& E,
  const Vec3d& q1, const Vec3d& q2,
  const Mat32d& nabla_q1, const Mat32d& nabla_q2,
  Vec2d& a, Vec2d& b)
{
  // See Eq. (2) and (3) of [1]
  a = nabla_q2.transpose() * E * q1;
  b = nabla_q1.transpose() * E.transpose() * q2;
}

void higher_order::stereo::RefineAffineCorrespondence(
  Mat2d& A,
  const Vec2d& a, const Vec2d& b)
{
  A -= (1 / a.dot(a)) * a * (a.transpose() * A + b.transpose());
}