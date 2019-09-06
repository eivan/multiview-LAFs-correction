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

#include <Eigen/Dense>

namespace higher_order {

  // an alias to indexing views
  using ViewID = uint32_t;

  // 2D vector with double elements
  using Vec2d = Eigen::Matrix<double, 2, 1>;
  // 3D vector with double elements
  using Vec3d = Eigen::Matrix<double, 3, 1>;
  // 2-by-2 matrix with double elements
  using Mat2d = Eigen::Matrix<double, 2, 2>;
  // 2-by-2 matrix with double elements
  using Mat3d = Eigen::Matrix<double, 3, 3>;
  // 2-by-2 matrix with double elements
  using Mat23d = Eigen::Matrix<double, 2, 3>;
  // 3-by-2 matrix with double elements
  using Mat32d = Eigen::Matrix<double, 3, 2>;
  // N-by-2 dynamic matrix with double elements
  using MatX2d = Eigen::Matrix<double, Eigen::Dynamic, 2>;

  /*! \class EpipolarConstraintOnLAFs
   * A struct that holds the vector-coefficients of a first-order epipolar constraint, and indices to the corresponding pair of view.
   */
  struct EpipolarConstraintOnLAFs {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /*! The first vector coefficient of the first-order epipolar constraint */
    Vec2d a;
    /*! The second vector coefficient of the first-order epipolar constraint */
    Vec2d b;

    /*! Indices to a pair of views */
    std::pair < ViewID, ViewID > view_pair;
  };

  /*! \class ObservedLAF
   * A struct that encapsulates the linear transformation part of a Local Affine Frame, and and index to the corresponding view.
   */
  struct ObservedLAF {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /*! The linear transformation part of a Local Affine Frame */
    Mat2d M;
    /*! The view-index of the corresponding Local Affine Frame */
    ViewID viewID;
  };

  // A vector containing first-order epipolar constraints (EpipolarConstraintOnLAFs)
  using EpipolarConstraints = std::vector<EpipolarConstraintOnLAFs, Eigen::aligned_allocator<EpipolarConstraintOnLAFs>>;

  // A vector constaining the linear transformation parts of Local Affine Frames (ObservedLAF)
  using ObservedLAFs = std::vector<ObservedLAF, Eigen::aligned_allocator<ObservedLAF>>;
}