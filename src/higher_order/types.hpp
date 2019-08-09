#pragma once

#include <Eigen/Dense>

namespace higher_order {

  using ViewID = uint32_t;

  using Vec2d = Eigen::Matrix<double, 2, 1>;
  using Vec3d = Eigen::Matrix<double, 3, 1>;

  using Mat2d = Eigen::Matrix<double, 2, 2>;
  using Mat3d = Eigen::Matrix<double, 3, 3>;

  using Mat23d = Eigen::Matrix<double, 2, 3>;
  using Mat32d = Eigen::Matrix<double, 3, 2>;
  using MatX2d = Eigen::Matrix<double, Eigen::Dynamic, 2>;

  struct EpipolarConstraintOnLAFs {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Vec2d a;
    Vec2d b;

    std::pair < ViewID, ViewID > view_pair;
  };

  struct ObservedLAF {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Mat2d M;
    ViewID viewID;
  };

  using EpipolarConstraints = std::vector<EpipolarConstraintOnLAFs, Eigen::aligned_allocator<EpipolarConstraintOnLAFs>>;

  using ObservedLAFs = std::vector<ObservedLAF, Eigen::aligned_allocator<ObservedLAF>>;
}