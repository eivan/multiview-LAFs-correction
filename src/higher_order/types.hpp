#pragma once

#include <Eigen/Dense>

namespace higher_order {
	
  using ViewID = size_t;

  using Mat23 = Eigen::Matrix<double, 2, 3>;
  using Mat32 = Eigen::Matrix<double, 3, 2>;
  using MatX2 = Eigen::Matrix<double, Eigen::Dynamic, 2>;

  struct EpipolarConstraintOnLAFs {
    std::pair < ViewID, ViewID > view_pair;
    Eigen::Vector2d a;
    Eigen::Vector2d b;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };

  using EpipolarConstraints = std::vector<EpipolarConstraintOnLAFs>;

  using ObservedLAFs = std::vector<std::pair<ViewID, Eigen::Matrix2d>>;
}