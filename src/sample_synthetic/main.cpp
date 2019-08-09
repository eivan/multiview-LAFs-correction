#define _USE_MATH_DEFINES

#include <time.h>
#include <iomanip>
#include <iostream>
#include <set>

#include <higher_order/refinement.hpp>
#include <higher_order/structure_estimation.hpp>

#include "helper.hpp"


using namespace std;
using namespace higher_order;

void StereoTest_BMVC2016();
void StereoTest_BMVC2019();
void MultiviewTest_BMVC2019();

int main() {
  auto seed = time(0);

  srand(seed);
  cout
    << "================================================================================" << endl
    << "  Testing BMVC2016 (stereo only)" << endl
    << "================================================================================" << endl;
  StereoTest_BMVC2016();

  srand(seed);
  cout
    << "================================================================================" << endl
    << "  Testing BMVC2019 (stereo)" << endl
    << "================================================================================" << endl;
  StereoTest_BMVC2019();

  cout
    << "================================================================================" << endl
    << "  Testing BMVC2019 (multi-view)" << endl
    << "================================================================================" << endl;
  MultiviewTest_BMVC2019();

  return 0;
}

void CreateStereoPair(Mat3& K1, Mat3& K2, Mat3& R1, Mat3& R2, Vec3& t1, Vec3& t2) {
  // Create camera matrices
  K1 = (Mat3() << 1500, 0, 750, 0, 1500, 750, 0, 0, 1).finished();
  K2 = (Mat3() << 1500, 0, 750, 0, 1500, 750, 0, 0, 1).finished();

  // Generate poses (R, t)
  Vec3
    eye1{ 4.0, 10.0, 7.0 },
    eye2{ 6, 0, 9 },
    target{ 0,0,0 };

  R1 = LookAt_target(eye1, target);
  R2 = LookAt_target(eye2, target);

  t1 = -R1 * eye1;
  t2 = -R2 * eye2;
}

auto CreateMultiviewScene(size_t num_views) {
  assert(num_views > 1);

  struct View {
    ViewID index;
    Mat3 K, R;
    Vec3 t;

    inline Vec2 operator()(const Vec3& X) const {
      return project(K, R * X + t);
    }

    inline Mat23 gradient(const Vec3& X) const {
      return project_gradient(K, R * X + t) * R;
    }
  };

  std::vector<View> views(num_views);

  Vec3 target{ 0,0,0 };

  for (ViewID i = 0; i < num_views; ++i) {

    Vec3 eye = Vec3::Random().normalized() * 10;

    views[i].index = i;
    views[i].K = (Mat3() << 1500, 0, 750, 0, 1500, 750, 0, 0, 1).finished();
    views[i].R = LookAt_target(eye, target);
    views[i].t = -views[i].R * eye;
  }

  return views;
}

void StereoTest_BMVC2016() {

  Mat3 K1, K2, R1, R2;
  Vec3 t1, t2;
  CreateStereoPair(K1, K2, R1, R2, t1, t2);

  // Compute GT essential and fundamental matrices
  Mat3 E12;
  EssentialFromRt(R1, t1, R2, t2, &E12);
  Mat3 F12 = K2.inverse().transpose() * E12 * K1.inverse();

  // Generate a surflet (tangential plane) at 3D point X with normal N
  Vec3 X, N;
  X.setRandom();
  N.setRandom();
  N.normalize();
  Mat32 nulls = Nullspace<double>(N);

  // Project the surflet to get Local Affine Frames LAF(x1, M1) and LAF(x2, M2)
  Vec2
    x1 = project(K1, R1 * X + t1),
    x2 = project(K2, R2 * X + t2);
  Mat23
    x1_gradient = project_gradient(K1, R1 * X + t1) * R1,
    x2_gradient = project_gradient(K2, R2 * X + t2) * R2;
  Mat2
    M1 = x1_gradient * nulls,
    M2 = x2_gradient * nulls;

  // Compute the affine correspondence AC(x1, x2, A)
  Mat2 A_gt = M2 * M1.inverse();

  // Compute the first order constraint on an affine correspondence
  Vec2 a_12, b_12;
  stereo::CreateFirstOrderConstraint(F12, x1, x2, a_12, b_12);

  // add noise 
  ViewID view_1{ 0 }, view_2{ 1 };
  ObservedLAFs LAFs_Noisy{
    {view_1, M1},
    {view_2, M2}
  };
  for (auto& obs : LAFs_Noisy) {
    obs.second += Mat2::Random() * 1.0;
  }
  Mat2
    M1_noisy = LAFs_Noisy[0].second,
    M2_noisy = LAFs_Noisy[1].second;
  Mat2 A_noisy = M2_noisy * M1_noisy.inverse();

  // Perform refinement
  Mat2 A_refined = A_noisy;
  stereo::RefineAffineCorrespondence(A_refined, a_12, b_12);


  // Measure errors
  double err_gt = (A_gt.transpose() * a_12 + b_12).norm();
  double err_noisy = (A_noisy.transpose() * a_12 + b_12).norm();
  double err_refined = (A_refined.transpose() * a_12 + b_12).norm();

  cout << std::setprecision(7) << std::scientific;
  cout << "Magnitude of the epipolar constraint on (ground truth) A:\n\t" << err_gt << endl;
  cout << "Magnitude of the epipolar constraint on A_noisy:\n\t" << err_noisy << endl;
  cout << "Magnitude to the epipolar constraint after refinement:\n\t" << err_refined << endl;

  cout << endl;

  Vec3 N_noisy = stereo::EstimateNormal_FNE(x1_gradient, x2_gradient, A_noisy);
  Vec3 N_refined = stereo::EstimateNormal_FNE(x1_gradient, x2_gradient, A_refined);

  cout << std::setprecision(4) << std::fixed
    << "Normal estimation:" << endl
    << "  ground truth: \t" << N.transpose() << endl
    << "  estim. (noisy): \t" << N_noisy.transpose() << "\t angular error (deg):\t" <<
    (180 / M_PI) * acos(abs(N.dot(N_noisy))) << endl
    << "  estim. (refined): \t" << N_refined.transpose() << "\t angular error (deg):\t" <<
    (180 / M_PI) * acos(abs(N.dot(N_refined))) << endl;

  cout << endl;
}

void StereoTest_BMVC2019() {
  Mat3 K1, K2, R1, R2;
  Vec3 t1, t2;
  CreateStereoPair(K1, K2, R1, R2, t1, t2);

  // Compute GT essential and fundamental matrices
  Mat3 E12;
  EssentialFromRt(R1, t1, R2, t2, &E12);
  Mat3 F12 = K2.inverse().transpose() * E12 * K1.inverse();

  // Generate a surflet (tangential plane) at 3D point X with normal N
  Vec3 X, N;
  X.setRandom();
  N.setRandom();
  N.normalize();
  Mat32 nulls = Nullspace<double>(N);

  // Project the surflet to get Local Affine Frames LAF(x1, M1) and LAF(x2, M2)
  Vec2
    x1 = project(K1, R1 * X + t1),
    x2 = project(K2, R2 * X + t2);
  Mat23
    x1_gradient = project_gradient(K1, R1 * X + t1) * R1,
    x2_gradient = project_gradient(K2, R2 * X + t2) * R2;

  // Compute the first order constraint on an affine correspondence
  Vec2 a_12, b_12;
  stereo::CreateFirstOrderConstraint(F12, x1, x2, a_12, b_12);

  ViewID view_1{ 0 }, view_2{ 1 };

  ObservedLAFs LAFs_GT{
    {view_1, /* M1 = */ x1_gradient * nulls },
    {view_2, /* M2 = */ x2_gradient * nulls }
  };

  // Test
  ObservedLAFs LAFs_Noisy = LAFs_GT;

  // add noise
  for (auto& obs : LAFs_Noisy) {
    obs.second += Mat2::Random() * 1.0;
  }

  ObservedLAFs LAFs_Refined = LAFs_Noisy;

  multiview::RefineLocalAffineFrames(
    EpipolarConstraints{
      { {view_1, view_2}, a_12, b_12 }
    },
    LAFs_Refined
  );

  // Measure errors

  Mat2
    M1_gt = LAFs_GT[0].second,
    M2_gt = LAFs_GT[1].second,
    M1_noisy = LAFs_Noisy[0].second,
    M2_noisy = LAFs_Noisy[1].second,
    M1_refined = LAFs_Refined[0].second,
    M2_refined = LAFs_Refined[1].second;

  double err_gt = (M2_gt.transpose() * a_12 + M1_gt.transpose() * b_12).norm();
  double err_noisy = (M2_noisy.transpose() * a_12 + M1_noisy.transpose() * b_12).norm();
  double err_refined = (M2_refined.transpose() * a_12 + M1_refined.transpose() * b_12).norm();

  cout << std::setprecision(7) << std::scientific;
  cout << "Magnitude of the epipolar constraint on (ground truth) A:\n\t" << err_gt << endl;
  cout << "Magnitude of the epipolar constraint on A_noisy:\n\t" << err_noisy << endl;
  cout << "Magnitude to the epipolar constraint after refinement:\n\t" << err_refined << endl;
  cout << endl;

  Mat2 A_noisy = M2_noisy * M1_noisy.inverse();
  Vec3 N_noisy = stereo::EstimateNormal_FNE(x1_gradient, x2_gradient, A_noisy);

  Mat2 A_refined = M2_refined * M1_refined.inverse();
  Vec3 N_refined = stereo::EstimateNormal_FNE(x1_gradient, x2_gradient, A_refined);

  cout << std::setprecision(4) << std::fixed
    << "Normal estimation:" << endl
    << "  ground truth: \t" << N.transpose() << endl
    << "  estim. (noisy): \t" << N_noisy.transpose() << "\t angular error (deg):\t" << (180 / M_PI) * acos(abs(N.dot(N_noisy))) << endl
    << "  estim. (refined): \t" << N_refined.transpose() << "\t angular error (deg):\t" << (180 / M_PI) * acos(abs(N.dot(N_refined))) << endl;
  cout << endl;
}

void MultiviewTest_BMVC2019() {
  const auto views = CreateMultiviewScene(5);

  Vec3 X, N;
  X.setRandom();
  N.setRandom();
  N.normalize();
  Mat32 nulls = Nullspace<double>(N);

  // project surface onto image
  vector<Vec2> xs(views.size());
  ObservedLAFs Ms(views.size());

  for (auto& view : views) {
    xs[view.index] = view(X);
    Ms[view.index].first = view.index;
    Ms[view.index].second = view.gradient(X) * nulls;
  }

  // construct epipolar constraints
  EpipolarConstraints constraints;

  for (ViewID i = 0; i < views.size()-1; ++i) {
    for (ViewID j = i+1; j < views.size(); ++j) {
      // Compute GT essential and fundamental matrices
      Mat3 E12;
      EssentialFromRt(views[i].R, views[i].t, views[j].R, views[j].t, &E12);
      Mat3 F12 = views[j].K.inverse().transpose() * E12 * views[i].K.inverse();

      // Compute the first order constraint on an affine correspondence
      Vec2 a_12, b_12;
      stereo::CreateFirstOrderConstraint(F12, xs[i], xs[j], a_12, b_12);

      EpipolarConstraintOnLAFs constraint{ { i, j }, a_12, b_12 };
      constraints.emplace_back(constraint);
    }
  }

  // add noise
  ObservedLAFs Ms_noisy = Ms;
  for (auto& obs : Ms_noisy) {
    obs.second += Mat2::Random() * 1.0;
  }

  // refine
  ObservedLAFs Ms_Refined = Ms_noisy;
  multiview::RefineLocalAffineFrames(
    constraints,
    Ms_Refined
  );

  // measure errors
  Mat2
    M1_gt = Ms[0].second,
    M2_gt = Ms[1].second,
    M1_noisy = Ms_noisy[0].second,
    M2_noisy = Ms_noisy[1].second,
    M1_refined = Ms_Refined[0].second,
    M2_refined = Ms_Refined[1].second;
  Mat23
    x1_gradient = views[0].gradient(X),
    x2_gradient = views[1].gradient(X);

  Mat2 A_gt = M2_gt * M1_gt.inverse();

  Mat2 A_noisy = M2_noisy * M1_noisy.inverse();
  Vec3 N_noisy = stereo::EstimateNormal_FNE(x1_gradient, x2_gradient, A_noisy);

  Mat2 A_refined = M2_refined * M1_refined.inverse();
  Vec3 N_refined = stereo::EstimateNormal_FNE(x1_gradient, x2_gradient, A_refined);

  cout << std::setprecision(4) << std::fixed
    << "Normal estimation:" << endl
    << "  ground truth: \t" << N.transpose() << endl
    << "  estim. (noisy): \t" << N_noisy.transpose() << "\t angular error (deg):\t" << (180 / M_PI) * acos(abs(N.dot(N_noisy))) << endl
    << "  estim. (refined): \t" << N_refined.transpose() << "\t angular error (deg):\t" << (180 / M_PI) * acos(abs(N.dot(N_refined))) << endl;
  cout << endl;
}