#include <iostream>
#include <unordered_map>

#include <Eigen/Dense>
#include <third_party/cmdLine/cmdLine.h>
#include <third_party/progress/progress_display.hpp>
#include <third_party/stlplus3/filesystemSimplified/file_system.hpp>

#include <openMVG/multiview/essential.hpp>
#include <openMVG/sfm/sfm_data.hpp>
#include <openMVG/sfm/sfm_data_io.hpp>

#include <higher_order/refinement.hpp>

using namespace openMVG;
using namespace openMVG::sfm;
using namespace openMVG::geometry;

void RefineLAFsForLandmark(
  const SfM_Data& sfm_data,
  Landmark& landmark,
  bool use_observed_region_centers = false);

void RefineAllLAFs(SfM_Data& sfm_data);

int main(int argc, char** argv) {
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sOutDir = "";

  cmd.add(make_option('i', sSfM_Data_Filename, "sfmdata"));
  cmd.add(make_option('o', sOutDir, "outdir"));

  try {
    if (argc == 1) {
      throw std::string("Invalid command line parameter.");
    }
    cmd.process(argc, argv);
  }
  catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata] filename, the SfM_Data file\n"
      << "[-o|--outdir] path where refined SfM_Data with refined LAFs\n"
      << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir)) {
    stlplus::folder_create(sOutDir);
  }

  // Read the input SfM scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \"" << sSfM_Data_Filename
      << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  // Perform refinement for all landmarks
  RefineAllLAFs(sfm_data);

  // Save the output scene with the refined observations
  std::string output_sSfM_Data_Filename = stlplus::create_filespec(
    sOutDir, "refined_" + stlplus::basename_part(sSfM_Data_Filename), "bin");
  if (!Save(sfm_data, output_sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "Cannot write to output SfM_Data file \"" << output_sSfM_Data_Filename
      << "\"." << std::endl;
    return EXIT_FAILURE;
  }
}

void RefineLAFsForLandmark(
  const SfM_Data& sfm_data,
  Landmark& landmark,
  bool use_observed_region_centers) {
  using namespace higher_order;

  auto& obs = landmark.obs;

  EpipolarConstraints constraints;
  constraints.reserve((obs.size() * (obs.size() - 1)) / 2);

  ObservedLAFs LAFs(obs.size());

  struct CachedView {
    openMVG::Vec3 bearing;
    openMVG::Mat32 gradient;
    openMVG::Mat2 distortion_gradient;

    inline Mat2 distort(const Mat2& M) const {
      return 
         distortion_gradient * M;
    }

    inline Mat2 undistort(const Mat2& M) const {
      return 
         distortion_gradient.inverse() * M;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };

  std::unordered_map<ViewID, CachedView> cached_views;

  // populate LAFs and pre-evaluate bearing vectors and gradients, etc
  auto index = 0;
  for (const auto& observation : obs) {
    auto viewID = observation.first;
    auto& ob = observation.second;
    const Mat2 M = ob.M;

    const auto view = sfm_data.views.at(viewID).get();
    const auto intrinsic = sfm_data.intrinsics.at(view->id_intrinsic).get();
    const auto& pose = sfm_data.GetPoseOrDie(view);

    const Vec2 x = use_observed_region_centers
      ? intrinsic->get_ud_pixel(ob.x)
      : intrinsic->project(pose(landmark.X), true);

    cached_views[viewID] = {
      /*bearing*/ intrinsic->operator()(x),
      /*gradient*/ intrinsic->gradient(x),
      intrinsic->get_d_pixel_gradient(x)
    };

    LAFs[index].viewID = viewID;
    LAFs[index].M = cached_views[viewID].undistort(M);
    ++index;
  }

  for (auto it = obs.begin(); it != obs.end(); ++it) {
    const auto& viewID_i = it->first;
    const auto view_i = sfm_data.views.at(viewID_i).get();
    const auto& pose_i = sfm_data.GetPoseOrDie(view_i);

    // bearing vector and gradient corresponding to viewID_i
    const auto& [q_i, q_grad_i, _3] = cached_views[viewID_i];

    auto tmp = it;
    std::advance(tmp, 1);
    for (auto jt = tmp; jt != obs.end(); ++jt) {
      const auto& viewID_j = jt->first;
      const auto view_j = sfm_data.views.at(viewID_j).get();
      const auto& pose_j = sfm_data.poses.at(view_j->id_pose);

      // bearing vector and gradient corresponding to viewID_j
      const auto& [q_j, q_grad_j, _3] = cached_views[viewID_j];

      Mat3 E;
      // TODO: this can be made a bit more efficient.
      EssentialFromRt(pose_i.rotation(), pose_i.translation(),
        pose_j.rotation(), pose_j.translation(), &E);

      EpipolarConstraintOnLAFs constraint;
      constraint.view_pair = { viewID_i, viewID_j };

      stereo::CreateFirstOrderConstraint(E,
        q_i, q_j, q_grad_i, q_grad_j,
        constraint.a, constraint.b);

      constraints.emplace_back(constraint);
    }
  }

  multiview::RefineLocalAffineFrames(
    constraints,
    LAFs
  );

  for (const auto& laf : LAFs) {
    auto& viewID = laf.viewID;
    auto& M = laf.M;
    obs[viewID].M = cached_views[viewID].distort(M);
  }
}

void RefineAllLAFs(SfM_Data& sfm_data) {
  C_Progress_display my_progress_bar(sfm_data.GetLandmarks().size(), std::cout,
    "\n- PROCESSING LANDMARKS (Refining Local Affine Frames) -\n");

  for (auto& landmark : sfm_data.structure) {
    auto& obs = landmark.second.obs;

    RefineLAFsForLandmark(sfm_data, landmark.second);

    ++my_progress_bar;
  }
}