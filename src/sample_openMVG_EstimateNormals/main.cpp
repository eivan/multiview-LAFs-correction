#include <iostream>
#include <unordered_map>

#include <Eigen/Dense>
#include <third_party/cmdLine/cmdLine.h>
#include <third_party/progress/progress_display.hpp>
#include <third_party/stlplus3/filesystemSimplified/file_system.hpp>

#include <openMVG/sfm/sfm_data.hpp>
#include <openMVG/sfm/sfm_data_io.hpp>

#include <higher_order/structure_estimation.hpp>

using namespace openMVG;
using namespace openMVG::sfm;
using namespace openMVG::geometry;

void EstimateNormalsForLandmark(
  const SfM_Data& sfm_data,
  Landmark& landmark,
  bool use_observed_region_centers = false);

void FlipNormals(SfM_Data& sfm_data);

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

  // Perform surface normal estimation
  C_Progress_display my_progress_bar(sfm_data.GetLandmarks().size(), std::cout,
    "\n- PROCESSING LANDMARKS (Estimating surface normals) -\n");
  for (auto& landmark : sfm_data.structure) {
    auto& obs = landmark.second.obs;

    EstimateNormalsForLandmark(sfm_data, landmark.second);

    ++my_progress_bar;
  }

  // Make the surface normals orient towards viewpoints
  FlipNormals(sfm_data);

  // Save the output scene as PLY file
  std::string output_Filename = stlplus::create_filespec(
    sOutDir, stlplus::basename_part(sSfM_Data_Filename), "ply");
  if (!Save(sfm_data, output_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "Cannot write to output PLY file \"" << output_Filename
      << "\"." << std::endl;
    return EXIT_FAILURE;
  }
}

void EstimateNormalsForLandmark(
  const SfM_Data& sfm_data,
  Landmark& landmark,
  bool use_observed_region_centers)
{
  using namespace higher_order;

  auto& obs = landmark.obs;

  std::vector<Mat23> gradients(obs.size());
  std::vector<Mat2> Ms(obs.size());

  // populate LAFs and pre-evaluate bearing vectors and gradients, etc
  auto index = 0;
  for (const auto& observation : obs) {
    auto viewID = observation.first;
    auto& ob = observation.second;
    const Mat2 M = ob.M;

    const auto view = sfm_data.views.at(viewID).get();
    const auto intrinsic = sfm_data.intrinsics.at(view->id_intrinsic).get();
    const auto& pose = sfm_data.GetPoseOrDie(view);

    const Vec3 Y = use_observed_region_centers
      // case1: without scale, using the bearing vector
      ? intrinsic->operator()(intrinsic->get_ud_pixel(ob.x))
      // case2: using correct and accurate scale (known depth)
      : pose(landmark.X);

    gradients[index] =
      intrinsic->project_gradient(Y, false) *
      pose.rotation();

    Ms[index] = ob.M;

    ++index;
  }

  // Estimate normal using FNE
  landmark.N = stereo::EstimateNormal_FNE(
    gradients[0], gradients[1], 
    Ms[1] * Ms[0].inverse());
}

void FlipNormals(SfM_Data& sfm_data) {
  for (auto& landmark : sfm_data.structure) {
    int num_inverted = 0;
    for (const auto& ob : landmark.second.obs) {
      const auto& view = sfm_data.GetViews().at(ob.first).get();
      const auto& pose = sfm_data.GetPoses().at(view->id_pose);
      const auto& intrinsic =
        sfm_data.GetIntrinsics().at(view->id_intrinsic).get();

      Vec3 v = pose.rotation().transpose() * intrinsic->operator()(ob.second.x);
      if (landmark.second.N.dot(v) > 0) ++num_inverted;
    }

    if (num_inverted > landmark.second.obs.size() / 2) landmark.second.N *= -1;
  }
}