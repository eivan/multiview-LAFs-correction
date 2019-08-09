#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <unordered_map>

#include <unsupported/Eigen/MatrixFunctions>

#include <third_party/cmdLine/cmdLine.h>
#include <third_party/progress/progress_display.hpp>
#include <third_party/stlplus3/filesystemSimplified/file_system.hpp>
#include <third_party/vectorGraphics/svgDrawer.hpp>

#include <openMVG/types.hpp>
#include <openMVG/sfm/sfm_data.hpp>
#include <openMVG/sfm/sfm_data_io.hpp>

using namespace openMVG;
using namespace openMVG::sfm;
using namespace openMVG::geometry;
using namespace std;

void DrawSVGs(const SfM_Data& sfm_data, const std::string& sOutDir);

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
      << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  // Read the input SfM scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \"" << sSfM_Data_Filename
      << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir)) {
    stlplus::folder_create(sOutDir);
  }

  DrawSVGs(sfm_data, sOutDir);

  return 0;
}

void DrawSVGs(const SfM_Data& sfm_data, const std::string& sOutDir) {

  // Generate random colors to draw regions with.
  // Matching regions are sharing colors.
  std::unordered_map<IndexT, std::string> randhexs;
  for (auto& landmark : sfm_data.structure) {
    const char* hex_digits = "0123456789ABCDEF";
    std::string str = "#FFFFFF";
    for (int i = 1; i < 7; ++i) {
      str[i] = hex_digits[rand() % 16];
    }
    randhexs[landmark.first] = str;
  }

  C_Progress_display my_progress_bar(sfm_data.GetViews().size(), std::cout,
    "\n- PROCESSING VIEWS (Writing SVGs) -\n");

  // Draw SVGs for all views.
  for (auto& v : sfm_data.views) {
    auto view = v.second.get();

    const std::string view_filename =
      stlplus::create_filespec(sfm_data.s_root_path, view->s_Img_path);

    svg::svgDrawer svgStream(view->ui_width, view->ui_height);
    svgStream.drawImage(view_filename, view->ui_width, view->ui_height);

    for (auto& lm : sfm_data.structure) {
      auto it = lm.second.obs.find(view->id_view);
      if (it != lm.second.obs.end()) {
        const auto& M = it->second.M;
        const auto& x = it->second.x;

        auto svd = M.jacobiSvd(Eigen::ComputeFullU);
        double rx = (svd.singularValues().x()), ry = (svd.singularValues().y());

        openMVG::Mat2 U = svd.matrixU();
        svgStream.drawAffine(x.x(), x.y(), rx, ry, U(0, 0), U(0, 1), U(1, 0),
          U(1, 1), svg::svgStyle().stroke(randhexs.at(lm.first), 2.0), true);

        svgStream.drawLine(x.x(), x.y(), x.x() + M(0, 0), x.y() + M(1, 0),
          svg::svgStyle().stroke("red"));
        svgStream.drawLine(x.x(), x.y(), x.x() - M(0, 1), x.y() - M(1, 1),
          svg::svgStyle().stroke("blue"));
      }

    }

    // Save the SVG file
    std::string output_filename = stlplus::create_filespec(
      sOutDir, "view_" + std::to_string(view->id_view), ".svg");
    
    std::ofstream svgFile(output_filename.c_str());
    if (svgFile.is_open()) {
      svgFile << svgStream.closeSvgFile().str();
      svgFile.close();
    }

    ++my_progress_bar;
  }

}
