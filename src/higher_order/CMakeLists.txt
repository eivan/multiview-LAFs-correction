if (OpenMVG_USE_INTERNAL_EIGEN)
  set(EIGEN_INCLUDE_DIR_HINTS ${PROJECT_SOURCE_DIR}/dependencies/openMVG/third_party/eigen)
endif (OpenMVG_USE_INTERNAL_EIGEN)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

# Set header and source files to be added to the library
file(GLOB_RECURSE HDRSSRCS_LIB 
  "*.hpp"
  "*.cpp"
)

add_library(higher_order STATIC
  ${HDRSSRCS_LIB}
)

target_link_libraries(higher_order
  Eigen3::Eigen
)