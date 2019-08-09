if (USE_INTERNAL_OPENMVG)
  set(EIGEN_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/dependencies/openMVG/src/third_party/eigen)
else (USE_INTERNAL_OPENMVG)
  find_package (Eigen3 3.3 REQUIRED NO_MODULE)
  set(EIGEN_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
endif (USE_INTERNAL_OPENMVG)

# Set header and source files to be added to the library
file(GLOB_RECURSE HDRSSRCS_LIB 
  "*.hpp"
  "*.cpp"
)

add_library(higher_order STATIC
  ${HDRSSRCS_LIB}
)

#target_link_libraries(higher_order
#  Eigen3::Eigen
#)
target_include_directories(higher_order
  PUBLIC
    ${EIGEN_INCLUDE_DIRS}
)