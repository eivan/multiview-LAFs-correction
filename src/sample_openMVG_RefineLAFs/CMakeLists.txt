file(GLOB_RECURSE HDRSSRCS_SAMPLE 
  "*.hpp"
  "*.cpp"
)

add_executable(Sample_OpenMVG_RefineLAFs
  ${HDRSSRCS_SAMPLE}
)

if (USE_INTERNAL_OPENMVG)
  target_link_libraries(Sample_OpenMVG_RefineLAFs
    higher_order
    openMVG_camera
    openMVG_exif
    openMVG_features
    openMVG_geodesy
    openMVG_geometry
    openMVG_graph
    openMVG_image
    openMVG_linearProgramming
    openMVG_matching
    openMVG_matching_image_collection
    openMVG_multiview
    openMVG_numeric
    openMVG_robust_estimation
    openMVG_sfm
    openMVG_system
  )
else (USE_INTERNAL_OPENMVG)
  target_link_libraries(Sample_OpenMVG_RefineLAFs
    higher_order
    OpenMVG::openMVG_camera
    OpenMVG::openMVG_exif
    OpenMVG::openMVG_features
    OpenMVG::openMVG_geodesy
    OpenMVG::openMVG_geometry
    OpenMVG::openMVG_graph
    OpenMVG::openMVG_image
    OpenMVG::openMVG_linearProgramming
    OpenMVG::openMVG_matching
    OpenMVG::openMVG_matching_image_collection
    OpenMVG::openMVG_multiview
    OpenMVG::openMVG_numeric
    OpenMVG::openMVG_robust_estimation
    OpenMVG::openMVG_sfm
    OpenMVG::openMVG_system
  )
endif (USE_INTERNAL_OPENMVG)

add_dependencies(Sample_OpenMVG_RefineLAFs
	higher_order
)