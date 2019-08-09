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
    openMVG_geometry
    openMVG_multiview
    openMVG_numeric
    openMVG_sfm
    openMVG_system
  )
else (USE_INTERNAL_OPENMVG)
  target_link_libraries(Sample_OpenMVG_RefineLAFs
    higher_order
    OpenMVG::openMVG_camera
    OpenMVG::openMVG_geometry
    OpenMVG::openMVG_multiview
    OpenMVG::openMVG_numeric
    OpenMVG::openMVG_sfm
    OpenMVG::openMVG_system
  )
endif (USE_INTERNAL_OPENMVG)

add_dependencies(Sample_OpenMVG_RefineLAFs
	higher_order
)