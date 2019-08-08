file(GLOB_RECURSE HDRSSRCS_SAMPLE 
  "*.hpp"
  "*.cpp"
)

add_executable(Sample_Synthetic
  ${HDRSSRCS_SAMPLE}
)

target_link_libraries(Sample_Synthetic
  higher_order
)

add_dependencies(Sample_Synthetic
	higher_order
)