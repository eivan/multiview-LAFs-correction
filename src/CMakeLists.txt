cmake_minimum_required(VERSION 3.8)

project(multiview-LAFs-correction LANGUAGES CXX)

# indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

# indicate if OPENMP should be enabled
option(USE_OPENMP "Use OPENMP" ON)

option(USE_INTERNAL_OPENMVG "Use eivan/OpenMVG" ON)
option(BUILD_OPENMVG_EXAMPLES "Build OpenMVG examples" ON)

# ==============================================================================
# Check C++17 support
# ==============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ==============================================================================
# OpenMP
# ==============================================================================
find_package(OpenMP)
if (USE_OPENMP)
	if(NOT OPENMP_FOUND)
		message(FATAL_
		
		"OPENMP not found.")
	endif()
	add_definitions(-DUSE_OPENMP)
	set(TRGT_LNK_LBS_ADDITIONAL OpenMP::OpenMP_CXX)
endif (USE_OPENMP)

# ==============================================================================
# OpenMVG
# ==============================================================================

if (USE_INTERNAL_OPENMVG)
  set(OpenMVG_BUILD_DOC           OFF CACHE BOOL "Do not build the OpenMVG documentation")
  set(OpenMVG_BUILD_EXAMPLES      OFF CACHE BOOL "Do not build the OpenMVG samples applications")
  set(OpenMVG_BUILD_GUI_SOFTWARES OFF CACHE BOOL "Do not build the OpenMVG GUI softwares (QT5)")
  set(OpenMVG_USE_AS_SUBMODULE    ON CACHE BOOL "The internal eivan/OpenMVG should loose its install capabilities")
  include_directories(SYSTEM dependencies/openMVG/src)
  add_subdirectory(dependencies/openMVG/src)
else (USE_INTERNAL_OPENMVG)
  find_package(OpenMVG REQUIRED)
endif (USE_INTERNAL_OPENMVG)

# ==============================================================================
# Includes
# ==============================================================================
include_directories (
	${PROJECT_SOURCE_DIR}
)

# ==============================================================================
# Structure: Library
# ==============================================================================
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

add_subdirectory(higher_order)

# ==============================================================================
# Structure: Applications
# ==============================================================================
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

add_subdirectory(sample_synthetic)

if (BUILD_OPENMVG_EXAMPLES)
  add_subdirectory(sample_openMVG_LAFsToSVG)
  add_subdirectory(sample_openMVG_RefineLAFs)
  add_subdirectory(sample_openMVG_EstimateNormals)
endif (BUILD_OPENMVG_EXAMPLES)

# ==============================================================================
# Documentation
# ==============================================================================

# build doc if option is enabled
if (BUILD_DOC)
	# check if Doxygen is installed
	find_package(Doxygen)
	if (DOXYGEN_FOUND)
		# set input and output files
		set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
		set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

		# request to configure the file
		configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
		message("Doxygen build started")

		# note the option ALL which allows to build the docs together with the application
		add_custom_target( doc_doxygen ALL
			COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			COMMENT "Generating API documentation with Doxygen"
			VERBATIM )
	else (DOXYGEN_FOUND)
	  message("Doxygen is needed to be installed to generate the doxygen documentation")
	endif (DOXYGEN_FOUND)
endif (BUILD_DOC)