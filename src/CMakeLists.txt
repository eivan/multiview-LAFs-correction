cmake_minimum_required(VERSION 3.8)

project(multiview-LAFs-correction LANGUAGES CXX)

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

add_subdirectory(higher_order)

# ==============================================================================
# Structure: Applications
# ==============================================================================

add_subdirectory(sample_synthetic)

if (BUILD_OPENMVG_EXAMPLES)
  add_subdirectory(sample_openMVG_LAFsToSVG)
  add_subdirectory(sample_openMVG_RefineLAFs)
  add_subdirectory(sample_openMVG_EstimateNormals)
endif (BUILD_OPENMVG_EXAMPLES)