# Optimal Multi-view Correction of Local Affine Frames

Build instructions
------------------

Required tools:

- CMake
- Git
- C/C++ compiler

Dependencies:

- Eigen3

Optional dependencies:

- OpenMVG ([a modified version](/eivan/openMVG/tree/develop))

Note:

- CMAKE variables you can configure:
<a name="cmakevariables"></a>

  - USE_INTERNAL_OPENMVG (ON(default)/OFF)
      - Build OpenMVG as a submodule
  - BUILD_OPENMVG_EXAMPLES (ON(default)/OFF)
      - Build sample applications that use openMVG
	  
Checking out the project and build it
--------------------------------------

- [Getting the project](#checkout)
- [Compiling](#compiling)

Getting the project
--------------------
<a name="checkout"></a>

Getting the sources (and the submodules):
```shell
$ git clone --recursive https://github.com/eivan/multiview-LAFs-correction.git
```
or
```shell
$ git clone https://github.com/eivan/multiview-LAFs-correction.git
$ cd multiview-LAFs-correction
$ git submodule init
$ git submodule update
```

Note that if you do not intend to build the samples using OpenMVG, you do not need GIT submodules and recursive checkout.
```shell
$ git clone https://github.com/eivan/multiview-LAFs-correction.git
```

Compiling
-------------------
<a name="compiling"></a>

1. Make a directory for the build files to be generated.
```shell
$ mkdir build_dir
$ cd build_dir
```

2. Configure CMAKE.
```shell
$ cmake-gui ../src
```
See the choice of [CMAKE options](#cmakevariables).

3. Compile.