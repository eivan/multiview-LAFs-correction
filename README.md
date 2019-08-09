# Optimal Multi-view Correction of Local Affine Frames
Code for our BMVC2019 paper: "Optimal Multi-view Correction of Local Affine Frames" by Ivan Eichhardt and Daniel Barath

Preprint available at [arXiv.org](https://arxiv.org/abs/1905.00519).

Cite it as
```
% TODO
```

Introduction
------------

A method is proposed for correcting the parameters of a sequence of detected local affine frames through multiple views. The technique requires the epipolar geometry to be pre-estimated between each image pair. It exploits the constraints which the camera movement implies, in order to apply a closed-form correction to the parameters of the input affinities. Also, it is shown that the rotations and scales obtained by partially affine-covariant detectors, e.g. AKAZE or SIFT, can be upgraded to be full affine frames by the proposed algorithm. It is validated both in synthetic experiments and on publicly available real-world datasets that the method almost always improves the output of the evaluated affine-covariant feature detectors. As a by-product, these detectors are compared and the ones obtaining the most accurate affine frames are reported. To demonstrate the applicability in real-world scenarios, we show that the proposed technique improves the accuracy of pose estimation for a camera rig, surface normal and homography estimation.

Building
--------

See [BUILD](https://github.com/eivan/multiview-LAFs-correction/blob/master/BUILD.md)

Dependencies:

- Eigen3
- OpenMVG ([a modified version](/eivan/openMVG/tree/develop)) (optional, for specific samples)