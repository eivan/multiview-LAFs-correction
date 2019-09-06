# Optimal Multi-view Correction of Local Affine Frames
Code for our BMVC2019 paper: "[Optimal Multi-view Correction of Local Affine Frames](https://bmvc2019.org/wp-content/uploads/papers/0816-paper.pdf)" by Ivan Eichhardt and Daniel Barath. Further available resources: [paper](https://bmvc2019.org/wp-content/uploads/papers/0816-paper.pdf), [supplementary material](https://bmvc2019.org/wp-content/uploads/papers/0816-supplementary.pdf).

Cite it as
```
@InProceedings{Eichhardt_Barath_2019_BMVC,
	author = {Eichhardt, Ivan and Barath, Daniel},
	title = {Optimal Multi-view Correction of Local Affine Frames},
	booktitle = {British Machine Vision Conference (BMVC)},
	month = {September},
	year = {2019}
}
```

Introduction
------------

A method is proposed for correcting the parameters of a sequence of detected local affine frames through multiple views. The technique requires the epipolar geometry to be pre-estimated between each image pair. It exploits the constraints which the camera movement implies, in order to apply a closed-form correction to the parameters of the input affinities. Also, it is shown that the rotations and scales obtained by partially affine-covariant detectors, e.g. AKAZE or SIFT, can be upgraded to be full affine frames by the proposed algorithm. It is validated both in synthetic experiments and on publicly available real-world datasets that the method almost always improves the output of the evaluated affine-covariant feature detectors. As a by-product, these detectors are compared and the ones obtaining the most accurate affine frames are reported. To demonstrate the applicability in real-world scenarios, we show that the proposed technique improves the accuracy of pose estimation for a camera rig, surface normal and homography estimation.

Building
--------

See [BUILD](https://github.com/eivan/multiview-LAFs-correction/blob/master/BUILD.md)

Dependencies:

- Eigen3
- OpenMVG ([a modified version](https://github.com/eivan/openMVG/tree/develop)) (optional, for specific samples)

OpenMVG examples (optional)
---------------------------

Follow the Wiki of the [OpenMVG examples](https://github.com/eivan/multiview-LAFs-correction/wiki/Running-examples-that-use-OpenMVG) for a description of the OpenMVG-based tools provided with this repository.


