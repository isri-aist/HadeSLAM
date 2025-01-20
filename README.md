# HadeSLAM = stella_vslam + Equirectangular RGB-D modality

Please check stella_vslam's [README.md](./README-stella_vslam.md) file for its specific matters on setting up, installing and licensing. 


## Overview

HadeSLAM is the Equirectangular RGB-D capability added to stella_vslam that is a monocular, stereo, and RGBD visual SLAM system.

### Features

The notable features are the same as those of stella_vslam with the notable addition to handling RGB and Depth equirectangular data together via the example code [run_tum_rgbd_slam.cc](./example/run_tum_rgbd_slam.cc) that as been extended a little for that purpose. 


### Acknowledgements

HadeSLAM leverages [StellaVSLAM](https://github.com/stella-cv/stella_vslam). 

### Examples

## Dataset

You can run the Equirectangular RGB-D visual SLAM with the [CD2F-JRL](https://extra.u-picardie.fr/nextcloud/index.php/s/6XZZtfokRamqR9N) benchmarking dataset.

## Community

Please see [GitHub Discussions](https://github.com/stella-cv/stella_vslam/discussions) for specific questions or notice any bugs about the StellaVSLAM. 

For specific matter related to the Equirectangular RGB-D SLAM, please open an issue. 

## Currently working on

- Rebase on latest StellaVSLAM fro asking a requesting a pull
- ROS integration
- Stereo mode
- Implementation of extra camera models

The higher up the list, the higher the priority.
Feedbacks and contributions are welcome!



## Citation of HadeSLAM


If HadeSLAM helps your research, please cite the paper for which we provide a BibTeX entry:

```
@article{hadeslam2025,
  author = {Caron, Guillaume and Sevillano, Grace and Paula, Matheus and Imamoglu, Nevrez},
  title = {{Hybrid Appearance-Depth 360-degree Simulatenous Localization And Mapping}},
  journal = {anonymous submission to the IEEE Trans. on Robotics, special section on Visual SLAM},
  year = {2025},
  pages = {1--10},
  numpages = {10}
}
```
