# P2C-Clustering
Clustering module in IEEE PAMI 2013: A Framework for Automatic Modeling from Point Cloud Data. A robust unsupervised clustering algorithm P2C, based on a hierarchical statistical analysis of the geometric properties of the data.

The entire clustering module in IEEE PAMI 2013: A Framework for Automatic Modeling from Point Cloud Data. 
A variant of the first stage of the algorithm was first introduced in IEEE CVPR 2009: Automatic reconstruction of cities from remote sensor data and can be easily reproduced by changing the distributions as described in the paper.

The copyright for the included bilateral filtering code belongs to Sylvain Paris and Frédo Durand.


References:
1. IEEE PAMI 2013: A Framework for Automatic Modeling from Point Cloud Data
2. IEEE CVPR 2009: Automatic reconstruction of cities from remote sensor data

More information about this work: www.poullis.org

Technical details:

- The project file is provided for Code::Blocks IDE. 
- It requires the libraries Image Magick, fftw3, and opencv2.
- A small sample file is provided in the bin folders [structured_pointcloud.pfm]
- Usage: P2C-Clustering structured_pointcloud.pfm
- Input files can be generated from LiDAR data using the StructurePointcloud code also provide in a separate repository.

A sample input file can be found in the data folder.

*IMPORTANT: To use this software, YOU MUST CITE the following in any resulting publication:*

@article{poullis2013framework,
  title={A Framework for Automatic Modeling from Point Cloud Data},
  author={Poullis, Charalambos},
  journal={Pattern Analysis and Machine Intelligence, IEEE Transactions on},
  volume={35},
  number={11},
  pages={2563--2575},
  year={2013},
  publisher={IEEE}
}

@inproceedings{poullis2009automatic,
  title={Automatic reconstruction of cities from remote sensor data},
  author={Poullis, Charalambos and You, Suya},
  booktitle={Computer Vision and Pattern Recognition, 2009. CVPR 2009. IEEE Conference on},
  pages={2775--2782},
  year={2009},
  organization={IEEE}
}

Disclaimer: This code is part of a larger software package. 
I have tried to extract and clean up only the code relating to the clustering however some unused code may still be present.
