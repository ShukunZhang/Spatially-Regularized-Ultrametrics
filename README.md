# Spatially-Regularized-Ultrametrics

This repository contains the Matlab code for the paper Hyperspectral Image Clustering withSpatially-Regularized Ultrametrics

To use our implementation of SRUSC, please cite our paper at: https://arxiv.org/abs/2004.05048 

This implementation is partially adapted from the code https://jmurphy.math.tufts.edu/Code/LLPD_SpectralClustering_V2.1.zip. Please also cite the paper: Little, A., Maggioni, M., Murphy, J.M. Path-Based Spectral Clustering: Guarantees, Robustness to Outliers, and Fast Algorithms. Journal of Machine Learning Research, 21(6), pp. 1-66. 2020. https://www.jmlr.org/papers/volume21/18-085/18-085.pdf

Source of the HSI data set: http://www.ehu.eus/ccwintco/index.php/Hyperspectral_Remote_Sensing_Scenes#Salinas_scene

If you only need our method, please refer to the file in the folder SRUSC/scripts. The two synthetic HSI are FourSpheres and ThreeCubes, and two real HSI are SalinasA and PaviaU.

If you need all the comparisions, please use the older in the folder SRUSC/RunAll. Note: To run all the comparisions, please download the code from the following website:

Diffusion Learning (DL):https://jmurphy.math.tufts.edu/Code/ J.M. Murphy and M. Maggioni. Unsupervised clustering and active learning of hyperspectral images with nonlinear diffusion.
IEEE Transactions on Geoscience and Remote Sensing, 57(3):1829–1845, 2019. M. Maggioni and J.M. Murphy. Learning by unsupervised nonlinear diffusion. Journal of Machine Learning Research, 20(160):1–56, 2019.

Clustering by Fast Search and Find of Density Peaks (FSFDPC): https://people.sissa.it/~laio/Research/Res_clustering.php A. Rodriguez and A. Laio. Clustering by fast search and find of densitypeaks.Science, 344(6191):1492–1496, 2014.

Nonnegative Matrix Factorization (NMF): https://sites.google.com/site/nicolasgillis/code N. Gillis, D. Kuang, and H. Park. Hierarchical clustering of hyperspectral images using rank-two nonnegative matrix factorization. IEEE Transactions on Geoscience and Remote Sensing, 53(4):2066–2078, 2015.

Local Covariance Matrix Representation: (LCMR): https://github.com/henanjun/LCMR L. Fang, N. He, S. Li, A.J. Plaza, and Javier J. Plaza. A new spatial–spectral feature extraction method for hyperspectral images using local covariance matrix representation. IEEE Transactions on Geoscience and Remote Sensing, 56(6):3534–3546, 2018.
