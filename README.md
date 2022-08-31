# DDM-UQ
Software package for differential dynamic microscopy with uncertainty quantification
MATLAB version 0.5.3

Description:
This package is the application of DDM-UQ method (uncertainty quantification in differential dynamic microscopy). It allows 2-D particle movement simulation; and fully automated and high-throughput estimation and analysis of image structure function, intermediate scattering function and mean squared displacement (MSD) using Gaussian process regression. 

Reference:
Gu, M., Luo, Y., He, Y., Helgeson, M. E., & Valentine, M. T. (2021). Uncertainty quantification and estimation in differential dynamic microscopy. Physical Review E, just accepted, arXiv preprint arXiv:2105.01200.

Installation:
To use this package, please make sure you have installed the optimization toolbox in MATLAB.

Contents:
DDM-UQ.html - documentation file with description, result of some examples and syntax. 
example.m - demonstration examples.
simulation.m - simulation module, allow user simulates 2-D particle movement.
processing.m - processing module, process Fourier transformation for intensity profiles.
analysis.m - analysis module, robustly estimate the image structure function, mean squared displacement and other quantities of interest. 
...(other functions: kernel and log likelihood)

Authors:
Yue He, Mengyang Gu

Department of Statistics and Applied Probability 
University of California, Santa Barbara

Email: mengyang@pstat.ucsb.edu

This research is supported by NSF BioPACIFIC (DMR-1933487) and DMS-2053423. The authors thank the contributions by Yimin Luo, Matthew Helgeson and Megan Valentine from UCSB. 
