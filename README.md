# Deconvolution and demixing of calcium imaging data

[![Join the chat at https://gitter.im/epnev/ca_source_extraction](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/epnev/ca_source_extraction?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

The code implements a method for simultaneous source extraction and spike inference from large scale calcium imaging movies. The code is suitable for the analysis of somatic imaging data. Implementation for the analysis of dendritic/axonal imaging data will be added in the future. 

The algorithm is presented in more detail in

Pnevmatikakis, E.A., Soudry, D., Gao, Y., Machado, T., Merel, J., ... & Paninski, L. (2016). Simultaneous denoising, deconvolution, and demixing of calcium imaging data. Neuron, in press, http://dx.doi.org/10.1016/j.neuron.2015.11.037

Pnevmatikakis, E.A., Gao, Y., Soudry, D., Pfau, D., Lacefield, C., ... & Paninski, L. (2014). A structured matrix factorization framework for large scale calcium imaging data analysis. arXiv preprint arXiv:1409.2903. http://arxiv.org/abs/1409.2903

Code description and related packages
=======

This repository contains a MATLAB implementation of the spatio-temporal demixing, i.e., (source extraction) code for large scale calcium imaging data. Related code can be found in the following links:

## Matlab 
- [Constrained deconvolution and source extraction with CNMF (this package)](https://github.com/epnev/ca_source_extraction)
- [MCMC spike inference](https://github.com/epnev/continuous_time_ca_sampler)
- [Group LASSO initialization and spatial CNMF](https://github.com/danielso/ROI_detect)

## Python
- [Constrained deconvolution for neural activity (spike) extraction](https://github.com/epnev/constrained_foopsi_python)
- [Source extraction with CNMF](https://github.com/agiovann/SOURCE_EXTRACTION_PYTHON)
- [Group LASSO initialization and spatial CNMF](https://github.com/danielso/ROI_detect)

## Integration with other libraries
- [SIMA](http://www.losonczylab.org/sima/1.3/): The [constrained deconvolution](https://github.com/losonczylab/sima/blob/master/sima/spikes.py) method has been integrated with SIMA, a Python based library for calcium imaging data analysis.
- [Thunder](http://thunder-project.org/): The [group LASSO initialization and spatial CNMF](https://github.com/j-friedrich/thunder/tree/LocalNMF) method has been integrated with Thunder, a library for large scale neural data analysis with Spark.


Usage and Documentation
=======
Check the demo scripts and documentation.pdf to get started.

Dependencies
========
The code uses the following libraries
- [NumPy](http://www.numpy.org/)
- [SciPy](http://www.scipy.org/)
- [Matplotlib](http://matplotlib.org/)
- [Scikit-Learn](http://scikit-learn.org/stable/)
- [Tifffile](https://pypi.python.org/pypi/tifffile) For reading tiff files. Other choices can work there too.

External Dependencies
============

The constrained deconvolution method (constrained_foopsi_python.py) can estimate with two different methods, each of which requires some additional packages:
1. 'spgl1'
  - [SPGL1](https://github.com/epnev/SPGL1_python_port) Python implementation of the SPGL1 method. Please use the "forked" repository linked here.
2. 'cvx'
  - [CVXOPT](http://cvxopt.org/) Required.
  - [PICOS](http://picos.zib.de/) Required.
  - [MOSEK](https://www.mosek.com/) Optional but strongly recommended for speed improvement, free for academic use.

In general 'spgl1' can be faster, but the python implementation is not as fast as in Matlab and not thoroughly tested.

Questions, comments, issues
=======
Please use the gitter chat room (use the button above) for questions and comments and create an issue for any bugs you might encounter.

License
=======

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


[//]: #  (function name                           | description )
[//]: #  (----------------------------------------|-----------------------------------)
[//]: #  (demo_script.m                           | wrapper code <br />)
[//]: #  (demoMovie.tif                           | Sample dataset for running the code (by W. Yang and D. Peterka) <br />)
[//]: #  (update_spatial_components.m             | update spatial components given temporal components and data <br />)
[//]: #  (update_temporal_components.m            | update temporal components given spatial components and data <br />)
[//]: #  (merge_ROIs.m                            | merge spatially overlapping components that are temporally correlated <br />)
[//]: #  (utilities/arpfit.m                      | estimation of noise level for every pixel and global time constants <br />)
[//]: #  (utilities/bigread2.m                    | read (parts of) large tiff stacks)
[//]: #  (utilities/com.m:                        | calculation of the center of mass of each component <br />)
[//]: #  (utilities/correlation_image.m           | calculates the correlation image of the movie <br />)
[//]: #  (utilities/extract_DF_F.m                | transforming the temporal components in the DF/F domain <br />)
[//]: #  (utilities/graph_connected_comp.m        | finds the connected components in a graph <br />)
[//]: #  (utilities/greedyROI2d.m                 | Greedy method for initializing the spatial and temporal components <br />)
[//]: #  (utilities/interp_missing_data.m         | Filling in missing data using linear interpolation <br />)
[//]: #  (utilities/lars_regression_noise.m       | solve a basis pursuit denoising problem using the LARS algorithm <br />)
[//]: #  (utilities/make_G_matrix.m               | construct a convolution/deconvolution matrix <br />)
[//]: #  (utilities/make_patch_video.m            | construct a video that displays the results of the algorithm <br />)
[//]: #  (utilities/order_ROIs.m                  | order found components based on their temporal activation and spatial size <br />)
[//]: #  (utilities/plain_foopsi.m                | projection of fluorescence onto the cone formed by the indicator dynamics )
[//]: #  (utilities/plot_contours.m               | contour plot of found components and creation of a json file <br />)
[//]: #  (utilities/tiff_reader.m                 | loading a tiff stack into matlab <br />)
[//]: #  (utilities/threshold_components.m        | mild post-processing of spatial components <br />)
[//]: #  (utilities/view_patches.m                | plotting of each found component and its temporal activation <br />)
