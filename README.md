[![Join the chat at https://gitter.im/epnev/ca_source_extraction](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/epnev/ca_source_extraction?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

CaImAn-MATLAB
======
<img src="https://github.com/simonsfoundation/CaImAn/blob/master/docs/LOGOS/Caiman_logo_FI.png" width="500" align="right">

A Computational toolbox for large scale **Ca**lcium **Im**aging data **An**alysis.
The code implements the CNMF algorithm [[1]](#neuron) for simultaneous source extraction and spike inference from large scale calcium imaging movies. Many more features are included (see below). The code is suitable for the analysis of somatic imaging data. Improved implementation for the analysis of dendritic/axonal imaging data will be added in the future. 


## Features and methods included

* **Source extraction** 

    * Separates different sources based on constrained nonnegative matrix factorization (CNMF) [[1-2]](#neuron)
    * Deals with heavily overlaping and neuropil contaminated movies     
    * Selection of inferred sources using a [pre-trained convolutional neural network classifier](https://github.com/flatironinstitute/CaImAn-MATLAB/wiki/Component-classification-with-a-convolutional-neural-network)
    * [Component registration across different sessions/days](https://github.com/flatironinstitute/CaImAn-MATLAB/wiki/Registering-ROIs-across-different-sessions-%5C--days) 

* **Denoising, deconvolution and spike extraction**

    * Constrained foopsi method for inferring neural activity from fluorescence traces [[1]](#neuron)
    * Near online implementation using the OASIS algorihtm [[3]](#oasis)
    * MCMC algorithm for Bayesian spike inference [[4]](#mcmc)
    
* **Handling of very large datasets**

    * [Memory mapping and parallel processing in patches](https://github.com/flatironinstitute/CaImAn-MATLAB/wiki/Processing-of-large-datasets)
    
* **Motion correction**

    * Fast parallelizable non-rigid motion correction using the NoRMCorre algorithm [[5]](#normcorre). Separate standalone package can be found [here](https://github.com/simonsfoundation/NoRMCorre). It will be included in this package in the future.
    
New: Renaming to CaImAn-MATLAB
======
We moved the code into the Flatiron Institute github account and renamed the repository to CaImAn-MATLAB to bring it more in touch with the [CaImAn](https://github.com/simonsfoundation/CaImAn) Python package. Everything else is the same. The old link ```https://github.com/epnev/ca_source_extraction``` redirects here.

# Citation

If you use this code please cite the corresponding papers where original methods appeared (see References below), as well as: 

<a name="caiman"></a>[1] Giovannucci A., Friedrich J., Gunn P., Kalfon J., Koay S.A., Taxidis J., Najafi F., Gauthier J.L., Zhou P., Tank D.W., Chklovskii D.B., Pnevmatikakis E.A. (2018). CaImAn: An open source tool for scalable Calcium Imaging data Analysis. bioarXiv preprint. [[paper]](https://doi.org/10.1101/339564)

# References

The following references provide the theoretical background and original code for the included methods. 

### Deconvolution and demixing of calcium imaging data

<a name="neuron"></a>[1] Pnevmatikakis, E.A., Soudry, D., Gao, Y., Machado, T., Merel, J., ... & Paninski, L. (2016). Simultaneous denoising, deconvolution, and demixing of calcium imaging data. Neuron 89(2):285-299, [[paper]](http://dx.doi.org/10.1016/j.neuron.2015.11.037). 

<a name="struct"></a>[2] Pnevmatikakis, E.A., Gao, Y., Soudry, D., Pfau, D., Lacefield, C., ... & Paninski, L. (2014). A structured matrix factorization framework for large scale calcium imaging data analysis. arXiv preprint arXiv:1409.2903. [[paper]](http://arxiv.org/abs/1409.2903). 

<a name="oasis"></a>[3] Friedrich J. and Paninski L. Fast active set methods for online spike inference from calcium imaging. NIPS, 29:1984-1992, 2016. [[paper]](https://papers.nips.cc/paper/6505-fast-active-set-methods-for-online-spike-inference-from-calcium-imaging), [[Github repository - Python]](https://github.com/j-friedrich/OASIS), [[Github repository - MATLAB]](https://github.com/zhoupc/OASIS_matlab).

<a name="mcmc"></a>[4] Pnevmatikakis, E. A., Merel, J., Pakman, A., & Paninski, L. Bayesian spike inference from calcium imaging data. In Signals, Systems and Computers, 2013 Asilomar Conference on (pp. 349-353). IEEE, 2013. [[paper]](https://arxiv.org/abs/1311.6864), [[Github repository - MATLAB]](https://github.com/epnev/continuous_time_ca_sampler).

### Motion Correction

<a name="normcorre"></a>[5] Pnevmatikakis, E.A., and Giovannucci A. (2017). NoRMCorre: An online algorithm for piecewise rigid motion correction of calcium imaging data. Journal of Neuroscience Methods, 291:83-92 [[paper]](https://doi.org/10.1016/j.jneumeth.2017.07.031), [[Github repository - MATLAB]](https://github.com/simonsfoundation/normcorre).

Code description
=======

The best way to start is by looking at the various demos.
- [demo_script.m](https://github.com/epnev/ca_source_extraction/blob/master/demo_script.m): A simple demo with a small dataset included in the repo to display the notation and basic operations
- [demo_script_class.m](https://github.com/flatironinstitute/CaImAn-MATLAB/blob/master/demo_script_class.m): Replicates the [demo_script.m](https://github.com/epnev/ca_source_extraction/blob/master/demo_script.m) file in a cleaner way using a CNMF object.
- [demo_patches.m](https://github.com/epnev/ca_source_extraction/blob/master/demo_patches.m): A larger demo displaying the process of memory mapping and spliting the field of view in patches to be processed in parallel and then combined.
- [demo_patches_class.m](https://github.com/epnev/ca_source_extraction/blob/master/demo_patches_class.m): Similar to [demo_patches.m](https://github.com/epnev/ca_source_extraction/blob/master/demo_patches.m) but using the CNMF object.
- [run_pipeline.m](https://github.com/epnev/ca_source_extraction/blob/master/run_pipeline.m): Demo for the complete pipeline of motion correction, source separation and spike extraction for large datasets. More details about the pipeline can be found [here](https://github.com/epnev/ca_source_extraction/wiki/Complete-analysis-pipeline).
- [3D/demo_3D.m](https://github.com/epnev/ca_source_extraction/blob/master/3D/demo_3D.m): Demo for processing of 3D volumetric imaging data.

# Python

A complete analysis Python pipeline including motion correction, source extraction and activity deconvolution is performed through the package [CaImAn](https://github.com/simonsfoundation/caiman). This package also includes method for online processing of calcium imaging data and elements of behavioral analysis in head fixed mice. 

Usage and Documentation
=======
Check the demo scripts and the [wiki](https://github.com/epnev/ca_source_extraction/wiki) to get started.

Dependencies
========
The following matlab toolboxes are needed for the default parameter settings:

- Statistics and Machine Learning Toolbox
- Image processing toolbox

Depending on the settings the following toolboxes may also be required

- Neural networks toolbox (required for component classifier)
- Signal processing toolbox (recommended but not required)
- Parallel computing toolbox (recommended for large datasets but not required)
- Optimization toolbox (not required)

Depending on the settings the following packages may also be required

- The CVX library which can be downloaded from http://cvxr.com/cvx/download/ (after unpacking CVX open Matlab and run cvx_setup from inside the CVX directory to properly install and add CVX to the Matlab path). **CVX is no longer required**.
- SPGL1 package from https://github.com/mpf/spgl1 (for solving constrained_foopsi using SPGL1)

# Developers

This package is mainly developed and maintained by [Eftychios A. Pnevmatikakis](https://github.com/epnev) (Flatiron Institute, Simons Foundation) with help from a lot of [contributors](https://github.com/flatironinstitute/CaImAn-MATLAB/graphs/contributors).

# Acknowledgements

Special thanks to the following people for letting us use their datasets for our various demo files:

* Weijian Yang, Darcy Peterka, Rafael Yuste, Columbia University
* Sue Ann Koay, David Tank, Princeton University
* Diego Pacheco Pinedo, Mala Murthy, Princeton University
* Clay Lacefied, Randy Bruno, Columbia University


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
