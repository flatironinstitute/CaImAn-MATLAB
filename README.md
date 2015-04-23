# ca_source_extraction

The code implements a method for simultaneous source extraction and spike inference from large scale calcium imaging movies. The code is suitable for the analysis of somatic imaging data. Implementation for the analysis of dendritic/axonal imaging data will be added in the future. The algorithm is presented in

Pnevmatikakis, E. A., Gao, Y., Soudry, D., Pfau, D., Lacefield, C., Poskanzer, K., ... & Paninski, L. (2014). A structured matrix factorization framework for large scale calcium imaging data analysis. arXiv preprint arXiv:1409.2903. http://arxiv.org/abs/1409.2903

function name                           | description 
----------------------------------------|-----------------------------------
demo_script.m                           | wrapper code <br />
update_spatial_components.m             | update spatial components given temporal components and data <br />
update_temporal_components.m            | update temporal components given spatial components and data <br />
merge_ROIs.m                            | merge spatially overlapping components that are temporally correlated <br />
utilities/arpfit.m                      | estimation of noise level for every pixel and global time constants <br />
utilities/com.m:                        | calculation of the center of mass of each component <br />
utilities/graph_connected_comp.m        | finds the connected components in a graph <br />
utilities/greedyROI2d.m                 | Greedy method for initializing the spatial and temporal components <br />
utilities/lars_regression_noise.m       | solve a basis pursuit denoising problem using the LARS algorithm <br />
utilities/make_G_matrix.m               | construct a convolution/deconvolution matrix <br />
utilities/make_patch_video.m            | construct a video that displays the results of the algorithm <br />
utilities/order_ROIs.m                  | order found components based on their temporal activation and spatial size <br />
utilities/plain_foopsi.m                | projection of fluorescence onto the cone formed by the indicator dynamics 
utilities/plot_contours.m               | contour plot of found components and creation of a json file <br />
utilities/tiff_reader.m                 | loading a tiff stack into matlab <br />
utilities/view_patches.m                | plotting of each found component and its temporal activation <br />
