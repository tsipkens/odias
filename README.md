# ODIAS
(***O***ne-***D***imensional ***I***nversion of ***A***erosol ***S***ize distributions)

**NOTE**: This is currently an experimental code that is being validated against experimental data. 

This MATLAB program is designed to invert aerosol size distributions for a range of devices, with a focus on methods in the Bayesian framework. 

## Setup

This program has three dependences that are included as git submodules: 

1. The **tfer_pma** submodule, available at https://github.com/tsipkens/mat-tfer-pma, contains MATLAB code to compute the transfer function of particle mass analyzers (including the centrifugal particle mass analyzer and aerosol particle mass analyzer) and to compute basic aerosol properties. Functions in this submodule are necessary to compute the kernel (the quantity that related aerosol measurements  by a range of instruments to their underlying particle size distributions). As such, while this package is primarily necessary if considering particle mass analyzer transfer functions, the package also includes basic functions for computing particle mobility necessary for computing the DMA transfer function. 

2. The **autils** submodule, available at https://github.com/tsipkens/autils, adds utilities for basic aerosol calculations, including size conversions and distribution moment calculation. 

3. The **cmap** submodule, available at https://github.com/tsipkens/cmap, adds perceptually uniform colormaps to the program. This submodule is optional in that one could also replace references in existing scripts to the colormaps that would otherwise be in that package. 

As a result, the folders corresponding to these submodules will initially be empty. Their are multiple routes to downloading these submodules. If using git, one can initially clone the repository using 

```shell
git clone git://github.com/tsipkens/odias --recurse-submodules
```

which will automatically download the submodules when downloading overall program. Alternatively, the submodules can be downloaded manually from the above sources and placed in the `cmap/` and `tfer_pma/` folders. In either case, to be used directly, these packages should then be added to the Matlab path at the beginning of any script using

```Matlab
addpath('cmap', 'tfer_pma', 'autils');
```

For **tfer_pma**, functions in the **kernel** package will add this folder to the path automatically, whenever necessary, such that it is not necessary to explicitly include the above command in high level scripts. 

----

#### Acknowledgements and credit

This code was primarily written by Timothy Sipkens at the University of British Columbia and the National Research Council of Canada. Some code is taken from a corresponding two-dimensional size distribution inversion code available [here][mat2d], with the corresponding acknowledgements. 

Two-step Tikhonov inversion follows from the work of [Huckle and Sedlacek (2012)][Huckle2012], as recently explored by [Petters (2021)][Petters2021]. The Twomey method follows from [Twomey (1975)][Twomey1975]. The Twomey-Markowski method follows from the adaptation to the Twomey method by [Markowski (1987)][Markowski1987]. Exponential distance method adapts the method presented by [Sipkens et al. (2020)][Sipkens2020] for two dimensions to the one-dimensional case here. Function evaluating the Fuchs unipolar charging model was originally written by Tyler Johnson, with performance improvements during inclusion in this codebase. 

We also wish to acknowledge a competing code by Petters available [here][PettersCode], which contains regularization tools for aerosol size distribution inversion in Julia. 



[Huckle2012]: https://onlinelibrary.wiley.com/doi/abs/10.1002/pamm.201210310
[PettersCode]: https://github.com/mdpetters/RegularizationTools.jl
[Twomey1975]: https://www.sciencedirect.com/science/article/pii/0021999175900285
[Markowski1987]: https://www.tandfonline.com/doi/abs/10.1080/02786828708959153
[Sipkens2020]: https://doi.org/10.1016/j.jaerosci.2020.105565
[Petters2021]: https://amt.copernicus.org/preprints/amt-2021-51/
[mat2d]: https://github.com/tsipkens/mat-2d-aerosol-inversion
