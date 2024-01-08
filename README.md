# ODIAS
(***O***ne-***D***imensional ***I***nversion of ***A***erosol ***S***ize distributions)

**NOTE**: This is currently an experimental code that is being validated against experimental data. 

This MATLAB program is designed to invert aerosol size distributions for a range of devices, with a focus on methods in the Bayesian framework. 

## Setup

This program has three dependencies that are included as git submodules: 

1. The **tfer** submodule, available at https://github.com/tsipkens/tfer, contains MATLAB code to compute transfer functions of various kinds, including DMAs, PMAs, and chargers. This is an essential submodule and itself contains the **mat-tfer-pma** submodule available at https://github.com/tsipkens/mat-tfer-pma. 

2. The **autils** submodule, available at https://github.com/tsipkens/autils, adds utilities for basic aerosol calculations, including size conversions and distribution moment calculations. 

3. The **cmap** submodule, available at https://github.com/tsipkens/cmap, adds perceptually uniform colormaps to the program. This submodule is optional in that one could also replace references in existing scripts to the colormaps that would otherwise be in that package. 

As a result, the folders corresponding to these submodules will initially be empty. Their are multiple routes to downloading these submodules. If using git, one can initially clone the repository using 

```shell
git clone git://github.com/tsipkens/odias --recurse-submodules
```

which will automatically download the submodules when downloading overall program. Alternatively, the submodules can be downloaded manually from the above sources and placed in the corresponding folders. In either case, to be used directly, these packages should then be added to the Matlab path at the beginning of any script using

```Matlab
addpath autils tfer cmap;
```

## Packages

Packages refer to folders prefaced with a **+** symbol and are used to group similar functions. Use of functions in these packages requires one to use the package name in the function call. For example, when calling the `tikhonov(...)` function in the **+invert** package, one must use `invert.tikhonov(...)`. Packages include: 

### +kernel

The **+kernel** package, which contains functions to couple transfer functions and charging fractions from the **+tfer** package to describe systems of classifiers. 

### +invert

The **+invert** package contains methods that can be used to invert data to find size distributions, often by adding prior information to stabalize the inversion (e.g., Twomey-Markowski ([Markowski, 1987][Markowski1987]) and Tikhonov regularization). 

### +ac

The **+ac** package contains methods for computing the average charge and particle size transmitted by classifiers at a given size-to-charge setpoint. This is particularily useful when the full size distribution is not of interest and when using unipolar chargers with a classifier. 

Functions in this package are closely associated with a draft manuscript summarizing average charge (AC) algorithms. The iterative-average charge (IAC) algorithm predates this publication and is first defined in [Corbin et al. (2022)][Corbin2022]. 

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
[Corbin2022]: https://doi.org/10.1016/j.carbon.2022.02.037
[mat2d]: https://github.com/tsipkens/mat-2d-aerosol-inversion
