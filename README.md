# ODIAS
(***O***ne-***D***imensional ***I***nversion of ***A***erosol ***S***ize distributions)

**NOTE**: This is experimental code and could change drastically at any time. 

This Matlab program is designed to invert aerosol size distributions for a range of devices, with a focus on methods in the Bayesian framework. 

## Setup

This program has two dependences that are included as git submodules: 

1. The **tfer_pma** submodule, available at https://github.com/tsipkens/mat-tfer-pma, contains Matlab code to compute the transfer function of particle mass analyzers (including the centrifugal particle mass analyzer and aerosol particle mass analyzer) and to compute basic aerosol properties. Functions in this submodule are necessary to compute the kernel (the quantity that related aerosol measurements  by a range of instruments to their underlying particle size distributions). As such, while this package is primarily necessary if considering particle mass analyzer transfer functions, the package also includes basic functions for computing particle mobility necessary for computing the DMA transfer function. 

2. The **cmap** submodule, available at https://github.com/tsipkens/cmap, adds perceptually uniform colormaps to the program. This submodule is optional in that one could also replace references in existing scripts to the colormaps that would otherwise be in that package. 

As a result, the folders corresponding to these submodules will initially be empty. Their are multiple routes to downloading these submodules. If using git, one can initially clone the repository using 

```shell
git clone git://github.com/tsipkens/odias --recurse-submodules
```

which will automatically download the submodules when downloading overall program. Alternatively, the submodules can be downloaded manually from the above sources and placed in the `cmap/` and `tfer_pma/` folders. In either case, to be used directly, these packages should then be added to the Matlab path at the beginning of any script using

```Matlab
addpath('tfer_pma', 'cmap');
```

For **tfer_pma**, functions in the **kernel** package will add this folder to the path automatically, whenever necessary, such that it is not necessary to explicitly include the above command in high level scripts. 
