# ODIAS
(***O***ne-***D***imensional ***I***nversions of ***A***erosol ***S***ize distributions)

This Matlab program is designed to invert aerosol size distributions for a range of devices, with a focus on methods in the Bayesian framework. 

## Setup

This program depends on the **tfer_pma** submodule, available at https://github.com/tsipkens/mat-tfer-pma. This submodule contains Matlab code to compute the transfer function of particle mass analyzers (including the centrifugal particle mass analyzer and aerosol particle mass analyzer) and to compute basic aerosol properties. Functions in this submodule are necessary to compute the kernel (the quantity that related aerosol measurements  by a range of instruments to their underlying particle size distributions). 

As a result, this folder will initially be empty. Their are multiple routes to downloading these submodules. If using git, one can initially clone the repository using 

```shell
git clone git://github.com/tsipkens/mat-2d-aerosol-inversion --recurse-submodules
```

which will automatically download the submodules when downloading overall program. Alternatively, the submodules can be downloaded manually from the above sources and placed in `tfer_pma/` folder. 

In either case, to be used directly, these packages should then be added to the Matlab path at the beginning of any script using

```Matlab
addpath('tfer_pma', 'cmap');
```