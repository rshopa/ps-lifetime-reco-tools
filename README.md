# Ps lifetime reconstruction tools

This is a set of non-linear minimisation tools, written in [R language](https://cran.r-project.org/). The main objective is to fit [positronium](https://en.wikipedia.org/wiki/Positronium) (Ps) lifetime spectra, either separate or as a multi-voxel array representing Positronium lifetime imaging (PLI). More on the latter [here](https://doi.org/10.1109/TMI.2024.3357659), [here](https://koza.if.uj.edu.pl/files/d0bf762aba4f68695caac0f1d5745279/09059856.pdf) or [here](https://arxiv.org/abs/2501.04145).


## Prerequisites
Tested on Linux Ubuntu 22.04 LTE with R version 4.1.3 and on CERN CentOS Linux release 7.9.2009 (Core) with R version 4.1.1. The scripts are designed to operate with OpenMP. You would also need a GCC compiler and the following R packages installed:
* [jsonlite](https://cran.r-project.org/web/packages/jsonlite/index.html)
* [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)
* [minpack.lm](https://cran.r-project.org/web/packages/minpack.lm/index.html)
* [ks](https://cran.r-project.org/web/packages/ks/index.html)

The fitting is based on [Levenberg-Marquardt Algorithm](https://www.sciencedirect.com/topics/engineering/levenberg-marquardt-algorithm), aka damped least-squares (DLS), which is used to solve non-linear least squares problems. In the directory ```/rsrc```,  there are executables to fit standalone histograms (```fitOneHistogramLMA2.R```, ```fitOneHistogramLMA3.R```) and multi-voxel (```runMultiVoxelFitLMA2.R```, ```runMultiVoxelFitLMA3.R```). 

The data format is ASCII for single histograms (first column - time delay in ns, second - histogram counts in a.u.). For multi-voxel, it is binary, requiring two files -- for voxel indices and for histograms counts as an array (each row is one histogram). The files can be stored either as an ```R``` ```.rds``` format, or binaries in little-endian (2-byte integers for voxel indices with columns -- for axes, and 4-byte floats -- for histograms). In ```R```, the structure of the two is as follows:

```
> str(voxel_IDs)
 int [1:22788, 1:3] -18 -18 -18 -18 -18 -18 -18 -18 -18 -18 ...
> str(histograms)
 num [1:22788, 1:137] 0 0 0 0 0 0 0 0 0 0 ...
```
Please refer to ```/examples``` for more information.

(c) 2025 Roman Shopa
