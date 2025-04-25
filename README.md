# Ps lifetime reconstruction tools

This is a set of non-linear minimisation tools, written in [R language](https://cran.r-project.org/), to fit [positronium](https://en.wikipedia.org/wiki/Positronium) (Ps) lifetime spectra, either separate or as a multi-voxel array representing Positronium lifetime imaging (PLI). More on the latter [here](https://doi.org/10.1109/TMI.2024.3357659), [here](https://koza.if.uj.edu.pl/files/d0bf762aba4f68695caac0f1d5745279/09059856.pdf) or [here](https://arxiv.org/abs/2501.04145) .


## Prerequisites
Tested on Linux Ubuntu 22.04 LTE with R version 4.1.3 and on CERN CentOS Linux release 7.9.2009 (Core) with R version 4.1.1. The scripts are designed to operate with OpenMP. You would also need a GCC compiler and the following R packages installed:
* [jsonlite](https://cran.r-project.org/web/packages/jsonlite/index.html)
* [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)
* [minpack.lm](https://cran.r-project.org/web/packages/minpack.lm/index.html)
* [ks](https://cran.r-project.org/web/packages/ks/index.html)

## The fitting of Ps lifetime spectra
The curve fitting is based on [Levenberg-Marquardt Algorithm](https://www.sciencedirect.com/topics/engineering/levenberg-marquardt-algorithm), aka damped least-squares, which is used to solve non-linear least squares problems. The minimisation of each spectrum -- a histogram of time delay between the prompt emission and electron-positron annihilation (a good description is given in [this article](https://doi.org/10.1038/s42005-020-00440-z)) -- is conducted in multiple stages:

1) the constant background is estimated from the tails of the spectrum.

2) with the background fixed, a multi-component decay model is used to fit the histogram in linear scale, extracting the smearing parameters and (optionally) the mean lifetime for direct annihilation.

3) having fixed the variables estimated at the previous stage, another minimisation is launched, but for the histogram in logarithmic scale. That allows for more precise assessment of ortho-positronium (o-Ps) lifetime.

Optionally, an additional stage can be run, with the fixed o-Ps lifetime -- to refine the parameters fit at the linear stage.

## Data requirements

The executables are in the directory ```/rsrc```. They are designed to fit either standalone histograms (```fitOneHistogramLMA2.R```, ```fitOneHistogramLMA3.R```) or multi-voxel arrays (```runMultiVoxelFitLMA2.R```, ```runMultiVoxelFitLMA3.R```).

The data format is ASCII for single histograms (first column - time delay in ns, second - histogram counts in a.u.). For multi-voxel, it is binary, requiring two files -- for voxel indices and for histograms counts as an array (each row is one histogram). The files can be stored either as an ```R``` ```.rds``` format, or binaries in little-endian (2-byte integers for voxel indices with columns -- for axes, and 4-byte floats -- for histograms). In ```R```, the structure of the two is as follows:

```
> str(voxel_IDs)
 int [1:22788, 1:3] -18 -18 -18 -18 -18 -18 -18 -18 -18 -18 ...
> str(histograms)
 num [1:22788, 1:137] 0 0 0 0 0 0 0 0 0 0 ...
```

The exemplary configuration JSON file is placed in ```/params_json```, along with the detailed description. Also, please refer to the examples in ```/examples``` directory.

(c) 2025 Roman Shopa
