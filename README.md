# Ps lifetime reconstruction tools

This is a set of non-linear minimisation tools, written in [R language](https://cran.r-project.org/). The main objective is to fit [positronium](https://en.wikipedia.org/wiki/Positronium) (Ps) lifetime spectra, either separate or as a multi-voxel array representing Positronium lifetime imaging (PLI). More on the latter [here](https://doi.org/10.1109/TMI.2024.3357659), [here](https://koza.if.uj.edu.pl/files/d0bf762aba4f68695caac0f1d5745279/09059856.pdf) or [here](https://arxiv.org/abs/2501.04145).


## Prerequisites
Tested on Linux Ubuntu 22.04 LTE with R version 4.1.3 and on CERN CentOS Linux release 7.9.2009 (Core) with R version 4.1.1. The scripts are designed to operate with OpenMP. You would also need a GCC compiler and the following R packages installed:
* [jsonlite](https://cran.r-project.org/web/packages/jsonlite/index.html)
* [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)
* [minpack.lm](https://cran.r-project.org/web/packages/minpack.lm/index.html)
* [ks](https://cran.r-project.org/web/packages/ks/index.html)
