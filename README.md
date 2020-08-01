**lotR**: **l**atent class analysis of **o**bservations organized by **t**ree in **R**

An R package for Latent Class Models for Observations Organized by Tree Structure (lotR)

zhenkewu badges:
[![Travis CI Build Status](https://travis-ci.org/zhenkewu/lotR.svg?branch=master)](https://travis-ci.org/zhenkewu/lotR)

**Maintainer**: Zhenke Wu, zhenkewu@umich.edu

## Table of content
- [1. Installation](#id-section1)
- [2. Overview](#id-section2)
- [2. Example](#id-section3)

<div id='id-section1'/>

Installation
--------------
```r
install.packages("devtools",repos="https://cloud.r-project.org")
devtools::install_github("zhenkewu/lotR")
```
<div id='id-section2'/>

Overview
----------
`lotR` is designed for analyzing multivariate binary observations while integrating 
additional sample related information represented by each observation's membership
in the leaves of a given tree. The observations that are closer in the tree are 
_a priori_ more likely to be grouped together and fitted by a LCM with identical 
LCM parameters. The model is built on spike-and-slab priors on the increments of 
a Gaussian diffusion process for each node of the tree. The model is self-adaptive 
in that it automatically choose the optimal grouping of observations to fit 
distinct latent class models. The posterior inferential algorithm is based on 
variational inference and can provide approximate posterior uncertainty quantifications.


**lotR** works for 

* multivariate binary responses
	-  known cut level, pre-specified grouping of observations
    -  unknown cut level, which requires self-adaptive grouping of observations


<div id='id-section3'/>

Examples 
---------

* _lotR_ is self-adaptive: leaves close in the tree means they are more likely to be grouped together.
![](inst/example_figure/lotR_self_adaptive.png)

* _lotR_ produces similar results as `BayesLCA` on fully collapsed tree (ignoring tree information)
![](inst/example_figure/comparison_with_std.png)




