---
title: "GMM-Lasso Replication Package"
author: Zhentao Shi
date: February 2015
output: pdf_document
---


This repository hosts functions that generate the simulation results in ["Estimation of Sparse Structural Parameters with Many Endogenous Variables"](http://www.tandfonline.com/doi/full/10.1080/07474938.2015.1092805) (2016).

## Functionality

* **master_DGP1.R** is the master file.
* **estimation.R** is the main function for the estimation and the calculation of the summary statistics such as the MSE, squared-bias and variance. It also has many small functions to do the background work.
* **DGP.R** only contains functions that generate the data.

