


This repository hosts functions that generate the simulation results in

Zhentao Shi (2016): ["Estimation of Sparse Structural Parameters with Many Endogenous Variables"](http://www.tandfonline.com/doi/full/10.1080/07474938.2015.1092805), *Econometric Reviews*, 35(8-10): 1582-1608

## Functionality of Scripts in v1.0

* **master_DGP1.R** is the master file.
* **estimation.R** is the main function for the estimation and the calculation of the summary statistics such as the MSE, squared-bias and variance. It also has many small functions to do the background work.
* **DGP.R** only contains functions that generate the data.

## Works for future research

* `solution_path_of_gmm_lasso.lyx` contains idea to transform LARS for GMM_Lasso. I am interested in developing it into a note.
* `super_master_inference.R` is for another project for the asymptotic distribution of GMM-Lasso. `func_inference.R` contains the workhorse files.
* `cv.lars.2sls` is a modified version to implement cross-validation to choose the tuning parameter.
