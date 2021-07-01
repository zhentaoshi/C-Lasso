# C-Lasso

This is the Matlab code for the empirical applications and simulations of

* [Liangjun Su](http://www.mysmu.edu/faculty/ljsu/), [Zhentao Shi](http://www.zhentaoshi.com/) and [Peter Phillips](http://korora.econ.yale.edu/phillips/): [“Identifying Latent Structures in Panel Data”](http://onlinelibrary.wiley.com/doi/10.3982/ECTA12560/full) (2016), *Econometrica*, Vol.84, No.6, 2215-2264.

Please contact Zhentao Shi ([zhentao.shi@cuhk.edu.hk](zhentao.shi@cuhk.edu.hk)) if you have any question about the code.

**R users please check [github.com/zhan-gao/classo](https://github.com/zhan-gao/classo).**

A follow-up paper is composed to further investigate the computational speed of C-Lasso. Please refer to:

* Zhan Gao and Zhentao Shi (2020): "[Implementing Convex Optimization in R: Two Econometric Examples](https://arxiv.org/abs/1806.10423)", arXiv:1806.10423

### Computation Environment

For the Matlab code, [CVX](http://cvxr.com/cvx/download/) must be installed to implement convex optimization.
[Mosek](https://www.mosek.com/resources/downloads) is recommended to facilitate CVX, but not necessary.

### Generic Functions

We add a folder `generic_functions` for the estimation procedures.
The functions are ready to take input and return output.

* `SSP_PLS_est.m` is a generic function to implement PLS.
* `PLS_example.m` is a minimum example of PLS.

### Development Plan after Publication

In response to demand, we may further consider

* provide user-friendly Matlab interface for general use (currently working under `generic_functions`)

We welcome interested researchers to develop the code with us.


## Note for v1.0: Replication Package
The results in the paper are generated under
* Matlab 8.5
* CVX 2.1 (http://cvxr.com/cvx/download/)  
* Mosek 7.1 (https://www.mosek.com/resources/downloads).

CVX must be installed and linked with Matlab, and Mosek is invoked as the solver through the command `cvx_solver mosek`. Without Mosek, a user can still run the code with CVX if he comments out this line.

The empirical applications can be exactly replicated by the commented `master.m` in folders
* `app_saving_PLS`: for Section 5.1
* `app_saving_PGMM`: for Section 5.1
* `app_civil_war`: for Section 5.2
* `app_democracy`: for Section S4.3

Data are also provided in each folder.

The workhorse scripts that execute the iterative algorithm in Section 3.1 of the Supplementary Material are
* `PLS_est.m`: for PLS estimation
* `PGMM_est.m`: for PGMM estimation
* `PNL_est.m`: for the PPL (Panel Probit) estimation

The scripts in folders `simulations` generate the simulation results. The master files are either `master_**` or `**_super`. Super parameters, such as `N`, `T` and `Rep`, should be provided outside of the main function or script.
