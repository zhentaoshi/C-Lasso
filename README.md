# C-Lasso Replication Package

#### Liangjun Su, Zhentao Shi and Peter Phillips

We provide all code for the empirical applications and simulations in [“Identifying Latent Structures in Panel Data”](http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2448189). Please contact Zhentao Shi ([zhentao.shi@cuhk.edu.hk](zhentao.shi@cuhk.edu.hk)) if you have any question about the code.

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
