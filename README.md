# cvxsem
Convex Formulation for Regularized Estimation of Structural Equation Models

This project include proximal algorithms (ADMM and PPXA) to solve the confirmatory and sparse SEM (structural equation modelling). The algorithms are based on the analysis in the manuscript

A. Pruttiakaravanich and J, Songsiri, Convex Formulation for Regularized Estimation of Structural Equation Models, Submitted to Signal Processing, 2019.

The main functions are 
1) confirmatory_sem_ppxa.m
2) sparse_sem_ppxa.m
3) main_code.m
4) modelsel_sem.m
4) other algorithms are tested in confirmatory_sem_admm, confirmatory_sem_cvx, sparse_sem_admm, sparse_sem_cvx

To see an example of solving SEM problems, just run 
main_code.m 

or if you would like to test exploratory SEM 
run  modelsel_sem.m

According to our paper, the codes used in running experiments are

4) gen_datasem.m
5) performance_sem.m
6) compare_sem_regsem.m
7) sem_abide_fmri.m  (application on learning brain connectivity from fMRI data)
8) sem_air.m (application on learning a network from air pollution data)

Other sub folders contain raw data of experimental results. 
