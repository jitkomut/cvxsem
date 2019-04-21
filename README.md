# cvxsem
Convex Formulation for Regularized Estimation of Structural Equation Models

This project include proximal algorithms (ADMM and PPXA) to solve the confirmatory and sparse SEM (structural equation modelling). The algorithms are based on the analysis in the manuscript

A. Pruttiakaravanich and J, Songsiri, Convex Formulation for Regularized Estimation of Structural Equation Models, Submitted to Signal Processing, 2019.

The three functions are 
1) confirmatory_sem_admm
2) sparse_sem_admm
3) sparse_sem_admm

The function 2) and 3) solve the same problem but with different algorithms. Both of them are proximal algorithms. Their convergences are slightly different under various conditions of problem parameter. PPXA seems to be faster when solving dense problem and ADMM is a bit faster when solving sparse problem. PPXA requires more auxillary variables.

run

main_run_cvxsem.m

to see a usage of the main algorithms.
