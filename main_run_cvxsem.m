%% ADMM consensus for sparse SEM
% the algorithm is ADMM based on Boyd book page 50 
% Distributed Optimization and Stat learning via ADMM


clear; clc; close all;

addpath('../');

n = 20; 
% generate random nonzero index indA 
density = 0.2;
R = eye(n)+sprandn(n,n,density);
indA = find(R); % index where P(A) = 0  or P(X2) = I
allind = (1:n^2)'; 
indnotA = setdiff(allind,indA,'rows');
I = eye(n);
subplot(121);plot_spy(indA,n,'s');title('index of A = 0');
subplot(122);plot_spy(indnotA,n,'s');title('index of not A');

% parameters for original problem

Y = randn(2*n,n); S0 = cov(Y);
lambdamin = min(eig(S0));
alpha0 = 1*lambdamin; 
tmpS = alpha0*eye(n)-S0;
gamma_max = (1/alpha0)*max(abs(tmpS(indnotA)));
gamma = 0.05*gamma_max;

% it requires log_det function by CVX

%% Sparse SEM comparison

[s1,h1] = sparse_sem_admm(S0,gamma,indA,alpha0); % currently used version
[s2,h2] = sparse_sem_ppxa(S0,gamma,indA,alpha0); % PPXA algorithm



%% Confirmatory SEM comparison

[sc, hc] = confirmatory_sem_admm(S0,indA,alpha0); % current version

