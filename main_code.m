% This file is a main file based on paper,
% 
% A. Pruttiakaravanich and J. Songsiri, "Convex Formulation for Regularized Estimation
% of Structural Equation Models"
% 
% Our main formulations are 1) confirmatory SEM (formulation (6)) and
% sparse SEM (formulation (9).
% 
% Our main algorithms for solving the two formulations are 
% 1) ADMM (Stephen Boyd book) and 2) PPXA (Combette book)
% 
% This file contains 3 sections as follows.
% 
% section 1:
% We provide two examples of solving our confirmatory sem by using ADMM and PPXA. 
% The correctness of results from two algorithm are
% compared with a generic convex program solver, CVX. We also plot the performance between
% these two algorithms.
% 
% section 2:
% We provide examples of solving our sparse sem by using ADMM and PPXA.
% The correctness of results from two algorithm are
% compare with a generic convex solver, CVX. We also plot the performance between
% these two algorithms.
% 
% section 3:
% We provide the use of exploratory SEM by using a function, modelsel_sem.m
% The output from this function is the best path matrix selected from 5
% model selection criteria as BIC, AIC, AICc, KIC and KICc.
% 
% Developer: Anupon Pruttiakaravanich, Dr. Jitkomut Songsiri
% Date: May 2, 2019


%% example of solving convex confirmatory sem
% For the problem setting of this section, we assume that
% we need to explore a relation structure among 20 observed variables (n = 20)
% with some prior assumptions about relation structure. The data we have is from the
% observed matrix (Y) which is measured from 100 time points (N=100). We then
% compute the sample covariance matrix (S) from Y and suggest to set alpha = min(eig(S))
% we solve our confirmaty sem by 2 algorithms, PPXA and ADMM. The relation
% structure among the observed variable can be seen from s1.A or s2.A.
% The correctness of results from two algorithm are
% compare with a generic solver, CVX. We also plot the performance between
% these two algorithm.

clear; clc; close all;
addpath('../');

%number of observed variable
n = 20; 

% generate random nonzero index indA 
density = 0.2;
R = eye(n)+sprandn(n,n,density);
indA = find(R);                             % index where Aij = 0 from prior assumption

% parameters for original problem
N = 100;                                    %number of observations (sampples)
Y = randn(N,n);                             %observed matrix
S0 = cov(Y);                             	%computing a sample covariance matrix

lambdamin = min(eig(S0))
alpha0 = lambdamin;                         %alpha is a problem parameter.
                                            %we suggest to set alpha = min(eig(S0)
                                            
                                            
X0 = [S0/alpha0 eye(n); eye(n) alpha0*(eye(n))];    %initial solution (optional)

%solving convex confirmatory sem by using generic solver(CVX)
[X1,pcvx] = confirmatory_sem_cvx(S0,alpha0,indA);

%solving convex confirmatory sem, formulation(6), by using ppxa algorithm

tic;
[s1, h1] = confirmatory_sem_ppxa(S0,indA,alpha0,'initial',X0,'tolfun',1e-6,'tolx',1e-6);
toc;

%solving convex confirmatory sem by using admm algorithm
tic;
[s2, h2] = confirmatory_sem_admm(S0,indA,alpha0,'initial',X0,'tolfun',1e-6,'tolx',1e-6);
toc;

%comparing results wit CVX
[pcvx h1.objval(end) h2.objval(end)]
[norm(s1.X-X1) norm(s2.X-X1)]

%ploting performance of ppxa and admm algorithm
figure;
n1 = length(h1.objval); n2 = length(h2.objval);
semilogy(1:n1,abs((h1.objval-pcvx)/pcvx),...
    1:n2,abs((h2.objval-pcvx)/pcvx));
legend('PPXA','ADMM');
xlabel('number of iterations');
ylabel('relative error in f');

%% example of solving sparse_sem
% For the problem setting of this section, we assume that
% we need to explore a relation structure among 20 observed variables (n = 20)
% with some prior assumptions about relation structure. The data we have is from the
% observed matrix (Y) which is measured from 100 time points (N=100). We then
% compute the sample covariance matrix (S) from Y and suggest to set alpha = min(eig(S))
% we compute gamma_max and heuristically choose gamma = 0.1*gamma_max.
% we solve our confirmaty sem by 2 algorithms, PPXA and ADMM. The relation
% structure among the observed variable can be seen from s1.A or s2.A.
% The correctness of results from two algorithms are
% compare with a generic solver, CVX. We also plot the performance between
% these two algorithm.
clear; clc; close all;
addpath('../');

%number of observed variable
n = 20;

% parameters for original problem
N = 100;                                    %number of observations (sampples)
Y = randn(N,n);                             %observed matrix
S0 = cov(Y);                             	%computing a sample covariance matrix

lambdamin = min(eig(S0))
alpha0 = lambdamin;                         %alpha is a problem parameter.
                                            %we suggest to set alpha = min(eig(S0)

% generate random nonzero index indA 
density = 0.2;
R = eye(n)+sprandn(n,n,density);
indA = find(R);                             % index where Aij = 0 from prior assumption
allind = (1:n^2)'; 
indnotA = setdiff(allind,indA,'rows');

% calculatring gamma_max, the value of gamma that forces all entries in returned path matrix to be zero. 
% the derivation of gamma_max is described in section 9.2 in our paper.
tmpS = alpha0*eye(n)-S0;
gamma_max = (1/alpha0)*max(abs(tmpS(indnotA)));
gamma = 0.1*gamma_max; % selecting gamma (regularization parameter) for our sparse_sem

%solving sparse sem using generic solver(CVX)
[X1,pcvx] = sparse_sem_cvx(S0,alpha0,gamma,indA,indnotA);

%solving sparse sem using ppxa algorithm
tic;
[s1, h1] = sparse_sem_ppxa(S0,gamma,indA,alpha0);
toc;

%solving sparse sem using admm algorithm
tic;
[s2, h2] = sparse_sem_admm(S0,gamma,indA,alpha0);
toc;

%example of using our admm solver with some modified options
X0 = zeros(2*n,2*n);
tic;
[s3, h3] = sparse_sem_admm(S0,gamma,indA,alpha0,'initial',X0,'criterion','res','maxiter',1000,'eabs',1e-5,'erel',1e-5,'freqprint',100);
toc;

%comparing results wit CVX
[pcvx h1.objval(end) h2.objval(end) h3.objval(end) ]
[norm(s1.X-X1) norm(s2.X-X1) norm(s3.X-X1)]

%ploting performance of ppxa and admm algorithm
figure;
n1 = length(h1.objval); n2 = length(h2.objval);
semilogy(1:n1,abs((h1.objval-pcvx)/pcvx),...
    1:n2,abs((h2.objval-pcvx)/pcvx));
legend('PPXA','ADMM');
xlabel('number of iterations');
ylabel('relative error in f');

%% 3) the use of exploratory SEM from section 3.4). 
% For the problem setting of this section, we assume that we have an
% observed matrix from 10 variables and 100 time points as Y with size
% 10x100 and we need to find the best path matrix that describes the relation
% among these 10 variables. We impose the number of candidate model to be 25 
% (25 value of gamma from 0 to gamma_max) and we run our exploratory sem
% via a funcion called modelsel_sem. The function returns the best path
% matrix that minimizes BIC, AIC, AICC, KIC and KICC scores.
% 
clear; clc; close all;

addpath('../');
n = 10;                             %number of observed variables    
N = 100;                            %number of observations (samples)    

%daga generating
Y = randn(n,N);                     %observed matrix (number of variable (n) x time point (N))
no_candidate = 25;

% generate random nonzero index indA 
density = 0.2;
R = eye(n)+sprandn(n,n,density);
indA = find(R);                     %index of A where P(A) = 0
    
[output] = modelsel_sem(no_candidate,Y,indA)
                                    %output contains a set of candidate A
                                    %and the best A from each model
                                    %selection criterion

% note that the user should adjust options in the algorithm (in
% modelsel_sem.m) when solving problems in high dimension