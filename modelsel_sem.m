% This file is a function used for Exploratory SEM. It concerns two step
% optimization, which is described in section3.4, Exploratory SEM,  from
% 
% A. Pruttiakaravanich and J. Songsiri, "Convex Formulation for Regularized Estimation
% of Structural Equation Models"
% 
% This function requires input as
%
%   1) no_candidate     : number of candidate model (or number of gamma)
%   2) Y                : observed matrix in a form of Y = [variables (n) x time points (N)]
%   3) indA0            : index of A that need to be set as 0 from prior assumption
%   4) alpha(optional)  : the default value of alpha is "alpha = min(eig(S))" 
%
% and it returns output as
%   1) A        : a set of path matrices from all candidate models
%   2) BICmin   : minimum BIC scores
%   3) A_bic    : the model selected from BIC scores
%   4) AICmin   : minimum AIC scores
%   5) A_aic    : the model selected from AIC scores
%   6) AICCmin 	: minimum AICC scores
%   7) A_aicc  	: the model selected from AICC scores
%   8) KICmin   : minimum BIC scores
%   9) A_kic    : the model selected from KIC scores
%   10) KICCmin	: minimum KICC scores
%   11) A_kicc  : the model selected from KICC scores
% 
% Developer: Anupon Pruttiakaravanich, Dr. Jitkomut Songsiri
% 
%% exploratory SEM
function [output] = modelsel_sem(no_candidate,Y,indA0,alpha)

if(no_candidate <= 0)
    display('WARNING!! number of candidate model must be positive')
    output = 0;
    return;
end

S = cov(Y');

size_Y = size(Y);
N = size_Y(2);                          %number of time points
n = size_Y(1);                          %number of variables (ROIs)   

if ~exist('alpha','var')
    % alpha does not exist, so default it to set alpha = min(eig(S))
    alpha = min(eig(S));      
end

%no assumption for path matrix
allind = (1:n^2)'; 
indnotA0 = setdiff(allind,indA0,'rows');


%gamma
tmpS = alpha*eye(n)-S;
gamma_max = (1/alpha)*max(abs(tmpS(indnotA0)));

%mode : gamma with log scale
gamma = (gamma_max)*[0 logspace(-4,0,no_candidate-1)];

for k = 1:no_candidate
    if (k == 1)
        %solve convex confirm SEM
        indA = indA0;
        [sol_confirm, hist] = confirmatory_sem_ppxa(S,indA,alpha,'maxiter',100000,'tolfun',1e-5,'tolx',1e-5);            
        Acon(:,:,k) = sol_confirm.A;           
        fprintf('DONE solving convex confirmatory sem: candidate model(%d)\n',k);
    else
        %solve sparse SEM
        [sol_sparse,~] = sparse_sem_ppxa(S,gamma(k),indA0,alpha,'maxiter',100000,'tolfun',1e-5,'tolx',1e-5);    
        indA = find(sol_sparse.A == 0);             	%#zero entries in A will be used in solving confirm SEM
        fprintf('DONE solving convex sparse sem: candidate model(%d)\n',k);

        %solve convex confirm SEM
        [sol_confirm, hist] = confirmatory_sem_ppxa(S,indA,alpha,'maxiter',100000,'tolfun',1e-5,'tolx',1e-5);
        Acon(:,:,k) = sol_confirm.A;

        fprintf('DONE solving convex confirmatory sem: candidate model(%d)\n',k);
    end 
    
    %compute BIC scores of each model
    %#effective parameter (#overall parameters - #parameter in A which are set to be 0 from the constraint indA = 0)
    d = 2*n - 2*length(indA);                   %number of effective parameters
    minus2L = N*(hist.objval(end));             %-2L

    %information criterion sccores
    BIC(k) = minus2L + d*log(N);
    AIC(k) = minus2L + 2*d;
    AICC(k) = minus2L + (2*d*N)/(N-d-1);
    KIC(k) = minus2L + 3*d;
    KICC(k) = minus2L + ((d+1)*((3*N)-d-2)/(N-d-2)) + (d/(N-d));  
end
output.A = Acon;

[BICmin, ind_BICmin] = min(BIC);
output.BIC = BICmin;
output.A_bic = Acon(:,:,ind_BICmin);

[AICmin, ind_AICmin] = min(AIC);
output.AIC = AICmin;
output.A_aic = Acon(:,:,ind_AICmin);

[AICcmin, ind_AICCmin] = min(AICC);
output.AICC = AICcmin;
output.A_aicc = Acon(:,:,ind_AICCmin);

[KICmin, ind_KICmin] = min(KIC);
output.KIC = KICmin;
output.A_kic = Acon(:,:,ind_KICmin);

[KICCmin, ind_KICCmin] = min(KICC);
output.KICC = KICCmin;
output.A_kicc = Acon(:,:,ind_KICCmin);
end