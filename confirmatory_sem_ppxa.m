
function [sol, history] = confirmatory_sem_ppxa(S0,indA,alpha0,varargin)
%% Parallel Proximal Algorithm (PPXA) for confirmatory SEM
% 
% SPARSE_SEM_PPXA solves the problem
% 
% minimize -logdet(X1) + Tr(SX1) 
% subject to  X >= 0, 0 <= X4 <= alpha*I, P(X2) = I
% 
% DEFAULT USAGE: 
% 
% [SOL, H] = confirmatory_sem_ppxa(S0,indA,alpha0)
% 
% USAGE WITH OPTIONS:
% 
%
% [SOL, H] = confirmatory_sem_ppxa(S0,indA,alpha0,'initial',X0,'maxiter',1000,'tolfun',1e-5,'tolx',1e-5,'freqprint',5)
% 
% Required input arguments: 
%   1) S0       : sample covariance matrix
%   2) indA     : linear indices that A_ij is zero (index set I_A)
%   3) alpha0 	: a problem parameter. normally we suggest to set alpha = min(eig(S0)).
% 
% Optional arguments:
% 
%   5) 'initial' : intitial solution (X0) must be a symmetric matrix with size 2n*2n
%   6) 'maxiter' : a maximum number of iterations  (positive integer) 
%   7) 'tolfun'  : tolerance of the relative change of cost function(f)
%   8) 'tolx'    : tolerance of the relative change of solution (X)
%   9) 'freqprint': frequency of printing algorithm iterations on screen
%
% Returned output arguments: 
% 
%   1) SOL: a structure variable with the fields
%      SOL.X : a symmtric X = [X1 X2 ; X2 X4] each Xk is n x n
%      SOL.A : our estimated path matrix
%   2) HISTORY: a struture variable containing historical algorithm values
%      objective values, relative change in X and objective function, residual norms
% 
% 
% We implement PPXA based on Combett book (Proximal Algorithm)
% 
% minimize f1(X1) + f2(X2)
% subject to X1 = X2
% where 
% f1(Y) = -logdet(Y1)+ Tr(S*Y1) + I{ P(Y2) = I } + I{ 0 <= X4 <= alpha*I }
% f2(Y) = I{ Y >= 0} 
% 
% I{ X in C} denotes an indicator function)
% 
% We also solve the scaled problem using (S,alpha) = (beta*S0,1) 
% where S0 is the original sample covariance and beta = 1/alpha0
% the original solution is obtained by scaling back (check Proposition 2 in the
% paper)
% 
% 
% The algorithm is based on PPXA page 203 in Combett book (ALG 10.27)
% described in 
% 
% A. Pruttiakaravanich and J. Songsiri, "Convex Formulation for Regularized Estimation
% of Structural Equation Models"% 
% 
% Additional note: it requires 'det_rootn' function by CVX to compute log_det of X > 0.
% 
% Author: Anupon Pruttiakaravanich and Jitkomut Songsiri
%  Date: May 2, 2019
% 
%%  Scaled problem
beta = 1/alpha0; % scaling factor
alpha = 1; % used in scaled problem
S = beta*S0; % use this S to solve the scaled problem
n = size(S0,1);
m = 2;  
lambda = 1.8; 
rho = 0.1; % PPXA algorithm parameters (tuned)


%setting optional parameters (X0, criterion, maxiter, tolfun, tolx, eabs, rel)
defaultMaxIter = 50000;     %maximum iteration
defaultTolFun = 1e-7;       %relative change of f
defaultTolX = 1e-7;         %relative change of X (solution)
defaultEabs = 1e-7;         %absolute tol for residual
defaultErel = 1e-7;         %relative tol for residual
defaultFreqPrint = 10;      %printing frequency
defaultCriterion = 'rel';   %default criterion : relative error

Y1 = S\eye(n); Y4 = alpha*eye(n); Y2 = zeros(n); Y0 = [Y1 Y2';Y2 Y4];
defaultX0 = Y0;

p = inputParser;

validmatrixinputS = @(x) ismatrix(x) && issymmetric(x) && (sum(size(x) == [n,n]) == 2);
validmatrixinputX = @(x) ismatrix(x) && issymmetric(x) && (sum(size(x) == [2*n,2*n]) == 2);
validScalarPosNumForMaxIter = @(x) isnumeric(x) && isscalar(x) && (x >= 10);
validScalarPosNumForCriterion = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x < 1);
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validScalarPosNumForPrint = @(x) isnumeric(x) && isscalar(x) && (x > 1);
validCriterion = {'rel','res'};
checkCriterion = @(x) any(validatestring(x,validCriterion));

addRequired(p,'S0',validmatrixinputS);
addRequired(p,'indA');
addRequired(p,'alpha0',validScalarPosNum);
addOptional(p,'criterion',defaultCriterion,checkCriterion);
addParameter(p,'initial',defaultX0,validmatrixinputX);
addParameter(p,'maxiter',defaultMaxIter,validScalarPosNumForMaxIter);
addParameter(p,'tolfun',defaultTolFun,validScalarPosNumForCriterion);
addParameter(p,'tolx',defaultTolX,validScalarPosNumForCriterion);
addParameter(p,'eabs',defaultEabs,validScalarPosNumForCriterion);    
addParameter(p,'erel',defaultErel,validScalarPosNumForCriterion);    
addParameter(p,'freqprint',defaultFreqPrint,validScalarPosNumForPrint);

parse(p,S0,indA,alpha0,varargin{:}); 

% numerical parameters
criterion = p.Results.criterion;      	% set stopping criterion
MAXITER = p.Results.maxiter;            % set number of maximum iterations
Y0 = p.Results.initial;             	% set initial guess X0
FREQ_PRINT = p.Results.freqprint; 
if(strcmp(criterion,'res'))
    E_abs = p.Results.eabs;             % absolute tolerance
    E_rel = p.Results.erel;             % relative tolerance                         
else
    TOL_RELF = p.Results.tolfun;        % relative change of cost objective
    TOL_RELX = p.Results.tolx;          % relative change of solution    
end      
PRINT_RESULT = 1;


% FOUR VARIABLES:  Y,P,Pbar (auxillary), X (main)
Y = repmat(Y0,1,1,m);
X = Y0; w = 1/m;
P = repmat(Y0,1,1,m); Pbar = Y0;

f = -n*log(det_rootn(X(1:n,1:n)))+trace(S0*X(1:n,1:n));

if(strcmp(criterion,'rel'))
    
   	%print information to user
    fprintf(['----------ALGORITHM PARAMETERS----------\n',...
            'stopping criterion : relative error\n',...
            'relative change of objective function (f): %0.4d\n',...
            'relative change of solution (X): %0.4d\n',...
            'maximum iteration : %d\n',...
            'printing frequency : %d\n',...
            '----------STARTING ALGORITHM----------\n'],TOL_RELF,TOL_RELX,MAXITER,FREQ_PRINT);
        
    if (PRINT_RESULT == 1)
        fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
        'r norm', 's norm', 'rel change in f', 'rel change in x', 'objective');
    end

    for ii=1:MAXITER
        fold = f;
        Xold = X;
        Pold = Pbar;

        P(:,:,1) = proxlogdet(Y(:,:,1),S,w/rho,alpha,indA);
        P(:,:,2) = proxpdf(Y(:,:,2));

        Pbar = mean(P,3);

        for k=1:m
            Y(:,:,k) = Y(:,:,k) + lambda*(2*Pbar - X - P(:,:,k));
        end
        X = X + lambda*(Pbar-X);

            % cost objective value of scaled problem
        history.objval1(ii) = -n*log(det_rootn(X(1:n,1:n)))+trace(S*X(1:n,1:n));

        W = X; % W is the solution of original problem
        W(1:n,1:n) = X(1:n,1:n)*beta; W(n+1:end,n+1:end) = X(n+1:end,n+1:end)/beta;

        % cost objective value of original problem
        history.objval(ii) = -n*log(det_rootn(W(1:n,1:n)))+trace(S0*W(1:n,1:n));

        f = history.objval(ii);

        history.relf(ii) = abs( (fold - f)/fold ) ; 
        history.relx(ii) = norm(X-Xold)/norm(Xold);

        history.r_norm(ii) = norm( reshape(P-Pbar,m*(2*n)^2,1));
        history.s_norm(ii) = sqrt(m)*rho*norm(Pbar-Pold);

        if (PRINT_RESULT && mod(ii,FREQ_PRINT) == 0)
            fprintf('%3d\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.2f\n', ii, ...
            history.r_norm(ii), history.s_norm(ii), history.relf(ii), history.relx(ii),history.objval(ii));
        end

        if ( history.relf(ii)  <= TOL_RELF) && ( history.relx(ii) <= TOL_RELX )
            break;
        end
    end

    if ii == MAXITER
        history.converge = 0; display('The algorithm hits the max number of iterations');
    else
        history.converge = 1; display('relative change of X and f are less than the desired tolerance');
    end    
    
else
    %use residual error as stopping criterion
    %we need to calculate eps pri and eps dual instead of relative error
    %but it has not been completed yet
    
    %print information to user
    fprintf(['----------ALGORITHM PARAMETERS----------\n',...
            'stopping criterion : residual error\n',...
            'absolute residual tolerance: %0.4d\n',...
            'relative residual tolerance: %0.4d\n',...
            'maximum iteration : %d\n',...
            'printing frequency : %d\n',...
            '----------STARTING ALGORITHM----------\n'],E_abs,E_rel,MAXITER,FREQ_PRINT)    
    
    if (PRINT_RESULT == 1)
        fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
        'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
    end    

    for ii=1:MAXITER
        fold = f;
        Xold = X;
        Pold = Pbar;

        P(:,:,1) = proxlogdet(Y(:,:,1),S,w/rho,alpha,indA);
        P(:,:,2) = proxpdf(Y(:,:,2));

        Pbar = mean(P,3);

        for k=1:m
            Y(:,:,k) = Y(:,:,k) + lambda*(2*Pbar - X - P(:,:,k));
        end
        X = X + lambda*(Pbar-X);

            % cost objective value of scaled problem
        history.objval1(ii) = -n*log(det_rootn(X(1:n,1:n)))+trace(S*X(1:n,1:n));

        W = X; % W is the solution of original problem
        W(1:n,1:n) = X(1:n,1:n)*beta; W(n+1:end,n+1:end) = X(n+1:end,n+1:end)/beta;

        % cost objective value of original problem
        history.objval(ii) = -n*log(det_rootn(W(1:n,1:n)))+trace(S0*W(1:n,1:n));

        f = history.objval(ii);
        history.relf(ii) = abs( (fold - f)/fold ) ; 
        history.relx(ii) = norm(X-Xold)/norm(Xold);

        history.r_norm(ii) = norm( reshape(P-Pbar,m*(2*n)^2,1));
        history.s_norm(ii) = sqrt(m)*rho*norm(Pbar-Pold);

        %we need to compute Erel and Eabs here!!
%         history.eps_pri(ii) = ??
%         history.eps_dual(ii) = ??

        if (PRINT_RESULT && mod(ii,FREQ_PRINT) == 0)
             fprintf('%3d\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.2f\n', ii, ...
            history.r_norm(ii), history.eps_pri(ii) , history.s_norm(ii), history.eps_dual(ii), history.objval(ii));
        end

        %break program
        if (history.r_norm(ii) < history.eps_pri(ii) && ...
           history.s_norm(ii) < history.eps_dual(ii))
             break;
        end

    end

    if ii == MAXITER
        history.converge = 0; display('the algorithm hits the max number of iterations');
    else
        history.converge = 1; display('residual norms are less than a tolerance');
    end
end
sol.X = X; 
sol.X(1:n,1:n) = X(1:n,1:n)*beta;
sol.X(n+1:end,n+1:end) = X(n+1:end,n+1:end)/beta; 
% use P(:,:,1) because P(X2) = I exactly
X2 = P(n+1:end,1:n,1);
sol.X(n+1:end,1:n) = X2; sol.X(1:n,n+1:end) = X2';
sol.A = eye(n)-sol.X(n+1:end,1:n); % A = I-X2;

%% Proximal operator of loget 
% f1(X) = -logdet(X1) + Tr(SX1) + I{ P(X2) = I} + I{ 0 <= X4 <= alpha*I } 
% 
% minimize_X -logdet(X1) + Tr(SX1) + (rho/2) || Y - X ||^2
% subject to P(X2) = I,  0 <= X4 <= alpha*I 

function[Z] = proxlogdet(Y,S,rho,alpha,indA)

n = size(S,1); 
Z = zeros(2*n,2*n); % allocate Z first
I = eye(n);

% Z1
[U,D] = eig(rho*Y(1:n,1:n)-S);
d = diag(D);
Dz1 = diag(0.5*(d+sqrt(d.^2+4*rho) )/rho);
Z(1:n,1:n) = U*Dz1*U';

% Z2
Z2 = Y(n+1:end,1:n);
Z2(indA) = I(indA); % P(Z2) = I
Z(n+1:end,1:n) = Z2; Z(1:n,n+1:end) = Z2';

% Z4
[U4,D4] = eig(Y(n+1:end,n+1:end));
Dz4 = diag(max(0,min(alpha,diag(D4))));
Z(n+1:end,n+1:end) = U4*Dz4*U4';

end

%% Proximal operator of cone constraint
% minimize_{ X >= 0 } (1/2) || Y - X ||^2

function[Z] = proxpdf(Y)

[U,D] = eig(Y);
Dz = diag([max(0,diag(D))]);
Z = U*Dz*U';
end

end