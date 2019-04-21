function [sol, history] = confirmatory_sem_admm(S0,indA,alpha0)
% CONFIRMATORY_SEM_ADMM solve the problem
% 
% minimize -logdet(X1) + Tr(SX1) 
% subject to  X >= 0, 0 <= X4 <= alpha*I, P(X2) = I
% 
% using ADMM algorithm based on the paper
% A. Pruttiakaravanich and J. Songsiri, "Convex Formulation for Regularized Estimation
% of Structural Equation Models"
%
% USAGE:
% [SOL, H] = confirmatory_sem_admm(S0,gamma,indA,alpha0,x0)
% 
%
% it requires input as follows.
%   1) S0 : sample covariance matrix
%   2) gamma : penalty parameter
%   3) indA : linear indices that A_ij is zero (index set I_A)
%   4) alpha0 (optional) : normaly we use alpha0 = min(eig(S))
%   5) X0 : initial solution (optional). 
%
% and it returns outputs: 'sol' as a structure variable with fields
% 
%   1) SOL.X : a symmtric X = [X1 X2 ; X2 X4] each Xk is n x n
%   2) SOL.A : our estimated path matrix
%   3) H: historical algorithm values: objective values, iterations, ...
%      objective values, primal and dual residual norms
% 


% We implement ADMM based on the global consensus problem, page 50 of 
% Stephen Boyd book on Distributed Optimization 
% 
% minimize f1(X) + f2(Z)
% subject to X - Z = 0
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
% note: it requires 'det_rootn' function by CVX to compute log_det of X > 0
% 
%% confirmatory SEM

beta = 1/alpha0; % scaling factor
alpha = 1; % used in scaled problem
S = beta*S0; % use this S to solve the scaled problem
m = 2;  rho = 1*beta; % heuristic choice
n = size(S0,1);
allind = (1:n^2)'; 
indnotA = setdiff(allind,indA,'rows');

% numerical parameters
MAXITER = 100000; FREQ_PRINT = 1000; PRINT_RESULT = 1;  
E_abs = 1e-5;               %absolute tolerance
E_rel = 1e-5;               %relative tolerance

X = zeros(2*n,2*n); Z = zeros(2*n,2*n); Y = zeros(2*n,2*n); 

if (PRINT_RESULT == 1)
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
    'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end
            
for ii=1:MAXITER
    Zold = Z;
    X = proxlogdet(Z-Y/rho,S,rho,alpha,indA);
    Z = proxpdf(X+Y/rho);
    Y = Y + rho*(X-Z);
      
    % cost objective value of scaled problem
    history.objval1(ii) = -n*log(det_rootn(X(1:n,1:n)))+trace(S*X(1:n,1:n));
    
    W = X; % W is the solution of original problem
    W(1:n,1:n) = X(1:n,1:n)*beta; W(n+1:end,n+1:end) = X(n+1:end,n+1:end)/beta;
    
    % cost objective value of original problem
    history.objval(ii) = -n*log(det_rootn(W(1:n,1:n)))+trace(S0*W(1:n,1:n));
    
    history.eps_pri(ii) = sqrt((2*n)^2)*E_abs + ...
        E_rel*max([norm(X(:)) norm(Z(:))]) ; 
    history.eps_dual(ii) = sqrt((2*n)^2)*E_abs + ...
        E_rel*norm(vec(Y));

    history.r_norm(ii) = norm( vec(X-Z));
    history.s_norm(ii) = norm(Z-Zold);

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

sol.X = X; 
sol.X(1:n,1:n) = X(1:n,1:n)*beta;
sol.X(n+1:end,n+1:end) = X(n+1:end,n+1:end)/beta;

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