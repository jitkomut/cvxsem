%% Parallel Proximal Algorithm (PPXA) for sparse SEM
% the algorithm is based on PPXA page 203 in Combett book (ALG 10.27)
% STILL DIVERGE

function [sol, history] = sparse_sem_ppxa(S0,gamma,indA,alpha0,varargin)
% SPARSE_SEM_PPXA solve the problem
% 
% minimize -loget(X1) + Tr(SX1) + 2*gamma \sum_{(i,j) not in I_A} | (X2)_ij | 
% subject to X >= 0, 0 <= X4 <= alpha*I , P(X2) = I
% We implement PPXA based on Combett book (Proximal Algorithm)
% 
% minimize f1(X1) + f2(X2) + f3(X3) 
% C = {(X1,X2,X3) | X1= X2 = X3  } 
% where 
% f1(Y) = -logdet(Y1)+ Tr(SY1) + I{ 0 <= Y4 <= alpha*I } 
% f2(Y) = 2*gamma sum_{(i,j) not in indA} |(Y2)ij| + I{ P(Y2) = I} 
% f3(Y) = I{ Y >= 0} 
% 
% I{ X in C} denotes an indicator function)
% 
% We also solve the scaled problem using (S,alpha) = (beta*S0,1) 
% where S0 is the original sample covariance and beta = 1/alpha0
% the original solution is obtained by scaling back (check Proposition 2 in the
% paper)
% 
% note: it requires 'det_rootn' function by CVX to compute log_det of X > 0
% 
% 
%  Scaled problem
beta = 1/alpha0; % scaling factor
alpha = 1; % used in scaled problem
S = beta*S0; % use this S to solve the scaled problem

n = size(S0,1);
allind = (1:n^2)'; 
indnotA = setdiff(allind,indA,'rows');


% numerical parameters
MAXITER = 5000; FREQ_PRINT = 100; PRINT_RESULT = 1;
E_abs = 1e-5;               %absolute tolerance
E_rel = 1e-5;               %relative tolerance
TOL_RELF = 1e-8; % relative change of cost objective
TOL_RELX = 1e-6; % relative change of solution

m = 3;  
lambda = 1.8; rho = 0.1; % PPXA algorithm parameters (tuned)

% FOUR VARIABLES:  Y,P,Pbar (auxillary), X (main)

Y1 = S\eye(n); Y4 = alpha*eye(n); Y2 = zeros(n); Y0 = [Y1 Y2';Y2 Y4];
Y = repmat(Y0,1,1,m);
X = Y0; w = 1/m;
P = repmat(Y0,1,1,m); Pbar = Y0;
f = objval(X,S0,n,gamma,indnotA);


if (PRINT_RESULT == 1)
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
    'r norm', 's norm', 'rel change in f', 'rel change in x', 'objective');
end

for ii=1:MAXITER
    fold = f;
    Xold = X;
    Pold = Pbar;
    P(:,:,1) = proxlogdet(Y(:,:,1),S,w/rho,alpha); 
    P(:,:,2) = proxl1(Y(:,:,2),n,gamma*rho/w,indA,indnotA); 
    P(:,:,3) = proxpdf(Y(:,:,3)); 
    
    Pbar = mean(P,3);

    for k=1:m
        Y(:,:,k) = Y(:,:,k) + lambda*(2*Pbar - X - P(:,:,k));
    end
    X = X + lambda*(Pbar-X);
    
    % cost objective value of scaled problem
    history.objval1(ii) = objval(X,S,n,gamma,indnotA);
    
    W = X; % W is the solution of original problem
    W(1:n,1:n) = X(1:n,1:n)*beta; W(n+1:end,n+1:end) = X(n+1:end,n+1:end)/beta;
    
    history.objval(ii) = objval(W,S0,n,gamma,indnotA);
    
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

sol.X = X; 
% use the solution X2 (from soft thresholding to get exact zeros)
sol.X(n+1:end,1:n) = P(n+1:end,1:n,2); sol.X(1:n,n+1:end) = P(1:n,n+1:end,2);
sol.X(1:n,1:n) = X(1:n,1:n)*beta;
sol.X(n+1:end,n+1:end) = X(n+1:end,n+1:end)/beta;

sol.A = eye(n)-sol.X(n+1:end,1:n); % A = I-X2;

%% Cost objective
function[f] = objval(X,S,n,gamma,indnotA)

X2 = X(n+1:end,1:n);
% log_det function (by CVX) cannot be used within parfor
% f = -log_det(X(1:n,1:n))+trace(S*X(1:n,1:n))+2*gamma*sum(abs(X2(indnotA)));

% det function (built-in matlab) cannot be directly used as log(det(X))
f = -n*log(det_rootn(X(1:n,1:n)))+trace(S*X(1:n,1:n))+2*gamma*sum(abs(X2(indnotA))); 

if imag(f) ~= 0
    f = Inf; % if block X1 is not >= 0 then det < 0 and logdet is complex
end
end


%% Proximal operator of loget 
% f(X) = -logdet(X1) + Tr(S*X1) + I{ 0 <= X4 <= alpha*I } 
% minimize_X -logdet(X1) + Tr(SX1) + (rho/2) || Y - X ||^2 
% subject to X1 >= 0,  0 <= X4 <= alpha*I

function[Z] = proxlogdet(Y,S,rho,alpha)

n = size(S,1); 
Z = Y; % pre-assign 

[U,D] = eig(rho*Y(1:n,1:n)-S);
d = diag(D);
Dz1 = diag(0.5*(d+sqrt(d.^2+4*rho) )/rho);
Z(1:n,1:n) = U*Dz1*U'; % modify on Z1

[U4,D4] = eig(Y(n+1:end,n+1:end));
Dz4 = diag(max(0,min(alpha,diag(D4))));
Z(n+1:end,n+1:end) = U4*Dz4*U4'; % modify on Z4

end

%% Proximal operator of l1
% minimize_X  2a sum_{(i,j) not IA} | (X_2)ij | + (1/2) || Y - X ||^2
% subject to  P(X2) = I

function[Z] = proxl1(Y,n,a,indA,indnotA)

I = eye(n);
Z = Y;
Y2 = Y(n+1:end,1:n); Z2 = zeros(n);
Z2(indnotA) = max(abs(Y2(indnotA))-a,0).*sign(Y2(indnotA)); % factor 2g/2
Z2(indA) = I(indA); % P(Z2) = I
Z(n+1:end,1:n) = Z2; Z(1:n,n+1:end) = Z2';
end

%% Proximal operator of cone constraint
% minimize_{ X >= 0 } (1/2) || Y - X ||^2

function[Z] = proxpdf(Y)

[U,D] = eig(Y);
Dz = diag([max(0,diag(D))]);
Z = U*Dz*U';
end

end
