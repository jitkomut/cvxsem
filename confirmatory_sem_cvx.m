%% confirmatory SEM by CVX

function[X1,cvx_optval] = confirmatory_sem_cvx(S,alpha,indA)

n = size(S,1);
I = eye(n);
cvx_begin sdp
    variable X1(2*n,2*n) symmetric
    variable X2(n,n)
    minimize -log_det(X1(1:n,1:n))+trace(S*X1(1:n,1:n))
        subject to
            X1 >= 0;
            X2 == X1(n+1:end,1:n);
            X2(indA) == I(indA);
            X1(n+1:end,n+1:end) <= alpha*eye(n);
cvx_end
end