% This file is to prepare the result of experiment described in section
% 5.1, based on 
% 
% A. Pruttiakaravanich and J. Songsiri, "Convex Formulation for Regularized Estimation
% of Structural Equation Models"
% 
% In this experiment, we explore the accuracy of our sparse sem in
% detecting nonzero pattern of a true model. We would like to observe what
% factor can affect to the performance of our formulation.
% 
% Firstly, we set n = 10 and choose Atrue having random sparsity patterns 
% and generate measurements from Y = (inv(I - Atrue))*e where e is normally distributed with
% variance of 0.1. The matrix S is then computed as the sample covariance of Y.
% 
% In data generating process, we vary 2 factors as
% 1) the number of samples (N)
% 2) the density of non-zero elements in the true model (Atrue)
% and in estimation process, we also vary the condition of known zero
% locations as required in P(A) = 0.
% 
% To examine the performance,
% we apply typical measures: True positives (TP), True negative (TN), False positives (FP), False
% negatives (FN), TP rate (TPR), FP rate (FPR), where positives are non-zero entries in A and
% negatives are zeros in A. These measeure are computed for each model
% (value of gamma)
% 
% Finally, the accuracy (or the error) is computed from the model (the path matrix, A)
% which has been selected from each model selection criterion. This
% accuracy are averaged over 40 trials.
%   
% Developer: Anupon Pruttiakaravanich, Dr. Jitkomut Songsiri
% 
%% This file is a code for experment 5.1 (performance of sem)

clear;  clc;

%data generating
n = 10;                     %number of observed variables
num_trial = 50;           	%number of trial;
density = 0.2;           	%density of nonzero entries in Atrue (Atrue is sparse)
num_candidate = 40;         %number of candidate model (#point of gamma along x-axis)
N = 100;                    %number of observation (samples)

%running with various number of trials
for k=1:num_trial
 
    input_data = gen_datasem(n,N,density);
    
    Atrue = input_data.Atrue;
    indz_Atrue = find(Atrue == 0);
    indz_Atrue_off_diagonal = setdiff(indz_Atrue,(1:n+1:n^2));

    S = input_data.S;
    alpha = min(eig(S));
   
    %no assumption about zero index
    %indexing of Ahat from assumption
    inddiag = find(eye(n));    
    indA0 = inddiag;           %index of I_A for a_ij = 0; (no assumption)
    allind = (1:n^2)'; 
    indnotA0 = setdiff(allind,indA0,'rows');
    
    %[USING THE BELOW COMMENTED LINES WHEN WE ASSUME TO KNOW TRUE ZERO LOCATIONS IN ATRUE]
    %redefine the indexing of Ahat from assumption (known location of zero)
%     known_percentage = 0.2;             %know zero location 20%
%     indices = randperm(length(indz_Atrue_off_diagonal),round(known_percentage*length(indz_Atrue_off_diagonal)));
%     indassump = indAtrue(indices);
%     inddiag = find(eye(n));    
%     indA0 = [indassump';inddiag];        %index of I_A for a_ij = 0; (has assumption)
%     allind = (1:n^2)'; 
%     indnotA0 = setdiff(allind,indA0,'rows'); 
    
    %gamma
    tmpS = alpha*eye(n)-S;
    gamma_max = (1/alpha)*max(abs(tmpS(indnotA0)));

    %mode : gamma with log scale (40 points)
	gamma = (gamma_max)*[0 logspace(-4,0,num_candidate-1)];
    
    for j = 1:num_candidate
        if (j == 1)
            %solve convex confirm SEM
            indA = indA0;
            [sol_confirm, hist] = confirmatory_sem_ppxa(S,indA,alpha,'tolfun',1e-5,'tolx',1e-5);           
            Acon(:,:,j) = sol_confirm.A;            
            fprintf('DONE solving convex confirmatory sem: candidate model(%d) of trials(%d)\n',j,k);
        else
            %solve sparse SEM
            [sol_sparse,~] = sparse_sem_ppxa(S,gamma(j),indA0,alpha,'tolfun',1e-5,'tolx',1e-5);    
            indA = find(sol_sparse.A == 0);             	%#zero entries in A will be used in solving confirm SEM
            fprintf('DONE solving convex sparse sem: candidate model(%d) of trials(%d)\n',j,k);

            %solve convex confirm SEM
            [sol_confirm, hist] = confirmatory_sem_ppxa(S,indA,alpha,'tolfun',1e-5,'tolx',1e-5);
            Acon(:,:,j) = sol_confirm.A;
            
            fprintf('DONE solving convex confirmatory sem: candidate model(%d)\n',k);
        end 

        %compute BIC scores of each model
        %#effective parameter (#overall parameters - #parameter in A which are set to be 0 from the constraint indA = 0)
        d = 2*n - 2*length(indA);                   %number of effective parameters
        minus2L = N*(hist.objval(end));             %-2L

        BIC(j) = minus2L + d*log(N);
        AIC(j) = minus2L + 2*d;
        AICC(j) = minus2L + (2*d*N)/(N-d-1);
        KIC(j) = minus2L + 3*d;
        KICC(j) = minus2L + ((d+1)*((3*N)-d-2)/(N-d-2)) + (d/(N-d));  
    end
    
    [BICmin, ind_BICmin] = min(BIC);
    A_bic = Acon(:,:,ind_BICmin);

    [AICmin, ind_AICmin] = min(AIC);
    A_aic = Acon(:,:,ind_AICmin);

    [AICCmin, ind_AICCmin] = min(AICC);
    A_aicc = Acon(:,:,ind_AICCmin);

    [KICmin, ind_KICmin] = min(KIC);
    A_kic = Acon(:,:,ind_KICmin);

    [KICCmin, ind_KICCmin] = min(KICC);  
    A_kicc = Acon(:,:,ind_KICCmin);
    
    %computing error of each criterions
    [TP_bic,TN_bic,FP_bic,FN_bic] = comp_error(Atrue,Acon(:,:,ind_BICmin),n);
    [TP_aic,TN_aic,FP_aic,FN_aic] = comp_error(Atrue,Acon(:,:,ind_AICmin),n);
   	[TP_aicc,TN_aicc,FP_aicc,FN_aicc] = comp_error(Atrue,Acon(:,:,ind_AICCmin),n);
	[TP_kic,TN_kic,FP_kic,FN_kic] = comp_error(Atrue,Acon(:,:,ind_KICmin),n);
   	[TP_kicc,TN_kicc,FP_kicc,FN_kicc] = comp_error(Atrue,Acon(:,:,ind_KICCmin),n);

    %save output and result
    output{k} = struct('S',S,'Atrue',Atrue,...
                        'A_bic',A_bic,'A_aic',A_aic,'A_aicc',A_aicc,'A_kic',A_kic,'A_kicc',A_kicc,...
                        'TP_bic',TP_bic,'TN_bic',TN_bic,'FP_bic',FP_bic,'FN_bic',FN_bic,...
                        'TP_aic',TP_aic,'TN_aic',TN_aic,'FP_aic',FP_aic,'FN_aic',FN_aic,...
                        'TP_aicc',TP_aicc,'TN_aicc',TN_aicc,'FP_aicc',FP_aicc,'FN_aicc',FN_aicc,...
                        'TP_kic',TP_kic,'TN_kic',TN_kic,'FP_kic',FP_kic,'FN_kic',FN_kic,...
                        'TP_kicc',TP_kicc,'TN_kicc',TN_kicc,'FP_kicc',FP_kicc,'FN_kicc',FN_kicc);
    
    fprintf('-----------------------------\n');
    fprintf('DONE trials(%d)\n',k);
    fprintf('-----------------------------\n');
    
%     save('experiment51_result.mat','output')
end


