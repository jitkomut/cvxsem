% This file is a code for expenriment described is in section 5.2, 
% comparison with an existing method, based on,
% 
% A. Pruttiakaravanich and J. Songsiri, "Convex Formulation for Regularized Estimation
% of Structural Equation Models"
% 
% In this experiment, we compare our result with the result from REGSEM.
% Regsem package in R was developed to solve structural equation modeling 
% with regularization term including both ridge and the least absolute 
% shrinkage and selection operator (lasso).
% 
% Firstly, we generate Atrue corresponding to the structure described in
% section 5.2, and then generate the obaserved matrix by Y = (inv(I -
% Atrue))*e, where e is normally distributed with
% variance of 0.1. In the data generating process, we vary the number of
% observation between N=100 and N=1000. We repeat this process to have
% 100 sets of Y (100 trials).
% 
% We use REGSEM package in R to estimate the path matrix. We then import
% results of REGSEM to matlab, compute averaged TP,TN,FP,FN over 100 trials
% and finally plot ROC curve.
% 
% Developer: Anupon Pruttiakaravanich, Dr. Jitkomut Songsiri
% 
%% comparing with REGSEM
clear; clc;
addpath('../regsem_experiment')
addpath('../regsem_experiment/result_from_matlab')
addpath('../regsem_experiment/result_from_regsem')

load Atrue.mat
load regsem_N1000_lin.mat

q = 2;    r=4;    s=2;    t=4;
n = q+r+s+t;
no_of_estimated_para = (q*r) + (r*s) + (s*t);
num_trial = 100;
num_candidate = 50;
index = [0 logspace(-5, 0, num_candidate - 1)];

%for calculating TP, TN, FP, FN
Htrue = Atrue(1:2, 3:6);
Jtrue = Atrue(3:6, 7:8);
Ktrue = Atrue(7:8, 9:12);

% reindexing the result from REGSEM
for j=1:num_trial
    for i=1:num_candidate
    myA = zeros(12,12);
    myA(1,3) = R(1,i,j);
    myA(1,4) = R(2,i,j); 
    myA(1,5) = R(3,i,j);
    myA(1,6) = R(4,i,j);
    myA(2,3) = R(5,i,j);
    myA(2,4) = R(6,i,j);
    myA(2,5) = R(7,i,j);
    myA(2,6) = R(8,i,j);
    myA(3,7) = R(9,i,j);
    myA(3,8) = R(10,i,j);
    myA(4,7) = R(11,i,j);
    myA(4,8) = R(12,i,j);
    myA(5,7) = R(13,i,j);
    myA(5,8) = R(14,i,j);
    myA(6,7) = R(15,i,j);
    myA(6,8) = R(16,i,j);
    myA(7,9) = R(17,i,j);
    myA(7,10) = R(18,i,j);
    myA(7,11) = R(19,i,j);
    myA(7,12) = R(20,i,j);
    myA(8,9) = R(21,i,j);
    myA(8,10) = R(22,i,j);
    myA(8,11) = R(23,i,j);
    myA(8,12) = R(24,i,j);
  
    
    %calculate TP, TN, FP, FN
    Hhat = myA(1:2, 3:6);
    Jhat = myA(3:6, 7:8);
    Khat = myA(7:8, 9:12);
    
    [TP_H,TN_H,FP_H,FN_H] = comp_error_input_not_symmetric(Htrue,Hhat);
    [TP_J,TN_J,FP_J,FN_J] = comp_error_input_not_symmetric(Jtrue,Jhat);
    [TP_K,TN_K,FP_K,FN_K] = comp_error_input_not_symmetric(Ktrue,Khat);
    TP(j,i) = TP_H + TP_J + TP_K;
    TN(j,i) = TN_H + TN_J + TN_K;
    FP(j,i) = FP_H + FP_J + FP_K;
    FN(j,i) = FN_H + FN_J + FN_K;
    
    fprintf('computing TP,TN,FP,FN of result from REGSEM: trial no.%d, candidate model no.%d\n',j,i);
    end
end

%plot ROC of result from REGSEM
TPmean = mean(TP);
TNmean = mean(TN);
FPmean = mean(FP);
FNmean = mean(FN);

TPR = TPmean./(TPmean + FNmean);
FPR = FPmean./(FPmean + TNmean);
% figure;
plot(FPR,TPR,'-o')
axis square;
xlim([0 1])
ylim([0 1])
set(gca,'FontSize',15)
hold on;

%
clear; clc;
load rec_path_ex_1varnoise_N1000_linscale_100trial.mat

num_trial = 100;
num_candidate = 50;

%for calculating TP, TN, FP, FN
Htrue = Atrue(1:2, 3:6);
Jtrue = Atrue(3:6, 7:8);
Ktrue = Atrue(7:8, 9:12);

for j=1:num_trial
    for i=1:num_candidate
    
    %calculate TP, TN, FP, FN
    Hhat = A_saved{j,i}(1:2, 3:6);
    Jhat = A_saved{j,i}(3:6, 7:8);
    Khat = A_saved{j,i}(7:8, 9:12);
    
    [TP_H,TN_H,FP_H,FN_H] = comp_error_input_not_symmetric(Htrue,Hhat);
    [TP_J,TN_J,FP_J,FN_J] = comp_error_input_not_symmetric(Jtrue,Jhat);
    [TP_K,TN_K,FP_K,FN_K] = comp_error_input_not_symmetric(Ktrue,Khat);
    TP(j,i) = TP_H + TP_J + TP_K;
    TN(j,i) = TN_H + TN_J + TN_K;
    FP(j,i) = FP_H + FP_J + FP_K;
    FN(j,i) = FN_H + FN_J + FN_K;
    fprintf('computing TP,TN,FP,FN of result from MATLAB: trial no.%d, candidate model no.%d\n',j,i);
    end
end

% plot ROC of result from MATLAB
TPmean = mean(TP);
TNmean = mean(TN);
FPmean = mean(FP);
FNmean = mean(FN);

TPR = TPmean./(TPmean + FNmean);
FPR = FPmean./(FPmean + TNmean);
% figure;
plot(FPR,TPR,'-o')
axis square;
xlim([0 1])
ylim([0 1])
legend('regsem','sparse sem')
xlabel('FPR')
ylabel('TPR')
title('ROC')
set(gca,'FontSize',15)
hold off;

