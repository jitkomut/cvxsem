clear; clc;
addpath('../regsem_experiment')
addpath('../regsem_experiment/result_from_matlab')
addpath('../regsem_experiment/result_from_regsem')


%% Plot graph

load result2plot_regsemN100

% figure;
plot(FPRregsem,TPRregsem,'-+b',FPRsem,TPRsem,'-or','MarkerSize',12,'Linewidth',3)
xlim([0 1]); ylim([0 1]);axis square
xlabel('FPR');ylabel('TPR');title('ROC');
legend('regsem','sparse SEM','location','southeast');
ax = gca; ax.FontSize = 22; 
% set(gcf, 'Position', get(0, 'Screensize'));
ax = gca; [ax] = nowhitespace(ax);

print -depsc compare_result_regsem_sparsesem_N100.eps
savefig('compare_result_regsem_sparsesem_N100')

%% Plotgraph N=1000
load result2plot_regsemN1000

% figure;
plot(FPRregsem,TPRregsem,'-+b',FPRsem,TPRsem,'-or','MarkerSize',12,'Linewidth',3)
xlim([0 1]); ylim([0 1]);axis square
xlabel('FPR');ylabel('TPR');title('ROC');
legend('regsem','sparse SEM','location','southeast');
ax = gca; ax.FontSize = 22; 
% set(gcf, 'Position', get(0, 'Screensize'));
ax = gca; [ax] = nowhitespace(ax);

print -depsc compare_result_regsem_sparsesem_N1000.eps
savefig('compare_result_regsem_sparsesem_N1000')