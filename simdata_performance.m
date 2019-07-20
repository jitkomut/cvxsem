
%% TOTAL ERROR Atrue DENSE
% totalerr = 5 methods x 3 cases of known zeros x 2 cases of N


totalerr(:,:,1) = [0.1669 0.1536 0.1889 0.1573 0.1802;...
                    0.1644 0.1518 0.1858 0.1562 0.1784;...
                    0.1602 0.1467 0.1804 0.1511 0.1729];

totalerr(:,:,2) = [0.1342 0.1331 0.1331 0.1324 0.1324;...
                    0.1307 0.1289 0.1289 0.1289 0.1289;
                     0.1260 0.1242 0.1242 0.1242 0.1242];

close all;                
[fig1] = plotbarchart(totalerr,'N=100','N=100,000',0:0.02:0.2,0:0.02:0.14);

print('-painters','-depsc','simulated_data_averaged_error_Atrue_dense')
savefig('simulated_data_averaged_error_Atrue_dense');

%% TOTAL ERROR Atrue SPARSE
% totalerr = 5 methods x 3 cases of known zeros x 2 cases of N

totalerr(:,:,1) = [0.5 0.5722 0.5322 0.5478 0.4744; ...
                    0.3944 0.4756 0.4289 0.4433 0.4011;...
                    0.2000 0.2333 0.1756 0.2267 0.1611];
                
totalerr(:,:,2) = [0.5622 0.7133 0.6944 0.6944 0.6944;...
                    0.4500 0.5500 0.5500 0.5378 0.5378;...
                    0.2544 0.2822 0.2822 0.2811 0.2811];

close all;                
[fig1] = plotbarchart(totalerr,'N=100','N=100,000',0:0.1:0.6,0:0.1:0.8);

print('-painters','-depsc','simulated_data_averaged_error_Atrue_sparse') 
savefig('simulated_data_averaged_error_Atrue_sparse');

%% FN Atrue DENSE
% FN = 5 methods x 3 cases of known zeros x 2 cases of N

FN(:,:,1) = [0.0573 0.0356 0.0876 0.0431 0.0773;...
               0.0582 0.0360 0.0882 0.0433 0.0778;...
               0.0568 0.0356 0.0887 0.0429 0.0778];
                
FN(:,:,2) = [0.0022 0.0002 0.0002 0.0002 0.0002;...
             0.0022 0.0002 0.0002 0.0002 0.0002;...
             0.0022 0.0002 0.0002 0.0002 0.0002]*(1e-3);

close all;                
[fig1] = plotbarchart(FN,'N=100','N=100,000',0:0.02:0.09,[]);

print('-painters','-depsc','simulated_data_averaged_FN_Atrue_dense') 
savefig('simulated_data_averaged_FN_Atrue_dense');

%% FN Atrue SPARSE
% FN = 5 methods x 3 cases of known zeros x 2 cases of N

FN(:,:,1) = [0.0256 0.0156 0.0378 0.0189 0.0433;...
             0.0244 0.0122 0.0256 0.0144 0.0400;...
             0.0178 0.0089 0.0556 0.0100 0.0578];
                
FN(:,:,2) = [0.0044 0.0022 0.0022 0.0022 0.0022;...
             0.0044 0.0011 0.0022 0.0022 0.0022;...
             0.0033 0.0011 0.0011 0.0011 0.0011]*(1e-3);

close all;                
[fig1] = plotbarchart(FN,'N=100',...
    'N=100,000',0:0.02:0.09,[]);

print('-painters','-depsc','simulated_data_averaged_FN_Atrue_sparse') 
savefig('simulated_data_averaged_FN_Atrue_sparse');

%% FP Atrue DENSE
% FP = 5 methods x 3 cases of known zeros x 2 cases of N

FP(:,:,1) = [0.1096 0.1180 0.1013 0.1142 0.1029;...
             0.1062 0.1158 0.0976 0.1129 0.1007;...
             0.1027 0.1111 0.0918 0.1082 0.0951];
                
FP(:,:,2) = [0.1320 0.1329 0.1329 0.1322 0.1322;...
             0.1284 0.1287 0.1287 0.1287 0.1287;...
             0.1238 0.1240 0.1240 0.1240 0.1240];

close all;                
[fig1] = plotbarchart(FP,'N=100',...
    'N=100,000',0:0.02:0.12,0:0.02:0.14);
% 
print('-painters','-depsc','simulated_data_averaged_FP_Atrue_dense') 
savefig('simulated_data_averaged_FP_Atrue_dense');


%% FP Atrue SPARSE
% FP = 5 methods x 3 cases of known zeros x 2 cases of N

FP(:,:,1) = [0.4744 0.5567 0.4944 0.5289 0.4311;...
             0.3700 0.4633 0.4033 0.4289 0.3611;...
             0.1822 0.2244 0.1200 0.2167 0.1033];
                
FP(:,:,2) = [0.5578 0.7111 0.7111 0.6922 0.6922;...
             0.4456 0.5489 0.5489 0.5356 0.5356;...
             0.2511 0.2811 0.2811 0.2800 0.2800];

close all;                
[fig1] = plotbarchart(FP,'N=100',...
    'N=100,000',0:0.1:0.6,0:0.1:0.8);
% 
print('-painters','-depsc','simulated_data_averaged_FP_Atrue_sparse') 
savefig('simulated_data_averaged_FP_Atrue_sparse');
