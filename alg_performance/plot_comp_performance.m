function [avgtime1,sdtime1,avgiter1,sditer2] = plot_comp_performance(textfile1,textfile2)
%to plot comparative performance of two algorithms, all files that need to be loaded must be
%in the same folder of this function.
%
%example :
% 
% clear; clc; close all;
% 
%
% plot_performance('data_ppxa.txt','data_admm.txt');

addpath('../');

    
fileID1 = fopen(textfile1);
filename1 = textscan(fileID1,'%s');

fileID2 = fopen(textfile2);
filename2 = textscan(fileID2,'%s');

if size(filename1{1},1) ~= size(filename2{1},1)
    error('data from two algorithms do not have the same size');
end

num_file= size(filename1{1},1);


trial = 50;

% we compare two algorithms under the cases of small lambda_min(S) 
% we make plots of CPU times and Number of iterations

x = zeros(num_file,1);


avgtime1 = zeros(2,num_file); % vary gamma x num_file
avgiter1 = zeros(2,num_file); % vary gamma x num_file
avgtime2 = zeros(2,num_file); % vary gamma x num_file
avgiter2 = zeros(2,num_file); % vary gamma x num_file

sdtime1 = zeros(2,num_file); % vary gamma x TWO ALGs x num_file
sditer1 = zeros(2,num_file); % vary gamma x TWO ALGs x num_file
sdtime2 = zeros(2,num_file); % vary gamma x TWO ALGs x num_file
sditer2 = zeros(2,num_file); % vary gamma x TWO ALGs x num_file


for ii=1:num_file

        load(filename1{1}{ii});
        fprintf('loading %s \n',filename1{1}{ii});
        
        x(ii) = result.n;
        if find(result.iter >= 5000)
            display('some trial hit maximum number of iterations');
        end
                
        % these measures are 2 x 2 matrices.
%         after squeeze it's 2 x num_trials

        avgtime1(:,ii) = mean(squeeze(result.time(1,:,:)),2);    %mean of CPU time (seconds)
        avgiter1(:,ii) = ceil(mean(squeeze(result.iter(1,:,:)),2));
        
        sdtime1(:,ii) = std(squeeze(result.time(1,:,:)),0,2);    %sd of CPU time (seconds)
        sditer1(:,ii) = ceil(std(squeeze(result.iter(1,:,:)),0,2));
        
        load(filename2{1}{ii});
        fprintf('loading %s \n',filename2{1}{ii});
                     
        avgtime2(:,ii) = mean(squeeze(result.time(1,:,:)),2);    %mean of CPU time (seconds)
        avgiter2(:,ii) = ceil(mean(squeeze(result.iter(1,:,:)),2));
        
        sdtime2(:,ii) = std(squeeze(result.time(1,:,:)),0,2);    %sd of CPU time (seconds)
        sditer2(:,ii) = ceil(std(squeeze(result.iter(1,:,:)),0,2));

end

%%

%     figure1 : plot CPU time when gamma = 0.05gamma_max
    figure;
    subplot(221);
    errorbar(x,avgtime1(1,:),sdtime1(1,:),'--x','MarkerSize',8,...
        'MarkerEdgeColor','red','MarkerFaceColor','red','linewidth',1.5);
    hold on;
    
    errorbar(x,avgtime2(1,:),sdtime2(1,:),'-.o','MarkerSize',8,...
        'MarkerEdgeColor','blue','MarkerFaceColor','blue','linewidth',1.5);
    hold off;
    lgd = legend('ADMM','PPXA'); title('\gamma = 0.05\gamma_{max}','FontSize',18);
    xlabel('n','FontSize',18);
    ylabel('CPU time (seconds)','FontSize',18);
    set(lgd,'FontSize',16,'Location','SouthEast');
  	set(gca,'XTick',x);
    set(gca,'XTickLabel',x);  
    
    subplot(222);
    errorbar(x,avgtime1(2,:),sdtime1(2,:),'--x','MarkerSize',8,...
        'MarkerEdgeColor','red','MarkerFaceColor','red','linewidth',1.5);
    hold on;
    
    errorbar(x,avgtime2(2,:),sdtime2(2,:),'-.o','MarkerSize',8,...
        'MarkerEdgeColor','blue','MarkerFaceColor','blue','linewidth',1.5);
    hold off;
    lgd = legend('ADMM','PPXA'); title('\gamma = 0.8\gamma_{max}','FontSize',18);
    xlabel('n','FontSize',18);
    ylabel('CPU time (seconds)','FontSize',18);
    set(lgd,'FontSize',16,'Location','SouthEast');
  	set(gca,'XTick',x);
    set(gca,'XTickLabel',x);  
    
    
    %figure2 : plot average number of iteration
%     figure;
    subplot(223);
    errorbar(x,avgiter1(1,:),sditer1(1,:),'--x','MarkerSize',8,...
        'MarkerEdgeColor','red','MarkerFaceColor','red','linewidth',1.5);
    hold on;
    
    errorbar(x,avgiter2(1,:),sditer2(1,:),'-.o','MarkerSize',8,...
        'MarkerEdgeColor','blue','MarkerFaceColor','blue','linewidth',1.5);
    hold off;

    lgd = legend('ADMM','PPXA'); title('\gamma = 0.05\gamma_{max}','FontSize',18);
    xlabel('n','FontSize',18);
    ylabel('Number of iterations','FontSize',18);
    set(lgd,'FontSize',16,'Location','SouthEast');
  	set(gca,'XTick',x);
    set(gca,'XTickLabel',x); 
    
    subplot(224);
    errorbar(x,avgiter1(2,:),sditer1(2,:),'--x','MarkerSize',8,...
        'MarkerEdgeColor','red','MarkerFaceColor','red','linewidth',1.5);
    hold on;
    
    errorbar(x,avgiter2(2,:),sditer2(2,:),'-.o','MarkerSize',8,...
        'MarkerEdgeColor','blue','MarkerFaceColor','blue','linewidth',1.5);
    hold off;

    lgd = legend('ADMM','PPXA');  title('\gamma = 0.8\gamma_{max}','FontSize',18);
    xlabel('n','FontSize',18);
    ylabel('Number of iterations','FontSize',18);
    set(lgd,'FontSize',16,'Location','SouthEast');
  	set(gca,'XTick',x);
    set(gca,'XTickLabel',x); 
    
    
    %% Plot only this case (eigen of S is small and gamma = 0.05)
        figure;
    subplot(121);
    errorbar(x,avgtime1(1,:),sdtime1(1,:),'--x','MarkerSize',8,...
        'MarkerEdgeColor','red','MarkerFaceColor','red','linewidth',1.5);
    hold on;
    
    errorbar(x,avgtime2(1,:),sdtime2(1,:),'-.o','MarkerSize',8,...
        'MarkerEdgeColor','blue','MarkerFaceColor','blue','linewidth',1.5);
    hold off;
    lgd = legend('ADMM','PPXA'); 
    %title('\gamma = 0.05\gamma_{max}','FontSize',16);
    xlabel('n','FontSize',14);
    ylabel('CPU time (seconds)','FontSize',14);
    set(lgd,'FontSize',12,'Location','SouthEast');
  	set(gca,'XTick',x);
    set(gca,'XTickLabel',x);  
    
    subplot(122);
        errorbar(x,avgiter1(1,:),sditer1(1,:),'--x','MarkerSize',8,...
        'MarkerEdgeColor','red','MarkerFaceColor','red','linewidth',1.5);
    hold on;
    
    errorbar(x,avgiter2(1,:),sditer2(1,:),'-.o','MarkerSize',8,...
        'MarkerEdgeColor','blue','MarkerFaceColor','blue','linewidth',1.5);
    hold off;

    lgd = legend('ADMM','PPXA'); 
    %title('\gamma = 0.05\gamma_{max}','FontSize',18);
    xlabel('n','FontSize',14);
    ylabel('Number of iterations','FontSize',14);
    set(lgd,'FontSize',12,'Location','SouthEast');
  	set(gca,'XTick',x);
    set(gca,'XTickLabel',x); 
    
        
