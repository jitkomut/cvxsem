function [] = plot_performance(textfile,varargin)
%to plot performance of algorithm, all files that need to be loaded must be
%in the same folder of this function.
%
%example :
% 
% clear; clc; close all;
% 
%
% plot_performance('data_ppxa.txt',1);


% textfile = 'datarho1.txt'
% textfile = 'dataepsminus6.txt'
addpath('../');

    
fileID = fopen(textfile);
filename = textscan(fileID,'%s');

num_file = size(filename{1},1);
trial = 50;
x = zeros(num_file,1);
avgtime = zeros(2,2,num_file);
avgiter = zeros(2,2,num_file);
sdtime = zeros(2,2,num_file);
sditer = zeros(2,2,num_file);
avgmineig = zeros(2,num_file);
sdmineig = zeros(2,num_file);


for ii=1:num_file

        load(filename{1}{ii});
        
        fprintf('loading %s \n',filename{1}{ii});
         
        avgmineig(:,ii) = mean(result.eigen(:,1,:),3);
        sdmineig(:,ii) = std(result.eigen(:,1,:),0,3);
        
        x(ii) = result.n;
        
        if find(result.iter >= 5000)
            display('some trial hit maximum number of iterations');
        end
        
        % these measures are 2 x 2 matrices.

        avgtime(:,:,ii) = mean((result.time),3);    %mean of CPU time (seconds)
        avgiter(:,:,ii) = ceil(mean(result.iter,3));
        
        sdtime(:,:,ii) = std((result.time),0,3);    %sd of CPU time (seconds)
        sditer(:,:,ii) = ceil(std(result.iter,0,3));
end

%%

%     figure1 : plot CPU time when varying gamma
    figure;
    errorbar(x,squeeze(avgtime(1,1,:)),squeeze(sdtime(1,1,:)),'--x','MarkerSize',8,...
        'MarkerEdgeColor','red','MarkerFaceColor','red','linewidth',1.5)
    hold on;
    errorbar(x,squeeze(avgtime(1,2,:)),squeeze(sdtime(1,2,:)),'-.o','MarkerSize',5,...
        'MarkerEdgeColor','blue','MarkerFaceColor','blue','linewidth',1.5)
    hold off;
    lgd = legend('\gamma = 0.05\gamma_{max}','\gamma = 0.8\gamma_{max}');
    xlabel('n','FontSize',18);
    ylabel('CPU time (seconds)','FontSize',18);
    set(lgd,'FontSize',16,'Location','SouthEast');
  	set(gca,'XTick',x);
    set(gca,'XTickLabel',x);  
    
    
    %figure2 : plot CPU time when varying S0 factor
    figure;
    errorbar(x,squeeze(avgtime(1,1,:)),squeeze(sdtime(1,1,:)),'--x','MarkerSize',8,...
        'MarkerEdgeColor','red','MarkerFaceColor','red','linewidth',1.5)
    hold on;
    errorbar(x,squeeze(avgtime(2,1,:)),squeeze(sdtime(2,1,:)),'-.o','MarkerSize',5,...
        'MarkerEdgeColor','blue','MarkerFaceColor','blue','linewidth',1.5)
    hold off;
    lgd = legend('small \lambda_{min}(S)','moderate \lambda_{min}(S)');
    xlabel('n','FontSize',18);
    ylabel('CPU time seconds)','FontSize',18);
    set(lgd,'FontSize',16,'Location','SouthEast');
    set(gca,'XTick',x);
    set(gca,'XTickLabel',x);    
    
    %figure3 : plot number of iterations when varying gamma
    figure;
    errorbar(x,squeeze(avgiter(1,1,:)),squeeze(sditer(1,1,:)),'--x','MarkerSize',8,...
        'MarkerEdgeColor','red','MarkerFaceColor','red','linewidth',1.5)
    hold on;
    errorbar(x,squeeze(avgiter(1,2,:)),squeeze(sditer(1,2,:)),'-.o','MarkerSize',5,...
        'MarkerEdgeColor','blue','MarkerFaceColor','blue','linewidth',1.5)
    hold off;
    lgd = legend('\gamma = 0.05\gamma_{max}','\gamma = 0.8\gamma_{max}');
    xlabel('n','FontSize',18);
    ylabel('number of iterations','FontSize',18);
    set(lgd,'FontSize',16);
  	set(gca,'XTick',x);
    set(gca,'XTickLabel',x);  
    
    
    %figure4 : plot number of iterations when varying S0 factor
    figure;
    errorbar(x,squeeze(avgiter(1,1,:)),squeeze(sditer(1,1,:)),'--x','MarkerSize',8,...
        'MarkerEdgeColor','red','MarkerFaceColor','red','linewidth',1.5)
    hold on;
    errorbar(x,squeeze(avgiter(2,1,:)),squeeze(sditer(2,1,:)),'-.o','MarkerSize',5,...
        'MarkerEdgeColor','blue','MarkerFaceColor','blue','linewidth',1.5)
    hold off
    lgd = legend('small \lambda_{min}(S)','moderate \lambda_{min}(S)');
    xlabel('n','FontSize',18);
    ylabel('number of iterations','FontSize',18);
    set(lgd,'FontSize',16);
    set(gca,'XTick',x);
    set(gca,'XTickLabel',x);   
    
    
    fprintf('DONE!! PLOTTING \n');   
