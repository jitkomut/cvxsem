% This file is to prepare the result of experiment described in section
% 6.2, based on 
% 
% A. Pruttiakaravanich and J. Songsiri, "Convex Formulation for Regularized Estimation
% of Structural Equation Models"
% 
% In this experiment, we explore a common brain network between control and 
% autism groups under resting-state condition learned from 
% Autism Brain Imaging Data Exchange (ABIDE) dataset.
% The data used in the experiment has been preprocessed from Preprocessed 
% Connectomes Project (PCP). See more information in a below link
% 
% http://preprocessed-connectomes-project.org/abide/
% 
% The preprocessing step has been done by Configurable Pipeline for the Analysis 
% of Connectomes (CPAC). In the detail, this pipeline performs a structural
% preprocessing of skull-stripping using AFNIs 3dSkullStrip, three-issue type brain segmentation
% using FAST in FSL, skull-stripped brain normalization to MNI152 with linear and non-linear
% registrations using ANTs. Then it performs a slice timing correction and motion realignment,
% respectively. The image intensity was normalized by 4D global mean and a band-pass filtering
% in the range of 0.01-0.1 Hz is applied. All images from every subject are transformed from the
% original to MNI152 template. To reduce the data dimension, we average the time series over a
% region of interest (ROI) using Automated Anatomical Labeling (AAL) template.
% 
% The preprocessed data of two group, AUTISM and CONTROL, is located abide_raw_data_input folder.
% Each data file provides a struct variable with field
%   1) data         : an observed matrix with size 296x116 (296 time points x 116 ROIs)
%   2) colheaders   : the label of each ROI (or node)
% 
% During estimation process, we reduce the number of ROIs from 116 to 90 due to
% the limitation of time point. In conclusion, we use our method to expore 
% the relationship among 90 ROIs of AUTISM and CONTROL group.
% 
% user can find more info about the label of each ROI (or node) from file "Node_AAL90.node". 
% This file provides the detail as follows.
%   column 1)   :   node's position on X-AXIS
%   column 2)   :   node's position on Y-AXIS
%   column 3)   :   node's position on Z-AXIS
%   column 4)   :   node's color's code
%   column 5)   :   node's size
%   column 6)   :   node's name
% 
% NOTE: 
%   1) The ROI's names in "Node_AAL90.node" have been sorted corresponding to
% colheaders already. For example, 
%       #2001	=	PreCG.L
%       #2002	=	PreCG.R
%       #2101	=	SFGdor.L
%    	#2102	=	SFGdor.R
% 
%   2) The file "Node_AAL90.node" is also used for plotting the brain
% network using BrainNetViewer.
% 
% https://www.nitrc.org/projects/bnv/
% 
% Developer: Anupon Pruttiakaravanich, Dr. Jitkomut Songsiri
% 
%% experiment for autism group
clear; clc; close all;
addpath ..\codes\abide_raw_data_input\autism-1
num_candidate = 25;         %number of candidate model (#point of gamma along x-axis)

 
for k=272:326
    filename_open = sprintf('UM_1_0050%d_rois_aal.1D',k)
    file = importdata(filename_open);

    Y = file.data;
    n = 116;            	%# of variables (ROIs)
    
    %standardize Y                       
    sd = std(Y,0,1);        %sd
    md = mean(Y);           %mean
    for i=1:n
    Y(:,i) = (Y(:,i)-md(i))/sd(i);
    end
    
    %reduce 116 ROIs to 90 ROIs [we use 1st version of AAL atlas due to the limit of number of observations (N)]
    Y(:,91:end)=[];   
 	n = 90;                 %# of variables is now reduced to 90
    indA = (1:n+1:n^2)';  	%index of A where P(A) = 0
                            %we do not have any prior assumptions for this experiment
                            
    [output] = modelsel_sem(num_candidate,Y',indA);

    filename_save = sprintf('result_UM_1_0050%d_rois_aal',k);
%     save (filename_save)
end

%% experiment for control group
clear; clc; close all;
addpath ..\codes\abide_raw_data_input\control-2
num_candidate = 25;         %number of candidate model (#point of gamma along x-axis)

 
for k=329:381
    filename_open = sprintf('UM_1_0050%d_rois_aal.1D',k);
    file = importdata(filename_open);

    Y = file.data;
    n = 116;            	%# of variables (ROIs)
    
    %standardize Y                       
    sd = std(Y,0,1);        %sd
    md = mean(Y);           %mean
    for i=1:n
    Y(:,i) = (Y(:,i)-md(i))/sd(i);
    end
    
    %reduce 116 ROIs to 90 ROIs [we use 1st version of AAL atlas due to the limit of number of observations (N)]
    Y(:,91:end)=[];   
 	n = 90;                 %# of variables is now reduced to 90
    indA = (1:n+1:n^2)';  	%index of A where P(A) = 0
                            %we do not have any prior assumptions for this experiment
                            
    [output] = modelsel_sem(num_candidate,Y',indA);

    filename_save = sprintf('result_UM_1_0050%d_rois_aal',k);
%     save (filename_save)
end