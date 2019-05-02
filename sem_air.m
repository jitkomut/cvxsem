% This file is to prepare the result of experiment described in section
% 6.1, based on 
% 
% A. Pruttiakaravanich and J. Songsiri, "Convex Formulation for Regularized Estimation
% of Structural Equation Models"
% 
% In this experiment, we explore a relation structure among 11 climate variables 
% including greenhouse gas (SO2, NO2, O3, CO), solar radiation, relative humidity, 
% temperature, particulate matter (PM10), pressure and wind speed, in residential sites of
% Bangkok, the capital of Thailand, during February 15 to May 15, 2007-2009. The hourly data were 
% hourly measured from eight weather stations during 10.00 A.M. to 3.00 P.M.
% and they are standardized to have zero mean and unit variance. We split data into training 
% and validation sets with the ratio of 2 : 1. 
% 
% In order to obtain a reasonable constraint, we perform a
% partial correlation analysis of 11 observed variables using partialcorr()
% in MATLAB.
% 
% The preprocessed data is located in "climate_raw_data_input" folder.
% The data. It contains 3 variables as
%   y               : observed matrix [N x n]
%   n               : number of observed variable 
%   indz_assump_A   : a list of zero index of A obtatined from partial
%                     correlation analysis
% 
% Developer: Anupon Pruttiakaravanich, Dr. Jitkomut Songsiri
% 
%% experiment of air pollution and weather data
clear; clc; close all;

addpath ..\codes\climate_raw_data_input

load nonsi_data.mat

num_candidate = 25;             %number of candidate model (#point of gamma along x-axis)
N = length(y_nonsi(:,1));       %number of observation    
N = ceil((2/3)*N);          
Y = y_nonsi(1:N,:);             %using 2/3 for training
                                %leaving 1/3 for validation
                                                    
%standardize Y                        
sd = std(Y,0,1);        %sd
md = mean(Y);           %mean
for i=1:n
Y(:,i) = (Y(:,i)-md(i))/sd(i);
end

indA = indz_assump_A          	%assumption for path matrix from patrial correlation analysis
                                %indz_assump_A is obtained from partial anylysis during the cleaning data process

[output] = modelsel_sem(num_candidate,Y',indA);

