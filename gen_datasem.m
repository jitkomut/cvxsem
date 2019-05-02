% This function is to generate data. It requires in put as follows.
%
% n : number of variables
% N : number of observations
% density : density of nonzero entries
% 
% And it returns out put as follows.
% 
% S : sample covariance matrix
% Y : observed matrix variable
% Atrue : the true path matrix
% 
% example
% 
% n = 10;
% N= 100;
% density = 0.3;
% [output] = gen_datasem(n,N,density)

function [output] = gen_datasem(n,N,density)
mean = 0;                               %mean of data;
var = 0.1;                              %variance of data
e = mean + sqrt(var) * randn(n,N );     %noise

Atrue = sprandn(n,n,density);
Atrue = Atrue - diag(diag(Atrue));
Y = (eye(n)-Atrue)\e;
S = cov(Y');

field1 = 'S';               value1 = S;
field2 = 'Y';               value2 = Y;
field3 = 'Atrue';         	value3 = Atrue;
output = struct(field1,value1,field2,value2,field3,value3);