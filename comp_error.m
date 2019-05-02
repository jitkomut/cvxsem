% Compute error of true sparsity pattern
% ind : the linear index of the first pattern
% indx : the linear index of the second pattern
% Both can be generated from sub2ind([n n],i,j) if (i,j) are indices of zero
% n : dimension

function [true_pos,true_neg,false_pos,false_neg] = comp_error(Atrue,Ahat,n)

ind1 = find(Atrue ~=0);                     %index of nonzero entries in Atrue
ind2 = find(Atrue == 0);                    %index of zero entries in Atrue
indx = find(Ahat ~=0);
indz = find(Ahat == 0);

corr_ind1 = intersect(ind1,indx,'rows');    %intersect the true nonzero entries index
corr_ind2 = intersect(ind2,indz,'rows');    %intersect the true zero intries index
wrong_ind = setdiff(indx,ind1,'rows');
miss_ind = setdiff(ind1,indx,'rows');

% correct nonzero entries
A = zeros(n);
A(corr_ind1) = 1;
% spy(A,'s');hold on;
true_pos = sum(sum(A));

% correct zero entries
A = zeros(n);
A(corr_ind2) = 1;
% spy(A,'s');hold on;
true_neg = sum(sum(A)) - n;          % minus with diagonal entries since it is always zero

% plot wrongly-identified sparsity
A = zeros(n);
A(wrong_ind) = 1;
% spy(A,'or');
false_pos = sum(sum(A));

% plot missed sparsity
A = zeros(n);
A(miss_ind) = 1;
% spy(A,'+k');
false_neg = sum(sum(A));

% error = (false_pos + false_neg)/((n^2)-n);