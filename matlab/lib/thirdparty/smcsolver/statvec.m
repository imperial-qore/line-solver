function K = statvec(A)
% STAT(A) calculates a stationary distribution of a stochatic matrix A, i.e. a
% left eigenvector corresponding to eigenvalue 1 with nonnegative elements
% summing to one.
 
% Last modification: 7th March 1997


S = size(A,1);
e = ones(S,1);
B = [A - eye(S),e];
y = [zeros(1,S),1];
K = y / B;
