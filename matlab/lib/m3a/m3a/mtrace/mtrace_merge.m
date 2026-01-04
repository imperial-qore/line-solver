function [T L] = mtrace_merge(t1, t2)
% Merges two traces in a single marked (multiclass) trace
% Input
%   t1: inter-arrival times of the first trace
%   t2: inter-arrival times of the second trace
% Output
%   T: inter-arrival times of the marked process
%   L: labels of the marked process

[TCUM,IDX] = sort([0; cumsum(t1); cumsum(t2)]);
T = diff(TCUM);
IDX = IDX(2:end);
L = zeros(length(T),1);
L(IDX >= 2 & IDX <= 1+length(t1)) = 1;
L(IDX >= 2+length(t1) & IDX <= 1+length(t1)+length(t2)) = 2;