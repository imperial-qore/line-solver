function [S,Q11,Q12,Q21,Q22,T]=ctmc_stochcomp(Q,I)
% [S,Q11,Q12,Q21,Q22,T] = CTMC_STOCHCOMP(Q,I)
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if nargin==1
    I=1:ceil(length(Q)/2);
end
Ic = setdiff(1:length(Q),I);
Q11 = Q(I,I);
Q12 = Q(I,Ic);
Q21 = Q(Ic,I);
Q22 = Q(Ic,Ic);
%I = eye(size(Q22));
T = (-Q22) \ Q21;
T = Q12*T;
S = Q11+T;
end
