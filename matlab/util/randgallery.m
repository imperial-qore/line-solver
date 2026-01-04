function [L,it]=randgallery(M,R,it,maxint)
% L=RANDGALLERY(M,R,s,maxint)
% Return integer matrix of size MxR with seed s and entries in [0,maxint]
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if ~exist('it','var')
    it=randi(10000);
end
if ~exist('maxint','var')
    maxint=100;
end
L = gallery('integerdata', maxint, [M, R], it);
end
