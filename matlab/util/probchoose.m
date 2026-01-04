function [pos,f] = probchoose(p)
% pos = PROBCHOOSE(P)
% Choose an element according to probability vector P
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

f = cumsum(p);
r = rand;
pos = maxpos(r<=f);
end