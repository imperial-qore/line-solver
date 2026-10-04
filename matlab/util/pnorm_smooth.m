function ghat = pnorm_smooth(x, c, p)
% GHAT = PNORM_SMOOTH(X, C, P)
%
% Smoothed processor-share constraint approximation using p-norm.
% Based on Ruuskanen et al., PEVA 151 (2021).
%
% @param x Total queue length at station
% @param c Number of servers (capacity)
% @param p Smoothing parameter (larger = closer to min function)
%
% @return ghat Smoothed capacity constraint value in [0, 1]
%
% The p-norm smooth function approximates min(x,c)/x more smoothly than
% the softmin function, which helps ODE stability.
%
% Formula: ghat = 1 / (1 + (x/c)^p)^(1/p)
%
% As p -> infinity, ghat -> min(1, c/x) = min(x,c)/x
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if x <= 0 || c <= 0
    ghat = 0;
    return;
end

ratio = x / c;
if p <= 0
    ghat = min(1, c/x); % fallback to hard min
else
    ghat = 1 / (1 + ratio^p)^(1/p);
end

if isnan(ghat)
    ghat = 0;
end
end
