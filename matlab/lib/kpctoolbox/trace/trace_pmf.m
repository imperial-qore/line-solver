function [pmf,px] = trace_pmf(X)
% MAX(X) equally spaced bins spanning the data range, which is the layout
% HIST(X, MAX(X)) used. Naming the edges keeps that layout explicit and is not
% subject to the 65536-bin ceiling histcounts applies to a scalar bin count.
nbins = max(X);
[pmf] = histcounts(X, linspace(min(X), nbins, nbins+1))' ./ numel(X); 
px = unique(X);
end
