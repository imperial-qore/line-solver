function [m, mall] = getFirstPassTMoments(self, A, B, nmax)
% [M, MALL] = GETFIRSTPASSTMOMENTS(A, B, NMAX)
%
% Moments of order 1..NMAX of the first passage time from the state set A into
% the state set B. M is the (1 x NMAX) moment vector for a passage started
% uniformly in A; MALL is (nstates x NMAX), one row per starting state, zero on
% B and Inf where B cannot be reached.
%
% NO TRANSFORM INVERSION AND NO TIME GRID ARE INVOLVED. The moments come from
% Eq. 3 of P. G. Harrison and W. J. Knottenbelt, "Passage Time Distributions in
% Large Markov Chains", 2002 -- one linear solve per order -- so they are exact
% and are not limited by the horizon a CDF would have to be truncated at. This
% is the cheapest way to get the variance or the skewness of a passage time in
% LINE.
%
% A and B name states as in GETCDFFIRSTPASST.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4 || isempty(nmax)
    nmax = 3;
end
% lang='cpp' serves this getter from line-cli's -a firstpasstmom arm, added
% alongside the getter itself. The state sets cross the boundary as STATE ROWS
% rather than as row indices, because the two enumerations need not order (or
% even purge) the space identically -- the same convention getCdfFirstPassT
% uses, and the reason local_idx below resolves them here first.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    n = size(self.getStateSpace(), 1);
    space = self.getStateSpace();
    Bidx = local_idx(self, B, n);
    Aidx = local_idx(self, A, n);
    if isempty(Bidx)
        line_error(mfilename, 'The target state set B is empty: a first passage time into no state is undefined.');
    end
    Arows = space(Aidx, :);
    Brows = space(Bidx, :);
    [m, mall] = CPPLINE.firstPassTMoments(self.name, self.model, self.options, Arows, Brows, nmax);
    return
end

Q = self.getGenerator();
n = size(Q,1);
if isempty(B)
    line_error(mfilename, 'The target state set B is empty: a first passage time into no state is undefined.');
end

Bidx = local_idx(self, B, n);
Aidx = local_idx(self, A, n);
if isempty(Aidx)
    pi0 = [];
else
    pi0 = zeros(1,n);
    pi0(Aidx) = 1/numel(Aidx);
end

[mall, m] = ctmc_passage_moments(Q, pi0, Bidx, nmax);
end

function idx = local_idx(self, S, n)
if isempty(S)
    idx = [];
    return
end
if isvector(S) && all(S == round(S)) && all(S >= 1) && all(S <= n)
    idx = unique(reshape(S,1,[]));
    return
end
space = self.getStateSpace();
idx = zeros(1, size(S,1));
for i = 1:size(S,1)
    r = matchrow(space, reshape(S(i,:),1,[]));
    if r <= 0
        line_error(mfilename, 'A state given to getFirstPassTMoments is not in the state space.');
    end
    idx(i) = r;
end
idx = unique(idx);
end
