function [RD, out] = getCdfFirstPassT(self, A, B)
% [RD, OUT] = GETCDFFIRSTPASST(A, B)
%
% Distribution of the FIRST PASSAGE TIME from the state set A into the state
% set B, on the CTMC underlying this model. RD is an [n x 2] matrix whose
% first column is F(t) and whose second is t, the column order every other
% CDF getter in LINE uses.
%
% A and B name states either as 1-based ROW INDICES into the state space
% returned by getStateSpace, or as matrices of state rows, which are resolved
% against that space. An empty A starts from the conditional stationary law on
% the complement of B.
%
% THIS IS NOT getCdfRespT. That getter times a tagged job between an arrival
% at a station and its departure, through the event filtration; this one times
% the chain between two sets of states the caller names, and answers questions
% the filtration cannot express -- the writer cycle time of a readers-writers
% model, the time to fill a buffer, the time to leave a degraded region.
%
% Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time Distributions
% in Large Markov Chains", 2002.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% lang='cpp' serves this getter from line-cli's -a firstpasst arm. The state
% sets cross the boundary as STATE ROWS, not row indices: the C++ engine does
% not purge vanishing states, so the two enumerations need not agree row for
% row, while a state row resolves by content in either space.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    T0 = tic;
    space = self.getStateSpace();
    n = size(space,1);
    [Aidx, Bidx] = local_resolve(self, A, B, n);
    [RD, out] = CPPLINE.cdfFirstPassT(self.name, self.model, self.options, ...
        space(Aidx,:), space(Bidx,:));
    out.source = Aidx;
    out.target = Bidx;
    out.runtime = toc(T0);
    return
end

T0 = tic;
sn = self.getStruct;
Q = self.getGenerator();
n = size(Q,1);

[Aidx, Bidx] = local_resolve(self, A, B, n);

config = self.getOptions.config;
if isfield(config,'passage_method') && ~isempty(config.passage_method)
    method = config.passage_method;
else
    method = 'expm';
end

if isempty(Aidx)
    pi0 = [];
else
    pi0 = zeros(1,n);
    pi0(Aidx) = 1/numel(Aidx);
end

% The horizon is chosen the way the response-time getter chooses it: 100
% events at the slowest rate in the chain.
nonZeroRates = abs(Q(Q~=0));
nonZeroRates = nonZeroRates(nonZeroRates > GlobalConstants.FineTol);
Thor = abs(100/min(nonZeroRates));
tset = linspace(0, Thor, 1000);

opt = struct('method', method);
[F, f, out] = ctmc_passage_time(Q, pi0, Bidx, tset, opt);
RD = [F(:), tset(:)];
out.tset = tset;
out.density = f;
out.source = Aidx;
out.target = Bidx;
out.runtime = toc(T0);
end

function [Aidx, Bidx] = local_resolve(self, A, B, n)
% A state set may be given as row indices or as state rows; both resolve to
% row indices here, and an unrecognised row is an error rather than a silent
% drop, since a passage into a state that is not in the space is not a slow
% passage but an undefined one.
Aidx = local_one(self, A, n, 'A');
Bidx = local_one(self, B, n, 'B');
if isempty(Bidx)
    line_error(mfilename, 'The target state set B is empty: a first passage time into no state is undefined.');
end
end

function idx = local_one(self, S, n, name)
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
        line_error(mfilename, sprintf('A state given in set %s is not in the state space.', name));
    end
    idx(i) = r;
end
idx = unique(idx);
end
