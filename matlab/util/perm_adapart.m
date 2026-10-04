function [val, info] = perm_adapart(A, options)
% [VAL, INFO] = PERM_ADAPART(A)
% [VAL, INFO] = PERM_ADAPART(A, OPTIONS)
%
% Adaptive partitioning (AdaPart) sampler for the permanent of a nonnegative
% matrix. The space of the permutations is recursively partitioned, each part
% is bounded by the Soules column bound, and a part is drawn with probability
% proportional to its bound; the acceptance ratio of the resulting rejection
% sampler scales the root bound into an unbiased estimate of the permanent.
%
% Twin of jline.lib.perm.AdaPartSampler and of the python
% line_solver.api.perm.AdaPartSampler. The three codebases agree in
% distribution but not sample by sample, since each draws from its own
% generator.
%
% Input:
%   A       - nonnegative square matrix
%   OPTIONS - optional struct with fields
%             mode        'classic' (default), 'time' or 'sample'
%             maxAccepted acceptance budget of 'classic' mode (default 100)
%             maxTime     time budget in ms of 'time' mode (default 30000)
%             maxSamples  draw budget of 'sample' mode (default 450)
%             seed        seed of a private mt19937ar stream (default: global)
%
% Output:
%   VAL  - estimate of the permanent
%   INFO - struct with fields accepted (0/1 per draw), time (ms per draw),
%          permStep (running estimate) and zub (root Soules bound)
%
% This is a Monte Carlo estimator: at the default budget expect a relative
% error near 1e-2, falling as 1/sqrt(budget). It is exact when the Soules
% bound is tight, which includes constant and column-constant matrices.
%
% Reference:
%   J. Kuck, T. Dao, H. Rezatofighi, A. Sabharwal, S. Ermon, "Approximating
%   the Permanent by Sampling from Adaptive Partitions", NeurIPS, 2019.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    options = struct();
end
if ~isfield(options, 'mode'), options.mode = 'classic'; end
if ~isfield(options, 'maxAccepted'), options.maxAccepted = 100; end
if ~isfield(options, 'maxTime'), options.maxTime = 30000; end
if ~isfield(options, 'maxSamples'), options.maxSamples = 450; end
if ~isfield(options, 'maxDraws'), options.maxDraws = 1e6; end
if ~isfield(options, 'seed'), options.seed = []; end

if any(A(:) < 0)
    line_error(mfilename, 'Matrix must be non-negative.');
end
if size(A,1) ~= size(A,2)
    line_error(mfilename, 'Matrix must be square.');
end

n = size(A, 1);
if n == 0
    val = 1; info = struct('accepted', [], 'time', [], 'permStep', [], 'zub', 1);
    return
end

perm_require_support(A, mfilename);

if isempty(options.seed)
    stream = RandStream.getGlobalStream();
else
    stream = RandStream('mt19937ar', 'Seed', options.seed);
end

zub = soules_bound(A, n);
accepted = 0;
total = 0;
acceptedList = [];
timeList = [];
t0 = tic;

switch options.mode
    case 'time'
        while toc(t0)*1000 < options.maxTime
            s = draw_sample(A, n, stream);
            accepted = accepted + s; total = total + 1;
            acceptedList(end+1) = s; %#ok<AGROW>
            timeList(end+1) = toc(t0)*1000; %#ok<AGROW>
        end
    case 'sample'
        while numel(acceptedList) < options.maxSamples
            s = draw_sample(A, n, stream);
            accepted = accepted + s; total = total + 1;
            acceptedList(end+1) = s; %#ok<AGROW>
            timeList(end+1) = toc(t0)*1000; %#ok<AGROW>
        end
    otherwise
        % Bounded independently of the scaling: with perm(A)=0 the acceptance
        % probability is 0 and this loop would never terminate.
        while accepted < options.maxAccepted
            if total >= options.maxDraws
                line_error(mfilename, ['Only %d of the %d required acceptances ' ...
                    'were obtained in %d draws. Raise options.maxDraws or use ' ...
                    'the exact engine.'], accepted, options.maxAccepted, total);
            end
            s = draw_sample(A, n, stream);
            accepted = accepted + s; total = total + 1;
            acceptedList(end+1) = s; %#ok<AGROW>
            timeList(end+1) = toc(t0)*1000; %#ok<AGROW>
        end
end

if total > 0
    val = zub * accepted / total;
else
    val = 0;
end

info.accepted = acceptedList;
info.time = timeList;
info.permStep = zub * cumsum(acceptedList) ./ (1:numel(acceptedList));
info.zub = zub;
end


function accept = draw_sample(A, n, stream)
% Draw one partition path, returning 1 if accepted and 0 if rejected.
% Column j of the assignment holds the row bound to it, 0 marking a free column.
S = zeros(1, n);

while any(S(:) == 0)
    sInit = S(1, :);
    ub = soules_bound(modify_matrix(A, sInit, n), n);
    zubS = ub;
    init = true;

    while (ub >= zubS) || init
        init = false;
        % Only elements with a free column can be refined; expanding a complete
        % assignment yields no children and stalls the sampler.
        refinable = find(any(S == 0, 2));
        if isempty(refinable)
            break
        end
        pick = refinable(randi(stream, numel(refinable)));
        sSub = S(pick, :);
        S(pick, :) = [];

        subMatrix = modify_matrix(A, sSub, n);
        % Discount the bound of the element actually removed, not the root
        % bound; see _kb/03-api-layer.md.
        subUb = soules_bound(subMatrix, n);
        [newUb, j] = select_column(subMatrix, subUb, ub, sSub, n);

        for i = 1:n
            if ~any(sSub == i)
                sAdd = sSub;
                sAdd(j) = i;
                S(end+1, :) = sAdd; %#ok<AGROW>
            end
        end
        S = unique(S, 'rows', 'stable');

        noProgress = newUb >= ub;
        ub = newUb;
        % The Soules bound is tight on matrices with equal entries, so
        % refinement cannot improve it and "refine until improved" would never
        % exit. Stop on the first non-improving expansion instead; a tight
        % bound means the draw is accepted with probability 1.
        if noProgress
            break
        end
    end

    c = compute_probabilities(S, zubS, A, n, stream);
    if c == size(S, 1) + 1
        accept = 0;
        return
    end
    S = pick_subset(S, c, n);
end
accept = 1;
end


function [newUb, j] = select_column(sMatrix, removedUb, ub, sSub, n)
% Pick the column whose expansion minimizes the summed Soules bound.
% Only columns still free in sSub are candidates: scoring an already-assigned
% column just re-derives its own constraint, which always looks cheapest, so
% the sampler would re-split the same column forever.
ubi = inf(1, n);
for i = 1:n
    if sSub(i) ~= 0
        continue
    end
    tot = 0;
    for jj = 1:n
        assignment = zeros(1, n);
        assignment(i) = jj;
        tot = tot + soules_bound(modify_matrix(sMatrix, assignment, n), n);
    end
    ubi(i) = tot;
end
[~, j] = min(ubi);
newUb = ub - removedUb + ubi(j);
end


function c = compute_probabilities(S, zubS, A, n, stream)
% Draw a partition element, or the slack index size(S,1)+1 meaning rejection.
nS = size(S, 1);
p = zeros(1, nS + 1);
for k = 1:nS
    p(k) = soules_bound(modify_matrix(A, S(k, :), n), n);
end
if sum(p(1:nS)) > 0
    p(1:nS) = p(1:nS) / zubS;
end
p(nS+1) = 1 - sum(p(1:nS));
tot = sum(p);
if tot > 0
    p = abs(p) / tot;
end
c = find(cumsum(p) >= rand(stream), 1, 'first');
if isempty(c)
    c = nS + 1;
end
end


function S = pick_subset(S, c, n)
% Keep the drawn element, completing it when a single free column is left.
sInter = S(c, :);
if sum(sInter == 0) == 1
    missing = setdiff(1:n, sInter(sInter ~= 0));
    sInter(find(sInter == 0, 1)) = missing(1);
end
S = sInter;
end


function out = modify_matrix(m, t, n)
% Zero out the entries excluded by the partial assignment t.
mask = false(n, n);
for j = 1:n
    if t(j) ~= 0
        mask(t(j), j) = true;
    end
end
for i = 1:n
    if ~any(t == i)
        for j = 1:n
            if t(j) == 0
                mask(i, j) = true;
            end
        end
    end
end
out = zeros(n, n);
out(mask) = m(mask);
end


function b = soules_bound(m, n)
% Soules upper bound of the permanent, a product of column bounds.
g = zeros(1, n+1);
f = 1;
for k = 1:n
    f = f * k;
    g(k+1) = f^(1/k);
end
delta = g(n+1:-1:2) - g(n:-1:1);
b = prod(delta * sort(m, 1));
end
