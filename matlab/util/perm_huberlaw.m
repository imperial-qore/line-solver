function [val, info] = perm_huberlaw(A, options)
% [VAL, INFO] = PERM_HUBERLAW(A)
% [VAL, INFO] = PERM_HUBERLAW(A, OPTIONS)
%
% Huber-Law acceptance-rejection sampler for the permanent of a nonnegative
% matrix. The matrix is rescaled to be doubly stochastic, a permutation is
% drawn column by column under the Huber-Law upper bound on the permanent of
% what is left, and the acceptance ratio times the rescaling constant is an
% unbiased estimate of the permanent.
%
% Twin of jline.lib.perm.HuberLawSampler and of the python
% line_solver.api.perm.HuberLawSampler. The three codebases agree in
% distribution but not sample by sample, since each draws from its own
% generator.
%
% Input:
%   A       - nonnegative square matrix
%   OPTIONS - optional struct with fields
%             mode       'classic' (default), 'time' or 'sample'
%             delta      relative accuracy target, sets the acceptance budget
%                        K = 14*delta^-2*log(2/epsilon) (default 0.1)
%             epsilon    failure probability target (default 0.1)
%             alpha2     convergence threshold of the doubly stochastic
%                        rescaling (default 1e-6)
%             maxSamples draw budget of 'sample' mode (default 1000)
%             maxTime    time budget in ms of 'time' mode (default 30000)
%             seed       seed of a private mt19937ar stream (default: global)
%
% Output:
%   VAL  - estimate of the permanent
%   INFO - struct with fields accepted (0/1 per draw), time (ms per draw),
%          permStep (running estimate) and cst (rescaling constant)
%
% ALPHA3, the greedy assignment lower bound that sets the flooring level, is
% taken on the SCALED matrix. Computing it on A while flooring A/max(A) floors
% nearly every entry and overestimates the permanent by orders of magnitude on
% any matrix whose entries are not already O(1); see _kb/03-api-layer.md.
%
% Reference:
%   M. Huber, J. Law, "Fast approximation of the permanent for very dense
%   problems", SODA, 2008.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    options = struct();
end
if ~isfield(options, 'mode'), options.mode = 'classic'; end
if ~isfield(options, 'delta'), options.delta = 0.1; end
if ~isfield(options, 'epsilon'), options.epsilon = 0.1; end
if ~isfield(options, 'alpha2'), options.alpha2 = 1e-6; end
if ~isfield(options, 'maxSamples'), options.maxSamples = 1000; end
if ~isfield(options, 'maxTime'), options.maxTime = 30000; end
if ~isfield(options, 'maxSinkhorn'), options.maxSinkhorn = 10000; end
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
    val = 1;
    info = struct('accepted', [], 'time', [], 'permStep', [], 'cst', 1);
    return
end

perm_require_support(A, mfilename);

% A private stream keeps the global one untouched, as in perm_adapart.
if isempty(options.seed)
    stream = [];
else
    stream = RandStream('mt19937ar', 'Seed', options.seed);
end

[C, cst] = sub_rescale(A, n, options);

t0 = tic;
accepted = [];
elapsed = [];
switch options.mode
    case 'time'
        while toc(t0)*1000 < options.maxTime
            [accepted, elapsed] = sub_draw(C, n, stream, accepted, elapsed, t0);
        end
    case 'sample'
        while numel(accepted) < options.maxSamples
            [accepted, elapsed] = sub_draw(C, n, stream, accepted, elapsed, t0);
        end
    otherwise
        K = floor(14 * options.delta^(-2) * log(2/options.epsilon));
        % Bounded independently of the scaling: with perm(A)=0 the acceptance
        % probability is 0 and this loop would never terminate. A cap that
        % RETURNS a number would be a workaround, so it raises.
        while sum(accepted) < K
            if numel(accepted) >= options.maxDraws
                line_error(mfilename, ['Only %d of the %d required acceptances were ' ...
                    'obtained in %d draws. The acceptance probability is too low ' ...
                    'for this budget; raise options.maxDraws, relax options.delta, ' ...
                    'or use the exact engine.'], sum(accepted), K, numel(accepted));
            end
            [accepted, elapsed] = sub_draw(C, n, stream, accepted, elapsed, t0);
        end
end

if isempty(accepted)
    val = 0;
else
    val = sum(accepted)/numel(accepted) * cst;
end
info = struct('accepted', accepted, 'time', elapsed, ...
    'permStep', cst * cumsum(accepted) ./ (1:numel(accepted)), 'cst', cst);
end

function [accepted, elapsed] = sub_draw(C, n, stream, accepted, elapsed, t0)
sigma = sub_sample(C, n, stream);
accepted(end+1) = double(~any(sigma == n+1)); %#ok<AGROW>
elapsed(end+1) = toc(t0)*1000; %#ok<AGROW>
end

function sigma = sub_sample(C, n, stream)
% Draw a permutation column by column. A rejected draw is flagged by returning
% the out-of-range row index n+1 in every position (the JAR and python use n,
% which is out of range there because they index from 0).
M = C;
sigma = zeros(1, n);
for j = 1:n
    ub = prod(sub_h(sum(M, 2))) / exp(n);
    p = sub_precomputing(M, j, n) / ub;
    prob = [p(:)', 1 - sum(p)];
    if prob(end) < 0
        prob(1:n) = prob(1:n) / sum(prob(1:n));
        prob(end) = 0;
    end
    if isempty(stream)
        u = rand();
    else
        u = rand(stream);
    end
    sel = n + 1;
    csum = 0;
    for i = 1:(n+1)
        csum = csum + prob(i);
        if u <= csum
            sel = i;
            break
        end
    end
    if sel == n + 1
        sigma = (n+1) * ones(1, n);
        return
    end
    sigma(j) = sel;
    keep = true(n, n);
    keep(sel, :) = false;
    keep(:, j) = false;
    Mnew = zeros(n, n);
    Mnew(keep) = M(keep);
    Mnew(sel, j) = M(sel, j);
    M = Mnew;
end
end

function p = sub_precomputing(M, j, n)
% Unnormalized selection weight of each row for column J.
c = M(:, j);
r = sum(M, 2) - c;
hr = sub_h(r);
p = prod(hr) ./ hr .* c / exp(n-1);
end

function hr = sub_h(r)
% Huber-Law bound factor of a row holding remaining mass R.
hr = zeros(size(r));
big = r >= 1;
hr(big) = r(big) + 0.5*log(max(r(big),1)) + exp(1) - 1;
hr(~big) = 1 + (exp(1)-1)*r(~big);
end

function [C, cst] = sub_rescale(A, n, options)
% Rescale to doubly stochastic and return the constant that undoes it.
% A is strictly positive here (perm_require_support above), so no log floor is
% needed, and the assignment is the real maximum-weight one rather than the
% row-by-row greedy that used to stand in for it.
logA = log(A);
assignment = perm_maxweight_assignment(logA);

maxel = max(A(:));
Ms = A / maxel;

% alpha3 is a permanent lower bound OF THE SCALED MATRIX, which is the one
% floored on the next line.
alpha3 = 1;
for i = 1:n
    alpha3 = alpha3 * Ms(i, assignment(i));
end
alpha1 = alpha3 * options.delta / 3 / factorial(n);
Ms = max(Ms, alpha1);

[B, X, Y] = sub_doubly_stochastic(Ms, n, options.alpha2, options.maxSinkhorn);

Z = eye(n);
for i = 1:n
    Z(i,i) = 1 / max(B(i,:));
end
C = Z * B;

hprod = prod(sub_h(sum(C, 2)) / exp(1));
dprod = prod(diag(X) .* diag(Y) .* diag(Z));
cst = hprod / dprod * maxel^n;
end

function [B, X, Y] = sub_doubly_stochastic(M, n, alpha2, maxSweeps)
% Alternate column and row normalization until both margins are within ALPHA2.
% Capped: a row that sums to zero leaves ROWERR at 1 forever, and the guarded
% normalization below skips it, so the loop used to spin without terminating.
B = M;
X = eye(n);
Y = eye(n);
rowErr = Inf;
colErr = Inf;
sweeps = 0;
while rowErr > alpha2 || colErr > alpha2
    sweeps = sweeps + 1;
    if sweeps > maxSweeps
        line_error('perm_huberlaw', ['The doubly stochastic rescaling did not ' ...
            'converge in %d sweeps (row error %g, column error %g against a ' ...
            'tolerance of %g). The usual cause is a matrix without total ' ...
            'support.'], maxSweeps, rowErr, colErr, alpha2);
    end
    csum = sum(B, 1);
    for j = 1:n
        if csum(j) > 0
            B(:,j) = B(:,j) / csum(j);
            Y(:,j) = Y(:,j) / csum(j);
        end
    end
    rsum = sum(B, 2);
    for i = 1:n
        if rsum(i) > 0
            B(i,:) = B(i,:) / rsum(i);
            X(i,:) = X(i,:) / rsum(i);
        end
    end
    colErr = max(abs(sum(B, 1) - 1));
    rowErr = max(abs(sum(B, 2) - 1));
end
end
