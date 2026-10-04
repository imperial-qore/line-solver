function [QLen, Loss, lG, niter] = lossn_rec(nu, A, C)
% [QLEN, LOSS, LG, NITER] = LOSSN_REC(NU, A, C)
% Exact analysis of a loss network by MDD-rec: the normalising constant is the
% sum of a product form over the admissible set {n >= 0 : A n <= C}, which is
% what a decision diagram holding that set computes in one memoised walk.
%
% A Kelly loss network carries offered load nu_r on route r and admits a call
% only while the resource constraint A n <= C still holds after it. The
% stationary law is the truncation of independent Poisson counts to that set,
%
%   P(n) = (1/G) prod_r nu_r^{n_r} / n_r!,   G = sum_{A n <= C} prod_r ...,
%
% so g_r(k) = nu_r^k/k! and MDD_REC returns G. By PASTA the acceptance
% probability of a class-r call is the ratio of two such constants,
%
%   1 - B_r = G(C - A e_r) / G(C),
%
% which is one further diagram per class.
%
% WHY THIS EXISTS ALONGSIDE LOSSN_MANJUNATH. The Manjunath-Sikdar transform
% evaluates G exactly as a multidimensional residue, and the residue argument
% counts WHOLE UNITS: it needs an integral A and C. On a region declaring a
% fractional class size or capacity the analyzer had no exact route at all and
% fell back to the Erlang fixed point, an approximation. MDD-rec needs only
% that the admissible set be finite and bounded coordinate by coordinate, which
% a fractional constraint still is, so it is exact there too. It is also an
% exact alternative to the Monte Carlo summation LOSSN_MCI estimates.
%
% -- Input
% NU   : 1 x K offered load per class
% A    : J x K non-negative resource requirement matrix
% C    : J x 1 capacity vector
% -- Output
% QLEN : 1 x K mean number of class-r calls in progress (the carried load)
% LOSS : 1 x K blocking probability
% LG   : log of the normalising constant G(C)
% NITER: number of diagram walks performed, K + 1
%
% -- Reference
% F. P. Kelly, "Loss networks", Annals of Applied Probability 1(3), 1991.
% S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant for
% product-form models of distributed systems with synchronisation", Future
% Generation Computer Systems 111 (2020) 475-490.
%
% See also LOSSN_MANJUNATH, LOSSN_ERLANGFP, LOSSN_MCI, MDD_REC.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

nu = nu(:)';
C = C(:)';
K = numel(nu);
if size(A, 2) ~= K
    line_error(mfilename, sprintf('A has %d columns but there are %d classes', size(A, 2), K));
end
if size(A, 1) ~= numel(C)
    line_error(mfilename, sprintf('A has %d rows but C has %d entries', size(A, 1), numel(C)));
end
if any(A(:) < 0)
    line_error(mfilename, 'the resource matrix A must be non-negative');
end

% ---- per-class bound: the most calls the tightest constraint alone admits
bound = zeros(1, K);
for r = 1:K
    j = find(A(:, r) > 0);
    if isempty(j)
        line_error(mfilename, sprintf(['class %d consumes no resource, so the admissible set is ' ...
            'unbounded in that coordinate and its normalising constant diverges'], r));
    end
    bound(r) = floor(min(C(j) ./ A(j, r)'));
    if bound(r) < 0, bound(r) = 0; end
end

g = cell(1, K);
for r = 1:K
    k = 0:bound(r);
    g{r} = (nu(r) .^ k) ./ gamma(k + 1);
end

lGfull = i_logG(A, C, bound, g);
if ~isfinite(lGfull)
    line_error(mfilename, 'the admissible set is empty: no call of any class fits within C');
end
lG = lGfull;

% ---- carried load per class, from the marginals of the same diagram
mdds = i_diagram(A, C, bound);
G = exp(lG);
QLen = zeros(1, K);
for r = 1:K
    pk = mdd_rec_marginal(mdds, g, r) / G;
    QLen(r) = (0:numel(pk) - 1) * pk(:);
end

% ---- blocking: 1 - B_r = G(C - A e_r)/G(C), Kelly's ratio, by PASTA
Loss = zeros(1, K);
for r = 1:K
    Cr = C - A(:, r)';
    if any(Cr < 0)
        Loss(r) = 1;                            % the call never fits
        continue
    end
    lGr = i_logG(A, Cr, bound, g);
    if ~isfinite(lGr)
        Loss(r) = 1;
    else
        Loss(r) = 1 - exp(lGr - lG);
    end
    Loss(r) = min(1, max(0, Loss(r)));
end

niter = K + 1;
end

% ------------------------------------------------------------------------
function mdds = i_diagram(A, C, bound)
% The admissible set {n >= 0 : A n <= C}, generated one call at a time from the
% empty network. Adding a call is the only move, so the breadth-first closure
% visits exactly the admissible vectors.
K = numel(bound);
domain = bound + 1;
nextfun = @(s) i_next(s, A, C, bound, K);
mdd = mdd_reachset(domain, zeros(1, K), nextfun);
mdds = mdd.toStruct();
end

% ------------------------------------------------------------------------
function T = i_next(s, A, C, bound, K)
T = zeros(0, K);
for r = 1:K
    if s(r) >= bound(r), continue; end
    t = s; t(r) = t(r) + 1;
    if all(A * t(:) <= C(:) + 1e-12)
        T(end + 1, :) = t; %#ok<AGROW>
    end
end
end

% ------------------------------------------------------------------------
function lG = i_logG(A, C, bound, g)
% log G over the admissible set at capacity C, keeping the per-class domains of
% the FULL problem so that one set of factors g serves every reduced capacity.
if any(C < 0), lG = -Inf; return; end
mdds = i_diagram(A, C, bound);
G = mdd_rec(mdds, g);
if ~(G > 0), lG = -Inf; else, lG = log(G); end
end
