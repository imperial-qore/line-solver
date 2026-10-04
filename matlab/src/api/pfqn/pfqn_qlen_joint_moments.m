function out = pfqn_qlen_joint_moments(L, N, Z, pairs, route, lGsrc, options)
% out = pfqn_qlen_joint_moments(L, N, Z, pairs, route, lGsrc, options)
%
% Joint moments of the queue-length vector of a closed product-form network,
% obtained from normalizing constants.
%
% The coordinates are (station,class) pairs. Two pairs sharing a class give the
% cross-station covariance of that class; two pairs sharing a station give the
% cross-class covariance at that station, which is what a class-oriented method
% of moments (pfqn_comomrm and its relatives) is positioned to deliver. Two
% exact routes reach the joint survival array, and both end in the same
% conversion, the tail edge of the house of moments (api/moment) followed by
% the joint central-moment and cumulant conversions:
%
%   - SINGLE CLASS (R = 1), route 'tail'. The survival probabilities are ratios
%     of normalizing constants of the network itself,
%
%       P(n_i >= k_i for all i) = (prod_i L_i^k_i) * G(N - sum_i k_i) / G(N)
%
%     which holds because a load-independent single-class station has the
%     geometric occupancy L_i^n. Only N+1 constants of the ORIGINAL model are
%     needed, which is why any normalizing-constant algorithm serves it.
%
%   - MULTICLASS, route 'pmf'. The geometric factorization fails, since a
%     multiclass load-independent station carries the multinomial occupancy
%     f_i(n_i) = (|n_i|)! prod_r L_ir^n_ir / n_ir!. What holds instead is the
%     joint law of the selected stations in terms of the COMPLEMENTARY network,
%     the model with those stations deleted and the think times kept,
%
%       P(n_i = m_i, i in S) = prod_i f_i(m_i) * G_(S^c)(N - sum_i m_i) / G(N)
%
%     The survival array is the reverse cumulative sum of that array, exactly,
%     since the box covers the support.
%
% Neither the factorial nor the raw moments have a one-constant closed form;
% the survival array is the queue-length functional that does. The
% normalizing-constant algorithm is INJECTED rather than called at a fixed
% site: the whole set of populations is known before any evaluation, so it is
% emitted in one batch and an algorithm that produces several constants in one
% pass serves it without recomputation.
%
% Input:
%   L: service demand matrix of the QUEUEING stations (MxR). Delay stations
%      belong in Z, their marginals following a different law. Load-dependent
%      and multiserver stations are out of scope for both routes
%   N: population vector (1xR)
%   Z: think time vector (1xR), zeros if empty
%   pairs: Px2 matrix of 1-based (station,class) pairs, one per dimension of
%          the returned arrays. Defaults to every class of every station
%   route: 'auto' (default), 'tail' or 'pmf'
%   lGsrc: where log G comes from. Empty calls pfqn_nc. A function handle is
%          invoked ONCE per network as lGsrc(Lsub, pops), pops being a PxR
%          matrix of populations, and must return P values of log G with NaN
%          where it cannot serve; those are filled in by pfqn_nc. A numeric
%          array is read as a table indexed by population, which is what a
%          convolution sweep produces for free. The 'pmf' route queries the
%          COMPLEMENTARY network, so a table must be its table
%   options: options struct passed to pfqn_nc. Defaults to
%            SolverNC.defaultOptions with method 'exact', since an approximate
%            normalizing constant would silently make the moments approximate
%
% Output:
%   out: struct with the joint arrays over the selected coordinates (tail,
%        binomial, factorial, raw, central, cumulant), the mean vector, the
%        covariance matrix cov, and info holding route, points, served, evals
%        and the pairs used
%
% Example:
%   out = pfqn_qlen_joint_moments([2 1], [4 3], [0.5 0.8], [1 1; 1 2]);
%   cov12 = out.cov(1,2);
%
% Reference:
% M. Reiser and S. S. Lavenberg. Mean-value analysis of closed multichain
% queuing networks. Journal of the ACM, 27(2):313-322, 1980.

[M,R] = size(L);
N = N(:).';
if numel(N) ~= R
    line_error(mfilename,'The population vector N must have one entry per class.');
end
if nargin < 3 || isempty(Z)
    Z = zeros(1,R);
end
Z = Z(:).';
if numel(Z) ~= R
    line_error(mfilename,'The think time vector Z must have one entry per class.');
end
if nargin < 4 || isempty(pairs)
    pairs = zeros(M*R,2);
    t = 0;
    for i = 1:M
        for r = 1:R
            t = t + 1;
            pairs(t,:) = [i, r];
        end
    end
end
if size(pairs,2) ~= 2 || isempty(pairs)
    line_error(mfilename,'The pairs must be given as a Px2 matrix of (station,class) indices.');
end
if any(pairs(:,1) < 1) || any(pairs(:,1) > M) || any(pairs(:,2) < 1) || any(pairs(:,2) > R)
    line_error(mfilename,'A (station,class) pair is out of range.');
end
if size(unique(pairs,'rows'),1) ~= size(pairs,1)
    line_error(mfilename,'The (station,class) pairs must be distinct.');
end
if nargin < 5 || isempty(route)
    route = 'auto';
end
if nargin < 6
    lGsrc = [];
end
if nargin < 7 || isempty(options)
    % the point of this routine is an EXACT moment array, so the default is the
    % exact normalizing constant rather than the adaptive dispatch of pfqn_nc.
    % parseOptions replaces the defaults wholesale, so the struct is built from
    % SolverNC.defaultOptions and only the method is overridden
    options = SolverNC.defaultOptions;
    options.method = 'exact';
end
if strcmp(route,'auto')
    if R == 1
        route = 'tail';
    else
        route = 'pmf';
    end
end
if strcmp(route,'tail') && R > 1
    line_error(mfilename,'The tail route needs the geometric occupancy of a single-class load-independent station; with several classes the multinomial factor breaks the survival identity, so use the pmf route.');
end
if ~any(strcmp(route,{'tail','pmf'}))
    line_error(mfilename,'The route must be auto, tail or pmf.');
end

d = size(pairs,1);
dims = zeros(1,d);
for j = 1:d
    dims(j) = N(pairs(j,2)) + 1;
end

if strcmp(route,'tail')
    % the whole population set is known up front: N minus the total order
    nel = prod(dims);
    need = [];
    a = ones(1,d);
    for ia = 1:nel
        s = sum(a-1);
        if N - s >= 0
            need(end+1,:) = N - s; %#ok<AGROW>
        end
        a = local_odometer(a, dims);
    end
    need = unique([need; N], 'rows');
    [lg, served, evals] = local_batch_lg(L, need, Z, lGsrc, options);
    lgN = lg(local_findrow(need, N));
    tail = zeros([dims 1]);
    a = ones(1,d);
    for ia = 1:nel
        s = sum(a-1);
        if N - s >= 0
            acc = 0;
            ok = true;
            for j = 1:d
                if a(j) > 1
                    if L(pairs(j,1),pairs(j,2)) <= 0
                        ok = false;
                        break
                    end
                    acc = acc + (a(j)-1)*log(L(pairs(j,1),pairs(j,2)));
                end
            end
            if ok
                tail(ia) = exp(acc + lg(local_findrow(need, N - s)) - lgN);
            end
        end
        a = local_odometer(a, dims);
    end
else
    % the joint law of the selected stations needs every class of those
    % stations, so the internal box runs over (station,class) and the requested
    % pairs are marginalized out of it afterwards
    stations = unique(pairs(:,1)).';
    ns = numel(stations);
    coords = zeros(ns*R,2);
    t = 0;
    for i = stations
        for r = 1:R
            t = t + 1;
            coords(t,:) = [i, r];
        end
    end
    dc = size(coords,1);
    cdims = zeros(1,dc);
    for j = 1:dc
        cdims(j) = N(coords(j,2)) + 1;
    end
    Lsub = L;
    Lsub(stations,:) = [];
    ncel = prod(cdims);
    need = [];
    a = ones(1,dc);
    for ia = 1:ncel
        n = N;
        for j = 1:dc
            n(coords(j,2)) = n(coords(j,2)) - (a(j)-1);
        end
        if all(n >= 0)
            need(end+1,:) = n; %#ok<AGROW>
        end
        a = local_odometer(a, cdims);
    end
    need = unique(need, 'rows');
    if isempty(Lsub)
        lgc = zeros(size(need,1),1);
        for p = 1:size(need,1)
            lgc(p) = local_delay_lg(Z, need(p,:));
        end
        served = size(need,1);
        evals = 0;
    else
        [lgc, served, evals] = local_batch_lg(Lsub, need, Z, lGsrc, options);
    end
    [lgNv, ~, evals0] = local_batch_lg(L, N, Z, [], options);
    lgN = lgNv(1);
    evals = evals + evals0;

    marg = zeros([dims 1]);
    a = ones(1,dc);
    for ia = 1:ncel
        n = N;
        for j = 1:dc
            n(coords(j,2)) = n(coords(j,2)) - (a(j)-1);
        end
        if all(n >= 0)
            gc = lgc(local_findrow(need, n));
            if isfinite(gc)
                acc = gc - lgN;
                ok = true;
                for i = stations
                    tot = 0;
                    for j = 1:dc
                        if coords(j,1) == i
                            tot = tot + (a(j)-1);
                        end
                    end
                    acc = acc + gammaln(tot+1);
                    for j = 1:dc
                        if coords(j,1) == i && a(j) > 1
                            if L(i,coords(j,2)) <= 0
                                ok = false;
                                break
                            end
                            acc = acc + (a(j)-1)*log(L(i,coords(j,2))) - gammaln(a(j));
                        end
                    end
                    if ~ok
                        break
                    end
                end
                if ok
                    sub = ones(1,d);
                    for j = 1:d
                        for jc = 1:dc
                            if coords(jc,1) == pairs(j,1) && coords(jc,2) == pairs(j,2)
                                sub(j) = a(jc);
                            end
                        end
                    end
                    subc = num2cell(sub);
                    marg(subc{:}) = marg(subc{:}) + exp(acc);
                end
            end
        end
        a = local_odometer(a, cdims);
    end
    tail = marg;
    for mode = 1:d
        tail = flip(cumsum(flip(tail,mode),mode),mode);
    end
end

b = moment_joint_binomial_from_tail(tail);
f = moment_joint_factorial_from_binomial(b);
m = moment_joint_raw_from_factorial(f);
mc = moment_joint_central_from_raw(m);
kap = moment_joint_cumulant_from_raw(m);

meanv = zeros(1,d);
covm = zeros(d,d);
for j = 1:d
    e = ones(1,d);
    e(j) = 2;
    ec = num2cell(e);
    meanv(j) = m(ec{:});
    for l = 1:d
        aa = ones(1,d);
        aa(j) = aa(j) + 1;
        aa(l) = aa(l) + 1;
        ac = num2cell(aa);
        covm(j,l) = kap(ac{:});
    end
end

out = struct();
out.tail = tail;
out.binomial = b;
out.factorial = f;
out.raw = m;
out.central = mc;
out.cumulant = kap;
out.mean = meanv;
out.cov = covm;
out.info = struct('route', route, 'points', size(need,1), 'served', served, ...
    'evals', evals, 'exact', true, 'pairs', pairs, 'dims', dims);
end

function a = local_odometer(a, dims)
% Advance a 1-based multi-index, first dimension fastest.
for l = 1:numel(dims)
    a(l) = a(l) + 1;
    if a(l) <= dims(l)
        return
    end
    a(l) = 1;
end
end

function p = local_findrow(rows, key)
% Index of a population vector in the deduplicated request matrix.
p = find(all(rows == repmat(key, size(rows,1), 1), 2), 1);
end

function lg = local_delay_lg(Z, n)
% Log normalizing constant of a pure-delay network, prod_r Z_r^n_r / n_r!.
lg = 0;
for r = 1:numel(n)
    if n(r) == 0
        continue
    end
    if Z(r) <= 0
        lg = -Inf;
        return
    end
    lg = lg + n(r)*log(Z(r)) - gammaln(n(r)+1);
end
end

function [lg, served, evals] = local_batch_lg(Lsub, pops, Z, lGsrc, options)
% Evaluate log G at a batch of populations, honouring the injected source.
P = size(pops,1);
lg = nan(P,1);
served = 0;
if ~isempty(lGsrc)
    if isa(lGsrc,'function_handle')
        got = lGsrc(Lsub, pops);
        got = got(:);
        if numel(got) ~= P
            line_error(mfilename,'The lGsrc handle must return one value per requested population.');
        end
        lg = got;
        served = sum(isfinite(lg));
    elseif isnumeric(lGsrc)
        for p = 1:P
            sub = num2cell(pops(p,:) + 1);
            lg(p) = lGsrc(sub{:});
        end
        served = sum(isfinite(lg));
    else
        line_error(mfilename,'lGsrc must be empty, a function handle or a numeric table.');
    end
end
evals = 0;
R = size(pops,2);
for p = 1:P
    if isfinite(lg(p))
        continue
    end
    lg(p) = pfqn_nc(zeros(1,R), Lsub, pops(p,:), Z, options);
    evals = evals + 1;
end
end
