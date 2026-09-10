function result = qsys_mgisrgi_whitt(lambda, mu, s, r, patience, varargin)
% QSYS_MGISRGI_WHITT Engineering solution of the M/GI/s/r+GI queue.
%
% RESULT = QSYS_MGISRGI_WHITT(LAMBDA, MU, S, R, PATIENCE) computes every
% standard steady-state measure of the call-center model M/GI/s/r+GI: Poisson
% arrivals at rate LAMBDA, iid general service times of mean 1/MU, S servers, R
% extra waiting spaces and iid patience times with a general distribution.
%
% THE TWO APPROXIMATIONS. The general patience law is replaced by STATE-
% DEPENDENT Markovian abandonment: a customer who is jth from the end of a queue
% abandons at rate delta_j = h(j/lambda), where h = f/(1-F) is the patience
% hazard rate, because a customer in that position has been waiting for about
% j/lambda (eq. 3.3). The total abandonment rate with k waiting is then
% Delta_k = sum_{j<=k} delta_j (eq. 3.4). The general service law is replaced by
% an exponential of the same mean (Section 5), which is accurate here because
% with many servers and non-negligible abandonment the model behaves like a loss
% system, where the service law is insensitive beyond its mean. What is left is
% the Markovian M/M/s/r+M(n) model, solved as a birth-and-death process.
%
% WHAT THE PATIENCE LAW CONTRIBUTES. Only the hazard function NEAR THE ORIGIN,
% not the mean and not the tail: waits are O(1/sqrt(s)) in the many-server
% regime, so a customer either abandons early or never. That is the paper's main
% modelling insight and the reason a two-moment patience fit is not enough.
%
% PATIENCE accepts three forms:
%   scalar THETA     - exponential patience of rate THETA, h(t) = THETA. The
%                      approximations are then EXACT and the model is Erlang A
%                      (see QSYS_ERLANGA).
%   function handle  - the hazard rate h(t), used as in eq. (3.3).
%   struct('ccdf',G) - the complementary cdf G(t) = 1-F(t), used through the
%                      integrated form Delta_k = -log G(k/lambda) of eq. (3.6),
%                      which is the variant to use when the density is not
%                      smooth.
%
% RESULT = QSYS_MGISRGI_WHITT(..., 'wPoints', T) also returns the waiting-time
% cdfs at the times T, obtained by numerically inverting the transforms of
% eqs. (7.22)-(7.23) and (7.32)-(7.33) with the Abate-Whitt EULER algorithm.
% Other options: 'maxQueue' (truncation level used when R is Inf, default
% 100000), 'tol' (relative tail tolerance for that truncation, default 1e-14),
% 'invMethod' and 'invN' (inversion method and node count, default 'euler', 41).
%
% Returns a struct with fields:
%   queueLengthDist  - P(N = k) for k = 0..s+r, N the number in system
%   probLoss         - P(an arrival is blocked) = p_{s+r}, zero when R is Inf
%   probNoWait       - P(W = 0) among entering customers
%   probServed       - P(S), an entering customer is eventually served
%   probAbandon      - P(A) = 1 - P(S)
%   meanNumber       - E[N]         varNumber        - Var[N]
%   meanQueueLength  - E[Q]         varQueueLength   - Var[Q], Q = (N-s)^+
%   meanWaitServed   - E[W|S]       varWaitServed    - Var[W|S]
%   meanWaitAbandon  - E[W|A]       varWaitAbandon   - Var[W|A]
%   meanWait         - E[W] over entering customers, zeros included
%   secondMomentWait - E[W^2] over entering customers
%   utilization      - E[min(N,s)]/s, the fraction of servers busy
%   throughput       - rate of served customers, lambda(1-P_loss)P(S)
%   abandonRate      - rate of abandoning customers, lambda(1-P_loss)P(A)
%   abandonRates     - delta_j, j = 1..r    totalAbandonRates - Delta_k, k = 0..r
%   waitPoints       - the times T, when requested
%   cdfWaitServed    - P(W <= t | S)        cdfWaitAbandon   - P(W <= t | A)
%   cdfWait          - P(W <= t) over entering customers
%
% ACCURACY. Exact for M/M/s/r+M. Elsewhere the paper reports errors of a few
% percent against simulation, degrading as the service SCV moves away from 1.
%
% Example:
%   % M/M/100/200+M with mean patience 1, from Table 1 of the paper
%   res = qsys_mgisrgi_whitt(102, 1/10, 100, 200, 1);
%
% Reference: W. Whitt (2005). Engineering solution of a basic call-center model.
% Management Science 51(2), 221-235.
%
% See also QSYS_ERLANGA, QSYS_GGSGI_FLUID, LAPLACE_INVERT_EULER.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

options = qsys_mgisrgi_whitt_options(varargin);

if lambda <= 0
    line_error(mfilename, 'The arrival rate lambda must be positive.');
end
if mu <= 0
    line_error(mfilename, 'The service rate mu must be positive.');
end
s = round(s);
if s < 1
    line_error(mfilename, 'The number of servers s must be at least 1.');
end
if r < 0
    line_error(mfilename, 'The number of extra waiting spaces r must be non-negative.');
end

[hazardFun, ccdfFun, isExponential, theta] = qsys_mgisrgi_whitt_patience(patience);

% The abandonment rates delta_j and their partial sums Delta_k. Eq. (3.3)-(3.4)
% use the hazard directly; eq. (3.6) integrates it, which is what the ccdf form
% supplies exactly. Delta is stored with a leading Delta_0 = 0, so Dlt(k+1)=Delta_k.
finiteR = isfinite(r);
if finiteR
    rr = round(r);
else
    rr = options.maxQueue;
end

% The birth-death recursion of eqs. (7.4)-(7.7). x is unnormalized with x_s = 1;
% when r is Inf the upward leg stops as soon as the tail is negligible.
xUp = zeros(1, rr+1);       % xUp(k+1) = x_{s+k}, k = 0..rr
xUp(1) = 1;
Dlt = zeros(1, rr+1);       % Dlt(k+1) = Delta_k
delta = zeros(1, rr);       % delta(j) = delta_j
smu = s*mu;
kUsed = rr;
for k = 0:rr-1
    j = k + 1;
    [delta(j), Dlt(j+1)] = qsys_mgisrgi_whitt_rates(j, lambda, Dlt(j), hazardFun, ccdfFun);
    xUp(k+2) = lambda * xUp(k+1) / (smu + Dlt(j+1));
    if ~finiteR && xUp(k+2) < options.tol * max(xUp(1:k+2)) && k >= 1
        kUsed = j;
        break
    end
end
if ~finiteR
    if kUsed == rr && rr > 0
        line_error(mfilename, sprintf(['the queue-length tail is still %g of its peak at the ' ...
            'truncation level %d; with r = Inf the patience law must make the chain ergodic ' ...
            '(raise maxQueue if the model is genuinely that large)'], ...
            xUp(rr+1)/max(xUp), rr));
    end
    xUp = xUp(1:kUsed+1);
    Dlt = Dlt(1:kUsed+1);
    delta = delta(1:kUsed);
    rr = kUsed;
end

% The downward leg, eq. (7.5), over the states 0..s-1 where not all servers are busy.
xDown = zeros(1, s);        % xDown(k) = x_{k-1}, k = 1..s
xk = 1;                     % x_s
for k = s:-1:1
    xk = k * mu * xk / lambda;
    xDown(k) = xk;
end

x = [xDown, xUp];           % x(k+1) = x_k for k = 0..s+rr
p = x / sum(x);
pa = p;
if finiteR
    probLoss = p(end);
else
    probLoss = 0;
end
pa = pa / (1 - probLoss);   % eq. (7.8), the state seen by an ENTERING customer

% Queue-length and occupancy moments, taken directly from the distribution.
kAll = 0:(s+rr);
qAll = max(0, kAll - s);
meanNumber = sum(kAll .* p);
varNumber = sum(((kAll - meanNumber).^2) .* p);
meanQueue = sum(qAll .* p);
varQueue = sum(((qAll - meanQueue).^2) .* p);
utilization = sum(min(kAll, s) .* p) / s;

% Customer experience. Everything below conditions on the state seen at arrival
% and averages over it, exactly as in Sections 7.2-7.4.
probNoWait = sum(pa(1:s));  % eq. (7.9), states 0..s-1

sigma = zeros(1, rr);       % sigma(k) = P(a customer arriving into position k is served)
Msum = zeros(1, rr);        % M_k     = sum_j m_k(j)
Vsum = zeros(1, rr);        % V_k     = sum_j m_k(j)^2
EWa1 = zeros(1, rr);        % E[W_k 1_A]
EWa2 = zeros(1, rr);        % E[W_k^2 1_A]
for k = 1:rr
    [phik, mk] = qsys_mgisrgi_whitt_kernel(k, smu, Dlt, delta);
    sigma(k) = prod(1 - phik);
    Msum(k) = sum(mk);
    Vsum(k) = sum(mk.^2);
    % Eqs. (7.28)-(7.29): abandoning at the jth departure epoch costs the sum of
    % the first j interdeparture times, whose mean and second moment accumulate.
    surv = 1;
    cumM = 0;
    cumV = 0;
    for j = 1:k
        cumM = cumM + mk(j);
        cumV = cumV + mk(j)^2;
        w = surv * phik(j);
        EWa1(k) = EWa1(k) + w * cumM;
        EWa2(k) = EWa2(k) + w * (cumV + cumM^2);
        surv = surv * (1 - phik(j));
    end
end

% Averaging over the arrival state: finding s+k in system puts the arrival in
% position k+1, so the weights are pa_{s+k} for k = 0..r-1.
wArr = pa(s+1:s+rr);        % wArr(k+1) = pa_{s+k}
probServed = probNoWait + sum(wArr .* sigma);
probAbandon = 1 - probServed;
EWS1 = sum(wArr .* sigma .* Msum);                      % eq. (7.16)
EWS2 = sum(wArr .* sigma .* (Vsum + Msum.^2));          % eq. (7.17)
EWA1 = sum(wArr .* EWa1);                               % eq. (7.26)
EWA2 = sum(wArr .* EWa2);                               % eq. (7.27)

result.queueLengthDist = p;
result.probLoss = probLoss;
result.probNoWait = probNoWait;
result.probServed = probServed;
result.probAbandon = probAbandon;
result.meanNumber = meanNumber;
result.varNumber = varNumber;
result.meanQueueLength = meanQueue;
result.varQueueLength = varQueue;
result.utilization = utilization;
result.throughput = lambda * (1 - probLoss) * probServed;
result.abandonRate = lambda * (1 - probLoss) * probAbandon;
result.meanWaitServed = qsys_mgisrgi_whitt_ratio(EWS1, probServed);
result.varWaitServed = max(0, qsys_mgisrgi_whitt_ratio(EWS2, probServed) - result.meanWaitServed^2);
result.meanWaitAbandon = qsys_mgisrgi_whitt_ratio(EWA1, probAbandon);
result.varWaitAbandon = max(0, qsys_mgisrgi_whitt_ratio(EWA2, probAbandon) - result.meanWaitAbandon^2);
result.meanWait = EWS1 + EWA1;
result.secondMomentWait = EWS2 + EWA2;
result.abandonRates = delta;
result.totalAbandonRates = Dlt;
result.numWaitingSpaces = rr;
result.isExponentialPatience = isExponential;
result.patienceRate = theta;

if ~isempty(options.wPoints)
    t = options.wPoints(:).';
    ws = @(z) qsys_mgisrgi_whitt_transform(z, wArr, sigma, smu, Dlt, delta, true);
    ab = @(z) qsys_mgisrgi_whitt_transform(z, wArr, sigma, smu, Dlt, delta, false);
    FS = zeros(1, numel(t));
    FA = zeros(1, numel(t));
    for i = 1:numel(t)
        FS(i) = laplace_invert(@(z) ws(z)/z, t(i), options.invMethod, options.invN);
        FA(i) = laplace_invert(@(z) ab(z)/z, t(i), options.invMethod, options.invN);
    end
    FS = min(max(FS, 0), probServed - probNoWait);
    FA = min(max(FA, 0), probAbandon);
    result.waitPoints = t;
    result.cdfWaitServed = (probNoWait + FS) / max(probServed, realmin);   % eq. (7.24)
    result.cdfWaitAbandon = FA / max(probAbandon, realmin);                % eq. (7.34)
    result.cdfWait = probNoWait + FS + FA;                                 % eq. (7.35)
end
end

function opt = qsys_mgisrgi_whitt_options(args)
% Name-value options with the defaults documented in the header.
opt = struct('wPoints', [], 'maxQueue', 100000, 'tol', 1e-14, 'invMethod', 'euler', 'invN', 41);
for i = 1:2:numel(args)
    name = args{i};
    if i+1 > numel(args)
        line_error(mfilename, sprintf('option %s has no value', char(name)));
    end
    switch lower(char(name))
        case 'wpoints'
            opt.wPoints = args{i+1};
        case 'maxqueue'
            opt.maxQueue = round(args{i+1});
        case 'tol'
            opt.tol = args{i+1};
        case 'invmethod'
            opt.invMethod = args{i+1};
        case 'invn'
            opt.invN = round(args{i+1});
        otherwise
            line_error(mfilename, sprintf('unknown option %s', char(name)));
    end
end
end

function [hazardFun, ccdfFun, isExponential, theta] = qsys_mgisrgi_whitt_patience(patience)
% Resolve the three accepted forms of the patience argument.
hazardFun = [];
ccdfFun = [];
isExponential = false;
theta = NaN;
if isnumeric(patience) && isscalar(patience)
    if patience < 0
        line_error(mfilename, 'the patience rate theta must be non-negative');
    end
    theta = patience;
    isExponential = true;
    hazardFun = @(t) theta * ones(size(t));
elseif isa(patience, 'function_handle')
    hazardFun = patience;
elseif isstruct(patience) && isfield(patience, 'ccdf')
    ccdfFun = patience.ccdf;
elseif isstruct(patience) && isfield(patience, 'hazard')
    hazardFun = patience.hazard;
else
    line_error(mfilename, ['patience must be a scalar rate, a hazard-rate function handle ' ...
        'or a struct with field ccdf or hazard']);
end
end

function [deltaj, Deltaj] = qsys_mgisrgi_whitt_rates(j, lambda, DeltaPrev, hazardFun, ccdfFun)
% One step of eqs. (3.3)-(3.4) (hazard form) or (3.5)-(3.6) (ccdf form).
%
% DIVERGENCE from the printed eqs. (3.5)-(3.6): they read
% delta_j = int_{(j-1)/lambda}^{j/lambda} h(t) dt and Delta_k = -log F^c(k/lambda),
% which are cumulative hazards, i.e. dimensionless, while delta and Delta are
% rates everywhere else in the paper. They are the AVERAGE hazard over an
% interval of length 1/lambda, so the factor lambda is missing. Restoring it
% makes the ccdf form reduce to the exact Erlang A rates under exponential
% patience, which the paper states this approximation does (eq. 7.12); the
% literal form gives theta/lambda instead of theta and is wrong by that factor.
if isempty(ccdfFun)
    deltaj = hazardFun(j/lambda);
    Deltaj = DeltaPrev + deltaj;
else
    g = ccdfFun(j/lambda);
    if g <= 0
        line_error(mfilename, sprintf(['the patience ccdf vanishes at t = %g, so every ' ...
            'customer has abandoned by then; supply a hazard handle instead'], j/lambda));
    end
    Deltaj = -lambda*log(g);
    deltaj = Deltaj - DeltaPrev;
end
if deltaj < 0
    line_error(mfilename, 'the patience law produced a negative abandonment rate');
end
end

function [phik, mk] = qsys_mgisrgi_whitt_kernel(k, smu, Dlt, delta)
% Eqs. (7.10)-(7.11): with k waiting, the total departure rate before the jth
% departure epoch is s*mu + Delta_k - Delta_{j-1}, of which delta_j is the share
% belonging to the customer of interest.
j = 1:k;
rate = smu + Dlt(k+1) - Dlt(j);
mk = 1 ./ rate;
phik = delta(j) ./ rate;
end

function v = qsys_mgisrgi_whitt_ratio(num, den)
% A conditional moment is 0/0 when the conditioning event cannot happen.
if den <= 0
    v = 0;
else
    v = num / den;
end
end

function val = qsys_mgisrgi_whitt_transform(z, wArr, sigma, smu, Dlt, delta, served)
% Eq. (7.22)-(7.23) when SERVED, eq. (7.32)-(7.33) otherwise. Both fold the same
% per-position kernel: the wait is a sum of exponentials with rates
% 1/m_k(j), truncated at the departure epoch that serves or loses the customer.
val = 0;
for k = 1:numel(wArr)
    [phik, mk] = qsys_mgisrgi_whitt_kernel(k, smu, Dlt, delta);
    rate = 1 ./ mk;
    if served
        val = val + wArr(k) * sigma(k) * prod(rate ./ (rate + z));
    else
        surv = 1;
        chain = 1;
        for j = 1:k
            chain = chain * rate(j) / (rate(j) + z);
            val = val + wArr(k) * surv * phik(j) * chain;
            surv = surv * (1 - phik(j));
        end
    end
end
end
