function result = qsys_mapg1k(D0, D1, svc, K, varargin)
% RESULT = QSYS_MAPG1K(D0, D1, SVC, K)
%
% Exact analysis of a MAP/G/1/K queue with tail drop: Markovian arrivals,
% arbitrary service time distribution F, and a finite buffer of K packets
% (the position held by the packet in transmission included).
%
% Unlike QSYS_MAPG1, the service time is NOT fitted to a phase-type
% distribution: F enters exactly, through the functionals A_m and Q_m
% evaluated by uniformization of the arrival MAP. Unlike QSYS_MG1K_LOSS,
% which embeds the same way but assumes Poisson input, arrivals may be a
% general MAP, so flows with equal rate but different interarrival
% variability or autocorrelation are told apart.
%
%   D0, D1 - MAP parameter matrices (M x M), D0 + D1 an irreducible generator
%   SVC    - service time descriptor, a struct with field 'type':
%              'gamma'   : fields alpha (shape), theta (scale). Covers Exp
%                          (alpha=1) and Erlang (alpha integer).
%              'det'     : field d (constant service time)
%              'ph'      : fields alpha (1 x p row), T (p x p subgenerator)
%              'density' : field pdf (handle), optional field tmax
%   K      - buffer size in packets, K >= 1
%
% RESULT = QSYS_MAPG1K(..., 'tol', TOL) sets the uniformization truncation
% tolerance (default 1e-12). RESULT = QSYS_MAPG1K(..., 'nmax', N) caps the
% uniformization order.
%
% Returns a struct with fields:
%   p0              - stationary probability of an empty buffer
%   pK              - stationary probability of a full buffer
%   lossProbability - loss ratio of the aggregate arrival stream, 1-T/lambda
%   throughput      - aggregate throughput [pkts/s]
%   lambda          - aggregate arrival rate of the MAP
%   meanServiceTime - S = E[service time]
%   utilization     - 1 - p0
%   rho             - offered load lambda*S
%   nmax            - uniformization order used
%   sigma           - stationary law of the embedded chain, K*M entries
%   pKvec           - 1 x M, P(buffer full, phase j), summing to pK
%   p0vec           - 1 x M, P(buffer empty, phase j), summing to p0
%   plevel          - 1 x (K+1), time-stationary P(level = l), l = 0..K
%   meanQueueLength - E[number in system], sum_l l*plevel(l+1)
%
% Method. The chain embedded at departure epochs is used, in the state
% (n,j): n = 0..K-1 packets left behind by a departure, j = MAP phase. With
% A_m the matrix of "m arrivals during a service, phase i -> j",
%   n >= 1:  n' = n-1+min(m, K-n),  via A_m, overflow sum_{m>=K-n} A_m
%   n == 0:  the phase first jumps by (-D0)^{-1}*D1 (the idle period ends at
%            an arrival), and the service then proceeds as from n = 1.
% Its stationary law sigma gives, by Markov renewal reward,
%   E[cycle] = S + sum_j sigma(0,j)*idle_j,   idle = (-D0)^{-1}*e
%   T  = 1/E[cycle],   p0 = (sum_j sigma(0,j)*idle_j)/E[cycle] = 1 - T*S
%   pK = E[time at level K per cycle]/E[cycle], from Q_m,
% where Q_m is the expected time within a service with exactly m arrivals so
% far. Time-stationary p0 and pK follow, so no PASTA assumption is needed on
% the MAP side.
%
% This is not the transform solution of Theorem 1 of [1], which is stated in
% terms of a sequence R_m obeying R(z) = z*(A(z)-z*I)^{-1}. That sequence
% grows geometrically, at a rate set by the smallest zero of det(A(z)-z*I),
% while the quantity extracted from it stays O(s) as s -> 0+; cond(G(s))
% therefore grows like that ratio^K and crosses the double-precision ceiling
% near K = 20 for the flows of [1] under gamma service with CV = 2, and near
% K = 10 under constant service, where A_0 = exp(D0*d) has entries O(1e-8).
% Reference [1] evaluates its formulae in arbitrary precision, so the
% restriction is invisible there. The embedded chain used here has every
% entry a probability or a time and is stable for any K and any F.
%
% TEST (M/M/1/5, exact loss 0.04812030):
%  r=qsys_mapg1k(-2,2,struct('type','gamma','alpha',1,'theta',1/3),5); r.pK
%
% References:
% [1] Chydzinski, A. Per-Flow Throughput of a FIFO Buffer. Applied System
%     Innovation 2026, 9, 112.
% [2] Niu, Z.; Cooper, R.B. Transform-Free Analysis of M/G/1/K and Related
%     Queues. Mathematics of Operations Research 1993, 18, 486-510.
%
% See also QSYS_MAPG1K_PERFLOW, QSYS_MAPG1, QSYS_MG1K_LOSS, QSYS_MMCK.

p = inputParser;
addParameter(p, 'tol', 1e-12);
addParameter(p, 'nmax', 200000);
parse(p, varargin{:});
tol = p.Results.tol;
nmaxCap = p.Results.nmax;

M = size(D0, 1);
if size(D0, 2) ~= M || any(size(D1) ~= [M M])
    line_error(mfilename, 'D0 and D1 must be square matrices of equal size.');
end
if K < 1 || K ~= round(K)
    line_error(mfilename, 'Buffer size K must be a positive integer.');
end
beta = -diag(D0);
if any(beta <= 0)
    line_error(mfilename, 'D0 must have strictly negative diagonal entries.');
end

% Uniformization constant: theta >= max_i beta_i keeps I+D0/theta substochastic
theta = max(beta);

% c_n = E[exp(-theta*S)*(theta*S)^n/n!], summing to f(0) = 1
[cn, Smean] = i_service(svc, theta, tol, nmaxCap);
nmax = numel(cn) - 1;

% see _kb/03-api-layer.md (qsys/ family) for rationale
tailc = flipud(cumsum(flipud(cn(:))));
dn = [tailc(2:end); 0]/theta;

% A_m and Q_m for m = 0..K-1, plus B0 = sum_m A_m = E[exp((D0+D1)*S)]
mmax = max(K-1, 0);
A = zeros(M, M, mmax+1);
Q = zeros(M, M, mmax+1);
Sn = zeros(M, M, mmax+1);
Sn(:,:,1) = eye(M);
B0 = zeros(M);
Qtot = zeros(M);
Pn = eye(M);
Pt0 = eye(M) + D0/theta;
Pt1 = D1/theta;
PD = eye(M) + (D0+D1)/theta;
for n = 0:nmax
    for m = 0:min(n, mmax)
        A(:,:,m+1) = A(:,:,m+1) + Sn(:,:,m+1)*cn(n+1);
        Q(:,:,m+1) = Q(:,:,m+1) + Sn(:,:,m+1)*dn(n+1);
    end
    B0 = B0 + Pn*cn(n+1);
    Qtot = Qtot + Pn*dn(n+1);       % sum_m Q_m = int_0^inf exp(D*x)*(1-F(x))dx
    if n < nmax
        Snew = zeros(M, M, mmax+1);
        for m = 0:min(n+1, mmax)
            acc = zeros(M);
            if m <= n
                acc = acc + Sn(:,:,m+1)*Pt0;
            end
            if m >= 1 && m-1 <= n
                acc = acc + Sn(:,:,m)*Pt1;
            end
            Snew(:,:,m+1) = acc;
        end
        Sn = Snew;
        Pn = Pn*PD;
    end
end

e = ones(M, 1);
negD0inv = inv(-D0);
Psi = negD0inv*D1;            % phase at the arrival that ends an idle period
idle = negD0inv*e;            % expected idle time from each phase

% Embedded chain at departure epochs, state (n,j) -> index n*M+j
P = zeros(K*M, K*M);
lastblk = (K-1)*M + (1:M);
for n = 1:K-1
    rows = n*M + (1:M);
    Bacc = B0;
    for m = 0:K-n-1
        P(rows, (n-1+m)*M + (1:M)) = P(rows, (n-1+m)*M + (1:M)) + A(:,:,m+1);
        Bacc = Bacc - A(:,:,m+1);
    end
    % Bacc = sum_{m>=K-n} A_m: every further arrival overflows the buffer
    P(rows, lastblk) = P(rows, lastblk) + Bacc;
end
rows = 1:M;
Bacc = B0;
for m = 0:K-2
    P(rows, m*M + (1:M)) = P(rows, m*M + (1:M)) + Psi*A(:,:,m+1);
    Bacc = Bacc - A(:,:,m+1);
end
P(rows, lastblk) = P(rows, lastblk) + Psi*Bacc;

rowdev = max(abs(sum(P, 2) - 1));
if rowdev > 1e-8
    line_error(mfilename, sprintf(['Embedded chain rows deviate from 1 by %.2e. ' ...
        'The uniformization series for A_m has not converged; raise ''nmax''.'], rowdev));
end

sigma = dtmc_solve(P);
sigma = sigma(:).';
sigma0 = sigma(1:M);

% Markov renewal reward over the interval between successive departures
idleTime = sigma0*idle;
Ecyc = Smean + idleTime;
T = 1/Ecyc;
p0 = idleTime/Ecyc;

% see _kb/03-api-layer.md (qsys/ family) for rationale
Qcum = zeros(M, M, mmax+1);
acc = zeros(M);
for m = 0:mmax
    acc = acc + Q(:,:,m+1);
    Qcum(:,:,m+1) = acc;
end
timeKvec = zeros(1, M);
for n = 1:K-1
    r = K-n-1;                       % Qcum(:,:,r+1) = sum_{m=0}^{r} Q_m
    timeKvec = timeKvec + sigma(n*M + (1:M))*(Qtot - Qcum(:,:,r+1));
end
if K >= 2
    timeKvec = timeKvec + sigma0*Psi*(Qtot - Qcum(:,:,K-1));
else
    timeKvec = timeKvec + sigma0*Psi*Qtot;
end
pKvec = timeKvec/Ecyc;
pK = sum(pKvec);

% see _kb/03-api-layer.md (qsys/ family) for rationale
timeL = zeros(1, K+1);
timeL(1) = idleTime;
for n = 1:K-1
    sn_row = sigma(n*M + (1:M));
    for l = n:K-1
        timeL(l+1) = timeL(l+1) + sn_row*Q(:,:,l-n+1)*e;
    end
end
s0Psi = sigma0*Psi;
for l = 1:K-1
    timeL(l+1) = timeL(l+1) + s0Psi*Q(:,:,l)*e;
end
timeL(K+1) = sum(timeKvec);
plevel = timeL/Ecyc;
massdev = abs(sum(plevel) - 1);
if massdev > 1e-8
    line_error(mfilename, sprintf(['Level distribution has mass %.12f. The Q_m ' ...
        'series has not converged; raise ''nmax''.'], sum(plevel)));
end
meanQ = (0:K)*plevel(:);
% Time at level 0 is the idle period alone, whose phase law is (-D0)^{-1}
p0vec = (sigma0*negD0inv)/Ecyc;

lambda = map_lambda({D0, D1});

result = struct();
result.p0 = p0;
result.pK = pK;
result.throughput = T;
result.lossProbability = 1 - T/lambda;
result.lambda = lambda;
result.meanServiceTime = Smean;
result.utilization = 1 - p0;
result.rho = lambda*Smean;
result.nmax = nmax;
result.sigma = sigma;
result.pKvec = pKvec;
result.p0vec = p0vec;
result.plevel = plevel;
result.meanQueueLength = meanQ;
result.analyzer = 'qsys_mapg1k';
end

% -------------------------------------------------------------------------
function Smean = i_svcmean(svc)
% Mean service time S of the descriptor.
switch lower(svc.type)
    case 'gamma'
        Smean = svc.alpha*svc.theta;
    case 'det'
        Smean = svc.d;
    case 'ph'
        Smean = -svc.alpha(:).'*(svc.T\ones(size(svc.T,1),1));
    case 'density'
        if isfield(svc, 'tmax'); tmax = svc.tmax; else; tmax = Inf; end
        Smean = i_dquad(@(x) x, svc.pdf, tmax);
    otherwise
        line_error(mfilename, sprintf('Unsupported service type ''%s''.', svc.type));
end
end

% -------------------------------------------------------------------------
function [cn, Smean] = i_service(svc, theta, tol, nmaxCap)
% c_n = E[exp(-theta*S)*(theta*S)^n/n!] for n = 0..nmax, and the mean S.
% sum_{n>=0} c_n = E[exp(-theta*S)*exp(theta*S)] = 1 exactly, which both
% sets the truncation order and certifies it.
if ~isfield(svc, 'type')
    line_error(mfilename, 'Service descriptor must have a ''type'' field.');
end
Smean = i_svcmean(svc);
switch lower(svc.type)
    case 'gamma'
        al = svc.alpha; th = svc.theta;
        fn = @(nn) exp(nn*log(theta*th) - gammaln(nn+1) + gammaln(al+nn) ...
            - gammaln(al) - (al+nn)*log1p(th*theta));
    case 'det'
        d = svc.d;
        fn = @(nn) exp(-theta*d + nn*log(theta*d) - gammaln(nn+1));
    case 'ph'
        alv = svc.alpha(:).'; T = svc.T;
        tv = -T*ones(size(T,1), 1);
        % c_n = theta^n * alpha * (theta*I-T)^{-(n+1)} * t
        Minv = inv(theta*eye(size(T)) - T);
        fn = @(nn) i_phblk(nn, alv, Minv, tv, theta);
    case 'density'
        if isfield(svc, 'tmax'); tmax = svc.tmax; else; tmax = Inf; end
        fn = @(nn) i_cquad(svc.pdf, theta, nn, tmax);
    otherwise
        line_error(mfilename, sprintf('Unsupported service type ''%s''.', svc.type));
end
cn = i_grow(fn, tol, nmaxCap, i_guess(theta*Smean, nmaxCap));
if abs(1 - sum(cn)) > 1e-6
    line_warning(mfilename, sprintf(['Uniformization series truncated at n=%d with ' ...
        'residual %g; increase ''nmax''.'], numel(cn)-1, abs(1 - sum(cn))));
end
end

% -------------------------------------------------------------------------
function cn = i_grow(fn, tol, cap, n0)
% Build c_0..c_N in blocks, stopping when the series sums to 1 within tol or
% when a whole block adds nothing in floating point, i.e. the representable
% series is exhausted. The second criterion terminates paths whose terms are
% known only to quadrature accuracy, where the first can never be met.
cn = fn((0:n0).');
while numel(cn) < cap
    if abs(1 - sum(cn)) <= tol
        break
    end
    n = numel(cn);
    add = fn((n:min(cap-1, n + 63)).');
    cn = [cn; add];
    if sum(add) <= eps*sum(cn)
        break
    end
end
cn = cn(:);
end

% -------------------------------------------------------------------------
function v = i_phblk(nn, alv, Minv, tv, theta)
v = zeros(numel(nn), 1);
for i = 1:numel(nn)
    v(i) = theta^nn(i)*(alv*(Minv^(nn(i)+1))*tv);
end
end

% -------------------------------------------------------------------------
function v = i_cquad(pdf, theta, nn, tmax)
v = zeros(numel(nn), 1);
for i = 1:numel(nn)
    n = nn(i);
    v(i) = i_dquad(@(x) exp(-theta*x + n*log(theta*x) - gammaln(n+1)), pdf, tmax);
end
end

% -------------------------------------------------------------------------
function v = i_dquad(w, pdf, tmax)
% E[w(S)] for a service law given by a density, under the substitution
% x = exp(u). An integrable density may diverge at the origin (the gamma
% density behaves as x^(alpha-1), i.e. x^(-0.75) at the CV=2 shape used in
% [1]), which caps adaptive quadrature on [0,tmax] at a few digits. The
% Jacobian exp(u) turns x^(alpha-1)dx into exp(alpha*u)du, which decays
% smoothly as u -> -Inf for any alpha > 0, so the singularity disappears
% rather than being resolved.
%
% The transformed integrand tends to 0 at both ends: w is bounded and
% integrability of pdf forces x*pdf(x) -> 0 as x -> 0 and as x -> Inf. In
% floating point those limits are reached as 0*Inf, since exp(u) underflows
% while pdf(exp(u)) overflows, so the NaN produced there is an artifact of
% the substitution and is replaced by the analytic limit.
if isfinite(tmax)
    ulim = log(tmax);
else
    ulim = Inf;
end
v = integral(@(u) i_finite(w(exp(u)).*pdf(exp(u)).*exp(u)), -Inf, ulim, ...
    'AbsTol', 1e-300, 'RelTol', 1e-13);
end

% -------------------------------------------------------------------------
function y = i_finite(y)
y(~isfinite(y)) = 0;
end

% -------------------------------------------------------------------------
function n = i_guess(m, cap)
% Initial uniformization order: mean plus a generous deviation allowance.
n = min(cap, max(32, ceil(m + 10*sqrt(max(m, 1)) + 32)));
end
