function result = qsys_mmapg1k(D0, D1c, svc, K, varargin)
% RESULT = QSYS_MMAPG1K(D0, D1C, SVC, K)
%
% Exact per-class throughput and loss ratio of an MMAP[K]/G/1/K queue with
% tail drop: marked Markovian arrivals, arbitrary service time distribution
% F common to all classes, and a finite buffer of K packets (the position
% held by the packet in transmission included).
%
% Two classes of equal arrival rate but different interarrival variability
% or autocorrelation receive different loss ratios, which is the effect that
% motivates [1]. Aggregate-only finite-buffer analyses cannot express it:
% they return a single blocking probability p and set T_k = lambda_k*(1-p),
% making the loss ratio identical across classes by construction.
%
%   D0   - M x M hidden transition matrix of the arrival MMAP
%   D1C  - 1 x R cell array, D1C{k} = M x M arrival matrix of class k.
%          D0 + sum_k D1C{k} must be an irreducible generator.
%   SVC  - service time descriptor, see QSYS_MAPG1K
%   K    - buffer size in packets, K >= 1
%
% Options are forwarded to QSYS_MAPG1K ('tol', 'nmax').
%
% Returns a struct with fields:
%   throughput(k)      - throughput of class k [pkts/s]
%   lossRatio(k)       - loss ratio of class k, in [0,1]
%   lambda(k)          - arrival rate of class k [pkts/s]
%   lambdaAggregate, throughputAggregate, lossAggregate
%   p0, pK             - empty/full buffer probabilities
%   pKvec              - 1 x M, P(buffer full, phase j)
%   plevel             - 1 x (K+1), time-stationary P(level = l)
%   meanQueueLength    - E[number in system]
%   meanServiceTime, rho, utilization
%
% Method. The aggregate MAP {D0, sum_k D1C{k}} drives QSYS_MAPG1K, whose
% embedded chain returns the joint law of buffer level and MAP phase. A
% class-k arrival leaves phase i at rate (D1C{k}*e)_i, so the rate of class-k
% arrivals that meet a full buffer is pKvec*D1C{k}*e, and
%   lambda_k = pi*D1C{k}*e,   L_k = (pKvec*D1C{k}*e)/lambda_k,
% with pi the stationary phase law of the aggregate. This is exact: no
% independence between classes is assumed and no PASTA argument is used, the
% phase resolution of pKvec doing the work instead.
%
% Relation to [1]. Reference [1] instead keeps one flow exact and replaces
% the rest by a Poisson stream of the same rate, invoking Palm-Khinchin, and
% repeats that once per flow. That approximation is needed only when the
% joint arrival process is unavailable, its exact form costing prod_n M_n
% phases for N independent flows. When the joint MMAP is already at hand, as
% it is inside a solver that propagates MMAPs between stations, the phase
% resolution above is both exact and cheaper. Use QSYS_MAPG1K_PERFLOW for
% the setting of [1], where flows are given as N separate MAPs.
%
% Assumes a single server and a service law that is iid and independent of
% class. Per-class service makes the departure rate depend on which class
% holds the server, which this model does not represent.
%
% References:
% [1] Chydzinski, A. Per-Flow Throughput of a FIFO Buffer. Applied System
%     Innovation 2026, 9, 112.
%
% See also QSYS_MAPG1K, QSYS_MAPG1K_PERFLOW, QSYS_MG1K_LOSS.

if ~iscell(D1c)
    line_error(mfilename, 'D1C must be a cell array of per-class D1 matrices.');
end
R = numel(D1c);
M = size(D0, 1);
D1 = zeros(M);
for k = 1:R
    if any(size(D1c{k}) ~= [M M])
        line_error(mfilename, sprintf('D1C{%d} must be %dx%d.', k, M, M));
    end
    D1 = D1 + D1c{k};
end

r = qsys_mapg1k(D0, D1, svc, K, varargin{:});

e = ones(M, 1);
pit = map_prob({D0, D1});
pit = pit(:).';

lam = zeros(1, R);
Lk = zeros(1, R);
Tk = zeros(1, R);
for k = 1:R
    lam(k) = pit*D1c{k}*e;
    if lam(k) > 0
        Lk(k) = (r.pKvec*D1c{k}*e)/lam(k);
    else
        Lk(k) = 0;
    end
    Tk(k) = lam(k)*(1 - Lk(k));
end

result = struct();
result.throughput = Tk;
result.lossRatio = Lk;
result.lambda = lam;
result.lambdaAggregate = sum(lam);
result.throughputAggregate = sum(Tk);
result.lossAggregate = r.lossProbability;
result.p0 = r.p0;
result.pK = r.pK;
result.pKvec = r.pKvec;
result.plevel = r.plevel;
result.meanQueueLength = r.meanQueueLength;
result.meanServiceTime = r.meanServiceTime;
result.utilization = r.utilization;
result.rho = r.rho;
result.analyzer = 'qsys_mmapg1k';
end
