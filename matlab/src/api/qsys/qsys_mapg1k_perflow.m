function result = qsys_mapg1k_perflow(MAPS, svc, K, varargin)
% RESULT = QSYS_MAPG1K_PERFLOW(MAPS, SVC, K)
%
% Per-flow throughput and loss ratio of a FIFO buffer with tail drop that is
% fed by N flows of arbitrary, mutually different statistical character.
%
% Flow n is described by its own MAP, so two flows may share an arrival rate
% and still differ in the shape and autocorrelation of their interarrival
% times. The buffer holds K packets including the one in transmission, and
% the transmission time follows an arbitrary distribution F.
%
%   MAPS - 1 x N cell array, MAPS{n} = {D0n, D1n}, the MAP of flow n. The
%          modulating orders M_n may differ from flow to flow.
%   SVC  - service time descriptor, see QSYS_MAPG1K
%   K    - buffer size in packets, K >= 1
%
% Options are forwarded to QSYS_MAPG1K ('tol', 'nmax').
%
% Returns a struct with fields:
%   throughput(n)      - throughput of flow n [pkts/s]
%   lossRatio(n)       - loss ratio of flow n, in [0,1]
%   lambda(n)          - arrival rate of flow n [pkts/s]
%   lambdaAggregate    - sum_n lambda(n)
%   throughputAggregate- sum_n throughput(n)
%   lossAggregate      - aggregate loss ratio, sum_n L(n)*lambda(n)/lambda
%   p0(n), pK(n)       - empty/full buffer probabilities of the n-th model
%   rho                - offered load lambdaAggregate * E[S]
%
% Method. The exact model of N flows would need a Markov chain tracking the
% modulating state of every flow jointly with the buffer occupancy, hence
% prod_n M_n * (K+1) states; [1] notes this is already out of reach at N=10,
% M_n=3, K=10 (over 4e11 transition matrix entries). Instead, one model per
% flow is solved: flow n is kept exactly as MAP_n, while the other N-1 flows
% are replaced by a single Poisson stream of rate lambda - lambda_n. That
% substitution is justified by the Palm-Khinchin limiting theorem on the
% superposition of many point processes, so it is an approximation that
% improves as N grows, and it is applied N times, once per flow, so no flow
% is ever the one being Poissonized when its own throughput is computed.
% Superposing MAP_n with the Poisson background yields the MAP
%   D0 = D0n - lambdaBar_n*I,   D1 = D1n + lambdaBar_n*I,      ([1], eq. 5)
% which is passed to QSYS_MAPG1K. Each solve costs O((K*M)^3) for the
% embedded chain, so the whole sweep is O(N*(K*M)^3) against the O(M^(3N)*K^3)
% of the exact joint model: linear rather than exponential in the flow count.
%
% Accuracy. [1] reports errors against simulation of the exact model below
% about 8% for N >= 9 with K >= 20, falling to 2.1% at K=50 and 0.5% at
% K=100, and to 1.2% at N=900. Errors are largest when flows are few, highly
% variable, and the buffer is small.
%
% TEST (Table 2 of [1], K=20, rho=1, gamma service with CV=2; the nine flows
% are defined in eq. (55)-(63)). Tables 2-9 of [1] all reproduce to within
% 0.07 pkts/s, the precision to which they are published:
%   throughput ~ [89.8 179.6 269.4 85.2 169.3 253.4 56.9 110.8 166.6] pkts/s
%
% References:
% [1] Chydzinski, A. Per-Flow Throughput of a FIFO Buffer. Applied System
%     Innovation 2026, 9, 112. Theorem 1.
%
% See also QSYS_MAPG1K, QSYS_MG1K_LOSS, QSYS_MAPG1.

if ~iscell(MAPS)
    line_error(mfilename, 'MAPS must be a cell array of {D0,D1} pairs.');
end
N = numel(MAPS);
if N < 1
    line_error(mfilename, 'At least one flow is required.');
end

lam = zeros(1, N);
for n = 1:N
    if ~iscell(MAPS{n}) || numel(MAPS{n}) < 2
        line_error(mfilename, sprintf('MAPS{%d} must be a {D0,D1} cell.', n));
    end
    lam(n) = map_lambda(MAPS{n});
end
lamTot = sum(lam);

Tn = zeros(1, N);
Ln = zeros(1, N);
p0 = zeros(1, N);
pK = zeros(1, N);
Smean = NaN;
for n = 1:N
    D0n = MAPS{n}{1};
    D1n = MAPS{n}{2};
    lamBar = lamTot - lam(n);
    % Superposition of MAP_n with a Poisson background of rate lamBar, eq. (5)
    Mn = size(D0n, 1);
    D0 = D0n - lamBar*eye(Mn);
    D1 = D1n + lamBar*eye(Mn);
    r = qsys_mapg1k(D0, D1, svc, K, varargin{:});
    Smean = r.meanServiceTime;
    p0(n) = r.p0;
    pK(n) = r.pK;
    % Throughput of flow n, eq. (20): the aggregate departure rate (1-p0)/S
    % less the background throughput lamBar*(1-pK), the background loss ratio
    % being pK by PASTA since the background is Poisson.
    Tn(n) = (1 - r.p0)/Smean + r.pK*lamBar - lamBar;
    Ln(n) = 1 - Tn(n)/lam(n);
end

result = struct();
result.throughput = Tn;
result.lossRatio = Ln;
result.lambda = lam;
result.lambdaAggregate = lamTot;
result.throughputAggregate = sum(Tn);
result.lossAggregate = sum(Ln.*lam)/lamTot;
result.p0 = p0;
result.pK = pK;
result.meanServiceTime = Smean;
result.rho = lamTot*Smean;
result.analyzer = 'qsys_mapg1k_perflow';
end
