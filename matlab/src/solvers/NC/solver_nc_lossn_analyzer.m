function [Q,U,R,T,C,X,lG,runtime,iter,method] = solver_nc_lossn_analyzer(sn, options)
% SOLVER_NC_LOSSN_ANALYZER Analyzes open loss networks with FCR
%
% This analyzer handles open queueing networks with a single multiclass
% Delay node inside a Finite Capacity Region (FCR) with DROP policy.
%
% Method selection (options.method):
%   'erlangfp' (default) - Erlang fixed-point (reduced-load) approximation.
%   'mci'                - Monte Carlo importance-sampling summation
%                          (Ross-Wang 1992): estimates the normalization
%                          constant g(C) and class blocking with confidence
%                          intervals. Set options.samples and options.seed.

Tstart = tic;
K = sn.nclasses;  % number of classes
M = sn.nstations;

line_debug('NC loss network analyzer starting: method=%s, nstations=%d, nclasses=%d', options.method, M, K);

% 1. Extract arrival rates from Source
nu = zeros(1, K);
for r = 1:K
    sourceIdx = sn.refstat(r);
    nu(r) = sn.rates(sourceIdx, r);  % arrival rate
end

% 2. Find delay station in FCR and extract constraints
regionMatrix = sn.region{1};
stationsInFCR = find(any(regionMatrix(:,1:end-1) >= 0, 2) | regionMatrix(:,end) >= 0);
delayIdx = stationsInFCR(1);

globalMax = regionMatrix(delayIdx, K+1);
classMax = regionMatrix(delayIdx, 1:K);

% Handle unbounded constraints (replace -1 with large value for Erlang)
if globalMax < 0
    globalMax = 1e6;
end
classMax(classMax < 0) = 1e6;

% 3. Build A matrix (J x K) where J = K+1 links
% Link 1: global constraint (all classes contribute)
% Links 2..K+1: per-class constraints (only class r contributes to link r+1)
J = K + 1;
A = zeros(J, K);
A(1, :) = 1;  % global link: all classes contribute
for r = 1:K
    A(r+1, r) = 1;  % per-class link: only class r contributes
end

% 4. Build C vector (J x 1)
C_vec = zeros(J, 1);
C_vec(1) = globalMax;
C_vec(2:end) = classMax(:);

% 5. Select method and solve
tokens = strsplit(lower(options.method), '.');
useMCI = any(strcmp(tokens, 'mci'));

lG = NaN;  % normalization constant (finite only for mci)
if useMCI
    mciopt = struct();
    if isfield(options, 'samples') && ~isempty(options.samples) && ~isinf(options.samples)
        mciopt.samples = options.samples;
    end
    if isfield(options, 'seed') && ~isempty(options.seed)
        mciopt.seed = options.seed;
    end
    [QLen, Loss, lG, ~, niter] = lossn_mci(nu, A, C_vec, mciopt);
    method = 'mci';
else
    [QLen, Loss, E, niter] = lossn_erlangfp(nu, A, C_vec);
    method = 'erlangfp';
end

% 6. Convert to standard outputs
Q = zeros(M, K);
U = zeros(M, K);
T = zeros(M, K);
R = zeros(M, K);

% At delay node: QLen is effective throughput (after loss)
for r = 1:K
    T(delayIdx, r) = QLen(r);  % effective throughput = arrival rate * (1-loss)
    % Source emits the accepted (post-drop) rate into the network so that the
    % routing-based arrival-rate computation (sn_get_arvr_from_tput) yields a
    % non-zero ArvR at the delay, consistent with the flow-conserving departure
    % throughput and with the simulated rate reported by SolverJMT.
    sourceIdx = sn.refstat(r);
    T(sourceIdx, r) = QLen(r);
    mu_r = sn.rates(delayIdx, r);  % service rate at delay
    Q(delayIdx, r) = QLen(r) / mu_r;  % Little's law: Q = X * S
    R(delayIdx, r) = 1 / mu_r;  % response time = service time (infinite server)
    U(delayIdx, r) = T(delayIdx, r) / mu_r;  % "utilization" for delay
end

X = QLen(:)';  % system throughput per class (row vector)
C = zeros(1, K);  % cycle time not applicable
iter = niter;
runtime = toc(Tstart);
end
