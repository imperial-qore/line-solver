function [Q,U,R,T,C,X,lG,runtime,iter,method] = solver_nc_lossn_analyzer(sn, options)
% SOLVER_NC_LOSSN_ANALYZER Analyzes open loss networks with FCR
%
% This analyzer handles open queueing networks with a single multiclass
% Delay node inside a Finite Capacity Region (FCR) with DROP policy.
%
% The FCR admission rule is A n <= C on the per-class occupancy vector n of
% the region, where the rows of A are assembled from every constraint the
% region declares: the global job cap, the memory budget weighted by the
% per-class sizes, the per-class job caps, and any explicit linear constraint
% set with FiniteCapacityRegion.setConstraint. Rows left unbounded are
% dropped rather than given a surrogate capacity.
%
% Method selection (options.method):
%   'exact' (default)    - Manjunath-Sikdar transform (lossn_manjunath): the
%                          normalization constant is obtained exactly as a
%                          multidimensional contour integral evaluated by
%                          residues. Requires integer A and C.
%   'rec'                - MDD-rec (lossn_rec): the same constant as the exact
%                          sum over the admissible set, obtained by one
%                          memoised walk of the decision diagram holding it.
%                          Places no integrality demand on A or C.
%   'erlangfp'           - Erlang fixed-point (reduced-load) approximation.
%   'mci'                - Monte Carlo importance-sampling summation
%                          (Ross-Wang 1992): estimates the normalization
%                          constant g(C) and class blocking with confidence
%                          intervals. Set options.samples and options.seed.
%
% The default is the residue transform on an integral region and MDD-REC on a
% fractional one. It used to fall back to 'erlangfp' there, an approximation,
% because the residue argument counts whole units; MDD-rec needs only that the
% admissible set be finite and bounded per coordinate, which it still is, so
% the fractional case is now exact as well.

Tstart = tic;
K = sn.nclasses;  % number of classes
M = sn.nstations;

line_debug('NC loss network analyzer starting: method=%s, nstations=%d, nclasses=%d', options.method, M, K);

% 1. Locate the delay station inside the region
regionMatrix = sn.region{1};
stationsInFCR = find(any(regionMatrix(:,1:end-1) >= 0, 2) | regionMatrix(:,end) >= 0);
delayIdx = stationsInFCR(1);

% 2. Offered load per class. A route carries load nu_r = arrival rate times
% mean holding time inside the region, i.e. the visit ratio at the delay
% divided by its service rate. Passing the bare arrival rate would be correct
% only for unit mean service times.
nu = zeros(1, K);
for r = 1:K
    sourceIdx = sn.refstat(r);
    lambda_r = sn.rates(sourceIdx, r);
    c = find(sn.chains(:,r));
    V_r = 1;
    if ~isempty(c)
        % visits is indexed by STATEFUL node, refstat and delayIdx by station
        vref = sn.visits{c}(sn.stationToStateful(sn.refstat(r)), r);
        if vref > 0
            V_r = sn.visits{c}(sn.stationToStateful(delayIdx), r) / vref;
        end
    end
    mu_r = sn.rates(delayIdx, r);
    nu(r) = lambda_r * V_r / mu_r;
end

% 3. Assemble the admission constraints A n <= C_vec of the region
[A, C_vec] = lossn_region_constraints(sn, 1, delayIdx, K);

% 4. Select method and solve
tokens = strsplit(lower(options.method), '.');
isIntegral = all(abs(A(:) - round(A(:))) < 1e-9) && all(abs(C_vec - round(C_vec)) < 1e-9);

if any(strcmp(tokens, 'mci'))
    chosen = 'mci';
elseif any(strcmp(tokens, 'erlangfp'))
    chosen = 'erlangfp';
elseif any(strcmp(tokens, 'rec'))
    chosen = 'rec';
elseif any(strcmp(tokens, 'exact')) || any(strcmp(tokens, 'manjunath')) || any(strcmp(tokens, 'ms'))
    chosen = 'exact';
else
    chosen = 'exact';
    if ~isIntegral
        % the residue argument counts whole units; MDD-rec does not
        chosen = 'rec';
    end
end

lG = NaN;  % normalization constant (finite for exact and mci)
switch chosen
    case 'mci'
        mciopt = struct();
        if isfield(options, 'samples') && ~isempty(options.samples) && ~isinf(options.samples)
            mciopt.samples = options.samples;
        end
        if isfield(options, 'seed') && ~isempty(options.seed)
            mciopt.seed = options.seed;
        end
        [QLen, Loss, lG, ~, niter] = lossn_mci(nu, A, C_vec, mciopt);
        method = 'lossn.mci';
    case 'rec'
        [QLen, Loss, lG, niter] = lossn_rec(nu, A, C_vec);
        method = 'lossn.rec';
    case 'exact'
        [QLen, Loss, lG, niter] = lossn_manjunath(nu, A, C_vec);
        method = 'lossn.exact';
    otherwise
        [QLen, Loss, ~, niter] = lossn_erlangfp(nu, A, C_vec);
        method = 'lossn.erlangfp';
end

% 5. Convert to standard outputs. QLen is the carried load E[n_r], so the
% carried throughput follows from Little's law at the infinite server.
Q = zeros(M, K);
U = zeros(M, K);
T = zeros(M, K);
R = zeros(M, K);
Xc = zeros(1, K);

for r = 1:K
    sourceIdx = sn.refstat(r);
    lambda_r = sn.rates(sourceIdx, r);
    Xc(r) = lambda_r * (1 - Loss(r));       % carried (accepted) rate
    mu_r = sn.rates(delayIdx, r);
    T(delayIdx, r) = Xc(r);
    % The source emits the accepted (post-drop) rate so that the routing-based
    % arrival-rate computation (sn_get_arvr_from_tput) yields a non-zero ArvR
    % at the delay, consistent with the flow-conserving departure throughput
    % and with the simulated rate reported by SolverJMT.
    T(sourceIdx, r) = Xc(r);
    Q(delayIdx, r) = QLen(r);               % mean number in the region
    R(delayIdx, r) = 1 / mu_r;              % response time = service time (IS)
    U(delayIdx, r) = QLen(r);               % IS "utilization" is the mean busy servers
end

X = Xc;                 % system throughput per class (row vector)
C = zeros(1, K);        % cycle time not applicable
iter = niter;
runtime = toc(Tstart);
end

% ------------------------------------------------------------------------

function [A, C_vec] = lossn_region_constraints(sn, f, stationIdx, K)
% Rows of the admission rule A n <= C for region f, in the order in which the
% simulation engines test them. Every row is a function of the per-class
% occupancy of the region only, which is what the FiniteCapacityRegion API
% can express, so no row can distinguish stations inside the region.

regionMatrix = sn.region{f};
A = zeros(0, K);
C_vec = zeros(0, 1);

% Global job cap: sum_r n_r <= globalMaxJobs
globalMax = regionMatrix(stationIdx, K+1);
if globalMax >= 0
    A(end+1, :) = ones(1, K);
    C_vec(end+1, 1) = globalMax;
end

% Memory budget: sum_r classSize_r n_r <= globalMaxMemory
if isfield(sn, 'regionmaxmem') && numel(sn.regionmaxmem) >= f && ~isempty(sn.regionmaxmem{f})
    memMatrix = sn.regionmaxmem{f};
    maxmem = memMatrix(stationIdx);
    if maxmem >= 0
        sz = ones(1, K);
        if isfield(sn, 'regionsz') && ~isempty(sn.regionsz)
            sz = sn.regionsz(f, 1:K);
        end
        A(end+1, :) = sz;
        C_vec(end+1, 1) = maxmem;
    end
end

% Per-class job caps: n_r <= classMaxJobs_r (already folded with classMaxMemory)
classMax = regionMatrix(stationIdx, 1:K);
for r = 1:K
    if classMax(r) >= 0
        row = zeros(1, K);
        row(r) = 1;
        A(end+1, :) = row;
        C_vec(end+1, 1) = classMax(r);
    end
end

% Explicit linear constraints from FiniteCapacityRegion.setConstraint
if isfield(sn, 'regionlincon') && size(sn.regionlincon, 1) >= f && ~isempty(sn.regionlincon{f,1})
    linA = sn.regionlincon{f,1};
    linB = sn.regionlincon{f,2};
    for k = 1:size(linA, 1)
        A(end+1, :) = linA(k, 1:K);
        C_vec(end+1, 1) = linB(k);
    end
end

if isempty(C_vec)
    line_error(mfilename, 'Finite capacity region declares no bounded constraint.');
end
end
