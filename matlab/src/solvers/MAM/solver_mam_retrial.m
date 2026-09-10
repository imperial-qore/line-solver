function [QN,UN,RN,TN,CN,XN,totiter,perf] = solver_mam_retrial(sn, options)
% [QN,UN,RN,TN,CN,XN,TOTITER] = SOLVER_MAM_RETRIAL(SN, OPTIONS)
%
% Solves queueing models with customer impatience:
%
% 1. RETRIAL (orbit impatience): BMAP/PH/N/N bufferless retrial queues
%    Reference: Dudin et al., "Analysis of BMAP/PH/N-Type Queueing System
%    with Flexible Retrials Admission Control", Mathematics 2025, 13(9), 1434.
%
% 2. RENEGING (queue abandonment): MAP/M/s+G queues with patience
%    Reference: O. Gursoy, K. A. Mehr, N. Akar, "The MAP/M/s + G Call Center
%    Model with Generally Distributed Patience Times"
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Check for reneging (queue abandonment) first. The shape predicate is the
% one MAM_RETRIAL_APPLICABLE and SolverMAM.supportsModelMethod ask as well.
[isReneging, renegingInfo] = mam_reneging_applicable(sn);
if isReneging
    [QN,UN,RN,TN,CN,XN,totiter] = solveReneging(sn, options, renegingInfo);
    perf = struct('analyzer','LINE:solver_mam_reneging');
    return;
end

% Check for retrial topology
[isRetrial, retInfo] = qsys_is_retrial(sn);
if ~isRetrial
    line_error(mfilename, 'No valid impatience configuration detected (retrial or reneging).');
end

% Initialize output arrays
M = sn.nstations;
K = sn.nclasses;
QN = zeros(M, K);
UN = zeros(M, K);
RN = zeros(M, K);
TN = zeros(M, K);
CN = zeros(1, K);
XN = zeros(1, K);

% Extract indices
sourceIdx = retInfo.sourceIdx;
queueIdx = retInfo.stationIdx;
classIdx = retInfo.classIdx;
N = retInfo.N;

% The matrix-analytic retrial engine is defined only for Markovian arrival
% and phase-type service. Reject anything else here, naming the station and
% class, rather than letting a malformed representation reach the generator.
assertMarkovian(sn, sourceIdx, classIdx, 'arrival');
assertMarkovian(sn, queueIdx, classIdx, 'service');
warnIfApproximated(sn, queueIdx, classIdx);

% Extract arrival process (BMAP) from source
% sn.proc{sourceIdx}{classIdx} contains the arrival process
arrivalProc = sn.proc{sourceIdx}{classIdx};

% Convert arrival process to BMAP matrices D = {D0, D1, ...}
D = extractBMAPMatrices(arrivalProc);

if isempty(D)
    line_error(mfilename, 'Could not extract BMAP matrices from arrival process.');
end

% Extract PH service distribution from queue
% sn.proc{queueIdx}{classIdx} contains {alpha, A} for PH
serviceProc = sn.proc{queueIdx}{classIdx};
[beta, S] = extractPHParams(serviceProc);

if isempty(beta) || isempty(S)
    line_error(mfilename, 'Could not extract PH parameters from service process.');
end

% Extract retrial rate alpha from the canonical sn.retrialProc {D0,D1} field
alpha = 0.1;  % Default retrial rate
if isfield(sn, 'retrialProc') && ~isempty(sn.retrialProc)
    if size(sn.retrialProc, 1) >= queueIdx && size(sn.retrialProc, 2) >= classIdx
        retrialDist = sn.retrialProc{queueIdx, classIdx};
        if ~isempty(retrialDist) && iscell(retrialDist)
            % For Exp(alpha), D0 = -alpha
            alpha = -retrialDist{1}(1,1);
        end
    end
end

% Extract orbit impatience gamma (default 0)
gamma = retInfo.gamma;
if isfield(sn, 'orbitImpatience') && ~isempty(sn.orbitImpatience)
    if size(sn.orbitImpatience, 1) >= queueIdx && size(sn.orbitImpatience, 2) >= classIdx
        impatienceDist = sn.orbitImpatience{queueIdx, classIdx};
        if ~isempty(impatienceDist) && iscell(impatienceDist)
            % For Exp(gamma), D0 = -gamma
            gamma = -impatienceDist{1}(1,1);
        end
    end
end

% Extract batch rejection probability p (default 0)
p = retInfo.p;
if isfield(sn, 'batchRejectProb') && ~isempty(sn.batchRejectProb)
    if size(sn.batchRejectProb, 1) >= queueIdx && size(sn.batchRejectProb, 2) >= classIdx
        p = sn.batchRejectProb(queueIdx, classIdx);
    end
end

% Extract admission threshold R (from FCR or default N-1)
R = retInfo.R;

% Orbit truncation is adaptive, not iter_max; see _kb/06-solver-catalog.md for rationale
maxLevel = [];
if isfield(options, 'config') && isstruct(options.config) ...
        && isfield(options.config, 'orbit_maxlevel') && ~isempty(options.config.orbit_maxlevel)
    maxLevel = options.config.orbit_maxlevel;
end
% setOrbit finite orbit is a MODEL property (real loss), passed straight through;
% see _kb/06-solver-catalog.md for rationale
if isfield(sn, 'orbitMaxJobs') && ~isempty(sn.orbitMaxJobs) ...
        && size(sn.orbitMaxJobs,1) >= queueIdx && size(sn.orbitMaxJobs,2) >= classIdx ...
        && sn.orbitMaxJobs(queueIdx, classIdx) >= 0
    maxLevel = sn.orbitMaxJobs(queueIdx, classIdx);
end

retrialPolicy = RetrialPolicy.LINEAR;
if isfield(sn, 'retrialPolicy') && ~isempty(sn.retrialPolicy) ...
        && size(sn.retrialPolicy,1) >= queueIdx && size(sn.retrialPolicy,2) >= classIdx ...
        && sn.retrialPolicy(queueIdx, classIdx) > 0
    retrialPolicy = sn.retrialPolicy(queueIdx, classIdx);
end

tailTol = 1e-6;
if isfield(options, 'config') && isstruct(options.config) ...
        && isfield(options.config, 'orbit_tailtol') && ~isempty(options.config.orbit_tailtol)
    tailTol = options.config.orbit_tailtol;
end

tol = 1e-10;
if isfield(options, 'tol') && ~isempty(options.tol)
    tol = options.tol;
end

verbose = false;
if isfield(options, 'verbose') && options.verbose
    verbose = true;
end

% Call qsys BMAP/PH/N/N retrial solver
perf = qsys_bmapphnn_retrial(D, beta, S, N, alpha, gamma, p, R, ...
    'MaxLevel', maxLevel, 'Tolerance', tol, 'TailTolerance', tailTol, ...
    'RetrialPolicy', retrialPolicy, 'Verbose', verbose);

% Map to LINE output format
% Queue length includes both orbit and servers
QN(queueIdx, classIdx) = perf.L_orbit + perf.N_server;
UN(queueIdx, classIdx) = perf.Utilization;
TN(queueIdx, classIdx) = perf.Throughput;

% Response time via Little's law
if perf.Throughput > 0
    RN(queueIdx, classIdx) = QN(queueIdx, classIdx) / perf.Throughput;
else
    RN(queueIdx, classIdx) = Inf;
end

% System-level metrics
XN(classIdx) = perf.Throughput;
CN(classIdx) = RN(queueIdx, classIdx);

% Return iteration count (truncation level used)
totiter = perf.truncLevel;

end

%% Helper functions

function assertMarkovian(sn, ist, r, role)
% Reject a station-class process that is not a valid Markovian (D0,D1,...)
% representation. Non-Markovian distributions (Det, traces, NHPP rate
% schedules) and disabled classes reach here as scalars or as NaN-filled
% matrices; they have no place in a matrix-analytic generator.
stationName = sn.nodenames{sn.stationToNode(ist)};
className = sn.classnames{r};
proc = sn.proc{ist}{r};

if isempty(proc) || ~iscell(proc) || numel(proc) < 2
    line_error(mfilename, sprintf(['The %s process of class ''%s'' at station ''%s'' is not a Markovian ' ...
        '(D0,D1) representation. The matrix-analytic retrial solver requires phase-type or ' ...
        'Markovian-arrival distributions; use SolverCTMC, SolverSSA or SolverLDES instead.'], ...
        role, className, stationName));
end

n0 = size(proc{1}, 1);
for e = 1:numel(proc)
    De = proc{e};
    if ~isnumeric(De) || ~ismatrix(De) || size(De,1) ~= size(De,2) || size(De,1) ~= n0
        line_error(mfilename, sprintf(['The %s process of class ''%s'' at station ''%s'' has ' ...
            'inconsistent (D0,D1) block dimensions.'], role, className, stationName));
    end
    if any(~isfinite(De(:)))
        line_error(mfilename, sprintf(['The %s process of class ''%s'' at station ''%s'' contains ' ...
            'NaN or Inf entries: the class is disabled at this station or the distribution has no ' ...
            'Markovian representation.'], role, className, stationName));
    end
end
end

function warnIfApproximated(sn, ist, r)
% Non-Markovian service distributions are replaced by an Erlang MAP of
% matching mean (see MNetwork.refreshProcessRepresentations). That is a
% documented approximation, not the requested distribution, so say so.
if ~isfield(sn, 'procid') || isempty(sn.procid)
    return;
end
if size(sn.procid, 1) < ist || size(sn.procid, 2) < r
    return;
end
switch sn.procid(ist, r)
    case {ProcessType.DET, ProcessType.REPLAYER, ProcessType.UNIFORM, ...
            ProcessType.GAMMA, ProcessType.PARETO, ProcessType.WEIBULL, ...
            ProcessType.LOGNORMAL}
        stationName = sn.nodenames{sn.stationToNode(ist)};
        line_warning(mfilename, sprintf(['Service distribution %s at station ''%s'' is not phase-type. ' ...
            'The matrix-analytic retrial solver uses an Erlang approximation matching its mean and ' ...
            'squared coefficient of variation (%d phases); results are approximate.'], ...
            ProcessType.toText(sn.procid(ist, r)), stationName, sn.phases(ist, r)));
end
end

function D = extractBMAPMatrices(proc)
% Extract BMAP matrices {D0, D1, ...} from process representation
% LINE stores arrival processes in MAP format: {D0, D1, D2, ...}
% where D0 is the "hidden" generator and D1, D2, ... are arrival matrices

D = {};

if isempty(proc) || ~iscell(proc)
    return;
end

% Check if already in BMAP format (cell of cell arrays - unlikely)
if iscell(proc{1})
    D = proc;
    return;
end

% For LINE, arrival processes are stored as {D0, D1, D2, ...}
% D0 is a square matrix with negative diagonal entries (subgenerator)
% D1, D2, ... are non-negative matrices for arrivals

if length(proc) >= 2 && isnumeric(proc{1}) && isnumeric(proc{2})
    D0 = proc{1};
    D1 = proc{2};

    n = size(D0, 1);

    % Check if D0 looks like a generator (negative diagonal)
    % For scalar Exp(rate), D0 = -rate (negative)
    if n == 1
        % Scalar case (exponential)
        if D0 < 0 && D1 > 0
            % Already in MAP format {D0, D1}
            D = proc;
            return;
        end
    else
        % Matrix case - check if D0 has negative diagonal
        diagD0 = diag(D0);
        if all(diagD0 < 0) || all(diagD0 <= 0)
            % Looks like MAP format {D0, D1, ...}
            D = proc;
            return;
        end
    end

    % If we get here, try to interpret as PH format {alpha, T}
    % and convert to MAP
    alpha = D0;  % Initial probability vector
    T = D1;      % Subgenerator

    if isrow(alpha) || (numel(alpha) == size(T, 1) && size(T, 1) == size(T, 2))
        alpha = alpha(:)';
        n = length(alpha);
        if size(T, 1) == n && size(T, 2) == n
            % Convert PH to MAP: D0 = T, D1 = (-T*e)*alpha
            D0_new = T;
            D1_new = (-T * ones(n, 1)) * alpha;
            D = {D0_new, D1_new};
            return;
        end
    end
end

% Fallback: assume it's already in MAP format
D = proc;

end

function [beta, S] = extractPHParams(proc)
% Extract PH parameters (beta, S) from process representation
% LINE stores PH as MAP format: {D0, D1} where D0=T (subgenerator), D1=S0*alpha

beta = [];
S = [];

if isempty(proc) || ~iscell(proc)
    return;
end

% LINE stores PH in MAP format: {D0, D1}
% D0 = T (subgenerator matrix)
% D1 = S0 * alpha (exit rate times initial prob)
if length(proc) >= 2 && isnumeric(proc{1}) && isnumeric(proc{2})
    D0 = proc{1};
    D1 = proc{2};

    % Validate dimensions
    n = size(D0, 1);
    if size(D0, 2) ~= n || size(D1, 1) ~= n || size(D1, 2) ~= n
        return;
    end

    % S = D0 (subgenerator)
    S = D0;

    % S0 = -S * ones (exit rates)
    S0 = -S * ones(n, 1);

    % Extract alpha from D1 = S0 * alpha
    % Find a row with non-zero exit rate
    idx = find(S0 > 1e-10, 1);
    if isempty(idx)
        % All rows have zero exit rate - use uniform
        beta = ones(1, n) / n;
    else
        % alpha = D1(idx,:) / S0(idx)
        beta = D1(idx, :) / S0(idx);
    end

    % Ensure beta is row vector and sums to 1
    beta = beta(:)';
    if abs(sum(beta) - 1) > 1e-6
        % Normalize
        if sum(beta) > 0
            beta = beta / sum(beta);
        else
            beta = ones(1, n) / n;
        end
    end
end

end

%% Reneging (MAPMsG) helper functions

function [QN,UN,RN,TN,CN,XN,totiter] = solveReneging(sn, options, info)
% SOLVERENEGING Solve MAP/M/s+G queue using MAPMsG library
%
% Reference: O. Gursoy, K. A. Mehr, N. Akar, "The MAP/M/s + G Call Center
% Model with Generally Distributed Patience Times"

M = sn.nstations;
K = sn.nclasses;
QN = zeros(M, K);
UN = zeros(M, K);
RN = zeros(M, K);
TN = zeros(M, K);
CN = zeros(1, K);
XN = zeros(1, K);

sourceIdx = info.sourceIdx;
queueIdx = info.queueIdx;
classIdx = info.classIdx;

% Extract MAP arrival process matrices (C, D in MAPMsG notation)
arrivalProc = sn.proc{sourceIdx}{classIdx};
C = arrivalProc{1};  % D0 (subgenerator)
D = arrivalProc{2};  % D1 (arrival transitions)
MAPSIZE = size(C, 1);

% Extract service rate and server count
mu = info.serviceRate;
SERVERSIZE = info.nServers;

% Extract patience distribution and convert to regimes
patienceProc = sn.patienceProc{queueIdx, classIdx};
[BoundaryLevels, ga, QUANTIZATION] = convertPatienceToRegimes(patienceProc, options);

% Build MRMFQ matrices following MAPMsGCompiler.m logic
I = eye(MAPSIZE);
em = ones(MAPSIZE, 1);
lmap = D * em;  % Arrival rate vector

% Initialize Qy0 (boundary generator at level 0)
Qy0 = zeros((SERVERSIZE+1)*MAPSIZE, (SERVERSIZE+1)*MAPSIZE);
for row = 1:SERVERSIZE+1
    if row == 1
        Qy0((row-1)*MAPSIZE+1:(row-1)*MAPSIZE+MAPSIZE, (row-1)*MAPSIZE+1:(row-1)*MAPSIZE+MAPSIZE) = C;
        Qy0((row-1)*MAPSIZE+1:(row-1)*MAPSIZE+MAPSIZE, row*MAPSIZE+1:row*MAPSIZE+MAPSIZE) = D;
    elseif row == SERVERSIZE+1
        Qy0((row-1)*MAPSIZE+1:(row-1)*MAPSIZE+MAPSIZE, (row-1)*MAPSIZE+1:(row-1)*MAPSIZE+MAPSIZE) = -(row-1)*mu*I;
        Qy0((row-1)*MAPSIZE+1:(row-1)*MAPSIZE+MAPSIZE, (row-2)*MAPSIZE+1:(row-2)*MAPSIZE+MAPSIZE) = (row-1)*mu*I;
    else
        Qy0((row-1)*MAPSIZE+1:(row-1)*MAPSIZE+MAPSIZE, (row-1)*MAPSIZE+1:(row-1)*MAPSIZE+MAPSIZE) = C - (row-1)*mu*I;
        Qy0((row-1)*MAPSIZE+1:(row-1)*MAPSIZE+MAPSIZE, (row-2)*MAPSIZE+1:(row-2)*MAPSIZE+MAPSIZE) = (row-1)*mu*I;
        Qy0((row-1)*MAPSIZE+1:(row-1)*MAPSIZE+MAPSIZE, row*MAPSIZE+1:row*MAPSIZE+MAPSIZE) = D;
    end
end

% Initialize Qy for each regime (with abandonment)
Qy = zeros((SERVERSIZE+1)*MAPSIZE, (SERVERSIZE+1)*MAPSIZE, QUANTIZATION);
for regimecount = 1:QUANTIZATION
    for row = SERVERSIZE:SERVERSIZE+1
        if row == SERVERSIZE
            Qy((row-1)*MAPSIZE+1:(row-1)*MAPSIZE+MAPSIZE, (row-1)*MAPSIZE+1:(row-1)*MAPSIZE+MAPSIZE, regimecount) = ga(regimecount+1)*D + C;
            Qy((row-1)*MAPSIZE+1:(row-1)*MAPSIZE+MAPSIZE, row*MAPSIZE+1:row*MAPSIZE+MAPSIZE, regimecount) = (1-ga(regimecount+1))*D;
        elseif row == SERVERSIZE+1
            Qy((row-1)*MAPSIZE+1:(row-1)*MAPSIZE+MAPSIZE, (row-1)*MAPSIZE+1:(row-1)*MAPSIZE+MAPSIZE, regimecount) = -(row-1)*mu*I;
            Qy((row-1)*MAPSIZE+1:(row-1)*MAPSIZE+MAPSIZE, (row-2)*MAPSIZE+1:(row-2)*MAPSIZE+MAPSIZE, regimecount) = (row-1)*mu*I;
        end
    end
end

% Build drift matrices
Rydiag = -ones(1, size(Qy0, 1));
Rydiag(size(Qy0,1)-MAPSIZE+1:size(Qy0,1)) = -Rydiag(size(Qy0,1)-MAPSIZE+1:size(Qy0,1));
Ry = diag(Rydiag);

ydriftregimes = zeros(QUANTIZATION, length(Rydiag));
Ryregimes = zeros(size(Ry,1), size(Ry,2), QUANTIZATION);
for regimecount = 1:QUANTIZATION
    ydriftregimes(regimecount,:) = Rydiag;
    Ryregimes(:,:,regimecount) = Ry;
end

% Combine boundaries and regimes
Qybounds = cat(3, Qy0, Qy);
ydriftbounds = cat(1, Rydiag, ydriftregimes);

% Prepare boundary levels (remove first, add large value at end)
B = BoundaryLevels;
B(1) = [];
B(end+1) = 10000000;

% Call MRMFQ solver
[coefficients, boundaries, Lzeromulti, Lnegmulti, Lposmulti, Anegmulti, Aposmulti] = ...
    MRMFQSolver(Qy, Qybounds, ydriftregimes, ydriftbounds, B);

% Compute steady-state results following MAPMsGCompiler.m
zeromass = boundaries{1};

% Compute integrals for each regime
integral = zeros(QUANTIZATION, length(zeromass));
waitintegral = zeros(QUANTIZATION, length(zeromass));
abandonintegral = zeros(QUANTIZATION, length(zeromass));
successfulintegral = zeros(QUANTIZATION, length(zeromass));

for d = 1:QUANTIZATION
    if d == 1
        prevB = 0;
    else
        prevB = B(d-1);
    end

    Lz = Lzeromulti{d};
    Ln = Lnegmulti{d};
    Lp = Lposmulti{d};
    An = Anegmulti{d};
    Ap = Aposmulti{d};
    coef = coefficients{d};

    deltaB = B(d) - prevB;

    integrand = [Lz * deltaB; ...
                 An \ (expm(An * deltaB) - eye(size(An))) * Ln; ...
                 Ap \ (eye(size(Ap)) - expm(-Ap * deltaB)) * Lp];

    integral(d,:) = coef * integrand;
    waitintegral(d,:) = ((B(d) + prevB)/2) * (1 - ga(d+1)) * coef * integrand;
    abandonintegral(d,:) = ga(d+1) * coef * integrand;
    successfulintegral(d,:) = (1 - ga(d+1)) * coef * integrand;
end

% Map integrals to arrival rates
normalization = SERVERSIZE * MAPSIZE;
IntegralMapped = zeros(size(integral));
AbandonIntegralMapped = zeros(size(integral));
WaitIntegralMapped = zeros(size(integral));
ZeroMassMapped = zeros(1, length(zeromass));

for r = 1:length(zeromass)
    lmapIdx = mod(r-1, MAPSIZE) + 1;
    IntegralMapped(:,r) = integral(:,r) * lmap(lmapIdx);
    AbandonIntegralMapped(:,r) = abandonintegral(:,r) * lmap(lmapIdx);
    WaitIntegralMapped(:,r) = waitintegral(:,r) * lmap(lmapIdx);
    ZeroMassMapped(r) = zeromass(r) * lmap(lmapIdx);
end

% Compute performance metrics
totalMass = sum(ZeroMassMapped) + sum(sum(IntegralMapped(:,1:normalization)));
AbandonProb = sum(sum(AbandonIntegralMapped(:,1:normalization))) / totalMass;
ExpectedWait = sum(sum(WaitIntegralMapped(:,1:normalization))) / (totalMass * (1 - AbandonProb));

% Compute arrival rate (lambda)
lambda = sum(lmap);

% Map to LINE output format
% Throughput = arrival rate * (1 - abandonment probability)
throughput = lambda * (1 - AbandonProb);
TN(queueIdx, classIdx) = throughput;

% Utilization
UN(queueIdx, classIdx) = throughput / (mu * SERVERSIZE);

% Response time = expected wait + expected service time
RN(queueIdx, classIdx) = ExpectedWait + 1/mu;

% Queue length via Little's law
QN(queueIdx, classIdx) = throughput * RN(queueIdx, classIdx);

% System-level metrics
XN(classIdx) = throughput;
CN(classIdx) = RN(queueIdx, classIdx);

% Iteration count
totiter = QUANTIZATION;

end

function [BoundaryLevels, ga, QUANTIZATION] = convertPatienceToRegimes(patienceProc, options)
% CONVERTPATIENCETOREGIMES Convert patience distribution to MAPMsG regimes
%
% Converts a patience distribution (in MAP/PH format) to piecewise-constant
% abandonment function for MAPMsG.

% Default quantization
QUANTIZATION = 11;
if isfield(options, 'config') && isfield(options.config, 'mapmsg_quantization')
    QUANTIZATION = options.config.mapmsg_quantization;
end

% Extract patience rate from the distribution
% For exponential patience Exp(gamma): patienceProc = {-gamma, gamma}
if iscell(patienceProc) && length(patienceProc) >= 2
    % Extract rate from subgenerator
    gamma = -patienceProc{1}(1,1);
else
    % Default patience rate
    gamma = 0.1;
end

% Generate boundary levels and abandonment probabilities
% For exponential patience with rate gamma:
% F(t) = 1 - exp(-gamma*t) is the CDF (abandonment probability by time t)
maxTime = 10;  % Maximum time horizon for regimes
BoundaryLevels = linspace(0, maxTime, QUANTIZATION);

% Compute abandonment probability at each boundary
% ga(k) = F(BoundaryLevels(k)) for piecewise constant approximation
ga = zeros(1, QUANTIZATION + 1);
ga(1) = 0;  % No abandonment at time 0
for k = 1:QUANTIZATION
    % Abandonment probability at midpoint of regime
    midpoint = (BoundaryLevels(k) + BoundaryLevels(min(k+1, QUANTIZATION))) / 2;
    if k == 1
        midpoint = BoundaryLevels(k) / 2;
    end
    ga(k+1) = 1 - exp(-gamma * midpoint);
end

end
