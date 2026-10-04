%{ @file cache_retrieval_inputs.m
 %  @brief Extract delayed-hit retrieval algorithm inputs from a NetworkStruct
 %
 %  @author LINE Development Team
%}

%{
 % @brief Builds the inputs of the retrieval (delayed-hit) analytic algorithms from a model
 %
 % @details
 % Given the NetworkStruct of a cache equipped with a retrieval system
 % (Cache.setRetrievalSystem), reconstructs the inputs required by the
 % retrieval_* algorithms (retrieval_nc, retrieval_metrics, retrieval_fpi,
 % retrieval_fpi_latency):
 %
 %   m            cache list capacities (1 x h)
 %   lambda       per-item arrival rates (1 x n)         = sourceRate * pread
 %   gamma        access factors gamma_{i,j} (n x h)     via cache_gamma_lp
 %   eta          fetching demands eta_{s,i} (n x (r+1)) col 1 = IS, cols 2.. = PS
 %   alpha        cell(1,S); alpha{s}(1,:,i) PH entry vector of item i at station s
 %   T            cell(1,S); T{s}(:,:,i) PH subgenerator of item i at station s
 %   R            (S+1) x (S+1) x n routing matrices (index 1 = cache/outside,
 %                2..S+1 = retrieval stations), R(a,b,i) probability a->b for item i
 %   station_type (1 x S) string per retrieval station, one of "IS", "PS",
 %                "SIRO", "FCFS", "LCFSPR" (mapped from sn.sched)
 %
 % The retrieval system is single-class (IRM): exactly one read class may route
 % into the retrieval system. Supported retrieval scheduling policies are IS,
 % PS, SIRO, FCFS and LCFSPR. IS is independent; PS, SIRO, FCFS and LCFSPR use
 % the mean-field sharing slowdown. PS and LCFSPR (symmetric/insensitive BCMP
 % disciplines) admit general phase-type service and class-dependent rates;
 % SIRO and FCFS require exponential service with identical per-class rates.
 % Any other scheduling policy raises an error, matching retrieval_fpi_latency.
 %
 % @par Syntax:
 % @code
 % [m,lambda,gamma,eta,alpha,T,R,station_type] = cache_retrieval_inputs(sn)
 % @endcode
%}
function [m,lambda,gamma,eta,alpha,T,R,station_type] = cache_retrieval_inputs(sn, lambdaOverride)
% [...] = CACHE_RETRIEVAL_INPUTS(SN) builds the retrieval-algorithm inputs from
% an OPEN cache model (per-item rate lambda = sourceRate * pread).
%
% [...] = CACHE_RETRIEVAL_INPUTS(SN, LAMBDAOVERRIDE) uses the supplied read-class
% arrival rate LAMBDAOVERRIDE (a scalar rate for the single read class) instead of
% the Source throughput, so the same inputs can be built for a CLOSED integrated
% cache-queueing sublayer where the read rate comes from the network solution
% (da_cacheqn_retrieval). All other inputs (gamma, eta, alpha, T, R) are unchanged.

if nargin < 2
    lambdaOverride = [];
end

ci = find(sn.nodetype == NodeType.Cache);
if numel(ci) ~= 1
    line_error(mfilename, 'Retrieval analysis requires exactly one Cache node.');
end
ch = sn.nodeparam{ci};
if ~isfield(ch, 'retrievalSystemCapacity') || ch.retrievalSystemCapacity <= 0
    line_error(mfilename, 'The Cache node has no retrieval system (call setRetrievalSystem).');
end

m = ch.itemcap(:).';
n = ch.nitems;
h = numel(m);
K = sn.nclasses;

% --- read class (single-class IRM) ---
rk = keys(ch.retrievalSystemQueueIndices);
if numel(rk) ~= 1
    line_error(mfilename, 'Retrieval analysis supports a single read class.');
end
jobinClass = double(rk(1)) + 1;               % stored 0-indexed
queueNodes = ch.retrievalSystemQueueIndices{rk(1)};
queueNodes = double(queueNodes(:).');
S = numel(queueNodes);
if S == 0
    line_error(mfilename, 'The retrieval system has no stations.');
end

% --- per-item arrival rates lambda(i) = readRate * pread(i) ---
% readRate is the Source throughput for an open model, or the caller-supplied
% closed read-class arrival rate (from the network solution) when overridden.
pread = ch.pread{jobinClass};
if ~isempty(lambdaOverride)
    readRate = lambdaOverride;
else
    source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
    if isempty(source_ist)
        line_error(mfilename, ['Retrieval analysis of a closed model requires an explicit read ' ...
            'rate (call cache_retrieval_inputs(sn, lambdaOverride)); no Source node found.']);
    end
    readRate = sn.rates(source_ist, jobinClass);
end
if isnan(readRate), readRate = 0; end
lambda = readRate * pread(:).';             % 1 x n

% --- gamma via the existing plain-cache utility (n x h) ---
lambda3d = zeros(1, n, h);
for k = 1:n
    for l = 1:(h+1)
        lambda3d(1, k, l) = lambda(k);
    end
end
Rcost = ch.accost;
if isempty(Rcost)
    % Default linear cache routing: item flows from list l to list l+1.
    Rcost = cell(1, n);
    for k = 1:n
        Rmat = diag(ones(1, h), 1);
        Rmat(h+1, h+1) = 1;
        Rcost{1, k} = Rmat;
    end
end
gamma = cache_gamma_lp(lambda3d, Rcost);      % n x h

% --- station types (IS / PS / SIRO / FCFS) ---
% see _kb/09-ldes-and-cache.md (delayed-hit latency: station-type equivalences)
station_type = strings(1, S);
for s = 1:S
    sst = sn.nodeToStation(queueNodes(s));
    if sn.sched(sst) == SchedStrategy.INF
        station_type(s) = "IS";
    elseif sn.sched(sst) == SchedStrategy.PS
        station_type(s) = "PS";
    elseif sn.sched(sst) == SchedStrategy.SIRO
        station_type(s) = "SIRO";
    elseif sn.sched(sst) == SchedStrategy.FCFS
        station_type(s) = "FCFS";
    elseif sn.sched(sst) == SchedStrategy.LCFSPR
        station_type(s) = "LCFSPR";
    else
        line_error(mfilename, ['Retrieval analysis supports only IS, PS, SIRO, FCFS and LCFSPR ' ...
            'retrieval stations; station %d uses an unsupported scheduling policy.'], queueNodes(s));
    end
end

% --- per-item PH service (alpha,T) per station, routing R ---
alpha = cell(1, S);
T = cell(1, S);
fsz = zeros(1, S);
for s = 1:S
    sst = sn.nodeToStation(queueNodes(s));
    rcls0 = ch.retrievalClasses(1, jobinClass);
    fsz(s) = size(sn.proc{sst}{rcls0}{1}, 1);
    if (station_type(s) == "SIRO" || station_type(s) == "FCFS") && fsz(s) > 1
        line_error(mfilename, ['Retrieval analysis supports SIRO/FCFS retrieval stations only with ' ...
            'exponential (single-phase) service; station %d has phase-type service.'], queueNodes(s));
    end
    alpha{s} = zeros(1, fsz(s), n);
    T{s} = zeros(fsz(s), fsz(s), n);
end
R = zeros(S+1, S+1, n);
lin = @(node, cls) (node-1)*K + cls;          % rtnodes flat index
for i = 1:n
    rcls = ch.retrievalClasses(i, jobinClass);
    for s = 1:S
        sst = sn.nodeToStation(queueNodes(s));
        alpha{s}(1, :, i) = sn.pie{sst}{rcls}(:).';
        T{s}(:, :, i) = sn.proc{sst}{rcls}{1};
    end
    % routing: index 1 = cache (outside), 2..S+1 = retrieval stations
    for s = 1:S
        R(1, s+1, i) = sn.rtnodes(lin(ci, rcls), lin(queueNodes(s), rcls));        % cache -> queue s
        R(s+1, 1, i) = sn.rtnodes(lin(queueNodes(s), rcls), lin(ci, rcls));        % queue s -> cache
        for sp = 1:S
            R(s+1, sp+1, i) = sn.rtnodes(lin(queueNodes(s), rcls), lin(queueNodes(sp), rcls));
        end
    end
end

% --- eta(i,1) = sum of IS visits*mean; eta(i,1+p) = PS station p ---
% SIRO and FCFS reduce to PS only with class-independent service rates; LCFSPR
% is insensitive and exempt.
for s = find(station_type == "FCFS" | station_type == "SIRO")
    taus = zeros(1, n);
    for i = 1:n
        taus(i) = -alpha{s}(:, :, i) / T{s}(:, :, i) * ones(fsz(s), 1);
    end
    if max(taus) - min(taus) > 1e-9 * max(taus)
        line_error(mfilename, ['Retrieval analysis requires class-independent (identical) mean ' ...
            'service rates at SIRO/FCFS station %d.'], queueNodes(s));
    end
end

isIdx = find(station_type == "IS");
psIdx = find(station_type == "PS" | station_type == "SIRO" | station_type == "FCFS" | station_type == "LCFSPR");   % SIRO/FCFS/LCFSPR as PS
r = numel(psIdx);
eta = zeros(n, r+1);
for i = 1:n
    Ri = R(:, :, i);
    a = Ri(1, 2:S+1);                         % outside -> station entry probs
    Pmat = Ri(2:S+1, 2:S+1);                  % station -> station
    visits = a / (eye(S) - Pmat);             % expected visits per fetch
    tau = zeros(1, S);
    for s = 1:S
        tau(s) = -alpha{s}(:, :, i) / T{s}(:, :, i) * ones(fsz(s), 1);
    end
    eta_s = visits .* tau;
    eta(i, 1) = sum(eta_s(isIdx));
    for p = 1:r
        eta(i, 1+p) = eta_s(psIdx(p));
    end
end
end
