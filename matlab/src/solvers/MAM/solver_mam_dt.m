function [QN,UN,RN,TN,CN,XN,totiter,method,pql] = solver_mam_dt(sn, options, slotLength)
% [QN,UN,RN,TN,CN,XN,TOTITER,METHOD,PQL] = SOLVER_MAM_DT(SN, OPTIONS, SLOTLENGTH)
%
% Discrete-time (slotted) analysis of an open network whose interarrival and
% service laws all live on the slot lattice. A single queueing station is
% solved EXACTLY by the Q-MAM discrete-time algorithms, Q_DT_PH_PH_1 when both
% laws are renewal discrete phase-type and Q_DT_MAP_MAP_1 when either side is a
% DMAP. Several stations are solved by a discrete-time parametric
% decomposition, which is an approximation: see solver_mam_dt_network below.
%
% Time is measured in slots internally and converted back on exit, so QN and
% UN are dimensionless, TN is per time unit and RN is in time units.
%
% Convention: late arrival system with delayed access (LAS-DA). Within a slot
% a completion resolves before the arrival, an arrival cannot enter service in
% the slot it arrives in, and every metric is read after both events. This is
% the convention of the Q-MAM discrete-time queues and of the LDES slotted
% engine, so the two are directly comparable.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
K = sn.nclasses;
QN = zeros(M,K); UN = zeros(M,K); RN = zeros(M,K); TN = zeros(M,K);
CN = zeros(1,K); XN = zeros(1,K);
totiter = 1;
pql = [];

% The scope rules live in MAM_DT_SUPPORTS, which SolverMAM.supportsModelMethod
% asks for every method on a slotted model: one body, two callers.
[dtOk, dtWhy] = mam_dt_supports(sn);
if ~dtOk
    line_error(mfilename, dtWhy);
end

% Per-station discrete-time processes, in slots
DMAPst = cell(M,1);
for ist = 1:M
    for r = 1:K
        if sn.procid(ist,r) == ProcessType.DISABLED || isnan(sn.rates(ist,r))
            continue;
        end
        DMAPst{ist} = solver_mam_dt_law(sn, ist, r, slotLength);
    end
end

sourceIdx = find(sn.sched == SchedStrategy.EXT);
queueIdx = find(sn.sched ~= SchedStrategy.EXT);

if numel(queueIdx) == 1
    [QNq, UNq, TNq, pql] = solver_mam_dt_qsys(DMAPst{sourceIdx}, DMAPst{queueIdx}, sn, sourceIdx, queueIdx, options);
    method = 'dt.qmam';
else
    [QNs, UNs, TNs, totiter] = solver_mam_dt_network(sn, DMAPst, sourceIdx, queueIdx, options);
    QNq = QNs; UNq = UNs; TNq = TNs;
    method = 'dt.dec';
end

% Station metrics, back in time units
lambdaSlot = dmap_lambda(DMAPst{sourceIdx});
for idx = 1:numel(queueIdx)
    ist = queueIdx(idx);
    QN(ist,1) = QNq(idx);
    UN(ist,1) = UNq(idx);
    TN(ist,1) = TNq(idx) / slotLength;
    if TNq(idx) > 0
        RN(ist,1) = QNq(idx) / TNq(idx) * slotLength;
    end
end
TN(sourceIdx,1) = lambdaSlot / slotLength;

XN(1) = lambdaSlot / slotLength;
CN(1) = sum(QN(:,1)) / XN(1);

end

%% Discrete-time law of station ist, class r, expressed in slots
function DMAP = solver_mam_dt_law(sn, ist, r, slotLength)
procType = sn.procid(ist,r);
meanSlots = 1 / (sn.rates(ist,r) * slotLength);
if procType == ProcessType.DMAP
    % the (D0,D1) pair survives verbatim in sn.proc: it already has MAP shape
    DMAP = sn.proc{ist}{r};
    if slotLength ~= 1
        line_error(mfilename, ['A DMAP is defined on its own slot, so it cannot be combined with ' ...
            'options.config.slotlength=%g. Set the slot length to 1 or rescale the DMAP.'], slotLength);
    end
else
    [alpha, A] = dph_from_dist(procType, meanSlots, sn.scv(ist,r));
    DMAP = dph_to_dmap(alpha, A);
end
end

%% Exact single-station analysis through the Q-MAM discrete-time queues
function [QN, UN, TN, pql] = solver_mam_dt_qsys(ARV, SVC, sn, sourceIdx, queueIdx, options)
maxNumComp = 1000;
if isfield(options,'config') && isfield(options.config,'dt_maxlevel') && ~isempty(options.config.dt_maxlevel)
    maxNumComp = options.config.dt_maxlevel;
end

arvIsRenewal = dmap_is_renewal(ARV);
svcIsRenewal = dmap_is_renewal(SVC);

if arvIsRenewal && svcIsRenewal
    % both sides are renewal discrete phase-type: the PH/PH/1 entry point
    [alpha, T] = dmap_to_dph(ARV);
    [beta, S] = dmap_to_dph(SVC);
    ql = Q_DT_PH_PH_1(alpha, T, beta, S, 'MaxNumComp', maxNumComp);
else
    ql = Q_DT_MAP_MAP_1(ARV{1}, ARV{2}, SVC{1}, SVC{2}, 'MaxNumComp', maxNumComp);
end

ql = ql / sum(ql);
pql = ql;
QN = (0:(length(ql)-1)) * ql(:);
UN = 1 - ql(1);
TN = dmap_lambda(ARV);
end

%% Discrete-time parametric decomposition over several stations
function [QN, UN, TN, totiter] = solver_mam_dt_network(sn, DMAPst, sourceIdx, queueIdx, options)
% Each station is solved as a DBMAP/DMAP/1 queue given its arrival stream, and
% its departure stream is extracted from the truncated stationary chain and
% split by the routing probabilities. Superposing discrete streams produces
% batches, which is why the station solve is M/G/1-type rather than a QBD.
% Unlike the single-station entry this is an APPROXIMATION: the departure
% process is truncated at a finite level and compressed back to a bounded
% phase dimension, and correlations between the streams entering a station are
% not preserved.

nq = numel(queueIdx);
P = solver_mam_dt_routing(sn, sourceIdx, queueIdx);
V = cellsum(sn.visits);

lambda = dmap_lambda(DMAPst{sourceIdx});
spaceMax = 128;
if isfield(options,'config') && isfield(options.config,'space_max') && ~isempty(options.config.space_max)
    spaceMax = options.config.space_max;
end
iterMax = 100;
if isfield(options,'iter_max') && ~isempty(options.iter_max)
    iterMax = options.iter_max;
end
iterTol = GlobalConstants.CoarseTol;
if isfield(options,'iter_tol') && ~isempty(options.iter_tol)
    iterTol = options.iter_tol;
end

% Initial departure streams: Bernoulli of the exact station throughput, the
% slotted counterpart of seeding a decomposition with Poisson streams
DEP = cell(nq,1);
for idx = 1:nq
    p = lambda * V(queueIdx(idx),1);
    DEP{idx} = {1-p, p};
end
SRC = DMAPst{sourceIdx};

QN = zeros(nq,1); UN = zeros(nq,1); TN = zeros(nq,1);
QNprev = QN; QNprev2 = QN;
UNprev = UN; TNprev = TN;
totiter = 0;
for it = 1:iterMax
    totiter = it;
    for idx = 1:nq
        ARV = solver_mam_dt_arrivals(SRC, DEP, P, idx, spaceMax);
        [QN(idx), UN(idx), TN(idx), ~, DEPidx] = mg1_dt_queue(ARV, DMAPst{queueIdx(idx)}, options);
        DEP{idx} = dmap_compress(DEPidx, spaceMax);
    end
    if it > 1 && max(abs(QN - QNprev) ./ max(QNprev, GlobalConstants.Zero)) < iterTol
        break;
    end
    if it > 2 && max(abs(QN - QNprev2) ./ max(QNprev2, GlobalConstants.Zero)) < iterTol
        % Feedback loops settle into a period-two cycle rather than a point:
        % re-solving a station with the departure process it just produced
        % moves it back. The cycle amplitude sits far below the error of the
        % decomposition itself, so the midpoint is reported instead of
        % burning iter_max sweeps on an orbit that will not close.
        QN = (QN + QNprev) / 2;
        UN = (UN + UNprev) / 2;
        TN = (TN + TNprev) / 2;
        break;
    end
    QNprev2 = QNprev;
    QNprev = QN; UNprev = UN; TNprev = TN;
end
end

%% Arrival stream of queue idx: source share plus thinned upstream departures
function ARV = solver_mam_dt_arrivals(SRC, DEP, P, idx, spaceMax)
ARV = [];
if P(1, idx) > 0
    ARV = dmap_thin(SRC, P(1, idx));
end
for j = 1:numel(DEP)
    p = P(1 + j, idx);
    if p <= 0
        continue;
    end
    contrib = dmap_thin(DEP{j}, p);
    if isempty(ARV)
        ARV = contrib;
    else
        ARV = dmap_super(ARV, contrib);
        ARV = dmap_compress_batch(ARV, spaceMax);
    end
end
if isempty(ARV)
    line_error(mfilename, 'Queue %d receives no arrivals in the discrete-time routing matrix.', idx);
end
end

%% Station-to-station routing probabilities, source first, sink stripped
function P = solver_mam_dt_routing(sn, sourceIdx, queueIdx)
I = sn.nnodes;
K = sn.nclasses;
Pn = zeros(I, I);
for ind = 1:I
    for jnd = 1:I
        Pn(ind, jnd) = sn.rtnodes((ind-1)*K + 1, (jnd-1)*K + 1);
    end
end
% the Sink feeds back into the Source to keep sn.rt stochastic; an open
% traffic equation must not see that edge
sinkNodes = find(sn.nodetype == NodeType.Sink);
Pn(sinkNodes, :) = 0;

stationNodes = [sn.stationToNode(sourceIdx), arrayfun(@(x) sn.stationToNode(x), queueIdx(:)')];
interNodes = setdiff(1:I, [stationNodes, sinkNodes(:)']);

Pss = Pn(stationNodes, stationNodes);
if isempty(interNodes)
    P = Pss;
else
    Psn = Pn(stationNodes, interNodes);
    Pnn = Pn(interNodes, interNodes);
    Pns = Pn(interNodes, stationNodes);
    P = Pss + Psn * ((eye(size(Pnn)) - Pnn) \ Pns);
end
% column 1 is the source, which receives nothing
P = P(:, 2:end);
end
