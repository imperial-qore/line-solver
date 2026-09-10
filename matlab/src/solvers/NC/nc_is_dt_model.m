function dt = nc_is_dt_model(sn, options)
% DT = NC_IS_DT_MODEL(SN, OPTIONS)
%
% Classify SN against the two families of discrete-time product-form models
% of Daduna (2001) that SOLVER_NC_DT_ANALYZER can solve exactly. The route is
% requested with OPTIONS.CONFIG.SLOTTED, the same knob SolverLDES uses to run
% on a discrete time scale, and is never auto-detected: a Geometric service
% time is a perfectly ordinary continuous-time model unless the caller says
% the model lives on a slot lattice.
%
% Returns a struct with field KIND:
%   'bernoulli1' - one Bernoulli server fed by a state dependent Bernoulli
%                  arrival stream (chapter 2). Fields STATION, ARRIVALPROB,
%                  SERVICEPROB, CAPACITY.
%   'cycle'      - closed cycle of Bernoulli servers (chapter 3). Fields
%                  ORDER (stations in cycle order), SERVICEPROB (J-by-N
%                  matrix p_j(n)), POPULATION.
%   'none'       - not a discrete-time product-form model. Field REASON says
%                  why, and the analyzer turns it into an error rather than
%                  falling back to a continuous-time approximation.
%
% The admissible feature set is narrow because the discrete-time product form
% is narrow. Beyond the geometric service requirement:
%   - a cycle is the only topology; Daduna, section 4.1, records that general
%     discrete-time topologies of FCFS Bernoulli servers have no product form;
%   - every station must be a single server. Pestien and Ramakrishnan, quoted
%     in the same reference before example 2.10, proved that a multiserver
%     node inside a cycle of geometrical queues destroys the product form for
%     any finite server count;
%   - class switching is rejected, and a multichain cycle is admitted only
%     through the aggregate population (see SOLVER_NC_DT_ANALYZER).
%
% See also SOLVER_NC_DT_ANALYZER, DPFQN_NC, DPFQN_NCLD, DQSYS_BERNOULLI1
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

dt = struct('kind', 'none', 'reason', '');

M = sn.nstations;
R = sn.nclasses;

% ---- geometric service everywhere a class is actually served --------------
served = false(M, R);
for ist = 1:M
    for r = 1:R
        if isfinite(sn.rates(ist, r)) && sn.rates(ist, r) > 0
            served(ist, r) = true;
        end
    end
end
for ist = 1:M
    for r = 1:R
        if served(ist, r) && sn.procid(ist, r) ~= ProcessType.GEOMETRIC
            dt.reason = sprintf('station %d serves class %d with a %s process; a discrete-time model needs Geometric service and interarrival times', ...
                ist, r, ProcessType.toText(sn.procid(ist, r)));
            return
        end
    end
end

% ---- class switching is out of scope --------------------------------------
if isfield(sn, 'csmask') && ~isempty(sn.csmask)
    cs = sn.csmask;
    cs(1:R+1:end) = false;
    if any(cs(:))
        dt.reason = 'class switching is not covered by the discrete-time product form';
        return
    end
end
if ~isempty(sn.cdscaling) && any(~cellfun(@isempty, sn.cdscaling))
    dt.reason = 'class-dependent scaling is not covered by the discrete-time product form';
    return
end
if ~isempty(sn.jdscaling) && any(~cellfun(@isempty, sn.jdscaling))
    dt.reason = 'joint-dependent scaling is not covered by the discrete-time product form';
    return
end

isOpen = any(isinf(sn.njobs));
if isOpen
    dt = dt_open_single(sn, dt);
else
    dt = dt_closed_cycle(sn, dt);
end
end

% ==========================================================================
function dt = dt_open_single(sn, dt)
% Chapter 2: Source -> Bernoulli server -> Sink, one open class.
M = sn.nstations;
R = sn.nclasses;

if R ~= 1
    dt.reason = 'the discrete-time single-node route handles one open class';
    return
end
srcs = find(sn.sched == SchedStrategy.EXT);
if numel(srcs) ~= 1
    dt.reason = 'an open discrete-time model needs exactly one Source';
    return
end
queues = setdiff(1:M, srcs);
if numel(queues) ~= 1
    dt.reason = sprintf('the discrete-time single-node route handles one queueing station, found %d', numel(queues));
    return
end
ist = queues;
if sn.sched(ist) ~= SchedStrategy.FCFS
    dt.reason = 'a Bernoulli server is a FCFS station';
    return
end
if isfinite(sn.nservers(ist)) && sn.nservers(ist) ~= 1
    dt.reason = 'a Bernoulli server is a single-server station; use load dependence for the multiserver approximation of example 2.10';
    return
end

b = sn.rates(srcs, 1);                 % Geometric(b) interarrivals
p = sn.rates(ist, 1);                  % Geometric(p) service
if ~isfinite(b) || b <= 0 || b > 1
    dt.reason = 'the source arrival probability must lie in (0,1]';
    return
end
if ~isfinite(p) || p <= 0 || p > 1
    dt.reason = 'the service probability must lie in (0,1]';
    return
end

cap = dt_capacity(sn, ist, 1);
lld = dt_lld(sn, ist);
if isinf(cap) && ~isempty(lld)
    dt.reason = 'a load-dependent Bernoulli server needs a finite capacity to bound the state space';
    return
end
if isinf(cap) && b >= p
    dt.reason = 'an unbounded discrete-time queue needs an arrival probability below the service probability';
    return
end

if isinf(cap)
    svc = p;
else
    svc = p * dt_lld_vector(lld, cap);
    if any(svc <= 0) || any(svc > 1)
        dt.reason = 'load dependence must keep the service probability inside (0,1]';
        return
    end
end

dt.kind = 'bernoulli1';
dt.station = ist;
dt.source = srcs;
dt.arrivalProb = b;
dt.serviceProb = svc;
dt.capacity = cap;
end

% ==========================================================================
function dt = dt_closed_cycle(sn, dt)
% Chapter 3: a closed cycle of state dependent Bernoulli servers.
M = sn.nstations;
R = sn.nclasses;
N = sum(sn.njobs(isfinite(sn.njobs)));
if N <= 0 || floor(N) ~= N
    dt.reason = 'the closed population must be a positive integer';
    return
end

for ist = 1:M
    if sn.sched(ist) ~= SchedStrategy.FCFS
        dt.reason = sprintf('station %d is not FCFS; a cycle of Bernoulli servers is FCFS throughout', ist);
        return
    end
    if isfinite(sn.nservers(ist)) && sn.nservers(ist) ~= 1
        dt.reason = sprintf('station %d has %d servers; a multiserver node inside a cycle of geometrical queues has no product form', ist, sn.nservers(ist));
        return
    end
end

% One service probability per station: the discrete-time cycle serves every
% class at the node rate, so a class-dependent rate is outside the model.
p = zeros(1, M);
for ist = 1:M
    rr = sn.rates(ist, :);
    rr = rr(isfinite(rr) & rr > 0);
    if isempty(rr)
        dt.reason = sprintf('station %d serves no class', ist);
        return
    end
    if (max(rr) - min(rr)) > GlobalConstants.FineTol * max(rr)
        dt.reason = sprintf('station %d has a class-dependent service probability; the discrete-time cycle needs one Bernoulli server per node', ist);
        return
    end
    p(ist) = rr(1);
    if p(ist) <= 0 || p(ist) >= 1
        dt.reason = sprintf('station %d has service probability %g; the product form of theorem 3.2 needs p in (0,1)', ist, p(ist));
        return
    end
end

order = dt_cycle_order(sn);
if isempty(order)
    dt.reason = 'the stations do not form a single deterministic cycle; discrete-time FCFS networks of other topologies have no product form';
    return
end

P = zeros(M, max(N,1));
for k = 1:M
    ist = order(k);
    alpha = dt_lld_vector(dt_lld(sn, ist), N);
    P(k, :) = p(ist) * alpha;
end
if any(P(:) <= 0) || any(P(:) > 1)
    dt.reason = 'load dependence must keep every service probability inside (0,1]';
    return
end

dt.kind = 'cycle';
dt.order = order;
dt.serviceProb = P;
dt.population = N;
end

% ==========================================================================
function order = dt_cycle_order(sn)
% Station order along the cycle, or empty when the routing is not a single
% deterministic cycle visiting every station exactly once.
order = [];
M = sn.nstations;
rtst = sn_rt_stations(sn);
R = sn.nclasses;

succ = zeros(1, M);
for i = 1:M
    tgt = 0;
    for j = 1:M
        w = 0;
        for r = 1:R
            for s = 1:R
                w = w + rtst((i-1)*R+r, (j-1)*R+s);
            end
        end
        if w > GlobalConstants.FineTol
            if abs(w - R) > GlobalConstants.CoarseTol && abs(w - 1) > GlobalConstants.CoarseTol
                return   % fractional routing out of station i
            end
            if tgt ~= 0
                return   % more than one successor
            end
            tgt = j;
        end
    end
    if tgt == 0 || tgt == i
        return
    end
    succ(i) = tgt;
end

visited = false(1, M);
order = zeros(1, M);
cur = 1;
for k = 1:M
    if visited(cur)
        order = [];
        return
    end
    visited(cur) = true;
    order(k) = cur;
    cur = succ(cur);
end
if cur ~= 1 || ~all(visited)
    order = [];
end
end

% ==========================================================================
function cap = dt_capacity(sn, ist, r)
% Effective buffer capacity of station ist for class r, Inf when unbounded.
cap = Inf;
if isfield(sn, 'cap') && ~isempty(sn.cap) && numel(sn.cap) >= ist
    cap = min(cap, sn.cap(ist));
end
if isfield(sn, 'classcap') && ~isempty(sn.classcap) && size(sn.classcap,1) >= ist
    cap = min(cap, sn.classcap(ist, r));
end
end

% ==========================================================================
function lld = dt_lld(sn, ist)
% Load-dependent scaling vector of station ist, empty when undeclared.
lld = [];
if isfield(sn, 'lldscaling') && ~isempty(sn.lldscaling) && size(sn.lldscaling,1) >= ist
    row = sn.lldscaling(ist, :);
    if ~isempty(row) && any(abs(row - 1) > GlobalConstants.FineTol)
        lld = row;
    end
end
end

% ==========================================================================
function alpha = dt_lld_vector(lld, N)
% Expand a load-dependent scaling to alpha(1..N), holding the last declared
% value beyond the tabulated range, as the load-dependent NC solvers do.
alpha = ones(1, N);
if isempty(lld)
    return
end
lld = lld(:)';
n = min(N, numel(lld));
alpha(1:n) = lld(1:n);
if N > numel(lld)
    alpha(numel(lld)+1:N) = lld(end);
end
end
