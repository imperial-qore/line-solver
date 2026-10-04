function terms = fluid_moment_terms(sn, options)
% TERMS = FLUID_MOMENT_TERMS(SN, OPTIONS)
%
% Event-based representation of the fluid population process, as required by
% the moment-closure methods of SolverFLD.
%
% The closing ODEs are a density-dependent Markov population process
%
%   dx/dt = F(x) = D * r(x),     r_e(x) = rateBase(e) * g_e(x)
%
% with D the jump matrix of ODE_JUMPS_NEW and g the rate-factor vector of
% ODE_RATES_CLOSING. SOLVER_FLUID_ODES discards D and r once it has composed
% the right-hand side, but the covariance equation of the linear noise
% approximation and the 1/N refinement need them separately: the diffusion
% matrix is D*diag(r(x))*D', which cannot be recovered from F alone. This
% function rebuilds that representation from SN and returns it together with
% drift, rate and Jacobian handles that all take the closure variance as an
% explicit second argument.
%
% Parameters:
%   sn      - NetworkStruct
%   options - solver options
%
% Returns:
%   terms - struct with fields D, rateBase, eventIdx, q_indices, Kic,
%           enabled, w, sched, S, M, K, nstate, stationBlock, classBlock,
%           minExact, ratesFcn, driftFcn, jacFcn
%
% See also SOLVER_FLUID_MOMENTS, FLUID_DRIFT_JACOBIAN, ODE_JUMPS_NEW.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
K = sn.nclasses;
N = sn.nclosedjobs;
Mu = sn.mu;
Phi = sn.phi;
PH = sn.proc;
sched = sn.sched;
% station-major routing: sn.rt is indexed by stateful node, and the jump set
% below indexes it by station (see SN_RT_STATIONS)
rt = sn_rt_stations(sn);
S = sn.nservers;

% NHPP service is read as a flow, mirroring SOLVER_FLUID
if isfield(sn,'procid')
    for ist = 1:M
        for k = 1:K
            if sn.procid(ist,k) == ProcessType.NHPP && ~isempty(Mu{ist}{k}) ...
                    && ~isnan(Mu{ist}{k}(1))
                lam = Mu{ist}{k}(1);
                PH{ist}{k} = {-lam, lam};
                Phi{ist}{k} = 1;
            end
        end
    end
end

for ist = 1:M
    for k = 1:K
        if isnan(Mu{ist}{k})
            Mu{ist}{k} = [];
            Phi{ist}{k} = [];
        end
    end
    if isinf(S(ist))
        % a pure open model has no closed population, so N = 0 would make the
        % utilization divisor vanish
        S(ist) = max(N, 1);
    end
end

% state layout, mirroring SOLVER_FLUID_ODES
w = ones(M,K);
enabled = false(M,K);
q_indices = zeros(M,K);
Kic = zeros(M,K);
cumidx = 1;
for i = 1:M
    for c = 1:K
        if isempty(Mu{i}{c})
            numphases = 0;
            enabled(i,c) = false;
        else
            numphases = length(Mu{i}{c});
            enabled(i,c) = true;
        end
        q_indices(i,c) = cumidx;
        Kic(i,c) = numphases;
        cumidx = cumidx + numphases;
    end
end
for i = 1:M
    if sched(i) == SchedStrategy.DPS || sched(i) == SchedStrategy.GPS
        w(i,:) = sn.schedparam(i,:);
    end
end

lldscaling = [];
if isfield(sn,'lldscaling') && ~isempty(sn.lldscaling)
    lldscaling = sn.lldscaling;
end

D = ode_jumps_new(M, K, enabled, q_indices, rt, Kic);
[rateBase, eventIdx] = ode_rate_base(sn, Phi, Mu, PH, M, K, enabled, q_indices, rt, Kic, sched, D);

% THE MOMENT CLOSURE READS THE SAME REDUCED EVENT SET AS EVERY OTHER ROUTE.
% It used to refuse the reduction, on the grounds that it needs the untransformed
% event set; what it actually needs is to be able to say which (station,class)
% each event is a completion of, and EMAP carries exactly that across the
% composition -- an event folded through an immediate coordinate keeps a row with
% weight on every original event it stands for, including the two completions a
% pass-through realises at once. The diffusion D*diag(r)*D' is then the diffusion
% of the reduced process, which is the right one: the eliminated coordinate holds
% O(1/InfRate) mass and contributes noise of the same order.
% Kept at their pre-reduction values: the event ATTRIBUTES below (which
% (station,class) an event is a completion of) are indexed by original event,
% and EMAP is what carries them onto the reduced ones.
nevents0 = numel(rateBase);
eventIdx0 = eventIdx;
Emap = [];
absorb = [];
if fluid_hide_immediate(sn, options)
    [D, rateBase, eventIdx, ~, Emap, absorb] = ...
        ode_eliminate_immediate(D, rateBase, eventIdx, sn, options);
end

% a time-varying rate multiplier makes the drift non-autonomous, so the
% process has no stationary covariance to solve for
[~, rt_Mmat] = solver_fluid_ratemult(numel(rateBase), M, K, enabled, q_indices, Kic, Mu, eventIdx, options);
if ~isempty(rt_Mmat)
    line_error(mfilename,'Moment-closure methods require an autonomous drift, but options.config.rate_traj/nhpp_sched make the rates time-varying. Use options.method=''closing'' or ''matrix''.');
end

nstate = sum(Kic(:));

% Open and mixed models: the covariance lives on the QUEUE coordinates only.
%
% The closing representation models a Source as an EXT pseudo-station holding
% unit mass, so its coordinate is a normalisation constant, not a job count.
% Building D*diag(r)*D' over it would invent noise for a direction that has no
% population. Projecting those coordinates away leaves exactly the right open
% event set, because the closing form already emits the correct events: with a
% single-phase source the EXT rate factor is 1 - sum(of nothing) = 1
% identically, so an arrival is a CONSTANT-rate event whose jump, once the
% source row is dropped, is a lone +1 into the destination queue -- the
% canonical exogenous Poisson arrival, with diffusion intensity lambda. The
% return leg that LINE routes Sink -> Source becomes a lone -1, the departure.
% The EXT row of the Jacobian is identically zero for a single-phase source
% (FLUID_DRIFT_JACOBIAN, EXT branch), so A(covIdx,covIdx) is exactly the
% Jacobian of the projected drift, not an approximation of it.
%
% A MULTI-PHASE source is refused: those coordinates are the phase of ONE
% arrival process, a single Markov chain rather than a population, so its
% fluctuations are O(1) and the linear noise approximation does not apply to
% them at any scale.
isExt = (sched == SchedStrategy.EXT);
covMask = true(nstate,1);
for i = 1:M
    if ~isExt(i)
        continue
    end
    for c = 1:K
        if Kic(i,c) > 0
            if Kic(i,c) > 1
                line_error(mfilename,sprintf(['The moment-closure methods need a Poisson arrival stream, but the ' ...
                    'source of class %d is a %d-phase process. Those coordinates track the phase of a single ' ...
                    'arrival process rather than a population, so they carry no linear noise approximation. Use ' ...
                    'an exponential inter-arrival time, or options.method=''matrix''.'], c, Kic(i,c)));
            end
            covMask(q_indices(i,c):(q_indices(i,c)+Kic(i,c)-1)) = false;
        end
    end
end
covIdx = find(covMask);

stationBlock = cell(M,1);
classBlock = cell(M,K);
for i = 1:M
    lo = q_indices(i,1);
    hi = q_indices(i,K) + Kic(i,K) - 1;
    stationBlock{i} = lo:hi;
    for c = 1:K
        classBlock{i,c} = q_indices(i,c):(q_indices(i,c)+Kic(i,c)-1);
    end
end

% A STATION THAT CANNOT FILL ITS SERVERS HAS NOTHING TO CLOSE. min(n_i,c_i) is
% the identity on the whole support whenever the occupancy of station i is
% bounded above by its server count, and there the Gaussian closure is not an
% improvement on the first-order one, it is an ERROR: it spreads a normal
% marginal over n_i > c_i, mass the station can never hold, and returns
% E[min(n_i,c_i)] < n_i. On a closed model with one job per chain the exact
% answer is R = D at every queue (a job cannot queue behind itself), which the
% first-order closure reproduces to machine precision while the closure reads
% 0.4758 against 0.5 on the queue length and 1.1017 against 1 on the response
% time. The bound is the total population of every chain with a class enabled
% here; an open chain contributes Inf and never qualifies. SOLVER_FLUID_MOMENTS
% holds the drift variance of these stations at zero, exactly as it does for
% the delay stations, whose min() is likewise absent.
minExact = false(M,1);
if isfield(sn,'chains') && ~isempty(sn.chains) && isfield(sn,'njobs')
    nchains = size(sn.chains,1);
    for i = 1:M
        if isExt(i) || sched(i) == SchedStrategy.INF || ~isfinite(sn.nservers(i))
            continue
        end
        bound = 0;
        covered = false(1,K);
        for ch = 1:nchains
            inch = find(sn.chains(ch,:));
            covered(inch) = true;
            % VISITS, not `enabled`: a station may declare a service time for
            % every class while the routing never sends most of them there, and
            % counting those chains inflates the bound past the server count.
            % sn.visits is indexed by STATEFUL node, not by station.
            if numel(sn.visits) >= ch && ~isempty(sn.visits{ch})
                here = any(sn.visits{ch}(sn.stationToStateful(i),inch) > 0);
            else
                here = any(enabled(i,inch));
            end
            if here
                bound = bound + sum(sn.njobs(inch));
            end
        end
        if any(enabled(i,:) & ~covered)
            continue % a class outside every chain carries no population bound
        end
        minExact(i) = isfinite(bound) && bound <= sn.nservers(i) + GlobalConstants.FineTol;
    end
end

% classify events so throughputs can be read off the rate vector: ODE_RATE_BASE
% emits every service completion first, then every intra-PH phase change, and
% EVENTIDX records the source coordinate of each event. Summing the completion
% rates sourced at (i,c) gives the class-c throughput at station i exactly,
% because the routing probabilities and the PH entry vector each sum to one
% over the destinations enumerated there.
nDeparture = 0;
for i = 1:M
    for c = 1:K
        if enabled(i,c)
            for j = 1:M
                for l = 1:K
                    if rt((i-1)*K+c,(j-1)*K+l) > 0
                        nDeparture = nDeparture + Kic(i,c)*Kic(j,l);
                    end
                end
            end
        end
    end
end
% EVISDEPARTURE, EVSTATION and EVCLASS are indexed by ORIGINAL event, which is
% what EMAP maps onto. Without a reduction EMAP is the identity and the two
% indexings coincide, exactly as before.
nevents = nevents0;
evIsDeparture = false(nevents,1);
evIsDeparture(1:nDeparture) = true;
coordStation = zeros(nstate,1);
coordClass = zeros(nstate,1);
for i = 1:M
    for c = 1:K
        idx = q_indices(i,c):(q_indices(i,c)+Kic(i,c)-1);
        coordStation(idx) = i;
        coordClass(idx) = c;
    end
end
evStation = coordStation(eventIdx0);
evClass = coordClass(eventIdx0);

terms = struct();
terms.M = M; terms.K = K; terms.nstate = nstate;
terms.evIsDeparture = evIsDeparture;
% Identity when nothing was eliminated, so a caller writes one expression for
% both cases: TN = r' * (Emap * indicator_over_original_events).
if isempty(Emap)
    Emap = speye(numel(rateBase), nevents);
end
terms.Emap = Emap;
terms.absorb = absorb;
terms.evStation = evStation;
terms.evClass = evClass;
terms.D = D;
terms.rateBase = rateBase;
terms.eventIdx = eventIdx;
terms.q_indices = q_indices;
terms.Kic = Kic;
terms.enabled = enabled;
terms.w = w;
terms.sched = sched;
terms.S = S;
terms.stationBlock = stationBlock;
terms.classBlock = classBlock;
terms.minExact = minExact;
terms.covIdx = covIdx;
terms.isExt = isExt;

terms.lldscaling = lldscaling;
% every handle takes (x, sigma2, covblk): sigma2 closes the min(), covblk
% closes the DPS share ratio, and both are supplied by SOLVER_FLUID_MOMENTS
terms.factorFcn = @(x,s2,cb) ode_rates_closing_factors(x(:), M, K, enabled, q_indices, Kic, S, w, sched, s2, lldscaling, cb);
terms.ratesFcn = @(x,s2,cb) ode_rates_closing(x(:), M, K, enabled, q_indices, Kic, S, w, sched, rateBase, eventIdx, s2, lldscaling, cb);
terms.driftFcn = @(x,s2,cb) D * terms.ratesFcn(x,s2,cb);
terms.jacFcn = @(x,s2,cb) fluid_drift_jacobian(x(:), M, K, enabled, q_indices, Kic, S, w, sched, rateBase, eventIdx, D, s2, lldscaling, cb);
end
