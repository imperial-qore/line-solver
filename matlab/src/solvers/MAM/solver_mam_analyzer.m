function [QN,UN,RN,TN,CN,XN,runtime,method,totiter,percResults] = solver_mam_analyzer(sn, options)
% [QN,UN,RN,TN,CN,XN,RUNTIME,METHOD] = SOLVER_MAM_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart = tic;

% Discrete-time (slotted) models are recognized from the distributions and
% routed to the Q-MAM discrete-time algorithms. The test must run BEFORE
% sn_nonmarkov_toph, which would fit a continuous CME to a Geometric and erase
% the lattice; see _kb/06-solver-catalog.md for the LAS-DA convention.
[isDiscreteTime, slotLength] = sn_is_discrete_time(sn, options);
if isDiscreteTime
    [QN,UN,RN,TN,CN,XN,totiter,method] = solver_mam_dt(sn, options, slotLength);
    percResults = [];
    runtime = toc(Tstart);
    return;
end

% Exact fast-path: correlated single-class MAP/MAP/1; see _kb/06-solver-catalog.md for rationale
% A finite buffer a closed class can fill is refused for every method: the
% analyzers below carry a buffer as a LOSS buffer and the closed routes read
% none, see MAM_BUFFER_REFUSAL (the predicate SolverMAM.supportsModelMethod
% reports, so this raise is reached only with enableChecks off).
bufferWhy = mam_buffer_refusal(sn);
if ~isempty(bufferWhy)
    line_error(mfilename, bufferWhy);
end
[isMapMap1, QN, UN, RN, TN, CN, XN] = solver_mam_mapmap1_exact(sn);
if isMapMap1
    method = 'exact.mapmap1';
    totiter = 1;
    percResults = [];
    runtime = toc(Tstart);
    return;
end

% Preserve deterministic distributions for exact MAP/D/c analysis
if ~isfield(options, 'config')
    options.config = struct();
end
if ~isfield(options.config, 'preserveDet')
    % The bgchain method builds a phase-type mixture of the open service laws
    % and has no deterministic-service branch, so it needs the Det fitted to PH
    % rather than preserved.
    options.config.preserveDet = ~strcmpi(options.method, 'bgchain');
end
% The conversion RETAGS procid to APH/ME/MAP, so afterwards nothing names the law
% the user declared. MMAP[K]/G[K]/1 works off that law's transform (sn.lst) rather
% than off the surrogate, so its gate needs the tags as they stand here; only DET
% survives the retagging.
procidDeclared = sn.procid;

% Convert non-Markovian distributions to PH (Det preserved if preserveDet=true)
sn = sn_nonmarkov_toph(sn, options);
sn.procidDeclared = procidDeclared;

% Check if the model is mixed (has both open and closed classes)
isOpen = sn_is_open_model(sn);
isClosed = sn_is_closed_model(sn);
% isOpen means EVERY class is open and isClosed that every class is closed, so
% a mixed model is neither, not both
isMixed = sn_has_mixed_classes(sn);
% Setup/delay-off is read by solver_mam_basic AND, since 2026-09, by
% solver_mam_ldqbd, which solves the closed regime exactly through
% qbd_setupdelayoff_closed. bgchain and mna still drop it, so the default must
% not hand a setup model to those two; the ldqbd branch decides for itself in
% ISCLOSEDDELAYQUEUE, which admits the regime that analysis covers (single
% server, exponential service, no load dependence) and refuses the rest.
hasSetup = isfield(sn, 'hassetup') && any(sn.hassetup);

line_debug('MAM analyzer starting: method=%s, isOpen=%d, isClosed=%d', options.method, isOpen, isClosed);

% Mixed models are supported by the dec.source method

if nargin<2 || isempty(options.config) || ~isfield(options.config,'merge')
    % Preserve existing config fields (like preserveDet) while setting defaults
    if ~isfield(options, 'config') || isempty(options.config)
        options.config = struct();
    end
    if ~isfield(options.config, 'merge')
        options.config.merge = 'super';
    end
    if ~isfield(options.config, 'compress')
        options.config.compress = 'mixture.order1';
    end
    if ~isfield(options.config, 'space_max')
        options.config.space_max = 128;
    end
    if ~isfield(options.config, 'etaqa_trunc')
        options.config.etaqa_trunc = 8;
    end
end

method = options.method;
percResults = []; % Initialize as empty

switch options.method
    case 'dec.mmap'
        % service distribution per class scaled by utilization used as
        % departure process
        line_debug('Using dec.mmap method, calling solver_mam');
        [QN,UN,RN,TN,CN,XN,totiter] = solver_mam(sn, options);
    case {'default', 'dec.source'}
        % Check if network is a valid Fork-Join topology for FJ_codes
        [isHomogeneous, fjInfo] = fj_is_homogeneous(sn);

        if isHomogeneous
            % Use FJ_codes for Fork-Join analysis
            if strcmpi(options.method, 'default')
                line_debug('Default method: using FJ_codes for Fork-Join topology\n');
            end
            line_debug('Detected Fork-Join topology, using FJ_codes method');
            [QN,UN,RN,TN,CN,XN,totiter,percResults] = solver_mam_fj(sn, options);
            method = 'qiu';
        elseif sn_has_fork_join(sn) && sn_is_open_model(sn)
            % Use MMAP-based FJ decomposition with mmap_max synchronization
            line_debug('Detected general Fork-Join topology, using dec.source.mmap method');
            [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_basic_mmap(sn, options);
            method = 'dec.source.mmap';
        else
            % Check if network is a valid BMAP/PH/N/N bufferless retrial queue
            [isRetrial, ~] = qsys_is_retrial(sn);

            % Check if network has reneging (queue abandonment) for MAPMsG
            isReneging = mam_has_reneging_patience(sn);

            if isRetrial
                % Use BMAP/PH/N/N retrial solver
                if strcmpi(options.method, 'default')
                    line_debug('Default method: using retrial for BMAP/PH/N/N bufferless topology\n');
                end
                line_debug('Detected BMAP/PH/N/N retrial topology, using retrial method');
                [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_retrial(sn, options);
                method = 'retrial';
            elseif isReneging
                % Use MAP/M/s+G reneging solver (MAPMsG)
                if strcmpi(options.method, 'default')
                    line_debug('Default method: using MAPMsG for MAP/M/s+G with reneging\n');
                end
                line_debug('Detected reneging topology, using MAPMsG method');
                [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_retrial(sn, options);
                method = 'reneging';
            elseif strcmpi(options.method, 'default') && isClosedDelayQueue(sn)
                % Single-class closed Delay+Queue: LDQBD is exact; see _kb/06-solver-catalog.md for rationale
                line_debug('Default method: using LDQBD for single-class closed Delay/Queue\n');
                [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_ldqbd(sn, options);
                method = 'ldqbd';
            elseif strcmpi(options.method, 'default') && isClosed && ~hasSetup ...
                    && mam_bgchain_applicable(sn, options) && bgchainClosedExact(sn)
                % A closed model is the degenerate case of the background chain:
                % with no open work to take a share of the servers the chain is
                % the EXACT closed CTMC at chain granularity, so it dominates the
                % mna fixed point wherever its state space fits. BGCHAINAPPLIES
                % sizes that state space and BGCHAINCLOSEDEXACT checks the chain
                % is built from a service law it represents exactly.
                line_debug('Default method: using bgchain for a closed model\n');
                [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_bgchain(sn, options);
                method = 'bgchain';
            elseif strcmpi(options.method, 'default') && isClosed && ~hasSetup && mnaApplies(sn)
                % A closed model has no arrival stream for dec.source to build
                % its Poisson surrogate from: it replaces each closed chain by a
                % source at the current throughput iterate and never enforces the
                % population, so the answer neither conserves N nor separates the
                % classes. MNA closes the same traffic equations by bisecting the
                % per-class throughput against N; see _kb/06-solver-catalog.md.
                line_debug('Default method: using mna for a closed model\n');
                [QN,UN,RN,TN,CN,XN,~,totiter] = solver_mna_closed(sn, options);
                method = 'mna';
                if ~mnaConserves(sn, QN, RN, TN)
                    % The outer bisection did not close on the population: the
                    % last step rescales each chain onto N regardless, so the
                    % failure is invisible in QN alone and only Little's law on
                    % the unrescaled R and T still shows it.
                    line_debug('mna did not close on the population, falling back to dec.source\n');
                    [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_basic(sn, options);
                    method = 'dec.source';
                end
            elseif strcmpi(options.method, 'default') && isMixed && ~hasSetup ...
                    && mam_bgchain_applicable(sn, options)
                % Mixed: the closed classes are solved exactly as a background
                % chain and the open ones as QBDs driven by it, which is 4-5
                % significant digits against CTMC where dec.source is 10-24% out
                line_debug('Default method: using bgchain for a mixed model\n');
                [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_bgchain(sn, options);
                method = 'bgchain';
            else
                % arrival process per chain rescaled by visits at each node
                if strcmpi(options.method, 'default')
                    line_debug('Default method: using dec.source\n');
                end
                line_debug('Using default/dec.source method, calling solver_mam_basic');
                [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_basic(sn, options);
                method = 'dec.source';
            end
        end
    case 'dec.poisson'
        % analyze the network with Poisson streams
        line_debug('Using dec.poisson method with space_max=1, calling solver_mam_basic');
        options.config.space_max = 1;
        [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_basic(sn, options);
    case 'mna'
        % The predicate SolverMAM.supportsModelMethod asks, so the gate that
        % decides whether to OFFER 'mna' and this run cannot drift apart.
        [mnaOk, mnaWhy] = mam_mna_applicable(sn);
        if ~mnaOk
            line_error(mfilename, mnaWhy);
        end
        if sn_is_open_model(sn)
            line_debug('Using MNA method for open model, calling solver_mna_open');
            [QN,UN,RN,TN,CN,XN,~,totiter] = solver_mna_open(sn, options);
        else
            line_debug('Using MNA method for closed model, calling solver_mna_closed');
            [QN,UN,RN,TN,CN,XN,~,totiter] = solver_mna_closed(sn, options);
        end
    case 'bgchain'
        % Mixed networks: the closed classes are a background modulating chain,
        % the open classes are QBDs driven by it. With several closed chains an
        % outer iteration tags one chain at a time and aggregates the rest, so
        % the background chain always carries two classes.
        line_debug('Using bgchain method, calling solver_mam_bgchain');
        [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_bgchain(sn, options);
    case 'ldqbd'
        % Level-Dependent QBD method for single-class Delay/Queue networks:
        % closed (finite population) or open (Poisson arrivals, truncated).
        % The shape is MAM_LDQBD_APPLICABLE's, which the gate and
        % solver_mam_ldqbd ask as well.
        [ldOk, ldWhy] = mam_ldqbd_applicable(sn);
        if ~ldOk
            line_error(mfilename, ldWhy);
        end
        line_debug('Using LDQBD method, calling solver_mam_ldqbd');
        [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_ldqbd(sn, options);
    case 'dec.source.mmap'
        % MMAP-based FJ decomposition with mmap_max synchronization
        line_debug('Using dec.source.mmap method, calling solver_mam_basic_mmap');
        [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_basic_mmap(sn, options);
    case 'retrial'
        % The BMAP/PH/N/N retrial analyzer, which 'default' also resolves to on
        % a retrial topology and on a reneging patience. Advertised by
        % listValidMethods, so it must be reachable by name; refused by name
        % off those topologies rather than answering a different model.
        % The predicate SolverMAM.supportsModelMethod asks, so the gate that
        % decides whether to OFFER 'retrial' and this run cannot drift apart.
        [retrialOk, retrialWhy] = mam_retrial_applicable(sn);
        if ~retrialOk
            line_error(mfilename, retrialWhy);
        end
        line_debug('Using retrial method, calling solver_mam_retrial');
        [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_retrial(sn, options);
    otherwise
        line_error(mfilename,'Unknown method: %s', options.method);
end

for i=1:sn.nstations
    switch sn.sched(i)
        case SchedStrategy.EXT
            TN(i,:) = sn.rates(i,:);
    end
end

QN(isnan(QN))=0;
CN(isnan(CN))=0;
RN(isnan(RN))=0;
UN(isnan(UN))=0;
XN(isnan(XN))=0;
TN(isnan(TN))=0;

runtime = toc(Tstart);
end

function isLdqbd = isClosedDelayQueue(sn)
% ISCLOSEDDELAYQUEUE Check if model is a single-class closed Delay+Queue
%
% The CLOSED regime of MAM_LDQBD_APPLICABLE, the predicate the gate and
% solver_mam_ldqbd ask. The open Source+Queue regime is deliberately excluded
% here: it truncates the level space at options.cutoff, so it is not
% unconditionally preferable to dec.source and stays opt-in.

[ok, ~, regime] = mam_ldqbd_applicable(sn);
isLdqbd = ok && strcmp(regime, 'closed');

% A SETUP/DELAY-OFF QUEUE is exact here only at a single server with
% exponential service and no load dependence, which is what
% qbd_setupdelayoff_closed models. Outside that, LDQBD refuses by name, so the
% default has to fall through to the decomposition rather than reach the
% refusal; dec.source carries the same analysis approximately.
if isLdqbd && isfield(sn,'hassetup') && any(sn.hassetup)
    qIdx = find(sn.sched == SchedStrategy.FCFS, 1);
    if sn.hassetup(qIdx)
        if sn.nservers(qIdx) > 1 || sn_has_load_dependence(sn)
            isLdqbd = false;
        else
            PHq = sn.proc{qIdx}{1};
            isLdqbd = numel(PHq{1}) == 1;   % exponential service only
        end
    end
end

end

function tf = mnaApplies(sn)
% MNAAPPLIES Check that the closed MNA analyzer covers this model as a DEFAULT
%
% MAM_MNA_APPLICABLE carries the correctness rules a by-name 'mna' must clear
% (mixed model, closed round-robin, self-looping class, fork-join, closed
% class switching, a discipline the sweep never updates). The two tests below
% are PREFERENCES: they do not refuse a caller who asks for mna, they only
% keep the default off shapes where dec.source is cheaper or better.

tf = false;

if ~mam_mna_applicable(sn)
    return;
end

for ist = 1:sn.nstations
    % The PS branch of solver_mna_closed forms U = S*T and the geometric bound
    % from it WITHOUT dividing by the number of servers, so a multiserver PS
    % station is misrepresented; dec.source is exact on the non-queueing regime
    % (c >> N) that shape usually stands for.
    if sn.sched(ist) == SchedStrategy.PS && sn.nservers(ist) > 1
        return;
    end
    % A multiclass FCFS station makes the flow sweep superpose one MMAP per
    % class and then solve MMAPPH1FCFS at level sum(N)+1: measured against an
    % exact CTMC that costs 12-61s where dec.source costs 0.03s and is not
    % more accurate (mean relative error 0.11-0.38 against 0.04-0.19). The
    % single-class case is both cheap and better, so keep only that one.
    if sn.sched(ist) == SchedStrategy.FCFS && sn.nclasses > 1
        return;
    end
end

tf = true;

end

function tf = mnaConserves(sn, QN, RN, TN)
% MNACONSERVES Check that the closed MNA outer bisection closed on N
%
% solver_mna_closed rescales each chain onto its population as a last step, so
% a diverged bisection still returns queue lengths that sum to N and the
% failure is invisible in QN. R and T are NOT rescaled, so Little's law over
% the whole network, sum_i T(i,k)*R(i,k) = N_k, still reads the raw iterate: a
% converged run lands within 5e-4 of N and a diverged one is orders of
% magnitude out, or negative.

tf = false;

if any(~isfinite(QN(:))) || any(~isfinite(RN(:))) || any(~isfinite(TN(:)))
    return;
end

Npred = sum(TN .* RN, 1);
for k = 1:sn.nclasses
    if ~isfinite(sn.njobs(k)) || sn.njobs(k) <= 0
        continue;
    end
    if abs(Npred(k) - sn.njobs(k)) > 0.01 * sn.njobs(k)
        return;
    end
end

tf = true;

end

function tf = bgchainClosedExact(sn)
% BGCHAINCLOSEDEXACT Check that the background chain represents this closed
% model's service laws exactly
%
% A station that is not PS or INF must satisfy BOTH conditions below, because the
% background chain makes two separate first-moment substitutions there.
%
%   1. MAM_BGCHAIN_CTMC builds its generator from the MEAN service time alone.
%      That is exact at a PS or INF station, which is insensitive to the service
%      law beyond its first moment, and exact under any discipline when the law
%      IS exponential. Measured on a closed Delay+FCFS cycle with Erlang-3
%      service, the mean-only chain reads 2.8% off SolverCTMC.
%   2. The capacity a station's closed jobs hold is split over the background
%      classes in proportion to their COUNTS, which is service in random order.
%      That is exact under PS, and exact under FCFS only when the classes are
%      served at the SAME rate -- an FCFS station with class-dependent rates
%      reads 25.2% off SolverCTMC on a two-chain closed cycle, against 3.5e-16
%      when the two rates are made equal.
%
% solver_mna_closed carries the phase-type representation instead, so neither
% surrogate may be chosen as the DEFAULT. Asking for bgchain by name still gets
% it, with both approximations documented in SOLVER_MAM_BGCHAIN.

tf = false;

for ist = 1:sn.nstations
    if sn.sched(ist) == SchedStrategy.INF || sn.sched(ist) == SchedStrategy.PS ...
            || sn.sched(ist) == SchedStrategy.EXT
        continue;
    end
    rate_here = NaN;
    for r = 1:sn.nclasses
        if ~isfinite(sn.rates(ist, r)) || sn.rates(ist, r) <= 0
            continue;
        end
        if sn.procid(ist, r) ~= ProcessType.EXP
            return;
        end
        if isnan(rate_here)
            rate_here = sn.rates(ist, r);
        elseif abs(sn.rates(ist, r) - rate_here) > GlobalConstants.CoarseTol * rate_here
            return;
        end
    end
end

tf = true;

end
