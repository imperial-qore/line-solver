function [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t,iters,runtime] = solver_fluid_matrix(sn, options)

% [QN,UN,RN,TN,CN,RUNTIME] = SOLVER_FLUID_MATRIX(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;    %number of stations
K = sn.nclasses;    %number of classes
pie = sn.pie;
PH = sn.proc;
% see _kb/06-solver-catalog.md for rationale
if isfield(sn,'procid')
    for ist = 1:M
        for r = 1:K
            if sn.procid(ist,r) == ProcessType.NHPP
                lam = sn.rates(ist,r);
                if isfinite(lam) && lam > 0
                    PH{ist}{r} = {-lam, lam};
                    pie{ist}{r} = 1;
                end
            elseif sn.procid(ist,r) == ProcessType.MAPT || sn.procid(ist,r) == ProcessType.PHT
                % see solver_fluid.m: the nominal preserves the phase count
                [D0bar, D1bar] = sn_schedule_nominal(sn, ist, r);
                PH{ist}{r} = {D0bar, D1bar};
                pie{ist}{r} = map_pie({D0bar, D1bar});
            end
        end
    end
end
NK = sn.njobs';  %initial population
S = sn.nservers;
infServers = isinf(S);
S(infServers) = sum(NK);
nphases = sn.phases;
%refstat = sn.refstat; % reference station
weights = ones(M,K);

% Extract station-to-station routing matrix from stateful-to-stateful matrix
% using stochastic complementation to resolve routing through non-station
% stateful nodes (e.g., Router nodes). The closing family needs the same
% reduction, so it lives in one place.
P = sn_rt_stations(sn);

% see _kb/06-solver-catalog.md for rationale
for src_ist = 1:M
    if sn.sched(src_ist) == SchedStrategy.EXT
        % This is a Source station - remove feedback routing TO it for open classes
        for r = 1:K
            % Check if this is an open class with external arrivals
            if ~isnan(sn.rates(src_ist, r)) && sn.rates(src_ist, r) > 0
                % Zero out routing TO this Source for this class from all other stations
                src_col = (src_ist - 1) * K + r;  % Column index in P for (Source, class r)
                for from_ist = 1:M
                    if from_ist ~= src_ist  % Don't modify Source's own outgoing routing
                        for from_r = 1:K
                            from_row = (from_ist - 1) * K + from_r;
                            P(from_row, src_col) = 0;  % Remove feedback to Source
                        end
                    end
                end
            end
        end
    end
end

% ODE building as per Ruuskanen et al., PEVA 151 (2021).
Psi = [];
A = [];
B = [];
for ist=1:M
    for r=1:K
        if nphases(ist,r)==0
            Psi = blkdiag(Psi,0);
            B = blkdiag(B,0);
            A = blkdiag(A,NaN);
        else
            Psi = blkdiag(Psi,PH{ist}{r}{1});
            B = blkdiag(B,sum(PH{ist}{r}{2},2));
            A = blkdiag(A,pie{ist}{r}');
        end
    end
end
W = Psi + B*P*A';

% see _kb/06-solver-catalog.md for rationale

% First, identify Source stations and their arrival rates per class
source_arrivals = zeros(M, K);  % source_arrivals(src, r) = arrival rate at source src for class r
for src_ist = 1:M
    if sn.sched(src_ist) == SchedStrategy.EXT
        for r = 1:K
            if ~isnan(sn.rates(src_ist, r)) && sn.rates(src_ist, r) > 0
                source_arrivals(src_ist, r) = sn.rates(src_ist, r);
            end
        end
    end
end

% Build Alambda: arrivals go to QUEUE phases (where jobs route from Source), not Source phases
Alambda_full = zeros(size(A,1), 1);
state = 0;
for ist=1:M
    for r=1:K
        if nphases(ist,r) > 0
            if sn.sched(ist) == SchedStrategy.EXT
                % Source station: do NOT add arrivals here (arrivals go to downstream queues)
                state = state + nphases(ist,r);
            else
                % Queue station: check if it receives arrivals from any Source
                arrival_rate_to_queue = 0;
                for src_ist = 1:M
                    if source_arrivals(src_ist, r) > 0
                        % Get routing probability from Source to this queue for this class
                        src_row = (src_ist - 1) * K + r;
                        queue_col = (ist - 1) * K + r;
                        routing_prob = P(src_row, queue_col);
                        arrival_rate_to_queue = arrival_rate_to_queue + source_arrivals(src_ist, r) * routing_prob;
                    end
                end

                if arrival_rate_to_queue > 0
                    % Apply arrivals according to entrance probability ζ^{i,r}
                    for k=1:nphases(ist,r)
                        state = state + 1;
                        Alambda_full(state) = pie{ist}{r}(k) * arrival_rate_to_queue;
                    end
                else
                    state = state + nphases(ist,r);
                end
            end
        else
            % Add placeholder for disabled class to match W matrix structure
            state = state + 1;
            % Alambda_full(state) stays 0 since there are no arrivals
        end
    end
end

% remove disabled transitions
keep = find(~isnan(sum(W,1)));
W = W(keep,:);
W = W(:,keep);
Alambda = Alambda_full(keep);  % Also filter arrival vector

% see _kb/06-solver-catalog.md for rationale
hide_imm_requested = fluid_hide_immediate(sn, options);

% Eliminate immediate transitions if requested or auto-detected
state_map_imm = [];
npre_elim = [];
W_pre_elim = [];
Alambda_pre_elim = [];
if hide_imm_requested
    W_pre_elim = W;  % Save for Alambda correction
    Alambda_pre_elim = Alambda;
    [W, state_map_imm] = eliminate_immediate_matrix(W, sn, options);
end

Qa = []; % state mapping to queues (called Q(a) in Ruuskanen et al.)
SQC = zeros(M*K,0); % to compute per-class queue length at the end
SUC = zeros(M*K,0); % to compute per-class utilizations at the end
STC = zeros(M*K,0); % to compute per-class throughput at the end
x0_build = []; % Build x0 with same structure as W matrix
%x0 = []; % initial state
state = 0;
init_sol_idx = 0;  % Index into options.init_sol for enabled classes
for ist=1:M
    for r=1:K
        if nphases(ist,r)==0
            % Add placeholder for disabled transition (matching W matrix structure)
            state = state + 1;
            Qa(1,state) = ist; %#ok<*AGROW>
            SQC((ist-1)*K+r,state) = 0;  % No queue contribution
            SUC((ist-1)*K+r,state) = 0;  % No utilization contribution
            STC((ist-1)*K+r,state) = 0;  % No throughput contribution
            x0_build(state,1) = 0;  % No initial population
        else
            for k=1:nphases(ist,r)
                state = state + 1;
                Qa(1,state) = ist;
                if isnan(sn.rates(ist,r))
                    % Class has phases but is disabled (NaN rate) -
                    % solver_fluid_initsol skips these, so do not
                    % increment init_sol_idx
                    SQC((ist-1)*K+r,state) = 0;
                    SUC((ist-1)*K+r,state) = 0;
                    STC((ist-1)*K+r,state) = 0;
                    x0_build(state,1) = 0;
                else
                    init_sol_idx = init_sol_idx + 1;
                    SQC((ist-1)*K+r,state) = 1;
                    SUC((ist-1)*K+r,state) = 1/S(ist);
                    STC((ist-1)*K+r,state) = sum(sn.proc{ist}{r}{2}(k,:));
                    x0_build(state,1) = options.init_sol(init_sol_idx);
                end
            end
        end
        % code to initialize all jobs at ref station
        %if i == refstat(r)
        %    x0 = [x0; NK(r)*pie{i}{r}']; % initial probability of PH
        %else
        %    x0 = [x0; zeros(nphases(i,r),1)];
        %end
    end
end
x0 = x0_build;

% Apply keep filtering for disabled transitions (NaN in W matrix)
Qa = Qa(keep);
SQC = SQC(:, keep);
SUC = SUC(:, keep);
STC = STC(:, keep);
x0 = x0(keep);

% Save full matrices before immediate elimination for metric reconstruction
Qa_full = Qa;
STC_full = STC;
imm_states_in_keep = [];  % Indices of immediate states within keep-filtered space

% Apply state mapping if immediate transitions were eliminated
if ~isempty(state_map_imm)
    % Identify eliminated (immediate) states
    imm_states_in_keep = setdiff(1:length(Qa), state_map_imm);
    npre_elim = length(Qa);

    Qa = Qa(state_map_imm);

    % THE READ-OFF MATRICES ARE CORRECTED, NOT TRUNCATED. An eliminated
    % coordinate holds O(1/InfRate) mass and yet carries a FINITE throughput,
    % because the rate read off it is InfRate itself; dropping its column would
    % silently delete every completion the instantaneous phase makes and the
    % station's throughput would stop balancing against its neighbours. The mass
    % it holds per unit mass of the timed states is x_I = x_T*Q_TI*(-Q_II)^-1,
    % so folding that in is exact.
    Gsoj = W_pre_elim(state_map_imm, imm_states_in_keep) / ...
        (-W_pre_elim(imm_states_in_keep, imm_states_in_keep));
    SQC = SQC(:, state_map_imm) + SQC(:, imm_states_in_keep) * Gsoj';
    SUC = SUC(:, state_map_imm) + SUC(:, imm_states_in_keep) * Gsoj';
    STC = STC(:, state_map_imm) + STC(:, imm_states_in_keep) * Gsoj';

    % THE INITIAL POINT IS PROJECTED, NOT TRUNCATED. Dropping the rows would
    % delete whatever mass the initial condition parked on an eliminated
    % coordinate -- zero on a cold start, which puts every job in phase 1, but
    % not on a warm start from an earlier LN iterate -- and those are jobs, so
    % the chain would lose population before the first step. (-Q_II)\Q_IT is the
    % absorption distribution of the immediate block, the same complement the
    % arrival correction below uses.
    if ~isempty(imm_states_in_keep) && any(x0(imm_states_in_keep) ~= 0)
        Aimm = (-W_pre_elim(imm_states_in_keep, imm_states_in_keep)) \ ...
            W_pre_elim(imm_states_in_keep, state_map_imm);
        x0 = x0(state_map_imm) + Aimm' * x0(imm_states_in_keep);
    else
        x0 = x0(state_map_imm);
    end

    % see _kb/06-solver-catalog.md for rationale
    lambda_T = Alambda_pre_elim(state_map_imm);
    lambda_I = Alambda_pre_elim(imm_states_in_keep);
    if any(lambda_I ~= 0)
        Q_IT = W_pre_elim(imm_states_in_keep, state_map_imm);
        Q_II = W_pre_elim(imm_states_in_keep, imm_states_in_keep);
        Alambda = lambda_T - Q_IT' * (Q_II' \ lambda_I);
    else
        Alambda = lambda_T;
    end
end

% Build SQ matrix to compute total queue length per station in ODEs
% SQ(s,:) sums all states at the same station as state s
nstates = length(x0);
SQ = zeros(nstates, nstates);
for s = 1:nstates
    ist = Qa(s);
    SQ(s, Qa == ist) = 1;  % Sum all states at the same station
end

% Identify Source station states (EXT scheduler)
% For Source stations, theta should be 0.0 to effectively bypass Source in dynamics.
% This matches the ground truth implementation where Source is excluded from state space.
% Arrivals are injected directly into queue phases via Alambda.
isSourceState = false(nstates, 1);
for s = 1:nstates
    ist = Qa(s);
    if sn.sched(ist) == SchedStrategy.EXT
        isSourceState(s) = true;
        x0(s) = 0;  % Initialize Source phases to 0 (no mass at Source)
    end
end

%x0

tol = options.tol;
timespan = options.timespan;
itermax = options.iter_max;
odeopt = odeset('AbsTol', tol, 'RelTol', tol, 'NonNegative', 1:length(x0));
nonZeroRates = abs(W(abs(W)>0)); nonZeroRates=nonZeroRates(:);
if isempty(nonZeroRates)
    trange = [timespan(1), timespan(2)];
    if ~isfinite(trange(2))
        trange(2) = 1;
    end
else
    trange = [timespan(1),min(timespan(2),abs(10*itermax/min(nonZeroRates)))];
end

% Whether the p-norm smoothing is in force, and with which exponent. FLUID_PSTAR
% is the single rule: options.config.pstar switches it on under any name, and the
% method 'pnorm' defaults it to 20, as SolverFluid.java and the C++ FluidOptions
% do. Selecting it on the option alone made 'pnorm' run the unsmoothed drift.
[use_pnorm, pstar_val] = fluid_pstar(options.method, options, M);

if use_pnorm
    % Create per-phase pstar array (pQa) using the filtered Qa mapping
    pQa = pstar_val(Qa(:));  % Use Qa which is already filtered by keep and state_map_imm
    Sa_pnorm = S(Qa(:));  % Column vector for pnorm_theta
    % An INF station has k = infinity, so its share is 1 with no min() to
    % smooth; S holds sum(NK) there, which would smooth it as a k = N queue
    isInfState = infServers(Qa(:));
end

% One theta for the drift and for the metrics: eq. (23) reads the utilization
% off the SAME share the ODE integrated, so U and T must not revert to min()
if use_pnorm
    % p-norm smoothing as per Ruuskanen et al., PEVA 151 (2021), eq. (26)-(27)
    % ghat = 1 / (1 + (x/c)^p)^(1/p) where x is queue length, c is servers, p is pstar
    theta_func = @(x) pnorm_theta(x, SQ, Sa_pnorm, pQa, isSourceState, isInfState);
else
    % Standard matrix method without smoothing, eq. (12)
    Sa_ode = S(Qa(:));  % Column vector for element-wise operations
    theta_func = @(x) compute_theta(x, SQ, Sa_ode, isSourceState);
end

T0 = tic;
iters = 1;
ode_failed = false;
try
    % dx/dt = W^T * theta(x) + A*lambda
    ode_func = @(t,x) W'*theta_func(x) + Alambda;
    if options.stiff
        [t, xvec_t] = ode_solve_stiff(ode_func, trange, x0, odeopt, options);
    else
        [t, xvec_t] = ode_solve(ode_func, trange, x0, odeopt, options);
    end
catch me
    if contains(me.identifier, 'lsoda')
        ode_failed = true;
    else
        rethrow(me);
    end
end

% On LSODA failure, retry with hide_immediate toggled (only once)
if ode_failed
    is_retry = isfield(options.config, 'lsoda_retry') && options.config.lsoda_retry;
    if ~is_retry
        options_retry = options;
        options_retry.config.lsoda_retry = true;
        if hide_imm_requested
            % hide_immediate was active (auto-detected or explicit) and failed — retry without it
            options_retry.config.hide_immediate = false;
        else
            % Normal solve failed — retry with hide_immediate
            options_retry.config.hide_immediate = true;
        end
        [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t,iters,runtime] = solver_fluid_matrix(sn, options_retry);
        return;
    else
        % Both attempts failed — return empty lastSol to signal failure
        warning('lsoda:failed', 'LSODA failed on both normal and hide_immediate attempts');
        QN = NaN(M, K);
        UN = NaN(M, K);
        RN = NaN(M, K);
        TN = NaN(M, K);
        xvec_it = {};  % empty signals failure to caller (runAnalyzer checks isempty)
        QNt = cell(M, K);
        UNt = cell(M, K);
        TNt = cell(M, K);
        xvec_t = zeros(1, sum(nphases(:)));
        t = 0;
        iters = 0;
        runtime = toc(T0);
        return;
    end
end
runtime = toc(T0);

% DEGENERATE DRIFT: re-integrate with a closed saturation term, do not touch the
% answer that came back. min(E[n],c) is FLAT above the server count, so a network
% of saturated stations has a CONTINUUM of fixed points and this method returns
% whichever one the integrator stopped at -- [9 1] against an exact [5 5] on two
% identical saturated stations in a closed cycle, and [8 2] with two servers each.
% The repair is applied to the DRIFT, not to the point: the same trajectory is
% integrated again with E[min(n,c)] in place of min(E[n],c), which is strictly
% increasing and so isolates one fixed point. A selection rule imposed after the
% fact would not be a solution of anything.
%
% WHY A CLOSURE AND NOT A SMOOTHED min: any smoothing sharp enough to stay
% faithful to min away from the kink is numerically FLAT far from it. The
% Boltzmann softmin at alpha=20 carries a restoring force of exp(-160) at the
% [9 1] point, and the p-norm trades the two off directly (pstar=2 recovers
% [5 5], pstar=8 gives [7.64 2.36], pstar=128 gives [8.94 1.06]). The closure
% escapes the trade-off because its slope comes from the VARIANCE of the
% marginal rather than from a smoothing width.
%
% Only a model that is ACTUALLY degenerate pays for it: the test is a
% null-direction probe at the returned point, so a well-posed model integrates
% once and is unchanged. See BUGS.md and _kb/06-solver-catalog.md.
if ~use_pnorm && ~ode_failed && ~isempty(xvec_t)
    isInfState_ode = infServers(Qa(:));
    if fluid_degenerate_fixed_point(xvec_t(end,:)', ode_func, SQC, K, ...
            isSourceState, isInfState_ode, max(abs(W(:))))
        theta_func = @(x) compute_theta_closed(x, SQ, Sa_ode, isSourceState, isInfState_ode);
        ode_func_c = @(t,x) W'*theta_func(x) + Alambda;
        try
            if options.stiff
                [t_c, xvec_c] = ode_solve_stiff(ode_func_c, trange, x0, odeopt, options);
            else
                [t_c, xvec_c] = ode_solve(ode_func_c, trange, x0, odeopt, options);
            end
            if ~isempty(xvec_c) && all(isfinite(xvec_c(end,:)))
                t = t_c;
                xvec_t = xvec_c;
                line_printf('Fluid: the first-order fixed point is not isolated (two or more saturated stations), so the drift was re-integrated with a closed saturation term.\n');
            end
        catch
            % A failed repair leaves the unrepaired answer standing rather than
            % turning a wrong number into no number.
        end
    end
end

Tmax = size(xvec_t,1);
QNtmp = cell(1,Tmax);
UNtmp = cell(1,Tmax);
RNtmp = cell(1,Tmax);
TNtmp = cell(1,Tmax);
S = repmat(S,1,K)'; S=S(:);
for j=1:Tmax
    x = xvec_t(j,:)';
    QNtmp{j} = zeros(K,M);
    TNtmp{j} = zeros(K,M);
    UNtmp{j} = zeros(K,M);
    RNtmp{j} = zeros(K,M);

    QNtmp{j}(:) = SQC*x;
    % The same share the drift used, smoothed or not
    theta_j = theta_func(x);
    TNtmp{j}(:) = STC*theta_j;
    UNtmp{j}(:) = SUC*theta_j;
    % Little's law is invalid in transient so this vector is not returned
    % except the last element as an approximation of the actual RN
    RNtmp{j}(:) = QNtmp{j}(:)./TNtmp{j}(:);

    QNtmp{j} = QNtmp{j}';
    UNtmp{j} = UNtmp{j}';
    RNtmp{j} = RNtmp{j}';
    TNtmp{j} = TNtmp{j}';
end
% steady state metrics
for j=1:Tmax
    QNtmp{j} = QNtmp{j}(:);
    UNtmp{j} = UNtmp{j}(:);
    RNtmp{j} = RNtmp{j}(:);
    TNtmp{j} = TNtmp{j}(:);
end

% compute cell array with time-varying metrics for stations and classes
QNtmp = cell2mat(QNtmp)';
UNtmp = cell2mat(UNtmp)';
RNtmp = cell2mat(RNtmp)';
TNtmp = cell2mat(TNtmp)';
QNt = cell(M,K);
UNt = cell(M,K);
RNt = cell(M,K);
TNt = cell(M,K);
for ist=1:M
    for r=1:K
        QNt{ist,r} = QNtmp(:,(r-1)*M+ist);
        UNt{ist,r} = UNtmp(:,(r-1)*M+ist);
        RNt{ist,r} = RNtmp(:,(r-1)*M+ist);
        TNt{ist,r} = TNtmp(:,(r-1)*M+ist);
    end
end
QN = reshape(QNtmp(end,:),M,K);
UN = reshape(UNtmp(end,:),M,K);
RN = reshape(RNtmp(end,:),M,K);
TN = reshape(TNtmp(end,:),M,K);

% Set throughput at Source stations to arrival rates for open classes
% Source stations have theta = 0 in the ODE (to bypass Source in dynamics),
% but their throughput should equal the external arrival rate.
for ist=1:M
    if sn.sched(ist) == SchedStrategy.EXT
        for r=1:K
            if ~isnan(sn.rates(ist, r)) && sn.rates(ist, r) > 0
                TN(ist, r) = sn.rates(ist, r);
            end
        end
    end
end

% Reconstruct throughput for eliminated immediate states using flow conservation
if ~isempty(imm_states_in_keep)
    % For each eliminated state, throughput = sum of incoming throughputs via routing
    for idx = 1:length(imm_states_in_keep)
        s = imm_states_in_keep(idx);
        ist = Qa_full(s);  % Station of this eliminated state
        % Find which class this state belongs to
        r = 0;
        state_count = 0;
        for rr = 1:K
            state_count = state_count + nphases(ist, rr);
            if s <= sum(Qa_full == ist & (1:length(Qa_full)) <= state_count)
                r = rr;
                break;
            end
        end
        if r == 0
            % Fallback: find class by examining STC_full
            for rr = 1:K
                if STC_full((ist-1)*K+rr, s) > 0
                    r = rr;
                    break;
                end
            end
        end
        if r > 0
            % Compute incoming throughput using routing matrix P
            % T(ist, r) = sum over all (j, l) of T(j, l) * P((j,l) -> (ist,r))
            incoming_tput = 0;
            for j = 1:M
                for l = 1:K
                    p_jl_ir = P((j-1)*K+l, (ist-1)*K+r);
                    if p_jl_ir > 0
                        incoming_tput = incoming_tput + TN(j, l) * p_jl_ir;
                    end
                end
            end
            % Also add external arrivals if this is a source station
            if sn.sched(ist) == SchedStrategy.EXT && ~isnan(sn.rates(ist, r))
                incoming_tput = incoming_tput + sn.rates(ist, r);
            end
            TN(ist, r) = incoming_tput;
        end
    end
end

% XVEC_IT IS HANDED BACK IN THE PRE-ELIMINATION LAYOUT. The caller stores it as
% options.init_sol for the next FCFS iterate, which rebuilds the state over every
% phase, so a vector shortened by the immediate elimination would run that iterate
% off the end of init_sol. The eliminated coordinates come back empty, which is
% what they hold.
if ~isempty(state_map_imm) && numel(state_map_imm) == size(xvec_t,2)
    xend = zeros(1, npre_elim);
    xend(state_map_imm) = xvec_t(end,:);
    xvec_it = {xend};
else
    xvec_it = {xvec_t(end,:)};
end
end

function theta = pnorm_theta(x, SQ, Sa, pQa, isSourceState, isInfState)
% PNORM_THETA - Mass in service under p-norm smoothing
% As per Ruuskanen et al., PEVA 151 (2021), eq. (26)-(27)
% theta = x .* ghat, ghat = 1 / (1 + (sumXQa/Sa)^pQa)^(1/pQa)
% An INF station carries no min() to smooth, so its share stays 1.

sumXQa = GlobalConstants.FineTol + SQ * x;
ghat = ones(size(x));
for i = 1:length(x)
    if isInfState(i)
        continue
    end
    xVal = sumXQa(i);
    cVal = Sa(i);
    pVal = pQa(i);
    if cVal > 0 && pVal > 0
        ghatVal = 1.0 / (1 + (xVal / cVal)^pVal)^(1/pVal);
        if isnan(ghatVal)
            ghat(i) = 0;
        else
            ghat(i) = ghatVal;
        end
    end
end

theta = x .* ghat;
% For Source stations, override to 0.0 to bypass Source in dynamics
% (matching ground truth where Source is excluded from state space)
theta(isSourceState) = 0.0;
end

function theta = compute_theta_closed(x, SQ, Sa, isSourceState, isInfState)
% COMPUTE_THETA_CLOSED - Mass in service with a VARIANCE-CARRYING saturation term
%
% E[min(n,c)] under the station's equilibrium geometric marginal, rather than
% min(E[n],c):
%
%   n ~ Geometric(mean m)  =>  E[min(n,c)] = sum_{k=1..c} p^k = m*(1 - p^c),
%                              p = m/(1+m).
%
% It has the two properties the hard min lacks and the degeneracy repair needs:
% strictly increasing in m everywhere (slope 1/(1+m)^2 at c=1, so still 1e-2 at
% m=9, a restoring force the integrator can follow inside its horizon), and the
% same asymptotes, -> c as m -> inf and -> m as m -> 0. It is the first-order
% face of what the 'dae' rung does by seeding the variance positive, which is
% why both isolate the same fixed point.

sumXQa = GlobalConstants.FineTol + SQ * x;
p = sumXQa ./ (1 + sumXQa);
emin = sumXQa .* (1 - p .^ max(Sa, 0));
emin(~isfinite(emin)) = 0;
% An INF station has a server per job: there is no min() to close, and Sa holds
% the whole population there, which the closure would read as a finite queue.
emin(isInfState) = sumXQa(isInfState);
emin = min(emin, sumXQa);
theta = x ./ sumXQa .* emin;
theta(isSourceState) = 0.0;
end

function tf = fluid_degenerate_fixed_point(x, ode_func, SQC, K, isSourceState, isInfState, rateScale)
% FLUID_DEGENERATE_FIXED_POINT - Is the returned point one of a CONTINUUM?
%
% A station whose queue exceeds its server count has theta pinned at the server
% count: min(S,sum_x) stops depending on sum_x, so the drift cannot tell one
% split of the mass between two such stations from another. The test is direct:
% move a little mass of ONE CLASS from one station to another along a
% population-conserving direction and see whether the drift moves at all. Both
% directions are tried, because the integrator typically stops on the BOUNDARY
% of the degenerate set, where one of the two does change the drift.
%
% PER CLASS, NOT PER STATION. A direction that moves a station's mass across ALL
% its classes is not one the model can take: a SelfLoopingClass is pinned at one
% station and can never leave, so the direction is infeasible, the drift is
% trivially unchanged along it, and a well-posed model reads as degenerate. That
% is what it did to sanity_CQN_2q_psfcfs_1class_1slcateachqueue, whose two queues
% each hold a self-looping job: RespT came back 1.4336 against a baseline of
% 0.726303. Moving ONE class between two stations it occupies IS feasible, and a
% self-looping class occupies exactly one station, so no pair exists for it.
%
% The DIRECTIONAL DERIVATIVE is the scale-free quantity to threshold: a live
% direction moves the drift at the station's own service rate and a null one only
% by the FineTol the share carries, four orders apart.

tf = false;
d0 = ode_func(0, x);
if isempty(d0) || max(abs(d0)) > 1e-6 * max(1, max(abs(x)))
    return
end
n = numel(x);
M = size(SQC,1)/K;
step = 1e-3 * max(1, max(abs(x)));
rateScale = max(rateScale, 1e-12);
for r = 1:K
    groups = {};
    for i = 1:M
        members = find(SQC((i-1)*K+r,:) > 0);
        % Source and INF states are dropped: a Source carries theta = 0 by
        % construction and an INF station carries theta = x with no min() to pin.
        members = members(~isSourceState(members) & ~isInfState(members));
        if ~isempty(members)
            groups{end+1} = members(:);
        end
    end
    if numel(groups) < 2
        continue
    end
    mass = cellfun(@(g) sum(x(g)), groups);
    for a = 1:numel(groups)
        if mass(a) <= step
            continue
        end
        for b = 1:numel(groups)
            if a == b
                continue
            end
            ga = groups{a}; gb = groups{b};
            d = zeros(n,1);
            d(ga) = d(ga) - x(ga)/mass(a);          % take, proportionally
            if mass(b) > 0
                d(gb) = d(gb) + x(gb)/mass(b);      % give, proportionally
            else
                d(gb) = d(gb) + 1/numel(gb);
            end
            dd = ode_func(0, x + step*d) - d0;
            if max(abs(dd))/step <= 1e-4 * rateScale
                tf = true;
                return
            end
        end
    end
end
end

function theta = compute_theta(x, SQ, Sa, isSourceState)
% COMPUTE_THETA - Compute theta vector for fluid ODE
% For regular stations: theta = x./(SQ*x) .* min(Sa, SQ*x)
% For Source stations (EXT scheduler): theta = 0.0
%   - Source is conceptually excluded from state space (matching ground truth)
%   - Arrivals are injected directly into queues via Alambda
%   - Setting theta = 0 prevents Source from contributing to W' * theta

sumXQa = GlobalConstants.FineTol + SQ * x;
theta = x ./ sumXQa .* min(Sa, sumXQa);

% Override theta for Source stations
% Source stations have theta = 0.0 to bypass Source in dynamics
theta(isSourceState) = 0.0;
end