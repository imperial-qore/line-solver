function [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t,iters,runtime] = solver_fluid_matrix(sn, options)

% [QN,UN,RN,TN,CN,RUNTIME] = SOLVER_FLUID_MATRIX(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;    %number of stations
K = sn.nclasses;    %number of classes
pie = sn.pie;
PH = sn.proc;
P_full = sn.rt;  % Full routing matrix (stateful nodes)
NK = sn.njobs';  %initial population
S = sn.nservers;
infServers = isinf(S);
S(infServers) = sum(NK);
nphases = sn.phases;
%refstat = sn.refstat; % reference station
weights = ones(M,K);

% Extract station-to-station routing matrix from stateful-to-stateful matrix
% Build mapping of station indices in the class-blocked routing matrix
station_indices = [];
for ist = 1:M
    isf = sn.stationToStateful(ist);
    for r = 1:K
        station_indices = [station_indices, (isf-1)*K + r];
    end
end
P = P_full(station_indices, station_indices);

% Remove Sink->Source feedback routing for open classes
% In open networks, jobs exit at Sink and should not recirculate back to Source.
% The routing matrix includes this feedback (added by getRoutingMatrix.m for
% pseudo-closed network analysis), but it causes incorrect flow balance in the
% fluid ODE formulation because arrivals are already accounted for via Alambda.
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

% Build arrival rate vector Aλ for mixed/open networks (Ruuskanen et al., PEVA 151 (2021), Eq. 7)
% Following the ground truth implementation: Source is conceptually excluded from state space.
% Arrivals go directly to queue phases (not Source phases) weighted by routing probabilities.
% This matches dx = W' * theta + A * lambda where lambda represents arrivals INTO queues.

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

% Eliminate immediate transitions if requested
state_map_imm = [];
if isfield(options.config, 'hide_immediate') && options.config.hide_immediate
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
                init_sol_idx = init_sol_idx + 1;
                Qa(1,state) = ist;
                SQC((ist-1)*K+r,state) = 1;
                SUC((ist-1)*K+r,state) = 1/S(ist);
                STC((ist-1)*K+r,state) = sum(sn.proc{ist}{r}{2}(k,:));
                x0_build(state,1) = options.init_sol(init_sol_idx);
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

    Qa = Qa(state_map_imm);
    SQC = SQC(:, state_map_imm);
    SUC = SUC(:, state_map_imm);
    STC = STC(:, state_map_imm);
    x0 = x0(state_map_imm);
    Alambda = Alambda(state_map_imm);  % Filter arrival vector for eliminated states
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
trange = [timespan(1),min(timespan(2),abs(10*itermax/min(nonZeroRates)))];

% Check if p-norm smoothing should be used (pstar parameter)
use_pnorm = isfield(options, 'pstar') && ~isempty(options.pstar) || ...
            (isfield(options.config, 'pstar') && ~isempty(options.config.pstar));

if use_pnorm
    % Get pstar values - expand scalar to per-station array
    if isfield(options, 'pstar') && ~isempty(options.pstar)
        pstar_val = options.pstar;
    else
        pstar_val = options.config.pstar;
    end
    if isscalar(pstar_val)
        pstar_val = pstar_val * ones(M, 1);
    end
    % Create per-phase pstar array (pQa) using the filtered Qa mapping
    pQa = pstar_val(Qa(:));  % Use Qa which is already filtered by keep and state_map_imm
    Sa_pnorm = S(Qa(:));  % Column vector for pnorm_ode
end

T0 = tic;
iters = 1;
if use_pnorm
    % p-norm smoothing ODE as per Ruuskanen et al., PEVA 151 (2021)
    % ghat = 1 / (1 + (x/c)^p)^(1/p) where x is queue length, c is servers, p is pstar
    % dx/dt = W^T * θ̂(x,p) + Aλ  (Eq. 27 for mixed networks)
    ode_pnorm_func = @(t,x) pnorm_ode(x, W, SQ, Sa_pnorm, pQa, Alambda, isSourceState);
    if options.stiff
        [t, xvec_t] = ode_solve_stiff(ode_pnorm_func, trange, x0, odeopt, options);
    else
        [t, xvec_t] = ode_solve(ode_pnorm_func, trange, x0, odeopt, options);
    end
else
    % Standard matrix method without smoothing
    % dx/dt = W^T * θ(x) + Aλ  (Eq. 12 for mixed networks)
    Sa_ode = S(Qa(:));  % Column vector for element-wise operations
    % Define theta function with special handling for Source stations
    theta_func = @(x) compute_theta(x, SQ, Sa_ode, isSourceState);
    if options.stiff
        [t, xvec_t] = ode_solve_stiff(@(t,x) W'*theta_func(x) + Alambda, trange, x0, odeopt, options);
    else
        [t, xvec_t] = ode_solve(@(t,x) W'*theta_func(x) + Alambda, trange, x0, odeopt, options);
    end
end
runtime = toc(T0);

Tmax = size(xvec_t,1);
QNtmp = cell(1,Tmax);
UNtmp = cell(1,Tmax);
RNtmp = cell(1,Tmax);
TNtmp = cell(1,Tmax);
Sa = S(Qa(:));  % Column vector for element-wise operations
S = repmat(S,1,K)'; S=S(:);
for j=1:Tmax
    x = xvec_t(j,:)';
    QNtmp{j} = zeros(K,M);
    TNtmp{j} = zeros(K,M);
    UNtmp{j} = zeros(K,M);
    RNtmp{j} = zeros(K,M);

    QNtmp{j}(:) = SQC*x;
    % Use compute_theta for consistent handling of Source stations
    theta_j = compute_theta(x, SQ, Sa, isSourceState);
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

xvec_it = {xvec_t(end,:)};
end

function dxdt = pnorm_ode(x, W, SQ, Sa, pQa, Alambda, isSourceState)
% PNORM_ODE - ODE derivative using p-norm smoothing
% As per Ruuskanen et al., PEVA 151 (2021)
% dxdt = W' * (x .* ghat) + Aλ  (Eq. 27)
% where ghat = 1 / (1 + (sumXQa/Sa)^pQa)^(1/pQa)

sumXQa = GlobalConstants.FineTol + SQ * x;
ghat = zeros(size(x));
for i = 1:length(x)
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
    else
        ghat(i) = 1;
    end
end

% Compute effective rate (x .* ghat), with special handling for Source stations
theta_eff = x .* ghat;
% For Source stations, override to 0.0 to bypass Source in dynamics
% (matching ground truth where Source is excluded from state space)
theta_eff(isSourceState) = 0.0;

dxdt = W' * theta_eff + Alambda;
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