function sys = solver_fluid_symodes(sn, options)
% SYS = SOLVER_FLUID_SYMODES(SN, OPTIONS)
% Build a symbolic description of the mean-field ODE system integrated by
% SolverFLD for the ODE-based methods. Two representations are produced,
% mirroring the numerical solver code paths:
%
%   form = 'W': dx/dt = W'*theta(x) + lambda   (methods: default, matrix, pnorm)
%       W        (nstates x nstates) phase-transition rate matrix
%       Alambda  (nstates x 1) exogenous arrival rate vector
%       theta_s  = x_s * min(n_i, S_i)/n_i  ('min' smoothing), or
%                = x_s * (1+(n_i/S_i)^p_i)^(-1/p_i)  ('pnorm' smoothing)
%       with theta_s = 0 for states of Source (EXT) stations, and n_i the
%       total mass at the station i owning state s.
%
%   form = 'J': dx/dt = J*r(x)                 (methods: closing, statedep, softmin)
%       J        (nstates x nevents) stoichiometry (jump) matrix
%       r_e(x)   = coeff(e) * factor_e(x), where factor_e depends on the
%                scheduling strategy of the station owning the state
%                variable that drives event e.
%
% Factor types (form 'J'):
%   'lin'      x_v
%   'min'      x_v * min(n_i, S_i)/n_i
%   'ext1'     1 - sum of x over phases 2..end of the class at the source
%   'dps'      x_v / (c0 + ntilde_i), scaled weights folded into coeff;
%              ntilde_i = sum_r w_ir n_ir (normalized DPS weights)
%   'dpspw'    piecewise: x_v if n_i <= S_i, else S_i*w_ir*x_v/ntilde_i
%   'fcfsw'    x_v * min(n_i, S_i)/nhat_i, phase weight folded into coeff;
%              nhat_i = sum_u w_u x_u with w_u = -1/D0(k,k) (mean phase time)
%   'fcfsws'   x_v * softmin(n_i, S_i; alpha)/nhat_i, phase weight in coeff
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
K = sn.nclasses;

method = options.method;
method = strrep(method, 'fluid.', '');
switch method
    case {'default','matrix','pnorm'}
        form = 'W';
        if strcmp(method,'default')
            method = 'matrix';
        end
    case {'closing','statedep','softmin'}
        form = 'J';
    otherwise
        line_error(mfilename, sprintf('Symbolic ODE export is unsupported for method ''%s''. Supported methods: default, matrix, pnorm, closing, statedep, softmin.', method));
end

sys = struct();
sys.form = form;
sys.method = method;
sys.stationNames = cell(M,1);
sys.classNames = cell(K,1);
sys.schedNames = cell(M,1);
for i = 1:M
    sys.stationNames{i} = sn.nodenames{sn.stationToNode(i)};
    sys.schedNames{i} = SchedStrategy.toText(sn.sched(i));
end
for r = 1:K
    sys.classNames{r} = sn.classnames{r};
end
sys.sched = sn.sched;

switch form
    case 'W'
        sys = symodes_wform(sys, sn, options, M, K);
    case 'J'
        sys = symodes_jform(sys, sn, M, K, method);
end
end

function sys = symodes_wform(sys, sn, options, M, K)
% Mirror the ODE construction of solver_fluid_matrix (Ruuskanen et al.,
% PEVA 151 (2021)): dx/dt = W'*theta(x) + Alambda.

pie = sn.pie;
PH = sn.proc;
NK = sn.njobs';
S = sn.nservers;
S(isinf(S)) = sum(NK);
nphases = sn.phases;

% station-to-station routing matrix via stochastic complementation
station_indices = [];
for ist = 1:M
    isf = sn.stationToStateful(ist);
    for r = 1:K
        station_indices = [station_indices, (isf-1)*K + r]; %#ok<AGROW>
    end
end
P = dtmc_stochcomp(sn.rt, station_indices);

% remove Sink->Source feedback routing for open classes
for src_ist = 1:M
    if sn.sched(src_ist) == SchedStrategy.EXT
        for r = 1:K
            if ~isnan(sn.rates(src_ist, r)) && sn.rates(src_ist, r) > 0
                src_col = (src_ist - 1) * K + r;
                for from_ist = 1:M
                    if from_ist ~= src_ist
                        for from_r = 1:K
                            P((from_ist - 1) * K + from_r, src_col) = 0;
                        end
                    end
                end
            end
        end
    end
end

% W = Psi + B*P*A'
Psi = [];
A = [];
B = [];
for ist = 1:M
    for r = 1:K
        if nphases(ist,r) == 0
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

% exogenous arrival rates into queue phases
source_arrivals = zeros(M, K);
for src_ist = 1:M
    if sn.sched(src_ist) == SchedStrategy.EXT
        for r = 1:K
            if ~isnan(sn.rates(src_ist, r)) && sn.rates(src_ist, r) > 0
                source_arrivals(src_ist, r) = sn.rates(src_ist, r);
            end
        end
    end
end
Alambda_full = zeros(size(A,1), 1);
state = 0;
for ist = 1:M
    for r = 1:K
        if nphases(ist,r) > 0
            if sn.sched(ist) == SchedStrategy.EXT
                state = state + nphases(ist,r);
            else
                arrival_rate_to_queue = 0;
                for src_ist = 1:M
                    if source_arrivals(src_ist, r) > 0
                        arrival_rate_to_queue = arrival_rate_to_queue + source_arrivals(src_ist, r) * P((src_ist - 1) * K + r, (ist - 1) * K + r);
                    end
                end
                if arrival_rate_to_queue > 0
                    for k = 1:nphases(ist,r)
                        state = state + 1;
                        Alambda_full(state) = pie{ist}{r}(k) * arrival_rate_to_queue;
                    end
                else
                    state = state + nphases(ist,r);
                end
            end
        else
            state = state + 1;
        end
    end
end

% state metadata over the same enumeration, then keep-filter
stateStation = zeros(size(A,1),1);
stateClass = zeros(size(A,1),1);
statePhase = zeros(size(A,1),1);
state = 0;
for ist = 1:M
    for r = 1:K
        if nphases(ist,r) == 0
            state = state + 1;
            stateStation(state) = ist;
            stateClass(state) = r;
            statePhase(state) = 0; % placeholder, removed by keep filter
        else
            for k = 1:nphases(ist,r)
                state = state + 1;
                stateStation(state) = ist;
                stateClass(state) = r;
                statePhase(state) = k;
            end
        end
    end
end

keep = find(~isnan(sum(W,1)));
W = W(keep,:);
W = W(:,keep);
sys.keep = keep;
sys.W = W;
sys.Alambda = Alambda_full(keep);
sys.stateStation = stateStation(keep);
sys.stateClass = stateClass(keep);
sys.statePhase = statePhase(keep);
sys.nstates = length(keep);
sys.S = S;
sys.isSource = (sn.sched(sys.stateStation) == SchedStrategy.EXT);

% smoothing selection mirrors use_pnorm in solver_fluid_matrix
use_pnorm = isfield(options, 'pstar') && ~isempty(options.pstar) || ...
    (isfield(options.config, 'pstar') && ~isempty(options.config.pstar));
if use_pnorm
    if isfield(options, 'pstar') && ~isempty(options.pstar)
        pstar_val = options.pstar;
    else
        pstar_val = options.config.pstar;
    end
    if isscalar(pstar_val)
        pstar_val = pstar_val * ones(M, 1);
    end
    sys.smoothing = 'pnorm';
    sys.pstar = pstar_val(:);
else
    sys.smoothing = 'min';
    sys.pstar = [];
end
end

function sys = symodes_jform(sys, sn, M, K, method)
% Mirror the event enumeration of ode_jumps_new/ode_rate_base and the rate
% scaling of ode_rates_closing/ode_statedep/ode_softmin: dx/dt = J*r(x).

Mu = sn.mu;
Phi = sn.phi;
PH = sn.proc;
sched = sn.sched;
rt = sn.rt;
S = sn.nservers;
N = sn.nclosedjobs;
for ist = 1:M
    for k = 1:K
        if isnan(Mu{ist}{k})
            Mu{ist}{k} = [];
            Phi{ist}{k} = [];
        end
    end
    if isinf(S(ist))
        S(ist) = N;
    end
end

if any(sched == SchedStrategy.EXT) && any(strcmp(method, {'statedep','softmin'}))
    line_error(mfilename, sprintf('The ''%s'' ODE method does not support open models, so their ODE system cannot be exported. Use the ''matrix'' or ''closing'' method instead.', method));
end

% state indexing as in solver_fluid_odes
q_indices = zeros(M,K);
Kic = zeros(M,K);
enabled = false(M,K);
cs = 1;
for i = 1:M
    for c = 1:K
        numphases = length(Mu{i}{c});
        q_indices(i,c) = cs;
        enabled(i,c) = numphases > 0;
        Kic(i,c) = numphases;
        cs = cs + numphases;
    end
end
nstates = cs - 1;

stateStation = zeros(nstates,1);
stateClass = zeros(nstates,1);
statePhase = zeros(nstates,1);
for i = 1:M
    for c = 1:K
        for k = 1:Kic(i,c)
            stateStation(q_indices(i,c)+k-1) = i;
            stateClass(q_indices(i,c)+k-1) = c;
            statePhase(q_indices(i,c)+k-1) = k;
        end
    end
end

% normalized DPS weights (as in ode_rates_closing/ode_statedep/ode_softmin)
dpsw = ones(M,K);
for i = 1:M
    if sched(i) == SchedStrategy.DPS
        dpsw(i,:) = sn.schedparam(i,:);
        dpsw(i,:) = dpsw(i,:)/sum(dpsw(i,:));
    end
end

% FCFS phase weights used by statedep/softmin: w = -1/D0(k,k)
fcfsPhaseW = zeros(nstates,1);
if any(strcmp(method, {'statedep','softmin'}))
    for i = 1:M
        if sched(i) == SchedStrategy.FCFS
            for c = 1:K
                if enabled(i,c)
                    for k = 1:Kic(i,c)
                        fcfsPhaseW(q_indices(i,c)+k-1) = -1/PH{i}{c}{1}(k,k);
                    end
                end
            end
        end
    end
end

% under statedep/softmin, stations with an unmatched scheduling strategy
% contribute no outgoing events (their switch has no otherwise branch)
handledSched = [SchedStrategy.INF, SchedStrategy.EXT, SchedStrategy.PS, ...
    SchedStrategy.FCFS, SchedStrategy.DPS];

% event enumeration (departures first, then phase changes) as in
% ode_jumps_new/ode_rate_base; zero-rate events are dropped
J = zeros(nstates,0);
coeff = [];
eventVar = [];
factorType = {};
factorData = {};
ne = 0;
for i = 1:M
    for c = 1:K
        if enabled(i,c)
            xic = q_indices(i,c);
            for j = 1:M
                for l = 1:K
                    if rt((i-1)*K+c,(j-1)*K+l) > 0
                        if isempty(PH{j}{l})
                            pievec = 1;
                        else
                            pievec = map_pie(PH{j}{l});
                        end
                        xjl = q_indices(j,l);
                        for ki = 1:Kic(i,c)
                            for kj = 1:Kic(j,l)
                                if any(strcmp(method, {'statedep','softmin'}))
                                    if sched(i) == SchedStrategy.INF && j == i
                                        continue % self-loop departures skipped at INF stations
                                    end
                                    if ~any(sched(i) == handledSched)
                                        continue
                                    end
                                end
                                base = Phi{i}{c}(ki) * Mu{i}{c}(ki) * rt((i-1)*K+c,(j-1)*K+l) * pievec(kj);
                                if base > 0
                                    ne = ne + 1;
                                    col = zeros(nstates,1);
                                    col(xic+ki-1) = col(xic+ki-1) - 1;
                                    col(xjl+kj-1) = col(xjl+kj-1) + 1;
                                    J(:,ne) = col;
                                    coeff(ne,1) = base; %#ok<AGROW>
                                    eventVar(ne,1) = xic+ki-1; %#ok<AGROW>
                                    [factorType{ne,1}, factorData{ne,1}, coeff(ne,1)] = ...
                                        symodes_factor(method, sched(i), i, c, ki, coeff(ne,1), ...
                                        q_indices, Kic, S, dpsw, fcfsPhaseW, xic+ki-1); %#ok<AGROW>
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
for i = 1:M
    for c = 1:K
        if enabled(i,c)
            if any(strcmp(method, {'statedep','softmin'})) && ~any(sched(i) == handledSched)
                continue
            end
            xic = q_indices(i,c);
            for ki = 1:(Kic(i,c) - 1)
                for kip = 1:Kic(i,c)
                    if ki ~= kip
                        base = PH{i}{c}{1}(ki,kip);
                        if base > 0
                            ne = ne + 1;
                            col = zeros(nstates,1);
                            col(xic+ki-1) = -1;
                            col(xic+kip-1) = 1;
                            J(:,ne) = col;
                            coeff(ne,1) = base;
                            eventVar(ne,1) = xic+ki-1;
                            [factorType{ne,1}, factorData{ne,1}, coeff(ne,1)] = ...
                                symodes_factor(method, sched(i), i, c, ki, coeff(ne,1), ...
                                q_indices, Kic, S, dpsw, fcfsPhaseW, xic+ki-1);
                        end
                    end
                end
            end
        end
    end
end

sys.J = J;
sys.coeff = coeff;
sys.eventVar = eventVar;
sys.factorType = factorType;
sys.factorData = factorData;
sys.nevents = ne;
sys.nstates = nstates;
sys.stateStation = stateStation;
sys.stateClass = stateClass;
sys.statePhase = statePhase;
sys.S = S;
sys.dpsw = dpsw;
sys.fcfsPhaseW = fcfsPhaseW;
if strcmp(method,'softmin')
    sys.alpha = 20; % softmin parameter, as in solver_fluid_odes
end
end

function [ftype, fdata, coeff] = symodes_factor(method, schedi, i, c, ki, coeff, q_indices, Kic, S, dpsw, fcfsPhaseW, var)
% Scaling factor applied to the driving state variable of an event whose
% source station is i (class c, phase ki), matching the per-strategy rate
% scaling in ode_rates_closing (closing), ode_statedep and ode_softmin.
fdata = struct('station', i, 'class', c);
switch method
    case 'closing'
        switch schedi
            case SchedStrategy.INF
                ftype = 'lin';
            case SchedStrategy.EXT
                if ki == 1
                    ftype = 'ext1';
                    fdata.others = (q_indices(i,c)+1):(q_indices(i,c)+Kic(i,c)-1);
                else
                    ftype = 'lin';
                end
            case {SchedStrategy.PS, SchedStrategy.FCFS}
                ftype = 'min';
            case SchedStrategy.DPS
                % ode_rates_closing uses denominator mean(w) + sum_r w_r*n_r
                ftype = 'dps';
                fdata.c0 = mean(dpsw(i,:));
                coeff = coeff * S(i) * dpsw(i,c);
            otherwise
                % strategies without a case in ode_rates_closing keep rates = x
                ftype = 'lin';
        end
    case {'statedep','softmin'}
        switch schedi
            case SchedStrategy.INF
                ftype = 'lin';
            case SchedStrategy.PS
                ftype = 'min';
            case SchedStrategy.FCFS
                if strcmp(method,'softmin')
                    ftype = 'fcfsws';
                else
                    ftype = 'fcfsw';
                end
                coeff = coeff * fcfsPhaseW(var);
            case SchedStrategy.DPS
                ftype = 'dpspw';
        end
end
end
