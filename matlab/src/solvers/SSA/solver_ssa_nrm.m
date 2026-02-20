function [QN, UN, RN, TN, CN, XN, lG, sn] = solver_ssa_nrm(sn, options)
% SOLVER_SSA_NRM   Steady‑state analysis via the Next‑Reaction Method (SSA)
%
%   [QN, UN, RN, TN, CN, XN, LG, SN] = SOLVER_SSA_NRM(SN, OPTIONS)
%   runs a stochastic simulation of the queueing network described in SN
%   for OPTIONS.samples reaction firings using Gibson & Bruck's
%   Next‑Reaction Method.  During the run it:
%     • computes performance metrics directly during simulation;
%     • returns standard queueing performance measures.
%
%   Outputs
%     QN        – M×K matrix of mean queue lengths
%     UN        – M×K matrix of utilizations
%     RN        – M×K matrix of response times
%     TN        – M×K matrix of throughputs
%     CN        – 1×K vector of cycle times
%     XN        – 1×K vector of system throughputs
%     LG        – Logarithm of normalizing constant (not computed)
%     SN        – (Possibly updated) network structure.
%
%   See also SOLVER_SSA_NRM_SPACE, NEXT_REACTION_METHOD_DIRECT.

% ---------------------------------------------------------------------
% Parameters & shorthands
% ---------------------------------------------------------------------
samples = options.samples;
R  = sn.nclasses;
I  = sn.nnodes;
M  = sn.nstations;
K  = sn.nclasses;
state = sn.state;

% ---------------------------------------------------------------------
% Stoichiometry & reaction mapping (self‑loops included) ----------------
% ---------------------------------------------------------------------
S       = zeros(0, I*R);   % will transpose at the end
fromIdx = [];
toIdx   = [];
fromIR  = [];

% this is currently M^2*R^2, it can be lowered to M*R decoupling the
% routing
k = 0;
for ind = 1:I
    for r = 1:R
        k = k + 1;
        fromIR(k,:) = [ind, r];
        fromIdx(k) = (ind-1)*R + r;
        probIR{k} = [];
        toIdx{k} = [];
        Srow = zeros(1, I*R); % build stoichiometry row
        if sn.isslc(r)
            Srow(fromIdx(k))   = -Inf;
        else
            Srow(fromIdx(k)) = -1;
            for jnd = 1:I
                for s = 1:R
                    if sn.rtnodes((ind-1)*R+r, (jnd-1)*R+s) > 0
                        toIdx{k}(end+1) = (jnd-1)*R + s;
                        p = sn.rtnodes((ind-1)*R+r, (jnd-1)*R+s);
                        probIR{k}(end+1) = p;
                        Srow((jnd-1)*R + s)   = Srow((jnd-1)*R + s) + p;
                    end
                end
            end
        end
        S(k,:) = Srow;
    end
end
S = S.';   % states × reactions

% ---------------------------------------------------------------------
% Initial state vector --------------------------------------------------
% ---------------------------------------------------------------------
nvec0 = zeros(I*R,1); % initial state (aggregate state)
for ind=1:I
    if sn.isstateful(ind)
        [~,nir] = State.toMarginalAggr(sn, ind, state{sn.nodeToStateful(ind)});
        for r = 1:R
            if isinf(nir(r))
                if sn.nodetype(ind) == NodeType.Source
                    nir(r) = 1;
                else
                    line_error(mfilename, 'Infinite population error.');
                end
            end
            nvec0((ind-1)*R + r,1) = nir(r);
        end
    end
end

mi    = zeros(I,1);
rates = zeros(I,R);
for ind=1:I
    if sn.isstation(ind)
        for r=1:R
            ist = sn.nodeToStation(ind);
            muir = sn.rates(ist,r);
            if ~isnan(muir)
                rates(ind,r) = muir;
            end
            mi(ind,1) = sn.nservers(ist);
        end
    else
        for r=1:R
            rates(ind,r) = GlobalConstants.Immediate;
            mi(ind,1) = GlobalConstants.MaxInt;
        end
    end
    mi(isinf(mi)) = GlobalConstants.MaxInt;
end

% Propensity function ---------------------------------------------------
epstol = GlobalConstants.Zero;
a = {};
for j=1:length(fromIdx)
    if sn.isstation(fromIR(j,1))
        switch sn.sched(sn.nodeToStation(fromIR(j,1)))
            case SchedStrategy.EXT
                a{j} = @(X) rates(fromIR(j,1), fromIR(j,2));
            case SchedStrategy.INF
                a{j} = @(X) rates(fromIR(j,1), fromIR(j,2)) * X(fromIdx(j));
            case SchedStrategy.PS
                if R == 1 % single class
                    a{j} = @(X) rates(fromIR(j,1), fromIR(j,2)) * min( mi(fromIR(j,1)), X(fromIdx(j)));
                else
                    a{j} = @(X) rates(fromIR(j,1), fromIR(j,2)) * ( X(fromIdx(j)) ./ ...
                        (epstol+sum( X(((fromIR(j,1)-1)*R + 1):((fromIR(j,1)-1)*R + R)) ) )) * ...
                        min( mi(fromIR(j,1)), ...
                        (epstol+sum( X(((fromIR(j,1)-1)*R + 1):((fromIR(j,1)-1)*R + R)) ) ));
                end
        end
    else
        a{j} = @(X) rates(fromIR(j,1), fromIR(j,2)) * min(1, X(fromIdx(j)));
    end
end

% Propensity functions dependencies -----------------------------------
D = cell(1,size(S,2));
for k=1:size(D,2)
    J = find(S(:,k))'; % set of state variables affected by reaction k
    vecd = [];
    for j=1:length(J)
        % (ind-1)*R + r
        pos = J(j);
        r   = mod(pos-1, R) + 1;
        ind = ((pos-r)/R) + 1;
        vecd(end+1:end+R) = ((ind-1)*R + 1) : (ind*R);
    end
    % vecd now contains all state variables affected by the firing of
    % reaction k. We now find the propensity functions that depend
    % on those variables
    if ~isempty(vecd)
        vecd = unique(vecd);
        vecs = [];
        for j=1:length(vecd)
            vecs = [vecs,find(S(vecd(j),:)<0)];
        end
        D{k} = unique(vecs);
    else
        D{k} = [];
    end
end

% Having accounted for them in D, we can now remove self-loops markings
S(isinf(S))=0;

% ---------------------------------------------------------------------
% Initialize performance metric matrices
% ---------------------------------------------------------------------
lG = 0; % Not computed in SSA

% ---------------------------------------------------------------------
% Run SSA/NRM with direct metric computation
% ---------------------------------------------------------------------
[QN, UN, RN, TN, CN, XN] = next_reaction_method_direct(S, D, a, nvec0, samples, options, sn, fromIdx, fromIR);

end  % solver_ssa_nrm

% ======================================================================
% Next-Reaction Method with direct metric computation
% ======================================================================
function [QN, UN, RN, TN, CN, XN] = next_reaction_method_direct(S, D, a, nvec0, samples, options, sn, fromIdx, fromIR)

numReactions = size(S,2);
R = sn.nclasses;
I = sn.nnodes;
M = sn.nstations;
K = sn.nclasses;
QN = zeros(M, K);
UN = zeros(M, K);
RN = zeros(M, K);
TN = zeros(M, K);
CN = zeros(1, K);
XN = zeros(1, K);

% when a reaction fires, this matrix helps selecting the probability that a
% particular routing or phase is selected as a result ------------------
P = S; P(P<0)=P(P<0)+1';
fromIdxCell = cell(numReactions,1);
toIdxCell = cell(numReactions,1);
cdfVec = cell(numReactions,1);
for r=1:numReactions
    nnzP(r) = nnz(P(:,r));
    if nnzP(r)>1
        fromIdxCell{r} = find(S(:,r)<0);
        toIdxCell{r} = find(P(:,r));
        cdfVec{r} = cumsum(P(toIdxCell{r},r));
    end
end

% initialise Gillespie clocks ------------------------------------------
t   = 0;
for k=1:size(S,2)
    Ak(k) = a{k}(nvec0);
end
nvec   = nvec0;
Pk  = -log(rand(1,numReactions));
Tk  = zeros(1,numReactions);

tau = (Pk - Tk) ./ Ak;

% Performance tracking variables
totalTime = 0;
NK = sn.njobs'; % Jobs per class
servers = sn.nservers;

n = 1;
while n <= samples
    [dt, kfire] = min(tau);
    if isinf(dt), line_error(mfilename,'Deadlock. Quitting nrm method.'); end

    totalTime = totalTime + dt;

    % Accumulate state-dependent metrics during this time interval
    for ist = 1:M
        ind = sn.stationToNode(ist);
        for k = 1:K
            currentPop = nvec((ind-1)*R + k);

            % Accumulate queue length (QN)
            QN(ist, k) = QN(ist, k) + currentPop * dt;

            % Compute throughput contribution from departures
            irIndex = (ind-1)*R + k;
            depRate = 0;
            idx = find(fromIdx == irIndex);
            if ~isempty(idx)
                depRate = Ak(idx(1));
            end
            TN(ist, k) = TN(ist, k) + depRate * dt;

            % Compute utilization based on scheduling policy
            switch sn.sched(ist)
                case {SchedStrategy.INF, SchedStrategy.EXT}
                    UN(ist, k) = UN(ist, k) + currentPop * dt;
                case SchedStrategy.PS
                    totalPop = sum(nvec(((ind-1)*R + 1):(ind*R)));
                    if totalPop > 0
                        utilization = (currentPop / totalPop) * min(servers(ist), totalPop) / servers(ist);
                    else
                        utilization = 0;
                    end
                    UN(ist, k) = UN(ist, k) + utilization * dt;
            end
        end
    end

    t  = t + dt;

    % update aggregate state
    if nnzP(kfire)>1
        r = 1+find(rand>=cdfVec{kfire},1);
        if isempty(r)
            r = 1;
        end
        nvec(fromIdxCell{kfire}) = nvec(fromIdxCell{kfire}) - 1;
        nvec(toIdxCell{kfire}(r)) = nvec(toIdxCell{kfire}(r)) + 1;
    else
        nvec  = nvec + S(:,kfire);  % zero change for self-loops
    end

    Tk = Tk + Ak * dt;

    % update rates for all reactions dependent on the last fired reaction
    for k=D{kfire}
        Ak(k) = a{k}(nvec);
    end

    % update clocks
    Pk(kfire)  = Pk(kfire) - log(rand);
    tau        = (Pk - Tk) ./ Ak;
    tau(Ak==0) = inf;

    % do not count immediate events
    n = n + 1;
    print_progress(options, n);
end % while
% Print newline after progress counter
if isfield(options,'verbose') && options.verbose
    line_printf('\n');
end

% Normalize metrics by total time
if totalTime > 0
    for ist = 1:M
        for k = 1:K
            QN(ist, k) = QN(ist, k) / totalTime;
            UN(ist, k) = UN(ist, k) / totalTime;
            TN(ist, k) = TN(ist, k) / totalTime;
        end
    end
end

% Compute derived metrics
for k = 1:K
    % System throughput at reference station
    XN(1, k) = TN(sn.refstat(k), k);

    % Response times
    for ist = 1:M
        if TN(ist, k) > 0
            RN(ist, k) = QN(ist, k) / TN(ist, k);
        else
            RN(ist, k) = 0;
        end
    end

    % Cycle times
    if XN(1, k) > 0
        CN(1, k) = NK(k) / XN(1, k);
    end
end

% Handle NaN values
QN(isnan(QN)) = 0;
UN(isnan(UN)) = 0;
RN(isnan(RN)) = 0;
XN(isnan(XN)) = 0;
TN(isnan(TN)) = 0;
CN(isnan(CN)) = 0;

    function print_progress(opt, samples_collected)
        if ~isfield(opt,'verbose') || ~opt.verbose, return; end
        if samples_collected == 1e3
            line_printf('\bSSA samples: %8d', samples_collected);
        elseif opt.verbose == 2
            if samples_collected == 0
                line_printf('\bSSA samples: %9d', samples_collected);
            else
                line_printf('\b\b\b\b\b\b\b\b\b%9d', samples_collected);
            end
        elseif mod(samples_collected,1e3)==0 || opt.verbose == 2
            line_printf('\b\b\b\b\b\b\b\b\b%9d', samples_collected);
        end
    end
end  % next_reaction_method_direct