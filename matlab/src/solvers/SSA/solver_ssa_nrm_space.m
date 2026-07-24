function [pi, outspace, depRates, sn] = solver_ssa_nrm_space(sn, options)
% SOLVER_SSA_NRM_SPACE   Steady‑state analysis via the Next‑Reaction Method (SSA)
%
%   [PI, SSQ, ARVRATES, DEPRATES, SN] = SOLVER_SSA_NRM_SPACE(SN, OPTIONS)
%   runs a stochastic simulation of the queueing network described in SN
%   for OPTIONS.samples reaction firings using Gibson & Bruck's
%   Next‑Reaction Method.  During the run it:
%     • observes every distinct global state visited and the time spent in it;
%     • accumulates the per‑state propensities of **all** enabled reactions,
%       including *self‑loops* (service completions routed back to the same
%       queue).
%
%   Outputs
%     PI        – 1×S vector of empirical steady‑state probabilities
%                 (sojourn‑time fractions) for the S unique states;
%     SSQ       – S×(M·R) matrix listing those states row‑by‑row in the
%                 flattened (station, class) order;
%     DEPRATES  – S×(M·R) matrix of total departure rates from each queue
%                 in the corresponding state;
%     SN        – (Possibly updated) network structure.
%
%   See also NEXT_REACTION_METHOD.

% ---------------------------------------------------------------------
% Parameters & shorthands
% ---------------------------------------------------------------------
samples = options.samples;
R  = sn.nclasses;
I  = sn.nnodes;
state = sn.state;
% ---------------------------------------------------------------------
% Stoichiometry & reaction mapping (self‑loops included) ----------------
% ---------------------------------------------------------------------
S       = zeros(0, I*R);   % will transpose at the end
fromIdx = [];
toIdx   = [];
fromIR  = [];

% see _kb/06-solver-catalog.md for rationale (SSA NRM reaction mapping O(M^2*R^2))
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
bufferedSched = [SchedStrategy.FCFS, SchedStrategy.LCFS];
buffers0 = cell(I,1); % per-node ordered buffer of waiting job classes (FCFS/LCFS)
for ind=1:I
    buffers0{ind} = [];
end
for ind=1:I
    if sn.isstateful(ind)
        state_i = state{sn.nodeToStateful(ind)};
        [~,nir] = State.toMarginalAggr(sn, ind, state_i);
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

        % Populate buffers for buffered nodes from the raw state vector
        % (only stations have FCFS/LCFS scheduling; skip non-station
        % stateful nodes such as RROBIN dispatchers/Routers and Caches)
        ist = sn.nodeToStation(ind);
        if ist >= 1 && any(sn.sched(ist) == bufferedSched)
            sumK = sum(sn.phasessz(ist,:));
            sumNvars = sum(sn.nvars(ind,:));
            bufCols = size(state_i,2) - sumK - sumNvars;
            for pos = 1:bufCols
                classId = state_i(1,pos);
                if classId >= 1 && classId <= R
                    buffers0{ind}(end+1) = classId; % addLast
                end
                % classId == 0 means empty position, skip
            end
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
                a{j} = @(X, bufs) rates(fromIR(j,1), fromIR(j,2));
            case SchedStrategy.INF
                a{j} = @(X, bufs) rates(fromIR(j,1), fromIR(j,2)) * X(fromIdx(j));
            case SchedStrategy.PS
                if R == 1 % single class
                    a{j} = @(X, bufs) rates(fromIR(j,1), fromIR(j,2)) * min( mi(fromIR(j,1)), X(fromIdx(j)));
                else
                    a{j} = @(X, bufs) rates(fromIR(j,1), fromIR(j,2)) * ( X(fromIdx(j)) ./ ...
                        (epstol+sum( X(((fromIR(j,1)-1)*R + 1):((fromIR(j,1)-1)*R + R)) ) )) * ...
                        min( mi(fromIR(j,1)), ...
                        (epstol+sum( X(((fromIR(j,1)-1)*R + 1):((fromIR(j,1)-1)*R + R)) ) ));
                end
            case {SchedStrategy.FCFS, SchedStrategy.LCFS}
                % Rate proportional to the jobs actually being served, i.e. the
                % class-r population minus the class-r jobs waiting in buffer.
                a{j} = @(X, bufs) rates(fromIR(j,1), fromIR(j,2)) * ...
                    max(0, X(fromIdx(j)) - sum(bufs{fromIR(j,1)} == fromIR(j,2)));
        end
    else
        a{j} = @(X, bufs) rates(fromIR(j,1), fromIR(j,2)) * min(1, X(fromIdx(j)));
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
% Run SSA/NRM -----------------------------------------------------------
% ---------------------------------------------------------------------
if false %snIsClosedModel(sn)
    % mixed-radix hashing
    reactCache = containers.Map('KeyType','uint64','ValueType','any');
    njobs = sn.njobs;
    mixedradix = [cumprod(repmat(1+njobs,1,I))];
    mixedradix = [1,mixedradix(1:end-1)];
    hashfun = @(v, bufs) uint64(mixedradix*v(:));
else
    % buffer size unbounded so use string; the key combines the aggregate
    % state vector with the ordered contents of every node buffer, so that
    % FCFS/LCFS states differing only in queueing order remain distinct.
    reactCache = containers.Map('KeyType','char','ValueType','any');
    hashfun = @(v, bufs) [mat2str(v(:)'), '|', bufferHashAll(bufs)];
end
[t, nvecsim, bufferStates, ~, ~] = next_reaction_method(S, D, a, nvec0, buffers0, samples, options, reactCache, hashfun, fromIR, mi, R, sn);

% ---------------------------------------------------------------------
% Empirical state probabilities ----------------------------------------
% ---------------------------------------------------------------------
dt  = diff(t);
% Find unique states keyed on both the aggregate state and the buffers
numIntervals = size(nvecsim,2) - 1;
stateKeys = cell(numIntervals,1);
for i = 1:numIntervals
    stateKeys{i} = hashfun(nvecsim(:,i), bufferStates{i});
end
[~, ia, ic] = unique(stateKeys, 'stable');
outspace        = nvecsim(:, ia).';
outspaceBuffers = bufferStates(ia);
timeAccum    = accumarray(ic, dt(:));
pi           = timeAccum / sum(timeAccum);

% ---------------------------------------------------------------------
% Per‑state arrival / departure rates (self‑loops counted) -------------
% ---------------------------------------------------------------------
numStates = size(outspace,1);
depRates  = zeros(numStates, I*R);

for st = 1:numStates
    a_state = reactCache(hashfun(outspace(st,:)', outspaceBuffers{st}));
    for j = 1:length(fromIdx)
        depRates(st, fromIdx(j)) = depRates(st, fromIdx(j)) + a_state(j);
    end
end
end  % solver_ssa_nrm_space

% ======================================================================
% Next-Reaction Method core --------------------------------------------
% ======================================================================
function [t, nvec, bufferStates, kfires, rfires] = next_reaction_method(S, D, a, nvec0, buffers0, samples, options, reactcache, hashfun, fromIR, mi, R, sn)
numReactions = size(S,2);
rand_pool_size = 1e7;
buffers = buffers0; % working copy of the per-node ordered buffers

% when a reaction fires, this matrix helps selecting the probability that a
% particular routing or phase is selected as a result ------------------
P = S; P(P<0)=P(P<0)+1';
fromIdx = cell(numReactions,1);
toIdx = cell(numReactions,1);
cdfVec = cell(numReactions,1);
for r=1:numReactions
    nnzP(r) = nnz(P(:,r));
    if nnzP(r)>1
        fromIdx{r} = find(S(:,r)<0);
        toIdx{r} = find(P(:,r));
        cdfVec{r} = cumsum(P(toIdx{r},r));
    end
end

% initialise Gillespie clocks ------------------------------------------
t   = 0;
for k=1:size(S,2)
    Ak(k) = a{k}(nvec0, buffers);
end
nvec   = nvec0;
key = hashfun(nvec, buffers);
reactcache(key) = Ak;         % cache first state's propensities
Pk  = -log(rand(1,numReactions));
Tk  = zeros(1,numReactions);

tau = (Pk - Tk) ./ Ak;

% logs -----------------------------------------------------------------
tout = zeros(samples,1);
nvecout = zeros(length(nvec0),samples);
bufferStates = cell(1, samples+1);
bufferStates{1} = buffers; % buffers in the initial state
kfires = zeros(samples,1);
rfires = zeros(samples,1);
n = 1;
while n <= samples
    [dt, kfire] = min(tau);
    kfires(n) = kfire;
    if isinf(dt), line_error(mfilename,'Deadlock. Quitting nrm method.'); end

    t  = t + dt;

    % update aggregate state
    destPos = [];
    if nnzP(kfire)>1
        r = 1+find(rand>=cdfVec{kfire},1);
        if isempty(r)
            r = 1;
        end
        rfires(n) = r;
        nvec(fromIdx{kfire}) = nvec(fromIdx{kfire}) - 1;
        nvec(toIdx{kfire}(r)) = nvec(toIdx{kfire}(r)) + 1;
        destPos = toIdx{kfire}(r);
    else
        nvec  = nvec + S(:,kfire);  % zero change for self-loops
        dpos = find(S(:,kfire) > 0); % deterministic destination (single move)
        if ~isempty(dpos)
            destPos = dpos(1);
        end
    end

    % maintain FCFS/LCFS buffers given the source/destination of this firing
    buffers = updateBuffers(kfire, nvec, buffers, fromIR, destPos, mi, R, sn);

    Tk = Tk + Ak * dt;

    % update rates for all reactions dependent on the last fired reaction
    for k=D{kfire}
        Ak(k) = a{k}(nvec, buffers);
    end

    key = hashfun(nvec, buffers);
    if ~isKey(reactcache,key)
        reactcache(key) = Ak; % store propensities of new state
    end

    % maintain random number pool
    n_mod = mod(n,rand_pool_size);
    if n_mod == 1
        rand_pool = rand(1+min(rand_pool_size, samples-n),1);
    end

    % update clocks
    Pk(kfire)  = Pk(kfire) - log(rand_pool(n_mod));
    tau        = (Pk - Tk) ./ Ak;
    tau(Ak==0) = inf;

    % update measures
    tout(n) = t;
    nvecout(:,n) = nvec;
    bufferStates{n+1} = buffers;

    % do not count immediate events
    n = n + 1;
    print_progress(options, n);
end
% Print newline after progress counter
if isfield(options,'verbose') && options.verbose
    line_printf('\n');
end

t = [0; tout];
nvec = [nvec0, nvecout];

    function print_progress(opt, samples_collected)
        if ~isfield(opt,'verbose') || ~opt.verbose || batchStartupOptionUsed, return; end
        if samples_collected == 1e3
            line_printf('\nSSA samples: %8d', samples_collected);
        elseif opt.verbose == 2
            if samples_collected == 0
                line_printf('\nSSA samples: %9d', samples_collected);
            else
                line_printf('\b\b\b\b\b\b\b\b\b%9d', samples_collected);
            end
        elseif mod(samples_collected,1e3)==0 || opt.verbose == 2
            line_printf('\b\b\b\b\b\b\b\b\b%9d', samples_collected);
        end
    end
end  % next_reaction_method

% ======================================================================
% Buffer maintenance and hashing helpers for FCFS/LCFS nodes
% ======================================================================
function buffers = updateBuffers(kfire, nvec, buffers, fromIR, destPos, mi, R, sn)
% Maintain the ordered per-node buffers when reaction KFIRE fires. A
% departure removes the head/tail job of the source buffer; an arrival at a
% buffered destination whose servers are all busy joins the buffer head.
ind = fromIR(kfire,1); % source node of the firing

% Handle departure from FCFS/LCFS source node
if isFCFS(ind, sn)
    if ~isempty(buffers{ind})
        buffers{ind}(end) = []; % pollLast
    end
elseif isLCFS(ind, sn)
    if ~isempty(buffers{ind})
        buffers{ind}(1) = []; % pollFirst
    end
end

% Handle arrival at a buffered destination node
if ~isempty(destPos) && destPos > 0
    jnd = floor((destPos-1)/R) + 1;
    s   = mod(destPos-1, R) + 1;
    if isFCFS(jnd, sn) || isLCFS(jnd, sn)
        totalAtDest = sum(nvec(((jnd-1)*R + 1):(jnd*R)));
        if totalAtDest > mi(jnd)
            % All servers busy - arriving job joins back of buffer
            buffers{jnd} = [s, buffers{jnd}]; % addFirst
        end
        % Otherwise job went straight into service, buffer unchanged
    end
end
end

function s = bufferHashAll(bufs)
% Buffer hash over all node indices (used to distinguish unique states).
parts = cell(1, numel(bufs));
for ind = 1:numel(bufs)
    parts{ind} = [num2str(ind), ':', mat2str(bufs{ind}(:)')];
end
s = strjoin(parts, '|');
end

function tf = isFCFS(ind, sn)
tf = false;
if sn.isstation(ind)
    ist = sn.nodeToStation(ind);
    tf = (sn.sched(ist) == SchedStrategy.FCFS);
end
end

function tf = isLCFS(ind, sn)
tf = false;
if sn.isstation(ind)
    ist = sn.nodeToStation(ind);
    tf = (sn.sched(ist) == SchedStrategy.LCFS);
end
end
