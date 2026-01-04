function [pi, SSq, depRates, sn] = solver_ssa_nrm_OM2R2(sn, options)
% SOLVER_SSA_NRM   Steady‑state analysis via Gibson & Bruck's Next‑Reaction 
% Method (SSA) as modified in Anderson's, THE JOURNAL OF CHEMICAL PHYSICS
% 127, 214107, 2007.

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
rateIR  = [];

% This is currently M^2*R^2, it may be lowered to M*R decoupling the
% routing from the departure rate, but then the caching for the rates in 
% each state becomes less efficient as a given state in M*R can yield 
% multiple outgoing rates  depending on the routing target. Also, this 
% would still require a M^2*R^2 loop at the end to compute arvRates.
k = 0;
for ind = 1:I
    for r = 1:R
        for jnd = 1:I
            for s = 1:R
                if sn.rtnodes((ind-1)*R+r, (jnd-1)*R+s) == 0, continue; end
                k = k + 1;
                fromIdx(k) = (ind-1)*R + r;
                toIdx(k)   = (jnd-1)*R + s;
                rateIR(k,:) = [ind, r];
                probIR(k) = sn.rtnodes((ind-1)*R+r, (jnd-1)*R+s);

                % build stoichiometry row
                Srow = zeros(1, I*R);
                if fromIdx(k) ~= toIdx(k)
                    Srow(fromIdx(k)) = -1;
                    Srow(toIdx(k))   =  1;
                else % mark self-loops
                    Srow(fromIdx(k)) = -Inf;
                end
                S(k,:) = Srow;
            end
        end
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
                    nir(r) = GlobalConstants.MaxInt;
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
            rir = sn.rates(ist,r);
            if ~isnan(rir)
                rates(ind,r) = rir;
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
for k=1:length(fromIdx)
    if sn.isstation(rateIR(k,1))
        switch sn.sched(sn.nodeToStation(rateIR(k,1)))
            case SchedStrategy.EXT
                a{k} = @(X) rates(rateIR(k,1), rateIR(k,2)) * ...
                    probIR(k);
            case SchedStrategy.INF
                a{k} = @(X) rates(rateIR(k,1), rateIR(k,2)) * ...
                    probIR(k) * X(fromIdx(k));
            case SchedStrategy.PS
                if R == 1 % single class
                    a{k} = @(X) rates(rateIR(k,1), rateIR(k,2)) * ...
                        probIR(k) * ...
                        min( mi(rateIR(k,1)), X(fromIdx(k)));
                else
                    a{k} = @(X) rates(rateIR(k,1), rateIR(k,2)) * ...
                        probIR(k) * ...
                        ( X(fromIdx(k)) ./ ...
                        (epstol+sum( X(((rateIR(k,1)-1)*R + 1):((rateIR(k,1)-1)*R + R)) ) )) * ...
                        min( mi(rateIR(k,1)), ...
                        (epstol+sum( X(((rateIR(k,1)-1)*R + 1):((rateIR(k,1)-1)*R + R)) ) ));
                end
        end
    else
        a{k} = @(X) rates(rateIR(k,1), rateIR(k,2)) * ...
            probIR(k) * min(1, X(fromIdx(k)));
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
            vecs = [vecs,find(S(vecd(j),:)<=-1)];
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
    hashfun = @(v) uint64(mixedradix*v(:));
else
    % buffer size unbounded so use string
    reactCache = containers.Map('KeyType','char','ValueType','any');
    hashfun = @(v) mat2str(v(:)');
end
[t, X] = next_reaction_method(S, D, a, nvec0, samples, options, reactCache, hashfun);

% ---------------------------------------------------------------------
% Empirical state probabilities ----------------------------------------
% ---------------------------------------------------------------------
dt  = diff(t);
[SSq, ~, ic] = unique(X(:,1:end-1).', 'rows', 'stable');
timeAccum    = accumarray(ic, dt);
pi           = timeAccum / sum(timeAccum);

% ---------------------------------------------------------------------
% Per‑state arrival / departure rates (self‑loops counted) -------------
% ---------------------------------------------------------------------
numStates = size(SSq,1);
arvRates  = zeros(numStates, I*R);
depRates  = zeros(numStates, I*R);

for k = 1:numStates
    a_state = reactCache(hashfun(SSq(k,:)'));

    for ia = 1:length(fromIdx)
        depRates(k, fromIdx(ia)) = depRates(k, fromIdx(ia)) + a_state(ia);
        %arvRates(k, toIdx(ia))  = arvRates(k, toIdx(ia))  + a_state(ia);
    end
end
end  % solver_ssa_nrm

% ======================================================================
% Next-Reaction Method core --------------------------------------------
% ======================================================================
function [t, nvecout, kfire] = next_reaction_method(S, D, a, nvec0, samples, options, reactcache, hashfun)
numReactions = size(S,2);
rand_pool_size = 1e7;

% initialise Gillespie clocks -----------------------------------------
t  = 0;
Ak = zeros(1,size(S,2));
for k=1:size(S,2)
    Ak(k) = a{k}(nvec0);
end
nvec   = nvec0;
key = hashfun(nvec);
reactcache(key) = Ak;         % cache first state's propensities
Pk  = -log(rand(1,numReactions));
Tk  = zeros(1,numReactions);

tau = (Pk - Tk) ./ Ak;

% logs -----------------------------------------------------------------
tout = zeros(samples,1);
nvecsim = zeros(length(nvec0),samples);
kfires = zeros(samples,1);

n = 1;
while n <= samples
    [dt, kfire] = min(tau);
    kfires(n) = kfire; % record reaction that fired
    if isinf(dt), line_error(mfilename,'Deadlock. Quitting nrm method.'); end

    t  = t + dt;
    % update aggregate state
    nvec  = nvec + S(:,kfire);  % zero change for self-loops
    Tk = Tk + Ak * dt;

    % update rates for all reactions dependent on the last fired one
    for k=D{kfire}
        Ak(k) = a{k}(nvec);
    end

    key = hashfun(nvec);
    if ~isKey(reactcache,key)
        reactcache(key) = Ak; % store rates of new state
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
    nvecsim(:,n) = nvec;

    % do not count immediate events
    n = n + 1;
    print_progress(options, n);
end
% Print newline after progress counter
if isfield(options,'verbose') && options.verbose
    line_printf('\n');
end

t = [0; tout];
nvecout = [nvec0, nvecsim];

    function print_progress(opt, samples_collected)
        if ~isfield(opt,'verbose') || ~opt.verbose, return; end
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
