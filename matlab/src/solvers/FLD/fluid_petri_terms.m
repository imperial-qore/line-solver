function terms = fluid_petri_terms(sn, options)
% TERMS = FLUID_PETRI_TERMS(SN, OPTIONS)
%
% Event-based representation of the fluid marking process of a stochastic
% Petri net, in the contract SOLVER_FLUID_PETRI solves and FLUID_MOMENT_TERMS
% supplies for a queueing network.
%
% A GSPN IS ALREADY A DENSITY-DEPENDENT MARKOV POPULATION PROCESS, which is the
% object the moment-closure family of SolverFLD is built on: the marking is the
% population, a transition mode is a reaction, its incidence column is the jump,
% and the rate law lambda*min(enabling degree, servers) is the same min()
% non-linearity the min-normal closure exists to smooth. Nothing about the
% closure changes here; only where the drift comes from.
%
%   dx/dt = D * r(x, Sigma, phi, mu)
%
% THE STATE, x = [ m ; z ].
%
%   m(p,k)    token mass of class k at place p. One coordinate per (place,
%             class) pair that some arc touches, that the initial marking
%             loads, or that a Source feeds; a pair nothing reaches is dropped
%             rather than carried as a null direction of the Newton system.
%   y(j,h)    the number of mode-j servers running in phase h, for a mode whose
%             firing time has more than one phase. Their SUM is not free: the
%             ENABLE synchronization latches it instantaneously to
%             min(enabling degree, servers), so the latch is an ALGEBRAIC row
%             with one free-sign unknown mu_j -- the net rate at which servers
%             start or stop -- and the phase split evolves differentially.
%
%             CARRYING THE DISTRIBUTION INSTEAD OF THE COUNT LOOKS TIDIER AND IS
%             WRONG. Writing y = theta(m)*z and integrating z removes mu, but
%             the resulting z equation is missing the (dtheta/dt)/theta*(pie-z)
%             term and carries a spurious factor theta, so it agrees with the
%             count form only AT a fixed point. It also hides the coupling from
%             the marking into the firing rate, which the covariance needs: on a
%             closed cycle with Erlang(2) firing the count form reads 1.4564
%             against the exact 1.4531, the distribution form 1.4357.
%
% There is no Source coordinate: an exogenous arrival is a CONSTANT-propensity
% event depositing one token, as it is in the NRM SPN runner, which is the
% correct exogenous Poisson event once the source pool is projected away (see
% FLUID_MOMENT_TERMS for the same projection in the queueing case). There is no
% Sink coordinate either: a firing arc into a sink is mass leaving the net.
%
% THE EVENTS, one column of D each:
%
%   kind 1  a firing of mode j out of phase h into phase h'. The marking moves
%           by the incidence column post-pre; for a multi-phase mode the phase
%           moves by D1(h,h'), which is the PH restart alpha and the MAP's
%           landing phase in the same expression, since a PH's process pair has
%           D1 = s*alpha.
%   kind 2  an internal phase change of mode j, D0(h,h'): the marking is
%           unchanged, so this column contributes nothing to the diffusion of
%           the marking.
%   kind 3  an exogenous arrival, constant rate lambda*p, +1 token.
%   kind 4  a firing of an IMMEDIATE mode, at the algebraic flow phi_j.
%
% THE RATES:
%
%   theta_j(x) = g_inh_j(x) * E[ min_a( m_a/w_a ), c_j ]     enabled bindings
%   kind 1     = theta_j * z(j,h) * D1_j(h,h') * g_dep_j(m)
%   kind 2     = theta_j * z(j,h) * D0_j(h,h')
%   kind 3     = lambda * p
%   kind 4     = phi_j
%
% with the many-argument min closed by FLUID_MINMULTI_CLOSURE and the inhibitor
% indicator 1{m_b < thr} closed by the normal CDF Phi((thr-m_b)/sigma_b), which
% collapses to the hard indicator at sigma = 0 exactly as FLUID_MIN_CLOSURE
% does, so the first-order limit needs no second code path. For a single-phase
% mode z is absent and theta_j is the running-server count itself.
%
% THE COVARIANCE LIVES ON THE MARKING ONLY. z is a distribution, so its
% fluctuations are O(1) rather than O(sqrt(N)) and the linear noise
% approximation does not apply to it at any scale -- the same reason
% FLUID_MOMENT_TERMS projects the EXT source pool away. TERMS.covIdx therefore
% holds the marking coordinates and nothing else.
%
% Parameters:
%   sn      - NetworkStruct holding Place and Transition nodes
%   options - solver options
%
% Returns:
%   terms - struct; see the field list assembled at the end of this function
%
% See also SOLVER_FLUID_PETRI, FLUID_MINMULTI_CLOSURE, FLUID_MOMENT_TERMS.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
K = sn.nclasses;
I = sn.nnodes;

places = find(sn.nodetype == NodeType.Place);
places = places(:)';
transitions = find(sn.nodetype == NodeType.Transition);
transitions = transitions(:)';

% ---- the initial marking, read off the model state exactly as the NRM does
m0full = zeros(I, K);
for ind = places
    isf = sn.nodeToStateful(ind);
    [~, nir] = State.toMarginalAggr(sn, ind, sn.state{isf});
    for k = 1:K
        if isinf(nir(k))
            line_error(mfilename, sprintf('Place %s holds an infinite initial marking of class %d.', sn.nodenames{ind}, k));
        end
        m0full(ind, k) = nir(k);
    end
end

% ---- which (place, class) pairs carry a coordinate
touched = false(I, K);
for ind = transitions
    np = sn.nodeparam{ind};
    for m = 1:np.nmodes
        en = i_pad(np.enabling{m}, I, K, 0);
        fir = i_pad(np.firing{m}, I, K, 0);
        inh = i_pad(np.inhibiting{m}, I, K, Inf);
        touched = touched | (en > 0) | (fir ~= 0) | (isfinite(inh) & inh > 0);
    end
end
% an arrival makes its target a coordinate even when no arc mentions it
srcArr = zeros(0, 4); % [source node, source class, place node, place class]
for ind = 1:I
    if sn.nodetype(ind) ~= NodeType.Source
        continue
    end
    ist = sn.nodeToStation(ind);
    for r = 1:K
        lambda = sn.rates(ist, r);
        if isnan(lambda) || lambda <= 0
            continue
        end
        for jnd = places
            for s = 1:K
                p = sn.rtnodes((ind-1)*K + r, (jnd-1)*K + s);
                if p > 0
                    touched(jnd, s) = true;
                    srcArr(end+1, :) = [ind, r, jnd, s]; %#ok<AGROW>
                end
            end
        end
    end
end
keep = false(I, K);
for ind = places
    keep(ind, :) = touched(ind, :) | (m0full(ind, :) > 0);
end

pidx = zeros(I, K);
coordNode = []; coordClass = []; coordStation = [];
nm = 0;
for ind = places
    for k = 1:K
        if ~keep(ind, k)
            continue
        end
        nm = nm + 1;
        pidx(ind, k) = nm;
        coordNode(nm,1) = ind; %#ok<AGROW>
        coordClass(nm,1) = k; %#ok<AGROW>
        coordStation(nm,1) = sn.nodeToStation(ind); %#ok<AGROW>
    end
end

% ---- the modes, and the phase coordinates of the multi-phase ones
modes = struct('node',{},'mode',{},'timing',{},'arcSlot',{},'arcW',{},'inhSlot',{}, ...
    'inhThr',{},'c',{},'nph',{},'D0',{},'D1',{},'d1',{},'pie',{},'dep',{},'prio',{}, ...
    'weight',{},'Cvec',{},'zblk',{},'closable',{},'label',{});
nstate = nm;
for ind = transitions
    np = sn.nodeparam{ind};
    for m = 1:np.nmodes
        rec = i_buildmode(sn, np, ind, m, pidx, nm, I, K);
        if rec.nph > 1
            rec.zblk = nstate + (1:rec.nph);
            nstate = nstate + rec.nph;
        end
        modes(end+1) = rec; %#ok<AGROW>
    end
end
% The phase coordinates are appended after every marking coordinate, so a jump
% column built at nm width has to be grown once the total is known. GUARD THE
% GROWTH: `Cvec(nstate,1) = 0` on a vector that is ALREADY nstate long does not
% grow it, it ZEROES ITS LAST ENTRY -- which is every mode's arc at the last
% place on a net with no multi-phase mode, so every such net answered with an
% empty marking and a conservation row that said so.
for j = 1:numel(modes)
    if numel(modes(j).Cvec) < nstate
        modes(j).Cvec(nstate,1) = 0;
    end
end

timedIdx = find([modes.timing] == TimingStrategy.TIMED);
immIdx = find([modes.timing] == TimingStrategy.IMMEDIATE);

% ---- the event columns
D = zeros(nstate, 0);
rateBase = zeros(0,1);
evKind = zeros(0,1); evMode = zeros(0,1); evPhase = zeros(0,1); evTo = zeros(0,1);
evStation = zeros(0,1); evClass = zeros(0,1);
for j = timedIdx
    md = modes(j);
    if md.nph == 1
        col = md.Cvec;
        D(:,end+1) = col; %#ok<AGROW>
        rateBase(end+1,1) = md.d1(1); %#ok<AGROW>
        evKind(end+1,1) = 1; evMode(end+1,1) = j; evPhase(end+1,1) = 1; evTo(end+1,1) = 1; %#ok<AGROW>
        evStation(end+1,1) = 0; evClass(end+1,1) = 0; %#ok<AGROW>
    else
        for h = 1:md.nph
            for hp = 1:md.nph
                w = md.D1(h,hp);
                if w <= 0
                    continue
                end
                col = md.Cvec;
                col(md.zblk(hp)) = col(md.zblk(hp)) + 1;
                col(md.zblk(h)) = col(md.zblk(h)) - 1;
                D(:,end+1) = col; %#ok<AGROW>
                rateBase(end+1,1) = w; %#ok<AGROW>
                evKind(end+1,1) = 1; evMode(end+1,1) = j; evPhase(end+1,1) = h; evTo(end+1,1) = hp; %#ok<AGROW>
                evStation(end+1,1) = 0; evClass(end+1,1) = 0; %#ok<AGROW>
            end
        end
        for h = 1:md.nph
            for hp = 1:md.nph
                if hp == h
                    continue
                end
                w = md.D0(h,hp);
                if w <= 0
                    continue
                end
                col = zeros(nstate,1);
                col(md.zblk(hp)) = 1;
                col(md.zblk(h)) = -1;
                D(:,end+1) = col; %#ok<AGROW>
                rateBase(end+1,1) = w; %#ok<AGROW>
                evKind(end+1,1) = 2; evMode(end+1,1) = j; evPhase(end+1,1) = h; evTo(end+1,1) = hp; %#ok<AGROW>
                evStation(end+1,1) = 0; evClass(end+1,1) = 0; %#ok<AGROW>
            end
        end
    end
end
% One latch column per multi-phase mode: mu_j servers per unit time enter at the
% firing process's own entry distribution. The rate is FREE IN SIGN -- a mode
% whose enabling degree drops stops servers rather than starting them -- and it
% is zero at any fixed point, which is what the latch row enforces.
for j = timedIdx
    md = modes(j);
    if md.nph <= 1
        continue
    end
    col = zeros(nstate,1);
    col(md.zblk) = md.pie(:);
    D(:,end+1) = col; %#ok<AGROW>
    rateBase(end+1,1) = 1; %#ok<AGROW>
    evKind(end+1,1) = 5; evMode(end+1,1) = j; evPhase(end+1,1) = 0; evTo(end+1,1) = 0; %#ok<AGROW>
    evStation(end+1,1) = 0; evClass(end+1,1) = 0; %#ok<AGROW>
end
for j = immIdx
    D(:,end+1) = modes(j).Cvec; %#ok<AGROW>
    rateBase(end+1,1) = 1; %#ok<AGROW>
    evKind(end+1,1) = 4; evMode(end+1,1) = j; evPhase(end+1,1) = 0; evTo(end+1,1) = 0; %#ok<AGROW>
    evStation(end+1,1) = 0; evClass(end+1,1) = 0; %#ok<AGROW>
end
for a = 1:size(srcArr,1)
    snd = srcArr(a,1); r = srcArr(a,2); qnd = srcArr(a,3); l = srcArr(a,4);
    ist = sn.nodeToStation(snd);
    if isfield(sn,'procid') && sn.procid(ist, r) ~= ProcessType.EXP
        line_error(mfilename, sprintf(['Source %s has a non-exponential arrival for class %d. The fluid ' ...
            'Petri route models an arrival as a constant-propensity event, which a renewal stream with ' ...
            'memory is not; use SolverCTMC, SolverJMT or SolverSSA.'], sn.nodenames{snd}, r));
    end
    col = zeros(nstate,1);
    col(pidx(qnd,l)) = 1;
    D(:,end+1) = col; %#ok<AGROW>
    rateBase(end+1,1) = sn.rates(ist,r) * sn.rtnodes((snd-1)*K + r, (qnd-1)*K + l); %#ok<AGROW>
    evKind(end+1,1) = 3; evMode(end+1,1) = 0; evPhase(end+1,1) = 0; evTo(end+1,1) = 0; %#ok<AGROW>
    evStation(end+1,1) = ist; evClass(end+1,1) = r; %#ok<AGROW>
end
nev = size(D,2);

% ---- which Sigma entries the closure reads
% One unknown per (a,b) the drift actually looks at, and no more: the input
% coordinates of a closable mode pairwise, plus the variance of every inhibitor
% coordinate whose gate is smoothed. This is the direct analogue of the one
% sigma2 per closable station the queueing DAE carries, and it is what keeps
% Sigma out of the Newton vector.
pairKey = containers.Map('KeyType','char','ValueType','double');
covPairs = zeros(0,2);
for j = 1:numel(modes)
    md = modes(j);
    if md.closable
        for a = 1:numel(md.arcSlot)
            for b = a:numel(md.arcSlot)
                [covPairs, pairKey] = i_addpair(covPairs, pairKey, md.arcSlot(a), md.arcSlot(b));
            end
        end
    end
    if md.timing == TimingStrategy.TIMED
        for b = 1:numel(md.inhSlot)
            [covPairs, pairKey] = i_addpair(covPairs, pairKey, md.inhSlot(b), md.inhSlot(b));
        end
    end
end
npair = size(covPairs,1);

% ---- per-place consumption and production, for the metric reader
% A Place's throughput is the rate at which TOKENS leave it, so each consuming
% mode contributes its firing rate times the multiplicity of the arc it takes
% them through. That is the convention SolverCTMC reports and the one Little's
% law needs: on a place drained by one mode through an arc of weight two, the
% exact CTMC reports 3.158 where the firing rate is 1.579.
%
% NOTE THAT SOLVER_SSA'S NRM REPORTS THE OTHER ONE -- its CONSUMERS accumulator
% sums the propensity once per consuming mode, unweighted -- so the two solvers
% disagree at any place whose input arc has a multiplicity above one. The exact
% engine is the reference here.
consumers = cell(M, K);
consumerW = cell(M, K);
producers = cell(M, K);
for e = 1:nev
    if evKind(e) == 3
        producers{evStation(e), evClass(e)}(end+1) = e;
        continue
    end
    if evKind(e) ~= 1 && evKind(e) ~= 4
        continue
    end
    md = modes(evMode(e));
    for a = 1:numel(md.arcSlot)
        s = md.arcSlot(a);
        consumers{coordStation(s), coordClass(s)}(end+1) = e;
        consumerW{coordStation(s), coordClass(s)}(end+1) = md.arcW(a);
    end
end

% ---- the initial state
x0 = zeros(nstate,1);
for s = 1:nm
    x0(s) = m0full(coordNode(s), coordClass(s));
end
for j = 1:numel(modes)
    md = modes(j);
    if md.nph > 1
        % the servers the INITIAL marking latches, in the phase they start in
        e = Inf;
        for a = 1:numel(md.arcSlot)
            e = min(e, x0(md.arcSlot(a)) / md.arcW(a));
        end
        if isempty(md.arcSlot)
            e = 1;
        end
        x0(md.zblk) = min(e, md.c) * md.pie(:);
    end
end

terms = struct();
terms.M = M; terms.K = K; terms.I = I;
terms.places = places; terms.transitions = transitions;
terms.namesNode = sn.nodenames;
terms.nstate = nstate; terms.nm = nm;
terms.pidx = pidx;
terms.coordNode = coordNode; terms.coordClass = coordClass; terms.coordStation = coordStation;
terms.modes = modes; terms.timedIdx = timedIdx; terms.immIdx = immIdx;
terms.D = D; terms.rateBase = rateBase; terms.nev = nev;
terms.evKind = evKind; terms.evMode = evMode; terms.evPhase = evPhase; terms.evTo = evTo;
terms.evStation = evStation; terms.evClass = evClass;
terms.immCol = find(evKind == 4);
terms.latchCol = find(evKind == 5);
% The DIFFUSION counts the stochastic events only: an immediate flow and a
% server latch are both the limit of an infinitely fast mechanism whose
% fluctuation is slaved, not a Poisson stream with an intensity.
terms.stochCol = find(evKind ~= 4 & evKind ~= 5);
terms.latchMode = zeros(1,0);
for j = timedIdx
    if modes(j).nph > 1
        terms.latchMode(end+1) = j;
    end
end
% THE COVARIANCE COVERS EVERY COORDINATE, phases included. The in-flight firings
% of a mode are a population like any other, and dropping them from the Lyapunov
% solve severs the only path by which the marking feeds back into a multi-phase
% mode's firing rate -- the rate is linear in y and does not read the marking at
% all, so the reduced generator would lose that mode's whole restoring force. On
% a closed cycle with Erlang(2) firing that reads 1.4224 against the exact
% 1.4531, where carrying the phases reads 1.4564. What the latch removes is one
% DIRECTION per multi-phase mode, and LOCAL_CLAMP_TANGENT removes exactly that.
terms.covIdx = (1:nstate)';
terms.covPairs = covPairs; terms.npair = npair;
% (a,b) -> unknown index, symmetric, so the closure reads one Sigma entry in O(1)
terms.pairIndex = sparse(nm, nm);
for t = 1:npair
    terms.pairIndex(covPairs(t,1), covPairs(t,2)) = t;
    terms.pairIndex(covPairs(t,2), covPairs(t,1)) = t;
end
% the D column carrying each immediate mode's algebraic flow
terms.immColOf = zeros(numel(modes),1);
for t = 1:numel(terms.immCol)
    terms.immColOf(evMode(terms.immCol(t))) = terms.immCol(t);
end
terms.consumers = consumers; terms.consumerW = consumerW; terms.producers = producers;
terms.x0 = x0;
terms.m0full = m0full;
terms.options = options;
end

% -------------------------------------------------------------------------
function rec = i_buildmode(sn, np, ind, m, pidx, nm, I, K)
% One transition mode as a reaction record: its input arcs and their weights,
% its inhibitor arcs and their thresholds, its incidence column, and the firing
% process that times it.
rec = struct('node',ind,'mode',m,'timing',np.timing(m),'arcSlot',[],'arcW',[], ...
    'inhSlot',[],'inhThr',[],'c',1,'nph',1,'D0',[],'D1',[],'d1',[],'pie',1, ...
    'dep',[],'prio',1,'weight',1,'Cvec',zeros(nm,1),'zblk',[],'closable',false,'label','');
rec.label = sprintf('%s.%s', sn.nodenames{ind}, i_modename(np, m));

en = i_pad(np.enabling{m}, I, K, 0);
fir = i_pad(np.firing{m}, I, K, 0);
inh = i_pad(np.inhibiting{m}, I, K, Inf);

Cvec = zeros(nm,1);
[pr, cr] = find(en > 0);
for t = 1:numel(pr)
    p = pr(t); k = cr(t);
    w = en(p,k);
    if ~isfinite(w)
        line_error(mfilename, sprintf(['Mode %s has a non-finite enabling arc weight at %s. An arc that ' ...
            'no marking can satisfy disables the mode; declare a finite multiplicity.'], rec.label, sn.nodenames{p}));
    end
    if pidx(p,k) == 0
        line_error(mfilename, sprintf('Mode %s takes an enabling arc from %s, which is not a Place.', rec.label, sn.nodenames{p}));
    end
    rec.arcSlot(end+1) = pidx(p,k); %#ok<AGROW>
    rec.arcW(end+1) = w; %#ok<AGROW>
    Cvec(pidx(p,k)) = Cvec(pidx(p,k)) - w;
end
[pf, cf] = find(fir ~= 0);
for t = 1:numel(pf)
    p = pf(t); k = cf(t);
    if fir(p,k) <= 0
        continue % a negative entry marks an input place, whose token the PRE already removed
    end
    if pidx(p,k) == 0
        continue % a firing arc into a Sink is mass leaving the net
    end
    Cvec(pidx(p,k)) = Cvec(pidx(p,k)) + fir(p,k);
end
[ph, ch] = find(isfinite(inh));
for t = 1:numel(ph)
    p = ph(t); k = ch(t);
    if inh(p,k) < 0 || pidx(p,k) == 0
        continue % JMT writes a missing inhibitor arc as 0 or -1, never as a threshold
    end
    if inh(p,k) == 0
        continue
    end
    rec.inhSlot(end+1) = pidx(p,k); %#ok<AGROW>
    rec.inhThr(end+1) = inh(p,k); %#ok<AGROW>
end
rec.Cvec = Cvec;

c = np.nmodeservers(m);
if isempty(c) || isnan(c)
    c = 1;
end
rec.c = c;
rec.prio = np.firingprio(m);
rec.weight = np.fireweight(m);
if isfield(np,'firingdep') && numel(np.firingdep) >= m && ~isempty(np.firingdep{m})
    rec.dep = np.firingdep{m};
end

if rec.timing == TimingStrategy.IMMEDIATE
    rec.nph = 0;
    rec.d1 = 0;
    % An immediate mode carries no firing process: its flow is an algebraic
    % unknown of the DAE, not a rate.
    rec.closable = false;
    return
end

if isempty(np.firingproc{m}) || isnan(np.firingphases(m))
    line_error(mfilename, sprintf(['Mode %s has no Markovian firing process. SN_NONMARKOV_TOPH converts the ' ...
        'renewal families to phase type before the solver runs, so this is a distribution the fluid ' ...
        'Petri route cannot time; use SolverCTMC or SolverLDES.'], rec.label));
end
D0 = np.firingproc{m}{1};
D1 = np.firingproc{m}{2};
rec.nph = size(D0,1);
rec.D0 = D0; rec.D1 = D1;
rec.d1 = sum(D1,2);
if rec.nph > 1
    pie = np.firingpie{m};
    if isempty(pie)
        pie = [1, zeros(1, rec.nph-1)];
    end
    rec.pie = pie(:)'/sum(pie);
end

% A mode whose enabling degree cannot reach its server count has min() exact on
% its whole support, so closing it is not an improvement but an error, exactly
% as a station that cannot fill its servers is held first order in
% FLUID_MOMENT_TERMS. A single input arc with an unbounded server count is the
% common case, and it makes the whole drift LINEAR in that mode.
rec.closable = ~(numel(rec.arcSlot) <= 1 && ~isfinite(rec.c));
end

% -------------------------------------------------------------------------
function s = i_modename(np, m)
if isfield(np,'modenames') && numel(np.modenames) >= m
    s = char(np.modenames{m});
else
    s = sprintf('Mode%d', m);
end
end

% -------------------------------------------------------------------------
function A = i_pad(A, I, K, fillv)
% Pad an arc matrix to (nnodes x nclasses): ADDMODE sizes them at creation
% time, so a node or class added later leaves them short.
if size(A,1) < I || size(A,2) < K
    B = fillv*ones(I,K);
    B(1:size(A,1), 1:size(A,2)) = A;
    A = B;
end
A = A(1:I, 1:K);
end

% -------------------------------------------------------------------------
function [covPairs, pairKey] = i_addpair(covPairs, pairKey, a, b)
lo = min(a,b); hi = max(a,b);
key = sprintf('%d_%d', lo, hi);
if ~isKey(pairKey, key)
    covPairs(end+1,:) = [lo, hi];
    pairKey(key) = size(covPairs,1);
end
end
