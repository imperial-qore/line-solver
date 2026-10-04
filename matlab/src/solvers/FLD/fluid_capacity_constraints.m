function con = fluid_capacity_constraints(sn, terms)
% CON = FLUID_CAPACITY_CONSTRAINTS(SN, TERMS)
%
% Every capacity limit in the model as a linear constraint on the fluid state.
%
% TWO DECLARATIONS, ONE FAMILY OF ROWS. A finite capacity region caps a SET of
% stations jointly (addRegion); a station cap limits one station's own buffer
% (setCapacity/setClassCapacity). LINE stores the region form four ways -- a
% region-global job cap, a per-class job cap, a memory budget with per-class
% sizes, and an arbitrary linear constraint pair (A,b) -- and the station form
% two, a total and a per-class buffer. All six are the same object once written
% against the state:
%
%     Arow * x <= b,    Arow(s) = weight of the class of coordinate s,
%                                 zero outside the stations the limit covers
%
% so the solver carries one mechanism rather than six. The global cap is that
% row with every weight one, a per-class cap is the row with one class weighted,
% the memory budget weights each class by its size, and the linear constraint
% pair supplies the weights directly.
%
% WHAT DIFFERS IS NOT THE ROW BUT WHERE THE BLOCKED JOB GOES, which is decided
% per admission event in FLUID_CAPACITY_GATES rather than here.
%
% WHAT IS REFUSED, AND WHY IT IS NOT A CONSTRAINT. Only a waiting queue or a
% drop is a constraint on this drift. Under WAITQ a blocked job waits and is
% admitted later, so the population is conserved and only the admission FLOW is
% throttled; under DROP the flow the cap will not take is discarded, which is
% what LINE does with an open arrival at a full buffer. BAS/BBS/RSRD give the
% upstream station a blocked-server state and the retrial rules add an orbit.
% Each of those changes the event set itself, so each needs a different drift
% rather than a constraint on this one, and is refused by name.
%
% Parameters:
%   sn    - NetworkStruct, carrying nregions/region/regionmembers/regionrule/
%           regionweight/regionsz/regionmaxmem/regionlincon and cap/classcap/
%           droprule
%   terms - the event representation from FLUID_MOMENT_TERMS
%
% Returns:
%   con - struct with fields
%           A       (ncon x nstate) constraint rows in state space
%           As      (ncon x nstg)   the same rows over the staging coordinates,
%                                   filled in by FLUID_CAPACITY_EXTEND
%           b       (ncon x 1)      right-hand sides
%           region  (ncon x 1)      which region each row came from, 0 for a
%                                   station cap
%           station (ncon x 1)      which station, 0 for a region cap
%           class   (ncon x 1)      which class it limits, 0 when several
%           staged  (ncon x 1)      true where the blocked job waits in a room
%           label   (ncon x 1) cell human-readable origin, for messages
%           member  (nregions x nstate) region membership, for the rooms
%           coordClass (1 x nstate) the class of each coordinate
%         empty A/b when the model has no limit
%
% See also SOLVER_FLUID_DAE, FLUID_CAPACITY_GATES, FLUID_CAPACITY_STAGING,
% FLUID_CAPACITY_EXTEND, FLUID_MOMENT_TERMS.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = terms.M;
K = terms.K;
nstate = terms.nstate;

nregions = 0;
if isfield(sn,'nregions') && ~isempty(sn.nregions)
    nregions = sn.nregions;
end

% state -> (station, class), rebuilt from the layout TERMS indexes
coordStation = zeros(nstate,1);
coordClass = zeros(nstate,1);
for i = 1:M
    for c = 1:K
        if terms.Kic(i,c) > 0
            idx = terms.q_indices(i,c):(terms.q_indices(i,c)+terms.Kic(i,c)-1);
            coordStation(idx) = i;
            coordClass(idx) = c;
        end
    end
end

con = struct('A', zeros(0,nstate), 'As', zeros(0,0), 'b', zeros(0,1), ...
    'region', zeros(0,1), 'station', zeros(0,1), 'class', zeros(0,1), ...
    'staged', false(0,1), 'label', {cell(0,1)}, ...
    'member', false(max(nregions,1),nstate), 'coordClass', coordClass(:)', ...
    'nregions', nregions);

UNB = FiniteCapacityRegion.UNBOUNDED;
A = zeros(0,nstate); b = zeros(0,1);
rg = zeros(0,1); st = zeros(0,1); kl = zeros(0,1); sd = false(0,1);
lb = cell(0,1);
member = false(max(nregions,1), nstate);

% -- 1..4 the region forms ---------------------------------------------------
for f = 1:nregions
    members = false(M,1);
    if isfield(sn,'regionmembers') && numel(sn.regionmembers) >= f && ~isempty(sn.regionmembers{f})
        members = logical(sn.regionmembers{f}(:));
    end
    if ~any(members)
        continue
    end
    inRegion = members(max(coordStation,1)) & coordStation > 0;
    member(f,:) = inRegion(:)';

    % -- refusals, named before anything is built -----------------------------
    for r = 1:K
        if isfield(sn,'regionrule') && size(sn.regionrule,1) >= f && size(sn.regionrule,2) >= r
            rule = sn.regionrule(f,r);
            if rule ~= DropStrategy.WAITQ
                line_error(mfilename, sprintf(['Region %d applies %s to class %d. Only a waiting queue is a ' ...
                    'constraint on the fluid drift: it conserves the population and throttles the admission ' ...
                    'flow. %s changes the event set instead, so it needs a different drift rather than an ' ...
                    'algebraic equation on this one.'], f, DropStrategy.toText(rule), r, DropStrategy.toText(rule)));
            end
        end
    end
    if isfield(sn,'regionweight') && size(sn.regionweight,1) >= f
        wf = sn.regionweight(f,:);
        if any(abs(wf(:) - 1) > GlobalConstants.Zero)
            line_error(mfilename, sprintf(['Region %d sets per-class admission weights, which decide WHICH ' ...
                'blocked class enters when capacity frees up. The constraint form throttles the admission ' ...
                'flow in proportion to its own rate and carries no such priority, so the weights would be ' ...
                'ignored silently. Use SolverCTMC, SolverJMT or SolverSSA.'], f));
        end
    end

    Rmat = [];
    if isfield(sn,'region') && numel(sn.region) >= f
        Rmat = sn.region{f};
    end
    memberRow = find(members, 1);

    % -- 1. region-global job cap --------------------------------------------
    if ~isempty(Rmat) && size(Rmat,2) >= K+1
        gcap = Rmat(memberRow, K+1);
        if isfinite(gcap) && gcap ~= UNB && gcap >= 0
            row = zeros(1,nstate); row(inRegion) = 1;
            A(end+1,:) = row; b(end+1,1) = gcap; rg(end+1,1) = f; st(end+1,1) = 0; %#ok<AGROW>
            kl(end+1,1) = 0; sd(end+1,1) = true; %#ok<AGROW>
            lb{end+1,1} = sprintf('region %d global job cap', f); %#ok<AGROW>
        end
    end

    % -- 2. per-class job caps ------------------------------------------------
    if ~isempty(Rmat)
        for r = 1:K
            if size(Rmat,2) < r, continue, end
            ccap = Rmat(memberRow, r);
            if isfinite(ccap) && ccap ~= UNB && ccap >= 0
                sel = inRegion & (coordClass == r);
                if ~any(sel), continue, end
                row = zeros(1,nstate); row(sel) = 1;
                A(end+1,:) = row; b(end+1,1) = ccap; rg(end+1,1) = f; st(end+1,1) = 0; %#ok<AGROW>
                kl(end+1,1) = r; sd(end+1,1) = true; %#ok<AGROW>
                lb{end+1,1} = sprintf('region %d class %d job cap', f, r); %#ok<AGROW>
            end
        end
    end

    % -- 3. region-global memory budget --------------------------------------
    % REFRESHREGIONS has already folded a per-class memory limit into the
    % per-class job cap above (memjobs = floor(classMaxMemory/classSize)), so
    % only the region-global budget is left, and it is the job-count row with
    % each class weighted by its size.
    if isfield(sn,'regionmaxmem') && numel(sn.regionmaxmem) >= f && ~isempty(sn.regionmaxmem{f})
        mem = sn.regionmaxmem{f};
        gmem = mem(memberRow);
        if isfinite(gmem) && gmem ~= UNB && gmem >= 0
            sz = ones(1,K);
            if isfield(sn,'regionsz') && size(sn.regionsz,1) >= f
                sz = sn.regionsz(f,:);
            end
            row = zeros(1,nstate);
            for r = 1:K
                sel = inRegion & (coordClass == r);
                row(sel) = sz(r);
            end
            if any(row ~= 0)
                A(end+1,:) = row; b(end+1,1) = gmem; rg(end+1,1) = f; st(end+1,1) = 0; %#ok<AGROW>
                kl(end+1,1) = 0; sd(end+1,1) = true; %#ok<AGROW>
                lb{end+1,1} = sprintf('region %d memory budget', f); %#ok<AGROW>
            end
        end
    end

    % -- 4. explicit linear constraints --------------------------------------
    if isfield(sn,'regionlincon') && size(sn.regionlincon,1) >= f && ~isempty(sn.regionlincon{f,1})
        Alin = sn.regionlincon{f,1};
        blin = sn.regionlincon{f,2};
        for c = 1:size(Alin,1)
            row = zeros(1,nstate);
            for r = 1:min(K,size(Alin,2))
                if Alin(c,r) ~= 0
                    sel = inRegion & (coordClass == r);
                    row(sel) = Alin(c,r);
                end
            end
            if any(row ~= 0)
                A(end+1,:) = row; b(end+1,1) = blin(c); rg(end+1,1) = f; st(end+1,1) = 0; %#ok<AGROW>
                kl(end+1,1) = 0; sd(end+1,1) = true; %#ok<AGROW>
                lb{end+1,1} = sprintf('region %d linear constraint %d', f, c); %#ok<AGROW>
            end
        end
    end
end

% -- 5. per-station buffers --------------------------------------------------
% A STATION CAP IS THE ONE-STATION CASE OF THE SAME ROW, and every fluid method
% other than this one ignores it outright: nothing in the FLD tree reads sn.cap
% or sn.classcap, so a capped station was integrated as an unbounded one and the
% table reported more jobs in the buffer than the buffer holds (a closed
% Delay->Queue(cap 2) model returned 2.19 jobs in a buffer of 2).
cap = [];
if isfield(sn,'cap') && ~isempty(sn.cap)
    cap = sn.cap(:);
end
classcap = [];
if isfield(sn,'classcap') && ~isempty(sn.classcap)
    classcap = sn.classcap;
end
for i = 1:M
    if terms.isExt(i)
        continue % a source holds no jobs, so its cap caps nothing
    end
    atStation = coordStation == i;
    if ~any(atStation)
        continue
    end
    if numel(cap) >= i
        gcap = cap(i);
        if isfinite(gcap) && gcap >= 0
            row = zeros(1,nstate); row(atStation) = 1;
            A(end+1,:) = row; b(end+1,1) = gcap; rg(end+1,1) = 0; st(end+1,1) = i; %#ok<AGROW>
            kl(end+1,1) = 0; sd(end+1,1) = false; %#ok<AGROW>
            lb{end+1,1} = sprintf('station %d buffer', i); %#ok<AGROW>
        end
    end
    if size(classcap,1) >= i
        for r = 1:min(K,size(classcap,2))
            ccap = classcap(i,r);
            if ~isfinite(ccap) || ccap < 0
                continue
            end
            sel = atStation & (coordClass == r);
            if ~any(sel), continue, end
            row = zeros(1,nstate); row(sel) = 1;
            A(end+1,:) = row; b(end+1,1) = ccap; rg(end+1,1) = 0; st(end+1,1) = i; %#ok<AGROW>
            kl(end+1,1) = r; sd(end+1,1) = false; %#ok<AGROW>
            lb{end+1,1} = sprintf('station %d class %d buffer', i, r); %#ok<AGROW>
        end
    end
end

con.member = member;
if isempty(A)
    con.As = zeros(0,0);
    return
end

% A CAP THE POPULATION CANNOT REACH IS NOT A CAP, and dropping it here keeps it
% out of the active-set loop, where it would be tested on every pass and never
% bind while its multiplier stayed an unknown with no equation to pin it.
%
% THE BOUND IS PER CLASS, and it has to be: the heaviest weight times the whole
% population, which is what the region-only version used, never prunes a
% PER-CLASS cap set to its own class population -- exactly the row
% refreshCapacity derives at every station of every closed model. An open class
% makes the reachable total unbounded, so nothing on a row it weighs is pruned.
keep = true(size(b));
for c = 1:numel(b)
    if local_reach(A(c,:), coordClass, coordStation, sn.njobs, K) <= b(c) + GlobalConstants.Zero
        keep(c) = false;
    end
end

% A CAP ON COORDINATES THE IMMEDIATE REDUCTION FOLDED AWAY IS VACUOUS, NOT
% MALFORMED. options.config.hide_immediate stochastic-complements an
% Immediate-rate coordinate out of the event set -- the MMT transform's
% zero-service Join is one, until the fork-join fixed point gives it a
% synchronisation delay -- and the reduced drift then holds no mass there and has
% no event landing on it. Left in, such a row reaches FLUID_CAPACITY_GATES with no
% gating event and is reported as a limit the model cannot approach, which is the
% right message for a station that really does hold jobs and the wrong one here:
% this station holds none, so its buffer is satisfied identically.
gone = local_eliminated(terms, size(A,2));
for c = 1:numel(b)
    if keep(c) && ~any(A(c,:) ~= 0 & ~gone)
        keep(c) = false;
    end
end

% THE SAME ROW TWICE IS A SINGULAR NEWTON SYSTEM, not a redundancy the least
% squares absorbs: two identical rows both bind, each takes a multiplier, and
% nothing distinguishes them. refreshCapacity derives classcap from cap, so a
% single-class model declares the station total and the class buffer as the same
% row -- the common case, not a corner one. The TIGHTER bound survives; on a tie
% the region row does, because a region is an explicit construct with a waiting
% room of its own while a station cap is a buffer length.
for c = 1:numel(b)
    if ~keep(c), continue, end
    for d = c+1:numel(b)
        if ~keep(d) || any(abs(A(c,:) - A(d,:)) > GlobalConstants.Zero)
            continue
        end
        takeover = b(d) < b(c) - GlobalConstants.Zero || ...
            (abs(b(d) - b(c)) <= GlobalConstants.Zero && sd(d) && ~sd(c));
        if takeover
            b(c) = b(d); rg(c) = rg(d); st(c) = st(d); kl(c) = kl(d); sd(c) = sd(d);
            lb{c} = lb{d};
        end
        keep(d) = false;
    end
end

% A ROW ITS OWN PER-CLASS ROWS ALREADY IMPLY IS RANK, NOT INFORMATION. LINE
% derives the station total from the per-class buffers (cap = sum classcap), so a
% two-class station capped 3 and 3 declares a total of 6 as well -- and the total
% row is EXACTLY the sum of the two class rows, with a bound exactly their sum.
% All three then bind together, the constraint block has rank 2 with 3
% multipliers, and nothing distinguishes them: the least-squares step still lands
% on the right STATE, but the multipliers it reports are one arbitrary point of a
% line of solutions. Only an exact implication is dropped -- the per-class rows
% must cover every class the candidate weighs, at no smaller a weight, and their
% bounds must sum to no more than its own -- so a total TIGHTER than the sum of
% its parts, which is a genuine extra constraint, survives.
for c = 1:numel(b)
    if ~keep(c) || kl(c) ~= 0
        continue
    end
    classes = [];
    for r = 1:K
        if any((coordClass == r) & (coordStation > 0) & (A(c,:)' > 0))
            classes(end+1) = r; %#ok<AGROW>
        end
    end
    if isempty(classes)
        continue
    end
    implied = true; budget = 0;
    for r = classes
        sel = (coordClass == r) & (coordStation > 0) & (A(c,:)' > 0);
        part = 0;
        for d = 1:numel(b)
            if ~keep(d) || d == c || kl(d) ~= r || rg(d) ~= rg(c) || st(d) ~= st(c)
                continue
            end
            if all(A(d,sel)' >= A(c,sel)' - GlobalConstants.Zero)
                part = d;
                break
            end
        end
        if part == 0
            implied = false;
            break
        end
        budget = budget + b(part);
    end
    if implied && budget <= b(c) + GlobalConstants.Zero
        keep(c) = false;
    end
end

A = A(keep,:); b = b(keep); rg = rg(keep); st = st(keep); kl = kl(keep);
sd = sd(keep); lb = lb(keep);

% WHICH STATION RULES ARE A CONSTRAINT ON THIS DRIFT, and which are a different
% event set. THE RULE IS NOT WHAT DECIDES THE SEMANTICS -- the class type is,
% exactly as State.arrivalIsLost decides it for every other solver: a closed job
% is never lost (the upstream departure is disabled instead) and an open one is
% never held (it simply never enters). So a waiting queue and a drop declaration
% are BOTH constraints here and differ only in which of those two the class
% already implies; refreshCapacity in fact declares DROP by default at a capped
% station reached by an open class, so refusing DROP would refuse every open loss
% model. What is refused is the rules that add STATE.
%
% A rule is only a contradiction where the cap can BIND, so this runs on the
% surviving rows: a declaration on a buffer the population can never fill
% describes nothing, and refusing it would reject models that behave identically
% with and without it.
if isfield(sn,'droprule') && ~isempty(sn.droprule)
    for c = 1:numel(b)
        if st(c) == 0 || size(sn.droprule,1) < st(c)
            continue
        end
        for r = 1:min(K,size(sn.droprule,2))
            sel = (coordStation == st(c)) & (coordClass == r);
            if ~any(sel) || max(A(c,sel)) <= 0
                continue
            end
            rule = sn.droprule(st(c), r);
            if rule ~= DropStrategy.WAITQ && rule ~= DropStrategy.DROP
                line_error(mfilename, sprintf(['Station %d applies %s to class %d, and its buffer binds. ' ...
                    'Only a waiting queue or a drop is a constraint on this drift: the first conserves the ' ...
                    'population and throttles the admission flow, the second discards the flow the cap will ' ...
                    'not take. BAS/BBS/RSRD add a blocked-server state to the upstream station and the ' ...
                    'retrial rules add an orbit, so each needs a different drift rather than an algebraic ' ...
                    'equation on this one. Use SolverCTMC, SolverJMT, SolverSSA or SolverLDES.'], ...
                    st(c), DropStrategy.toText(rule), r));
            end
        end
    end
end

con.A = A;
con.As = zeros(numel(b),0);
con.b = b;
con.region = rg;
con.station = st;
con.class = kl;
con.staged = sd;
con.label = lb;
end

% ---------------------------------------------------------------------------
function gone = local_eliminated(terms, nstate)
% GONE(s) is true where the immediate reduction folded coordinate s away, so the
% reduced drift holds no mass there and no event lands on it.
%
% TERMS.ABSORB is the projector the reduction returns: the identity on a
% surviving coordinate and the absorption distribution on an eliminated one, so a
% zero diagonal is exactly the eliminated case. It is empty when nothing was
% eliminated, where every coordinate survives.
gone = false(1, nstate);
if ~isfield(terms,'absorb') || isempty(terms.absorb)
    return
end
n = min(nstate, size(terms.absorb,1));
gone(1:n) = full(diag(terms.absorb(1:n,1:n)))' == 0;
end

function total = local_reach(row, coordClass, coordStation, njobs, K)
% The largest ROW*x the population can produce, ignoring the coupling. An upper
% bound is what is wanted: too loose only costs a constraint that stays
% inactive, while too tight would discard a cap that does bind.
total = 0;
for r = 1:K
    sel = (coordClass == r) & (coordStation > 0);
    if ~any(sel)
        continue
    end
    w = max(row(sel));
    if w <= 0
        continue
    end
    nr = Inf;
    if numel(njobs) >= r
        nr = njobs(r);
    end
    if ~isfinite(nr)
        total = Inf;
        return
    end
    total = total + w * nr;
end
end
