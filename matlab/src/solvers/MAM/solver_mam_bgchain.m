function [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_bgchain(sn, options)
% [QN,UN,RN,TN,CN,XN,TOTITER] = SOLVER_MAM_BGCHAIN(SN, OPTIONS)
%
% Mixed-network solver that treats the CLOSED classes as a background
% modulating chain and the OPEN classes as matrix-analytic queues driven by it.
%
% The closed population vector of a mixed network is a finite continuous-time
% Markov chain in its own right: it is the only part of the model whose state
% space is bounded. This method solves it exactly (given the mean open
% occupancy) and hands the open classes a station-local Markovian ENVIRONMENT
% read off that chain, so each open station becomes a level-dependent QBD whose
% phase carries the number of closed jobs competing for its server. The two
% halves meet at a fixed point on the mean open occupancy: the background chain
% is built from it, and the QBDs recompute it.
%
%   1. background chain    closed population vector over the stations the
%                          closed classes visit, with the open classes present
%                          only through their mean occupancy (MAM_BGCHAIN_CTMC)
%   2. environment         the chain lumped onto the closed occupancy of one
%                          station (MAM_BGCHAIN_ENV)
%   3. open station        MAP/PH/c queue modulated by that environment, solved
%                          as a level-dependent QBD (MAM_BGCHAIN_STATION)
%   4. fixed point         the mean open occupancy feeds step 1 and closes
%
% A PURELY CLOSED MODEL IS THE DEGENERATE CASE OF THE SAME CONSTRUCTION, not a
% separate algorithm. With no open class steps 2-4 have nothing to do: CSHARE
% never leaves its initial min(e,c), which is exactly the capacity the closed
% jobs hold when no open work competes for it, so step 1 alone answers and it
% answers with the EXACT closed CTMC at chain granularity. The same limit is
% reached from the mixed side by letting the open arrival rate go to zero, and
% it agrees with exact MVA to the O(lambda) of that perturbation. What it costs
% is the state space below, which is why SOLVER_MAM_ANALYZER sizes the chain
% before choosing this method as the closed default.
%
% TAGGED-CLASS ITERATION. Step 1 is a population process of dimension (number
% of closed chains) x (number of stations), so its state space is exponential
% in the number of closed chains. With R > 1 closed chains the method keeps ONE
% chain free at a time: the tagged chain r is carried exactly, the other R-1
% are replaced by flow-equivalent aggregate classes whose population is their
% total and whose service time and routing at each station are their
% throughput-weighted means (Chandy-Herzog-Woo aggregation). Every chain takes
% its turn as the tagged one and reads its own metrics off the chain it is
% exact in; the open-class results are averaged over the passes, all of which
% estimate the same quantity.
%
% HOW MUCH TO AGGREGATE is options.config.bgaggr, the number G of aggregate
% classes; the background chain then carries 1 + G classes.
%
%   G = 1        the classic tagged/aggregate pair, and the default: the chain
%                stays two-class whatever R is, which is the cheapest option
%   1 < G < R-1  the untagged chains are split into G groups
%   G >= R-1     nothing is aggregated. Every closed chain is a class of its
%                own, the background chain is EXACT in the closed classes, and
%                ONE pass answers for all of them -- the tagged loop would
%                otherwise solve the same chain R times over. Passing R reaches
%                this, so asking for "no aggregation" needs no magic value
%
% What it costs is the state space, the product over the 1 + G classes of
% nchoosek(N_b + Mc - 1, Mc - 1), bounded by options.config.bgstates_max.
%
% WHICH CHAINS SHARE A GROUP is decided by similarity of per-station SERVICE
% DEMAND, in MAM_BGCHAIN_GROUPS. An aggregate carries the flow-weighted mean of
% its members' service times and routing, so it is exact when they place the
% same demand at every station and distorts in proportion to how far apart they
% are. Grouping the demand-similar chains together is what keeps the
% aggregation where it is harmless and away from the chains it would
% misrepresent.
%
% EXACTNESS, as measured against SolverCTMC and exact MVA on mixed models of two
% to four stations:
%
%   PS or INF, ANY service law (exponential, Erlang, HyperExp, Coxian), any
%     number of servers, Poisson or MAP arrivals, one to four closed chains
%                                                agreement to 4-5 significant
%                                                digits (<= 0.06% on every
%                                                queue length), where dec.source
%                                                is 10-24% out
%   FCFS, class-INDEPENDENT service rates        same, 5 significant digits
%   FCFS, class-DEPENDENT service rates          the closed queue lengths stay
%                                                within ~1%, the open queue
%                                                length reads 14-20% low
%
% PS is INSENSITIVE to the service law beyond its mean, and the method honours
% that rather than approximating it: at a PS station the open service is
% replaced by the exponential of the same mean before the QBD is built. That is
% not a shortcut. This QBD tracks ONE service phase for the whole station, so
% carrying the phase-type there makes the open queue length inherit the
% SCV-sensitivity of an M/PH/1 FCFS queue -- measured, a HyperExp of SCV 4 read
% 21% high where the exact answer is the exponential one to five digits.
%
% The remaining residual is FCFS-only and is one named approximation. A station
% serving classes at DIFFERENT rates under FCFS is served here in random order
% -- the server is held by an open job with probability k/(k+e) -- which is
% exact under PS but loses the head-of-line order FCFS actually imposes, and
% with it part of the variability of the wait. At an FCFS station the service
% law IS carried, collapsed into one phase-type process scaled by the share, the
% same collapse SOLVER_MAM_LDQBD documents. The background chain reads only the
% MEAN closed service time, which is exact under PS by the same insensitivity
% and a first-moment surrogate under FCFS. The level space is truncated at
% options.cutoff, which costs nothing at the tail probabilities the default
% chooses.
%
% OPTIONS.CONFIG.BGENV selects what the QBD carries on its environment axis:
%
%   'lump' (default) : the background chain aggregated onto the number of closed
%                      jobs THIS station holds (MAM_BGCHAIN_ENV)
%   'full'           : the background chain itself (MAM_BGCHAIN_ENVFULL), at a
%                      phase count bounded by OPTIONS.CONFIG.QBDPHASES_MAX
%
% 'full' IS AN ORACLE AND RETURNS THE SAME NUMBERS. The background chain is
% product-form by construction, so the lump is Norton-exact rather than
% approximate; the two agree to 6.7e-16 on every model measured. It is MATLAB
% ONLY, and it is here to re-check that property, not to improve on it. See
% MAM_BGCHAIN_ENVFULL for the argument and the measurements.
%
% THE NAME IS 'full', NOT 'exact', AND DELIBERATELY SO. It says what the QBD
% carries on its environment axis, and NOTHING about the answer: this method
% stays approximate under it, for the reasons listed above -- the mean-only
% closed service time, the mean-field coupling to the other stations, the FCFS
% random-order surrogate, the multiserver collapse, the level truncation and,
% at G < R-1, the Chandy-Herzog-Woo aggregation. Measured with 'full' selected,
% a two-chain FCFS model with class-dependent rates is still 83% out on a queue
% length against SolverCTMC.
%
% See also SOLVER_MAM_ANALYZER, MAM_BGCHAIN_CTMC, MAM_BGCHAIN_ENV,
% MAM_BGCHAIN_ENVFULL, MAM_BGCHAIN_STATION.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
K = sn.nclasses;
C = sn.nchains;
N = sn.njobs';
PH = sn.proc;
tol = options.tol;

if ~isfield(options, 'config') || isempty(options.config)
    options.config = struct();
end
if ~isfield(options.config, 'space_max')
    options.config.space_max = 128;
end
if ~isfield(options.config, 'qbdphases_max')
    options.config.qbdphases_max = 500;
end
if ~isfield(options.config, 'bgaggr')
    options.config.bgaggr = 1;
end
if ~isfield(options.config, 'bgenv') || isempty(options.config.bgenv)
    options.config.bgenv = 'lump';
end
options.config.bgenv = lower(char(options.config.bgenv));
if ~any(strcmp(options.config.bgenv, {'lump', 'full'}))
    line_error(mfilename, sprintf(['Unknown options.config.bgenv ''%s''. Use ''lump'' (the ' ...
        'default: the background chain aggregated onto the closed occupancy of the station) or ' ...
        '''full'' (the whole chain as the QBD environment, no lumping).'], options.config.bgenv));
end

% The model-class rules live in MAM_BGCHAIN_APPLICABLE, which
% SolverMAM.supportsModelMethod and the default chooser ask as well.
[bgOk, bgWhy] = mam_bgchain_applicable(sn, options);
if ~bgOk
    line_error(mfilename, bgWhy);
end

S = 1 ./ sn.rates;
S(~isfinite(S)) = 0;

[rtst, V] = sn_rt_stations(sn);
[Lchain, STchain, Vchain, alpha, Nchain] = sn_get_demands_chain(sn);

isopenchain = false(1, C);
for c = 1:C
    isopenchain(c) = any(isinf(N(sn.inchain{c})));
end
openChains = find(isopenchain);
closedChains = find(~isopenchain & Nchain > 0);
R = numel(closedChains); % R >= 1: MAM_BGCHAIN_APPLICABLE refused a purely open model above

%% Stations the closed chains visit: the support of the background chain
inbg = false(1, M);
for c = closedChains
    inbg = inbg | (Vchain(:, c)' > GlobalConstants.Zero);
end
cst = find(inbg);
Mc = numel(cst); % Mc >= 1, by the same predicate

%% Chain-level station routing, obtained by folding the class axis of sn.rt
Pchain = cell(1, C);
for c = 1:C
    inchain = sn.inchain{c}(:)';
    P = zeros(M, M);
    for i = 1:M
        for k = inchain
            if alpha(i, k) <= 0
                continue;
            end
            row = rtst((i-1)*K + k, :);
            acc = zeros(1, M);
            for kp = inchain
                acc = acc + row((0:M-1)*K + kp);
            end
            P(i, :) = P(i, :) + alpha(i, k) * acc;
        end
    end
    Pchain{c} = P;
end

%% Open arrival streams
lambdaChain = zeros(1, C);
chainArrival = cell(1, C);      % {D0,D1} of the chain source process
for c = openChains
    inchain = sn.inchain{c}(:)';
    isrc = sn.refstat(inchain(1));
    rates_c = sn.rates(isrc, inchain);
    lambdaChain(c) = sum(rates_c(isfinite(rates_c)));
    acc = [];
    for k = inchain
        if ~isfinite(sn.rates(isrc, k)) || sn.rates(isrc, k) <= 0
            continue;
        end
        Mk = PH{isrc}{k};
        if isempty(Mk) || any(isnan(Mk{1}(:)))
            continue;
        end
        if isempty(acc)
            acc = {Mk{1}, Mk{2}, Mk{2}};
        else
            acc = mmap_super_safe({acc, {Mk{1}, Mk{2}, Mk{2}}}, options.config.space_max, 'default');
        end
    end
    if isempty(acc)
        acc = mmap_exponential(GlobalConstants.Zero, 1);
    end
    chainArrival{c} = {acc{1}, acc{2}};
end

% per-class open arrival rate at each station
lambdaOpen = zeros(M, K);
for c = openChains
    for k = sn.inchain{c}(:)'
        lambdaOpen(:, k) = lambdaChain(c) * V(:, k);
    end
end
isopenclass = false(1, K);
for c = openChains
    isopenclass(sn.inchain{c}) = true;
end

%% Fixed-point state
QN = zeros(M, K); UN = zeros(M, K); RN = zeros(M, K);
TN = zeros(M, K); CN = zeros(1, K); XN = zeros(1, K);
% cshare(i,e+1) is the mean number of servers of station i that its e closed
% jobs hold once the open work has taken its share: the one quantity the two
% halves exchange. It starts at min(e,c), the value with no open work at all.
Ntot = sum(Nchain(closedChains));
cshare = zeros(M, Ntot + 1);
for i = 1:M
    cshare(i, :) = min(0:Ntot, sn.nservers(i));
end
Xclosed = zeros(1, C);
for c = closedChains
    denom = sum(Lchain(:, c));
    if denom > 0
        Xclosed(c) = Nchain(c) / denom;      % Bard-Schweitzer style lower bound
    end
end

%% How many aggregate classes the background chain carries, and which chains
% share each of them. options.config.bgaggr is the number of AGGREGATE classes
% G: G = 1 is the classic tagged/aggregate pair, G = R-1 leaves every untagged
% chain on its own and aggregates nothing. Anything above R-1 (R itself, say)
% clamps there, so asking for "no aggregation" needs no magic value.
naggr = min(max(round(options.config.bgaggr), 1), max(R - 1, 1));
% With nothing left to aggregate, ONE background chain carries every closed
% chain exactly and the tagged loop would solve the same chain R times over.
noAggr = (R == 1) || (naggr >= R - 1);

% The grouping is a property of the demands, not of the iterate, so it is fixed
% once here rather than recomputed inside the fixed point.
grpOf = cell(1, R);
othersOf = cell(1, R);
if ~noAggr
    for ridx = 1:R
        others = setdiff(closedChains, closedChains(ridx));
        othersOf{ridx} = others;
        grpOf{ridx} = mam_bgchain_groups(Lchain(cst, others), naggr);
    end
end

if noAggr
    passes = closedChains(1);
else
    passes = closedChains;
end
npass = numel(passes);

TN_1 = Inf(M, K);
totiter = 0;
relax = 0.5;

while max(max(abs(TN - TN_1))) > tol && totiter < options.iter_max
    totiter = totiter + 1;
    TN_1 = TN;

    QopenAcc = zeros(M, 1);
    UopenAcc = zeros(M, 1);
    cshareAcc = zeros(M, Ntot + 1);

    for pidx = 1:npass
        r = passes(pidx);
        ridx = find(closedChains == r, 1);

        %% Background classes: class 1 is the tagged chain, classes 2..1+G are
        % the flow-equivalent aggregates of the demand-similar groups. With no
        % aggregation every closed chain is a class of its own, in chain order.
        if noAggr
            members = num2cell(closedChains);
        else
            others = othersOf{ridx};
            grp = grpOf{ridx};
            members = cell(1, 1 + naggr);
            members{1} = r;
            for g = 1:naggr
                members{1 + g} = others(grp == g);
            end
        end
        B = numel(members);
        Nb = zeros(1, B);
        STb = zeros(Mc, B);
        Pb = cell(1, B);
        % A class can only hold jobs at the stations its members visit; see
        % MAM_BGCHAIN_CTMC on why enumerating the union instead is not merely
        % wasteful but makes the generator reducible.
        suppb = false(Mc, B);
        for b = 1:B
            mem = members{b};
            Nb(b) = sum(Nchain(mem));
            for oi = 1:numel(mem)
                suppb(:, b) = suppb(:, b) | (Vchain(cst, mem(oi)) > GlobalConstants.Zero);
            end
            if isscalar(mem)
                % a group of one is carried exactly: no mean to take
                STb(:, b) = STchain(cst, mem);
                Pb{b} = Pchain{mem}(cst, cst);
            else
                w = zeros(Mc, numel(mem));
                for oi = 1:numel(mem)
                    w(:, oi) = Xclosed(mem(oi)) * Vchain(cst, mem(oi));
                end
                rowsum = sum(w, 2);
                for ii = 1:Mc
                    if rowsum(ii) > 0
                        w(ii, :) = w(ii, :) / rowsum(ii);
                    else
                        w(ii, :) = 1 / numel(mem);
                    end
                end
                STb(:, b) = sum(w .* STchain(cst, mem), 2);
                Pagg = zeros(Mc, Mc);
                for oi = 1:numel(mem)
                    Pagg = Pagg + w(:, oi) .* Pchain{mem(oi)}(cst, cst);
                end
                Pb{b} = Pagg;
            end
            Pb{b} = bgchain_rownorm(Pb{b});
        end

        bg = mam_bgchain_ctmc(Nb, STb, Pb, sn.sched(cst), sn.nservers(cst), cshare(cst, :), ...
            suppb, options);

        %% Closed-class metrics of every chain this pass carries EXACTLY: the
        % tagged one always, and every chain when nothing was aggregated.
        if noAggr
            bext = 1:B;
        else
            bext = 1;
        end
        for b = bext
            rb = members{b}(1);
            inchain = sn.inchain{rb}(:)';
            for k = inchain
                QN(:, k) = 0; UN(:, k) = 0; RN(:, k) = 0; TN(:, k) = 0;
            end
            for ii = 1:Mc
                i = cst(ii);
                for k = inchain
                    if alpha(i, k) <= 0
                        continue;
                    end
                    % THROUGHPUT splits by VISIT share, OCCUPANCY by DEMAND
                    % share. A chain queue divided by alpha alone gives every
                    % class of a station the same response time, which is
                    % impossible at a Delay, where R must be the class service
                    % time; the weight is alpha*ST/STchain, the rule
                    % sn_deaggregate_chain_results applies. On
                    % cqn_twoclass_hyperl the visit-share split was 27% off the
                    % exact product-form QLen while the chain total was exact.
                    if STchain(i, rb) > GlobalConstants.Zero
                        w = alpha(i, k) * S(i, k) / STchain(i, rb);
                    else
                        w = alpha(i, k);
                    end
                    QN(i, k) = bg.QLen(ii, b) * w;
                    TN(i, k) = bg.Tput(ii, b) * alpha(i, k);
                    if sn.sched(i) == SchedStrategy.INF
                        UN(i, k) = QN(i, k);
                    else
                        UN(i, k) = bg.Ubusy(ii, b) * w;
                    end
                    if TN(i, k) > GlobalConstants.Zero
                        RN(i, k) = QN(i, k) / TN(i, k);
                    end
                end
            end
            iref = sn.refstat(inchain(1));
            vref = sum(Vchain(iref, rb));
            if vref > GlobalConstants.Zero
                Xclosed(rb) = sum(TN(iref, inchain)) / vref;
            else
                Xclosed(rb) = sum(TN(iref, inchain));
            end
        end

        %% Open classes: one modulated QBD per queueing station
        Uclosed = zeros(M, 1);
        Uclosed(cst) = sum(bg.Ubusy, 2);
        [Qo, Uo, gc] = bgchain_open_pass(sn, options, bg, cst, S, V, PH, ...
            lambdaOpen, lambdaChain, chainArrival, openChains, Uclosed, cshare);
        QopenAcc = QopenAcc + Qo;
        UopenAcc = UopenAcc + Uo;
        cshareAcc = cshareAcc + gc;
    end

    Qopen = QopenAcc / npass;
    Uopen = UopenAcc / npass; %#ok<NASGU> reported through the utilization law below
    cshare = (1 - relax) * cshare + relax * (cshareAcc / npass);

    %% Open-class metrics from the aggregate station results
    for i = 1:M
        kopen = find(isopenclass & lambdaOpen(i, :) > GlobalConstants.Zero);
        if isempty(kopen)
            for k = find(isopenclass)
                TN(i, k) = lambdaOpen(i, k);
                QN(i, k) = 0; UN(i, k) = 0; RN(i, k) = 0;
            end
            continue;
        end
        lam = lambdaOpen(i, kopen);
        lamtot = sum(lam);
        Smix = sum(lam .* S(i, kopen)) / lamtot;
        switch sn.sched(i)
            case SchedStrategy.EXT
                for k = kopen
                    TN(i, k) = lambdaOpen(i, k);
                    QN(i, k) = 0; UN(i, k) = 0; RN(i, k) = 0;
                end
            case SchedStrategy.INF
                for k = kopen
                    TN(i, k) = lambdaOpen(i, k);
                    RN(i, k) = S(i, k);
                    QN(i, k) = TN(i, k) * S(i, k);
                    UN(i, k) = QN(i, k);
                end
            otherwise
                Rtot = Qopen(i) / lamtot;
                for k = kopen
                    TN(i, k) = lambdaOpen(i, k);
                    if sn.sched(i) == SchedStrategy.PS
                        % processor sharing: residence scales with the demand
                        RN(i, k) = Rtot * S(i, k) / Smix;
                    else
                        % FCFS and its variants: the wait is class-blind, the
                        % service time is not
                        RN(i, k) = max(S(i, k), Rtot - Smix + S(i, k));
                    end
                    QN(i, k) = TN(i, k) * RN(i, k);
                    % Utilization Law: a c-server station holds TN*S/c of its
                    % capacity, the same identity solver_mam_basic reports.
                    UN(i, k) = TN(i, k) * S(i, k) / sn.nservers(i);
                end
        end
    end
end

%% System measures
for c = 1:C
    inchain = sn.inchain{c}(:)';
    if isopenchain(c)
        XN(inchain) = lambdaChain(c);
    else
        XN(inchain) = Xclosed(c);
    end
end
CN = sum(RN, 1);

QN(~isfinite(QN)) = 0;
UN(~isfinite(UN)) = 0;
RN(~isfinite(RN)) = 0;
TN(~isfinite(TN)) = 0;
CN(~isfinite(CN)) = 0;
XN(~isfinite(XN)) = 0;
end

function P = bgchain_rownorm(P)
% Row-normalizes a routing matrix, leaving an all-zero row as a self-loop so
% the background chain stays a proper Markov chain on its support.
for i = 1:size(P, 1)
    s = sum(P(i, :));
    if s > GlobalConstants.Zero
        P(i, :) = P(i, :) / s;
    else
        P(i, :) = 0;
        P(i, i) = 1;
    end
end
end

function [Qopen, Uopen, gnew] = bgchain_open_pass(sn, options, bg, cst, S, V, PH, ...
    lambdaOpen, lambdaChain, chainArrival, openChains, Uclosed, gref)
% One pass of the open side: for each station, build the environment the closed
% jobs present to it -- the background chain lumped onto its closed occupancy,
% or the whole chain under options.config.bgenv='full' -- and solve the
% resulting modulated level-dependent QBD.
M = sn.nstations;
Qopen = zeros(M, 1);
Uopen = zeros(M, 1);
gnew = gref;   % a station with no open queue keeps the share it had

for i = 1:M
    if sn.sched(i) == SchedStrategy.EXT || sn.sched(i) == SchedStrategy.INF
        continue;
    end
    kopen = find(lambdaOpen(i, :) > GlobalConstants.Zero);
    if isempty(kopen)
        continue;
    end

    %% Aggregate open arrival process at the station
    Da = [];
    for c = openChains
        inchain = sn.inchain{c}(:)';
        rate_ic = lambdaChain(c) * sum(V(i, inchain));
        if rate_ic <= GlobalConstants.Zero || lambdaChain(c) <= GlobalConstants.Zero
            continue;
        end
        Mc_arr = map_scale(chainArrival{c}, 1 / rate_ic);
        if isempty(Da)
            Da = Mc_arr;
        else
            sup = mmap_super_safe({{Da{1}, Da{2}, Da{2}}, {Mc_arr{1}, Mc_arr{2}, Mc_arr{2}}}, ...
                options.config.space_max, 'default');
            Da = {sup{1}, sup{2}};
        end
    end
    if isempty(Da)
        continue;
    end

    %% Arrival-weighted phase-type mixture of the open service laws
    lam = lambdaOpen(i, kopen);
    w = lam / sum(lam);
    alphas = [];
    Tblk = [];
    isPSstation = (sn.sched(i) == SchedStrategy.PS);
    for kk = 1:numel(kopen)
        k = kopen(kk);
        PHk = PH{i}{k};
        if isPSstation || isempty(PHk) || any(isnan(PHk{1}(:)))
            % PROCESSOR SHARING IS INSENSITIVE to the service law beyond its
            % mean, so carrying the phase-type representation at a PS station is
            % not merely unnecessary, it is WRONG. This QBD tracks ONE service
            % phase for the whole station, which makes the open queue length
            % inherit the SCV-sensitivity of an M/PH/1 FCFS queue; measured
            % against SolverCTMC, a HyperExp of SCV 4 then read 21% high where
            % the exact answer is the exponential one to five digits. The
            % exponential of the same mean is exact here, and it shrinks the
            % QBD's phase count as a side effect.
            PHk = map_exponential(max(S(i, k), GlobalConstants.Zero));
        else
            PHk = map_scale(PHk, S(i, k));
        end
        pik = map_pie(PHk);
        alphas = [alphas, w(kk) * pik(:)']; %#ok<AGROW>
        Tblk = blkdiag(Tblk, PHk{1});
    end

    %% Environment seen by the open classes at this station
    if any(cst == i)
        ii = find(cst == i, 1);
        switch options.config.bgenv
            case 'full'
                [A, ~, esup] = mam_bgchain_envfull(bg, ii);
            otherwise
                [A, ~, esup] = mam_bgchain_env(bg, ii);
        end
    else
        A = 0; esup = 0;
    end

    nphases = size(Da{1}, 1) * numel(esup) * numel(alphas);
    if nphases > options.config.qbdphases_max
        if strcmp(options.config.bgenv, 'full')
            envtxt = ['axis is the WHOLE background chain (options.config.bgenv=''full''), so it ' ...
                'grows as the chain does. Set options.config.bgenv=''lump'' to aggregate it onto ' ...
                'the closed occupancy of the station, raise options.config.qbdphases_max'];
        else
            envtxt = ['axis is the closed population held by the station, so it grows with the ' ...
                'closed population. Raise options.config.qbdphases_max'];
        end
        line_error(mfilename, sprintf(['The modulated QBD of station %d needs %d phases ' ...
            '(%d arrival x %d environment x %d service), above the limit of %d. The environment ' ...
            '%s, or reduce the closed population or the order of the arrival and service ' ...
            'processes.'], i, nphases, size(Da{1}, 1), numel(esup), numel(alphas), ...
            options.config.qbdphases_max, envtxt));
    end

    Kmax = bgchain_cutoff(options, sum(lam), sum(lam .* S(i, kopen)) / sum(lam), ...
        sn.nservers(i), Uclosed(i));
    ecap = min(esup, size(gref, 2) - 1);
    res = mam_bgchain_station(Da{1}, Da{2}, alphas, Tblk, A, esup, sn.nservers(i), ...
        gref(i, ecap + 1), Kmax, options);

    Qopen(i) = res.QLen;
    Uopen(i) = res.Util;
    % The QBD only saw the environment states the background chain reaches;
    % interpolate the rest so the next background chain has a share wherever it
    % may go, clipped to the capacity a closed job can physically hold.
    egrid = 0:(size(gref, 2) - 1);
    if numel(res.esup) > 1
        g = interp1(res.esup, res.cshare, egrid, 'linear', 'extrap');
    else
        g = res.cshare * ones(1, numel(egrid));
    end
    gnew(i, :) = min(max(g, 0), min(egrid, sn.nservers(i)));
end
end

function Kmax = bgchain_cutoff(options, lambda, Smix, nservers, Uclosed)
% Truncation level of the open queue: the explicit cutoff when the user set
% one, else enough levels for the geometric tail left by the closed traffic to
% be negligible.
if isfield(options, 'cutoff') && isscalar(options.cutoff) && isfinite(options.cutoff) ...
        && options.cutoff > 0
    Kmax = max(2, round(options.cutoff));
    return;
end
free = max(GlobalConstants.FineTol, 1 - Uclosed);
rho = lambda * Smix / (nservers * free);
rho = min(max(rho, 1e-3), 1 - 1e-3);
Kmax = ceil(log(1e-8) / log(rho));
Kmax = min(max(Kmax, 20), 200);
end
