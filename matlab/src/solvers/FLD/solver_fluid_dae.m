function [QN, UN, RN, TN, xvec_it, QNt, UNt, TNt, xvec_t, t, iters, runtime, moments] = solver_fluid_dae(sn, options)
% [QN, UN, RN, TN, XVEC_IT, QNT, UNT, TNT, XVEC_T, T, ITERS, RUNTIME, MOMENTS] = SOLVER_FLUID_DAE(SN, OPTIONS)
%
% Differential-algebraic formulation of the min-normal closure, backing
% options.method='dae'.
%
% SOLVER_FLUID_MOMENTS already solves a differential system (the mean)
% coupled to an algebraic one (the covariance). It solves them by SUCCESSIVE
% SUBSTITUTION: integrate the mean ODE to its fixed point at a held variance,
% solve the Lyapunov equation there, extract sigma2, repeat, up to 20 times
% and only to CoarseTol. This method states the same closure as one system
% and solves it as one system. Nothing about the closure changes -- the drift,
% the rate factors and the Lyapunov equation are the same functions, taken
% unmodified from FLUID_MOMENT_TERMS and FLUID_LYAPUNOV -- only the way the
% coupled equations are discharged.
%
% Two modes, chosen by the horizon, because the DAE is genuinely different in
% each:
%
%   STEADY STATE (options.timespan(2) infinite, the usual case)
%       Solve the algebraic system
%
%           0 = D*r(x, sigma2)              drift residual, nstate rows
%           0 = C*x - Nchain                population conservation, one row
%                                           per closed chain
%           0 = sigma2 - sigmaOf(x,sigma2)  closure consistency, one row per
%                                           closable station
%
%       simultaneously by damped Newton. Convergence is quadratic near the
%       root instead of the linear rate of substitution, and the answer is
%       converged to options.tol rather than to the outer CoarseTol.
%
%       THE FIXED POINT IS NOT FOUND BY INTEGRATING TO IT. Integrating a
%       stable ODE until it stops moving is a poor way to solve f(x)=0: the
%       cost is set by the slowest mode of the model rather than by the
%       accuracy wanted, which is exactly why the stiff models are expensive
%       here. A seed trajectory is still integrated once, cheaply and at the
%       first-order closure, because Newton needs a starting point inside the
%       basin; everything after that is algebraic.
%
%   TRANSIENT (finite horizon)
%       Integrate the index-1 DAE
%
%           d/dt x     = D*r(x, sigma2(t))        differential, kept rows
%           0          = C*x - Nchain             algebraic, one row per chain
%           d/dt Sigma = A*Sigma + Sigma*A' + Q   differential
%
%       with a singular mass matrix. This is the ONLY fluid method other than
%       'kp' that produces a time-varying second moment, and unlike
%       'minnormal' -- which evaluates the whole transient at the single
%       STATIONARY variance -- the variance here is the one the trajectory
%       actually had at each instant.
%
% WHAT THE ALGEBRAIC CONSTRAINT BUYS. Population conservation currently holds
% only to integrator tolerance: it is a consequence of the drift (the rows of
% D sum to zero on a closed chain), never an equation. Writing it as a
% constraint has two effects. It is enforced to solver tolerance rather than
% accumulated. And it is what makes the Newton system solvable at all: the
% drift Jacobian is singular along exactly the conserved directions -- the
% same singularity FLUID_LYAPUNOV works around by projecting onto range(D) --
% so the constraint rows supply the missing rank instead of a pseudo-inverse
% hiding it.
%
% WHY THE COVARIANCE IS NOT A NEWTON UNKNOWN. The obvious "solve everything
% at once" reading puts Sigma in the unknown vector. That is a trap: Sigma is
% nstate^2 entries, so the Jacobian is (nstate + nstate^2)^2 and the work is
% quartic -- strictly worse than the cubic Lyapunov solves it would replace,
% and it would tighten the state cap rather than relax it. Sigma is LINEAR in
% itself for a held x, so it is eliminated by one Lyapunov solve per residual
% evaluation and only sigma2, M numbers, joins x in the unknown vector.
%
% NO SIGMA2=0 SEED, SO NO KINK WORKAROUND. SOLVER_FLUID_MOMENTS must start its
% alternation at sigma2=0, where min(n,c) has no derivative, and a saturated
% model's first-order fixed point lands on that kink by construction; it
% carries a two-sided-Jacobian probe to decide hyperbolicity side-independently
% there. The simultaneous solve never adopts sigma2=0 as an iterate, so that
% probe is not needed. It does NOT rescue a fixed point that genuinely sits on
% the kink at the converged variance -- that is a real continuum of equilibria
% and no formulation removes it.
%
% Parameters:
%   sn      - NetworkStruct
%   options - solver options; options.tol sets the Newton tolerance,
%             options.config.dae_maxstate caps the simultaneous solve
%             (default 100) and options.config.dae_maxcov caps the transient
%             covariance (default 25)
%
% Returns:
%   as SOLVER_FLUID_MOMENTS, plus MOMENTS.Sigmat / MOMENTS.QVart / MOMENTS.tvar
%   on the transient route
%
% See also SOLVER_FLUID_MOMENTS, FLUID_MOMENT_TERMS, FLUID_LYAPUNOV.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
M = sn.nstations;
K = sn.nclasses;

terms = fluid_moment_terms(sn, options);
nstate = terms.nstate;

% The simultaneous solve carries one Lyapunov solve per residual evaluation
% and takes a finite-difference Jacobian over nstate + nclosable unknowns, so
% its cost is cubic per evaluation and quartic overall. That is affordable at
% the scale the moment methods already run at, but the crossover is lower than
% the 200 that SOLVER_FLUID_MOMENTS permits, so it gets its own limit.
maxstate = 100;
if isfield(options,'config') && isfield(options.config,'dae_maxstate') && ~isempty(options.config.dae_maxstate)
    maxstate = options.config.dae_maxstate;
end
if nstate > maxstate
    line_error(mfilename, sprintf(['The dae method solves a %d-unknown algebraic system with a ' ...
        'finite-difference Jacobian, above the limit of %d set by options.config.dae_maxstate. ' ...
        'Raise that limit, or use options.method=''minnormal'' for the same closure by ' ...
        'successive substitution.'], nstate, maxstate));
end

% The DPS and GPS shares close on the covariance BETWEEN station coordinates,
% not on the station total, so their closure state is a matrix block rather
% than the scalar sigma2. Carrying those blocks as Newton unknowns would put
% the quartic term back; carrying them as a chord would reintroduce the
% alternation this method exists to remove. Refuse instead of doing either
% quietly.
if any(terms.sched == SchedStrategy.DPS | terms.sched == SchedStrategy.GPS)
    line_error(mfilename, ['The dae method closes on the per-station variance only, but this model ' ...
        'has a DPS or GPS station whose share closes on the covariance BETWEEN its class ' ...
        'coordinates. Use options.method=''minnormal'', which carries those blocks through its ' ...
        'outer iteration.']);
end

% -- capacity limits, as admission constraints -------------------------------
% SolverFLD refuses these outright because the drift cannot express them. The
% DAE form can: a cap is a linear inequality on the state, and blocking is a
% throttle on the admission flow that keeps it satisfied. Both declarations are
% one family of rows -- a finite capacity region over a SET of stations, and a
% station's own buffer -- and what differs is where the blocked job goes, which
% is what GATES decides, per cap and per event.
con = fluid_capacity_constraints(sn, terms);
ncon = numel(con.b);
gates = fluid_capacity_gates(sn, terms, con);

% -- the waiting room outside each capped region ----------------------------
% Built before the conservation rows, which have to count the blocked mass. Only
% a STAGED cap has a room: a station buffer holds the job at the upstream station
% instead, where LINE's own semantics leaves it.
stg = fluid_capacity_staging(terms, con);
con = fluid_capacity_extend(con, stg, terms);

% -- population conservation, as equations ----------------------------------
% One row per closed chain. C*D is zero on those rows whenever the event set
% conserves the chain, which is the structural fact that makes the drift
% Jacobian singular; the check below is cheap and turns a malformed event set
% into a named error rather than a silently rank-deficient Newton system.
[C, Nvec] = local_conservation(sn, terms, stg);
if ~isempty(C)
    leak = max(abs(C(:,1:nstate) * terms.D), [], 'all');
    if leak > sqrt(GlobalConstants.Zero)
        line_error(mfilename, sprintf(['The event set does not conserve a closed chain: the largest ' ...
            'population leak per unit rate is %g, where it must be zero. The conservation constraint ' ...
            'would contradict the drift rather than complete it.'], leak));
    end
end

% -- which stations have a variance that enters the drift --------------------
% A delay or a source has no min() to close. A station that cannot fill its
% servers has min(n,c) = n on its whole support, so closing it is not an
% improvement but an error (see FLUID_MOMENT_TERMS). Those stations are held
% at zero and are not unknowns, which keeps the Newton system as small as the
% closure actually is.
closable = true(M,1);
for i = 1:M
    if terms.isExt(i) || terms.sched(i) == SchedStrategy.INF || ~isfinite(sn.nservers(i)) ...
            || terms.minExact(i) || isempty(terms.stationBlock{i})
        closable(i) = false;
    end
end
cidx = find(closable);
ncl = numel(cidx);
covblk = cell(M,1); % always empty here: DPS/GPS are refused above

% -- seed --------------------------------------------------------------------
% Newton needs a point in the basin, not an answer. One first-order solve
% supplies it, at the cost of the single integration this method is trying to
% avoid repeating twenty times.
seedopt = options;
seedopt.method = 'closing';
seedopt.config.moment_sigma2 = zeros(M,1);
seedopt.config.moment_cov = cell(M,1);
[~, xvec_seed, ~, ~, xvec_t, t, seed_iters] = solver_fluid(sn, seedopt);
x = xvec_seed{end}(:);

% THE VARIANCE IS SEEDED POSITIVE, which is the whole reason this route has no
% kink probe. sigma2 = 0 is where min(n,c) has no derivative, and a saturated
% model's first-order fixed point sits exactly there, so starting Newton at
% zero would reintroduce the degeneracy SOLVER_FLUID_MOMENTS spends two
% one-sided Jacobians ruling out. The station mean is an O(N) starting value in
% the right units -- it is the variance a Poisson population of that mean would
% have -- and it costs no Lyapunov solve to form.
sigma2 = zeros(M,1);
for j = 1:numel(cidx)
    i = cidx(j);
    blk = terms.stationBlock{i};
    if ~isempty(blk)
        sigma2(i) = max(GlobalConstants.FineTol, sum(x(blk)));
    end
end

% -- steady state: one simultaneous algebraic solve --------------------------
tol = options.tol;
if isempty(tol) || ~isfinite(tol)
    tol = GlobalConstants.FineTol;
end
newton_max = 50;
if isfield(options,'iter_max') && ~isempty(options.iter_max)
    newton_max = max(newton_max, options.iter_max);
end

% ACTIVE SET. A capacity constraint is an inequality, and an inequality has no
% residual to hand a Newton solver -- only the constraints that actually bind
% become equations. Each pass is a complete simultaneous solve, so the loop
% iterates over the COMBINATORICS of which caps bind, not over the closure.
%
% ONE MULTIPLIER PER ACTIVE ROW, not per region. A region's waiting room used to
% carry a single drain rate, so two caps of one region were two equalities
% against one control and the case was refused by name. A room gated by several
% active rows now drains at the harmonic composition of their rates and a held
% cap composes as a product of fractions, so each row keeps a control of its own;
% see FLUID_CAPACITY_LEGS.
%
% THE SET STARTS EMPTY, and the first pass therefore asks for the UNCONSTRAINED
% fixed point -- which is what makes the answer of a model whose caps bind
% independent of how far outside the seed happened to land. It is also not always
% solvable: an overloaded M/M/1/K has NO equilibrium without its cap, so that pass
% fails rather than converging and the caps the iterate violates are seeded into
% the set instead (below). Doing that up front for every model looked cheaper and
% changed answers: a region of 8 over two identical queues settled 7.43/0.57 from
% the clipped seed and 4/4 from the unconstrained one, and only the second is what
% every other codebase reports.
active = zeros(0,1);
seededFromFailure = false;
activeSolved = zeros(0,1);
iters = 0;
mult = zeros(0,1);
sg = zeros(stg.n,1);
aset_max = max(4, 2*ncon + 2);
converged = false; resnorm = Inf; clamped = false;
% the last iterate that WAS a fixed point, and the whole answer read off it. A
% pass that fails to converge leaves an iterate that is not a fixed point of
% anything, and reading the next active set off it is how a single bad pass turned
% into a walk through unrelated capacity combinations.
xOk = x; sigma2Ok = sigma2;
best = [];
for aset = 1:aset_max
    hasClamp = any(~con.staged(active));
    % Seed each waiting room with the mass that does not fit, and every
    % multiplier at unity -- an unthrottled fraction for a held or lost cap, and
    % the drain rate the region route has always started from.
    sg0 = zeros(stg.n,1);
    for k = 1:numel(active)
        c = active(k);
        if ~con.staged(c)
            continue
        end
        excess = max(0, con.A(c,:)*x - con.b(c));
        sel = find(stg.gatedBy(c,:));
        if ~isempty(sel)
            sg0(sel) = excess / numel(sel);
        end
    end
    u0 = [x; sg0; sigma2(cidx); ones(numel(active),1)];
    % THE CLAMPED COVARIANCE IS A FALLBACK, NOT THE DEFAULT. A cap that holds or
    % loses fixes its own combination of the state, so the honest linear-noise
    % approximation puts no fluctuation there (LOCAL_CLAMP_TANGENT) -- but the
    % truth is neither that nor the unprojected variance: the population under a
    % cap follows a TRUNCATED distribution, whose variance is smaller than the
    % unconstrained one and larger than zero. The unprojected solve is what every
    % other fluid method computes, so it is what runs first and what keeps the
    % answers of the models that already worked; the projection is tried only when
    % the unprojected system has no stationary covariance at all, which is exactly
    % the neutral case an overloaded loss station produces. Deciding once per PASS
    % rather than per residual matters: a projection that switched on and off
    % between iterates would give Newton a discontinuous system to converge on.
    if hasClamp
        clampAttempts = [false true];
    else
        clampAttempts = false;
    end
    u = u0; nit = 0;
    nohyp = [];
    for ia = 1:numel(clampAttempts)
        useClamp = clampAttempts(ia);
        ctx = struct('terms', terms, 'C', C, 'Nvec', Nvec, 'cidx', cidx, 'covblk', {covblk}, ...
            'nstate', nstate, 'conA', con.A, 'conAs', con.As, 'conb', con.b, ...
            'conStaged', con.staged, 'gates', gates, 'stg', stg, 'active', active, ...
            'clampT', []);
        if useClamp
            ctx.clampT = local_clamp_tangent(con, active, terms);
        end
        resid = @(uu, quiet) local_residual(uu, ctx, quiet);
        try
            % The first NSTATE unknowns (the means) are free and every one after
            % them -- variance, throttle -- is non-negative. Stated as a bound
            % VECTOR rather than as the count: a scalar bound is indistinguishable
            % from a one-unknown bound vector, and FLUID_DAE_PROJECT cannot tell
            % which was meant.
            lbDae = [-inf(nstate,1); zeros(numel(u0)-nstate,1)];
            [u, nit, converged, resnorm] = fluid_dae_newton(resid, u0, tol, newton_max, lbDae);
        catch me
            if ~strcmp(me.identifier,'LINE:FluidNonHyperbolic')
                rethrow(me);
            end
            if useClamp || ~hasClamp
                nohyp = me;
                break
            end
            continue
        end
        iters = iters + nit;
        clamped = useClamp;
        if converged
            break
        end
    end
    if ~isempty(nohyp)
        % THE UNCONSTRAINED FIXED POINT NEED NOT EXIST. An overloaded open station
        % has no equilibrium until its buffer bounds it, so the pass that asks for
        % one fails and the caps the iterate violates -- or left non-finite, which
        % no comparison catches -- are seeded into the active set instead. Once: a
        % second failure with caps already bound is the model's answer rather than a
        % starting point to improve.
        cand = zeros(0,1);
        for c = 1:ncon
            if any(active == c)
                continue
            end
            sel = con.A(c,:) > 0;
            if any(~isfinite(x(sel))) || con.A(c,:)*x > con.b(c) + max(1e-9, tol)
                cand(end+1,1) = c; %#ok<AGROW>
            end
        end
        if seededFromFailure || isempty(cand)
            rethrow(nohyp);
        end
        for k = 1:numel(cand)
            c = cand(k);
            sel = con.A(c,:) > 0;
            npos = sum(sel);
            if npos == 0 || ~(con.b(c) > 0)
                continue
            end
            val = con.A(c,:)*x;
            if isfinite(val) && val > con.b(c)
                x(sel) = x(sel) * (con.b(c)/val);
            elseif ~isfinite(val)
                x(sel) = con.b(c)/npos;
            end
        end
        x(~isfinite(x)) = 0;
        active = unique([active; cand]);
        seededFromFailure = true;
        continue
    end
    x = u(1:nstate);
    sg = u(nstate+(1:stg.n));
    sigma2 = zeros(M,1);
    sigma2(cidx) = max(0, u(nstate+stg.n+(1:numel(cidx))));
    mult = max(0, u(nstate+stg.n+numel(cidx)+1:end));
    % the set that produced THIS x and these multipliers, which is what the
    % metrics are read at; `active` below is the set to try NEXT
    activeSolved = active;
    if ncon == 0
        break
    end

    slack = con.b - con.A*x;
    if stg.n > 0
        slack = slack - con.As*sg;
    end
    if converged
        xOk = x; sigma2Ok = sigma2;
        best = struct('x', x, 'sg', sg, 'sigma2', sigma2, 'mult', mult, ...
            'active', activeSolved, 'resnorm', resnorm, 'clamped', clamped, ...
            'feasible', all(slack >= -max(1e-9, tol)));
    end
    violated = find(slack < -max(1e-9, tol) & ~ismember((1:ncon)', active));
    % THE RELEASE SIGNAL IS THE MULTIPLIER'S OWN UNITS, and the two kinds do not
    % share them. A held or lost cap throttles by a FRACTION, so a fraction that
    % came back above one was holding the flow down for no reason and the cap
    % never bound. A staged cap throttles by a RATE, which has no such scale --
    % there the signal is a waiting room with no blocked mass in it at all.
    released = zeros(0,1);
    for k = 1:numel(active)
        c = active(k);
        if con.staged(c)
            rooms = find(stg.gatedBy(c,:));
            if isempty(rooms) || sum(sg(rooms)) < 1e-9
                released(end+1,1) = c; %#ok<AGROW>
            end
        elseif mult(k) > 1 + max(1e-9, tol)
            released(end+1,1) = c; %#ok<AGROW>
        end
    end
    if ~converged
        % A FAILED PASS SAYS NOTHING ABOUT WHICH CAPS BIND. Its release signal is
        % still information -- a multiplier that ran above one was on its way out
        % of the active set -- but its state is not, so no row is ADDED from it and
        % the next pass restarts from the last point that was a fixed point.
        violated = zeros(0,1);
        x = xOk; sigma2 = sigma2Ok;
        if isempty(released)
            break
        end
    end
    if isempty(violated) && isempty(released)
        break
    end
    active = setdiff(union(active, violated), released);
    active = active(:);
end
if ncon > 0 && aset == aset_max
    line_warning(mfilename, ['The active set did not settle: the same capacity constraints kept ' ...
        'binding and releasing. The reported point satisfies the last set tried.']);
end

% A CONVERGED POINT BEATS THE LAST ITERATE. The loop can end on a pass that did
% not converge -- a cap the closure cannot hold at any multiplier will be added,
% fail and be released for as long as the loop runs -- and the last iterate of
% such a pass is not a fixed point of anything. Report the last converged one
% instead, and say plainly which caps it does not meet rather than presenting a
% cap-violating point as the answer.
if ~isempty(best) && ~converged
    x = best.x; sg = best.sg; sigma2 = best.sigma2; mult = best.mult;
    activeSolved = best.active; resnorm = best.resnorm; clamped = best.clamped;
    converged = true;
    if ~best.feasible
        val = con.A*x;
        if stg.n > 0
            val = val + con.As*sg;
        end
        over = find(val > con.b + max(1e-9, tol));
        names = strjoin(con.label(over), ', ');
        line_warning(mfilename, sprintf(['No fixed point of the closure satisfies %s. The closure wants ' ...
            'more jobs there than the cap allows and no admission multiplier holds it: the reported point ' ...
            'is the converged fixed point of the model without that cap, and it exceeds the cap. Use ' ...
            'SolverCTMC, SolverJMT, SolverSSA or SolverLDES for this model.'], names));
    end
end
if ~converged
    % Report rather than return a point that is not a fixed point. The
    % substitution route would simply have stopped at outer_max with no
    % indication, which is the failure mode this replaces.
    line_warning(mfilename, sprintf(['The simultaneous closure solve stopped at residual %.3e after %d ' ...
        'Newton steps without reaching %.3e. The reported point is the last iterate.'], ...
        resnorm, iters, tol));
end

% the variance as it entered the DRIFT, which is what any later solve on this
% fixed point must close at; identical to sigma2 here because the stations
% held at zero are excluded from the unknowns rather than zeroed afterwards
sigma2drift = sigma2;
covdrift = covblk;

% RMET is the flow that actually CROSSES each event, which is what the throughput
% table must report, and RFIRE the rate the event fires at, which is what the
% diffusion counts: a lost arrival fires and lands nowhere, so the two differ by
% exactly the loss.
ctxFinal = struct('terms', terms, 'C', C, 'Nvec', Nvec, 'cidx', cidx, 'covblk', {covblk}, ...
    'nstate', nstate, 'conA', con.A, 'conAs', con.As, 'conb', con.b, ...
    'conStaged', con.staged, 'gates', gates, 'stg', stg, 'active', activeSolved, ...
    'clampT', []);
if clamped
    ctxFinal.clampT = local_clamp_tangent(con, activeSolved, terms);
end
[~, rmet, rfire] = local_residual([x; sg; sigma2(cidx); mult], ctxFinal, false);
[~, Sigma] = local_sigma(x, sigma2, terms, cidx, covblk, rfire, ctxFinal.clampT);

% -- transient: the same closure, integrated as an index-1 DAE ---------------
Sigmat = [];
QVart = [];
tvar = [];
switches = cell(0,3);
if isfinite(options.timespan(2))
    % UNDER A CAP THIS IS A HYBRID DAE, integrated segment by segment with the
    % binding set updated at each located crossing. See LOCAL_TRANSIENT: what the
    % steady state settles once with an active-set loop, the trajectory settles
    % again at every fill and every drain.
    [t, xvec_t, Sigmat, tvar, switches] = local_transient(options, terms, C, Nvec, cidx, ...
        covblk, sigma2drift, t, xvec_t, con, gates, stg);
    if ~isempty(Sigmat)
        QVart = local_qvar_t(Sigmat, terms);
    end
end

% -- performance measures ----------------------------------------------------
% Read back from the same event representation that defines the drift, so
% throughput balances flow at the fixed point. Identical to the reader in
% SOLVER_FLUID_MOMENTS: the closure changed how it is solved, not what it is.
% The THROTTLED rates, so the throughput reported at a blocked region is the
% flow that actually crosses it rather than the nominal one. Reading the
% unthrottled vector here would balance flow everywhere except across the cap,
% which is the one place the model was asked about.
r = rmet;
gfac = terms.factorFcn(x, sigma2drift, covdrift);

Seff = terms.S;
if ~isempty(terms.lldscaling)
    for i = 1:min(M,size(terms.lldscaling,1))
        Seff(i) = max(terms.S(i), max(terms.lldscaling(i,:)));
    end
end

QN = zeros(M,K); UN = zeros(M,K); TN = zeros(M,K); RN = zeros(M,K);
for i = 1:M
    for k = 1:K
        blk = terms.classBlock{i,k};
        if isempty(blk)
            continue
        end
        QN(i,k) = sum(x(blk));
        if terms.sched(i) == SchedStrategy.INF
            UN(i,k) = QN(i,k);
        else
            UN(i,k) = sum(gfac(blk)) / Seff(i);
        end
        TN(i,k) = sum(r(terms.evIsDeparture & terms.evStation == i & terms.evClass == k));
    end
end

% UTILIZATION MUST BE READ FROM THE FLOW THAT ACTUALLY CROSSES once a cap binds.
% Away from a constraint the in-service fluid sum(gfac)/s and the carried
% utilization T/(mu*s) are the same number, because the drift balances; under an
% ACTIVE cap they are not -- the multiplier throttles the departures (TN) and
% leaves the in-service fluid (gfac) alone, so the two columns disagreed: a closed
% tandem capped at 1 reported Util 0.688 at a station whose own Tput/mu was 0.550,
% and the exact answer is neither. Every other solver reports the CARRIED
% utilization (Util = X*D, the utilization law), and SolverCTMC gives 0.444 for
% the same station, so that is the convention the throttled point has to keep.
% Unconstrained runs are bit-identical: the branch is entered only when the active
% set is non-empty. see _kb/06-solver-catalog.md
if ~isempty(activeSolved)
    for i = 1:M
        if terms.sched(i) == SchedStrategy.INF || terms.sched(i) == SchedStrategy.EXT
            continue
        end
        for k = 1:K
            if isempty(terms.classBlock{i,k})
                continue
            end
            mu = sn.rates(i,k);
            if isfinite(mu) && mu > 0
                UN(i,k) = TN(i,k) / (mu * Seff(i));
            end
        end
    end
end

% See solver_fluid_closing.m: TN is zero only to the integrator's accuracy.
RN(TN > GlobalConstants.Zero) = QN(TN > GlobalConstants.Zero) ./ TN(TN > GlobalConstants.Zero);

% -- transient measures ------------------------------------------------------
% On the DAE route sigma2 varies along the trajectory, so the rate factors are
% read at the variance the trajectory HAD at each instant. SOLVER_FLUID_MOMENTS
% reads its whole transient at the single stationary variance, which is the
% quasi-steady-state approximation of this.
nt = size(xvec_t,1);
QNt = cell(M,K); UNt = cell(M,K); TNt = cell(M,K);
Gt = zeros(nt, nstate);
Rt = zeros(nt, numel(terms.rateBase));
for s = 1:nt
    xs = xvec_t(s,:)';
    if ~isempty(Sigmat) && s <= size(Sigmat,3)
        s2s = local_sigma2_from(Sigmat(:,:,s), terms, cidx);
    else
        s2s = sigma2drift;
    end
    Gt(s,:) = terms.factorFcn(xs, s2s, covdrift)';
    Rt(s,:) = terms.ratesFcn(xs, s2s, covdrift)';
end
for i = 1:M
    for k = 1:K
        blk = terms.classBlock{i,k};
        if isempty(blk)
            QNt{i,k} = zeros(nt,1); UNt{i,k} = zeros(nt,1); TNt{i,k} = zeros(nt,1);
            continue
        end
        QNt{i,k} = sum(xvec_t(:,blk),2);
        if terms.sched(i) == SchedStrategy.INF
            UNt{i,k} = QNt{i,k};
        else
            UNt{i,k} = sum(Gt(:,blk),2) / Seff(i);
        end
        TNt{i,k} = sum(Rt(:, terms.evIsDeparture & terms.evStation == i & terms.evClass == k), 2);
    end
end

% -- moment report -----------------------------------------------------------
QVar = zeros(M,K);
for i = 1:M
    for k = 1:K
        blk = terms.classBlock{i,k};
        if ~isempty(blk)
            QVar(i,k) = max(0, sum(sum(Sigma(blk,blk))));
        end
    end
end
moments = struct();
moments.Sigma = Sigma;
moments.QVar = QVar;
moments.QStd = sqrt(QVar);
moments.sigma2 = sigma2;
moments.sigma2Drift = sigma2drift;
moments.refinement = [];
moments.outerIters = iters;
moments.stationBlock = terms.stationBlock;
moments.classBlock = terms.classBlock;
% the DAE-only outputs: a second moment along the trajectory rather than at
% the fixed point alone
moments.Sigmat = Sigmat;
moments.QVart = QVart;
moments.tvar = tvar;
moments.residual = resnorm;
moments.converged = converged;
moments.conservation = local_conservation_error(C, Nvec, [x; sg]);
% what the capacity constraints did: which bound, how hard each throttled the
% admission flow, and how much slack the rest had
% What the capacity constraints did: which bound, how much mass each is holding
% outside its region, and how fast that queue drains. STAGING IS REPORTED HERE
% AND NOT FOLDED INTO QN, because a blocked job is at no station -- the same
% choice LDES makes, whose station queues likewise sum to less than N.
% What the capacity limits did: which bound, how much mass each holds outside its
% region, and the multiplier each one settled at -- a fraction of the admissions
% allowed for a cap that holds the job upstream or loses it, a drain rate for one
% that stages it. STAGING IS REPORTED HERE AND NOT FOLDED INTO QN, because a
% STAGED job is at no station -- the same choice LDES makes, whose station queues
% likewise sum to less than N. A HELD job is not staged and IS folded in, at the
% upstream station where the reference counts it. SWITCHES records every time the
% trajectory made a cap start or stop binding.
% Assembled field by field rather than with STRUCT(...): a cell value there makes
% a struct ARRAY, one element per label, and two cell values of different sizes
% are an error rather than two fields.
capacity = struct();
capacity.label = con.label;
capacity.b = con.b;
capacity.value = local_capacity_value(con, x, sg);
capacity.active = activeSolved;
capacity.staged = con.staged;
capacity.region = con.region;
capacity.station = con.station;
capacity.staging = sg;
capacity.stagingRegion = stg.region;
capacity.stagingClass = stg.class;
capacity.blocked = sum(sg);
capacity.drain = mult;
capacity.multiplier = mult;
capacity.switches = switches;
moments.capacity = capacity;

xvec_it = {x(:)'};
iters = iters + seed_iters;
runtime = toc(T0);
end

% ---------------------------------------------------------------------------
function [C, Nvec] = local_conservation(sn, terms, stg)
% Population conservation, one row per CLOSED chain.
%
% An open chain has no conserved population and contributes nothing. The EXT
% coordinates are excluded because the closing representation holds unit mass
% there as a normalisation constant, not as a job count (FLUID_MOMENT_TERMS).
M = terms.M; K = terms.K; nstate = terms.nstate;
isExt = terms.isExt(:);
coordClass = zeros(nstate,1);
coordStation = zeros(nstate,1);
for i = 1:M
    for c = 1:K
        if terms.Kic(i,c) > 0
            idx = terms.q_indices(i,c):(terms.q_indices(i,c)+terms.Kic(i,c)-1);
            coordClass(idx) = c;
            coordStation(idx) = i;
        end
    end
end
nstg = 0;
if nargin >= 3 && ~isempty(stg), nstg = stg.n; end
C = zeros(0,nstate+nstg);
Nvec = zeros(0,1);
if ~isfield(sn,'chains') || isempty(sn.chains)
    return
end
for ch = 1:size(sn.chains,1)
    inch = find(sn.chains(ch,:));
    if isempty(inch)
        continue
    end
    Nch = sum(sn.njobs(inch));
    if ~isfinite(Nch) || Nch <= 0
        continue % open chain: nothing is conserved
    end
    sel = ismember(coordClass, inch) & coordStation > 0 & ~isExt(max(coordStation,1));
    if ~any(sel)
        continue
    end
    row = zeros(1,nstate+nstg);
    row([sel(:); false(nstg,1)]) = 1;
    % A JOB IN THE WAITING QUEUE IS STILL IN THE CHAIN. Leaving the staging
    % coordinates out of this row would let the constraint balance while the
    % blocked mass quietly left the model.
    if nstg > 0
        row(nstate + find(ismember(stg.class, inch))) = 1;
    end
    C(end+1,:) = row; %#ok<AGROW>
    Nvec(end+1,1) = Nch; %#ok<AGROW>
end
end

% ---------------------------------------------------------------------------
function v = local_capacity_value(con, x, sg)
if isempty(con.b)
    v = zeros(0,1);
else
    v = con.A * x;
    if ~isempty(sg) && size(con.As,2) == numel(sg)
        v = v + con.As * sg;
    end
end
end

% ---------------------------------------------------------------------------
function err = local_conservation_error(C, Nvec, x)
if isempty(C)
    err = 0;
else
    err = max(abs(C*x - Nvec));
end
end

% ---------------------------------------------------------------------------
function [sigma2, Sigma] = local_sigma(x, sigma2in, terms, cidx, covblk, r, clampT)
% One Lyapunov solve. Sigma is linear in itself for a held x, so this is a
% SOLVE and not an iteration -- which is why Sigma stays out of the Newton
% unknowns and only the M numbers it projects onto go in.
%
% R IS THE FIRING RATE VECTOR when a cap is active, so the diffusion matrix
% D*diag(r)*D' counts the events that actually happen. The drift Jacobian is NOT
% re-derived for the multiplier: it is constant within a mode, so it rescales the
% admission columns, and treating it as constant here is a first-order
% approximation of A rather than the exact one. It is flagged rather than hidden
% because the Gaussian closure over a capacity-truncated support is the larger
% approximation of the two.
%
% CLAMPT IS THE TANGENT SPACE OF THE CAPS THAT CLAMP. A cap that holds the job
% upstream or loses it fixes the constrained combination A*x at b for as long as
% it binds -- blocking answers the state instantaneously -- so that combination
% does not fluctuate and the noise lives on the subspace orthogonal to A. WITHOUT
% IT THE LYAPUNOV SOLVE HAS NO SOLUTION AT ALL, and not by accident: with the
% multiplier held constant the drift along the constrained direction is neutral
% (an overloaded M/M/1/K sits on a whole line of equilibria), so the Jacobian
% carries a zero eigenvalue there and 'LINE:FluidNonHyperbolic' is the correct
% verdict for the UNPROJECTED system. Projecting the jump directions is enough to
% state the reduced problem, because FLUID_LYAPUNOV restricts everything to
% range(D) already.
%
% A STAGED cap is NOT clamped and must not be projected: there the drain rate is
% constant and the room population is a state, so the region population genuinely
% fluctuates and its restoring force is the room.
A = terms.jacFcn(x, sigma2in, covblk);
if nargin < 6 || isempty(r)
    r = terms.ratesFcn(x, sigma2in, covblk);
end
if nargin < 7
    clampT = [];
end
idx = terms.covIdx;
Dc = terms.D(idx,:);
if ~isempty(clampT)
    Dc = clampT * Dc;
end
Qc = Dc * diag(r) * Dc';
Sc = fluid_lyapunov(A(idx,idx), Qc, Dc);
Sigma = zeros(terms.nstate);
Sigma(idx,idx) = Sc;
sigma2 = local_sigma2_from(Sigma, terms, cidx);
end

% ---------------------------------------------------------------------------
function T = local_clamp_tangent(con, active, terms)
% Orthogonal projector onto the subspace the CLAMPING caps leave free.
%
% One row per active cap that holds or loses, restricted to the covariance
% coordinates, and the projector that annihilates them: I - R'*(R*R')^-1*R. Empty
% when no active cap clamps, which is every model without a station buffer and
% every region model -- so the uncapped and the staged paths form the same
% Lyapunov system they always did.
T = [];
rows = active(~con.staged(active));
if isempty(rows)
    return
end
idx = terms.covIdx;
R = con.A(rows, idx);
if isempty(R) || ~any(abs(R(:)) > GlobalConstants.Zero)
    return
end
T = eye(numel(idx)) - R' * pinv(R * R') * R;
end

% ---------------------------------------------------------------------------
function sigma2 = local_sigma2_from(Sigma, terms, cidx)
sigma2 = zeros(terms.M,1);
for j = 1:numel(cidx)
    i = cidx(j);
    blk = terms.stationBlock{i};
    if ~isempty(blk)
        sigma2(i) = max(0, sum(sum(Sigma(blk,blk))));
    end
end
end

% ---------------------------------------------------------------------------
% ---------------------------------------------------------------------------
function [G, rmet, rfire] = local_residual(u, ctx, quiet)
% The coupled algebraic system, stacked. Returns [] when the closure cannot be
% evaluated at this iterate so that the line search can back off; the caller
% evaluates the seed with quiet=false, where a genuine failure must surface.
%
% THE UNKNOWNS ARE [x; staging; sigma2; mult]. Staging carries the blocked mass
% that a capped region will not admit yet; MULT is one multiplier per ACTIVE cap,
% and which multiplier it is depends on where that cap sends the mass it stops --
% a FRACTION of the admissions allowed for a cap that holds or loses, the RATE its
% waiting room drains at for one that stages. See FLUID_CAPACITY_LEGS.
terms = ctx.terms;
nstate = ctx.nstate;
cidx = ctx.cidx;
covblk = ctx.covblk;
stg = ctx.stg;
nstg = stg.n;
active = ctx.active;
nact = numel(active);
ncl = numel(cidx);
gates = ctx.gates;
staged = ctx.conStaged;

x  = u(1:nstate);
sg = u(nstate+(1:nstg));
s2 = zeros(terms.M,1);
% The Newton unknowns are unconstrained but a variance, a population and a
% multiplier are not. LOCAL_NEWTON projects the iterate into the feasible box, so
% these clamps are guards rather than the mechanism; clamping alone would flatten
% the Jacobian at the boundary and pin the unknown there for good.
s2(cidx) = max(0, u(nstate+nstg+(1:ncl)));
mult = max(0, u(nstate+nstg+ncl+1:end));

if quiet
    try
        r = terms.ratesFcn(x, s2, covblk);
        % THE FIRING RATE, not the nominal one: a cap that holds the job upstream
        % suppresses the event itself, so the upstream station keeps the mass and
        % reports it -- which is what LINE's own semantics does with a closed job
        % that finds no room (State.arrivalIsLost returns an empty successor and
        % the departure does not fire).
        [rup, rin, ds] = fluid_capacity_legs(terms, gates, stg, staged, active, x, sg, r, mult, false);
        s2new = local_sigma(x, s2, terms, cidx, covblk, rup, ctx.clampT);
    catch
        G = []; rmet = []; rfire = [];
        return
    end
else
    r = terms.ratesFcn(x, s2, covblk);
    [rup, rin, ds] = fluid_capacity_legs(terms, gates, stg, staged, active, x, sg, r, mult, false);
    s2new = local_sigma(x, s2, terms, cidx, covblk, rup, ctx.clampT);
end

if isempty(gates.Dn)
    drift = terms.D * rup;
else
    % THE THREE LEGS. The removal leaves at the rate the event FIRES and the
    % arrival lands at the rate mass actually ARRIVES; Dn + Dp = D, so this
    % collapses to D*rup wherever the two agree.
    %
    % A LOST ARRIVAL IS RETURNED TO THE SOURCE POOL, which is what DNEXT is for.
    % The closing representation holds unit mass on the EXT pseudo-station as a
    % NORMALISATION rather than as a job count, and its drift row is a real
    % equation: departures back to the pool balance arrivals out of it. Scaling
    % only the arrival leg of a lost event would leave that row unbalanced by
    % exactly the loss and the algebraic system would have no solution -- so on
    % the EXT rows the lost mass goes back where it came from, which is precisely
    % "the job never entered". A job lost on its way out of a REAL station is
    % destroyed instead: that station's row keeps the firing rate, because the job
    % did leave it.
    drift = gates.Dn * rup + gates.DnExt * rin + gates.Dp * rin;
end

G = [drift(:); ds(:)];
if ~isempty(ctx.C)
    G = [G; ctx.C*[x; sg] - ctx.Nvec];
end
G = [G; s2(cidx) - s2new(cidx)];
if nact > 0
    % THE UNDIFFERENTIATED CONSTRAINT, on purpose. At a fixed point every flow
    % already balances, so the differentiated form A*D*r = 0 is satisfied by
    % anything and pins no multiplier at all. A*x = b does pin it, through the
    % dependence of the fixed point on the multiplier.
    val = ctx.conA(active,:)*x;
    if nstg > 0
        val = val + ctx.conAs(active,:)*sg;
    end
    G = [G; val - ctx.conb(active)];
end
rmet = rin;
rfire = rup;
end

% ---------------------------------------------------------------------------
function [mult, ok] = local_hold_multipliers(terms, gates, stg, con, active, x, sg, s2, m0)
% The multipliers that hold the active caps at this state, by small Newton.
%
% Each active cap contributes one equation, d/dt(A*x + As*s) = 0: the constrained
% quantity is AT the cap, so holding it there means the flow across the cap
% balances. Each contributes one unknown too -- a fraction for a cap that holds or
% loses, an admitted flow for one that stages -- so the system is square and small
% (one entry per binding cap, never more than a handful), and a finite-difference
% Newton on it costs a few drift evaluations.
%
% This is what makes an event RESTART consistent: ode15s and RODAS both need the
% algebraic unknowns to satisfy their equations at the initial point of an index-1
% DAE, and the integrator's own Newton cannot be asked for them before the first
% step. It is also the FEASIBILITY test at an activation -- a fraction above one
% means the cap would have to admit more than arrives, so it is not binding after
% all and must not be activated.
mult = zeros(0,1);
ok = true;
active = active(:);
if isempty(active)
    return
end
mult = m0(:);
F = local_hold_resid(mult, terms, gates, stg, con, active, x, sg, s2);
for it = 1:40
    if max(abs(F)) < max(1e-12, 1e-10*(max(abs(mult))+1))
        break
    end
    J = zeros(numel(F), numel(mult));
    for j = 1:numel(mult)
        h = max(1e-7*abs(mult(j)), 1e-9);
        mp = mult; mp(j) = mp(j) + h;
        J(:,j) = (local_hold_resid(mp, terms, gates, stg, con, active, x, sg, s2) - F)/h;
    end
    warnstate = warning('off','MATLAB:rankDeficientMatrix');
    step = -(J \ F);
    warning(warnstate);
    if any(~isfinite(step))
        step = -(pinv(J)*F);
    end
    lam = 1; stepped = false;
    for ls = 1:20
        mn = max(0, mult + lam*step);
        Fn = local_hold_resid(mn, terms, gates, stg, con, active, x, sg, s2);
        if max(abs(Fn)) < max(abs(F))
            mult = mn; F = Fn; stepped = true;
            break
        end
        lam = lam/2;
    end
    if ~stepped
        break
    end
end
ok = max(abs(F)) < 1e-6;
end

% ---------------------------------------------------------------------------
function F = local_hold_resid(mult, terms, gates, stg, con, active, x, sg, s2)
r = terms.ratesFcn(x, s2, cell(terms.M,1));
[rup, rin, ds] = fluid_capacity_legs(terms, gates, stg, con.staged, active, x, sg, r, mult, true);
if isempty(gates.Dn)
    dx = terms.D * rup;
else
    dx = gates.Dn * rup + gates.DnExt * rin + gates.Dp * rin;
end
F = con.A(active,:)*dx;
if stg.n > 0
    F = F + con.As(active,:)*ds;
end
end

% ---------------------------------------------------------------------------
function [t, xvec_t, Sigmat, tvar, switches] = local_transient(options, terms, C, Nvec, cidx, ...
    covblk, sigma2ss, tfall, xfall, con, gates, stg)
% The transient closure as an index-1 DAE, integrated with a singular mass
% matrix.
%
% THE INTEGRATOR IS DISPATCHED, BUT THROUGH ITS OWN SLOT. ODE_SOLVE routes to
% options.odesolvers.accurateOdeSolver, which defaults to ode113 -- an explicit
% Adams method that cannot carry a singular mass matrix at all -- and
% ODE_SOLVE_STIFF may route to ode23s, which cannot carry a SINGULAR one. A DAE
% needs neither of those preferences, so it reads
% OPTIONS.ODESOLVERS.DAESOLVER, the fifth slot of the same struct, which
% defaults to @ode15s.
%
% THE OTHER CHOICE IS @RODAS, and it is why that slot exists. ode15s is MATLAB's
% and has no counterpart in the C++, native Python or JAR ports of this method:
% LSODA integrates y' = f and scipy's BDF/Radau take no mass matrix, so all three
% run the vendored Hairer-Wanner RODAS instead. Setting
% options.odesolvers.daeSolver = @rodas runs the same step sequence here, which is
% what makes the four codebases agree in the last digits; the default stays ode15s
% so an existing MATLAB answer does not move underneath a caller who did not ask.
% Both honour ODESET's EVENTS, which is what the hybrid route below needs.
%
% A SEPARATE SLOT RATHER THAN ACCURATESTIFFODESOLVER, which is also @ode15s
% today: the rest of the fluid solver sets NonNegative on its odeset, RODAS has
% no hook for holding a component at or above zero, and rodas.m refuses that
% field rather than ignoring it. Overloading the shared slot would therefore
% break every other fluid path the moment RODAS was selected for this one.
%
% UNDER A CAP THIS IS A HYBRID DAE, and it is integrated as one rather than
% refused. The system SWITCHES every time a cap starts or stops binding, so the
% horizon is covered by SEGMENTS: each segment is one index-1 DAE with a fixed set
% of binding caps, the segment ends at a located crossing, and the next one starts
% from that state with the set updated. Integrating with the set frozen instead
% would silently report the unconstrained trajectory through a cap the model
% declares, which is why this was refused before it was implemented.
%
% THE MULTIPLIER IS AN ALGEBRAIC UNKNOWN HERE, one per cap, carried in the state
% vector with a zero mass row, and its equation is the DIFFERENTIATED constraint
% d/dt(A*x + As*s) = 0. The undifferentiated form A*x = b would be an index-2 DAE
% -- the multiplier appears only after one differentiation -- and neither ode15s
% nor RODAS solves index 2. A cap that is not binding keeps its unknown pinned at
% the inert value (one for a fraction, zero for a flow) by an equation of its own,
% so the state layout is the same on every segment and a restart costs no
% re-indexing.
%
% A STAGED cap's unknown is the admitted FLOW and not the drain rate the steady
% state solves for; see FLUID_CAPACITY_LEGS for why a rate cannot start a constrained
% phase.
nstate = terms.nstate;
idx = terms.covIdx;
nc = numel(idx);
ncon = numel(con.b);
nstg = stg.n;
switches = cell(0,3);

maxcov = 25;
if isfield(options,'config') && isfield(options.config,'dae_maxcov') && ~isempty(options.config.dae_maxcov)
    maxcov = options.config.dae_maxcov;
end
% The covariance adds nc^2 differential states and ode15s forms its Jacobian by
% finite differences over all of them, so this grows as nc^4. Above the cap the
% mean is still integrated as a DAE -- conservation stays an equation -- but the
% variance is held at its stationary value, which is what 'minnormal' does for
% the whole transient anyway.
withcov = nc > 0 && nc <= maxcov;
if nc > maxcov
    line_warning(mfilename, sprintf(['The transient covariance would add %d differential states, above ' ...
        'the limit of %d set by options.config.dae_maxcov. Integrating the mean as a DAE with the ' ...
        'variance held at its stationary value; raise that limit for a time-varying second moment.'], ...
        nc*nc, maxcov*maxcov));
end

% The initial condition is taken from the seed trajectory rather than from
% options.init_sol so that it is certainly in the closing state layout that
% TERMS indexes; the two agree whenever every class carries a rate, and this
% does not depend on that.
x0 = xfall(1,:)';
if numel(x0) ~= nstate
    t = tfall; xvec_t = xfall; Sigmat = []; tvar = [];
    return
end

% One differential equation per closed chain is redundant -- the rows of D sum
% to zero there -- so one is replaced by the constraint rather than added to
% it. The row dropped is the one carrying the most mass at t=0, which keeps the
% algebraic equation away from a coordinate that is identically zero. Chains
% partition the classes, so no index is claimed twice.
algrow = zeros(size(C,1),1);
for r = 1:size(C,1)
    sup = find(C(r,1:nstate));
    [~, p] = max(x0(sup));
    algrow(r) = sup(p);
end

ncov = 0;
if withcov
    ncov = nc*nc;
end
nz = nstate + nstg + ncon + ncov;
oSg = nstate; oM = nstate + nstg; oCov = nstate + nstg + ncon;
inert = ones(ncon,1);
if ncon > 0
    inert(con.staged) = 0;   % a flow is inert at zero, a fraction at one
end

% A STATE ABOVE A CAP IS NOT A STATE THE MODEL CAN BE IN. The initial point comes
% from the caller -- the default for a closed model spreads the population over
% the stations -- so it may sit outside a cap; holding the cap from there would
% freeze the violation for the whole horizon, since the equation says the
% constrained quantity does not MOVE. Start on the cap instead, and say so.
%
% THE EXCESS IS MOVED, NOT DROPPED. Conservation is an algebraic row here, so an
% initial state that does not satisfy it is an INCONSISTENT initialisation and no
% index-1 solver may be handed one. The mass goes where the model would have put
% it: into the waiting room of a staged cap, and otherwise onto the coordinates
% that feed the capped stations, which is where a held job waits.
tol = options.tol;
if isempty(tol) || ~isfinite(tol)
    tol = GlobalConstants.FineTol;
end
sg0 = zeros(nstg,1);
excess = zeros(max(ncon,1),1);
over = zeros(0,1);
if ncon > 0
    over = find((con.A*x0) > con.b + max(1e-9, tol));
    for k = 1:numel(over)
        c = over(k);
        val = con.A(c,:)*x0;
        if ~(val > con.b(c) && con.b(c) > 0)
            continue
        end
        sel = con.A(c,:) > 0;
        excess(c) = sum(x0(sel)) * (1 - con.b(c)/val);
        x0(sel) = x0(sel) * (con.b(c)/val);
    end
end

% THE INITIAL MODE. A cap the initial state sits ON is already binding, so the
% first segment must carry it: starting inactive would integrate one step of the
% unconstrained system straight through the bound. Feasibility decides -- a
% multiplier above one means the drift is pulling the constrained quantity DOWN
% and the cap is not binding after all.
mInit = inert;
active = zeros(0,1);
if ~isempty(over)
    trial = sort(over(:));
    m0 = double(~con.staged(trial));
    [mm, ok] = local_hold_multipliers(terms, gates, stg, con, trial, x0, sg0, sigma2ss, m0);
    keep = trial(ok & (con.staged(trial) | mm <= 1 + 1e-9));
    if ~isempty(keep) && numel(keep) < numel(trial)
        [mm, ok] = local_hold_multipliers(terms, gates, stg, con, keep, x0, sg0, sigma2ss, ...
            double(~con.staged(keep)));
        trial = keep;
    end
    if ok
        active = trial(:);
        mInit(active) = mm;
    end
end
% and the mass the clip removed goes where the model would hold it
for k = 1:numel(over)
    c = over(k);
    if excess(c) <= 0
        continue
    end
    rooms = zeros(1,0);
    if nstg > 0 && con.staged(c) && any(active == c)
        rooms = find(stg.gatedBy(c,:));
    end
    if ~isempty(rooms)
        sg0(rooms) = sg0(rooms) + excess(c)/numel(rooms);
        continue
    end
    sel = con.A(c,:) > 0;
    pool = ~sel;
    if ~isempty(gates.Dn)
        feeders = false(nstate,1);
        for e = find(gates.gate(c,:))
            feeders = feeders | (gates.Dn(:,e) < 0) | (gates.DnExt(:,e) < 0);
        end
        if any(feeders(:) & pool(:))
            pool = feeders(:)' & pool;
        end
    end
    w = x0(pool);
    if sum(w) > GlobalConstants.FineTol
        x0(pool) = x0(pool) + excess(c) * w / sum(w);
    elseif any(pool)
        x0(pool) = x0(pool) + excess(c)/sum(pool);
    end
end
if ~isempty(over)
    line_warning(mfilename, sprintf(['The initial state holds more jobs than %s allows; the transient ' ...
        'starts on the cap, with the excess where the model would hold it -- in the waiting room, or at ' ...
        'the stations feeding the cap.'], strjoin(con.label(over), ', ')));
end

% Sigma(0) = 0 is the consistent initialisation, and it is also the physically
% right one: the population at t=0 is a known deterministic state, so it has no
% variance. C*x0 = N holds by construction, so the algebraic rows are satisfied
% at t=0 and no separate consistency solve is needed.
z0 = zeros(nz,1);
z0(1:nstate) = x0;
if nstg > 0
    z0(oSg+(1:nstg)) = sg0;
end
if ncon > 0
    z0(oM+(1:ncon)) = mInit;
end

daesolver = @ode15s;
if isfield(options,'odesolvers') && isfield(options.odesolvers,'daeSolver') ...
        && ~isempty(options.odesolvers.daeSolver)
    daesolver = options.odesolvers.daeSolver;
end

t = zeros(0,1);
Z = zeros(0,nz);
tcur = options.timespan(1);
z = z0;
armed = false(max(ncon,1),1);
if ncon > 0 && nstg > 0
    for k = 1:numel(active)
        c = active(k);
        if con.staged(c) && sum(sg0(stg.gatedBy(c,:))) > 1e-8
            armed(c) = true;
        end
    end
end
segMax = 4*ncon + 8;
frozen = false;
seg = 0;
while seg < segMax
    seg = seg + 1;
    gatedRooms = false(nstg,1);
    for k = 1:numel(active)
        if con.staged(active(k))
            gatedRooms = gatedRooms | stg.gatedBy(active(k),:)';
        end
    end
    diagM = ones(nz,1);
    diagM(algrow) = 0;
    if nstg > 0
        diagM(oSg+(1:nstg)) = double(gatedRooms);
    end
    if ncon > 0
        diagM(oM+(1:ncon)) = 0;
    end
    Mm = spdiags(diagM, 0, nz, nz);

    odefun = @(tt, zz) local_dae_rhs(zz, terms, C, Nvec, cidx, covblk, idx, algrow, ...
        withcov, sigma2ss, nstate, nc, con, gates, stg, active, gatedRooms, ...
        inert, oSg, oM, oCov);
    odeopt = odeset('Mass', Mm, 'MassSingular', 'yes', ...
        'RelTol', max(tol, 1e-10), 'AbsTol', max(tol, 1e-12));
    if ncon > 0
        odeopt = odeset(odeopt, 'Events', @evfun);
    end

    try
        if ncon > 0
            [ts, Zs, te, ze, ie] = feval(daesolver, odefun, ...
                [tcur, options.timespan(2)], z, odeopt);
        else
            [ts, Zs] = feval(daesolver, odefun, [tcur, options.timespan(2)], z, odeopt);
            te = []; ze = []; ie = [];
        end
    catch me
        line_warning(mfilename, sprintf(['The transient DAE failed (%s); reporting the seed trajectory ' ...
            'and the stationary variance instead.'], me.message));
        t = tfall; xvec_t = xfall; Sigmat = []; tvar = [];
        return
    end

    % THE OVERSHOOTING STEP IS NOT REPORTED, which is what makes the located
    % crossing worth locating: the reported path stops AT the cap and the next
    % segment starts there.
    if ~isempty(te)
        keepPts = ts < te(1) - eps(te(1));
        t = [t; ts(keepPts); te(1)]; %#ok<AGROW>
        Z = [Z; Zs(keepPts,:); ze(1,:)]; %#ok<AGROW>
    else
        t = [t; ts]; %#ok<AGROW>
        Z = [Z; Zs]; %#ok<AGROW>
    end

    if isempty(te)
        break % the horizon was reached with this set of caps binding
    end

    tcur = te(1);
    z = ze(1,:)';
    c = ie(1);
    xh = z(1:nstate);
    sgh = zeros(nstg,1);
    if nstg > 0
        sgh = z(oSg+(1:nstg));
    end
    s2h = sigma2ss;
    if withcov
        Sc = reshape(z(oCov+1:end), nc, nc);
        Sc = (Sc + Sc')/2;
        Sg = zeros(nstate); Sg(idx,idx) = Sc;
        s2h = local_sigma2_from(Sg, terms, cidx);
    end
    if any(active == c)
        active = active(active ~= c);
        armed(c) = false;
        z(oM+c) = inert(c);
        switches(end+1,:) = {tcur, c, 'release'}; %#ok<AGROW>
    else
        % ACTIVATING NEEDS THE MULTIPLIER THAT HOLDS THE CAP, both because the
        % restart must satisfy the algebraic rows and because that multiplier IS
        % the feasibility test: above one, the cap would have to admit more than
        % arrives, so it is not binding after all.
        trial = sort([active; c]);
        m0 = zeros(numel(trial),1);
        for k = 1:numel(trial)
            if any(active == trial(k))
                m0(k) = z(oM+trial(k));
            elseif con.staged(trial(k))
                m0(k) = 0;
            else
                m0(k) = 1;
            end
        end
        [mm, ok] = local_hold_multipliers(terms, gates, stg, con, trial, xh, sgh, s2h, m0);
        feasible = ok && all(con.staged(trial) | mm <= 1 + 1e-9);
        if feasible
            active = trial(:);
            z(oM+trial) = mm;
            switches(end+1,:) = {tcur, c, 'activate'}; %#ok<AGROW>
        else
            % a crossing whose cap cannot be held: stop rather than freeze a
            % violation into the algebraic rows
            switches(end+1,:) = {tcur, c, 'unheld'}; %#ok<AGROW>
            frozen = true;
            break
        end
    end
    if tcur >= options.timespan(2) - 1e-12
        break
    end
end
if seg >= segMax
    frozen = true;
end
if frozen
    line_warning(mfilename, sprintf(['The transient stopped switching after %d segments at t=%g of %g: ' ...
        'the reported trajectory ends there rather than continuing with a capacity constraint that does ' ...
        'not hold. Shorten options.timespan, or use SolverCTMC/SolverJMT/SolverSSA/SolverLDES.'], ...
        seg, tcur, options.timespan(2)));
end

xvec_t = Z(:,1:nstate);
if withcov
    Sigmat = zeros(nstate, nstate, numel(t));
    for s = 1:numel(t)
        Sc = reshape(Z(s, oCov+1:end), nc, nc);
        Sc = (Sc + Sc')/2; % the equation preserves symmetry; rounding does not
        Sg = zeros(nstate);
        Sg(idx,idx) = Sc;
        Sigmat(:,:,s) = Sg;
    end
    tvar = t;
else
    Sigmat = [];
    tvar = [];
end

    function [value, isterminal, direction] = evfun(~, zz)
        % The event functions, all of them, as one vector: a cap that is not
        % active is watched for REACHING its bound, one that is active for
        % stopping to bind -- a fraction above one, or a room that has emptied.
        %
        % ARMED IS A LATCH, and it has to be: a staged cap activates with an EMPTY
        % room, so "the room emptied" is true at the instant it starts binding and
        % the release would fire immediately. The latch only ever goes false ->
        % true, so a speculative evaluation of this function cannot corrupt it.
        value = ones(ncon,1);
        isterminal = ones(ncon,1);
        direction = -ones(ncon,1);
        xx = zz(1:nstate);
        ss = zeros(nstg,1);
        if nstg > 0
            ss = zz(oSg+(1:nstg));
        end
        val = con.A*xx;
        if nstg > 0
            val = val + con.As*ss;
        end
        for cc = 1:ncon
            if any(active == cc)
                if con.staged(cc)
                    rooms = find(stg.gatedBy(cc,:));
                    mass = 0;
                    if ~isempty(rooms)
                        mass = sum(ss(rooms));
                    end
                    if ~armed(cc) && mass > 1e-8
                        armed(cc) = true;
                    end
                    if armed(cc)
                        value(cc) = mass;
                    end
                else
                    value(cc) = 1 - zz(oM+cc);
                end
            else
                value(cc) = con.b(cc) - val(cc);
            end
        end
    end
end

% ---------------------------------------------------------------------------
function dz = local_dae_rhs(z, terms, C, Nvec, cidx, covblk, idx, algrow, withcov, ...
    sigma2ss, nstate, nc, con, gates, stg, active, gatedRooms, inert, oSg, oM, oCov)
% The right-hand side of one segment of the hybrid DAE.
%
% The mass matrix has zeroed the algebraic rows, so what is written there is the
% RESIDUAL and not a derivative: population conservation on one coordinate per
% closed chain, emptiness on a waiting room with no cap above it, and the
% DIFFERENTIATED constraint on every multiplier -- d/dt(A*x + As*s) = 0 where the
% cap binds, and the inert value where it does not.
ncon = numel(con.b);
nstg = stg.n;
x = z(1:nstate);
sg = zeros(nstg,1);
if nstg > 0
    sg = z(oSg+(1:nstg));
end
m = zeros(ncon,1);
if ncon > 0
    m = z(oM+(1:ncon));
end
if withcov
    Sc = reshape(z(oCov+1:end), nc, nc);
    Sc = (Sc + Sc')/2;
    Sigma = zeros(nstate);
    Sigma(idx,idx) = Sc;
    s2 = local_sigma2_from(Sigma, terms, cidx);
else
    Sc = [];
    s2 = sigma2ss;
end

r = terms.ratesFcn(x, s2, covblk);
if ncon > 0
    % THE TRANSIENT MULTIPLIER OF A STAGED CAP IS A FLOW, not the drain rate the
    % steady state solves for: at the instant a region fills its room is empty, so
    % theta*s is zero however large theta is. See FLUID_CAPACITY_LEGS.
    [rup, rin, ds] = fluid_capacity_legs(terms, gates, stg, con.staged, active, x, sg, r, m(active), true);
    dx = gates.Dn * rup + gates.DnExt * rin + gates.Dp * rin;
else
    rup = r; rin = r; ds = zeros(0,1);
    dx = terms.D * r;
end

dzx = dx;
for rr = 1:numel(algrow)
    val = C(rr,1:nstate)*x;
    if nstg > 0
        val = val + C(rr,nstate+1:end)*sg;
    end
    dzx(algrow(rr)) = val - Nvec(rr);
end

dz = zeros(numel(z),1);
dz(1:nstate) = dzx;
if nstg > 0
    dsg = ds;
    dsg(~gatedRooms) = sg(~gatedRooms);
    dz(oSg+(1:nstg)) = dsg;
end
if ncon > 0
    dm = zeros(ncon,1);
    for c = 1:ncon
        if any(active == c)
            val = con.A(c,:)*dx;
            if nstg > 0
                val = val + con.As(c,:)*ds;
            end
            dm(c) = val;
        else
            dm(c) = m(c) - inert(c);
        end
    end
    dz(oM+(1:ncon)) = dm;
end
if withcov
    A = terms.jacFcn(x, s2, covblk);
    Dc = terms.D(idx,:);
    Ac = A(idx,idx);
    dS = Ac*Sc + Sc*Ac' + Dc*diag(rup)*Dc';
    dz(oCov+1:end) = dS(:);
end
end

% ---------------------------------------------------------------------------
function QVart = local_qvar_t(Sigmat, terms)
M = terms.M; K = terms.K;
nt = size(Sigmat,3);
QVart = cell(M,K);
for i = 1:M
    for k = 1:K
        blk = terms.classBlock{i,k};
        if isempty(blk)
            QVart{i,k} = zeros(nt,1);
            continue
        end
        v = zeros(nt,1);
        for s = 1:nt
            v(s) = max(0, sum(sum(Sigmat(blk,blk,s))));
        end
        QVart{i,k} = v;
    end
end
end
