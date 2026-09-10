function [QN, UN, RN, TN, xvec_it, QNt, UNt, TNt, xvec_t, t, iters, runtime, moments] = solver_fluid_moments(sn, options)
% [QN, UN, RN, TN, XVEC_IT, QNT, UNT, TNT, XVEC_T, T, ITERS, RUNTIME, MOMENTS] = SOLVER_FLUID_MOMENTS(SN, OPTIONS)
%
% Second-order moment-closure fluid analysis, backing options.method
% 'minnormal' and 'refined'.
%
% The default fluid methods close the moment hierarchy at first order: the
% drift of the mean depends on E[min(X_i,c_i)], which they replace by
% min(E[X_i],c_i). No second moment ever enters, so no variance is produced
% and the mean itself is biased wherever min() is not locally linear. The
% two methods here reinstate the second moment:
%
%   'minnormal' Min-normal closure (Guenther, Stefanek, Bradley). The drift
%              uses E[min(X_i,c_i)] under a normal marginal whose variance is
%              produced by the covariance equation, so mean and covariance are
%              solved self-consistently by fixed-point iteration. This
%              corrects the mean, most visibly near rho = 1 where the
%              first-order closure is worst.
%
%   'refined'  Refined mean field (Gast). Adds the O(1/N) correction term of
%              the mean-field expansion to the MEAN-FIELD fixed point (not to
%              the 'minnormal' one, which already resums it), computed by
%              FLUID_REFINE_MEANFIELD.
%
% All performance measures are read back from the same event representation
% that defines the drift, so throughputs balance flow at the fixed point
% under whichever closure was used.
%
% Parameters:
%   sn      - NetworkStruct
%   options - solver options; options.method selects the closure
%
% Returns:
%   as SOLVER_FLUID_CLOSING, plus MOMENTS with fields Sigma (state-level
%   covariance), QVar (station-class queue-length variance), QStd, sigma2
%   (per-station population variance), refinement (1/N correction, 'refined'
%   only), outerIters
%
% See also FLUID_MOMENT_TERMS, FLUID_LYAPUNOV, FLUID_REFINE_MEANFIELD.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
M = sn.nstations;
K = sn.nclasses;
method = options.method;

terms = fluid_moment_terms(sn, options);

% the covariance is a dense nstate-by-nstate object and the Lyapunov solve is
% cubic in it, so refuse rather than silently crawl
maxstate = 200;
if isfield(options,'config') && isfield(options.config,'moment_maxstate') && ~isempty(options.config.moment_maxstate)
    maxstate = options.config.moment_maxstate;
end
if terms.nstate > maxstate
    line_error(mfilename,sprintf(['The moment-closure methods solve a %dx%d Lyapunov equation, above the ' ...
        'limit of %d set by options.config.moment_maxstate. Raise that limit or use ' ...
        'options.method=''closing''.'], terms.nstate, terms.nstate, maxstate));
end

sigma2 = zeros(M,1);
% The last variance whose Lyapunov solve SUCCEEDED, and the step taken toward
% the next one. See the damping in the outer loop for why an iterate that fails
% must retreat toward this rather than end the method.
sigma2ok = zeros(M,1);
covok = cell(M,1);
damp_min = 1/64;
outer_max = 20;
if isfield(options,'iter_max') && ~isempty(options.iter_max)
    outer_max = min(outer_max, max(2, options.iter_max));
end
% THE CLOSURE IS JUDGED AT 1e-8 RELATIVE, so it must not stop at 1e-3. This
% loop used to converge to CoarseTol, and that is not a fixed point to two
% hosts: measured on mqn_singleserver_ps, QLen at Queue1/OpenClass came back
%
%   iter_tol   picard03        picard06        cross-host relative
%   1e-3       42.396178997    42.818762547    9.9e-03  (picard03 declined)
%   1e-4       42.820613986    42.820611415    6.0e-08
%   1e-6       42.820599087    42.820598896    4.5e-09
%   1e-8       42.820599087    42.820598828    6.1e-09
%
% so 1e-6 is the loosest that lands inside the 1e-8 gate RUNTESTEXAMPLE holds
% SolverFLD to, and it costs ~20% over 1e-4 on that model. MIN, not plain
% assignment, so a caller asking for TIGHTER still gets it, exactly as
% OUTER_MAX above is bounded by ITER_MAX.
%
% It governs the INNER mean solve too, through MEANOPT below. That is not
% incidental: the transient non-hyperbolic iterate this method used to abort on
% was an artefact of a loosely converged inner solve, so tightening the outer
% loop alone would leave the abort in place.
mom_tol = 1e-6;
if isfield(options,'iter_tol') && ~isempty(options.iter_tol)
    mom_tol = min(mom_tol, options.iter_tol);
end
outer_tol = mom_tol;

meanopt = options;
meanopt.method = 'closing'; % the closure enters through config.moment_sigma2
meanopt.iter_tol = mom_tol;

Sigma = [];
iters = 0;
sigma2solve = sigma2; % the variance the mean solve actually used
% the DPS capacity share is a RATIO of populations, so closing it needs the
% covariance BETWEEN the station coordinates, not only the station total; the
% blocks are carried alongside sigma2 through the same fixed point
covblk = cell(M,1);
covsolve = covblk;
shareSched = terms.sched == SchedStrategy.PS | terms.sched == SchedStrategy.FCFS | ...
    terms.sched == SchedStrategy.DPS | terms.sched == SchedStrategy.GPS;
% A DELAY STATION HAS NO min() TO CLOSE, so its variance must never reach the
% drift -- only the report. This mask used to be applied to SIGMA2DRIFT after
% the loop and nowhere inside it, so every mean solve of the fixed point ran
% with the delay variance switched on. The rate factor there is mu*n, which the
% Gaussian correction turns into something that does not vanish with n: the
% coordinate is driven NEGATIVE, the drift is conservative so another
% coordinate grows to match, and the trajectory leaves the simplex for good. On
% CQN_Cox_CS_9 (Delay + PS + PS(c=5), N=6) the first window past sigma2 = 0
% moved 8.7e3 of mass and the drift norm reached 5.4e9; with NonNegative on,
% the integrator hid the blow-up by clamping and simply never returned.
noDriftVar = false(M,1);
for i = 1:M
    noDriftVar(i) = terms.sched(i) == SchedStrategy.INF || terms.sched(i) == SchedStrategy.EXT;
end
for outer = 1:outer_max
    % A TRANSIENT ITERATE MUST NOT VETO THE METHOD. The gate below asks whether
    % the LINEAR NOISE APPROXIMATION has a stationary covariance at the point
    % this iterate landed on; a fixed point that fails it is a model this
    % closure cannot answer. An INTERMEDIATE iterate that fails it is not --
    % it is a variance step that overshot, and the alternation has not finished.
    % Letting one abort the solve threw away answers the method can reach:
    % mqn_singleserver_ps ran iterate 1 at maxRealEig = -4.93e-03 (stable),
    % iterate 2 at +1.07e+01 (declined here), and the fixed point the fallback
    % then converged to was stable at -4.92e-03 -- one excursion between two
    % stable states. Worse, which side of it an iterate lands on is decided by
    % the host's last bits, so the SAME model answered 42.3962 on picard03 and
    % picard06 and 42.8207 on picard09, a 1e-2 relative spread against a 1e-8
    % gate. So an iterate that fails RETREATS toward the last variance that
    % succeeded, halving the step until the LNA is defined again; only a step
    % below DAMP_MIN, i.e. a failure arbitrarily close to an accepted point,
    % is the model's own non-hyperbolicity and still throws.
    step = 1;
    while true
        sigma2try = sigma2ok + step*(sigma2 - sigma2ok);
        covtry = local_blend_cov(covok, covblk, step);
        sigma2solve = sigma2try;
        covsolve = covtry;
        sigma2solve(noDriftVar) = 0;
        covsolve(noDriftVar) = {[]};
        % A station whose occupancy cannot reach its server count has min(n,c) = n
        % on the whole support, so the closure there must stay first order: see
        % FLUID_MOMENT_TERMS, which decides it from the chain populations. Its
        % share closure follows, because mu_r*(n_r/n)*min(n,c) collapses to mu_r*n_r
        % once the min is the identity. The covariance is still SOLVED for these
        % stations and still reported, it just does not enter the drift, exactly as
        % at the delay stations below.
        sigma2solve(terms.minExact) = 0;
        covsolve(terms.minExact) = {[]};
        meanopt.config.moment_sigma2 = sigma2solve;
        meanopt.config.moment_cov = covsolve;
        [~, xvec_it, ~, ~, xvec_t, t, inner_iters] = solver_fluid(sn, meanopt);
        iters = iters + inner_iters;

        x = xvec_it{end}(:);
        r = terms.ratesFcn(x, sigma2solve, covsolve);
        % A POINT ON A SATURATION KINK HAS NO JACOBIAN. The rate factor min(n_i,c_i)
        % has slope 1 below c_i and 0 above, and FLUID_DRIFT_JACOBIAN resolves the tie
        % onto the saturated side, so the verdict read off it would depend on which
        % side the integrator stopped. The VERDICT, not the point, has to be
        % side-independent: both one-sided Jacobians are ordinary matrices, so ASK
        % BOTH and decline only when a side fails. Refusing at every kink instead
        % throws away models the reference solves -- the first outer iterate runs at
        % sigma2 = 0, and a saturated model's first-order fixed point lands on the
        % kink by construction. Later iterates carry a positive sigma2 and are
        % smooth, so this costs two Jacobians on the seed and nothing after it. Twin
        % of FluidRateFactors.driftKinkStation/nudgedOffKink in the JAR.
        kink = local_kink_stations(terms, x, sigma2solve);
        if ~isempty(kink)
            for rel = [-1e-6, 1e-6]
                xs = local_nudge_off_kink(terms, x, kink, rel);
                try
                    local_lyapunov(terms.jacFcn(xs, sigma2solve, covsolve), ...
                        terms.ratesFcn(xs, sigma2solve, covsolve), terms);
                catch me
                    if ~strcmp(me.identifier, 'LINE:FluidNonHyperbolic')
                        rethrow(me);
                    end
                    throw(MException('LINE:FluidNonHyperbolic', ...
                        ['[%s.m] The fluid fixed point sits on the saturation kink of station %d ' ...
                        '(population equals its %g servers) and the two one-sided drift Jacobians ' ...
                        'there disagree on hyperbolicity, so which of them the linear noise ' ...
                        'approximation would use is decided by the integrator''s rounding residue ' ...
                        'rather than by the model. This is the saturated boundary of a continuum ' ...
                        'of equilibria; use options.method=''closing'' for the mean only. ' ...
                        'Underlying: %s'], mfilename, kink(1), terms.S(kink(1)), me.message));
                end
            end
        end
        A = terms.jacFcn(x, sigma2solve, covsolve);
        try
            Sigma = local_lyapunov(A, r, terms);
            break
        catch me
            % Nothing to retreat toward on the SEED iterate: sigma2 is still the
            % zero the alternation starts at, so a halved step re-solves the same
            % point. That failure IS the model's, and must cost one solve, not
            % seven -- every model that legitimately declines does so right here.
            atSeed = norm(sigma2 - sigma2ok, 1) == 0 && all(cellfun(@isempty, covblk));
            if ~strcmp(me.identifier, 'LINE:FluidNonHyperbolic') || step <= damp_min || atSeed
                rethrow(me);
            end
            step = step/2;
            line_debug(options, ['FLD minnormal: iterate %d is not hyperbolic (%s); ' ...
                'halving the variance step to %g and re-solving'], outer, me.message, step);
        end
    end
    sigma2ok = sigma2try;
    covok = covtry;

    sigma2new = zeros(M,1);
    covnew = cell(M,1);
    for i = 1:M
        blk = terms.stationBlock{i};
        if ~isempty(blk)
            sigma2new(i) = max(0, sum(sum(Sigma(blk,blk))));
            if shareSched(i)
                covnew{i} = Sigma(blk,blk);
            end
        end
    end

    delta = norm(sigma2new - sigma2try, 1) / max(1, norm(sigma2new, 1));
    for i = 1:M
        % sigma2 is the SUM of a block, so it can converge while the
        % off-diagonals the share closure reads are still moving
        if ~isempty(covnew{i})
            dc = covnew{i};
            if ~isempty(covtry{i})
                dc = dc - covtry{i};
            end
            delta = max(delta, norm(dc,1) / max(1, norm(covnew{i},1)));
        end
    end
    sigma2 = sigma2new;
    covblk = covnew;
    if delta < outer_tol
        break
    end
end
outerIters = outer;

% Metrics must be read at the SAME variance the mean solve used, not at the
% variance that solve produced. Using the latter evaluates the rate functions
% at a point that is not their fixed point, and throughput then fails to
% balance: with the variance held at zero for the mean solve, Tput came back
% 2.000000 at the delay against 1.949745 at the queue on Delay ->
% Queue(PS,c=2), N=6, a 2.5% gap in a closed cycle where the two must be
% equal. Input and output variance differ only within the outer tolerance once
% the fixed point has converged.
%
% The delay stations have no min() to close, so their variance must not enter
% the drift; keep it for reporting only. NODRIFTVAR already masked SIGMA2SOLVE
% on the way in, so this IS the variance the mean solve used, which is what the
% paragraph above requires.
sigma2drift = sigma2solve;
covdrift = covsolve;

refinement = [];
if any(strcmp(method, {'refined','fluid.refined'}))
    % The refinement is a truncated expansion about the MEAN-FIELD fixed
    % point, not about the Gaussian one. Adding it to the 'minnormal' fixed
    % point would count the same O(1/N) term twice: the Gaussian closure
    % already resums it, since expanding E[F(X)] to second order and setting
    % it to zero reproduces exactly the Gast correction equation. So the base
    % point is recomputed with the first-order closure, while the Hessian and
    % the Jacobian are taken from the smooth Gaussian drift, the hard min
    % being only piecewise linear and, at saturation, kinked exactly at the
    % fixed point.
    mfopt = options;
    mfopt.method = 'closing';
    mfopt.config.moment_sigma2 = zeros(M,1);
    mfopt.config.moment_cov = cell(M,1);
    sigma2mf = zeros(M,1);
    [~, xvec_mf, ~, ~, ~, ~, mf_iters] = solver_fluid(sn, mfopt);
    iters = iters + mf_iters;
    xmf = xvec_mf{end}(:);

    A = terms.jacFcn(xmf, sigma2drift, covdrift);
    Sigma = local_lyapunov(A, terms.ratesFcn(xmf, sigma2drift, covdrift), terms);
    % A LINEAR DRIFT NEEDS NO REFINEMENT, and that is not the degenerate call
    % FLUID_REFINE_MEANFIELD refuses. When every station is either an infinite
    % server or minExact -- min(n,c) is the identity on the reachable set, the
    % population bound never reaching c -- the drift is exactly affine there,
    % its Hessian vanishes and the O(1/N) correction is identically zero. The
    % mask above then zeroes all of sigma2drift, which the refinement reads as
    % "the caller handed me the first-order closure" and rejects. Settle it
    % here, where the reason for the zero is known: a null correction, not an
    % error. Delay + PS(c=2) at N=2 is the smallest case.
    if all(noDriftVar(:) | terms.minExact(:))
        refinement = zeros(size(xmf));
    else
        refinement = fluid_refine_meanfield(xmf, sigma2drift, Sigma, terms, covdrift);
    end
    x = xmf + refinement;
    x(x < 0) = 0;
    % the corrected point is a correction OF the mean-field fixed point, so
    % its rates are read with the mean-field (zero) variance
    sigma2drift = sigma2mf;
    covdrift = cell(M,1);
    r = terms.ratesFcn(x, sigma2drift, covdrift);
    for i = 1:M
        blk = terms.stationBlock{i};
        if ~isempty(blk)
            sigma2(i) = max(0, sum(sum(Sigma(blk,blk))));
        end
    end
else
    r = terms.ratesFcn(x, sigma2drift, covdrift);
end

% -- performance measures, read back from the event representation -----------
gfac = terms.factorFcn(x, sigma2drift, covdrift);

% a load-dependent station clears alpha(n) times the nominal work, so its
% utilization normalises by the peak scaling (T*S/peak, as in the CTMC)
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
        % Summed over ORIGINAL events through EMAP: a reduced event folded
        % through an immediate coordinate is a completion at more than one
        % (station,class), and its rate has to reach every one of them.
        TN(i,k) = r(:)' * (terms.Emap * double(terms.evIsDeparture & terms.evStation == i & terms.evClass == k));
    end
end
% See solver_fluid_closing.m: TN is zero only to the integrator's accuracy.
RN(TN > GlobalConstants.Zero) = QN(TN > GlobalConstants.Zero) ./ TN(TN > GlobalConstants.Zero);

% -- transients, evaluated on the trajectory the ODE actually integrated -----
nt = size(xvec_t,1);
QNt = cell(M,K); UNt = cell(M,K); TNt = cell(M,K);
Gt = zeros(nt, terms.nstate);
Rt = zeros(nt, numel(terms.rateBase));
for s = 1:nt
    xs = xvec_t(s,:)';
    Gt(s,:) = terms.factorFcn(xs, sigma2drift, covdrift)';
    Rt(s,:) = terms.ratesFcn(xs, sigma2drift, covdrift)';
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
        TNt{i,k} = Rt * (terms.Emap * double(terms.evIsDeparture & terms.evStation == i & terms.evClass == k));
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
% The same variance as it entered the DRIFT (zero at the delay stations). A
% later solve on this fixed point, i.e. the passage-time ODE of
% SOLVER_FLUID_PASSAGE_TIME, must close its capacity term at this one or it
% evaluates a first-order drift at a second-order fixed point.
moments.sigma2Drift = sigma2drift;
moments.refinement = refinement;
moments.outerIters = outerIters;
moments.stationBlock = terms.stationBlock;
moments.classBlock = terms.classBlock;

xvec_it{end} = x(:)';
runtime = toc(T0);
end

function Sigma = local_lyapunov(A, r, terms)
% Stationary covariance on the coordinates that carry a real population.
%
% For a closed model covIdx is every coordinate and this is the plain solve.
% For an open or mixed model it drops the EXT source pool, whose coordinate is
% a normalisation constant rather than a job count; see FLUID_MOMENT_TERMS for
% why the projection leaves exactly the right open event set. The result is
% scattered back to full size with zeros on the dropped rows so that
% stationBlock/classBlock indexing is unchanged downstream.
idx = terms.covIdx;
Dc = terms.D(idx,:);
Qc = Dc * diag(r) * Dc';
Sc = fluid_lyapunov(A(idx,idx), Qc, Dc);
Sigma = zeros(terms.nstate);
Sigma(idx,idx) = Sc;
end

function cov = local_blend_cov(covA, covB, step)
% COV = LOCAL_BLEND_COV(COVA, COVB, STEP)
%
% covA + step*(covB - covA), entry by entry, with an empty cell read as the
% zero matrix and an entry left empty when both sides are. The twin of the
% scalar sigma2 blend the damped step takes, so the two stay consistent.
cov = cell(size(covB));
for i = 1:numel(covB)
    a = [];
    if i <= numel(covA), a = covA{i}; end
    b = covB{i};
    if isempty(a) && isempty(b)
        cov{i} = [];
    elseif isempty(a)
        cov{i} = step*b;
    elseif isempty(b)
        cov{i} = (1-step)*a;
    else
        cov{i} = a + step*(b - a);
    end
end
end

function idx = local_kink_stations(terms, x, sigma2)
% IDX = LOCAL_KINK_STATIONS(TERMS, X, SIGMA2)
%
% Every station whose population sits ON the saturation kink n_i = c_i of the
% first-order rate factor, in increasing order, empty when none does. Only the
% branches that take the indicator derivative can sit on one: a positive sigma2
% or a load-dependent row makes the closure smooth, and an infinite server never
% saturates. Twin of FluidRateFactors.driftKinkStations in the JAR.
idx = [];
if any(sigma2 > 0)
    return % the Gaussian closure is smooth, it has no kink
end
tol = sqrt(eps);
for i = 1:terms.M
    if terms.sched(i) == SchedStrategy.INF || terms.sched(i) == SchedStrategy.EXT
        continue
    end
    if ~isempty(terms.lldscaling) && i <= size(terms.lldscaling,1) ...
            && any(terms.lldscaling(i,:) ~= 1)
        continue % psi is piecewise quadratic and closed smoothly
    end
    c = terms.S(i);
    if isinf(c) || c <= 0
        continue
    end
    blk = terms.stationBlock{i};
    if isempty(blk)
        continue
    end
    ni = sum(x(blk));
    if ni <= 0
        continue % g = x on an empty station, no saturation term
    end
    if abs(ni - c) <= tol*max(1, c)
        idx(end+1) = i; %#ok<AGROW>
    end
end
end

function y = local_nudge_off_kink(terms, x, kink, rel)
% Y = LOCAL_NUDGE_OFF_KINK(TERMS, X, KINK, REL)
%
% A copy of X with every station in KINK moved to c_i*(1+REL), i.e. strictly
% onto one side of its kink. The station's coordinates are scaled together, so
% the phase mix and every other station are untouched.
y = x;
for t = 1:numel(kink)
    i = kink(t);
    blk = terms.stationBlock{i};
    ni = sum(y(blk));
    if ~(ni > 0)
        continue
    end
    y(blk) = y(blk) * (terms.S(i)*(1 + rel)/ni);
end
end
