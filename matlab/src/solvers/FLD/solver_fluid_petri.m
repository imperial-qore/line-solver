function [QN, UN, RN, TN, xvec_it, QNt, UNt, TNt, xvec_t, t, iters, runtime, moments] = solver_fluid_petri(sn, options)
% [QN, UN, RN, TN, XVEC_IT, QNT, UNT, TNT, XVEC_T, T, ITERS, RUNTIME, MOMENTS] = SOLVER_FLUID_PETRI(SN, OPTIONS)
%
% Fluid (mean-field) analysis of a stochastic Petri net, as the differential-
% algebraic system of options.method='dae'.
%
% A GSPN IS A DENSITY-DEPENDENT MARKOV POPULATION PROCESS, so the min-normal
% closure SolverFLD already carries applies to it unchanged: the marking is the
% population, a transition mode is a reaction, its incidence column is the jump,
% and lambda*min(enabling degree, servers) is the same min() the closure exists
% to smooth. What this analyzer adds is not a new closure but the three things a
% Petri net needs that a queueing network does not, and each of them is a reason
% the DAE form is the only one that can host it:
%
%   CONSERVATION IS A P-INVARIANT, not a chain population. A net's conserved
%   quantities are the left null vectors of its jump matrix, and stating them as
%   EQUATIONS is what makes them hold to solver tolerance instead of to
%   integrator tolerance -- and what supplies the rank the drift Jacobian is
%   missing along exactly those directions.
%
%   AN IMMEDIATE TRANSITION IS AN ALGEBRAIC FLOW. In the fluid limit a zero-time
%   transition fires infinitely fast, so a marking that enables it cannot
%   persist and what survives is a flow phi >= 0 pinned by the constraint that
%   the binding input place holds no mass. An ODE has nowhere to put that; the
%   active-set loop this method already runs for capacity caps has exactly the
%   right shape for it.
%
%   A BOUNDED PLACE IS A LINEAR INEQUALITY, throttling the deposit leg of every
%   firing that would overfill it, which is the loss semantics the exact engines
%   apply (the NRM's APPLYPLACECAPS).
%
% THE UNKNOWNS are u = [x; s2; phi; zeta]:
%
%   x     the fluid state: the marking, followed by the servers each
%         multi-phase mode has running in each phase (see FLUID_PETRI_TERMS)
%   s2    the Sigma entries the closure reads, one per TERMS.covPairs row,
%         eliminated by one Lyapunov solve per residual evaluation so that the
%         full covariance never enters the Newton vector
%   phi   one firing flow per immediate mode
%   mu    one server-latch flow per multi-phase mode, free in sign
%   zeta  one admission fraction per binding place capacity
%
% and the stacked residual is
%
%   0 = Dneg*r + gain(zeta).*(Dpos*r)   drift, one row per coordinate
%   0 = C*x - N                         conservation, one row per invariant
%   0 = s2 - Sigma(covPairs)            closure consistency
%   0 = x_b            (pins)           immediate: the binding place is empty
%   0 = phi_j*w_l - phi_l*w_j (ratios)  immediate: the GSPN conflict rule
%   0 = sum_h y(j,h) - theta_j(x)       the server latch of a multi-phase mode
%   0 = A_c*x - b_c                     each binding place capacity
%
% solved by the same damped projected Newton the queueing DAE uses.
%
% WHAT THE SECOND MOMENT IS. The covariance is the linear noise approximation of
% the WHOLE state, dSigma = A Sigma + Sigma A' + D diag(r) D' over the STOCHASTIC
% event columns only -- the in-flight firings of a multi-phase mode are a
% population like any other, and dropping them severs the only path by which the
% marking reaches that mode's firing rate. An immediate flow and a server latch
% carry no noise of their own: both are the limit of an infinitely fast
% mechanism whose fluctuation is slaved, so their columns are excluded and the
% fluctuation is reduced onto the manifold their constraints define
% (LOCAL_CLAMP_TANGENT) -- ORTHOGONALLY for a latch and a capacity, and
% OBLIQUELY for an immediate pin, along the fast column itself, so that a
% deposit into the pinned place is forwarded rather than deleted.
%
% Parameters:
%   sn      - NetworkStruct holding Place and Transition nodes
%   options - solver options; options.tol sets the Newton tolerance and
%             options.config.dae_maxstate caps the simultaneous solve
%
% Returns:
%   as SOLVER_FLUID_DAE, with MOMENTS.petri carrying the marking covariance,
%   the per-mode firing flows and the conservation residuals
%
% See also SOLVER_FLUID_DAE, FLUID_PETRI_TERMS, FLUID_PETRI_IMMEDIATE,
% FLUID_PETRI_CONSERVATION, FLUID_PETRI_CONSTRAINTS.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
M = sn.nstations;
K = sn.nclasses;

[appl, why] = fluid_petri_applicable(sn, options);
if ~appl
    line_error(mfilename, sprintf('The fluid Petri route cannot solve this model: %s.', why));
end

terms = fluid_petri_terms(sn, options);
n = terms.nstate;
npair = terms.npair;

maxstate = 100;
if isfield(options,'config') && isfield(options.config,'dae_maxstate') && ~isempty(options.config.dae_maxstate)
    maxstate = options.config.dae_maxstate;
end
if n > maxstate
    line_error(mfilename, sprintf(['The fluid Petri route solves a %d-unknown algebraic system with a ' ...
        'finite-difference Jacobian, above the limit of %d set by options.config.dae_maxstate. Raise that ' ...
        'limit, or use SolverSSA for a net of this size.'], n, maxstate));
end

cons = fluid_petri_conservation(terms);
if cons.leak > sqrt(GlobalConstants.Zero)
    line_error(mfilename, sprintf(['The conserved directions and the jump matrix disagree: the largest leak ' ...
        'per unit rate is %g, where it must be zero.'], cons.leak));
end
con = fluid_petri_constraints(sn, terms);
ncon = numel(con.b);

% -- seed --------------------------------------------------------------------
% Newton needs a point in the basin, not an answer. One first-order integration
% supplies it: the same drift at zero variance, with the immediate modes given a
% large FINITE rate so that a vanishing place drains along the trajectory
% instead of by an algebraic pin. The stiffness that buys is confined to the
% seed -- the answer is the algebraic solve, which carries no large rate at all.
[tseed, xseed, seedLam] = local_seed(terms, options);
x = xseed(end,:)';

imm = fluid_petri_immediate(terms, x);
nimm = imm.n;

% THE VARIANCE IS SEEDED POSITIVE, for the reason SOLVER_FLUID_DAE seeds it
% positive: sigma2 = 0 is where min() has no derivative, and a saturated net's
% first-order fixed point sits exactly there. The coordinate mean is an O(N)
% starting value in the right units, the variance a Poisson population of that
% mean would have.
s2 = zeros(npair,1);
ondiag = terms.covPairs(:,1) == terms.covPairs(:,2);
for p = find(ondiag)'
    s2(p) = max(GlobalConstants.FineTol, x(terms.covPairs(p,1)));
end

% THE NEWTON TOLERANCE IS THE FINE ONE BY DEFAULT, not options.tol's general
% 1e-4. The conservation rows are LINEAR in the unknowns, so the residual norm
% IS the token-count error: stopping at 1e-4 would leave a closed net holding
% 4.0001 tokens and undercut the one property this formulation exists to
% guarantee. A caller who asks for something tighter still gets it.
tol = GlobalConstants.FineTol;
if isfield(options,'tol') && ~isempty(options.tol) && isfinite(options.tol) && options.tol < tol
    tol = options.tol;
end
newton_max = 50;
if isfield(options,'iter_max') && ~isempty(options.iter_max)
    newton_max = max(newton_max, options.iter_max);
end

% -- steady state: one simultaneous algebraic solve per active set -----------
active = zeros(0,1);
phi = zeros(nimm,1);
nlatch = numel(terms.latchMode);
mu = zeros(nlatch,1);
zeta = zeros(0,1);
iters = 0;
aset_max = max(6, 2*(ncon + nimm) + 2);
converged = false; resnorm = Inf;
best = [];
for aset = 1:aset_max
    ctx = local_context(terms, cons, con, imm, active);
    u0 = [x; s2; phi; mu; ones(numel(active),1)];
    lb = [-inf(n,1); -inf(npair,1); zeros(nimm,1); -inf(nlatch,1); zeros(numel(active),1)];
    lb(n + find(ondiag)) = 0;
    resid = @(uu, quiet) local_residual(uu, ctx, quiet);
    [u, nit, converged, resnorm] = fluid_dae_newton(resid, u0, tol, newton_max, lb);
    iters = iters + nit;
    [x, s2, phi, mu, zeta] = local_unpack(u, terms, imm, active);

    if isempty(best) || resnorm < best.resnorm
        best = struct('x',x,'s2',s2,'phi',phi,'mu',mu,'zeta',zeta,'active',active,'imm',imm, ...
            'resnorm',resnorm,'converged',converged);
    end

    % -- the active-set moves, in the order that a failure of each invalidates
    % the next: a negative flow means the mode does not fire at all, a negative
    % marking means the wrong coordinate was pinned, and only then is it worth
    % asking which capacity rows bind.
    moved = false;
    k = [];
    if nimm > 0
        k = find(phi < -max(tol, 1e-10)*max(1, max(abs(phi))), 1);
    end
    if ~isempty(k)
        imm.active(k) = false;
        imm = fluid_petri_immediate(terms, x, imm);
        moved = true;
    end
    if ~moved && nimm > 0
        s = find(x(1:terms.nm) < -max(tol, 1e-10), 1);
        if ~isempty(s)
            cand = find(imm.active(:).' & arrayfun(@(kk) any(terms.modes(terms.immIdx(kk)).arcSlot == s) ...
                && imm.bind(kk) ~= s, 1:nimm));
            if ~isempty(cand)
                imm.bind(cand(1)) = s;
                imm = fluid_petri_immediate(terms, x, imm);
                moved = true;
            end
        end
    end
    if ~moved && ncon > 0
        val = con.A * x;
        over = find(val > con.b + max(1e-9, tol));
        over = setdiff(over, active);
        if ~isempty(over)
            active = [active; over(:)]; %#ok<AGROW>
            moved = true;
        else
            rel = find(zeta > 1 + max(1e-9, tol));
            if ~isempty(rel)
                active(rel) = [];
                moved = true;
            end
        end
    end
    if ~moved
        break
    end
end

if ~isempty(best) && ~converged && best.converged
    x = best.x; s2 = best.s2; phi = best.phi; mu = best.mu; zeta = best.zeta;
    active = best.active; imm = best.imm;
    converged = true; resnorm = best.resnorm;
end
if ~converged
    line_warning(mfilename, sprintf(['The simultaneous closure solve stopped at residual %.3e after %d ' ...
        'Newton steps without reaching %.3e. The reported point is the last iterate.'], resnorm, iters, tol));
end
% A NEGATIVE MARKING IS NOT ROUNDING. It means the fixed point wanted mass a
% place cannot supply and no immediate mode could be rebound to pin it, so the
% answer is outside the model's own state space and is reported as such rather
% than clipped into it.
neg = find(x(1:terms.nm) < -max(tol, 1e-10));
if ~isempty(neg)
    line_warning(mfilename, sprintf(['The fixed point holds %g tokens at %s, which is negative: no ' ...
        'immediate transition could be rebound to pin that place at zero. Use SolverCTMC, SolverSSA or ' ...
        'SolverLDES for this net.'], x(neg(1)), terms.namesNode{terms.coordNode(neg(1))}));
end

ctx = local_context(terms, cons, con, imm, active);
[~, r, ~, Sigma] = local_residual([x; s2; phi; mu; zeta], ctx, false);

% -- transient ---------------------------------------------------------------
Sigmat = [];
tvar = [];
t = tseed;
xvec_t = xseed;
% The immediate flow ALONG the reported trajectory, which is not the
% steady-state one: on the seed path it is the large-finite-rate approximation
% the seed itself integrated, so the throughput table and the trajectory it is
% read off describe the same run.
phit = local_seed_flows(terms, xseed, seedLam);
if isfinite(options.timespan(2))
    [t, xvec_t, phit, Sigmat, tvar] = local_transient(terms, cons, con, imm, active, options, x, s2, phi, zeta);
end

% -- performance measures ----------------------------------------------------
[QN, UN, RN, TN] = local_metrics(terms, x, r, M, K);
[QNt, UNt, TNt] = local_metrics_t(terms, xvec_t, s2, phit, M, K);

xvec_it = {x};
moments = struct();
moments.Sigma = Sigma;
moments.petri = local_report(terms, cons, con, imm, active, x, r, phi, zeta, Sigma);
moments.QVar = local_qvar(terms, Sigma, M, K);
moments.QStd = sqrt(max(0, moments.QVar));
moments.Sigmat = Sigmat;
moments.tvar = tvar;
moments.method = 'dae';
moments.resnorm = resnorm;
runtime = toc(T0);
end

% ---------------------------------------------------------------------------
function ctx = local_context(terms, cons, con, imm, active)
% Everything the residual needs that does not change within one active set.
%
% THE CONSERVATION ROWS A BINDING CAP BREAKS ARE DROPPED. A capped place loses
% the tokens that do not fit, so a conserved quantity supported on it is not
% conserved while the cap binds; keeping its row would state an equation the
% drift contradicts, which is a singular Newton system rather than an
% inaccuracy. An open net has no such row to lose -- its arrivals already broke
% them -- which is where a bounded place is actually declared.
ctx = struct();
ctx.terms = terms;
ctx.imm = imm;
ctx.active = active(:);
ctx.con = con;
C = cons.C; N = cons.N;
if ~isempty(active) && ~isempty(C)
    hit = any(con.cover(active,:), 1);
    drop = any(C(:, hit) ~= 0, 2);
    C = C(~drop, :);
    N = N(~drop);
end
ctx.C = C; ctx.N = N;
ctx.Dp = max(terms.D, 0);
ctx.Dn = min(terms.D, 0);
end

% ---------------------------------------------------------------------------
function [x, s2, phi, mu, zeta] = local_unpack(u, terms, imm, active)
n = terms.nstate; np = terms.npair; ni = imm.n; na = numel(active);
nl = numel(terms.latchMode);
x = u(1:n);
s2 = u(n + (1:np));
phi = max(0, u(n + np + (1:ni)));
mu = u(n + np + ni + (1:nl));
zeta = max(0, u(n + np + ni + nl + (1:na)));
end

% ---------------------------------------------------------------------------
function T = local_clamp_tangent(terms, imm, con, active, th)
% The reduction of the fluctuation onto the manifold the fast and clamped
% directions leave free, over the whole state.
%
% AN IMMEDIATE PIN REDUCES OBLIQUELY, ALONG THE FAST REACTION ITSELF, and this
% is the one place where the orthogonal projector the queueing twin uses is
% wrong rather than merely different. A slow event that deposits into a pinned
% place is answered INSTANTLY by the immediate transition, so its effective jump
% is its own plus the immediate flow it triggers -- the token is forwarded, not
% lost. An orthogonal projection deletes the deposit instead, which destroys
% mass in the diffusion: on P1 -> P2 -> (immediate) -> P3 -> P1 with four tokens
% it reported Var[P1] = 3.88 against Var[P3] = 0.45 on a net whose two free
% places are exchangeable, and dragged the MEAN with it (3.26/0.74 against the
% 2.75/1.25 of the same net with the immediate transition eliminated by hand).
% With the oblique reduction the two answers agree to machine precision, which
% is the test that says the reduction is the right one.
%
%   P = I - Cf * G * E_B,   G = G0 * (E_B Cf G0)^-1
%
% with Cf the incidence columns of the active immediate modes, E_B the selector
% of the pinned coordinates, and G0 the weight-normalised group membership -- so
% the correction flow follows the SAME firing-weight split as the mean flow.
%
% A CAPACITY CAP REDUCES ORTHOGONALLY, as in the queueing twin: a place held at
% its capacity fixes that combination of the state and the mass that does not
% fit is genuinely lost, so there is nothing to forward it to.
%
% A SERVER LATCH REDUCES ORTHOGONALLY TOO, and it is why the phase coordinates
% can stay in the covariance at all: `sum_h y(j,h) - theta_j(m)` is identically
% zero in the exact process, because ENABLE answers the marking with no delay,
% so that combination carries no fluctuation. The row is the LINEARISED
% constraint, `[-dtheta_j/dm , 1 over the phase block]`, which is why this
% projector depends on the iterate and is rebuilt at every residual evaluation
% rather than once per active set.
idx = terms.covIdx;
nc = numel(idx);
T = eye(nc);

actk = find(imm.active(:).' & imm.bind(:).' > 0);
B = imm.pins(:).';
if ~isempty(actk) && ~isempty(B)
    Cf = zeros(nc, numel(actk));
    for a = 1:numel(actk)
        Cf(:,a) = terms.modes(terms.immIdx(actk(a))).Cvec(idx);
    end
    G0 = zeros(numel(actk), numel(B));
    for jb = 1:numel(B)
        grp = find(imm.bind(actk) == B(jb));
        w = zeros(numel(grp),1);
        for q = 1:numel(grp)
            w(q) = terms.modes(terms.immIdx(actk(grp(q)))).weight;
        end
        if sum(w) <= 0
            w = ones(numel(grp),1);
        end
        G0(grp, jb) = w/sum(w);
    end
    EB = zeros(numel(B), nc);
    for jb = 1:numel(B)
        EB(jb, B(jb)) = 1;
    end
    Mb = (EB*Cf)*G0;
    T = T - Cf*(G0*pinv(Mb))*EB;
end

R = zeros(0, nc);
for c = active(:).'
    R(end+1,:) = con.A(c, idx); %#ok<AGROW>
end
if nargin >= 5 && ~isempty(th)
    for j = terms.latchMode
        row = zeros(1, nc);
        row(terms.modes(j).zblk) = 1;
        row(th.dslot{j}) = row(th.dslot{j}) - th.dval{j}.';
        R(end+1,:) = row; %#ok<AGROW>
    end
end
if ~isempty(R) && any(abs(R(:)) > GlobalConstants.Zero)
    T = (eye(nc) - R' * pinv(R * R') * R) * T;
end
if norm(T - eye(nc), inf) <= GlobalConstants.Zero
    T = [];
end
end

% ---------------------------------------------------------------------------
function [G, rout, th, Sigma] = local_residual(u, ctx, quiet)
% The coupled algebraic system, stacked. Returns [] when the closure cannot be
% evaluated at this iterate so that the line search can back off; the first
% evaluation of a pass runs with QUIET false, where a genuine failure surfaces.
terms = ctx.terms;
imm = ctx.imm;
n = terms.nstate;
[x, s2, phi, mu, zeta] = local_unpack(u, terms, imm, ctx.active);

G = []; rout = []; Sigma = [];
th = fluid_petri_theta(terms, x, s2);
r = fluid_petri_rates(terms, x, s2, phi, mu, th);

% the deposit gate of every binding capacity, as a product of fractions
gain = ones(n,1);
for k = 1:numel(ctx.active)
    c = ctx.active(k);
    gain(ctx.con.cover(c,:)) = gain(ctx.con.cover(c,:)) * zeta(k);
end
drift = ctx.Dn*r + gain .* (ctx.Dp*r);

try
    A = fluid_petri_jacobian(terms, x, s2, phi, th);
    Sigma = local_sigma(terms, A, r, local_clamp_tangent(terms, imm, ctx.con, ctx.active, th));
catch me
    if ~quiet
        rethrow(me);
    end
    G = [];
    return
end
s2new = zeros(terms.npair,1);
for p = 1:terms.npair
    s2new(p) = Sigma(terms.covPairs(p,1), terms.covPairs(p,2));
end

G = drift;
if ~isempty(ctx.C)
    G = [G; ctx.C*x - ctx.N];
end
G = [G; s2 - s2new];
for q = 1:numel(imm.rows)
    row = imm.rows(q);
    switch row.kind
        case 'pin'
            G(end+1,1) = x(row.a); %#ok<AGROW>
        case 'ratio'
            G(end+1,1) = phi(row.a)*row.wb - phi(row.b)*row.wa; %#ok<AGROW>
        case 'zero'
            G(end+1,1) = phi(row.a); %#ok<AGROW>
    end
end
for q = 1:numel(terms.latchMode)
    j = terms.latchMode(q);
    G(end+1,1) = sum(x(terms.modes(j).zblk)) - th.theta(j); %#ok<AGROW>
end
for k = 1:numel(ctx.active)
    c = ctx.active(k);
    G(end+1,1) = ctx.con.A(c,:)*x - ctx.con.b(c); %#ok<AGROW>
end
rout = r;
end

% ---------------------------------------------------------------------------
function Sigma = local_sigma(terms, A, r, clampT)
% One Lyapunov solve, over the marking coordinates.
%
% THE DIFFUSION COUNTS THE STOCHASTIC EVENTS ONLY. An immediate flow is not a
% Poisson stream with an intensity; it is the limit of an infinitely fast one
% whose fluctuation is slaved to the marking it drains, and its pinned
% coordinate is projected out by CLAMPT rather than given a noise intensity.
%
% THE PHASE BLOCK IS OUTSIDE THE SOLVE, which is an approximation and is named
% as one: restricting A to the marking block freezes the phase distribution of a
% multi-phase mode at its stationary value, so the reported covariance is that
% of the marking driven by a renewal-blind firing stream. The MEAN is unaffected
% -- the phase coordinates are solved exactly, in the drift.
idx = terms.covIdx;
Dc = terms.D(idx, terms.stochCol);
Am = A(idx,idx);
if ~isempty(clampT)
    % BOTH the jump directions and the generator are reduced. FLUID_LYAPUNOV
    % restricts everything to range(Dc) already, so reducing Dc alone would fix
    % the subspace but leave the generator's ORTHOGONAL component on it, which
    % is not the reduced dynamics when the reduction is oblique.
    Dc = clampT * Dc;
    Am = clampT * Am;
end
Q = Dc * diag(r(terms.stochCol)) * Dc.';
Sc = fluid_lyapunov(Am, Q, Dc);
Sigma = zeros(terms.nstate);
Sigma(idx,idx) = Sc;
end

% ---------------------------------------------------------------------------
function [QN, UN, RN, TN] = local_metrics(terms, x, r, M, K)
% The station table, in the conventions the exact SPN engines report.
%
% A PLACE IS AN INF STATION: its queue length is its mean token count and its
% utilization is the same number, which is what SolverCTMC and the NRM both
% report for a Place. Its throughput is the rate at which TOKENS leave it, so a
% consuming mode contributes its firing rate times the multiplicity of the arc:
% that is what SolverCTMC reports and what Little's law needs (a place drained
% through an arc of weight two reports 3.158 where the firing rate is 1.579).
% SOLVER_SSA'S NRM SUMS THE UNWEIGHTED PROPENSITY instead, so the two disagree
% wherever an input arc has a multiplicity above one; the exact engine is the
% reference. A Source reports the arrival rate it injects, which is the
% reference-station throughput of its open class.
QN = zeros(M,K); UN = zeros(M,K); RN = zeros(M,K); TN = zeros(M,K);
for s = 1:terms.nm
    ist = terms.coordStation(s);
    k = terms.coordClass(s);
    QN(ist,k) = QN(ist,k) + x(s);
    UN(ist,k) = QN(ist,k);
end
for ist = 1:M
    for k = 1:K
        e = terms.consumers{ist,k};
        if ~isempty(e)
            w = terms.consumerW{ist,k};
            TN(ist,k) = TN(ist,k) + sum(w(:) .* r(e(:)));
        end
        e = terms.producers{ist,k};
        if ~isempty(e)
            TN(ist,k) = TN(ist,k) + sum(r(e));
        end
    end
end
RN(TN > 0) = QN(TN > 0) ./ TN(TN > 0);
end

% ---------------------------------------------------------------------------
function [QNt, UNt, TNt] = local_metrics_t(terms, xt, s2, phit, M, K)
% The same reader along a trajectory. The rates are evaluated at the STATIONARY
% closure variance, which is what every fluid method other than the transient
% covariance branch does, and at the immediate flow the trajectory itself
% carried at each instant.
nt = size(xt,1);
QNt = cell(M,K); UNt = cell(M,K); TNt = cell(M,K);
for ist = 1:M
    for k = 1:K
        QNt{ist,k} = zeros(nt,1); UNt{ist,k} = zeros(nt,1); TNt{ist,k} = zeros(nt,1);
    end
end
Rt = zeros(nt, terms.nev);
for a = 1:nt
    % the latch flow moves no marking, so it cannot change a place throughput
    Rt(a,:) = fluid_petri_rates(terms, xt(a,:)', s2, phit(a,:)', []).';
end
for s = 1:terms.nm
    ist = terms.coordStation(s);
    k = terms.coordClass(s);
    QNt{ist,k} = QNt{ist,k} + xt(:,s);
    UNt{ist,k} = QNt{ist,k};
end
for ist = 1:M
    for k = 1:K
        e = terms.consumers{ist,k};
        if ~isempty(e)
            TNt{ist,k} = TNt{ist,k} + Rt(:,e) * terms.consumerW{ist,k}(:);
        end
        e = terms.producers{ist,k};
        if ~isempty(e)
            TNt{ist,k} = TNt{ist,k} + sum(Rt(:,e), 2);
        end
    end
end
end

% ---------------------------------------------------------------------------
function QVar = local_qvar(terms, Sigma, M, K) %#ok<INUSD>
% The variance of each station-class token count, read off the marking
% covariance.
QVar = zeros(M,K);
for s = 1:terms.nm
    ist = terms.coordStation(s);
    k = terms.coordClass(s);
    QVar(ist,k) = max(0, Sigma(s,s));
end
end

% ---------------------------------------------------------------------------
function rep = local_report(terms, cons, con, imm, active, x, r, phi, zeta, Sigma)
% What the Petri route computes that the station table has no column for: the
% marking itself, the per-mode firing flows, the conserved quantities and how
% exactly they held, and which constraints ended up binding.
rep = struct();
rep.marking = zeros(terms.I, terms.K);
for s = 1:terms.nm
    rep.marking(terms.coordNode(s), terms.coordClass(s)) = x(s);
end
rep.markingVar = zeros(terms.I, terms.K);
for s = 1:terms.nm
    rep.markingVar(terms.coordNode(s), terms.coordClass(s)) = max(0, Sigma(s,s));
end
nmod = numel(terms.modes);
rep.modeLabel = cell(nmod,1);
rep.modeFlow = zeros(nmod,1);
for j = 1:nmod
    rep.modeLabel{j} = terms.modes(j).label;
    e = find(terms.evMode == j & (terms.evKind == 1 | terms.evKind == 4));
    if ~isempty(e)
        rep.modeFlow(j) = sum(r(e));
    end
end
rep.immediateFlow = phi;
rep.invariantLabel = cons.label;
rep.invariantValue = cons.N;
rep.invariantError = cons.C*x - cons.N;
rep.capacityLabel = con.label;
rep.capacityActive = active(:);
rep.capacityFraction = zeta(:);
rep.pinned = imm.pins(:);
rep.Sigma = Sigma;
end

% ---------------------------------------------------------------------------
function [t, xt, Lam] = local_seed(terms, options)
% One first-order trajectory, integrated from the initial marking.
%
% THE IMMEDIATE MODES GET A LARGE FINITE RATE HERE, and only here. The seed has
% to reach a marking in the basin of the fixed point, and the vanishing places
% have to be empty there for the active set to be read off it -- so the flow is
% approximated by a fast reaction rather than by an algebraic pin. The rate is
% scaled to the model's own timescale instead of taken from
% GlobalConstants.Immediate: 1e8 against a rate of order one is a stiffness the
% seed does not need, and the answer does not depend on the seed's accuracy.
n = terms.nstate;
rmax = 0;
for j = terms.timedIdx
    rmax = max(rmax, sum(terms.modes(j).d1));
end
for e = find(terms.evKind == 3).'
    rmax = max(rmax, terms.rateBase(e));
end
if ~(rmax > 0)
    rmax = 1;
end
Lam = min(GlobalConstants.Immediate, 1e4*rmax);

mass = sum(terms.x0(1:terms.nm));
T = 50*(mass + 1)/rmax;
odeopt = odeset('AbsTol', GlobalConstants.FineTol, 'RelTol', GlobalConstants.CoarseTol, ...
    'NonNegative', 1:n);
t = 0; xt = terms.x0.';
for attempt = 1:6
    [t, xt] = ode_solve_stiff(@(tt,xx) local_seed_drift(terms, xx, Lam), [0 T], terms.x0, odeopt, options);
    d = local_seed_drift(terms, xt(end,:)', Lam);
    if norm(d, inf) <= GlobalConstants.CoarseTol*max(1, rmax*(mass+1))
        break
    end
    T = 4*T;
end
end

% ---------------------------------------------------------------------------
function d = local_seed_drift(terms, x, Lam)
% The first-order drift: the same rates at zero variance, with an immediate mode
% firing at LAM times its enabling degree and its firing weight, and the server
% latch relaxed at LAM towards the enabling degree instead of solved. Both are
% approximations of an algebraic constraint by a fast reaction, and both are
% confined to the seed -- the answer is the algebraic solve.
x = max(0, x);
th = fluid_petri_theta(terms, x, zeros(terms.npair,1));
phi = zeros(numel(terms.immIdx),1);
for k = 1:numel(terms.immIdx)
    j = terms.immIdx(k);
    phi(k) = Lam * terms.modes(j).weight * th.theta(j);
end
mu = zeros(numel(terms.latchMode),1);
for q = 1:numel(terms.latchMode)
    j = terms.latchMode(q);
    mu(q) = Lam * (th.theta(j) - sum(x(terms.modes(j).zblk)));
end
r = fluid_petri_rates(terms, x, [], phi, mu, th);
d = terms.D * r;
end

% ---------------------------------------------------------------------------
function [t, xt, phit, Sigmat, tvar] = local_transient(terms, cons, con, imm, active, options, xss, s2ss, phiss, zetass)
% The transient closure as an index-1 DAE with a singular mass matrix.
%
% THE INITIAL MARKING JUMPS. An immediate transition fires in zero time, so a
% marking that enables one is not a state the trajectory ever occupies: the
% fluid analogue of the vanishing-marking collapse moves the initial marking
% along the immediate incidence columns until every pinned place is empty, and
% the trajectory starts THERE. That jump conserves every P-invariant by
% construction, since a conserved direction annihilates the columns it moves
% along.
%
% THE STRUCTURE, per segment:
%
%   d/dt x_s = drift_s                    a coordinate nothing pins
%   0        = x_b                        a coordinate an immediate mode pins
%   0        = drift_b                    the DIFFERENTIATED pin, which is what
%                                         determines that mode's flow
%   0        = phi_j*w_l - phi_l*w_j       the conflict rule, unchanged
%   0        = d/dt(sum_h y(j,h) - theta_j) the DIFFERENTIATED server latch,
%                                         which is what determines mu_j
%   0        = A_c * drift                 the differentiated capacity constraint
%   d/dt Sigma = A Sigma + Sigma A' + Q    the covariance, when it fits
%
% The undifferentiated pin x_b = 0 and its derivative drift_b = 0 are BOTH
% present, and they are not redundant: the first makes x_b an algebraic
% variable at its pinned value, the second is the equation that pins the flow.
% Using only the undifferentiated form would be index 2, which neither ode15s
% nor RODAS solves.
%
% THE INTEGRATOR IS OPTIONS.ODESOLVERS.DAESOLVER, the same slot SOLVER_FLUID_DAE
% reads: ODE_SOLVE routes to ode113, which carries no mass matrix at all, and
% ODE_SOLVE_STIFF may route to ode23s, which carries no SINGULAR one.
n = terms.nstate;
nimm = imm.n;
nlatch = numel(terms.latchMode);
nact = numel(active);
tspan = options.timespan;
if ~isfinite(tspan(2))
    t = 0; xt = xss.'; phit = phiss(:).'; Sigmat = []; tvar = [];
    return
end

maxcov = 25;
if isfield(options,'config') && isfield(options.config,'dae_maxcov') && ~isempty(options.config.dae_maxcov)
    maxcov = options.config.dae_maxcov;
end
nc = numel(terms.covIdx);
withcov = nc > 0 && nc <= maxcov;
if ~withcov && nc > maxcov
    line_warning(mfilename, sprintf(['The transient covariance would add %d differential states, above the ' ...
        'limit of %d set by options.config.dae_maxcov. Integrating the mean as a DAE with the variance held ' ...
        'at its stationary value.'], nc*nc, maxcov*maxcov));
end

x0 = local_collapse(terms, imm, terms.x0);
solver = options.odesolvers.daeSolver;

t = zeros(0,1); xt = zeros(0,n); phit = zeros(0,nimm);
Sigmat = []; tvar = zeros(0,1);
tcur = tspan(1);
xcur = x0;
phi = phiss; mu = zeros(nlatch,1); zeta = zetass;
Scur = zeros(nc);
seg_max = max(8, 4*(nimm + numel(con.b) + 1));
for seg = 1:seg_max
    % THE ACTIVE SET IS RE-READ EVERY SEGMENT, not carried from the top: a
    % located crossing may have added or released a capacity row, which changes
    % how many algebraic unknowns the state vector holds and where the
    % covariance block starts.
    nact = numel(active);
    ctx = local_context(terms, cons, con, imm, active);
    nvar = n + nimm + nlatch + nact;
    Mm = zeros(nvar);
    for s = 1:n
        Mm(s,s) = 1;
    end
    for p = imm.pins(:).'
        Mm(p,p) = 0;
    end
    if withcov
        Mm(nvar + nc*nc, nvar + nc*nc) = 0; % grow
        Mm(nvar+1:nvar+nc*nc, nvar+1:nvar+nc*nc) = eye(nc*nc);
        nvar = nvar + nc*nc;
    end
    z0 = [xcur; phi(:); mu(:); zeta(:)];
    if withcov
        z0 = [z0; Scur(:)]; %#ok<AGROW>
    end
    odeopt = odeset('Mass', Mm, 'MassSingular', 'yes', 'AbsTol', GlobalConstants.FineTol, ...
        'RelTol', GlobalConstants.CoarseTol, 'Events', @(tt,zz) local_events(tt, zz, ctx, withcov, nc));
    [tsg, zsg, te, ~, ie] = feval(solver, @(tt,zz) local_dae_rhs(tt, zz, ctx, withcov, nc), ...
        [tcur tspan(2)], z0, odeopt);
    t = [t; tsg]; %#ok<AGROW>
    xt = [xt; zsg(:,1:n)]; %#ok<AGROW>
    phit = [phit; zsg(:, n+(1:nimm))]; %#ok<AGROW>
    if withcov
        tvar = [tvar; tsg]; %#ok<AGROW>
        blk = reshape(zsg(:, n+nimm+nlatch+nact+(1:nc*nc)).', nc, nc, []);
        if isempty(Sigmat)
            Sigmat = blk;
        else
            Sigmat = cat(3, Sigmat, blk);
        end
    end
    tcur = tsg(end);
    zend = zsg(end,:).';
    xcur = zend(1:n);
    phi = zend(n+(1:nimm));
    mu = zend(n+nimm+(1:nlatch));
    zeta = zend(n+nimm+nlatch+(1:nact));
    if withcov
        Scur = reshape(zend(n+nimm+nlatch+nact+(1:nc*nc)), nc, nc);
    end
    if isempty(te) || tcur >= tspan(2) - GlobalConstants.FineTol
        break
    end
    [imm, active, changed, zeta] = local_switch(terms, con, imm, active, ie(end), xcur, phi, zeta);
    if ~changed
        break
    end
    % A COORDINATE THAT WAS JUST PINNED IS NOT AT ZERO YET. The crossing is
    % located where it reaches zero, so the restart is consistent to the
    % integrator's own event tolerance; ode15s and RODAS both absorb that much
    % inconsistency in an index-1 system, and the algebraic row drives it to
    % zero within the first step of the next segment.
end
if isempty(t)
    t = tspan(1); xt = x0.'; phit = phiss(:).';
end
end

% ---------------------------------------------------------------------------
function phit = local_seed_flows(terms, xt, Lam)
% The immediate firing flows the SEED trajectory carried, point by point: the
% large finite rate times the enabling degree and the firing weight, which is
% the approximation the seed integrated and therefore the one its throughput
% table has to be read at.
nt = size(xt,1);
ni = numel(terms.immIdx);
phit = zeros(nt, ni);
if ni == 0
    return
end
for a = 1:nt
    th = fluid_petri_theta(terms, xt(a,:)', zeros(terms.npair,1));
    for k = 1:ni
        j = terms.immIdx(k);
        phit(a,k) = Lam * terms.modes(j).weight * th.theta(j);
    end
end
end

% ---------------------------------------------------------------------------
function x = local_collapse(terms, imm, x0)
% The fluid vanishing-marking collapse: move the initial marking along the
% immediate incidence columns until every pinned place is empty.
%
% The move is s >= 0 along the columns of the active immediate modes, chosen so
% that the pinned coordinates land on zero. Every conserved quantity survives it
% for free, since a conserved direction annihilates those columns.
x = x0;
act = find(imm.active(:).');
pins = imm.pins(:);
if isempty(act) || isempty(pins)
    return
end
Cimm = zeros(terms.nstate, numel(act));
for a = 1:numel(act)
    Cimm(:,a) = terms.modes(terms.immIdx(act(a))).Cvec;
end
b = -x0(pins);
Aimm = Cimm(pins,:);
warnstate = warning('off','MATLAB:rankDeficientMatrix');
s = Aimm \ b;
warning(warnstate);
if any(~isfinite(s))
    return
end
s = max(0, s);
x = x0 + Cimm*s;
x(pins) = 0;
x(1:terms.nm) = max(0, x(1:terms.nm));
% The collapse moved the marking, so the servers a multi-phase mode has latched
% moved with it: re-read the enabling degree at the collapsed marking, or the
% trajectory starts off its own latch and the differentiated constraint has no
% way to bring it back.
th = fluid_petri_theta(terms, x, zeros(terms.npair,1));
for j = terms.latchMode
    x(terms.modes(j).zblk) = th.theta(j) * terms.modes(j).pie(:);
end
end

% ---------------------------------------------------------------------------
function dz = local_dae_rhs(tt, z, ctx, withcov, nc) %#ok<INUSL>
% The right-hand side of one segment, in the layout LOCAL_TRANSIENT documents.
terms = ctx.terms;
imm = ctx.imm;
n = terms.nstate; nimm = imm.n; nact = numel(ctx.active);
nlatch = numel(terms.latchMode);
x = z(1:n);
phi = z(n+(1:nimm));
mu = z(n+nimm+(1:nlatch));
zeta = z(n+nimm+nlatch+(1:nact));
th = fluid_petri_theta(terms, x, zeros(terms.npair,1));
r = fluid_petri_rates(terms, x, [], phi, mu, th);
gain = ones(n,1);
for k = 1:nact
    c = ctx.active(k);
    gain(ctx.con.cover(c,:)) = gain(ctx.con.cover(c,:)) * zeta(k);
end
drift = ctx.Dn*r + gain .* (ctx.Dp*r);

dz = drift;
for p = imm.pins(:).'
    dz(p) = x(p); % the algebraic level pin
end
for q = 1:numel(imm.rows)
    row = imm.rows(q);
    switch row.kind
        case 'pin'
            dz(end+1,1) = drift(row.a); %#ok<AGROW>
        case 'ratio'
            dz(end+1,1) = phi(row.a)*row.wb - phi(row.b)*row.wa; %#ok<AGROW>
        case 'zero'
            dz(end+1,1) = phi(row.a); %#ok<AGROW>
    end
end
% The DIFFERENTIATED latch, one row per multi-phase mode: the running-server
% count tracks the enabling degree exactly, so their rates of change agree and
% that is the equation MU solves.
for q = 1:nlatch
    j = terms.latchMode(q);
    acc = sum(drift(terms.modes(j).zblk));
    if ~isempty(th.dslot{j})
        acc = acc - th.dval{j}.' * drift(th.dslot{j});
    end
    dz(end+1,1) = acc; %#ok<AGROW>
end
for k = 1:nact
    c = ctx.active(k);
    dz(end+1,1) = ctx.con.A(c,:)*drift; %#ok<AGROW>
end
if withcov
    idx = terms.covIdx;
    A = fluid_petri_jacobian(terms, x, zeros(terms.npair,1), phi, th);
    Dc = terms.D(idx, terms.stochCol);
    Am = A(idx,idx);
    clampT = local_clamp_tangent(terms, imm, ctx.con, ctx.active, th);
    if ~isempty(clampT)
        Dc = clampT * Dc;
        Am = clampT * Am;
    end
    Q = Dc * diag(r(terms.stochCol)) * Dc.';
    S = reshape(z(n+nimm+nlatch+nact+(1:nc*nc)), nc, nc);
    dS = Am*S + S*Am.' + Q;
    dz = [dz; dS(:)];
end
end

% ---------------------------------------------------------------------------
function [value, isterminal, direction] = local_events(tt, z, ctx, withcov, nc) %#ok<INUSD>
% Where the segment ends: a flow that would go negative, a capacity that starts
% or stops binding, or a marking coordinate that would go negative.
terms = ctx.terms;
imm = ctx.imm;
n = terms.nstate; nimm = imm.n; nact = numel(ctx.active);
nlatch = numel(terms.latchMode);
x = z(1:n);
phi = z(n+(1:nimm));
zeta = z(n+nimm+nlatch+(1:nact));
value = zeros(0,1); direction = zeros(0,1);
for k = 1:nimm
    if imm.active(k)
        value(end+1,1) = phi(k); direction(end+1,1) = -1; %#ok<AGROW>
    else
        value(end+1,1) = 1; direction(end+1,1) = 0; %#ok<AGROW>
    end
end
for c = 1:numel(ctx.con.b)
    k = find(ctx.active == c, 1);
    if isempty(k)
        value(end+1,1) = ctx.con.b(c) - ctx.con.A(c,:)*x; direction(end+1,1) = -1; %#ok<AGROW>
    else
        value(end+1,1) = 1 - zeta(k); direction(end+1,1) = -1; %#ok<AGROW>
    end
end
for s = 1:terms.nm
    if any(imm.pins == s)
        value(end+1,1) = 1; direction(end+1,1) = 0; %#ok<AGROW>
    else
        value(end+1,1) = x(s) + GlobalConstants.FineTol; direction(end+1,1) = -1; %#ok<AGROW>
    end
end
isterminal = ones(numel(value),1);
end

% ---------------------------------------------------------------------------
function [imm, active, changed, zeta] = local_switch(terms, con, imm, active, ievent, x, phi, zeta) %#ok<INUSL>
% The active-set move a located crossing asks for, and the multiplier vector
% resized to match it: a newly bound capacity starts unthrottled, a released one
% loses its unknown.
changed = false;
nimm = imm.n;
ncon = numel(con.b);
if ievent <= nimm
    if imm.active(ievent)
        imm.active(ievent) = false;
        imm = fluid_petri_immediate(terms, x, imm);
        changed = true;
    end
    return
end
ievent = ievent - nimm;
if ievent <= ncon
    k = find(active == ievent, 1);
    if isempty(k)
        active = [active; ievent];
        zeta = [zeta(:); 1];
    else
        active(k) = [];
        zeta(k) = [];
    end
    changed = true;
    return
end
s = ievent - ncon;
cand = find(arrayfun(@(kk) imm.active(kk) && any(terms.modes(terms.immIdx(kk)).arcSlot == s) ...
    && imm.bind(kk) ~= s, 1:nimm));
if ~isempty(cand)
    imm.bind(cand(1)) = s;
    imm = fluid_petri_immediate(terms, x, imm);
    changed = true;
end
end
