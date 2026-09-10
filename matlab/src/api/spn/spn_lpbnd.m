function bnd = spn_lpbnd(sn, options)
% BND = SPN_LPBND(SN)
% BND = SPN_LPBND(SN, OPTIONS)
% Linear-programming bounds on the mean marking and the throughputs of a
% stochastic timed Petri net.
%
% The stationary chain is relaxed to a MOMENT POLYTOPE: the uniformized
% evolution equation is written for E[X_p], E[X_p^2] and E[X_p1 X_p2], which
% gives linear equalities among the mean marking x, the enabling probabilities
% q and the products y(p,t) = E[X_p e_t]; behavioural and probabilistic
% inequalities are added on top; and every reported measure is then obtained by
% minimising and maximising its linear form over that polytope. Any stationary
% point of the true chain satisfies every row, so the two optima BRACKET the
% exact value whatever the polytope leaves out.
%
% This is the Petri-net sibling of the QRF bounds in SolverBA: same technique,
% a different index space, and a LINEAR objective, so there is no stationary
% point to escape from and the answer is a property of the model alone.
%
% VARIABLES, over place levels l = 1..L and modes e = 1..E:
%
%   x(l)     E[X_l], mean tokens                            >= 0
%   q(e)     P(mode e enabled)                              in [0,1]
%   th(e)    throughput of mode e                           >= 0
%   u(e)     state-equation firing counts                   >= 0
%   y(l,e)   E[X_l e_e]                                     >= 0   (Markovian only)
%
% u is EXISTENTIAL and is not reported: E[X] is a convex combination of
% reachable markings, each of which is m0 + C h for some nonnegative integer h,
% so the mean satisfies m0 + C u for some nonnegative real u. It is not a mean
% firing count and has no steady-state value.
%
% LEVELS ARE (place, class) PAIRS, PLACE-MAJOR, level (pp-1)*R + k, the same
% coordinates SPN_MDD, SPN_SINVARIANTS and SPN_CONV use. MODES are (transition,
% mode) pairs in node order. Both index spaces are returned in BND.
%
% THE TOKEN COUNTS CREATED BY A FIRING ARE DETERMINISTIC IN LINE, which removes
% a whole branch of the reference: it allows sigma_{t,p}(n) to be random and
% splits the covariance family into an independent case (its eq. 7) and a
% selective one (its eq. 8). SETFIRINGOUTCOME takes an integer weight, so
% E[sigma^2] = sigma^2 and E[sigma_p1 sigma_p2] = sigma_p1 sigma_p2 hold
% exactly and eq. (7) is the correct form. Eq. (8) has no LINE model behind it
% and is deliberately absent.
%
% -- Input
% SN      : a NetworkStruct holding Places and Transitions
% OPTIONS : struct, all fields optional
%    markovian  true (default) uses the second-moment, covariance and Little's
%               law families, which need exponential firing times; false drops
%               them and the whole y block, leaving the operational bound,
%               which needs only a mean firing time and so admits any
%               phase-type law. A law that is not phase-type at all is refused
%               by both: SN carries such a mode as its own parameter list
%               rather than a (D0,D1) pair, and no mean can be read off it
%               without a per-distribution table
%    assumelive false (default). true adds the two liveness rows, which are
%               valid only on a live net; see the note below
%    init       initial tokens per place level, place-major; omitted reads the
%               declared marking through State.initialOccupancy
%    tol        LP feasibility slack added to the inequality sides, default 0
%    verbose    print the polytope size
%
% -- Output
% BND : struct with fields
%    places      place node indices, in level order
%    levelname   1 x L cell, "Place.Class"
%    modes       1 x E struct array with fields trans, mode, enab, inhib,
%                fire, rate
%    tokens      2 x L, row 1 the minimum and row 2 the maximum of x(l)
%    placeTput   2 x L, bracket on sum_e pi_e(l) mu_e q(e), tokens removed
%                from level l per unit time
%    modeTput    2 x E, bracket on th(e)
%    modeUtil    2 x E, bracket on q(e)
%    bound       1 x L, the a priori per-level bound B_l read off the
%                P-invariants, Inf where none constrains the level
%    nplacelevels, nclasses, markovian, nvars, nrows
%
% LIVENESS IS OFF BY DEFAULT, AND THAT IS DELIBERATE. The reference's two
% liveness rows (sum_t q_t >= 1 and x_p <= sum_t y_{p,t}) hold only on a live
% net, and liveness is not something this function can cheaply certify -- an
% inhibitor arc alone is enough to deadlock a net that looks well formed. A
% bound that silently assumed it would be wrong rather than loose on exactly
% the models where a bound is most wanted, so the rows are available under
% OPTIONS.assumelive and absent otherwise.
%
% WHAT THE ROWS ARE WORTH, MEASURED. They are the whole of the lower side. On
% the reference's own Table 2 (its Fig. 2b production line, five rate vectors)
% ASSUMELIVE reproduces its published l.b. column to four decimals -- 1.1653
% against 1.165, 1.8288 against 1.829, 1.5814 against 1.581, 1.3592 against
% 1.359, 1.3497 against 1.350 -- while without them the Markovian lower bound
% collapses onto the OPERATIONAL one (0.9302, 1.4815, 1.1111, 1.1110 against
% that column's 0.930, 1.481, 1.111, 1.111) on four of the five. The upper side
% needs neither row and matches the published u.b.2 either way.
%
% -- Reference
% Z. Liu, "Performance Analysis of Stochastic Timed Petri Nets Using Linear
% Programming Approach", IEEE Trans. Software Engineering 24(11), 1998,
% 1014-1030. The constraint families are its Table 1, p. 1022; the bracket
% statement is its Theorem 3, p. 1021.
%
% See also SPN_SINVARIANTS, SPN_METRICS, SPN_MDD, SOLVER_BA_SPNLP_ANALYZER.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(options), options = struct(); end
markovian = i_opt(options, 'markovian', true);
assumelive = i_opt(options, 'assumelive', false);
tol = i_opt(options, 'tol', 0);
verbose = i_opt(options, 'verbose', false);

places = find(sn.nodetype == NodeType.Place);
if isempty(places) || ~any(sn.nodetype == NodeType.Transition)
    line_error(mfilename, 'the model holds no Place or no Transition node');
end
places = places(:)';
R = sn.nclasses;
P = numel(places);
L = P * R;

md = i_modes(sn, places, R, markovian);
E = numel(md);

% Arc tables as E x L matrices: PI consumes, SG creates, ETA inhibits (Inf
% where there is no inhibitor arc), NET = SG - PI is the incidence.
PI = zeros(E,L); SG = zeros(E,L); ETA = inf(E,L);
for e = 1:E
    PI(e,:) = md(e).enab; SG(e,:) = md(e).fire; ETA(e,:) = md(e).inhib;
end
NET = SG - PI;
mu = [md.rate];

% Per-level a priori bounds and the conserved sums, both off the same
% minimal-support P-invariant basis. The reference writes its "cycle
% population" family for UNWEIGHTED cycles; SPN_SINVARIANTS returns the
% weighted invariants S m = V, which are equally linear and strictly tighter,
% so those are what is emitted.
if isfield(options, 'init') && ~isempty(options.init)
    init = options.init(:)';
else
    init = i_initmarking(sn, places, R);
end
inv = spn_sinvariants(sn, init);
S = inv.S;
V = inv.V;
m0 = inv.m0;
B = i_levelbounds(S, V, L);

% ---- variable layout
ix = 1:L;
iq = L + (1:E);
ith = L + E + (1:E);
iu = L + 2*E + (1:E);
nv = L + 3*E;
if markovian
    iy = reshape(nv + (1:(L*E)), L, E);   % iy(l,e)
    nv = nv + L*E;
else
    iy = [];
end

lb = zeros(nv,1);
ub = inf(nv,1);
ub(iq) = 1;
for l = 1:L
    if isfinite(B(l))
        ub(ix(l)) = B(l);
        if markovian, ub(iy(l,:)) = B(l); end
    end
end

eqr = []; eqc = []; eqv = []; beq = []; neq = 0;
inr = []; inc = []; inw = []; bin = []; nin = 0;

% ---- (1) throughput: th_e = mu_e q_e
for e = 1:E
    addrow([ith(e), iq(e)], [1, -mu(e)], 0, 'eq');
end

% ---- (2) flow balance: tokens are created at a level at the rate they are
% consumed there. Holds for any stable net, Markovian or not.
for l = 1:L
    addrow(iq, mu .* NET(:,l)', 0, 'eq');
end

% ---- (3)+(4) second moment and population covariance, from the stationarity
% of E[X_l1 X_l2] under the uniformized chain. Table 1 writes the q side as
% four sums over set intersections; since the memberships are exactly "sigma >
% 0" and "pi > 0", those collapse to
%     -(sigma_1 - pi_1)(sigma_2 - pi_2) = -net_1 net_2
% per mode, which also makes the l1 == l2 case reduce to the second-moment
% family with no separate derivation, as the reference notes it must.
if markovian
    for l1 = 1:L
        for l2 = l1:L
            % duplicate columns at l1 == l2 are summed by SPARSE, which is
            % what supplies the factor of two the reference's (6) carries
            addrow([iy(l1,:), iy(l2,:), iq], ...
                [mu .* NET(:,l2)', mu .* NET(:,l1)', mu .* (NET(:,l1) .* NET(:,l2))'], ...
                0, 'eq');
        end
    end
end

% ---- (5) liveness, only when the caller vouches for it
if assumelive
    addrow(iq, ones(1,E), 1, 'ge');
    if markovian
        for l = 1:L
            addrow([ix(l), iy(l,:)], [1, -ones(1,E)], 0, 'le');
        end
    end
end

% ---- (6) conflicting transitions: a mode that consumes no more and is
% inhibited no sooner is enabled whenever the other is
for e1 = 1:E
    for e2 = 1:E
        if e1 == e2, continue; end
        if all(PI(e1,:) <= PI(e2,:)) && all(ETA(e1,:) >= ETA(e2,:))
            addrow([iq(e1), iq(e2)], [1, -1], 0, 'ge');
        end
    end
end

% ---- (7) boundedness, per level; and (8) cycle population, as the weighted
% invariant equalities and their y companions
if markovian
    for l = 1:L
        if ~isfinite(B(l)), continue; end
        for e = 1:E
            addrow([iy(l,e), iq(e)], [1, -B(l)], 0, 'le');
            addrow([ix(l), iy(l,e), iq(e)], [1, -1, B(l)], B(l), 'le');
            if B(l) > 0
                addrow([ix(l), iy(l,e), iq(e)], [1 - 1/B(l), -1, 1], 0, 'ge');
            end
        end
    end
end
for i = 1:size(S,1)
    addrow(ix, S(i,:), V(i), 'eq');
    if markovian
        for e = 1:E
            addrow([iy(:,e)', iq(e)], [S(i,:), -V(i)], 0, 'eq');
        end
    end
end

% ---- (9) reachable marking: the mean lies in the state-equation cone
for l = 1:L
    addrow([ix(l), iu], [1, -NET(:,l)'], m0(l), 'eq');
end

% ---- (10) sample-path comparisons
if markovian
    mutot = sum(mu);
    for l = 1:L
        for e = 1:E
            addrow([iy(l,e), ix(l)], [1, -1], 0, 'le');
            if PI(e,l) > 0
                addrow([iy(l,e), iq(e)], [1, -PI(e,l)], 0, 'ge');
            end
            if isfinite(ETA(e,l))
                addrow([iy(l,e), iq(e)], [1, -(ETA(e,l) - 1)], 0, 'le');
            end
        end
        addrow([ix(l), iy(l,:)], [mutot, -mu], 0, 'ge');
    end
    for e = 1:E
        ent = find(PI(e,:) > 0);
        if isscalar(ent) && ~any(isfinite(ETA(e,:)))
            addrow([ix(ent), iy(ent,e)], [1, -1], PI(e,ent) - 1, 'le');
        end
    end
end

% ---- (11) enabling bounds, from Chernoff's inequality on the marking
for e = 1:E
    ent = find(PI(e,:) > 0);
    inh = find(isfinite(ETA(e,:)));
    D = numel(ent) + numel(inh);
    if D == 0, continue; end
    if ~isempty(ent) && all(isfinite(B(ent)))
        idx = iq(e); val = 1; rhs = 1; ok = true;
        for l = ent
            den = B(l) - PI(e,l) + 1;
            if den <= 0, ok = false; break; end
            idx = [idx, ix(l)]; val = [val, -1/den]; rhs = rhs - B(l)/den; %#ok<AGROW>
        end
        if ok
            for l = inh
                idx = [idx, ix(l)]; val = [val, 1/ETA(e,l)]; %#ok<AGROW>
            end
            addrow(idx, val, rhs, 'ge');
        end
    end
    % Upper side. Every term of the sum over input levels carries a "min"
    % operator, and Table 1's convention is that either operand may be taken;
    % each choice is a valid row and the whole set is the tightest linear
    % relaxation, so all of them are emitted while the count stays small.
    ok = true;
    for l = inh
        if ~isfinite(B(l)) || B(l) - ETA(e,l) + 1 <= 0, ok = false; break; end
    end
    if ok
        base = 0; bidx = []; bval = [];
        for l = inh
            den = B(l) - ETA(e,l) + 1;
            base = base + B(l)/den;
            bidx = [bidx, ix(l)]; bval = [bval, 1/den]; %#ok<AGROW>
        end
        nc = numel(ent);
        if nc <= 4, combos = 0:(2^nc - 1); else, combos = [0, 2^nc - 1]; end
        for c = combos
            idx = [iq(e), bidx]; val = [D, bval]; rhs = base;
            for j = 1:nc
                l = ent(j);
                if bitget(c, j) == 0
                    idx = [idx, ix(l)]; val = [val, -1/PI(e,l)]; %#ok<AGROW>
                else
                    rhs = rhs + 1;
                end
            end
            addrow(idx, val, rhs, 'le');
        end
    end
end

% ---- (12) Little's law at each level: the mean sojourn time of a token is at
% least the mean minimum firing time of the modes that can remove it
if markovian
    for l = 1:L
        out = sum(mu(PI(:,l) > 0));
        if out <= 0, continue; end
        addrow([ix(l), iq], [out, -(mu .* SG(:,l)')], 0, 'ge');
    end
end

Aeq = sparse(eqr, eqc, eqv, neq, nv);
Ain = sparse(inr, inc, inw, nin, nv);
if tol > 0, bin = bin + tol; end

if verbose
    line_printf(['\nSPN -> LP: %d place levels, %d modes, %d variables, ' ...
        '%d equalities, %d inequalities\n'], L, E, nv, neq, nin);
end
lpopt = optimoptions('linprog', 'Display', 'off');

tokens = nan(2, L);
placeTput = nan(2, L);
for l = 1:L
    [tokens(1,l), tokens(2,l)] = bracket(ix(l), 1);
    if any(PI(:,l) > 0)
        [placeTput(1,l), placeTput(2,l)] = bracket(iq, mu .* PI(:,l)');
    else
        placeTput(:,l) = 0;
    end
end
modeTput = nan(2, E);
modeUtil = nan(2, E);
for e = 1:E
    [modeTput(1,e), modeTput(2,e)] = bracket(ith(e), 1);
    [modeUtil(1,e), modeUtil(2,e)] = bracket(iq(e), 1);
end

levelname = cell(1, L);
for pp = 1:P
    for k = 1:R
        levelname{(pp-1)*R + k} = sprintf('%s.%s', sn.nodenames{places(pp)}, sn.classnames{k});
    end
end

bnd = struct('places', places, 'levelname', {levelname}, 'modes', md, ...
    'tokens', tokens, 'placeTput', placeTput, 'modeTput', modeTput, ...
    'modeUtil', modeUtil, 'bound', B, 'nplacelevels', L, 'nclasses', R, ...
    'markovian', markovian, 'nvars', nv, 'nrows', neq + nin);

% =========================================================================
    function addrow(idx, val, rhs, kind)
        % One row a'z (=,<=,>=) rhs. Duplicate column indices are left in
        % place: SPARSE sums them, which is how the l1 == l2 covariance row
        % picks up its factor of two.
        idx = idx(:)'; val = val(:)';
        keep = val ~= 0;
        idx = idx(keep); val = val(keep);
        if isempty(idx), return; end
        switch kind
            case 'eq'
                neq = neq + 1;
                eqr = [eqr, neq*ones(1,numel(idx))];
                eqc = [eqc, idx];
                eqv = [eqv, val];
                beq = [beq; rhs];
            case 'le'
                nin = nin + 1;
                inr = [inr, nin*ones(1,numel(idx))];
                inc = [inc, idx];
                inw = [inw, val];
                bin = [bin; rhs];
            case 'ge'
                nin = nin + 1;
                inr = [inr, nin*ones(1,numel(idx))];
                inc = [inc, idx];
                inw = [inw, -val];
                bin = [bin; -rhs];
        end
    end

% =========================================================================
    function [lo, hi] = bracket(idx, val)
        % Minimum and maximum of one linear form over the polytope.
        f = zeros(nv,1);
        f(idx) = val;
        lo = i_solve(f, Ain, bin, Aeq, beq, lb, ub, lpopt);
        hi = -i_solve(-f, Ain, bin, Aeq, beq, lb, ub, lpopt);
    end
end

% -------------------------------------------------------------------------
function v = i_opt(options, name, default)
if isfield(options, name) && ~isempty(options.(name))
    v = options.(name);
else
    v = default;
end
end

% -------------------------------------------------------------------------
function md = i_modes(sn, places, R, markovian)
% The (transition, mode) table over place-major levels. SPN_MDD builds the same
% table, but only as a step of reachable-set construction, which is the cost
% this bound exists to avoid; the arc conversion below is deliberately a copy
% of its I_ARCVEC rather than a call into it.
transitions = find(sn.nodetype == NodeType.Transition);
md = repmat(struct('trans',0,'mode',0,'enab',[],'inhib',[],'fire',[],'rate',0), 0, 1);
for ind = transitions(:)'
    param = sn.nodeparam{ind};
    for m = 1:param.nmodes
        if isfield(param,'timing') && numel(param.timing) >= m && ...
                param.timing(m) == TimingStrategy.IMMEDIATE
            line_error(mfilename, sprintf(['mode %d of node %d is IMMEDIATE; the moment ' ...
                'relaxation is written for a net whose transitions all have finite rates, ' ...
                'so vanishing states must be eliminated first'], m, ind));
        end
        if isfield(param,'firingdep') && numel(param.firingdep) >= m && ~isempty(param.firingdep{m})
            line_error(mfilename, sprintf(['mode %d of node %d has a marking-dependent firing ' ...
                'rate; the uniformization step needs one rate per mode'], m, ind));
        end
        if param.nmodeservers(m) ~= 1
            line_error(mfilename, sprintf(['mode %d of node %d has %g servers; the relaxation ' ...
                'is derived under single-server semantics, where the firing rate is mu*q. Its ' ...
                'infinite-server form needs the K-fold transition expansion of the reference''s ' ...
                'Section 7, which is not implemented'], m, ind, param.nmodeservers(m)));
        end
        % A PHASE-TYPE FIRING LAW IS WHERE THE TWO VARIANTS PART. The mean of a
        % (D0,D1) pair is pie*(-D0)^-1*1 and is all the operational bound needs;
        % the Markovian one needs the marking alone to be the state, which a
        % multi-phase mode breaks. A law that is not phase-type at all reaches
        % SN as its own parameter list (Pareto stores {shape, scale}), with no
        % mean recoverable without a per-distribution table, so BOTH variants
        % refuse it rather than guess.
        proc = param.firingproc{m};
        if ~iscell(proc) || numel(proc) < 2 || ~isnumeric(proc{1}) || ~ismatrix(proc{1}) || ...
                size(proc{1},1) ~= size(proc{1},2)
            line_error(mfilename, sprintf(['mode %d of node %d has a firing law that is not ' ...
                'phase-type; SN carries its own parameters rather than a (D0,D1) pair, so ' ...
                'neither the Markovian nor the operational bound can read a mean firing rate ' ...
                'from it. Use a phase-type law, or an exact solver'], m, ind));
        end
        D0 = full(proc{1}); D1 = full(proc{2});
        nph = size(D0,1);
        if markovian && nph > 1
            line_error(mfilename, sprintf(['mode %d of node %d has a phase-type firing time; ' ...
                'the relaxation is written over the marking alone, and a phase-type mode needs ' ...
                'the state-machine expansion of the reference''s Section 7, which is not ' ...
                'implemented. Use the operational bound, which needs only the mean'], m, ind));
        end
        if nph == 1
            rate = D1(1);
        else
            pv = param.firingpie{m};
            if isempty(pv), pv = [1, zeros(1, nph-1)]; end
            pv = pv(:)' / sum(pv);
            rate = 1 / (pv * ((-D0) \ ones(nph,1)));
        end
        if ~(rate > 0) || ~isfinite(rate)
            line_error(mfilename, sprintf(['mode %d of node %d has mean firing rate %g; a ' ...
                'bound needs a finite positive one'], m, ind, rate));
        end
        e = numel(md) + 1;
        md(e,1).trans = ind;
        md(e,1).mode = m;
        md(e,1).enab  = i_arcvec(param.enabling{m},   places, R, sn.nnodes, 0);
        md(e,1).inhib = i_arcvec(param.inhibiting{m}, places, R, sn.nnodes, Inf);
        md(e,1).fire  = i_arcvec(param.firing{m},     places, R, sn.nnodes, 0);
        md(e,1).rate = rate;
    end
end
if isempty(md)
    line_error(mfilename, 'the net has no firing mode');
end
md = md(:)';
end

% -------------------------------------------------------------------------
function v = i_arcvec(mat, places, R, nnodes, fillval)
% (nnodes x nclasses) arc matrix -> 1 x (P*R) place-major level vector. The row
% index is a NODE index; FILLVAL is 0 for enabling and firing and Inf for
% inhibiting, where Inf means "no arc" and 0 would mean "inhibited always".
m = reshape(mat, nnodes, R);
P = numel(places);
v = fillval * ones(1, P * R);
for pp = 1:P
    for k = 1:R
        x = m(places(pp), k);
        if isfinite(fillval) && fillval == 0, x = max(0, x); end
        v((pp - 1) * R + k) = x;
    end
end
end

% -------------------------------------------------------------------------
function init = i_initmarking(sn, places, R)
% Declared tokens per (place, class), read through the canonical decoder.
P = numel(places);
init = zeros(1, P * R);
anyset = false;
for pp = 1:P
    for k = 1:R
        n = State.initialOccupancy(sn, places(pp), k);
        init((pp - 1) * R + k) = n;
        anyset = anyset || n > 0;
    end
end
if ~anyset
    init = [];   % let SPN_SINVARIANTS take the reference-station default
end
end

% -------------------------------------------------------------------------
function B = i_levelbounds(S, V, L)
% Tightest a priori bound per level: an invariant with weight w on level l and
% conserved value V caps that level at floor(V/w).
B = inf(1, L);
for i = 1:size(S,1)
    for l = 1:L
        if S(i,l) > 0
            B(l) = min(B(l), floor(V(i) / S(i,l)));
        end
    end
end
end

% -------------------------------------------------------------------------
function val = i_solve(f, Ain, bin, Aeq, beq, lb, ub, lpopt)
% One LP. EXITFLAG is the wrong success predicate here (see
% _kb/11-conventions-and-gotchas.md); the finiteness of the answer is.
[x, fval] = linprog(f, Ain, bin, Aeq, beq, lb, ub, lpopt);
if isempty(x) || ~all(isfinite(x)) || ~isfinite(fval)
    val = NaN;
else
    val = fval;
end
end
