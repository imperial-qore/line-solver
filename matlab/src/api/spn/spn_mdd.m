function [mdds, desc, info] = spn_mdd(model, options)
% [MDDS, DESC, INFO] = SPN_MDD(MODEL)
% [MDDS, DESC, INFO] = SPN_MDD(MODEL, OPTIONS)
% Decision-diagram reachable set and Kronecker rate descriptor of a stochastic
% Petri net, so that MDD_MCD can analyse it.
%
% Levels are of two kinds:
%   place levels : one per (Place, class) pair, holding a token count. A
%                  multiclass net therefore has P*R of them, ordered
%                  place-major, so level (p-1)*R+k is class k in place p.
%   phase levels : one per mode whose firing time has more than one phase,
%                  holding the phase the running server occupies.
%
% The rate structure factorises exactly under single-server firing semantics:
% a mode fires at a constant rate whenever every input level holds its
% enabling multiplicity and no inhibitor level has reached its threshold, so
%
%   W_l^e[i, i + fire(l) - enab(l)] = 1     for enab(l) <= i < inhib(l)
%
% at every place level. A phase-type mode contributes two event families on
% its phase level, the internal phase changes D0 (marking unchanged) and the
% firings D1 (marking moved), each gated by the same per-level enabling
% indicators. Both are products of per-level terms, which is what Eq. 1 of the
% paper requires.
%
% -- Input
% MODEL   : a Network holding Places and Transitions
% OPTIONS : struct, fields
%             bound    - per-place-level token bound (scalar or 1 x P*R);
%                        inferred from a place invariant when omitted
%             phmemory - 'exact' (default) or 'resume', see the note below
%             descriptor - build the Kronecker rate descriptor (default true).
%                        Pass false for the MDD-rec route, which reads only the
%                        reachable set: see the note below
%             verbose  - print the net summary (default false)
% -- Output
% MDDS    : MDD.toStruct of the reachable set
% DESC    : Kronecker descriptor for MDD_MCD, with fields K, domain, events,
%           levelkind (1 place, 2 phase) and, when one exists,
%           invariant = struct('weights', w, 'value', v) with zero weight on
%           the phase levels
% INFO    : struct with fields places, placenames, classnames, levelname,
%           levelkind, modes, init, mdd, phaseof
%
% -- Phase-type firing and the memory policy
% LINE discards a running server's phase when its mode becomes disabled
% (State/afterGlobalEvent zeroes the server phases on the disable action), i.e.
% preemptive repeat. Resetting a mode's phase is then triggered by a JOINT
% condition on the place levels, which is not a product of per-level terms and
% has no Kronecker form. What this descriptor encodes is preemptive resume: a
% disabled mode's phase freezes and continues when the mode is re-enabled. The
% two policies coincide exactly when a mode is never disabled while running,
% so reachability records, for free, whether any phase-type mode was ever found
% disabled. OPTIONS.phmemory='exact' (the default) errors when one was;
% 'resume' proceeds deliberately with the resume semantics.
%
% -- Reachable-set-only mode (OPTIONS.descriptor = false)
% MDD_REC and SPN_METRICS need the reachable set and the metadata, never the
% rate descriptor: the product form supplies the rates. So the restrictions
% below that exist only because a Kronecker form must factorise per level --
% marking-dependent firing rates and multi-server modes drawing from several
% places -- are lifted, and DESC comes back carrying K, domain, levelkind and
% the invariant but no events. In exchange the firing times must be
% EXPONENTIAL: a phase level is a descriptor device, and a product-form
% marking process is memoryless in the marking alone. Immediate modes stay
% refused, since a vanishing marking carries no probability and would have to
% be eliminated from the reachable set first.
%
% -- Other restrictions (each is an error, never a silent approximation)
% No immediate transitions (they make vanishing states, which must be
% eliminated before a Kronecker rate descriptor exists) and no
% marking-dependent firing rates. A multi-server mode is accepted only when its
% enabling touches ONE level, because the enabling degree
% min_l floor(m(l)/enab(l)) is otherwise not a product of per-level terms.
%
% See also: mdd_mcd, mdd_descriptor, mdd_reachset.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(options), options = struct(); end
if ~isfield(options, 'verbose'),  options.verbose = false; end
if ~isfield(options, 'bound'),    options.bound = []; end
if ~isfield(options, 'phmemory'), options.phmemory = 'exact'; end
if ~isfield(options, 'descriptor'), options.descriptor = true; end

sn = model.getStruct();
R = sn.nclasses;
places = find(sn.nodetype == NodeType.Place);
transitions = find(sn.nodetype == NodeType.Transition);
if isempty(places) || isempty(transitions)
    line_error(mfilename, 'the model holds no Place or no Transition node');
end
P = numel(places);
L = P * R;                                  % place levels, place-major

% ---- collect the (transition, mode) pairs
md = struct('trans', {}, 'mode', {}, 'enab', {}, 'inhib', {}, 'fire', {}, ...
    'D0', {}, 'D1', {}, 'pie', {}, 'nph', {}, 'srv', {}, 'dep', {});
for tt = 1:numel(transitions)
    ind = transitions(tt);
    param = sn.nodeparam{ind};
    for m = 1:param.nmodes
        if isfield(param, 'timing') && numel(param.timing) >= m ...
                && param.timing(m) == TimingStrategy.IMMEDIATE
            line_error(mfilename, sprintf(['mode %d of node %d is IMMEDIATE; vanishing states ' ...
                'must be eliminated before the net has a Kronecker rate descriptor'], m, ind));
        end
        dep = [];
        if isfield(param, 'firingdep') && numel(param.firingdep) >= m
            dep = param.firingdep{m};
        end
        if ~isempty(dep) && options.descriptor
            line_error(mfilename, sprintf(['mode %d of node %d has a marking-dependent firing ' ...
                'rate; g(marking) is not a product of per-level terms'], m, ind));
        end
        proc = param.firingproc{m};
        if ~iscell(proc) || numel(proc) < 2
            line_error(mfilename, sprintf(['mode %d of node %d has no Markovian firing process; ' ...
                'a general distribution has no finite phase level'], m, ind));
        end
        D0 = full(proc{1}); D1 = full(proc{2});
        nph = size(D0, 1);
        pv = param.firingpie{m};
        if isempty(pv), pv = [1, zeros(1, nph - 1)]; end
        e = numel(md) + 1;
        md(e).trans = ind; md(e).mode = m;
        md(e).enab  = i_arcvec(param.enabling{m},   places, R, sn.nnodes, 0);
        md(e).inhib = i_arcvec(param.inhibiting{m}, places, R, sn.nnodes, Inf);
        md(e).fire  = i_arcvec(param.firing{m},     places, R, sn.nnodes, 0);
        md(e).D0 = D0; md(e).D1 = D1; md(e).pie = pv(:)' / sum(pv);
        md(e).nph = nph; md(e).srv = param.nmodeservers(m); md(e).dep = dep;
        if ~options.descriptor && nph > 1
            line_error(mfilename, sprintf(['mode %d of node %d has a phase-type firing time; ' ...
                'the reachable-set-only mode carries no phase level, and a product-form ' ...
                'marking process must be memoryless in the marking alone'], m, ind));
        end
    end
end
E = numel(md);
for e = 1:E
    if ~options.descriptor, break; end
    if md(e).srv ~= 1 && nnz(md(e).enab) > 1
        line_error(mfilename, sprintf(['mode %d of node %d has %g servers and draws from %d ' ...
            'levels; the enabling degree min_l floor(m(l)/enab(l)) is then not a product of ' ...
            'per-level terms and admits no Kronecker form'], ...
            md(e).mode, md(e).trans, md(e).srv, nnz(md(e).enab)));
    end
    if md(e).srv ~= 1 && nnz(md(e).enab) == 0
        line_error(mfilename, sprintf(['mode %d of node %d has %g servers but consumes from no ' ...
            'place, so its enabling degree is unbounded and its firing rate undefined'], ...
            md(e).mode, md(e).trans, md(e).srv));
    end
end

% ---- phase levels for the multi-phase modes
phaseof = zeros(1, E);                      % 0 = no phase level
if options.descriptor
    for e = 1:E
        if md(e).nph > 1
            phaseof(e) = L + nnz(phaseof) + 1;
        end
    end
end
Q = nnz(phaseof);
K = L + Q;

netm = zeros(E, L);
for e = 1:E, netm(e, :) = md(e).fire - md(e).enab; end

% ---- initial marking and per-level bounds
init0 = i_initmarking(model, sn, places, R);
[winv, vinv] = i_placeinvariant(netm, init0);
if ~isempty(options.bound)
    bound = options.bound(:)';
    if isscalar(bound), bound = bound * ones(1, L); end
elseif ~isempty(winv)
    bound = zeros(1, L);
    for l = 1:L
        if winv(l) > 0, bound(l) = floor(vinv / winv(l)); else, bound(l) = vinv; end
    end
else
    line_error(mfilename, ['the net has no place invariant with positive weights, so the ' ...
        'marking is not bounded a priori; pass OPTIONS.bound']);
end
domain = [bound + 1, zeros(1, Q)];
for e = 1:E
    if phaseof(e) > 0, domain(phaseof(e)) = md(e).nph; end
end

init = zeros(1, K);
init(1:L) = init0;
for e = 1:E
    if phaseof(e) > 0
        init(phaseof(e)) = find(md(e).pie > 0, 1) - 1;
    end
end

% ---- reachable set; the closure records which modes were ever disabled, so
% the phase-memory question is answered without a second pass over |S|
everDisabled = false(1, E);
nextfun = @(s) i_next(s, md, netm, phaseof, domain, L);
mdd = MDD(domain);
mdd.insert(init);
frontier = init; head = 1;
while head <= size(frontier, 1)
    s = frontier(head, :); head = head + 1;
    [T, dis] = nextfun(s);
    everDisabled = everDisabled | dis;
    for r = 1:size(T, 1)
        t = T(r, :);
        if ~mdd.member(t)
            mdd.insert(t);
            frontier(end + 1, :) = t; %#ok<AGROW>
        end
    end
    if head > 1024 && 2 * head > size(frontier, 1)
        frontier = frontier(head:end, :); head = 1;
    end
end
mdd.compact();
mdds = mdd.toStruct();

badph = find(phaseof > 0 & everDisabled);
if ~isempty(badph) && options.descriptor && ~strcmpi(options.phmemory, 'resume')
    line_error(mfilename, sprintf(['mode %d of node %d has a phase-type firing time AND is ' ...
        'disabled in some reachable marking. LINE discards the phase on disabling (preemptive ' ...
        'repeat) but that reset is a joint condition on the place levels and has no Kronecker ' ...
        'form, so this descriptor would encode preemptive resume instead and disagree with ' ...
        'SolverCTMC. Pass OPTIONS.phmemory=''resume'' to accept the resume semantics.'], ...
        md(badph(1)).mode, md(badph(1)).trans));
end

% ---- Kronecker event matrices
events = repmat(struct('a', 0, 'b', 0, 'lev', [], 'W', []), 0, 1);
for e = 1:E * double(options.descriptor)
    gate = find(md(e).enab > 0 | isfinite(md(e).inhib));      % levels gating enabling
    move = find(netm(e, :) ~= 0);                             % levels whose count moves
    touched = union(gate, move);
    degl = 0;
    if md(e).srv ~= 1, degl = find(md(e).enab > 0, 1); end     % level carrying the degree
    if phaseof(e) == 0
        % single-phase mode: one firing event, rate on the first touched level
        if isempty(touched), continue; end
        lev = []; W = {};
        for l = touched
            lev(end+1) = l;                                   %#ok<AGROW>
            W{end+1} = i_placemat(l, md(e), netm(e, l), domain(l), isequal(l, degl)); %#ok<AGROW>
        end
        W{1} = md(e).D1 * W{1};                               % scalar rate
        events(end+1, 1) = struct('a', md(e).trans, 'b', md(e).mode, 'lev', lev, 'W', {W}); %#ok<AGROW>
    else
        q = phaseof(e);
        % (1) internal phase changes: marking unchanged, gated by enabling
        D0off = md(e).D0 - diag(diag(md(e).D0));
        if nnz(D0off) > 0
            lev = []; W = {};
            for l = gate
                lev(end+1) = l;                               %#ok<AGROW>
                W{end+1} = i_placemat(l, md(e), 0, domain(l)); %#ok<AGROW>
            end
            lev(end+1) = q; W{end+1} = sparse(D0off);
            events(end+1, 1) = struct('a', md(e).trans, 'b', md(e).mode, 'lev', lev, 'W', {W}); %#ok<AGROW>
        end
        % (2) firings: marking moved, phase redrawn through D1
        lev = []; W = {};
        for l = touched
            lev(end+1) = l;                                   %#ok<AGROW>
            W{end+1} = i_placemat(l, md(e), netm(e, l), domain(l), isequal(l, degl)); %#ok<AGROW>
        end
        lev(end+1) = q; W{end+1} = sparse(md(e).D1);
        events(end+1, 1) = struct('a', md(e).trans, 'b', md(e).mode, 'lev', lev, 'W', {W}); %#ok<AGROW>
    end
end

desc.K = K;
desc.domain = domain;
desc.events = events;
desc.levelkind = [ones(1, L), 2 * ones(1, Q)];
if ~isempty(winv)
    desc.invariant = struct('weights', [winv, zeros(1, Q)], 'value', vinv);
end

% ---- descriptive information
info.places = places;
info.placenames = cell(1, P);
for pp = 1:P, info.placenames{pp} = sn.nodenames{places(pp)}; end
info.classnames = cell(1, R);
for k = 1:R, info.classnames{k} = sn.classnames{k}; end
info.levelkind = desc.levelkind;
info.levelname = cell(1, K);
for pp = 1:P
    for k = 1:R
        info.levelname{(pp - 1) * R + k} = sprintf('%s.%s', info.placenames{pp}, info.classnames{k});
    end
end
for e = 1:E
    if phaseof(e) > 0
        info.levelname{phaseof(e)} = sprintf('phase(%s.m%d)', sn.nodenames{md(e).trans}, md(e).mode);
    end
end
info.modes = md;
info.nnodes = sn.nnodes;
info.nclasses = R;
info.descriptor = options.descriptor;
info.init = init;
info.mdd = mdd;
info.phaseof = phaseof;
info.everDisabled = everDisabled;
info.nplacelevels = L;

if options.verbose
    line_printf('\nSPN -> MDD: %d places x %d classes = %d place levels, %d phase levels\n', P, R, L, Q);
    line_printf('  modes = %d, |S| = %d, bounds = %s\n', E, mdd.cardinality(), mat2str(bound));
end
end

% ------------------------------------------------------------------------
function v = i_arcvec(mat, places, R, nnodes, fillval)
% (nnodes x nclasses) arc matrix -> 1 x (P*R) place-major level vector
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

% ------------------------------------------------------------------------
function W = i_placemat(l, mde, net, d, applydegree)
% Local matrix of one mode at place level l: move the count by net, but only
% from local states that satisfy this level's enabling and inhibition.
% APPLYDEGREE scales each row by the enabling degree min(floor(m/enab), srv),
% i.e. the number of concurrently firing servers. It is carried by the single
% enabling level of a multi-server mode; with one server the degree is 1 and
% the indicator below is already the whole story.
i = 0:(d - 1);
ok = (i >= mde.enab(l)) & (i < mde.inhib(l));
j = i + net;
ok = ok & (j >= 0) & (j <= d - 1);
val = ones(1, d);
if nargin >= 5 && applydegree && mde.enab(l) > 0
    val = min(floor(i / mde.enab(l)), mde.srv);
end
W = sparse(i(ok) + 1, j(ok) + 1, val(ok), d, d);
end

% ------------------------------------------------------------------------
function init = i_initmarking(model, sn, places, R)
% Tokens per (place, class) at time zero, from the model state when set and
% from the reference station otherwise.
P = numel(places);
init = zeros(1, P * R);
anyset = false;
for pp = 1:P
    node = model.getNodeByIndex(places(pp));
    st = node.getState();
    if ~isempty(st)
        anyset = true;
        for k = 1:min(R, numel(st))
            init((pp - 1) * R + k) = st(k);
        end
    end
end
if ~anyset || all(init == 0)
    for k = 1:R
        ref = sn.refstat(k);
        for pp = 1:P
            if sn.nodeToStation(places(pp)) == ref
                init((pp - 1) * R + k) = sn.njobs(k);
            end
        end
    end
end
end

% ------------------------------------------------------------------------
function [w, val] = i_placeinvariant(netm, init)
% A place invariant w'*m = const is a NON-NEGATIVE right null vector of the
% incidence matrix over the place levels, and what bounds the marking is a
% STRICTLY POSITIVE one: the net is structurally bounded exactly when its places
% are covered by such invariants.
%
% Scanning a null-space BASIS for a positive vector is not enough, and the
% fork-join net is the counterexample: its two minimal-support invariants are
% (1,1,0,1) and (1,0,1,1), neither positive, while their sum (2,1,1,2) is. Which
% basis a codebase's null() returns then decides whether the net is accepted,
% which is how MATLAB and the JAR came to disagree on it. So the non-negative
% generators are computed directly, by Farkas' algorithm on [netm' | I] -- the
% same construction SPN_SINVARIANTS uses on the incidence matrix -- and summed.
w = []; val = [];
L = size(netm, 2);
if all(abs(sum(netm, 2)) < 1e-12)
    w = ones(1, L); val = sum(init); return
end
E = size(netm, 1);
M = netm';                                  % L x E, row l is level l's column
B = eye(L);                                 % the combination each row stands for
for e = 1:E
    keep = find(abs(M(:, e)) < 1e-12);
    pos = find(M(:, e) > 1e-12);
    neg = find(M(:, e) < -1e-12);
    Mn = M(keep, :); Bn = B(keep, :);
    for a = pos(:)'
        for b = neg(:)'
            c = -M(b, e) * M(a, :) + M(a, e) * M(b, :);
            d = -M(b, e) * B(a, :) + M(a, e) * B(b, :);
            g = i_gcdvec([c, d]);
            if g > 0, c = c / g; d = d / g; end
            Mn(end + 1, :) = c; Bn(end + 1, :) = d; %#ok<AGROW>
        end
    end
    [M, B] = i_minimalsupport(Mn, Bn);
end
if isempty(B), return; end
s = sum(B, 1);
if all(s > 1e-9)
    w = s / min(s(s > 1e-9));
    val = w * init(:);
end
end

% ------------------------------------------------------------------------
function g = i_gcdvec(v)
% gcd of the integral entries, 0 when any entry is not an integer.
g = 0;
for x = v(:)'
    if abs(x - round(x)) > 1e-9, g = 0; return; end
    g = gcd(g, abs(round(x)));
end
end

% ------------------------------------------------------------------------
function [M, B] = i_minimalsupport(M, B)
% Drop every row whose support strictly contains another's, which is what leaves
% the minimal supports and stops the pair expansion from blowing up.
n = size(B, 1);
supp = abs(B) > 1e-12;
drop = false(1, n);
for i = 1:n
    if drop(i), continue; end
    for j = 1:n
        if i == j || drop(j), continue; end
        if all(supp(j, :) <= supp(i, :)) && any(supp(i, :) > supp(j, :))
            drop(i) = true; break
        end
    end
end
M = M(~drop, :); B = B(~drop, :);
end

% ------------------------------------------------------------------------
function [T, disabled] = i_next(s, md, netm, phaseof, domain, L)
% Successor states of s, plus a flag per mode saying it was disabled here.
E = numel(md);
m = s(1:L);
T = zeros(0, numel(s));
disabled = false(1, E);
for e = 1:E
    if any(m < md(e).enab) || any(m >= md(e).inhib)
        disabled(e) = true;
        continue
    end
    if phaseof(e) == 0
        t = s; t(1:L) = m + netm(e, :);
        if all(t(1:L) >= 0) && all(t(1:L) <= domain(1:L) - 1)
            T(end + 1, :) = t; %#ok<AGROW>
        end
    else
        q = phaseof(e); ph = s(q) + 1;
        D0off = md(e).D0(ph, :); D0off(ph) = 0;
        for j = find(D0off ~= 0)
            t = s; t(q) = j - 1;
            T(end + 1, :) = t; %#ok<AGROW>
        end
        for j = find(md(e).D1(ph, :) ~= 0)
            t = s; t(1:L) = m + netm(e, :); t(q) = j - 1;
            if all(t(1:L) >= 0) && all(t(1:L) <= domain(1:L) - 1)
                T(end + 1, :) = t; %#ok<AGROW>
            end
        end
    end
end
end
