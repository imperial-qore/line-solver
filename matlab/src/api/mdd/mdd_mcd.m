function out = mdd_mcd(mdds, desc, options)
% OUT = MDD_MCD(MDDS, DESC)
% OUT = MDD_MCD(MDDS, DESC, OPTIONS)
% Miner-Ciardo-Donatelli approximate stationary analysis: solve a structured
% CTMC whose EXACT reachable state space is stored in a decision diagram, by
% building and iterating K level-CTMCs (a decision-diagram-guided aggregation),
% after A.S. Miner, G. Ciardo, S. Donatelli, "Using the exact state space of a
% Markov model to compute approximate stationary measures", SIGMETRICS 2000.
%
% The method never forms the |S|-state generator or probability vector. It
% keeps one CTMC per decision-diagram level k, over states M_k = {(p,i_k)} with
% p a level-k node and i_k a local state on a non-null arc, and iterates the
% coupled system to a fixed point. The single approximation (Eq. 5) is
% Pr{i_k | alpha} = Pr{i_k | p}: the local-state law at level k depends only on
% the node p, not the full path above it -- justified by the exact reachability
% the diagram encodes. For product-form models the method is EXACT (paper
% Sec. 5), so on a single-class closed QN it reproduces SolverCTMC.
%
% Orientation note: the paper indexes levels K (top/root) down to 1 (bottom/
% terminal); the MDD class uses level 1 as the root. This routine works in the
% paper's orientation, mapping paper level k to MDD level (K+1-k), i.e. to
% station (K+1-k).
%
% -- Input
% MDDS    : struct from MDD.toStruct (the reachable set, MDD orientation)
% DESC    : Kronecker rate descriptor from mdd_descriptor
% OPTIONS : struct, fields tol (1e-12), maxiter (500), verbose (false),
%           initpik (optional 1 x K cell of level warm-start vectors)
% -- Output
% OUT     : struct with fields
%             QLen   - 1 x K mean jobs per station
%             X      - 1 x K per-station throughput
%             U      - 1 x K utilisation
%             pik    - 1 x K cell, level-k stationary vectors over M_k
%             Mrows  - 1 x K cell, Mrows{k}(r,:) = [node local] of M_k row r
%             levelSizes - 1 x K, |M_k|
%             iters  - fixed-point iterations performed
%
% See also: mdd_descriptor, mdd_reachset, mdd_closedqn.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(options), options = struct(); end
if ~isfield(options, 'tol'),     options.tol = 1e-12;  end
if ~isfield(options, 'maxiter'), options.maxiter = 500; end
if ~isfield(options, 'verbose'), options.verbose = false; end

TERM_TRUE = -1;
K = mdds.K;

% ---- paper orientation: paper level k <-> MDD level (K+1-k) = station (K+1-k)
Pnode = cell(1, K); nn = zeros(1, K); dom = zeros(1, K); stationOf = zeros(1, K);
for k = 1:K
    oL = K + 1 - k;
    Pnode{k} = mdds.node{oL};
    nn(k) = mdds.nnodes(oL);
    dom(k) = mdds.domain(oL);
    stationOf(k) = oL;
end

% ---- per (event, paper level) local matrices W_k^e and enabling rates lambda
E = numel(desc.events);
Wall = cell(E, 1); lam = cell(E, 1);
for e = 1:E
    ev = desc.events(e);
    plev = K + 1 - ev.lev;                 % station levels -> paper levels
    Wall{e} = cell(1, K); lam{e} = cell(1, K);
    for k = 1:K
        idx = find(plev == k, 1);
        if isempty(idx)
            Wall{e}{k} = speye(dom(k));      % untouched level: identity
        else
            Wall{e}{k} = ev.W{idx};
        end
        lam{e}{k} = full(sum(Wall{e}{k}, 2));  % row sums = local enabling rate
    end
end

% ---- level-k CTMC state sets M_k = {(p, v) : arc p[v] non-null}
Mrows = cell(1, K); Midx = cell(1, K);
for k = 1:K
    tbl = Pnode{k};
    if k == 1
        [pp, cc] = find(tbl == TERM_TRUE);   % bottom level: TRUE arcs
    else
        [pp, cc] = find(tbl > 0);            % real child ids
    end
    pp = pp(:); cc = cc(:);                   % find() returns rows for 1-row tbl
    Mrows{k} = [pp, cc - 1];                 % [node, local(0-based)]
    Midx{k} = sparse(pp, cc, 1:numel(pp), nn(k), dom(k));
end
levelSizes = cellfun(@(m) size(m, 1), Mrows);

% ---- initialise level stationary vectors and node marginals
pik = cell(1, K); Prp = cell(1, K);
[above, below] = i_pathcounts(mdds, K);
haveInit = isfield(options, 'initpik') && ~isempty(options.initpik);
if ~haveInit
    pik = i_uniforminit(mdds, Mrows, K, above, below);
end
for k = 1:K
    if haveInit
        pik{k} = options.initpik{k}(:);
    end
    Prp{k} = accumarray(Mrows{k}(:, 1), pik{k}, [nn(k), 1]);
end

% ---- fixed-point iteration (Fig. 3, procedure Solve)
iters = 0; converged = false; delta = Inf;
for it = 1:options.maxiter
    iters = it;
    piold = pik;

    % ComputeBs, bottom-up: b_k^e[p] = sum_v Pr{v|p} * b_{k-1}^e[p[v]] * lambda
    bcell = cell(1, K);
    for k = 1:K
        bk = zeros(nn(k), E);
        rows = Mrows{k};
        for r = 1:size(rows, 1)
            p = rows(r, 1); v = rows(r, 2);
            if Prp{k}(p) <= 0, continue; end
            adjust = pik{k}(r) / Prp{k}(p);           % Pr{v|p}
            for e = 1:E
                le = lam{e}{k}(v + 1);
                if le == 0, continue; end             % e not locally enabled
                if k > 1
                    down = bcell{k - 1}(Pnode{k}(p, v + 1), e);
                else
                    down = 1;                          % terminal ONE
                end
                bk(p, e) = bk(p, e) + adjust * down * le;
            end
        end
        bcell{k} = bk;
    end

    % top-down: ComputeAs(k) then SolveLevel(k)
    Acell = cell(1, K);
    for e = 1:E, AcellK{e} = eye(nn(K)); end %#ok<AGROW>
    Acell{K} = AcellK;
    for k = K:-1:1
        if k < K
            Acell{k} = i_computeAs(k, Acell{k + 1}, Pnode, pik, Wall, Mrows, nn, E);
        end
        % SolveLevel(k): assemble R_k (Eq. 6), solve pi_k Q_k = 0, refresh Pr{p}
        Rk = i_computeMC(k, Acell{k}, bcell, Pnode, Wall, Mrows, Midx, levelSizes, E);
        Qk = Rk - diag(sum(Rk, 2));
        pik{k} = i_solvestat(Qk);
        Prp{k} = accumarray(Mrows{k}(:, 1), pik{k}, [nn(k), 1]);
    end

    % NaN-aware: MATLAB max() OMITS NaN, so max(delta,NaN) returns delta and a
    % diverged iterate would be reported as converged at the next sweep.
    delta = 0;
    for k = 1:K
        dk = max(abs(pik{k} - piold{k}));
        if ~isfinite(dk)
            line_error(mfilename, sprintf(['level %d iterate is not finite at iteration %d; ' ...
                'the level-%d CTMC did not yield a proper stationary vector.'], k, it, k));
        end
        delta = max(delta, dk);
    end
    if delta < options.tol, converged = true; break; end
end
if ~converged
    line_error(mfilename, sprintf(['the coupled level iteration did not converge in %d sweeps ' ...
        '(last change %.3e against tol %.3e); the level marginals returned would not be a ' ...
        'fixed point. Raise OPTIONS.maxiter or relax OPTIONS.tol.'], ...
        options.maxiter, delta, options.tol));
end

% ---- performance measures from the per-level marginals. QLen is the mean
% local value and is defined for any descriptor (jobs at a station, tokens in a
% place); X and U need the queueing parameters and are skipped without them.
isQN = isfield(desc, 'mu') && ~isempty(desc.mu) && isfield(desc, 'servers');
hasMap = isfield(desc, 'valuemap') && ~isempty(desc.valuemap);
QLen = zeros(1, K); X = []; U = [];
if isQN
    mu = desc.mu; servers = desc.servers;
    X = zeros(1, K); U = zeros(1, K);
end
for s = 1:K
    k = K + 1 - s;                       % paper level of station/place s
    % A level whose local state encodes more than a count (a station holding
    % both a population and a service phase) carries a map from local index to
    % the physical quantity; without one the index IS the quantity.
    if hasMap
        vmap = desc.valuemap{s}(:);
        v = vmap(Mrows{k}(:, 2) + 1);
    else
        v = Mrows{k}(:, 2);
    end
    pk = pik{k};
    QLen(s) = v' * pk;                    % E[occupancy of level s]
    if isQN
        busy = min(v, servers(s));
        X(s) = mu(s) * (busy' * pk);
        if isinf(servers(s))
            U(s) = v' * pk;
        else
            U(s) = (busy' * pk) / servers(s);
        end
    end
end

% The level chains are coupled only through rates, so nothing in the iteration
% forces the marginals to describe the same population; a fixed point that does
% not is a wrong answer, not an approximation, and must not be returned. The
% test is a conservation law of the model: the closed population for a QN, a
% place invariant w'*m = const for a Petri net.
[winv, vinv] = i_invariant(desc, K);
if ~isempty(winv)
    got = winv * QLen(:);
    if abs(got - vinv) > 1e-6 * max(1, abs(vinv))
        line_error(mfilename, sprintf(['the level marginals converged to an invariant value of ' ...
            '%.6g against the model value %g, so the fixed point reached is degenerate (the ' ...
            'level chains are mutually inconsistent). Supply OPTIONS.initpik with a consistent ' ...
            'starting law.'], got, vinv));
    end
end

out.QLen = QLen; out.X = X; out.U = U;
out.pik = pik; out.Mrows = Mrows; out.levelSizes = levelSizes; out.iters = iters;
[out.pathsPerLevel, out.noAggregation] = i_exactness(above, K);

if options.verbose
    line_printf('\nMDD-MCD approximate aggregation: %d levels, %d fixed-point iters\n', K, iters);
    line_printf('  level-CTMC sizes |M_k| = [%s] (max %d)\n', ...
        strtrim(sprintf('%d ', levelSizes)), max(levelSizes));
    line_printf('  mean queue lengths     = %s\n', mat2str(QLen, 5));
end
end

% ------------------------------------------------------------------------
function [w, val] = i_invariant(desc, K)
% Conservation law the converged marginals must satisfy, as w'*QLen = val.
w = []; val = [];
if isfield(desc, 'invariant') && ~isempty(desc.invariant)
    w = desc.invariant.weights(:)';
    val = desc.invariant.value;
elseif isfield(desc, 'N') && ~isempty(desc.N)
    w = ones(1, K);                       % closed QN: total population
    val = desc.N;
end
end

% ------------------------------------------------------------------------
function [above, below] = i_pathcounts(mdds, K)
% ABOVE{oL}(p) is the number of distinct root-to-p paths, |A(p)| in the paper's
% notation, and BELOW{oL}(p) the number of accepted states under p. Both are
% O(number of nodes) and serve the uniform initialisation and the exactness certificate.
TERM_TRUE = -1;
below = cell(1, K); above = cell(1, K);
for oL = K:-1:1
    nb = zeros(mdds.nnodes(oL), 1);
    tbl = mdds.node{oL};
    for v = 1:mdds.domain(oL)
        ch = tbl(:, v);
        if oL == K
            nb = nb + double(ch == TERM_TRUE);
        else
            nz = ch > 0;
            nb(nz) = nb(nz) + below{oL + 1}(ch(nz));
        end
    end
    below{oL} = nb;
end
for oL = 1:K
    above{oL} = zeros(mdds.nnodes(oL), 1);
end
above{1}(mdds.root) = 1;
for oL = 1:(K - 1)
    tbl = mdds.node{oL};
    for v = 1:mdds.domain(oL)
        ch = tbl(:, v);
        nz = ch > 0;
        if ~any(nz), continue; end
        above{oL + 1} = above{oL + 1} + ...
            accumarray(ch(nz), above{oL}(nz), [mdds.nnodes(oL + 1), 1]);
    end
end
end

% ------------------------------------------------------------------------
function [perLevel, noAggregation] = i_exactness(above, K)
% Structural certificate that the aggregation loses nothing.
%
% The single approximation is Pr{i_k | alpha} = Pr{i_k | p}: the local-state law
% is conditioned on the NODE rather than on the whole path above it. When a node
% is reached by exactly one path, |A(p)| = 1, conditioning on the node IS
% conditioning on the path and the identity is exact. If that holds at every
% node the fixed point is the exact stationary law, and no reference solve is
% needed to know it.
%
% SUFFICIENT, not necessary: a product-form model is exact too (paper Sec. 5)
% however much its diagram shares, and that is a property of the model rather
% than of the diagram. A false here means "not certified", never "approximate".
% Note also that max |A(p)| = 1 means no node is shared, i.e. the diagram
% compresses nothing, so exactness by this route and a useful saving are
% mutually exclusive.
perLevel = ones(1, K);
for k = 1:K
    oL = K + 1 - k;
    if ~isempty(above{oL}), perLevel(k) = max(above{oL}); end
end
noAggregation = all(perLevel <= 1 + 1e-12);
end

% ------------------------------------------------------------------------
function pik = i_uniforminit(mdds, Mrows, K, above, below)
% Uniform law over the EXACT reachable set, projected onto each level:
% Pr{(p,v)} = (paths root->p) * (states below arc p[v]) / |S|. A flat law over
% M_k instead treats level states as equiprobable irrespective of how many
% global states they stand for, which breaks the population invariant the
% diagram encodes; from about K=8 the coupled iteration then descends into the
% basin of the DEGENERATE empty-population fixed point (all mass on local state
% 0 at every level, a true fixed point since no station can then emit) and
% converges to it with zero residual. The projection below is consistent across
% levels by construction, so the iteration starts inside the physical simplex.
pik = cell(1, K);
for k = 1:K
    oL = K + 1 - k;                       % paper level k is MDD level K+1-k
    rows = Mrows{k};
    p = rows(:, 1); v = rows(:, 2) + 1;
    if oL == K
        w = above{oL}(p);                 % a TRUE arc stands for one state
    else
        ch = mdds.node{oL}(sub2ind(size(mdds.node{oL}), p, v));
        w = above{oL}(p) .* below{oL + 1}(ch);
    end
    pik{k} = w / sum(w);
end
end

% ------------------------------------------------------------------------
function Ak = i_computeAs(k, Aup, Pnode, pik, Wall, Mrows, nn, E)
% ComputeAs(k): A_k^e from A_{k+1}^e (the "from above" contribution), Fig. 3.
% The adjust denominator Pr{p[v]} is the FROM-ABOVE marginal of the child node,
% Pr{p} = sum_{alpha in A(p)} Pr{alpha} = sum over parents of pi_{k+1}. Using it
% (rather than the level-k CTMC marginal, which only equals it at convergence)
% makes adjust a proper conditional Pr{(parent,arc)|child} and pins the
% inter-level node marginals, removing the spurious fixed points.
rows1 = Mrows{k + 1};
PrAbove = zeros(nn(k), 1);
for r = 1:size(rows1, 1)
    child = Pnode{k + 1}(rows1(r, 1), rows1(r, 2) + 1);
    if child > 0, PrAbove(child) = PrAbove(child) + pik{k + 1}(r); end
end
Ak = cell(1, E);
for e = 1:E, Ak{e} = zeros(nn(k), nn(k)); end
for r = 1:size(rows1, 1)
    p = rows1(r, 1); v = rows1(r, 2);
    childp = Pnode{k + 1}(p, v + 1);            % p[v]: node at level k
    if childp <= 0 || PrAbove(childp) <= 0, continue; end
    adjust = pik{k + 1}(r) / PrAbove(childp);   % pi_{k+1}[(p,v)] / Pr{p[v]}
    for e = 1:E
        wrow = Wall{e}{k + 1}(v + 1, :);
        [~, wj, wval] = find(wrow);
        if isempty(wj), continue; end
        arow = Aup{e}(p, :);
        qcols = find(arow ~= 0);
        if isempty(qcols), continue; end
        for wi = 1:numel(wj)
            w = wj(wi) - 1; wv = wval(wi);
            for q = qcols
                childq = Pnode{k + 1}(q, w + 1);
                if childq <= 0, continue; end   % q[w] null
                Ak{e}(childp, childq) = Ak{e}(childp, childq) + arow(q) * wv * adjust;
            end
        end
    end
end
end

% ------------------------------------------------------------------------
function Rk = i_computeMC(k, Ak, bcell, Pnode, Wall, Mrows, Midx, levelSizes, E)
% ComputeMC(k): level-k rate matrix, R_k^e[(p,i),(q,j)] = A_k^e[p,q]*W_k^e[i,j]*b_{k-1}^e[p[i]] (Eq. 6)
nm = levelSizes(k);
rows = Mrows{k};
Ri = zeros(0, 1); Ci = zeros(0, 1); Vv = zeros(0, 1);
for r = 1:nm
    p = rows(r, 1); v = rows(r, 2);
    for e = 1:E
        wrow = Wall{e}{k}(v + 1, :);
        [~, wj, wval] = find(wrow);
        if isempty(wj), continue; end
        if k > 1
            bfac = bcell{k - 1}(Pnode{k}(p, v + 1), e);
        else
            bfac = 1;                            % terminal ONE
        end
        if bfac == 0, continue; end
        arow = Ak{e}(p, :);
        qcols = find(arow ~= 0);
        if isempty(qcols), continue; end
        for wi = 1:numel(wj)
            w = wj(wi) - 1; wv = wval(wi);
            for q = qcols
                di = Midx{k}(q, w + 1);
                if di == 0, continue; end        % destination (q,w) not in M_k
                Ri(end + 1, 1) = r;              %#ok<AGROW>
                Ci(end + 1, 1) = di;            %#ok<AGROW>
                Vv(end + 1, 1) = arow(q) * wv * bfac; %#ok<AGROW>
            end
        end
    end
end
Rk = sparse(Ri, Ci, Vv, nm, nm);
end

% ------------------------------------------------------------------------
function p = i_solvestat(Q)
% Stationary distribution of a small irreducible generator: p*Q = 0, sum p = 1.
% The normalisation is APPENDED rather than substituted for the last balance
% equation: overwriting a row discards a constraint and left Q' singular to
% working precision from about |M_k| = 325 upwards, so the solve returned NaN.
% The overdetermined system has full column rank whenever the level chain is
% irreducible, and QR least squares solves it stably at the same cost order.
n = size(Q, 1);
if n == 1, p = 1; return; end
A = [full(Q)'; ones(1, n)];
rhs = zeros(n + 1, 1); rhs(n + 1) = 1;
p = A \ rhs;
p(p < 0) = 0;
s = sum(p);
if ~isfinite(s) || s <= 0
    line_error(mfilename, sprintf(['level CTMC of order %d admits no proper stationary ' ...
        'distribution (the level generator is reducible or numerically degenerate).'], n));
end
p = p / s;
end
