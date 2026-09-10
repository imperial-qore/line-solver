function OUT = spn_bench_mdd(OPTIONS)
% OUT = SPN_BENCH_MDD()
% OUT = SPN_BENCH_MDD(OPTIONS)
% Cost profile of the decision-diagram aggregation on a SCALABLE stochastic
% Petri net, the flexible-manufacturing-system (FMS) family of
% G. Ciardo, K. Trivedi, "A decomposition approach for stochastic reward net
% models", Perform. Eval. 18(1), 1993, which is the standard scaling benchmark
% for decision-diagram state-space methods.
%
% The net has three machine lines, one job class each, that take a part from a
% common pool, process it and return it, plus an assembly stage that consumes
% one finished part of every class at once. The reachable set grows steeply in
% the part count N while the number of levels stays fixed, which is the regime
% the aggregation is built for: |S| grows polynomially of degree 3R in N and
% sum_k |M_k| does not.
%
% -- Input
% OPTIONS : struct, fields
%             parts    - vector of part counts to sweep (default 1:5)
%             maxstates- skip a row whose exact solve would exceed this many
%                        states (default 200000); the aggregation still runs
%             verbose  - print the table as it is produced (default true)
% -- Output
% OUT     : struct array, one entry per part count, with fields parts, nstates,
%           levels, sumMk, compression, iters, t_reach, t_solve, noAggregation
%           and, when the exact solve was run, t_exact and maxerr
%
% -- Remarks
% The exact reference is the explicit generator assembled through MDD/index,
% which is why MAXSTATES gates it: the aggregation is the cheap side and the
% point of the benchmark is to watch the two diverge.
%
% See also: spn_mdd, mdd_mcd, mdd_bench_mcd.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 1 || isempty(OPTIONS), OPTIONS = struct(); end
if ~isfield(OPTIONS, 'parts'),     OPTIONS.parts = 1:5; end
if ~isfield(OPTIONS, 'maxstates'), OPTIONS.maxstates = 200000; end
if ~isfield(OPTIONS, 'verbose'),   OPTIONS.verbose = true; end

if OPTIONS.verbose
    line_printf('\n%-6s %10s %7s %8s %7s %6s %9s %9s %9s %8s\n', ...
        'parts', '|S|', 'levels', 'sum|Mk|', 'compr', 'iters', ...
        't_reach', 't_solve', 't_exact', 'maxerr');
end

OUT = repmat(struct('parts', 0, 'nstates', 0, 'levels', 0, 'sumMk', 0, ...
    'compression', 0, 'iters', 0, 't_reach', 0, 't_solve', 0, ...
    'noAggregation', false, 't_exact', NaN, 'maxerr', NaN), 0, 1);

for np = OPTIONS.parts(:)'
    model = i_fms(np);
    t0 = tic;
    [mdds, desc, info] = spn_mdd(model, struct());
    t_reach = toc(t0);
    t0 = tic;
    out = mdd_mcd(mdds, desc, struct());
    t_solve = toc(t0);

    ns = info.mdd.cardinality();
    row = struct('parts', np, 'nstates', ns, 'levels', desc.K, ...
        'sumMk', sum(out.levelSizes), 'compression', ns / sum(out.levelSizes), ...
        'iters', out.iters, 't_reach', t_reach, 't_solve', t_solve, ...
        'noAggregation', out.noAggregation, 't_exact', NaN, 'maxerr', NaN);

    % The exact reference is the explicit generator over the SAME reachable set,
    % so any difference is the aggregation and not a different state space.
    if ns <= OPTIONS.maxstates
        t0 = tic;
        [Qex, ok] = i_exactmarking(info.mdd, desc);
        row.t_exact = toc(t0);
        if ok
            isplace = desc.levelkind == 1;
            row.maxerr = max(abs(out.QLen(isplace) - Qex(isplace)));
        end
    end
    OUT(end + 1, 1) = row; %#ok<AGROW>

    if OPTIONS.verbose
        line_printf('%-6d %10d %7d %8d %7.1fx %6d %9.3f %9.3f %9.3f %8.1e\n', ...
            row.parts, row.nstates, row.levels, row.sumMk, row.compression, ...
            row.iters, row.t_reach, row.t_solve, row.t_exact, row.maxerr);
    end
end
end

% ------------------------------------------------------------------------
function [QLen, ok] = i_exactmarking(mdd, desc)
% Mean tokens per level from the EXPLICIT generator over the reachable set,
% assembled through MDD/index so no state list is ever materialised.
QLen = []; ok = false;
n = mdd.cardinality();
S = mdd.enumerate();
K = desc.K;
ri = []; ci = []; vv = [];
for si = 1:n
    s = S(si, :);
    for e = 1:numel(desc.events)
        ev = desc.events(e);
        [tgt, rate] = i_fire(s, ev, K);
        for r = 1:size(tgt, 1)
            if rate(r) == 0, continue; end
            j = mdd.index(tgt(r, :));
            if j < 0 || j == si - 1, continue; end
            ri(end+1) = si; ci(end+1) = j + 1; vv(end+1) = rate(r); %#ok<AGROW>
        end
    end
end
Rk = sparse(ri, ci, vv, n, n);
Q = Rk - diag(sum(Rk, 2));
A = [full(Q)'; ones(1, n)];
rhs = zeros(n + 1, 1); rhs(n + 1) = 1;
p = A \ rhs;
p(p < 0) = 0;
if ~isfinite(sum(p)) || sum(p) <= 0, return; end
p = p / sum(p);
QLen = zeros(1, K);
for k = 1:K
    QLen(k) = sum(p(:)' .* S(:, k)');
end
ok = true;
end

% ------------------------------------------------------------------------
function [tgt, rate] = i_fire(s, ev, K)
% Every target of one Kronecker event from state s, with its rate. The event
% touches only ev.lev, so the product runs over those levels alone.
tgt = s; rate = 1;
for t = 1:numel(ev.lev)
    l = ev.lev(t);
    W = ev.W{t};
    row = W(s(l) + 1, :);
    [~, jj, vv] = find(row);
    if isempty(jj), tgt = zeros(0, K); rate = zeros(0, 1); return; end
    newtgt = zeros(0, K); newrate = zeros(0, 1);
    for a = 1:size(tgt, 1)
        for b = 1:numel(jj)
            cand = tgt(a, :); cand(l) = jj(b) - 1;
            newtgt(end+1, :) = cand;      %#ok<AGROW>
            newrate(end+1, 1) = rate(a) * vv(b); %#ok<AGROW>
        end
    end
    tgt = newtgt; rate = newrate;
end
end

% ------------------------------------------------------------------------
function model = i_fms(np)
% Flexible manufacturing system with three machine lines, NP parts each.
%
% One CLASS per line, so class c circulates Pool.c -> Busy.c -> Done.c and the
% assembly stage consumes one finished part of every class at once and returns
% one of each to the pool. A single-class version of the same net is NOT
% ergodic: with one token pool the parts can all pile into one line, leaving the
% assembly permanently starved and the chain absorbing, which mdd_mcd rightly
% refuses to solve. The class per line is what makes the join always eventually
% enabled, and it is also the only case exercising the multiclass (place,class)
% levels of spn_mdd.
model = Network(sprintf('FMS%d', np));
pool = Place(model, 'Pool');
busy = Place(model, 'Busy');
done = Place(model, 'Done');

R = 3;
jc = cell(1, R);
for c = 1:R
    jc{c} = ClosedClass(model, sprintf('Part%d', c), np, pool);
end

startT = cell(1, R); endT = cell(1, R);
rates = [1.0 1.3 1.7];
for c = 1:R
    startT{c} = Transition(model, sprintf('Start%d', c));
    m = startT{c}.addMode(sprintf('S%d', c));
    startT{c}.setDistribution(m, Exp(2.0));
    startT{c}.setNumberOfServers(m, 1);
    startT{c}.setEnablingConditions(m, jc{c}, pool, 1);
    startT{c}.setFiringOutcome(m, jc{c}, busy, 1);

    endT{c} = Transition(model, sprintf('End%d', c));
    m = endT{c}.addMode(sprintf('E%d', c));
    endT{c}.setDistribution(m, Exp(rates(c)));
    endT{c}.setNumberOfServers(m, 1);
    endT{c}.setEnablingConditions(m, jc{c}, busy, 1);
    endT{c}.setFiringOutcome(m, jc{c}, done, 1);
end

% the join: one finished part of every class in, one of every class back out,
% so the net stays conservative per class and keeps a place invariant per class
asm = Transition(model, 'Assemble');
m = asm.addMode('Asm');
asm.setDistribution(m, Exp(3.0));
asm.setNumberOfServers(m, 1);
for c = 1:R
    asm.setEnablingConditions(m, jc{c}, done, 1);
    asm.setFiringOutcome(m, jc{c}, pool, 1);
end

% Every (node, class) pair needs a routing strategy even when that class never
% visits the node: refreshStruct reads the strategy of each class at each node
% and rejects an unset one. The classes that never reach a transition are given
% the same successor as the class that does.
P = model.initRoutingMatrix();
for c = 1:R
    P.set(jc{c}, jc{c}, pool, startT{c}, 1.0);
    P.set(jc{c}, jc{c}, busy, endT{c}, 1.0);
    P.set(jc{c}, jc{c}, done, asm, 1.0);
    P.set(jc{c}, jc{c}, asm, pool, 1.0);
    for c2 = 1:R
        P.set(jc{c2}, jc{c2}, startT{c}, busy, 1.0);
        P.set(jc{c2}, jc{c2}, endT{c}, done, 1.0);
    end
end
model.link(P);

pool.setState(np * ones(1, R));
busy.setState(zeros(1, R));
done.setState(zeros(1, R));
end
