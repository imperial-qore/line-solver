%{
 % @brief Enumerate the conservation laws (the "physics") that any solution of a
 %        layered queueing network must satisfy, as symbolic relations over a
 %        LayeredNetworkStruct and, optionally, as numeric residuals.
 %
 % @details
 % A layered model is not free to report any tuple of throughputs, think times
 % and utilizations: five families of relations tie them together, and every one
 % of them is determined by the STRUCTURE of the model alone. This function
 % walks a LayeredNetworkStruct and emits them, one record per relation, with
 % the index sets to aggregate over, the constant coefficients and a printable
 % form. Nothing is solved here.
 %
 % The families, with `kind` as emitted:
 %
 % - `little`   Little's law on a task's THREAD POOL. The threads of task t form
 %              a closed cycle of one delay stage (the surrogate think time that
 %              SolverLN imputes to the task, plus the declared think time of a
 %              reference task) and one service stage (holding a request from
 %              above). Writing B(t,k) for the mean number of threads of t busy
 %              serving caller class k -- the per-class utilization expressed in
 %              JOB units -- the law reads
 %
 %                  X(t)*(Z(t) + z(t)) + sum_k B(t,k) = N(t)
 %
 %              which is the update SolverLN.updateThinkTimes iterates on. In
 %              the [0,1]-normalized utilization LINE reports for a queueing
 %              station, B(t,k) = N(t)*U(t,k) and the law is the familiar
 %              X(t)*(Z(t)+z(t)) = N(t)*(1 - sum_k U(t,k)); at an infinite
 %              server LINE's utilization is already a job count, so B = U.
 %              The caller classes k are the CALLS that target an entry of t --
 %              the in-edges of t in the call graph -- plus each entry of t that
 %              carries an OPEN ARRIVAL, a stream that holds a thread exactly as
 %              a call does and that a task can have alongside its callers. Both
 %              are structural neighbours of t, which makes the whole relation
 %              node-local. `termisentry` says which of the two a term is.
 %
 % - `callflow` Throughput conservation across one call: X(c) = X(src(c))*y(c),
 %              with y(c) the mean number of calls and src(c) the dispatching
 %              activity (the dispatching ENTRY for a forwarding call).
 %
 % - `entryflow` The requests an entry serves are the calls that reach it plus
 %              its open-arrival stream: X(e) = sum_c X(c) + lambda(e).
 %
 % - `actflow`  An activity executes v(a) times per invocation of its entry:
 %              X(a) = X(e)*v(a). The visit counts v come from the activity
 %              precedence graph and are returned in `out.visits`. An AND-JOIN is
 %              the one place where flow does not add up -- its target executes
 %              once per fork, not once per branch -- so the arcs into a join are
 %              scaled by 1/(number of joined branches); see joinScaledGraph.
 %
 % - `hostutil` The utilization law at a processor: the host demand executed
 %              there per unit time equals its busy servers,
 %              sum_a X(a)*D(a) = m(h)*U(h) (again m(h)*U(h) is a job count, and
 %              the factor m(h) drops at an infinite server).
 %
 % Together these close the system: `little` alone is one equation per task and
 % admits the all-zero solution, so a physics-informed loss built on it should
 % carry the flow and utilization families as well.
 %
 % @par Reference-task, forwarding and open-arrival tasks:
 % A reference task has no layer above, so its cycle closes on its OWN entries,
 % which enter K(t) as entry classes. A forwarding target and an entry-arrival
 % target are driven by a rate pinned outside the task, but the shape of the
 % relation is identical. All four cases therefore share one template, and the
 % `branch` field says which one produced the record; what it selects is only the
 % label and whether the utilization is a job count.
 %
 % @par Conventions (also returned in out.convention):
 % - Rates and populations in a `little` record are PER REPLICA, matching
 %   SolverLN.updateThinkTimes: X per replica is tput/repl and N is the declared
 %   multiplicity of one copy. Utilizations and throughputs elsewhere are as the
 %   solver reports them, i.e. totalled over replicas -- which is why the server
 %   count in `hostutil` is `mult` and not `mult*repl`.
 % - N(t) is `lqn.mult`. SolverLN instead iterates on `njobs`, which carries the
 %   interlocking corrections and, under replication, may be `maxmult`; both are
 %   returned per record so a caller can substitute.
 % - A task whose multiplicity is infinite has no finite thread pool, so its
 %   `little` record is emitted with const = Inf and `degenerate` set. The
 %   sustainable multiplicity `maxmult` is the finite surrogate SolverLN uses.
 % - S(k) is the entry SERVICE time servt (phase 1 plus phase 2), the time a
 %   thread is held, not the residence time residt the caller waits for. The
 %   difference is the phase-2 tail, and `corr.phase2` marks the tasks that have
 %   one. `corr.setup` marks a task that pays a setup charge, which lands on Z.
 %
 % @par Saturation:
 % SolverLN clamps the think time at zero, so the `little` equality is an
 % INEQUALITY at a saturated task: when sum_k U(t,k) -> 1 the right-hand side
 % reaches zero and Z can no longer absorb the imbalance. Records carry
 % `clamped` once instantiated, true when the equality is not attainable. A loss
 % built on these relations should use a one-sided (hinge) form there.
 %
 % @par Syntax:
 % @code
 % out = lqn_balance_equations(lqn)            % symbolic only
 % out = lqn_balance_equations(lqn, sol)       % also numeric residuals
 % out = lqn_balance_equations(lqn, sol, UN)   % including the host utilization law
 % lqn_balance_equations(lqn)                  % print the report
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lqn<td>LayeredNetworkStruct (from LayeredNetwork.getStruct())
 % <tr><td>sol<td>optional; a solved SolverLN, or any struct exposing the same
 %                fields (tput, util, thinkt, servt, residt), used to
 %                instantiate every relation and report its residual
 % <tr><td>UN<td>optional; the REPORTED utilization per element, the second
 %               return of getEnsembleAvg. Needed by the `hostutil` family, whose
 %               right-hand side is the reported processor utilization and not the
 %               `util` iterate (which a host leaves at zero). It is a separate
 %               argument because getEnsembleAvg re-enters iterate(), and a
 %               diagnostic must not re-run a fixed point on its own.
 % </table>
 %
 % @par Returns (struct out):
 % <table>
 % <tr><th>Field<th>Description
 % <tr><td>eqs<td>(1 x M) struct array of relations; see below
 % <tr><td>text<td>cellstr, the printable report
 % <tr><td>visits<td>(nidx x 1) activity visit count per invocation of its entry
 % <tr><td>A_little<td>sparse (ntasks x ncalls) 1 where call c is a caller class of task t
 % <tr><td>A_flow<td>sparse (ncalls x nidx) call-to-source incidence, weighted by y(c)
 % <tr><td>A_host<td>sparse (nhosts x nidx) host-to-activity incidence, weighted by D(a)
 % <tr><td>maxresidual<td>largest absolute residual over the non-degenerate records, NaN without sol
 % </table>
 %
 % Each element of `out.eqs` carries: `kind`, `branch`, `target` (absolute index
 % the relation is anchored on), `targetname`, `terms` (absolute indices to
 % aggregate over), `termisentry` (per term, true for an ENTRY class and false
 % for a CALL class, since the two read their rate from different vectors),
 % `coeff` (their constant coefficients), `const` (the
 % right-hand side constant), `mult`, `maxmult`, `repl`, `scaled` (true when the
 % per-class utilization must be multiplied by `mult` to reach job units),
 % `corr`, `degenerate`, `latex`, `text`, and -- once instantiated -- `lhs`,
 % `rhs`, `residual`, `relresidual`, `clamped`.
%}
function out = lqn_balance_equations(lqn, sol, un)
if nargin < 2
    sol = [];
end
if nargin < 3
    un = [];
end
hasSol = ~isempty(sol);

nidx = lqn.nidx;
visits = actVisits(lqn);

% Calls that target each task, i.e. the caller classes of its thread pool. The
% dispatching element of a forwarding call is an entry, not an activity.
callsInto = cell(lqn.tshift + lqn.ntasks, 1);
for cidx = 1:lqn.ncalls
    dst = lqn.callpair(cidx, 2);
    if dst < 1 || dst > nidx
        continue
    end
    tidx = lqn.parent(dst);
    if tidx >= 1 && tidx <= numel(callsInto)
        callsInto{tidx} = [callsInto{tidx}, cidx];
    end
end

if hasSol
    s = readSolution(lqn, sol, un, visits);
else
    s = [];
end

eqs = emptyRecord();
eqs(1) = [];

% ---------------- Little's law on each thread pool -----------------------
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    [branch, K, isEntryTerm] = threadPoolBranch(lqn, tidx, callsInto{tidx});
    if isempty(branch)
        continue
    end
    eqs(end+1) = littleRecord(lqn, tidx, branch, K, isEntryTerm, s); %#ok<AGROW>
end

% ---------------- flow conservation --------------------------------------
for cidx = 1:lqn.ncalls
    eqs(end+1) = callflowRecord(lqn, cidx, s); %#ok<AGROW>
end
for e = 1:lqn.nentries
    eidx = lqn.eshift + e;
    inc = incomingCalls(lqn, eidx);
    lam = arrivalRate(lqn, eidx);
    if isempty(inc) && lam == 0
        continue    % a reference entry is driven by its own cycle, not by flow
    end
    eqs(end+1) = entryflowRecord(lqn, eidx, inc, lam, s); %#ok<AGROW>
end
for a = 1:lqn.nacts
    aidx = lqn.ashift + a;
    eidx = entryOfActivity(lqn, aidx);
    if eidx == 0
        continue
    end
    eqs(end+1) = actflowRecord(lqn, aidx, eidx, visits(aidx), s); %#ok<AGROW>
end

% ---------------- utilization law at each processor ----------------------
for h = 1:lqn.nhosts
    hidx = lqn.hshift + h;
    eqs(end+1) = hostutilRecord(lqn, hidx, s); %#ok<AGROW>
end

out = struct();
out.eqs = eqs;
out.visits = visits;
out.convention = conventionText();
[out.A_little, out.A_flow, out.A_host] = incidences(lqn, eqs);
out.maxresidual = NaN;
if hasSol
    r = [eqs.residual];
    d = [eqs.degenerate];
    r = r(~d & ~isnan(r));
    if ~isempty(r)
        out.maxresidual = max(abs(r));
    end
end
out.text = report(lqn, out, hasSol);

if nargout == 0
    fprintf('%s\n', out.text{:});
    clear out
end
end

% ========================================================================
% record construction
% ========================================================================
function r = emptyRecord()
r = struct('kind', '', 'branch', '', 'target', 0, 'targetname', '', ...
    'terms', [], 'termisentry', logical([]), ...
    'coeff', [], 'const', 0, 'mult', NaN, 'maxmult', NaN, ...
    'repl', 1, 'scaled', false, 'corr', struct('phase2', false, 'setup', false), ...
    'degenerate', false, 'latex', '', 'text', '', ...
    'lhs', NaN, 'rhs', NaN, 'residual', NaN, 'relresidual', NaN, 'clamped', false, ...
    'perclassutil', []);
end

% Which thread-pool case task TIDX falls in, and the caller classes of its pool.
%
% The term set is uniform across the cases: every request stream that can hold a
% thread of the task contributes one class. That is each CALL targeting one of
% its entries, plus each entry carrying an OPEN ARRIVAL -- a stream that holds
% threads exactly as a call does, and that a task can have alongside its callers
% -- plus, for a reference task, its own entries, since the cycle of a reference
% task closes on itself and no layer above drives it.
%
% The branch label follows the case analysis of SolverLN.updateThinkTimes, which
% is what decides whether the utilization is a job count (infinite server) or is
% normalized to [0,1] (every other discipline).
function [branch, K, isEntryTerm] = threadPoolBranch(lqn, tidx, cin)
branch = '';
ents = lqn.entriesof{tidx};
arv = ents(arrayfun(@(e) arrivalRate(lqn, e) > 0, ents));
K = [cin(:)', arv(:)'];
isEntryTerm = [false(1, numel(cin)), true(1, numel(arv))];
if lqn.isref(tidx)
    branch = 'ref';
    rest = setdiff(ents(:)', arv(:)');
    K = [K, rest];
    isEntryTerm = [isEntryTerm, true(1, numel(rest))];
    return
end
blocking = cin(arrayfun(@(c) full(lqn.calltype(c)) ~= CallType.FWD, cin));
if ~isempty(blocking)
    if lqn.sched(tidx) == SchedStrategy.INF
        branch = 'inf';
    else
        branch = 'queueing';
    end
    return      % a forwarded request holds a thread too, so K already has it
end
if ~isempty(cin)
    branch = 'fwd';
    return
end
if ~isempty(arv)
    branch = 'arrival';
    return
end
% no caller, no arrival: the task is unreachable and carries no cycle
K = [];
isEntryTerm = logical([]);
end

function r = littleRecord(lqn, tidx, branch, K, isEntryTerm, s)
r = emptyRecord();
r.kind = 'little';
r.branch = branch;
r.target = tidx;
r.targetname = elemName(lqn, tidx);
r.terms = K;
r.termisentry = isEntryTerm;
r.coeff = ones(1, numel(K));
r.mult = lqn.mult(tidx);
r.maxmult = multField(lqn, 'maxmult', tidx);
r.repl = max(1, lqn.repl(tidx));
r.const = r.mult;
r.scaled = ~strcmp(branch, 'inf');
r.corr.phase2 = hasPhase2(lqn, tidx);
r.corr.setup = logical(full(lqn.hassetup(tidx)));
r.degenerate = ~isfinite(r.const) || isempty(K);

z = lqn_ref_thinktime(lqn, tidx);
Ksym = cell(1, numel(K));
for i = 1:numel(K)
    Ksym{i} = termName(lqn, isEntryTerm(i), K(i));
end
r.latex = littleLatex(r, z, Ksym);
r.text = r.latex;

if ~isempty(s)
    [X, B, U] = littleValues(lqn, tidx, K, isEntryTerm, s);
    Zt = s.thinkt(tidx);
    if isnan(Zt)
        Zt = 0;
    end
    r.lhs = X * (Zt + z) + sum(B);
    r.rhs = r.const;
    r.residual = r.lhs - r.rhs;
    r.relresidual = r.residual / max(abs(r.rhs), GlobalConstants.FineTol);
    r.clamped = isfinite(r.const) && (r.const - sum(B) - X * z) < 0;
    r.perclassutil = U;
end
end

function txt = littleLatex(r, z, Ksym)
X = sprintf('X(%s)', r.targetname);
Z = sprintf('Z(%s)', r.targetname);
zs = '';
if z > 0
    zs = sprintf(' + %g', z);
end
lhs = sprintf('%s*(%s%s)', X, Z, zs);
for i = 1:numel(Ksym)
    lhs = sprintf('%s + B(%s)', lhs, Ksym{i});
end
txt = sprintf('%s = %s', lhs, num2strc(r.const));
if r.scaled && ~isempty(Ksym)
    us = '';
    for i = 1:numel(Ksym)
        us = sprintf('%s - U(%s,%s)', us, r.targetname, Ksym{i});
    end
    txt = sprintf('%s\n           equivalently  %s*%s = %s*(1%s)', ...
        txt, X, Z, num2strc(r.const), us);
end
end

% Per-caller-class busy threads and normalized utilization of task TIDX.
%
% A class holds a thread for the SERVICE time of the entry it targets, at the
% rate it drives that entry: the call rate for a call class, the entry rate for
% an entry class (an open-arrival stream, or the self-driven cycle of a
% reference task).
function [X, B, U] = littleValues(lqn, tidx, K, isEntryTerm, s)
X = s.tput(tidx) / max(1, lqn.repl(tidx));
if isnan(X)
    X = 0;
end
B = zeros(1, numel(K));
for i = 1:numel(K)
    if isEntryTerm(i)
        B(i) = s.tput(K(i)) * s.servt(K(i));
    else
        B(i) = s.calltput(K(i)) * s.servt(lqn.callpair(K(i), 2));
    end
    if isnan(B(i))
        B(i) = 0;
    end
end
if lqn.sched(tidx) == SchedStrategy.INF || ~isfinite(lqn.mult(tidx)) || lqn.mult(tidx) <= 0
    U = B;
else
    U = B / lqn.mult(tidx);
end
end

function r = callflowRecord(lqn, cidx, s)
r = emptyRecord();
r.kind = 'callflow';
r.branch = callTypeName(lqn.calltype(cidx));
r.target = lqn.cshift + cidx;
r.targetname = callName(lqn, cidx);
src = lqn.callpair(cidx, 1);
y = lqn.callproc_mean(cidx);
if isnan(y)
    y = 0;
end
r.terms = src;
r.coeff = y;
r.const = 0;
r.degenerate = (y == 0);
r.text = sprintf('X(%s) = X(%s) * %g', r.targetname, elemName(lqn, src), y);
r.latex = r.text;
if ~isempty(s)
    r.lhs = s.calltput(cidx);
    r.rhs = s.tput(src) * y;
    r.residual = zeroIfNaN(r.lhs) - zeroIfNaN(r.rhs);
    r.relresidual = r.residual / max(abs(r.rhs), GlobalConstants.FineTol);
end
end

function r = entryflowRecord(lqn, eidx, inc, lam, s)
r = emptyRecord();
r.kind = 'entryflow';
r.target = eidx;
r.targetname = elemName(lqn, eidx);
r.terms = inc;
r.coeff = ones(1, numel(inc));
r.const = lam;
txt = sprintf('X(%s) =', r.targetname);
for i = 1:numel(inc)
    if i == 1
        txt = sprintf('%s X(%s)', txt, callName(lqn, inc(i)));
    else
        txt = sprintf('%s + X(%s)', txt, callName(lqn, inc(i)));
    end
end
if lam > 0
    if isempty(inc)
        txt = sprintf('%s %g', txt, lam);
    else
        txt = sprintf('%s + %g', txt, lam);
    end
    txt = sprintf('%s   (open arrival)', txt);
end
r.text = txt;
r.latex = txt;
if ~isempty(s)
    r.lhs = s.tput(eidx);
    r.rhs = lam;
    for i = 1:numel(inc)
        r.rhs = r.rhs + zeroIfNaN(s.calltput(inc(i)));
    end
    r.residual = zeroIfNaN(r.lhs) - r.rhs;
    r.relresidual = r.residual / max(abs(r.rhs), GlobalConstants.FineTol);
end
end

function r = actflowRecord(lqn, aidx, eidx, v, s)
r = emptyRecord();
r.kind = 'actflow';
r.target = aidx;
r.targetname = elemName(lqn, aidx);
r.terms = eidx;
r.coeff = v;
r.const = 0;
r.degenerate = ~isfinite(v);
r.text = sprintf('X(%s) = X(%s) * %g', r.targetname, elemName(lqn, eidx), v);
r.latex = r.text;
if ~isempty(s)
    r.lhs = s.tput(aidx);
    r.rhs = s.tput(eidx) * v;
    r.residual = zeroIfNaN(r.lhs) - zeroIfNaN(r.rhs);
    r.relresidual = r.residual / max(abs(r.rhs), GlobalConstants.FineTol);
end
end

function r = hostutilRecord(lqn, hidx, s)
r = emptyRecord();
r.kind = 'hostutil';
r.target = hidx;
r.targetname = elemName(lqn, hidx);
r.mult = lqn.mult(hidx);
r.repl = max(1, lqn.repl(hidx));
r.scaled = lqn.sched(hidx) ~= SchedStrategy.INF;
acts = [];
dem = [];
for tidx = lqn.tasksof{hidx}
    for aidx = lqn.actsof{tidx}
        d = lqn.hostdem_mean(aidx);
        if isnan(d) || d == 0
            continue
        end
        acts(end+1) = aidx; %#ok<AGROW>
        dem(end+1) = d; %#ok<AGROW>
    end
end
r.terms = acts;
r.coeff = dem;
% The server count is the declared multiplicity ALONE, not mult*repl: a
% replicated host reports its throughputs and its utilization as TOTALS over the
% copies, so the extra factor would double-count the replication. Checked on the
% two-replica processor of lqn_sockshop, where mult*repl overshoots the reported
% utilization by exactly the replication factor.
m = r.mult;
if ~r.scaled || ~isfinite(m)
    m = 1;
end
r.const = m;
r.branch = ternary(r.scaled, 'queueing', 'inf');
r.degenerate = isempty(acts);
lhs = '';
for i = 1:numel(acts)
    if i == 1
        lhs = sprintf('X(%s)*%g', elemName(lqn, acts(i)), dem(i));
    else
        lhs = sprintf('%s + X(%s)*%g', lhs, elemName(lqn, acts(i)), dem(i));
    end
end
if isempty(lhs)
    lhs = '0';
end
r.text = sprintf('%s = %s*U(%s)', lhs, num2strc(m), r.targetname);
r.latex = r.text;
if ~isempty(s) && ~isnan(s.un(hidx))
    r.lhs = 0;
    for i = 1:numel(acts)
        r.lhs = r.lhs + zeroIfNaN(s.tput(acts(i))) * dem(i);
    end
    r.rhs = m * s.un(hidx);
    r.residual = r.lhs - r.rhs;
    r.relresidual = r.residual / max(abs(r.rhs), GlobalConstants.FineTol);
end
end

% ========================================================================
% structural helpers
% ========================================================================

% Expected number of executions of every activity per invocation of its entry.
% The activity precedence arcs of a task are a transient Markov chain whose
% absorbing state is the reply, so the visit counts of the block solve
% v = e0*(I-P)^-1; a loop back-edge of weight 1-1/count returns count, and an
% AND-fork row summing above one returns the branching expectation, both of
% which are what a visit ratio means. Call arcs leave the block and drop out.
function v = actVisits(lqn)
v = zeros(lqn.nidx, 1);
G = joinScaledGraph(lqn);
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    A = lqn.actsof{tidx};
    if isempty(A)
        continue
    end
    A = unique(A(:))';
    pos = zeros(lqn.nidx, 1);
    pos(A) = 1:numel(A);
    P = full(G(A, A));
    M = eye(numel(A)) - P;
    for eidx = lqn.entriesof{tidx}
        bound = find(G(eidx, :));
        bound = bound(ismember(bound, A));
        if isempty(bound)
            continue
        end
        e0 = zeros(1, numel(A));
        e0(pos(bound)) = 1;
        vv = e0 / M;
        v(A) = v(A) + reshape(vv, size(v(A)));
    end
end
end

% The precedence graph with the arcs into every AND-join target divided by the
% number of branches the join waits for.
%
% An AND-JOIN is the one place where flow does not add up: its target executes
% ONCE per fork, not once per branch, so summing the inbound arcs would count it
% as many times as there are branches. Dividing recovers the rate of one branch
% exactly when the branches carry equal rate, the case for a well-formed
% fork/join block. `actpretype` marks the joined PREDECESSORS, so a join target
% is any successor of one.
function G = joinScaledGraph(lqn)
G = lqn.graph;
if ~isfield(lqn, 'actpretype') || isempty(lqn.actpretype)
    return
end
pre = full(lqn.actpretype);
andpre = find(pre(:)' == ActivityPrecedenceType.PRE_AND);
if isempty(andpre)
    return
end
andpre = andpre(andpre <= size(G, 1));
for j = 1:size(G, 2)
    joined = andpre(G(andpre, j) ~= 0);
    if numel(joined) > 1
        G(joined, j) = G(joined, j) / numel(joined);
    end
end
end

function inc = incomingCalls(lqn, eidx)
inc = [];
for cidx = 1:lqn.ncalls
    if lqn.callpair(cidx, 2) == eidx
        inc(end+1) = cidx; %#ok<AGROW>
    end
end
end

function lam = arrivalRate(lqn, eidx)
lam = 0;
if ~isfield(lqn, 'arrival_mean') || eidx > numel(lqn.arrival_mean)
    return
end
m = lqn.arrival_mean(eidx);
if isfinite(m) && m > GlobalConstants.FineTol
    lam = 1 / m;
end
end

% Entry an activity belongs to: the one whose activity block contains it. An
% activity shared by two entries is attributed to the first, matching how
% actVisits accumulates its visit counts.
function eidx = entryOfActivity(lqn, aidx)
eidx = 0;
tidx = lqn.parent(aidx);
if tidx < 1 || tidx > numel(lqn.entriesof)
    return
end
for e = lqn.entriesof{tidx}
    if any(lqn.actsof{e} == aidx)
        eidx = e;
        return
    end
end
end

function tf = hasPhase2(lqn, tidx)
tf = false;
if ~isfield(lqn, 'actphase')
    return
end
for eidx = lqn.entriesof{tidx}
    for aidx = lqn.actsof{eidx}
        a = aidx - lqn.ashift;
        if a >= 1 && a <= numel(lqn.actphase) && lqn.actphase(a) > 1
            tf = true;
            return
        end
    end
end
end

function [Al, Af, Ah] = incidences(lqn, eqs)
Al = sparse(lqn.ntasks, max(1, lqn.ncalls));
Af = sparse(max(1, lqn.ncalls), lqn.nidx);
Ah = sparse(lqn.nhosts, lqn.nidx);
for i = 1:numel(eqs)
    r = eqs(i);
    switch r.kind
        case 'little'
            % the call classes only; an entry class (open arrival, or the
            % self-driven cycle of a reference task) is not a call
            cols = r.terms(~r.termisentry);
            if ~isempty(cols)
                Al(r.target - lqn.tshift, cols) = 1;
            end
        case 'callflow'
            Af(r.target - lqn.cshift, r.terms) = r.coeff;
        case 'hostutil'
            if ~isempty(r.terms)
                Ah(r.target - lqn.hshift, r.terms) = r.coeff;
            end
    end
end
end

% ========================================================================
% solution adapter
% ========================================================================

% Read the iterate vectors a solved SolverLN carries. Any object or struct
% exposing tput/util/thinkt/servt/residt works, which is what makes this usable
% on an external prediction as well as on LINE's own fixed point.
function s = readSolution(lqn, sol, un, visits)
n = lqn.nidx;
s.tput = solvec(sol, 'tput', n);
s.util = solvec(sol, 'util', n);
s.thinkt = solvec(sol, 'thinkt', n);
s.servt = solvec(sol, 'servt', n);
s.residt = solvec(sol, 'residt', n);
s.un = reportedUtil(lqn, sol, un);
s.visits = visits;
% Per-call throughput is not an iterate of SolverLN: a call inherits the rate of
% its dispatching element scaled by the mean call count.
s.calltput = nan(lqn.ncalls, 1);
for cidx = 1:lqn.ncalls
    src = lqn.callpair(cidx, 1);
    y = lqn.callproc_mean(cidx);
    if isnan(y)
        y = 0;
    end
    if src >= 1 && src <= n
        s.calltput(cidx) = s.tput(src) * y;
    end
end
end

% The REPORTED utilization of every element, which is not the same quantity as
% the `util` iterate: the iterate holds a task's utilization as a SERVER in its
% own task layer -- the U that closes the thread-pool cycle -- and it is left at
% zero on a host. The utilization law at a processor is about the reported UN,
% which the caller supplies.
%
% It is NOT read off the solver here. getEnsembleAvg re-enters iterate() in every
% codebase, and a diagnostic must not re-run a fixed point as a side effect of
% being asked a question. Pass UN explicitly:
%
%   [QN,UN] = solver.getEnsembleAvg();
%   out = lqn_balance_equations(lqn, solver, UN);
function v = reportedUtil(lqn, sol, un)
n = lqn.nidx;
v = nan(n, 1);
if ~isempty(un)
    v = solvec(struct('un', un), 'un', n);
    return
end
for fname = {'un', 'UN'}
    if isstruct(sol) && isfield(sol, fname{1})
        v = solvec(sol, fname{1}, n);
        return
    end
end
end

function v = solvec(sol, fname, n)
v = nan(n, 1);
if isobject(sol)
    if ~isprop(sol, fname)
        return
    end
elseif ~isfield(sol, fname)
    return
end
w = sol.(fname);
if isempty(w)
    return
end
w = w(:);
m = min(numel(w), n);
v(1:m) = w(1:m);
end

% ========================================================================
% reporting
% ========================================================================
function txt = conventionText()
txt = {
    'Conventions:'
    '  little   : rates and populations PER REPLICA (X = tput/repl, N = mult of one copy).'
    '             B(t,k) is the per-class utilization in JOB units; for a queueing task'
    '             B = mult*U with U in [0,1], at an infinite server B = U directly.'
    '             S(k) is the entry SERVICE time (phase 1 + phase 2), the thread hold time.'
    '  other    : throughputs and utilizations as the solver reports them, totalled over replicas.'
    '  N(t)     : lqn.mult. SolverLN iterates on njobs (interlocking corrections, maxmult under'
    '             replication), reported per record as mult/maxmult.'
    };
end

function txt = report(lqn, out, hasSol)
txt = {};
txt{end+1} = sprintf('LQN balance equations: %d hosts, %d tasks, %d entries, %d activities, %d calls', ...
    lqn.nhosts, lqn.ntasks, lqn.nentries, lqn.nacts, lqn.ncalls);
txt = [txt, conventionText()'];
kinds = {'little', 'callflow', 'entryflow', 'actflow', 'hostutil'};
titles = {'thread-pool Little''s law', 'call-flow balance', ...
    'entry-flow balance', 'activity-flow balance', 'host utilization law'};
for k = 1:numel(kinds)
    sel = find(strcmp({out.eqs.kind}, kinds{k}));
    if isempty(sel)
        continue
    end
    txt{end+1} = ''; %#ok<AGROW>
    txt{end+1} = sprintf('--- %s  (kind=''%s'', %d relations) ---', titles{k}, kinds{k}, numel(sel)); %#ok<AGROW>
    for i = sel
        r = out.eqs(i);
        head = sprintf('[%3d] %-14s', i, r.targetname);
        if strcmp(r.kind, 'little')
            head = sprintf('%s %-9s mult=%s repl=%g', head, r.branch, num2strc(r.mult), r.repl);
        elseif strcmp(r.kind, 'hostutil')
            head = sprintf('%s %-9s mult=%s repl=%g', head, r.branch, num2strc(r.mult), r.repl);
        elseif ~isempty(r.branch)
            head = sprintf('%s %-9s', head, r.branch);
        end
        txt{end+1} = head; %#ok<AGROW>
        lines = strsplit(r.text, '\n');
        for j = 1:numel(lines)
            txt{end+1} = sprintf('       %s', lines{j}); %#ok<AGROW>
        end
        notes = {};
        if r.corr.phase2
            notes{end+1} = 'phase-2 tail on Z'; %#ok<AGROW>
        end
        if r.corr.setup
            notes{end+1} = 'setup charge on Z'; %#ok<AGROW>
        end
        if r.degenerate
            notes{end+1} = 'DEGENERATE (not usable as a residual)'; %#ok<AGROW>
        end
        if r.clamped
            notes{end+1} = 'SATURATED (equality unattainable, use a hinge)'; %#ok<AGROW>
        end
        if hasSol && ~r.degenerate && isnan(r.residual)
            notes{end+1} = 'not instantiated (pass UN for the host utilization law)'; %#ok<AGROW>
        end
        if ~isempty(notes)
            txt{end+1} = sprintf('       note: %s', strjoin(notes, ', ')); %#ok<AGROW>
        end
        if hasSol && ~r.degenerate && ~isnan(r.residual)
            txt{end+1} = sprintf('       lhs=%.6g  rhs=%.6g  residual=%.3e  rel=%.3e', ...
                r.lhs, r.rhs, r.residual, r.relresidual); %#ok<AGROW>
        end
    end
end
if hasSol
    txt{end+1} = ''; %#ok<AGROW>
    txt{end+1} = sprintf('max |residual| over %d non-degenerate relations: %.3e', ...
        sum(~[out.eqs.degenerate]), out.maxresidual); %#ok<AGROW>
end
end

% ========================================================================
% small utilities
% ========================================================================
function s = num2strc(x)
if isinf(x)
    s = 'Inf';
elseif x == fix(x)
    s = sprintf('%d', x);
else
    s = sprintf('%g', x);
end
end

function x = zeroIfNaN(x)
if isnan(x)
    x = 0;
end
end

function y = ternary(c, a, b)
if c
    y = a;
else
    y = b;
end
end

function v = multField(lqn, fname, tidx)
v = NaN;
if isfield(lqn, fname) && tidx <= numel(lqn.(fname))
    v = lqn.(fname)(tidx);
end
end

function nm = elemName(lqn, idx)
if idx >= 1 && idx <= numel(lqn.hashnames) && ~isempty(lqn.hashnames{idx})
    nm = lqn.hashnames{idx};
else
    nm = sprintf('#%d', idx);
end
end

function nm = callName(lqn, cidx)
if cidx >= 1 && cidx <= numel(lqn.callhashnames) && ~isempty(lqn.callhashnames{cidx})
    nm = lqn.callhashnames{cidx};
else
    nm = sprintf('C%d', cidx);
end
end

function nm = termName(lqn, isEntry, idx)
if isEntry
    nm = elemName(lqn, idx);
else
    nm = callName(lqn, idx);
end
end

function nm = callTypeName(ct)
switch ct
    case CallType.SYNC
        nm = 'sync';
    case CallType.ASYNC
        nm = 'async';
    case CallType.FWD
        nm = 'fwd';
    otherwise
        nm = '';
end
end
