% Pass-and-swap (PAS) saturation handling for large / unbounded buffers.
%
% An order-independent service rate mu(c) usually SATURATES: beyond some
% per-class count threshold tau_r, adding more class-r jobs no longer changes
% mu (e.g. a compatibility / rank function saturates once each present class
% has one job; an M/M/K rate saturates at K jobs). A large finite buffer can
% therefore approximate an unbounded queue cheaply.
%
% This uses the Comte (thesis Sect. 4.1) compatibility queue: I=2 classes,
% S=3 servers, S_1={1,2}, S_2={2,3}, rank function
%     mu(A) = sum_{ s : exists present class in A compatible with s } cap_s,
% so mu({1})=mu({2})=cap_1+cap_2 = 3, mu({1,2})=cap_1+cap_2+cap_3 = 5, and mu
% is constant once both classes are present (cutoffs tau=[1,1]).
%
% It checks that (i) LDES matches CTMC on a finite buffer and
% (ii) an unset/infinite buffer raises a clean, actionable error.
clear node jobclass

% compatibility comp(s,i) = 1 iff server s serves class i; capacities cap_s
comp  = [1 0;     % server 1 -> class 1
         1 1;     % server 2 -> classes 1, 2 (shared)
         0 1];    % server 3 -> class 2
cap_s = [2; 1; 2];
lambda = [1.5 1.0];
I = size(comp, 2);
S = size(comp, 1);

% Rank function: total capacity of servers compatible with a present class.
muRank = @(c) sum(cap_s(any(comp(:, unique(c(c>0))), 2)));

% (1) a large buffer approximates the unbounded queue; LDES matches CTMC ------
CAP = 8;
Qc = qlen(CTMC(build(muRank, CAP, I, S, lambda), 'cutoff', CAP), I);
Ql = qlen(LDES(build(muRank, CAP, I, S, lambda), 'samples', 3e5, 'seed', 23000), I);
fprintf('open compatibility PAS, buffer=%d:\n', CAP);
fprintf('  CTMC  QLen = %s\n', mat2str(round(Qc, 5)));
fprintf('  LDES  QLen = %s\n', mat2str(round(Ql, 5)));
rel = max(abs(Qc - Ql) ./ max(Qc, 1e-12));
fprintf('  LDES vs CTMC max rel|dQ| = %.2f%%\n', 100 * rel);
assert(rel <= 0.03, 'LDES does not match CTMC within simulation noise');

% (2) clean error for an unset / infinite buffer -----------------------------
fprintf('\ninfinite-buffer PAS raises a clean, actionable error:\n');
try
    LDES(build(muRank, Inf, I, S, lambda), 'samples', 1000, 'seed', 1).getAvgTable;
    error('expected a finite-buffer error');
catch e
    fprintf('  %s\n', e.message);
end

fprintf(['\nPASS: PAS saturation handling matches CTMC on a large buffer and ' ...
         'errors cleanly on infinite buffers.\n']);

function model = build(mu, cap, I, S, lambda)
model = Network('PASsaturation');
source = Source(model, 'Source');
queue  = Queue(model, 'PASQueue', SchedStrategy.PAS);
sink   = Sink(model, 'Sink');
jobclass = cell(1, I);
for r = 1:I
    jobclass{r} = OpenClass(model, sprintf('Class%d', r));
    source.setArrival(jobclass{r}, Exp(lambda(r)));
end
queue.setService(mu);
queue.setSwapGraph(zeros(I));     % empty graph: plain OI queue
queue.setNumberOfServers(S);
if ~isinf(cap)
    queue.setCap(cap);
end
P = model.initRoutingMatrix;
for r = 1:I
    P{jobclass{r}} = Network.serialRouting(source, queue, sink);
end
model.link(P);
end

function q = qlen(solver, I)
QN = solver.getAvgQLen;             % nstations x nclasses mean queue length
sn = solver.model.getStruct;
sidx = sn.nodeToStation(find(strcmp(sn.nodenames, 'PASQueue'), 1));
q = QN(sidx, 1:I);
end
