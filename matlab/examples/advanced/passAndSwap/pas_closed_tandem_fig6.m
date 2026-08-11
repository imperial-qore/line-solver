function pas_closed_tandem_fig6()
% PAS_CLOSED_TANDEM_FIG6  Closed tandem of two pass-and-swap (PAS) queues,
% reproducing Figures 5 and 6 of Comte & Dorsman, "Pass-and-Swap Queues"
% (2021, arXiv:2009.12299; Queueing Systems).
%
%   Topology (Fig. 6):   --> PASQueue1 --> PASQueue2 -->   (closed tandem)
%
% There are I=6 customer classes with a SINGLE customer of each class. Both
% queues share the swapping graph of Figure 5 (undirected; classes as nodes):
%
%        6                edges: 1-3, 1-4, 2-4, 2-5, 3-6, 4-6, 5-6
%      / | \              (class 6 at the top, classes 1,2 at the bottom).
%     3  4  5
%      \/ \/
%      1   2
%
% As in the paper's illustration, only the customer at the HEAD of each queue
% receives a positive service rate (single server): mu(c)=mu1 at queue 1 and
% nu(d)=mu2 at queue 2, independent of class so that the total rate depends
% only on the customer multiset (order-independence). Upon a service
% completion the departing customer is selected by the pass-and-swap scan over
% the swapping graph and is routed to the back of the OTHER queue as a
% customer of the same class.
%
% The model is solved exactly by CTMC and by the LDES simulator,
% and both are validated against an independent product-form reference that
% builds the Markov chain directly from the pass-and-swap algorithm (Theorem 5).
%
% NOTE: a closed PAS network with a non-empty swapping graph has a reducible
% generator (the pass-and-swap mechanism conserves the placement order, paper
% Prop. 2/4), so the stationary distribution is that of the recurrent
% component containing the initial state c=(1,2,3,4,5,6) of Figure 6a.

mu1 = 1.0;  mu2 = 1.3;          % head-only single-server rates at queue 1, 2
edges = [1 3; 1 4; 2 4; 2 5; 3 6; 4 6; 5 6];
G = zeros(6);
for k = 1:size(edges,1)
    G(edges(k,1),edges(k,2)) = 1;
    G(edges(k,2),edges(k,1)) = 1;
end

% ---- LINE model -------------------------------------------------------
model = Network('PASclosedTandem');
q1 = Queue(model, 'PASQueue1', SchedStrategy.PAS);
q2 = Queue(model, 'PASQueue2', SchedStrategy.PAS);
jobclass = cell(1,6);
for r = 1:6
    jobclass{r} = ClosedClass(model, sprintf('Class%d', r), 1, q1);
end
q1.setService(@(c) mu1);       % head-only: mu(prefix)=mu1 => Delta mu>0 only at head
q2.setService(@(c) mu2);
q1.setSwapGraph(G);  q1.setNumberOfServers(1);  q1.setCap(6);
q2.setSwapGraph(G);  q2.setNumberOfServers(1);  q2.setCap(6);
P = model.initRoutingMatrix;
for r = 1:6
    P{jobclass{r}}(q1, q2) = 1.0;
    P{jobclass{r}}(q2, q1) = 1.0;
end
model.link(P);

% A closed PAS network with a non-empty swapping graph is reducible (one
% recurrent component per placement order), so the initial job placement is a
% REQUIRED model input -- there is no valid default. We use the initial state
% of Figure 6a: queue 1 = (1,2,3,4,5,6) (class 1 oldest, at the head), queue 2
% empty. Solving with a reversed placement would select the mirror component.
q1.setState([1 2 3 4 5 6]);

fprintf('=== CTMC (exact) ===\n');
Tc = CTMC(model, 'cutoff', 6).getAvgTable;
disp(Tc);

fprintf('=== LDES (simulation, 6e5 samples) ===\n');
Tl = LDES(model, 'samples', 6e5, 'seed', 23000, 'verbose', false).getAvgTable;

% ---- independent product-form reference -------------------------------
[Qref1, Qref2] = reference_pas_tandem(G, [mu1 mu2]);
Qref = [Qref1, Qref2]';        % column order: q1 classes 1..6 then q2 classes 1..6

qc = Tc.QLen;  ql = Tl.QLen;
fprintf('\n  station    class   reference     CTMC       LDES\n');
for i = 1:height(Tc)
    fprintf('  %-9s  %-6s %9.5f %9.5f %9.5f\n', string(Tc.Station(i)), ...
        string(Tc.JobClass(i)), Qref(i), qc(i), ql(i));
end
errCTMC = max(abs(qc - Qref));
errLDES = max(abs(ql - Qref));
fprintf('\nCTMC vs reference:  max|dQ| = %.3e\n', errCTMC);
fprintf('LDES vs reference:  max|dQ| = %.3e (simulation noise)\n', errLDES);
assert(errCTMC <= 1e-9,  'CTMC does not match the pass-and-swap reference');
assert(errLDES <= 1e-2,  'LDES does not match the pass-and-swap reference');
fprintf('PASS: CTMC matches the pass-and-swap reference; LDES agrees within noise.\n');
end

% =======================================================================
function [Q1, Q2] = reference_pas_tandem(G, mu)
% Independent reference: build the CTMC of the closed PAS tandem directly from
% the pass-and-swap algorithm, by reachability from the initial state of Fig. 6a
% (queue 1 = (1,...,6) ascending, queue 2 empty), head-only service.
n = size(G,1);
init1 = 1:n;
states = {};  idx = containers.Map('KeyType','char','ValueType','double');
    function id = getid(l1, l2)
        key = [sprintf('%d,', l1), '|', sprintf('%d,', l2)];
        if isKey(idx, key)
            id = idx(key);
        else
            id = numel(states) + 1;  idx(key) = id;  states{id} = {l1, l2};
        end
    end
getid(init1, []);
fr = 1;  E = zeros(0,3);
while fr <= numel(states)
    l1 = states{fr}{1};  l2 = states{fr}{2};
    if ~isempty(l1)                                   % head of queue 1 completes
        [nl, dep] = psalgorithm(l1, 1, G);
        j = getid(nl, [l2, dep]);  E(end+1,:) = [fr, j, mu(1)]; %#ok<AGROW>
    end
    if ~isempty(l2)                                   % head of queue 2 completes
        [nl, dep] = psalgorithm(l2, 1, G);
        j = getid([l1, dep], nl);  E(end+1,:) = [fr, j, mu(2)]; %#ok<AGROW>
    end
    fr = fr + 1;
end
ns = numel(states);  Qm = zeros(ns);
for e = 1:size(E,1)
    Qm(E(e,1),E(e,2)) = Qm(E(e,1),E(e,2)) + E(e,3);
end
for i = 1:ns, Qm(i,i) = -sum(Qm(i,:)); end
pivec = ([Qm'; ones(1,ns)] \ [zeros(ns,1); 1])';
Q1 = zeros(1,n);  Q2 = zeros(1,n);
for i = 1:ns
    l1 = states{i}{1};  l2 = states{i}{2};
    for r = 1:n
        Q1(r) = Q1(r) + pivec(i)*sum(l1==r);
        Q2(r) = Q2(r) + pivec(i)*sum(l2==r);
    end
end
end

% =======================================================================
function [newlist, dep] = psalgorithm(list, p, G)
% Pass-and-swap algorithm: the customer at position p completes; it scans toward
% the back for the first swappable class, takes its place, the ejected customer
% repeats, and the last ejected customer departs. Classes shift one step along
% the chain and the served slot is removed.
nn = numel(list);  chain = p;  moving = list(p);  cur = p;
while true
    q = 0;
    for j = cur+1:nn
        if G(moving, list(j)), q = j; break; end
    end
    if q == 0, break; end
    chain(end+1) = q;  moving = list(q);  cur = q; %#ok<AGROW>
end
dep = list(chain(end));
newlist = list;
for i = 1:numel(chain)-1
    newlist(chain(i+1)) = list(chain(i));
end
newlist(chain(1)) = [];
end
