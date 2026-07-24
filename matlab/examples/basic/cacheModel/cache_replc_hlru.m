clear solver AvgTable;

% h-LRU / LRU(m) replacement: h LRU lists of capacities m(1..h); a miss
% inserts the item at the head of list 1, a hit in list l exchanges the item
% with the tail of list l+1. Exact in CTMC; simulated in SSA/LDES; MVA uses
% the characteristic-time (TTL) approximation of Gast and Van Houdt
% (SIGMETRICS 2015), which reduces to the Che approximation for h=1.

model = Network('model');

n = 6;        % number of items
m = [2 1];    % list capacities: list 1 holds 2 items, list 2 holds 1

source = Source(model, 'Source');
cacheNode = Cache(model, 'Cache', n, m, ReplacementStrategy.HLRU);
sink = Sink(model, 'Sink');

jobClass = OpenClass(model, 'InitClass', 0);
hitClass = OpenClass(model, 'HitClass', 0);
missClass = OpenClass(model, 'MissClass', 0);

source.setArrival(jobClass, Exp(1));

pAccess = Zipf(1.2, n);  % Zipf-like item references
cacheNode.setRead(jobClass, pAccess);

cacheNode.setHitClass(jobClass, hitClass);
cacheNode.setMissClass(jobClass, missClass);

P = model.initRoutingMatrix;
P{jobClass, jobClass}(source, cacheNode) = 1.0;
P{hitClass, hitClass}(cacheNode, sink) = 1.0;
P{missClass, missClass}(cacheNode, sink) = 1.0;
model.link(P);

solver{1} = CTMC(model,'keep',false,'cutoff',1);   % exact
AvgTable{1} = solver{1}.getAvgNodeTable; AvgTable{1}

model.reset;
solver{2} = MVA(model);                            % TTL approximation
AvgTable{2} = solver{2}.getAvgNodeTable; AvgTable{2}

model.reset;
solver{3} = SSA(model,'samples',1e4,'method','serial','seed',23000);
AvgTable{3} = solver{3}.getAvgNodeTable; AvgTable{3}
