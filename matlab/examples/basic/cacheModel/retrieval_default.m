% Retrieval system whose per-item miss routing AND service are taken, by
% default, from the read class: routing from the read class's edges among the
% retrieval queues in the top-level P matrix, and service from the read class's
% service distribution at each queue. Per-item overrides
% (setItemRoutingProb with the cache as source/dest / queue.setItemServiceRate)
% then reconfigure one item at finer granularity.
clc; clear solver AvgTable;

accessProb = [0.6, 0.3, 0.1];   % per-item access probabilities (3 items)

model = Network('DelayedHits');

n = numel(accessProb);          % number of items
capacity = [1];                 % per-level cache capacity

source    = Source(model, 'Source');
cacheNode = Cache(model, 'Cache', n, capacity, ReplacementStrategy.FIFO);
% PS retrieval stations: per-item (class-dependent) service rates are admissible
% in the analytical retrieval algorithm (FCFS/SIRO would require identical rates).
queue1    = Queue(model, 'Queue_1', SchedStrategy.PS);
queue2    = Queue(model, 'Queue_2', SchedStrategy.PS);
sink      = Sink(model, 'Sink');

jobClass  = OpenClass(model, 'InitClass', 0);
hitClass  = OpenClass(model, 'HitClass', 0);
missClass = OpenClass(model, 'MissClass', 0);

source.setArrival(jobClass, Exp(1));

% Read class service at each retrieval queue = default per-item fetch service.
queue1.setService(jobClass, Exp(2.0));
queue2.setService(jobClass, Exp(3.0));

pAccess = DiscreteSampler(accessProb);
cacheNode.setRead(jobClass, pAccess);
cacheNode.setHitClass(jobClass, hitClass);
cacheNode.setMissClass(jobClass, missClass);

% No serviceRates and no routingMatrices: both inherited from the read class.
cacheNode.setRetrievalSystem(jobClass, missClass, {queue1, queue2});

% Item-level overrides for item 1: skip Queue_2 and fetch faster at Queue_1.
cacheNode.setItemRoutingProb(jobClass, 1, queue1, queue2, 0.0);  % delete default edge
cacheNode.setItemRoutingProb(jobClass, 1, queue1, cacheNode, 1.0);  % exit after Queue_1
queue1.setItemServiceRate(cacheNode, jobClass, 1, 5.0);                 % faster item-1 fetch

P = model.initRoutingMatrix();
P{jobClass, jobClass}(source, cacheNode)  = 1.0;
% Default retrieval topology, drawn once for the read class: cache -> Q1 -> Q2 -> cache
P{jobClass, jobClass}(cacheNode, queue1)  = 1.0;   % entry into retrieval
P{jobClass, jobClass}(queue1, queue2)     = 1.0;   % Queue_1 -> Queue_2
P{jobClass, jobClass}(queue2, cacheNode)  = 1.0;   % exit back to cache
P{hitClass, hitClass}(cacheNode, sink)    = 1.0;
P{missClass, missClass}(cacheNode, sink)  = 1.0;
model.link(P);

% Simulation
SSA(model, 'samples', 100000, 'method', 'serial', 'seed', 1).getAvgCacheTable
LDES(model, 'samples', 1e6, 'seed', 1).getAvgCacheTable

% Analytical
MVA(model).getAvgCacheTable
NC(model).getAvgCacheTable

% Item-level cache occupancy
MVA(model).getAvgItemTable
NC(model).getAvgItemTable
