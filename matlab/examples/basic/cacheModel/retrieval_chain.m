% Cache with a chained retrieval system (delayed hits).
clc; clear solver AvgTable;

accessProb   = [0.6, 0.3, 0.1];   % per-item access probabilities (3 items)
nQueues      = 2;

model = Network('DelayedHits');

n = numel(accessProb);            % number of items
capacity = [1];                   % per-level cache capacity

source    = Source(model, 'Source');
cacheNode = Cache(model, 'Cache', n, capacity, ReplacementStrategy.FIFO);
queues = cell(1, nQueues);
for i = 1:nQueues
    queues{i} = Queue(model, sprintf('Queue_%d', i), SchedStrategy.FCFS);
end
sink      = Sink(model, 'Sink');

jobClass  = OpenClass(model, 'InitClass', 0);
hitClass  = OpenClass(model, 'HitClass', 0);
missClass = OpenClass(model, 'MissClass', 0);

source.setArrival(jobClass, Exp(1));

% Read class service at each retrieval queue = default per-item fetch service.
queues{1}.setService(jobClass, Exp(2.0));
queues{2}.setService(jobClass, Exp(3.0));

pAccess = DiscreteSampler(accessProb);
cacheNode.setRead(jobClass, pAccess);
cacheNode.setHitClass(jobClass, hitClass);
cacheNode.setMissClass(jobClass, missClass);

% No serviceRates and no routingMatrices: both inherited from the read class.
cacheNode.setRetrievalSystem(jobClass, missClass, queues);

P = model.initRoutingMatrix();
P{jobClass, jobClass}(source, cacheNode)     = 1.0;
% Retrieval chain drawn once for the read class: Cache -> Queue_1 -> Queue_2 -> Cache
P{jobClass, jobClass}(cacheNode, queues{1})  = 1.0;
P{jobClass, jobClass}(queues{1}, queues{2})  = 1.0;
P{jobClass, jobClass}(queues{2}, cacheNode)  = 1.0;
P{hitClass, hitClass}(cacheNode, sink)       = 1.0;
P{missClass, missClass}(cacheNode, sink)     = 1.0;
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
