% PS variant of retrieval_simple: a single processor-sharing retrieval station.
clc; clear solver AvgTable;

accessProb   = [49, 49, 49, 49, 7, 1, 1] / 205;   % per-item arrival probabilities (7 items)

model = Network('DelayedHits');

n = numel(accessProb);            % number of items
capacity = [6];                   % per-level cache capacity (total m = 6)

source    = Source(model, 'Source');
cacheNode = Cache(model, 'Cache', n, capacity, ReplacementStrategy.RR);
queue     = Queue(model, 'Queue', SchedStrategy.PS);
sink      = Sink(model, 'Sink');

jobClass  = OpenClass(model, 'InitClass', 0);
hitClass  = OpenClass(model, 'HitClass', 0);
missClass = OpenClass(model, 'MissClass', 0);

source.setArrival(jobClass, Exp(1));

% Read class service at the retrieval queue = default per-item fetch service.
queue.setService(jobClass, Exp(1.0));

pAccess = DiscreteSampler(accessProb);
cacheNode.setRead(jobClass, pAccess);
cacheNode.setHitClass(jobClass, hitClass);
cacheNode.setMissClass(jobClass, missClass);

% No serviceRates and no routingMatrices: both inherited from the read class.
cacheNode.setRetrievalSystem(jobClass, missClass, queue);

P = model.initRoutingMatrix();
P{jobClass, jobClass}(source, cacheNode)  = 1.0;
% Retrieval topology drawn once for the read class: cache <-> queue.
P{jobClass, jobClass}(cacheNode, queue)   = 1.0;
P{jobClass, jobClass}(queue, cacheNode)   = 1.0;
P{hitClass, hitClass}(cacheNode, sink)    = 1.0;
P{missClass, missClass}(cacheNode, sink)  = 1.0;
model.link(P);

% Simulation
SSA(model, 'samples', 5000, 'method', 'serial', 'seed', 1).getAvgCacheTable
LDES(model, 'samples', 1e6, 'seed', 1).getAvgCacheTable

% Analytical
MVA(model).getAvgCacheTable
NC(model).getAvgCacheTable

% Item-level cache occupancy
MVA(model).getAvgItemTable
NC(model).getAvgItemTable
