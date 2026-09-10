% Cache with probabilistic per-item routing through a retrieval system.
clc; clear solver AvgTable;

accessProb   = [0.6, 0.3, 0.1];   % per-item access probabilities (3 items)

model = Network('DelayedHits');

n = numel(accessProb);            % number of items
capacity = [2];                   % per-level cache capacity

source    = Source(model, 'Source');
cacheNode = Cache(model, 'Cache', n, capacity, ReplacementStrategy.FIFO);
isQueue   = Queue(model, 'IS Queue', SchedStrategy.INF);
queue1    = Queue(model, 'Queue 1', SchedStrategy.FCFS);
queue2    = Queue(model, 'Queue 2', SchedStrategy.FCFS);
queues    = {isQueue, queue1, queue2};
sink      = Sink(model, 'Sink');

jobClass  = OpenClass(model, 'InitClass', 0);
hitClass  = OpenClass(model, 'HitClass', 0);
missClass = OpenClass(model, 'MissClass', 0);

source.setArrival(jobClass, Exp(1));

% Read class service per queue = default per-item fetch service (uniform over items).
isQueue.setService(jobClass, Exp(2.0));
queue1.setService(jobClass, Exp(3.0));
queue2.setService(jobClass, Exp(3.0));

pAccess = DiscreteSampler(accessProb);
cacheNode.setRead(jobClass, pAccess);
cacheNode.setHitClass(jobClass, hitClass);
cacheNode.setMissClass(jobClass, missClass);

cacheNode.setRetrievalSystem(jobClass, missClass, queues);

% Per-item routing over [IS(1), Queue1(2), Queue2(3), Cache(4)], applied via the
% per-item routing methods. Index nQ+1 (last row/col) is the cache; row=from, col=to.
routingMatrices = { ...
    [0.00, 0.50, 0.00, 0.50;   % from IS
     0.00, 0.00, 0.70, 0.30;   % from Queue1
     0.00, 0.00, 0.00, 1.00;   % from Queue2
     0.70, 0.30, 0.00, 0.00];  % from Cache
    [0.00, 0.30, 0.00, 0.70;
     0.00, 0.00, 0.50, 0.50;
     0.00, 0.00, 0.00, 1.00;
     0.20, 0.80, 0.00, 0.00];
    [0.00, 0.60, 0.00, 0.40;
     0.00, 0.00, 0.40, 0.60;
     0.00, 0.00, 0.00, 1.00;
     0.50, 0.50, 0.00, 0.00]};

nQ = numel(queues);
for item = 1:n
    R = routingMatrices{item};
    for a = 1:nQ
        cacheNode.setItemRoutingProb(jobClass, item, cacheNode, queues{a}, R(nQ+1, a)); % cache -> queue a
        cacheNode.setItemRoutingProb(jobClass, item, queues{a}, cacheNode, R(a, nQ+1));  % queue a -> cache
        for b = 1:nQ
            cacheNode.setItemRoutingProb(jobClass, item, queues{a}, queues{b}, R(a, b));
        end
    end
end

P = model.initRoutingMatrix();
P{jobClass, jobClass}(source, cacheNode)  = 1.0;
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
