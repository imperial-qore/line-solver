% https://github.com/imperial-qore/line-solver/wiki/Getting-started#example-6-a-queueing-network-with-caching

% Block 1: nodes
model = Network('model');

% Block 1: nodes
clientDelay = Delay(model, 'Client');
cacheNode = Cache(model, 'Cache', 1000, 50, ReplacementStrategy.LRU);
cacheDelay = Delay(model, 'CacheDelay');

% Block 2: classes
clientClass = ClosedClass(model, 'ClientClass', 1, clientDelay, 0);
hitClass = ClosedClass(model, 'HitClass', 0, clientDelay, 0);
missClass = ClosedClass(model, 'MissClass', 0, clientDelay, 0);

clientDelay.setService(clientClass, Immediate());
cacheDelay.setService(hitClass, Exp.fitMean(0.2));
cacheDelay.setService(missClass, Exp.fitMean(1));

cacheNode.setRead(clientClass, Zipf(1.4,1000));
cacheNode.setHitClass(clientClass, hitClass);
cacheNode.setMissClass(clientClass, missClass);

% Block 3: topology
P = model.initRoutingMatrix;
% routing from client to cache
P{clientClass, clientClass}(clientDelay, cacheNode)=1;
% routing out of the cache
P{hitClass, hitClass}(cacheNode, cacheDelay)=1;
P{missClass, missClass}(cacheNode, cacheDelay)=1;
% return to the client
P{hitClass, clientClass}(cacheDelay, clientDelay)=1;
P{missClass, clientClass}(cacheDelay, clientDelay)=1;
% routing from cacheNode
model.link(P);

% Block 4: solution
ssaAvgTable = SolverSSA(model,'samples',2e4,'seed',1,'verbose',true).getAvgTable
if license('test','parallel_computing_toolbox')
    ssaAvgTablePara = SolverSSA(model,'samples',2e4,'seed',1,'verbose',true,'parallel').getAvgTable
end