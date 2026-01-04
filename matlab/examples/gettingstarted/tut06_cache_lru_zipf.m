% Example 6: A queueing network with caching

% Block 1: nodes
model = Network('Model');

nItems = 1000;
cacheSize = 50;
% Block 1: nodes
clientDelay = Delay(model, 'Client');
cacheNode = Cache(model, 'Cache', nItems, cacheSize, ReplacementStrategy.LRU());
cacheDelay = Delay(model, 'CacheDelay');

% Block 2: classes
clientClass = ClosedClass(model, 'ClientClass', 1, clientDelay, 0);
hitClass = ClosedClass(model, 'HitClass', 0, clientDelay, 0);
missClass = ClosedClass(model, 'MissClass', 0, clientDelay, 0);

clientDelay.setService(clientClass, Immediate());
cacheDelay.setService(hitClass, Exp.fitMean(0.2));
cacheDelay.setService(missClass, Exp.fitMean(1));

cacheNode.setRead(clientClass, Zipf(1.4,nItems));
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
ssaAvgTable = SSA(model,'samples',2e4,'seed',23000,'verbose',true).avgTable()
if license('test','parallel_computing_toolbox')
    ssaAvgTablePara = SSA(model,'samples',2e4,'seed',23000,'verbose',true,'parallel').avgTable()
end