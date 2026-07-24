# %%
from line_solver import *
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
model = Network('Model')
# Block 1: nodes
clientDelay = Delay(model, 'Client')
cacheNode = Cache(model, 'Cache', 1000, 50, ReplacementStrategy.LRU)
cacheDelay = Delay(model, 'CacheDelay')
# %%
# Block 2: classes
clientClass = ClosedClass(model, 'ClientClass', 1, clientDelay, 0)
hitClass = ClosedClass(model, 'HitClass', 0, clientDelay, 0)
missClass = ClosedClass(model, 'MissClass', 0, clientDelay, 0)

clientDelay.set_service(clientClass, Immediate())
cacheDelay.set_service(hitClass, Exp.fit_mean(0.2))
cacheDelay.set_service(missClass, Exp.fit_mean(1.0))

cacheNode.set_read(clientClass, Zipf(1.4, 1000))
cacheNode.set_hit_class(clientClass, hitClass)
cacheNode.set_miss_class(clientClass, missClass)
# %%
# Block 3: topology
P = model.init_routing_matrix()
# routing from client to cache
P.set(clientClass, clientClass, clientDelay, cacheNode, 1.0)
# routing out of the cache
P.set(hitClass, hitClass, cacheNode, cacheDelay, 1.0)
P.set(missClass, missClass, cacheNode, cacheDelay, 1.0)
# return to the client
P.set(hitClass, clientClass, cacheDelay, clientDelay, 1.0)
P.set(missClass, clientClass, cacheDelay, clientDelay, 1.0)
# routing from cacheNode
model.link(P)
# %%
ssaAvgTablePara = SSA(model, samples=20000, seed=1, verbose=True, method='serial').avg_table()
print(ssaAvgTablePara)