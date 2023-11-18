#!/usr/bin/env python
# coding: utf-8

# https://github.com/imperial-qore/line-solver/wiki/Getting-started#example-6-a-queueing-network-with-caching

# In[ ]:


from line_solver import *


# In[ ]:


# Block 1: nodes
model = Network('model')

# Block 1: nodes
clientDelay = Delay(model, 'Client')
cacheNode = Cache(model, 'Cache', 1000, 50, ReplacementStrategy.LRU)
cacheDelay = Delay(model, 'CacheDelay')

# Block 2: classes
clientClass = ClosedClass(model, 'ClientClass', 1, clientDelay, 0)
hitClass = ClosedClass(model, 'HitClass', 0, clientDelay, 0)
missClass = ClosedClass(model, 'MissClass', 0, clientDelay, 0)
clientDelay.setService(clientClass, Immediate())
cacheDelay.setService(hitClass, Exp.fitMean(0.2))
cacheDelay.setService(missClass, Exp.fitMean(1))
cacheNode.setRead(clientClass, Zipf(1.4, 1000))
cacheNode.setHitClass(clientClass, hitClass)
cacheNode.setMissClass(clientClass, missClass)

# Block 3: topology
P = model.initRoutingMatrix()
# routing from client to cache
P.set(clientClass, clientClass, clientDelay, cacheNode, 1.0)
# routing out of the cache
P.set(hitClass, hitClass, cacheNode, cacheDelay, 1.0)
P.set(missClass, missClass, cacheNode, cacheDelay, 1.0)
# # return to the client
P.set(hitClass, clientClass, cacheDelay, clientDelay, 1.0)
P.set(missClass, clientClass, cacheDelay, clientDelay, 1.0)
# # routing from cacheNode
model.link(P)

# SSA not working on this case in Java but MVA works fine
# Block 4: solution
# ssaAvgTable = SolverSSA(model,'samples',20000,'seed',1,'verbose',True).getAvgTable()
# print(ssaAvgTable)

# ssaAvgTablePara = SolverSSA(model,'samples',20000,'seed',1,'verbose',True,'parallel').getAvgTable()
# print(ssaAvgTable)

