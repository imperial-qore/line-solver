#!/usr/bin/env python3
"""Gallery Example: gallery_cache_routing"""

from line_solver import *

def gallery_cache_routing():
    """Open cache with hit/miss routed to distinct queues."""
    model = Network('Cache-Routing')
    n = 4
    m = 2
    source = Source(model, 'Source')
    cacheNode = Cache(model, 'Cache', n, m, ReplacementStrategy.LRU)
    hitQueue = Queue(model, 'HitQueue', SchedStrategy.FCFS)
    missQueue = Queue(model, 'MissQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    jobClass = OpenClass(model, 'InitClass', 0)
    hitClass = OpenClass(model, 'HitClass', 0)
    missClass = OpenClass(model, 'MissClass', 0)
    source.setArrival(jobClass, Exp(1))
    hitQueue.setService(hitClass, Exp(2.0))
    missQueue.setService(missClass, Exp(1.0))
    pAccess = DiscreteSampler([1.0 / n] * n)
    cacheNode.setRead(jobClass, pAccess)
    cacheNode.setHitClass(jobClass, hitClass)
    cacheNode.setMissClass(jobClass, missClass)
    P = model.init_routing_matrix()
    P.set(jobClass, jobClass, source, cacheNode, 1.0)
    P.set(hitClass, hitClass, cacheNode, hitQueue, 1.0)
    P.set(hitClass, hitClass, hitQueue, sink, 1.0)
    P.set(missClass, missClass, cacheNode, missQueue, 1.0)
    P.set(missClass, missClass, missQueue, sink, 1.0)
    model.link(P)
    return model



if __name__ == '__main__':
    model = gallery_cache_routing()
    print('Model built:', type(model).__name__)
