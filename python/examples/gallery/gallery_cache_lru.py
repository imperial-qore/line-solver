#!/usr/bin/env python3
"""Gallery Example: gallery_cache_lru"""

from line_solver import *

def gallery_cache_lru():
    """Closed cache model with LRU replacement (n=5 items, m=2 slots)."""
    model = Network('Cache-LRU')
    n = 5
    m = 2
    delay = Delay(model, 'Delay')
    cacheNode = Cache(model, 'Cache', n, m, ReplacementStrategy.LRU)
    jobClass = ClosedClass(model, 'JobClass', 1, delay, 0)
    hitClass = ClosedClass(model, 'HitClass', 0, delay, 0)
    missClass = ClosedClass(model, 'MissClass', 0, delay, 0)
    delay.setService(jobClass, Exp(1))
    pAccess = DiscreteSampler([1.0 / n] * n)
    cacheNode.setRead(jobClass, pAccess)
    cacheNode.setHitClass(jobClass, hitClass)
    cacheNode.setMissClass(jobClass, missClass)
    P = model.init_routing_matrix()
    P.set(jobClass, jobClass, delay, cacheNode, 1.0)
    P.set(hitClass, jobClass, cacheNode, delay, 1.0)
    P.set(missClass, jobClass, cacheNode, delay, 1.0)
    model.link(P)
    return model



if __name__ == '__main__':
    model = gallery_cache_lru()
    print('Model built:', type(model).__name__)
