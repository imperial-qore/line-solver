"""
Tandem of Two Caches Sharing One Item Set

A request reads Cache1; on a miss it reads Cache2 for THE SAME item; on a
second miss it is fetched from the origin. Item identity is carried by giving
each item its OWN job class, so no class switching happens on an arc into a
cache and the miss of item i at Cache1 IS the read class of item i at Cache2.

Item popularity is expressed by the per-class request rates rather than by a
popularity distribution the cache draws from, which is also how the TTL
cache-network analysis is parameterised.

Note on the capacities: with LRU(1) at BOTH levels, Cache2 can never hit. A
miss at Cache1 inserts the item there, so the next Cache1 miss is necessarily
for a different item and consecutive arrivals at Cache2 are always distinct.
Cache2 therefore needs room for more than one item to show any hit at all.

Twin of `matlab/examples/basic/cacheModel/cachenet_tandem.m` and of
`cachenet_tandem` in `cpp/examples/basic/cacheModel.cpp`.
"""

from line_solver import *


def cachenet_tandem():
    """Build the two-level cache network and return (model, cache1, cache2)."""
    n = 3                          # items, shared by both caches
    p_access = [0.5, 0.3, 0.2]     # item popularity
    lambda_ = 1.0                  # request rate scale
    cap1, cap2 = 1, 2              # items held by Cache1, Cache2

    model = Network('CacheTandem')

    think = Delay(model, 'Think')
    cache1 = Cache(model, 'Cache1', n, cap1, ReplacementStrategy.LRU)
    cache2 = Cache(model, 'Cache2', n, cap2, ReplacementStrategy.LRU)

    # One dedicated class per item for each role.
    read, hit1, hit2, miss = [], [], [], []
    for f in range(n):
        read.append(ClosedClass(model, f'Read{f + 1}', 1, think, 0))
        hit1.append(ClosedClass(model, f'Hit1_{f + 1}', 0, think, 0))
        hit2.append(ClosedClass(model, f'Hit2_{f + 1}', 0, think, 0))
        miss.append(ClosedClass(model, f'Miss_{f + 1}', 0, think, 0))
        think.set_service(read[f], Exp(lambda_ * p_access[f]))

    cache1.set_item_read_classes(read, hit1)       # read[f] reads item f at Cache1
    cache1.set_miss_cache(read[0], cache2, hit2)   # Cache1 miss of item f = Cache2 read of item f
    cache2.set_item_miss_class(read[0], miss)      # Cache2 miss leaves the network

    P = model.init_routing_matrix()
    for f in range(n):
        P.set(read[f], read[f], think, cache1, 1.0)
        P.set(hit1[f], read[f], cache1, think, 1.0)   # hit at Cache1
        P.set(hit2[f], read[f], cache2, think, 1.0)   # hit at Cache2
        P.set(miss[f], read[f], cache2, think, 1.0)   # miss at both
    # Cache1 -> Cache2 for each item is registered by set_miss_cache.
    model.link(P)

    return model, cache1, cache2


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model, cache1, cache2 = cachenet_tandem()

    print(CTMC(model).getAvgCacheTable())

    model.reset()
    print(SSA(model, samples=50000, seed=1).getAvgCacheTable())

    model.reset()
    print(LDES(model, samples=200000, seed=1).getAvgCacheTable())
