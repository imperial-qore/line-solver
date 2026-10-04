"""Tutorial 7: delayed-hit cache with a processor-sharing retrieval queue."""

import numpy as np

from line_solver import *


access_prob = list(np.array([49, 49, 49, 49, 7, 1, 1]) / 205)
model = Network("Retrieval")
source = Source(model, "Source")
cache = Cache(model, "Cache", len(access_prob), [6], ReplacementStrategy.RR)
queue = Queue(model, "Queue", SchedStrategy.PS)
sink = Sink(model, "Sink")

request = OpenClass(model, "InitClass", 0)
hit = OpenClass(model, "HitClass", 0)
miss = OpenClass(model, "MissClass", 0)
source.set_arrival(request, Exp(1))
cache.set_read(request, DiscreteSampler(access_prob))
cache.set_hit_class(request, hit)
cache.set_miss_class(request, miss)
queue.set_service(request, Exp(1))
cache.set_retrieval_system(request, miss, queue)

routing = model.init_routing_matrix()
routing.set(request, request, source, cache, 1.0)
routing.set(request, request, cache, queue, 1.0)
routing.set(request, request, queue, cache, 1.0)
routing.set(hit, hit, cache, sink, 1.0)
routing.set(miss, miss, cache, sink, 1.0)
model.link(routing)

NC(model).get_avg_table()
print("NC cache results")
print(cache.get_hit_ratio(), cache.get_miss_ratio())

MVA(model).get_avg_table()
print("MVA cache results")
print(cache.get_hit_ratio(), cache.get_miss_ratio(), cache.get_residt())
