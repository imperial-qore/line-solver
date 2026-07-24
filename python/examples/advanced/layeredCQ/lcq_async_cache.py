"""
Layered cache queueing model with asynchronous (non-blocking) cache access.

This example demonstrates async cache calls where:
- Client makes async call to cache (fire-and-forget, non-blocking)
- Client continues immediately without waiting for cache response
- Cache still uses POST_CACHE precedence for hit/miss determination
- Use for prefetching, cache warming, or non-critical operations

Compare with lcq_singlehost.ipynb which uses synchronous cache calls.
"""

from line_solver import *
import numpy as np

# Set verbose level
GlobalConstants.set_verbose(VerboseLevel.SILENT)

# Create layered network model
model = LayeredNetwork('AsyncCacheLQN')

# Client processor and task
P1 = Processor(model, 'P1', 1, SchedStrategy.PS)
T1 = Task(model, 'T1', 1, SchedStrategy.REF).on(P1)
E1 = Entry(model, 'E1').on(T1)

# Cache processor and task
totalitems = 4
cachecapacity = 2

# Create uniform access probabilities using Python list
access_probs = [1.0 / totalitems] * totalitems
pAccess = DiscreteSampler(access_probs)

PC = Processor(model, 'PC', 1, SchedStrategy.PS)
C2 = CacheTask(model, 'C2', totalitems, cachecapacity, ReplacementStrategy.LRU, 1).on(PC)
I2 = ItemEntry(model, 'I2', totalitems, pAccess).on(C2)

# Definition of activities
# Client activity with ASYNC call to cache (non-blocking, fire-and-forget)
A1 = Activity(model, 'A1', Immediate()).on(T1).bound_to(E1).asynch_call(I2, 1)

# Cache activities (unchanged from sync version)
AC2 = Activity(model, 'AC2', Immediate()).on(C2).bound_to(I2)
AC2h = Activity(model, 'AC2h', Exp(1.0)).on(C2).replies_to(I2)   # Cache hit
AC2m = Activity(model, 'AC2m', Exp(0.5)).on(C2).replies_to(I2)   # Cache miss

# Add cache access precedence (unchanged)
C2.add_precedence(ActivityPrecedence.cache_access(AC2, [AC2h, AC2m]))

# Solve the model
print("=== Solving Async Cache Model ===")
lnoptions = LN.default_options()
lnoptions.verbose = True
options = MVA.default_options()
options.verbose = False

solver = LN(model, lambda model: MVA(model, options), lnoptions)
AvgTable = solver.avg_table()

print("\n=== Average Performance Metrics ===")
print(AvgTable)

print("\n=== Compare with Synchronous Version ===")
print("To compare async vs sync cache access:")
print("1. Run this script (async): python lcq_async_cache.py")
print("2. Run sync notebook: lcq_singlehost.ipynb")
print("")
print("Expected differences:")
print("- Async: Lower client response time (no blocking)")
print("- Async: Higher client throughput (no cache bottleneck)")
print("- Similar cache hit/miss ratios (same cache logic)")
