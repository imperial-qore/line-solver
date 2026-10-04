# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.SILENT)
# %%
model = LayeredNetwork('cacheInLayeredNetwork')

# Client
P1 = Processor(model, 'P1', 1, SchedStrategy.PS)
T1 = Task(model, 'T1', 1, SchedStrategy.REF).on(P1)
E1 = Entry(model, 'E1').on(T1)

# Cache task
totalitems = 4
cachecapacity = 2
# Create uniform access probabilities using Python list
access_probs = [1.0 / totalitems] * totalitems
pAccess = DiscreteSampler(access_probs)
PC = Processor(model, 'PC', 1, SchedStrategy.PS)
C2 = CacheTask(model, 'C2', totalitems, cachecapacity, ReplacementStrategy.RR, 1).on(PC)
I2 = ItemEntry(model, 'I2', totalitems, pAccess).on(C2)

# Definition of activities
A1 = Activity(model, 'A1', Immediate()).on(T1).bound_to(E1).synch_call(I2, 1)
AC2 = Activity(model, 'AC2', Immediate()).on(C2).bound_to(I2)
AC2h = Activity(model, 'AC2h', Exp(1.0)).on(C2).replies_to(I2)  # Cache hit
AC2m = Activity(model, 'AC2m', Exp(0.5)).on(C2).replies_to(I2)  # Cache miss

# Add cache access precedence
C2.add_precedence(ActivityPrecedence.cache_access(AC2, [AC2h, AC2m]))
# %%
lnoptions = LN.default_options()
# lnoptions.iter_max = 1  # Note: iter_max option may not be available
lnoptions.verbose = True
options = MVA.default_options()
options.verbose = False

solver = LN(model, lambda model: MVA(model, options), lnoptions)
AvgTable = solver.avg_table()
print(AvgTable)