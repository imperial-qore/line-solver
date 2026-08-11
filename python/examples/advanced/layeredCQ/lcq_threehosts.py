# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.SILENT)
# %%
# Example of Layered Network with a multi-level cache
model = LayeredNetwork('LQNwithCaching')
nusers = 1
ntokens = 1

# Client
P1 = Processor(model, 'P1', 1, SchedStrategy.PS)
T1 = Task(model, 'T1', nusers, SchedStrategy.REF).on(P1)
E1 = Entry(model, 'E1').on(T1)

# Cache task
totalitems = 4
cachecapacity = [1, 1]
pAccess = DiscreteSampler([1.0/totalitems] * totalitems)
PC = Processor(model, 'Pc', 1, SchedStrategy.PS)
C2 = CacheTask(model, 'CT', totalitems, cachecapacity, ReplacementStrategy.RR, ntokens).on(PC)
I2 = ItemEntry(model, 'IE', totalitems, pAccess).on(C2)

# Server
P3 = Processor(model, 'P2', 1, SchedStrategy.PS)
T3 = Task(model, 'T2', 1, SchedStrategy.FCFS).on(P3)
E3 = Entry(model, 'E2').on(T3)
A3 = Activity(model, 'A2', Exp(5.0)).on(T3).bound_to(E3).replies_to(E3)

# Definition of activities
A1 = Activity(model, 'A1', Immediate()).on(T1).bound_to(E1).synch_call(I2, 1)
AC2 = Activity(model, 'Ac', Immediate()).on(C2).bound_to(I2)
AC2h = Activity(model, 'Ac_hit', Exp(1.0)).on(C2).replies_to(I2)
AC2m = Activity(model, 'Ac_miss', Exp(0.5)).on(C2).synch_call(E3, 1).replies_to(I2)

C2.add_precedence(ActivityPrecedence.cache_access(AC2, [AC2h, AC2m]))
# %%
lnoptions = LN.default_options()
lnoptions.verbose = 1
options = NC.default_options()
options.verbose = 0
solver1 = LN(model, lambda m: NC(m, options), lnoptions)
AvgTable1 = solver1.avg_table()
print(AvgTable1)
# %%
options2 = MVA.default_options()
options2.verbose = 0
solver2 = LN(model, lambda m: MVA(m, options2), lnoptions)
AvgTable2 = solver2.avg_table()
print(AvgTable2)