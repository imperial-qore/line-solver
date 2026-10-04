# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.SILENT)
# %%
print('This example illustrates a layered network with a loop.')

model = LayeredNetwork('myLayeredModel')

# Layer 1 - Reference task
P1 = Processor(model, 'P1', GlobalConstants.MaxInt, SchedStrategy.INF)
T1 = Task(model, 'T1', 1, SchedStrategy.REF).on(P1)
T1.set_think_time(Immediate())
E1 = Entry(model, 'Entry').on(T1)

# Layer 2 - Service task
P2 = Processor(model, 'P2', GlobalConstants.MaxInt, SchedStrategy.INF)
T2 = Task(model, 'T2', 1, SchedStrategy.INF).on(P2).set_think_time(Immediate())
E2 = Entry(model, 'E2').on(T2)

# Layer 3 - Backend task
P3 = Processor(model, 'P3', 5, SchedStrategy.PS)
T3 = Task(model, 'T3', 20, SchedStrategy.INF).on(P3)
T3.set_think_time(Exp.fit_mean(10))
E3 = Entry(model, 'E1').on(T3)  # Note: Entry named 'E1' in original

# Activities for T1 (with loop)
A1 = Activity(model, 'A1', Exp.fit_mean(1)).on(T1).bound_to(E1)
A2 = Activity(model, 'A2', Exp.fit_mean(2)).on(T1)
A3 = Activity(model, 'A3', Exp.fit_mean(3)).on(T1).synch_call(E2, 1)

# Activities for T2 (with fork-join)
B1 = Activity(model, 'B1', Exp.fit_mean(0.1)).on(T2).bound_to(E2)
B2 = Activity(model, 'B2', Exp.fit_mean(0.2)).on(T2)
B3 = Activity(model, 'B3', Exp.fit_mean(0.3)).on(T2)
B4 = Activity(model, 'B4', Exp.fit_mean(0.4)).on(T2)
B5 = Activity(model, 'B5', Exp.fit_mean(0.5)).on(T2)
B6 = Activity(model, 'B6', Exp.fit_mean(0.6)).on(T2).synch_call(E3, 1).replies_to(E2)

# Activities for T3 (with or-fork/join)
C1 = Activity(model, 'C1', Exp.fit_mean(0.1)).on(T3).bound_to(E3)
C2 = Activity(model, 'C2', Exp.fit_mean(0.2)).on(T3)
C3 = Activity(model, 'C3', Exp.fit_mean(0.3)).on(T3)
C4 = Activity(model, 'C4', Exp.fit_mean(0.4)).on(T3)
C5 = Activity(model, 'C5', Exp.fit_mean(0.5)).on(T3).replies_to(E3)
# %%
# Add precedence relationships

# T1: Loop with 3 iterations
T1.add_precedence(ActivityPrecedence.Loop(A1, [A2, A3], 3))

# T2: Serial connection and AND fork-join
T2.add_precedence(ActivityPrecedence.Serial(B4, B5))
T2.add_precedence(ActivityPrecedence.AndFork(B1, [B2, B3, B4]))
T2.add_precedence(ActivityPrecedence.AndJoin([B2, B3, B5], B6))

# T3: OR fork-join with probabilities
T3.add_precedence(ActivityPrecedence.OrFork(C1, [C2, C3, C4], [0.3, 0.3, 0.4]))
T3.add_precedence(ActivityPrecedence.OrJoin([C2, C3, C4], C5))
# %%
# Solve with different solvers

if LQNS.isAvailable():
    solver_lqns = LQNS(model)
    avg_table_lqns = solver_lqns.avg_table()
    print('\nLQNS Results:')
    print(avg_table_lqns)
# %%
solver_ln = LN(model)
avg_table_ln = solver_ln.avg_table()
print('\nLN Results:')
print(avg_table_ln)