# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.SILENT)
# %%
print('Example of layered model with a task paying setup and delay-off times')

model = LayeredNetwork('setup_example')

# definition of processors, tasks and entries
P1 = Processor(model, 'P1', GlobalConstants.MaxInt, SchedStrategy.INF)
T1 = Task(model, 'T1', 1, SchedStrategy.REF).on(P1)
E1 = Entry(model, 'E1').on(T1)

P2 = Processor(model, 'P2', 4, SchedStrategy.FCFS)
# Setup task: capacity, scheduling strategy, think time, setup time, delay-off time
T2 = SetupTask(model, 'F2', 6, SchedStrategy.FCFS).on(P2).set_think_time(Exp.fit_mean(8.0))
T2.setSetupTime(Exp(1.0))      # Activation time
T2.setDelayOffTime(Exp(2.0))   # Idle time before a server powers off

E2 = Entry(model, 'E2').on(T2)

# definition of activities
A1 = Activity(model, 'A1', Exp(1.0)).on(T1).bound_to(E1).synch_call(E2, 1)
A2 = Activity(model, 'A2', Exp(3.0)).on(T2).bound_to(E2).replies_to(E2)
# %%
# Define a helper function for solver selection (as in original MATLAB)
def myFunction(model):
    """Helper function to select appropriate solver for the model"""
    # This could implement logic to choose between MVA, MAM, etc.
    # For now, use MVA as default
    options = MVA.default_options()
    options.verbose = 0
    return MVA(model, options)

# Solve with LN solver using custom function
lnoptions = LN.default_options()
lnoptions.verbose = 0
lnoptions.seed = 23000

solver = LN(model, myFunction, lnoptions)
avg_table = solver.avg_table()
print('\nLN(MVA) Results:')
print(avg_table)