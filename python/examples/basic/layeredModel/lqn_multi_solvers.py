# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.SILENT)

# Clear variables (equivalent to MATLAB's clear solver AvgTable)
solver = None
AvgTable = None
# %%
# MATLAB: model = LayeredNetwork('LQN1');
model = LayeredNetwork('LQN1')

# definition of processors, tasks and entries
# MATLAB: P1 = Processor(model, 'P1', Inf, SchedStrategy.INF);
# MATLAB: T1 = Task(model, 'T1', 1, SchedStrategy.REF).on(P1);
# MATLAB: E1 = Entry(model, 'E1').on(T1);
P1 = Processor(model, 'P1', float('inf'), SchedStrategy.INF)
T1 = Task(model, 'T1', 1, SchedStrategy.REF).on(P1)
E1 = Entry(model, 'E1').on(T1)

# MATLAB: P2 = Processor(model, 'P2', Inf, SchedStrategy.INF);
# MATLAB: T2 = Task(model, 'T2', Inf, SchedStrategy.INF).on(P2);
# MATLAB: E2 = Entry(model, 'E2').on(T2);
P2 = Processor(model, 'P2', float('inf'), SchedStrategy.INF)
T2 = Task(model, 'T2', float('inf'), SchedStrategy.INF).on(P2)
E2 = Entry(model, 'E2').on(T2)

# definition of activities
# MATLAB: T1.set_think_time(Erlang.fit_mean_and_order(0.0001,2));
T1.set_think_time(Erlang.fit_mean_and_order(0.0001, 2))

# MATLAB: A1 = Activity(model, 'A1', Exp(1.0)).on(T1).boundTo(E1).synchCall(E2,3);
# MATLAB: A2 = Activity(model, 'A2', APH.fit_mean_and_scv(1,10)).on(T2).boundTo(E2).repliesTo(E2);
A1 = Activity(model, 'A1', Exp(1.0)).on(T1).bound_to(E1).synch_call(E2, 3)
A2 = Activity(model, 'A2', APH.fit_mean_and_scv(1, 10)).on(T2).bound_to(E2).replies_to(E2)
# %%
# instantiate solvers
# MATLAB: options = LQNS.defaultOptions;
# MATLAB: options.keep = true;
# MATLAB: options.verbose = 1;
# MATLAB: %options.method = 'lqsim';
# MATLAB: %options.samples = 1e4;
# MATLAB: lqnssolver = LQNS(model, options);
# MATLAB: AvgTableLQNS = lqnssolver.get_avg_table;
# MATLAB: AvgTableLQNS
if LQNS.isAvailable():
    options = LQNS.default_options()
    options.keep = True
    options.verbose = 1
    # options.method = 'lqsim'
    # options.samples = int(1e4)
    lqnssolver = LQNS(model, options)
    AvgTableLQNS = lqnssolver.avg_table()
    print('AvgTableLQNS:')
    print(AvgTableLQNS)
else:
    print("LQNS solver not available - skipping")

# this method runs the MVA solver in each layer
# MATLAB: lnoptions = LN.defaultOptions;
# MATLAB: lnoptions.verbose = 0;
# MATLAB: lnoptions.seed = 2300;  
# MATLAB: options = MVA.defaultOptions;
# MATLAB: options.verbose = 0;
# MATLAB: solver{1} = LN(model, @(model) MVA(model, options), lnoptions);
# MATLAB: AvgTable{1} = solver{1}.get_avg_table
# MATLAB: AvgTable{1}
solver = {}
AvgTable = {}

lnoptions = LN.default_options()
lnoptions.verbose = 0
lnoptions.seed = 2300
options = MVA.default_options()
options.verbose = 0
solver[1] = LN(model, lambda model_arg: MVA(model_arg, options), lnoptions)
AvgTable[1] = solver[1].avg_table()
print('AvgTable[1]:')
print(AvgTable[1])

# this method runs the NC solver in each layer
# MATLAB: lnoptions = LN.defaultOptions;
# MATLAB: lnoptions.verbose = 0;
# MATLAB: lnoptions.seed = 2300;
# MATLAB: options = NC.defaultOptions;
# MATLAB: options.verbose = 0;
# MATLAB: solver{2} = LN(model, @(model) NC(model, options), lnoptions);
# MATLAB: AvgTable{2} = solver{2}.get_avg_table
# MATLAB: AvgTable{2}
lnoptions = LN.default_options()
lnoptions.verbose = 0
lnoptions.seed = 2300
options = NC.default_options()
options.verbose = 0
solver[2] = LN(model, lambda model_arg: NC(model_arg, options), lnoptions)
AvgTable[2] = solver[2].avg_table()
print('AvgTable[2]:')
print(AvgTable[2])

# this method adapts with the features of each layer
# MATLAB: %solver{2} = LN(model, @(model) LINE(model, LINE.defaultOptions), lnoptions);
# MATLAB: %AvgTable{2} = solver{2}.get_avg_table
# MATLAB: %AvgTable{2}
# solver[3] = LN(model, lambda model_arg: LINE(model_arg, LINE.default_options()), lnoptions)
# AvgTable[3] = solver[3].avg_table()
# print('AvgTable[3]:')
# print(AvgTable[3])