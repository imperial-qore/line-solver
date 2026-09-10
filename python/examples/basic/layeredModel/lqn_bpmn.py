# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.SILENT)
# %%
model = LayeredNetwork('myLayeredModel')

# Definition of processors - matching MATLAB exactly
P = []
P.append(Processor(model, 'R1_Processor', 100, SchedStrategy.FCFS))     # P{1}
P.append(Processor(model, 'R2_Processor', float('inf'), SchedStrategy.INF))  # P{2}
P.append(Processor(model, 'R3_Processor', 2, SchedStrategy.FCFS))       # P{3}
P.append(Processor(model, 'R1A_Processor', 7, SchedStrategy.FCFS))      # P{4}
P.append(Processor(model, 'R1B_Processor', 3, SchedStrategy.FCFS))      # P{5}
P.append(Processor(model, 'R2A_Processor', 4, SchedStrategy.FCFS))      # P{6}
P.append(Processor(model, 'R2B_Processor', 5, SchedStrategy.FCFS))      # P{7}

# Definition of tasks - matching MATLAB exactly
T = []
T.append(Task(model, 'R1_Task', 100, SchedStrategy.REF).on(P[0]).set_think_time(Exp.fit_mean(20)))  # T{1}
T.append(Task(model, 'R2_Task', float('inf'), SchedStrategy.INF).on(P[1]).set_think_time(Immediate()))  # T{2}
T.append(Task(model, 'R3_Task', 2, SchedStrategy.FCFS).on(P[2]).set_think_time(Immediate()))       # T{3}
T.append(Task(model, 'R1A_Task', 7, SchedStrategy.FCFS).on(P[3]).set_think_time(Immediate()))      # T{4}
T.append(Task(model, 'R1B_Task', 3, SchedStrategy.FCFS).on(P[4]).set_think_time(Immediate()))      # T{5}
T.append(Task(model, 'R2A_Task', 4, SchedStrategy.FCFS).on(P[5]).set_think_time(Immediate()))      # T{6}
T.append(Task(model, 'R2B_Task', 5, SchedStrategy.FCFS).on(P[6]).set_think_time(Immediate()))      # T{7}

# Definition of entries - matching MATLAB exactly
E = []
E.append(Entry(model, 'R1_Ref_Entry').on(T[0]))          # E{1}
E.append(Entry(model, 'R2_Synch_A2_Entry').on(T[1]))     # E{2}
E.append(Entry(model, 'R2_Synch_A5_Entry').on(T[1]))     # E{3}
E.append(Entry(model, 'R3_Synch_A9_Entry').on(T[2]))     # E{4}
E.append(Entry(model, 'R1A_Synch_A1_Entry').on(T[3]))    # E{5}
E.append(Entry(model, 'R1A_Synch_A2_Entry').on(T[3]))    # E{6}
E.append(Entry(model, 'R1A_Synch_A3_Entry').on(T[3]))    # E{7}
E.append(Entry(model, 'R1B_Synch_A4_Entry').on(T[4]))    # E{8}
E.append(Entry(model, 'R1B_Synch_A5_Entry').on(T[4]))    # E{9}
E.append(Entry(model, 'R1B_Synch_A6_Entry').on(T[4]))    # E{10}
E.append(Entry(model, 'R2A_Synch_A7_Entry').on(T[5]))    # E{11}
E.append(Entry(model, 'R2A_Synch_A8_Entry').on(T[5]))    # E{12}
E.append(Entry(model, 'R2A_Synch_A11_Entry').on(T[5]))   # E{13}
E.append(Entry(model, 'R2B_Synch_A9_Entry').on(T[6]))    # E{14}
E.append(Entry(model, 'R2B_Synch_A10_Entry').on(T[6]))   # E{15}
E.append(Entry(model, 'R2B_Synch_A12_Entry').on(T[6]))   # E{16}
# %%
# Definition of activities - matching MATLAB exactly (all 29 activities)
A = []

# T1 activities (A{1}-A{6} in MATLAB)
A.append(Activity(model, 'A1_Empty', Immediate()).on(T[0]).bound_to(E[0]).synch_call(E[4], 1))  # A{1}
A.append(Activity(model, 'A2_Empty', Immediate()).on(T[0]).synch_call(E[5], 1))  # A{2}
A.append(Activity(model, 'A5_Empty', Immediate()).on(T[0]).synch_call(E[8], 1))  # A{3}
A.append(Activity(model, 'A6_Empty', Immediate()).on(T[0]).synch_call(E[9], 1))  # A{4}
A.append(Activity(model, 'A3_Empty', Immediate()).on(T[0]).synch_call(E[6], 1))  # A{5}
A.append(Activity(model, 'A4_Empty', Immediate()).on(T[0]).synch_call(E[7], 1))  # A{6}

# T2 activities (A{7}-A{13} in MATLAB)
A.append(Activity(model, 'E4_Empty', Immediate()).on(T[1]).bound_to(E[1]))  # A{7}
A.append(Activity(model, 'A7_Empty', Immediate()).on(T[1]).synch_call(E[10], 1))  # A{8}
A.append(Activity(model, 'A8_Empty', Immediate()).on(T[1]).synch_call(E[11], 1))  # A{9}
A.append(Activity(model, 'A9_Empty', Immediate()).on(T[1]).synch_call(E[13], 1))  # A{10}
A.append(Activity(model, 'A11_Empty', Immediate()).on(T[1]).synch_call(E[12], 1).replies_to(E[1]))  # A{11}
A.append(Activity(model, 'A12_Empty', Immediate()).on(T[1]).bound_to(E[2]).synch_call(E[15], 1).replies_to(E[2]))  # A{12}
A.append(Activity(model, 'A10_Empty', Immediate()).on(T[1]).synch_call(E[14], 1))  # A{13}

# T3 activity (A{14} in MATLAB)
A.append(Activity(model, 'A13', Exp.fit_mean(10)).on(T[2]).bound_to(E[3]).replies_to(E[3]))  # A{14}

# T4 activities (A{15}-A{18} in MATLAB)
A.append(Activity(model, 'A1', Exp.fit_mean(7)).on(T[3]).bound_to(E[4]).replies_to(E[4]))  # A{15}
A.append(Activity(model, 'A2', Exp.fit_mean(4)).on(T[3]).bound_to(E[5]))  # A{16}
A.append(Activity(model, 'A3', Exp.fit_mean(5)).on(T[3]).bound_to(E[6]).replies_to(E[6]))  # A{17}
A.append(Activity(model, 'A2_Res_Empty', Immediate()).on(T[3]).synch_call(E[1], 1).replies_to(E[5]))  # A{18}

# T5 activities (A{19}-A{22} in MATLAB)
A.append(Activity(model, 'A4', Exp.fit_mean(8)).on(T[4]).bound_to(E[7]).replies_to(E[7]))  # A{19}
A.append(Activity(model, 'A5', Exp.fit_mean(4)).on(T[4]).bound_to(E[8]))  # A{20}
A.append(Activity(model, 'A6', Exp.fit_mean(6)).on(T[4]).bound_to(E[9]).replies_to(E[9]))  # A{21}
A.append(Activity(model, 'A5_Res_Empty', Immediate()).on(T[4]).synch_call(E[2], 1).replies_to(E[8]))  # A{22}

# T6 activities (A{23}-A{25} in MATLAB)
A.append(Activity(model, 'A7', Exp.fit_mean(6)).on(T[5]).bound_to(E[10]).replies_to(E[10]))  # A{23}
A.append(Activity(model, 'A8', Exp.fit_mean(8)).on(T[5]).bound_to(E[11]).replies_to(E[11]))  # A{24}
A.append(Activity(model, 'A11', Exp.fit_mean(4)).on(T[5]).bound_to(E[12]).replies_to(E[12]))  # A{25}

# T7 activities (A{26}-A{29} in MATLAB)
A.append(Activity(model, 'A9', Exp.fit_mean(4)).on(T[6]).bound_to(E[13]))  # A{26}
A.append(Activity(model, 'A10', Exp.fit_mean(6)).on(T[6]).bound_to(E[14]).replies_to(E[14]))  # A{27}
A.append(Activity(model, 'A12', Exp.fit_mean(8)).on(T[6]).bound_to(E[15]).replies_to(E[15]))  # A{28}
A.append(Activity(model, 'A9_Res_Empty', Immediate()).on(T[6]).synch_call(E[3], 1).replies_to(E[13]))  # A{29}

# Add precedences - matching MATLAB lines 70-81
T[0].add_precedence(ActivityPrecedence.Serial(A[0], A[1]))  # Line 70: Serial(A{1}, A{2})
T[0].add_precedence(ActivityPrecedence.Serial(A[2], A[3]))  # Line 71: Serial(A{3}, A{4})
T[1].add_precedence(ActivityPrecedence.Serial(A[6], A[7]))  # Line 72: Serial(A{7}, A{8})
T[1].add_precedence(ActivityPrecedence.Serial(A[9], A[12]))  # Line 73: Serial(A{10}, A{13})
T[3].add_precedence(ActivityPrecedence.Serial(A[15], A[17]))  # Line 74: Serial(A{16}, A{18})
T[4].add_precedence(ActivityPrecedence.Serial(A[19], A[21]))  # Line 75: Serial(A{20}, A{22})
T[6].add_precedence(ActivityPrecedence.Serial(A[25], A[28]))  # Line 76: Serial(A{26}, A{29})
T[0].add_precedence(ActivityPrecedence.OrFork(A[1], [A[4], A[5]], [0.6, 0.4]))  # Line 77: OrFork(A{2}, {A{5}, A{6}}, [0.6,0.4])
T[1].add_precedence(ActivityPrecedence.AndFork(A[7], [A[8], A[9]]))  # Line 78: AndFork(A{8}, {A{9}, A{10}})
T[0].add_precedence(ActivityPrecedence.OrJoin([A[4], A[5]], A[2]))  # Line 79: OrJoin({A{5}, A{6}}, A{3})
T[1].add_precedence(ActivityPrecedence.AndJoin([A[8], A[12]], A[10]))  # Line 80: AndJoin({A{9}, A{13}}, A{11})
# %%
print('This example illustrates the solution of a complex layered queueing network extracted from a BPMN model.')

# Solve with LQNS solver
options = LQNS.default_options()
options.keep = True  # uncomment to keep the intermediate XML files generated while translating the model to LQNS

if LQNS.isAvailable():
    solver = LQNS(model)
    AvgTableLQNS = solver.avg_table()
    print('\nLQNS Results:')
    print(AvgTableLQNS)

# LN(MVA) solver
lnoptions = LN.default_options()
lnoptions.verbose = 0
options = MVA.default_options()
options.verbose = 0
solverLN = LN(model, lambda model_arg: MVA(model_arg, options), lnoptions)
AvgTableLN = solverLN.avg_table()
print('\nLN(MVA) Results:')
print(AvgTableLN)