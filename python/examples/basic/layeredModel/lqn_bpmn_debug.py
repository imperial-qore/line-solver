# Debug version of lqn_bpmn to trace LN solver iteration
from line_solver import *
import numpy as np

GlobalConstants.set_verbose(VerboseLevel.STD)

print('DEBUG: LN(MVA) iteration tracing for lqn_bpmn model\n')

model = LayeredNetwork('myLayeredModel')

# Definition of processors
P = []
P.append(Processor(model, 'R1_Processor', 100, SchedStrategy.FCFS))
P.append(Processor(model, 'R2_Processor', GlobalConstants.MaxInt, SchedStrategy.INF))
P.append(Processor(model, 'R3_Processor', 2, SchedStrategy.FCFS))
P.append(Processor(model, 'R1A_Processor', 7, SchedStrategy.FCFS))
P.append(Processor(model, 'R1B_Processor', 3, SchedStrategy.FCFS))
P.append(Processor(model, 'R2A_Processor', 4, SchedStrategy.FCFS))
P.append(Processor(model, 'R2B_Processor', 5, SchedStrategy.FCFS))

# Definition of tasks
T = []
T.append(Task(model, 'R1_Task', 100, SchedStrategy.REF).on(P[0]).set_think_time(Exp.fit_mean(20)))
T.append(Task(model, 'R2_Task', GlobalConstants.MaxInt, SchedStrategy.INF).on(P[1]).set_think_time(Immediate()))
T.append(Task(model, 'R3_Task', 2, SchedStrategy.FCFS).on(P[2]).set_think_time(Immediate()))
T.append(Task(model, 'R1A_Task', 7, SchedStrategy.FCFS).on(P[3]).set_think_time(Immediate()))
T.append(Task(model, 'R1B_Task', 3, SchedStrategy.FCFS).on(P[4]).set_think_time(Immediate()))
T.append(Task(model, 'R2A_Task', 4, SchedStrategy.FCFS).on(P[5]).set_think_time(Immediate()))
T.append(Task(model, 'R2B_Task', 5, SchedStrategy.FCFS).on(P[6]).set_think_time(Immediate()))

# Definition of entries
E = []
E.append(Entry(model, 'R1_Ref_Entry').on(T[0]))
E.append(Entry(model, 'R2_Synch_A2_Entry').on(T[1]))
E.append(Entry(model, 'R2_Synch_A5_Entry').on(T[1]))
E.append(Entry(model, 'R3_Synch_A9_Entry').on(T[2]))
E.append(Entry(model, 'R1A_Synch_A1_Entry').on(T[3]))
E.append(Entry(model, 'R1A_Synch_A2_Entry').on(T[3]))
E.append(Entry(model, 'R1A_Synch_A3_Entry').on(T[3]))
E.append(Entry(model, 'R1B_Synch_A4_Entry').on(T[4]))
E.append(Entry(model, 'R1B_Synch_A5_Entry').on(T[4]))
E.append(Entry(model, 'R1B_Synch_A6_Entry').on(T[4]))
E.append(Entry(model, 'R2A_Synch_A7_Entry').on(T[5]))
E.append(Entry(model, 'R2A_Synch_A8_Entry').on(T[5]))
E.append(Entry(model, 'R2A_Synch_A11_Entry').on(T[5]))
E.append(Entry(model, 'R2B_Synch_A9_Entry').on(T[6]))
E.append(Entry(model, 'R2B_Synch_A10_Entry').on(T[6]))
E.append(Entry(model, 'R2B_Synch_A12_Entry').on(T[6]))

# Definition of activities
A = []
A.append(Activity(model, 'A1_Empty', Immediate()).on(T[0]).bound_to(E[0]).synch_call(E[4], 1))
A.append(Activity(model, 'A2_Empty', Immediate()).on(T[0]).synch_call(E[5], 1))
A.append(Activity(model, 'A5_Empty', Immediate()).on(T[0]).synch_call(E[8], 1))
A.append(Activity(model, 'A6_Empty', Immediate()).on(T[0]).synch_call(E[9], 1))
A.append(Activity(model, 'A3_Empty', Immediate()).on(T[0]).synch_call(E[6], 1))
A.append(Activity(model, 'A4_Empty', Immediate()).on(T[0]).synch_call(E[7], 1))
A.append(Activity(model, 'E4_Empty', Immediate()).on(T[1]).bound_to(E[1]))
A.append(Activity(model, 'A7_Empty', Immediate()).on(T[1]).synch_call(E[10], 1))
A.append(Activity(model, 'A8_Empty', Immediate()).on(T[1]).synch_call(E[11], 1))
A.append(Activity(model, 'A9_Empty', Immediate()).on(T[1]).synch_call(E[13], 1))
A.append(Activity(model, 'A11_Empty', Immediate()).on(T[1]).synch_call(E[12], 1).replies_to(E[1]))
A.append(Activity(model, 'A12_Empty', Immediate()).on(T[1]).bound_to(E[2]).synch_call(E[15], 1).replies_to(E[2]))
A.append(Activity(model, 'A10_Empty', Immediate()).on(T[1]).synch_call(E[14], 1))
A.append(Activity(model, 'A13', Exp.fit_mean(10)).on(T[2]).bound_to(E[3]).replies_to(E[3]))
A.append(Activity(model, 'A1', Exp.fit_mean(7)).on(T[3]).bound_to(E[4]).replies_to(E[4]))
A.append(Activity(model, 'A2', Exp.fit_mean(4)).on(T[3]).bound_to(E[5]))
A.append(Activity(model, 'A3', Exp.fit_mean(5)).on(T[3]).bound_to(E[6]).replies_to(E[6]))
A.append(Activity(model, 'A2_Res_Empty', Immediate()).on(T[3]).synch_call(E[1], 1).replies_to(E[5]))
A.append(Activity(model, 'A4', Exp.fit_mean(8)).on(T[4]).bound_to(E[7]).replies_to(E[7]))
A.append(Activity(model, 'A5', Exp.fit_mean(4)).on(T[4]).bound_to(E[8]))
A.append(Activity(model, 'A6', Exp.fit_mean(6)).on(T[4]).bound_to(E[9]).replies_to(E[9]))
A.append(Activity(model, 'A5_Res_Empty', Immediate()).on(T[4]).synch_call(E[2], 1).replies_to(E[8]))
A.append(Activity(model, 'A7', Exp.fit_mean(6)).on(T[5]).bound_to(E[10]).replies_to(E[10]))
A.append(Activity(model, 'A8', Exp.fit_mean(8)).on(T[5]).bound_to(E[11]).replies_to(E[11]))
A.append(Activity(model, 'A11', Exp.fit_mean(4)).on(T[5]).bound_to(E[12]).replies_to(E[12]))
A.append(Activity(model, 'A9', Exp.fit_mean(4)).on(T[6]).bound_to(E[13]))
A.append(Activity(model, 'A10', Exp.fit_mean(6)).on(T[6]).bound_to(E[14]).replies_to(E[14]))
A.append(Activity(model, 'A12', Exp.fit_mean(8)).on(T[6]).bound_to(E[15]).replies_to(E[15]))
A.append(Activity(model, 'A9_Res_Empty', Immediate()).on(T[6]).synch_call(E[3], 1).replies_to(E[13]))

# Add precedences
T[0].add_precedence(ActivityPrecedence.Serial(A[0], A[1]))
T[0].add_precedence(ActivityPrecedence.Serial(A[2], A[3]))
T[1].add_precedence(ActivityPrecedence.Serial(A[6], A[7]))
T[1].add_precedence(ActivityPrecedence.Serial(A[9], A[12]))
T[3].add_precedence(ActivityPrecedence.Serial(A[15], A[17]))
T[4].add_precedence(ActivityPrecedence.Serial(A[19], A[21]))
T[6].add_precedence(ActivityPrecedence.Serial(A[25], A[28]))
T[0].add_precedence(ActivityPrecedence.OrFork(A[1], [A[4], A[5]], [0.6, 0.4]))
T[1].add_precedence(ActivityPrecedence.AndFork(A[7], [A[8], A[9]]))
T[0].add_precedence(ActivityPrecedence.OrJoin([A[4], A[5]], A[2]))
T[1].add_precedence(ActivityPrecedence.AndJoin([A[8], A[12]], A[10]))

# LN(MVA) solver with verbose output
print('\n=== Starting LN(MVA) Solver ===')

lnoptions = LN.default_options()
lnoptions.verbose = 1  # Enable verbose
lnoptions.iter_max = 50  # Limit iterations for debugging

options = MVA.default_options()
options.verbose = 0

solverLN = LN(model, lambda model_arg: MVA(model_arg, options), lnoptions)

# Print layer info
print(f'\nNumber of layers: {solverLN.nlayers}')

# Patch the iterate method to add debug output
original_analyze = solverLN.analyze

def debug_analyze(it, e):
    result, time = original_analyze(it, e)
    if result is not None:
        TN = result.get('TN', np.array([]))
        QN = result.get('QN', np.array([]))
        UN = result.get('UN', np.array([]))
        max_TN = np.nanmax(TN) if TN.size > 0 else 0
        max_QN = np.nanmax(QN) if QN.size > 0 else 0
        max_UN = np.nanmax(UN) if UN.size > 0 else 0
        print(f'  Iter {it}, Layer {e}: max(TN)={max_TN:.6f}, max(QN)={max_QN:.6f}, max(UN)={max_UN:.6f}')
    return result, time

solverLN.analyze = debug_analyze

# Run solver
AvgTable = solverLN.avg_table()

# Print iteration results summary
print('\n=== Iteration Summary ===')
if solverLN.results:
    numIters = len(solverLN.results)
    print(f'Total iterations: {numIters}')

    # Show last iteration details
    if numIters > 0:
        print(f'\n--- Last Iteration ({numIters}) ---')
        last_results = solverLN.results[-1]
        for e, r in enumerate(last_results):
            if r is not None:
                TN = r.get('TN', np.array([]))
                QN = r.get('QN', np.array([]))
                print(f'Layer {e}:')
                print(f'  TN: {TN.flatten()[:5]}...')  # First 5 values
                print(f'  QN: {QN.flatten()[:5]}...')

print('\n=== Final Results (Key Tasks) ===')
taskNames = ['R1_Task', 'R2_Task', 'R3_Task', 'R1A_Task', 'R1B_Task', 'R2A_Task', 'R2B_Task']
for name in taskNames:
    mask = AvgTable['Node'] == name
    if mask.any():
        row = AvgTable[mask].iloc[0]
        print(f"{name}: Tput={row['Tput']:.6f}, QLen={row['QLen']:.6f}, Util={row['Util']:.6f}")
