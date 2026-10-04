# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
N = 1
M = 2
E = 2
envModel = Environment('MyEnv', E)
envName = ['Stage1', 'Stage2']
envType = ['UP', 'DOWN']

rate = np.array([[2,1],[1,2]])
# %%
def renv_genqn(rate, N):
    # sn1
    qnet = Network('qn1')
    
    node = np.empty(2, dtype=object)
    node[0] = Delay(qnet, 'Queue1')
    node[1] = Queue(qnet, 'Queue2', SchedStrategy.PS)

    jobclass = np.empty(1, dtype=object)
    jobclass[0] = ClosedClass(qnet, 'Class1', N, node[0], 0)
    
    node[0].set_service(jobclass[0], Exp(rate[0]))
    node[1].set_service(jobclass[0], Exp(rate[1]))
    
    P = qnet.init_routing_matrix()
    P.set(jobclass[0],jobclass[0], [[0,1],[1,0]])
    qnet.link(P)
    return qnet
# %%
envSubModel = [renv_genqn(rate[:,0],N), renv_genqn(rate[:,1],N)]

for e in range(E):
    envModel.add_stage(e, envName[e], envType[e], envSubModel[e])
 
envRates = [[0,1], [0.5,0.5]]
for e in range(E):
    for h in range(E):
        if envRates[e][h]>0.0:
            envModel.add_transition(e, h, Exp(envRates[e][h]))
# %%
# Display stage table
print("Stage Table:")
envModel.stage_table()
# %%
# Configure solver options  
from line_solver import SolverType
options = SolverOptions(SolverType.ENV.value)
options.timespan = [0, float('inf')]
options.iter_max = 100
options.iter_tol = 0.01
options.method = 'default'
options.verbose = True

# Configure fluid solver options
sfoptions = SolverOptions(SolverType.FLUID.value)
sfoptions.timespan = [0, 1e3]
sfoptions.verbose = False

# Create solvers for each submodel
solvers = np.empty(E, dtype=object)
for e in range(E):
    solvers[e] = FLD(envSubModel[e], sfoptions)

# Create environment solver
envSolver = ENV(envModel, solvers, options)

# Get results
try:
    QN, UN, TN = envSolver.avg()
    print("Average queue lengths (QN):", QN)
    print("Average utilizations (UN):", UN)
    print("Average throughputs (TN):", TN)
    
    AvgTable = envSolver.avg_table()
    print("\nAverage Table:")
    print(AvgTable)
    
except Exception as e:
    print(f"Error during solving: {e}")
    print("Note: Environment solver features may not be fully implemented")
# %% [markdown]
# ## Alternative: CTMC Solver (Commented)
# The MATLAB version also shows an alternative using CTMC, which is commented out in the original.
# %%
# Alternative solver using CTMC (commented out as in MATLAB version)
# scoptions = SolverOptions()
# scoptions.timespan = [0, 1e3]
# scoptions.verbose = False
# envSolver = ENV(envModel, lambda model: CTMC(model, scoptions), options)
# QNc, UNc, TNc = envSolver.avg()
# AvgTableC = envSolver.avg_table()
# print("CTMC Results:")
# print("Average Table (CTMC):", AvgTableC)