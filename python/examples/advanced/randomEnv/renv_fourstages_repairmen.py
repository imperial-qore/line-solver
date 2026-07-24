# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
def renv_genqn(rate, N):
    """Helper function to generate a queueing network for random environment advanced."""
    qnet = Network('qn1')
    
    node = np.empty(2, dtype=object)
    node[0] = Delay(qnet, 'Queue1')
    node[1] = Queue(qnet, 'Queue2', SchedStrategy.PS)
    
    jobclass = np.empty(1, dtype=object)
    jobclass[0] = ClosedClass(qnet, 'Class1', N, node[0], 0)
    
    node[0].set_service(jobclass[0], Exp(rate[0]))
    node[1].set_service(jobclass[0], Exp(rate[1]))
    
    P = qnet.init_routing_matrix()
    P.set(jobclass[0], jobclass[0], [[0, 1], [1, 0]])
    qnet.link(P)
    
    return qnet
# %%
# Model parameters
N = 30  # Job population
M = 3   # Number of stations
E = 4   # Number of environment stages

# Create environment model
envModel = Environment('MyEnv', E)
envName = ['Stage1', 'Stage2', 'Stage3', 'Stage4']
envType = ['UP', 'DOWN', 'FAST', 'SLOW']

# Create rate matrix
rate = np.ones((M, E))
rate[M-1, :] = np.arange(1, E+1)  # rate(M,1:E)=(1:E)
rate[0, :] = np.arange(E, 0, -1)  # rate(1,1:E)=(E:-1:1)

print(f"Rate matrix:")
print(rate)
# %%
# Create queueing networks for each environment stage
qn1 = renv_genqn(rate[:, 0], N)
qn2 = renv_genqn(rate[:, 1], N)
qn3 = renv_genqn(rate[:, 2], N)
qn4 = renv_genqn(rate[:, 3], N)
envSubModel = [qn1, qn2, qn3, qn4]

# Add stages to environment model
for e in range(E):
    envModel.add_stage(e, envName[e], envType[e], envSubModel[e])
# %%
# Define environment transition rates
envRates = np.array([[0, 1, 0, 0],
                     [0, 0, 1, 1],
                     [1, 0, 0, 1],
                     [1, 1, 0, 0]]) / 2

print(f"Environment transition rates:")
print(envRates)

# Add transitions with APH distributions
for e in range(E):
    for h in range(E):
        if envRates[e, h] > 0:
            mean_time = 1.0 / envRates[e, h]
            envModel.add_transition(e, h, APH.fit_mean_and_scv(mean_time, 0.5))
# %% [markdown]
# The metasolver considers an environment with 4 stages and a queueing network with 3 stations.
# Every time the stage changes, the queueing network will modify the service rates of the stations.
# %%
print("The metasolver considers an environment with 4 stages and a queueing network with 3 stations.")
print("Every time the stage changes, the queueing network will modify the service rates of the stations.")

envModel.getStageTable()

# Create environment solver using lambda factory pattern
envSolver = ENV(envModel, lambda m: FLD(m))

# Get results
QN, UN, TN = envSolver.getAvg()
AvgTable = envSolver.getAvgTable()
print("AvgTable =")
print(AvgTable)