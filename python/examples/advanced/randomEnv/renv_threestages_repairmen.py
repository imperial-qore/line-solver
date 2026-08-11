# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
def circul(c):
    """Returns a circulant matrix of order c (MATLAB-compatible).

    Generates a circulant matrix using the same algorithm as MATLAB's circul function.
    For scalar input c, generates a c×c circulant matrix with last element of the
    generating vector set to 1.

    Args:
        c: Either an integer (order) or a vector for the first row

    Returns:
        numpy.ndarray: Circulant matrix matching MATLAB behavior
    """
    if np.isscalar(c):
        if c == 1:
            return np.array([[1]])
        else:
            v = np.zeros(c)
            v[-1] = 1  # Last element = 1
            return circul(v)

    # c is a vector - use MATLAB's algorithm: C = sum_t c(t) * R^t
    # where R is the right circular shift matrix
    n = len(c)
    I = np.eye(n)
    R = I[:, list(range(1, n)) + [0]]  # MATLAB: I(:, [2:n, 1])

    C = np.zeros((n, n))
    Rt = I.copy()
    for t in range(n):
        C = C + c[t] * Rt
        Rt = Rt @ R

    return C

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
N = 2  # Job population  
M = 2  # Number of stations
E = 3  # Number of environment stages

# Create environment model
envModel = Environment('MyEnv', E)
envName = ['Stage1', 'Stage2', 'Stage3']
envType = ['UP', 'DOWN', 'FAST']

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

envSubModel = [qn1, qn2, qn3]

# Add stages to environment model
for e in range(E):
    envModel.add_stage(e, envName[e], envType[e], envSubModel[e])
# %%
# Define environment transition rates using circulant matrix
envRates = circul(3)  # Creates a 3x3 circulant matrix

print(f"Environment transition rates (circulant matrix):")
print(envRates)

# Add transitions with Erlang distributions
for e in range(E):
    for h in range(E):
        if envRates[e, h] > 0:
            mean_time = 1.0 / envRates[e, h]
            order = e + h + 2  # Erlang order based on stage indices (e+h+2 for 0-indexed to match MATLAB's e+h for 1-indexed)
            envModel.add_transition(e, h, Erlang.fit_mean_and_order(mean_time, order))
# %% [markdown]
# The metasolver considers an environment with 3 stages and a queueing network with 2 stations.
# This example illustrates the computation of the infinitesimal generator of the system.
# %%
# Create solvers for each submodel using CTMC
solvers = np.empty(E, dtype=object)
for e in range(E):
    solvers[e] = CTMC(envSubModel[e])

# Create environment solver
envSolver = ENV(envModel, solvers)

# Get results
# Get the generator (infinitesimal generator matrix)
infGen, stageInfGen = envSolver.generator()

# Print in sparse matrix format similar to MATLAB
from scipy import sparse
infGen_sparse = sparse.csr_matrix(infGen)
print(f"infGen =")
print(f"  {infGen.shape[0]}×{infGen.shape[1]} sparse double matrix ({infGen_sparse.nnz} nonzeros)")

# Print first few non-zero entries in MATLAB format
nnz_count = 0
for i in range(infGen.shape[0]):
    for j in range(infGen.shape[1]):
        if abs(infGen[i, j]) > 1e-10:
            print(f"   ({i+1},{j+1})    {infGen[i, j]:10.5g}")
            nnz_count += 1
            if nnz_count >= 30:  # Limit output
                print("   ...")
                break
    if nnz_count >= 30:
        break

print()
print("stageInfGen =")
print(f"  1×{len(stageInfGen)} cell array")
for i, gen in enumerate(stageInfGen):
    print(f"    {{{gen.shape[0]}×{gen.shape[1]} double}}", end="")
print()

# Get environment-averaged performance metrics
avg_table = envSolver.getAvgTable()
print(avg_table)