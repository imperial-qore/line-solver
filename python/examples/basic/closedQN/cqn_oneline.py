# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Service times - S(i,r) - mean service time of class r at station i
D = np.array([[10, 5], [5, 9]])
# Number of jobs - N(r) - number of jobs of class r
N = np.array([1, 2])
# Think times - Z(r) - mean service time of class r at delay station i
Z = np.array([[91, 92], [93, 94]])
# %%
# Create cyclic model using the Network.cyclicPsInf method
model = Network.cyclicPsInf(N, D, Z)
# %%
# Solve using exact MVA
avgTable = MVA(model, method='exact').avg_table()
print(avgTable)