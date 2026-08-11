# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Create network
model = Network('model')

# Block 1: nodes
node1 = Delay(model, 'Delay')
node2 = Queue(model, 'Queue1', SchedStrategy.PS)
node3 = Queue(model, 'Queue2', SchedStrategy.PS)
node3.setNumberOfServers(2)  # Queue2 has 2 servers
# %%
# Block 2: classes
N = [1, 0, 4, 0]  # Population for each class
jobclass1 = ClosedClass(model, 'Class1', N[0], node1, 0)
jobclass2 = ClosedClass(model, 'Class2', N[1], node1, 0)
jobclass3 = ClosedClass(model, 'Class3', N[2], node1, 0)
jobclass4 = ClosedClass(model, 'Class4', N[3], node1, 0)
# %%
# Set service times for Delay
node1.setService(jobclass1, Exp.fit_mean(1.0))
node1.setService(jobclass2, Exp.fit_mean(1.0/2.0))  # Mean = 0.5
node1.setService(jobclass3, Exp.fit_mean(1.0))
node1.setService(jobclass4, Exp.fit_mean(1.0))

# Set service times for Queue1
node2.setService(jobclass1, Exp.fit_mean(1.0/3.0))  # Mean = 0.333
node2.setService(jobclass2, Exp.fit_mean(1.0/4.0))  # Mean = 0.25
node2.setService(jobclass3, Exp.fit_mean(1.0/5.0))  # Mean = 0.2
node2.setService(jobclass4, Exp.fit_mean(1.0))

# Set service times for Queue2
node3.setService(jobclass1, Exp.fit_mean(1.0))
node3.setService(jobclass2, Exp.fit_mean(1.0/3.0))  # Mean = 0.333
node3.setService(jobclass3, Exp.fit_mean(1.0/5.0))  # Mean = 0.2
node3.setService(jobclass4, Exp.fit_mean(1.0/2.0))  # Mean = 0.5
# %%
# Block 3: routing with class switching
# Create routing matrices for each class-to-class transition
K = 4  # Number of classes
P = {}

# P[(i,j)] represents routing from class i to class j
# Matrix dimensions: [from_node, to_node]

# Class1 routing
P[(jobclass1, jobclass1)] = np.array([[0,1,0], [0,0,1], [0,0,0]])
P[(jobclass1, jobclass2)] = np.array([[0,0,0], [0,0,0], [1,0,0]])  # Class switch at Queue2
P[(jobclass1, jobclass3)] = np.array([[0,0,0], [0,0,0], [0,0,0]])
P[(jobclass1, jobclass4)] = np.array([[0,0,0], [0,0,0], [0,0,0]])

# Class2 routing
P[(jobclass2, jobclass1)] = np.array([[0,0,0], [0,0,0], [1,0,0]])  # Class switch at Queue2
P[(jobclass2, jobclass2)] = np.array([[0,1,0], [0,0,1], [0,0,0]])
P[(jobclass2, jobclass3)] = np.array([[0,0,0], [0,0,0], [0,0,0]])
P[(jobclass2, jobclass4)] = np.array([[0,0,0], [0,0,0], [0,0,0]])

# Class3 routing
P[(jobclass3, jobclass1)] = np.array([[0,0,0], [0,0,0], [0,0,0]])
P[(jobclass3, jobclass2)] = np.array([[0,0,0], [0,0,0], [0,0,0]])
P[(jobclass3, jobclass3)] = np.array([[0,1,0], [0,0,1], [0,0,0]])
P[(jobclass3, jobclass4)] = np.array([[0,0,0], [0,0,0], [1,0,0]])  # Class switch at Queue2

# Class4 routing
P[(jobclass4, jobclass1)] = np.array([[0,0,0], [0,0,0], [0,0,0]])
P[(jobclass4, jobclass2)] = np.array([[0,0,0], [0,0,0], [0,0,0]])
P[(jobclass4, jobclass3)] = np.array([[0,0,0], [0,0,0], [1,0,0]])  # Class switch at Queue2
P[(jobclass4, jobclass4)] = np.array([[0,0,1], [0,0,0], [0,0,0]])  # Delay -> Queue2

model.link(P)
# %%
# Set initial state for probability calculation
# State format: [station][class] where -1 means ignored
n = np.array([[-1,-1,-1,-1],   # Delay state (ignored)
              [-1,-1,-1,-1],   # Queue1 state (ignored)
              [1, 0, 2, 1]])   # Queue2 state: 1 Class1, 0 Class2, 2 Class3, 1 Class4

# Set state for each node
nodes = [node1, node2, node3]
for i in range(len(nodes)):
    if not np.any(n[i] == -1):  # Only set state if not ignored
        nodes[i].setState(n[i])
# %%
# Solve with CTMC for exact state probabilities
options = {'verbose': 1, 'seed': 23000}
solver_ctmc = CTMC(model, options)
Pr_ctmc = solver_ctmc.getProbAggr(node3)
print(f'Station 3 is in state {n[2].tolist()} with probability {Pr_ctmc}')
print(f'Pr_ctmc = {Pr_ctmc}')

# %%
# Solve with NC (Normalizing Constant) method
print("\n=== NC Solution (Normalizing Constants) ===")
solver_nc = NC(model, options)
Pr_nc = solver_nc.getProbAggr(node3)
print(f'Station 3 is in state {n[2].tolist()} with probability {Pr_nc}')
print(f'Pr_nc = {Pr_nc}')