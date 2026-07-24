# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Create network
model = Network('model')

# Block 1: nodes with mixed scheduling
node1 = Delay(model, 'Delay')
node2 = Queue(model, 'Queue1', SchedStrategy.PS)   # Processor Sharing
node3 = Queue(model, 'Queue2', SchedStrategy.DPS)  # Discriminatory Processor Sharing

# Block 2: classes with different populations
jobclass1 = ClosedClass(model, 'Class1', 2, node1, 0)  # 2 jobs
jobclass2 = ClosedClass(model, 'Class2', 1, node1, 0)  # 1 job

print("Network created with mixed scheduling:")
print(f"  {node1.get_name()}: Delay node")
print(f"  {node2.get_name()}: PS (Processor Sharing)")
print(f"  {node3.get_name()}: DPS (Discriminatory Processor Sharing)")
print(f"\nClass populations: Class1={jobclass1.get_population()}, Class2={jobclass2.get_population()}")
# %%
# Set service distributions
node1.set_service(jobclass1, Exp.fit_mean(1.0/3.0))   # Delay, Class1: mean=1/3
node1.set_service(jobclass2, Exp.fit_mean(1.0/0.5))   # Delay, Class2: mean=1/0.5=2

# PS Queue - weights are ignored in PS scheduling
w1_ps = 5  # This weight is ignored since the node is PS
w2_ps = 1  # This weight is ignored since the node is PS
node2.set_service(jobclass1, Exp.fit_mean(1.0/0.1), w1_ps)
node2.set_service(jobclass2, Exp.fit_mean(1.0/1.0), w2_ps)

# DPS Queue - weights matter for service differentiation
w1_dps = 1  # Lower weight for Class1
w2_dps = 5  # Higher weight for Class2 (gets priority)
node3.set_service(jobclass1, Exp.fit_mean(1.0/0.1), w1_dps)
node3.set_service(jobclass2, Exp.fit_mean(1.0/1.0), w2_dps)

print("Service parameters configured:")
print(f"  Delay: Class1=Exp(3), Class2=Exp(0.5)")
print(f"  Queue1 (PS): Class1=Exp(0.1), Class2=Exp(1) - weights ignored")
print(f"  Queue2 (DPS): Class1=Exp(0.1,w={w1_dps}), Class2=Exp(1,w={w2_dps}) - Class2 has priority")
# %%
# Set up routing matrix
P = model.init_routing_matrix()

# Class1 routing: probabilistic from Delay (30% Queue1, 70% Queue2)
P.set(jobclass1, jobclass1, node1, node2, 0.3)  # Delay -> Queue1 (30%)
P.set(jobclass1, jobclass1, node1, node3, 0.7)  # Delay -> Queue2 (70%)
P.set(jobclass1, jobclass1, node2, node1, 1.0)  # Queue1 -> Delay (100%)
P.set(jobclass1, jobclass1, node3, node1, 1.0)  # Queue2 -> Delay (100%)

# Class2 routing: probabilistic from Delay (70% Queue1, 30% Queue2)
P.set(jobclass2, jobclass2, node1, node2, 0.7)  # Delay -> Queue1 (70%)
P.set(jobclass2, jobclass2, node1, node3, 0.3)  # Delay -> Queue2 (30%)
P.set(jobclass2, jobclass2, node2, node1, 1.0)  # Queue1 -> Delay (100%)
P.set(jobclass2, jobclass2, node3, node1, 1.0)  # Queue2 -> Delay (100%)

model.link(P)

print("Routing configured:")
print("  Class1: 30% Queue1, 70% Queue2 (more traffic to DPS queue)")
print("  Class2: 70% Queue1, 30% Queue2 (more traffic to PS queue)")
print("  All queues return 100% to Delay")
# %%
# Configure multiple solvers for comparison
print("=== Multi-Solver Analysis ===")
solver_list = []

# Add different solvers with their respective options
# CTMC
solver_list.append(CTMC(model, verbose=True))

# JMT
solver_list.append(JMT(model, verbose=True, samples=10000, seed=23000))

# Fluid
solver_list.append(FLD(model, verbose=True))

# MVA
solver_list.append(MVA(model, verbose=True))

# LDES (matches MATLAB)
solver_list.append(LDES(model, verbose=True, samples=10000, seed=23000))

print(f"Configured {len(solver_list)} solvers for comparison")
# %%
# Solve with each solver and compare results
avg_tables = []

for s, solver in enumerate(solver_list):
    print(f'\n=== SOLVER: {solver.get_name()} ===')
    try:
        avg_table = solver.avg_table()
        avg_tables.append(avg_table)
        print(avg_table)
    except RuntimeError as e:
        print(f"Solver not supported for this model: {e}")
        avg_tables.append(None)