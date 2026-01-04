"""
CDF of Response Time - Open Network with Two Classes

This example demonstrates response time CDF computation for an open network
with two job classes using both Fluid and JMT solvers.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create open network
model = Network('myModel')

# Block 1: nodes
source = Source(model, 'Source')
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')

# Block 2: classes
job_class1 = OpenClass(model, 'Class1', 0)
job_class2 = OpenClass(model, 'Class2', 0)

source.setArrival(job_class1, Exp.fit_mean(4.0))
source.setArrival(job_class2, Exp.fit_mean(4.0))

queue1.setService(job_class1, Exp.fit_mean(1.0))
queue1.setService(job_class2, Exp.fit_mean(1.0))

queue2.setService(job_class1, Exp.fit_mean(1.0))
queue2.setService(job_class2, Exp.fit_mean(1.0))

# Block 3: topology
P = model.initRoutingMatrix()
P[0][0][0, 1] = 1  # Class1: Source -> Queue1
P[0][1][1, 2] = 1  # Class1: Queue1 -> Queue2
P[1][0][2, 3] = 1  # Class1: Queue2 -> Sink

P[1][1][0, 1] = 1  # Class2: Source -> Queue1
P[1][0][1, 2] = 1  # Class2: Queue1 -> Queue2
P[0][1][2, 3] = 1  # Class2: Queue2 -> Sink

model.link(P)

print('=== CDF Response Time - Open Network with Two Classes ===\n')

# Solve with Fluid solver
print('Computing CDF with Fluid solver...')
rd_fluid = FLD(model, {'iter_max': 300}).get_cdf_resp_t()

# Solve with JMT simulation
print('Computing transient CDF with JMT simulation...')
rd_sim = JMT(model, samples=int(1e4), seed=23000).get_tran_cdf_resp_t()

# Compute average response time from CDFs
avg_resp_t_from_cdf_fluid = []
avg_resp_t_from_cdf_sim = []

for i in range(model.getNumberOfStations()):
    avg_resp_t_from_cdf_fluid.append([])
    avg_resp_t_from_cdf_sim.append([])
    for c in range(model.getNumberOfClasses()):
        # Fluid
        cdf_data = rd_fluid[i][c]
        if cdf_data is not None and len(cdf_data) > 1:
            mean_val = sum((cdf_data[j+1][0] - cdf_data[j][0]) * cdf_data[j+1][1]
                          for j in range(len(cdf_data)-1))
            avg_resp_t_from_cdf_fluid[i].append(mean_val)
        else:
            avg_resp_t_from_cdf_fluid[i].append(0)

        # Simulation
        cdf_data = rd_sim[i][c]
        if cdf_data is not None and len(cdf_data) > 1:
            mean_val = sum((cdf_data[j+1][0] - cdf_data[j][0]) * cdf_data[j+1][1]
                          for j in range(len(cdf_data)-1))
            avg_resp_t_from_cdf_sim[i].append(mean_val)
        else:
            avg_resp_t_from_cdf_sim[i].append(0)

print('\nAverage Response Time from CDF (Fluid):')
for i, station_avg in enumerate(avg_resp_t_from_cdf_fluid):
    if i > 0:  # Skip source
        print(f'Station {i}: {station_avg}')

print('\nAverage Response Time from CDF (Simulation):')
for i, station_avg in enumerate(avg_resp_t_from_cdf_sim):
    if i > 0:  # Skip source
        print(f'Station {i}: {station_avg}')

print('\nNote: Response time tail distributions are available in rd_fluid and rd_sim')
print('      for detailed analysis and plotting.')
