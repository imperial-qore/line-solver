"""
CDF of Response Time - Different Service Distributions

This example shows how response time CDFs vary with different service time
distributions (Exp, Erlang, HyperExp) in a closed queueing network.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create closed model with 2 classes
model = Network('model')

delay = Delay(model, 'Delay')
queue = Queue(model, 'Queue1', SchedStrategy.PS)

# Class 1: 1 job with exponential service
job_class1 = ClosedClass(model, 'Class1', 1, delay, 0)
delay.setService(job_class1, Exp.fit_mean(1.0))
queue.setService(job_class1, Exp.fit_mean(2.0))

# Class 2: 3 jobs with non-exponential service distributions
job_class2 = ClosedClass(model, 'Class2', 3, delay, 0)
delay.setService(job_class2, Erlang.fit_mean_and_order(4.0, 2))
queue.setService(job_class2, HyperExp.fit_mean_and_scv(5.0, 30.0))

# Routing matrix
P = [[None, None], [None, None]]
P[0][0] = [[0, 1], [1, 0]]  # circul(2)
P[0][1] = [[0, 0], [0, 0]]
P[1][0] = [[0, 0], [0, 0]]
P[1][1] = [[0, 1], [1, 0]]  # circul(2)

model.link(P)

print('=== CDF Response Time - Different Service Distributions ===\n')

# Solve with Fluid solver
print('Computing CDF with Fluid solver (steady-state)...')
rd_fluid = FLD(model).get_cdf_resp_t()

# Solve with JMT simulation
print('Computing transient CDF with JMT simulation...')
rd_sim = JMT(model, samples=int(1e5), seed=23000).get_tran_cdf_resp_t()

# Compute statistics from CDFs
avg_respt_from_cdf_fluid = []
avg_respt_from_cdf_sim = []

for i in range(model.getNumberOfStations()):
    avg_respt_from_cdf_fluid.append([])
    avg_respt_from_cdf_sim.append([])
    for c in range(model.getNumberOfClasses()):
        # Fluid
        cdf_data = rd_fluid[i][c]
        if cdf_data is not None and len(cdf_data) > 1:
            mean_val = sum((cdf_data[j+1][0] - cdf_data[j][0]) * cdf_data[j+1][1]
                          for j in range(len(cdf_data)-1))
            avg_respt_from_cdf_fluid[i].append(mean_val)
        else:
            avg_respt_from_cdf_fluid[i].append(0)

        # Simulation
        cdf_data = rd_sim[i][c]
        if cdf_data is not None and len(cdf_data) > 1:
            mean_val = sum((cdf_data[j+1][0] - cdf_data[j][0]) * cdf_data[j+1][1]
                          for j in range(len(cdf_data)-1))
            avg_respt_from_cdf_sim[i].append(mean_val)
        else:
            avg_respt_from_cdf_sim[i].append(0)

print('\nAverage Response Time from CDF (Fluid):')
print(avg_respt_from_cdf_fluid)

print('\nAverage Response Time from CDF (Simulation):')
print(avg_respt_from_cdf_sim)

# Get service process names
service_procs = []
stations = model.getStations()
for station in stations:
    service_procs.append([])
    for c in range(model.getNumberOfClasses()):
        try:
            proc = station.getServiceProcess(c)
            if proc is not None:
                service_procs[-1].append(proc.getName())
            else:
                service_procs[-1].append('N/A')
        except:
            service_procs[-1].append('N/A')

print('\nService Processes by Station and Class:')
for i, procs in enumerate(service_procs):
    print(f'Station {i}: {procs}')

print('\nNote: Class 1 uses Exponential service distributions.')
print('      Class 2 uses Erlang (delay) and HyperExponential (queue) distributions.')
print('      This shows how different service distributions affect response time CDFs.')
