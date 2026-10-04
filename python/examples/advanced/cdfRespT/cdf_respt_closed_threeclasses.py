"""
CDF of Response Time - Closed Network with Three Classes

This example demonstrates response time CDF computation for a closed network
with three job classes that have class-switching behavior.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create closed model with 3 classes
model = Network('model')

delay = Delay(model, 'Delay')
queue = Queue(model, 'Queue2', SchedStrategy.PS)

# Class definitions
job_class1 = ClosedClass(model, 'Class1', 1, delay, 0)
job_class2 = ClosedClass(model, 'Class2', 0, delay, 0)
job_class3 = ClosedClass(model, 'Class3', 0, delay, 0)

# Disable completion for class 1 to maintain population
job_class1.completes = False

# Set service times
delay.setService(job_class1, Exp(1/1))
delay.setService(job_class2, Exp(1/1))
delay.setService(job_class3, Exp(1/1))

queue.setService(job_class1, Exp(1/1))
queue.setService(job_class2, Erlang(1/2, 2))
queue.setService(job_class3, Exp(1/0.01))

M = model.getNumberOfStations()
K = model.getNumberOfClasses()

# Routing matrix with class switching
P = [[None for _ in range(K)] for _ in range(K)]

# Class 1 routing: Delay -> Queue (with class switch to Class2)
P[0][0] = [[0, 1], [0, 0]]
P[0][1] = [[0, 0], [1, 0]]
P[0][2] = [[0, 0], [0, 0]]

# Class 2 routing: Queue -> Delay (with class switch to Class1)
P[1][0] = [[0, 0], [1, 0]]
P[1][1] = [[0, 1], [0, 0]]
P[1][2] = [[0, 0], [0, 0]]

# Class 3 routing: simple circular (no class switching)
P[2][0] = [[0, 0], [0, 0]]
P[2][1] = [[0, 0], [0, 0]]
P[2][2] = [[0, 1], [1, 0]]  # circul(M)

model.link(P)

print('=== CDF Response Time - Closed Network with Three Classes ===\n')

# Solve with Fluid solver using state-dependent method
print('Computing with Fluid solver (state-dependent method)...')
solver = FLD(model, {'method': 'statedep', 'iter_max': 100})
avg_respt = solver.get_avg_respt()
print('\nAverage Response Time:')
print(avg_respt)

# Get CDF
rd_fluid = solver.get_cdf_resp_t()

# Compute statistics from CDF
avg_respt_from_cdf = []
sq_coeff_of_variation_resp_t_from_cdf = []

for i in range(model.getNumberOfStations()):
    avg_respt_from_cdf.append([])
    sq_coeff_of_variation_resp_t_from_cdf.append([])
    for c in range(model.getNumberOfClasses()):
        cdf_data = rd_fluid[i][c]
        if cdf_data is not None and len(cdf_data) > 1:
            # Mean
            mean_val = sum((cdf_data[j+1][0] - cdf_data[j][0]) * cdf_data[j+1][1]
                          for j in range(len(cdf_data)-1))
            avg_respt_from_cdf[i].append(mean_val)

            # Second moment and variance
            power_moment_2 = sum((cdf_data[j+1][0] - cdf_data[j][0]) * (cdf_data[j+1][1] ** 2)
                                for j in range(len(cdf_data)-1))
            variance = power_moment_2 - mean_val ** 2
            scv = variance / (mean_val ** 2) if mean_val > 0 else 0
            sq_coeff_of_variation_resp_t_from_cdf[i].append(scv)
        else:
            avg_respt_from_cdf[i].append(0)
            sq_coeff_of_variation_resp_t_from_cdf[i].append(0)

print('\nAverage Response Time from CDF:')
print(avg_respt_from_cdf)

print('\nSquared Coefficient of Variation from CDF:')
print(sq_coeff_of_variation_resp_t_from_cdf)

print('\nNote: This example demonstrates class-switching behavior:')
print('  - Class 1 jobs switch to Class 2 at the queue')
print('  - Class 2 jobs switch back to Class 1 at the delay')
print('  - Class 3 jobs follow a simple circular route without switching')
