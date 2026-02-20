"""
CDF of Response Time - Closed Network

This example demonstrates how to compute the cumulative distribution function (CDF)
of response times in a closed queueing network using both JMT simulation and Fluid solver.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create closed model with 2 stations
model = Network('model')

delay = Delay(model, 'Delay')
queue = Queue(model, 'Queue2', SchedStrategy.PS)

job_class = ClosedClass(model, 'Class1', 1, delay, 0)

# Set service distributions
serv_proc1 = Exp(1/0.1)
delay.setService(job_class, serv_proc1)
serv_proc2 = Erlang.fit_mean_and_scv(1, 1/3)
queue.setService(job_class, serv_proc2)

# Serial routing between stations
model.link(Network.serial_routing([delay, queue]))

print('=== CDF Response Time - Closed Network ===\n')

# Solve with JMT simulation
print('Computing CDF with JMT simulation...')
solver_jmt = JMT(model, seed=23000, samples=10000)
rd_sim = solver_jmt.get_cdf_resp_t()

# Compute statistics from CDF (simulation)
avg_respt_from_cdf_sim = []
sq_coeff_of_variation_resp_t_from_cdf_sim = []

for i in range(model.getNumberOfStations()):
    avg_respt_from_cdf_sim.append([])
    sq_coeff_of_variation_resp_t_from_cdf_sim.append([])
    for c in range(model.getNumberOfClasses()):
        cdf_data = rd_sim[i][c]
        if cdf_data is not None and len(cdf_data) > 1:
            # Compute mean from CDF
            mean_val = sum((cdf_data[j+1][0] - cdf_data[j][0]) * cdf_data[j+1][1]
                          for j in range(len(cdf_data)-1))
            avg_respt_from_cdf_sim[i].append(mean_val)

            # Compute second moment
            power_moment_2 = sum((cdf_data[j+1][0] - cdf_data[j][0]) * (cdf_data[j+1][1] ** 2)
                                for j in range(len(cdf_data)-1))
            variance = power_moment_2 - mean_val ** 2
            scv = variance / (mean_val ** 2) if mean_val > 0 else 0
            sq_coeff_of_variation_resp_t_from_cdf_sim[i].append(scv)
        else:
            avg_respt_from_cdf_sim[i].append(0)
            sq_coeff_of_variation_resp_t_from_cdf_sim[i].append(0)

# Solve with Fluid solver
print('Computing CDF with Fluid solver...')
solver_fluid = FLD(model)
rd_fluid = solver_fluid.get_cdf_resp_t()

# Compute statistics from CDF (fluid)
avg_respt_from_cdf_fluid = []
sq_coeff_of_variation_resp_t_from_cdf_fluid = []

for i in range(model.getNumberOfStations()):
    avg_respt_from_cdf_fluid.append([])
    sq_coeff_of_variation_resp_t_from_cdf_fluid.append([])
    for c in range(model.getNumberOfClasses()):
        cdf_data = rd_fluid[i][c]
        if cdf_data is not None and len(cdf_data) > 1:
            # Compute mean from CDF
            mean_val = sum((cdf_data[j+1][0] - cdf_data[j][0]) * cdf_data[j+1][1]
                          for j in range(len(cdf_data)-1))
            avg_respt_from_cdf_fluid[i].append(mean_val)

            # Compute second moment
            power_moment_2 = sum((cdf_data[j+1][0] - cdf_data[j][0]) * (cdf_data[j+1][1] ** 2)
                                for j in range(len(cdf_data)-1))
            variance = power_moment_2 - mean_val ** 2
            scv = variance / (mean_val ** 2) if mean_val > 0 else 0
            sq_coeff_of_variation_resp_t_from_cdf_fluid[i].append(scv)
        else:
            avg_respt_from_cdf_fluid[i].append(0)
            sq_coeff_of_variation_resp_t_from_cdf_fluid[i].append(0)

print('\nSince there is a single job, mean and squared coefficient of variation')
print('of response times are close, up to fluid approximation precision, to those')
print('of the service time distribution.\n')

avg_respt_from_theory = [serv_proc1.getMean(), serv_proc2.getMean()]
sq_coeff_of_variation_resp_t_from_theory = [serv_proc1.getSCV(), serv_proc2.getSCV()]

print('Average Response Time from Theory:')
print(avg_respt_from_theory)
print('\nAverage Response Time from CDF (Simulation):')
print(avg_respt_from_cdf_sim)
print('\nAverage Response Time from CDF (Fluid):')
print(avg_respt_from_cdf_fluid)
print('\nSquared Coefficient of Variation from Theory:')
print(sq_coeff_of_variation_resp_t_from_theory)
print('\nSquared Coefficient of Variation from CDF (Simulation):')
print(sq_coeff_of_variation_resp_t_from_cdf_sim)
print('\nSquared Coefficient of Variation from CDF (Fluid):')
print(sq_coeff_of_variation_resp_t_from_cdf_fluid)
