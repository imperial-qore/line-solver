"""
CDF of Response Time - Varying Population Sizes

This example demonstrates how response time CDFs change as the number
of jobs in a closed network increases.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Test with different population sizes
n_jobs_list = [1, 4, 8]

print('=== CDF Response Time - Varying Population Sizes ===\n')

results = []

for N in n_jobs_list:
    print(f'\n--- Population N = {N} jobs ---')

    model = Network('model')

    delay = Delay(model, 'Delay')
    queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.PS)

    job_class = ClosedClass(model, 'Class1', N, delay, 0)
    job_class.completes = False

    delay.setService(job_class, Exp(1/1))
    queue1.setService(job_class, Exp(1/2))
    queue2.setService(job_class, Exp(1/2))

    # Circular routing through all stations (circul(M))
    # circul(3) = [[0,1,0], [0,0,1], [1,0,0]]
    model.link([[0, 1, 0], [0, 0, 1], [1, 0, 0]])

    # Solve with Fluid solver
    solver = FLD(model, {'iter_max': 100})
    avg_respt = solver.get_avg_respt()
    rd_fluid = solver.get_cdf_resp_t()

    # Compute statistics from CDF
    avg_respt_from_cdf = []
    sq_coeff_of_variation_resp_t_from_cdf = []

    for c in range(model.getNumberOfClasses()):
        avg_respt_from_cdf.append([])
        sq_coeff_of_variation_resp_t_from_cdf.append([])
        for i in range(model.getNumberOfStations()):
            cdf_data = rd_fluid[i][c]
            if cdf_data is not None and len(cdf_data) > 1:
                # Mean
                mean_val = sum((cdf_data[j+1][0] - cdf_data[j][0]) * cdf_data[j+1][1]
                              for j in range(len(cdf_data)-1))
                avg_respt_from_cdf[c].append(mean_val)

                # Second moment and variance
                power_moment_2 = sum((cdf_data[j+1][0] - cdf_data[j][0]) * (cdf_data[j+1][1] ** 2)
                                    for j in range(len(cdf_data)-1))
                variance = power_moment_2 - mean_val ** 2
                scv = variance / (mean_val ** 2) if mean_val > 0 else 0
                sq_coeff_of_variation_resp_t_from_cdf[c].append(scv)
            else:
                avg_respt_from_cdf[c].append(0)
                sq_coeff_of_variation_resp_t_from_cdf[c].append(0)

    print(f'Average Response Time: {avg_respt}')
    print(f'Average Response Time from CDF: {avg_respt_from_cdf}')
    print(f'Squared Coefficient of Variation from CDF: {sq_coeff_of_variation_resp_t_from_cdf}')

    results.append({
        'N': N,
        'avg_respt': avg_respt,
        'avg_respt_from_cdf': avg_respt_from_cdf,
        'scv_from_cdf': sq_coeff_of_variation_resp_t_from_cdf,
        'cdf': rd_fluid
    })

print('\n\n=== Summary ===')
print('\nAs population increases, response times increase due to queueing effects.')
print('The CDF shifts to the right, showing higher probability of longer response times.')
print('\nResponse time statistics for each population:')
for result in results:
    print(f"\nN = {result['N']}: Mean Response Times = {result['avg_respt_from_cdf']}")
