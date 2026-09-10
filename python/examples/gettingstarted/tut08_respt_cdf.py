# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
import matplotlib.pyplot as plt
# %%
model = Network('Model')

# Block 1: nodes
node = np.empty(2, dtype=object)
node[0] = Delay(model, 'Delay')
node[1] = Queue(model, 'Queue1', SchedStrategy.PS)

# Block 2: classes
jobclass = np.empty(2, dtype=object)
jobclass[0] = ClosedClass(model, 'Class1', 5, node[0], 0)
node[0].set_service(jobclass[0], Exp(1.0))
node[1].set_service(jobclass[0], Exp(0.5))

# Block 3: topology
model.link(Network.serial_routing(node[0], node[1]))
# %%
# Block 4: solution
RDfluid = FLD(model).cdf_resp_t()
RDsim = JMT(model, seed=23000, samples=10000).cdf_resp_t()
# %%
# Plot results
import matplotlib.pyplot as plt

if RDsim[1][0] is not None and len(RDsim[1][0]) > 0:
    plt.figure(figsize=(10, 6))
    plt.semilogx(RDsim[1][0][:, 1], 1 - RDsim[1][0][:, 0], 'r-', label='jmt-transient')
    plt.semilogx(RDfluid[1][0][:, 1], 1 - RDfluid[1][0][:, 0], 'k--', label='fluid-steady')
    plt.legend(loc='best')
    plt.ylabel('Pr(T > t)')
    plt.xlabel('time t')
    plt.title('Response Time Distribution and Percentiles')
    plt.grid(True, alpha=0.3)
    plt.show(block=False)
    plt.close('all')
else:
    print("No simulation results available for plotting")
# %%
# Compute CDF-derived scalar statistics
M = model.getNumberOfStations()
K = model.getNumberOfClasses()

avg_respt_from_cdf_sim = []
sq_coeff_of_variation_resp_t_from_cdf_sim = []
for i in range(M):
    avg_respt_from_cdf_sim.append([])
    sq_coeff_of_variation_resp_t_from_cdf_sim.append([])
    for c in range(K):
        cdf_data = RDsim[i][c]
        if cdf_data is not None and len(cdf_data) > 1:
            mean_val = sum((cdf_data[j+1][0] - cdf_data[j][0]) * cdf_data[j+1][1]
                          for j in range(len(cdf_data)-1))
            power_moment_2 = sum((cdf_data[j+1][0] - cdf_data[j][0]) * (cdf_data[j+1][1] ** 2)
                                for j in range(len(cdf_data)-1))
            variance = power_moment_2 - mean_val ** 2
            scv = variance / (mean_val ** 2) if mean_val > 0 else 0
            avg_respt_from_cdf_sim[i].append(mean_val)
            sq_coeff_of_variation_resp_t_from_cdf_sim[i].append(scv)
        else:
            avg_respt_from_cdf_sim[i].append(0)
            sq_coeff_of_variation_resp_t_from_cdf_sim[i].append(0)

avg_respt_from_cdf_fluid = []
sq_coeff_of_variation_resp_t_from_cdf_fluid = []
for i in range(M):
    avg_respt_from_cdf_fluid.append([])
    sq_coeff_of_variation_resp_t_from_cdf_fluid.append([])
    for c in range(K):
        cdf_data = RDfluid[i][c]
        if cdf_data is not None and len(cdf_data) > 1:
            mean_val = sum((cdf_data[j+1][0] - cdf_data[j][0]) * cdf_data[j+1][1]
                          for j in range(len(cdf_data)-1))
            power_moment_2 = sum((cdf_data[j+1][0] - cdf_data[j][0]) * (cdf_data[j+1][1] ** 2)
                                for j in range(len(cdf_data)-1))
            variance = power_moment_2 - mean_val ** 2
            scv = variance / (mean_val ** 2) if mean_val > 0 else 0
            avg_respt_from_cdf_fluid[i].append(mean_val)
            sq_coeff_of_variation_resp_t_from_cdf_fluid[i].append(scv)
        else:
            avg_respt_from_cdf_fluid[i].append(0)
            sq_coeff_of_variation_resp_t_from_cdf_fluid[i].append(0)

print('\nAverage Response Time from CDF (Simulation):')
print(avg_respt_from_cdf_sim)
print('\nAverage Response Time from CDF (Fluid):')
print(avg_respt_from_cdf_fluid)
print('\nSquared Coefficient of Variation from CDF (Simulation):')
print(sq_coeff_of_variation_resp_t_from_cdf_sim)
print('\nSquared Coefficient of Variation from CDF (Fluid):')
print(sq_coeff_of_variation_resp_t_from_cdf_fluid)