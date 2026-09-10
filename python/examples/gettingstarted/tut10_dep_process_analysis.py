# %%
from line_solver import *
import numpy as np
# %%
model, source, queue, sink, oclass = gallery_merl1()
# %%
# Use cutoff=20 to limit state space size for phase-type distributions
# MATLAB auto-calculates cutoff as ceil(6000^(1/(M*K))) for open models
solver = CTMC(model, cutoff=150, seed=23000)
# %%
sa = solver.sample_sys_aggr(5000)
# %%
ind = model.get_node_index(queue) - 1  # Convert from 1-based to 0-based for event comparison

# Filter events for departures from the queue
dep_times = []
for event in sa.event:
    if event.node == ind and event.event == "DEP":
        dep_times.append(event.t)

print(f"Found {len(dep_times)} departure events from queue")

if len(dep_times) > 1:
    inter_dep_times = np.diff(dep_times)
    # Estimated squared coefficient of variation of departures
    scv_d_est = np.var(inter_dep_times) / np.mean(inter_dep_times)**2
    print(f"Simulated SCV of departures: {scv_d_est}")
else:
    print("Error: Insufficient departure events found")
    print(f"Total events generated: {len(sa.event)}")
    print(f"Sample of events: {[(e.node, e.event, e.t) for e in sa.event[:5]]}")
# %%
# Get queue utilization and waiting time
util = solver.avg_util()
util_queue = util[queue][0]
avg_wait_time = solver.avg_wait_t()
avg_wait_time_queue = avg_wait_time[ind]
# %%
# Marshall's exact formula for SCV of departures
scv_a = source.get_arrival(oclass).get_scv()
svc_rate = queue.get_service(oclass).get_rate()
scv_s = queue.get_service(oclass).get_scv()
scv_d = scv_a + 2*util_queue**2*scv_s - 2*util_queue*(1-util_queue)*svc_rate*avg_wait_time_queue
print(f"Theoretical SCV of departures (Marshall's formula): {scv_d}")

# Calculate relative error between simulated and theoretical SCV  
if 'scv_d_est' in locals():
    relative_error = abs(scv_d_est - scv_d[0]) / scv_d[0] * 100
    print(f"\n=== Departure Process Analysis Results ===")
    print(f"Simulated SCV of departures:   {scv_d_est:.6f}")
    print(f"Theoretical SCV (Marshall):    {scv_d[0]:.6f}")
    print(f"Relative error:                {relative_error:.2f}%")
else:
    print("\nCannot calculate relative error - simulation failed")