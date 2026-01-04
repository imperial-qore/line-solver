"""
Finite Capacity Region with Dropping vs M/M/1/K

This example shows that an FCR with dropping around a single queue
behaves like an M/M/1/K queue, where K is the FCR capacity.
Jobs arriving when the region is full are dropped (lost).
"""

from line_solver import *

# Parameters
arrival_rate = 0.8
service_rate = 1.0
K = 2  # Capacity limit

# Model 1: Queue with FCR (dropping)
model1 = Network('FCR Dropping')

source1 = Source(model1, 'Source')
queue1 = Queue(model1, 'Queue', SchedStrategy.FCFS)
sink1 = Sink(model1, 'Sink')

jobclass1 = OpenClass(model1, 'Class1', 0)

source1.set_arrival(jobclass1, Exp.fit_rate(arrival_rate))
queue1.set_service(jobclass1, Exp.fit_rate(service_rate))

P1 = model1.init_routing_matrix()
P1.set(jobclass1, jobclass1, source1, queue1, 1.0)
P1.set(jobclass1, jobclass1, queue1, sink1, 1.0)
model1.link(P1)

# Add FCR with dropping - jobs are lost when region is full
fcr = model1.add_region('FCR', queue1)
fcr.set_global_max_jobs(K)
fcr.set_drop_rule(jobclass1, True)  # True = drop jobs

# Model 2: M/M/1/K using queue capacity
model2 = Network('MM1K')

source2 = Source(model2, 'Source')
queue2 = Queue(model2, 'Queue', SchedStrategy.FCFS)
queue2.set_number_of_servers(1)
queue2.set_cap(K)  # Set queue capacity to K
sink2 = Sink(model2, 'Sink')

jobclass2 = OpenClass(model2, 'Class1', 0)

source2.set_arrival(jobclass2, Exp.fit_rate(arrival_rate))
queue2.set_service(jobclass2, Exp.fit_rate(service_rate))

P2 = model2.init_routing_matrix()
P2.set(jobclass2, jobclass2, source2, queue2, 1.0)
P2.set(jobclass2, jobclass2, queue2, sink2, 1.0)
model2.link(P2)

# Solve both models
solver1 = JMT(model1, seed=23000, samples=100000, verbose=VerboseLevel.SILENT)
solver2 = JMT(model2, seed=23000, samples=100000, verbose=VerboseLevel.SILENT)

# Compare results
print(f'\n=== Comparison: FCR Dropping vs M/M/1/K (K={K}) ===\n')

avg_table1 = solver1.get_avg_table()
avg_table2 = solver2.get_avg_table()

# Extract queue metrics (row 1 is the Queue node)
queue_idx = 1
print('Queue Metrics:')
print('                    FCR Dropping    M/M/1/K')
print(f'Queue Length:       {avg_table1.QLen[queue_idx]:.4f}          {avg_table2.QLen[queue_idx]:.4f}')
print(f'Utilization:        {avg_table1.Util[queue_idx]:.4f}          {avg_table2.Util[queue_idx]:.4f}')
print(f'Response Time:      {avg_table1.RespT[queue_idx]:.4f}          {avg_table2.RespT[queue_idx]:.4f}')
print(f'Throughput:         {avg_table1.Tput[queue_idx]:.4f}          {avg_table2.Tput[queue_idx]:.4f}')

# Theoretical M/M/1/K values
rho = arrival_rate / service_rate
if rho != 1:
    p0 = (1 - rho) / (1 - rho**(K + 1))
    pK = p0 * rho**K
    theoretical_Q = rho / (1 - rho) - (K + 1) * rho**(K + 1) / (1 - rho**(K + 1))
    lambda_eff = arrival_rate * (1 - pK)
    theoretical_U = lambda_eff / service_rate
    theoretical_X = lambda_eff
    theoretical_R = theoretical_Q / theoretical_X
    blocking_prob = pK
else:
    p0 = 1 / (K + 1)
    pK = p0
    theoretical_Q = K / 2
    lambda_eff = arrival_rate * (1 - pK)
    theoretical_U = lambda_eff / service_rate
    theoretical_X = lambda_eff
    theoretical_R = theoretical_Q / theoretical_X
    blocking_prob = pK

print('\nTheoretical M/M/1/K:')
print(f'Queue Length:       {theoretical_Q:.4f}')
print(f'Utilization:        {theoretical_U:.4f}')
print(f'Response Time:      {theoretical_R:.4f}')
print(f'Throughput:         {theoretical_X:.4f}')
print(f'Blocking Prob:      {blocking_prob:.4f}')

print('\n==> FCR with dropping matches M/M/1/K (jobs lost when full)')
