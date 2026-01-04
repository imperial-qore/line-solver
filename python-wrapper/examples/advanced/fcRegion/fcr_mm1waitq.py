"""
Finite Capacity Region with Blocking vs M/M/1

This example shows that an FCR with blocking (waitq) around a single
queue behaves identically to a standard M/M/1 queue, since jobs
simply wait when the region is "full" (no actual capacity limit effect).
"""

from line_solver import *

# Parameters
arrival_rate = 0.5
service_rate = 1.0

# Model 1: Queue with FCR (blocking)
model1 = Network('FCR Blocking')

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

# Add FCR with blocking - jobs wait when region is full
fcr = model1.add_region('FCR', queue1)
fcr.set_global_max_jobs(10)
fcr.set_drop_rule(jobclass1, False)  # False = block (wait)

# Model 2: Standard M/M/1 (no FCR)
model2 = Network('MM1')

source2 = Source(model2, 'Source')
queue2 = Queue(model2, 'Queue', SchedStrategy.FCFS)
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
print('\n=== Comparison: FCR Blocking vs M/M/1 ===\n')

avg_table1 = solver1.get_avg_table()
avg_table2 = solver2.get_avg_table()

# Extract queue metrics (row 1 is the Queue node)
queue_idx = 1
print('Queue Metrics:')
print('                    FCR Blocking    M/M/1')
print(f'Queue Length:       {avg_table1.QLen[queue_idx]:.4f}          {avg_table2.QLen[queue_idx]:.4f}')
print(f'Utilization:        {avg_table1.Util[queue_idx]:.4f}          {avg_table2.Util[queue_idx]:.4f}')
print(f'Response Time:      {avg_table1.RespT[queue_idx]:.4f}          {avg_table2.RespT[queue_idx]:.4f}')
print(f'Throughput:         {avg_table1.Tput[queue_idx]:.4f}          {avg_table2.Tput[queue_idx]:.4f}')

# Theoretical M/M/1 values
rho = arrival_rate / service_rate
theoretical_Q = rho / (1 - rho)
theoretical_U = rho
theoretical_R = 1 / (service_rate - arrival_rate)
theoretical_X = arrival_rate

print('\nTheoretical M/M/1:')
print(f'Queue Length:       {theoretical_Q:.4f}')
print(f'Utilization:        {theoretical_U:.4f}')
print(f'Response Time:      {theoretical_R:.4f}')
print(f'Throughput:         {theoretical_X:.4f}')

print('\n==> FCR with blocking matches M/M/1 (jobs wait, no loss)')
