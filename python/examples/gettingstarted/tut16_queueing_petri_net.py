"""Tutorial 16: closed queueing Petri net."""

from line_solver import *


population = 4
model = Network("QueueingPetriNet")
cpu = Place(model, "CPU", SchedStrategy.FCFS)
think = Place(model, "Think", SchedStrategy.INF)
jobs = ClosedClass(model, "Jobs", population, think, 0)
cpu.set_service(jobs, Exp(1.5))
think.set_service(jobs, Exp(0.5))

to_cpu = Transition(model, "toCPU")
mode1 = to_cpu.add_mode("m1")
to_cpu.set_timing_strategy(mode1, TimingStrategy.IMMEDIATE)
to_cpu.set_enabling_conditions(mode1, jobs, think, 1)
to_cpu.set_firing_outcome(mode1, jobs, cpu, 1)

to_think = Transition(model, "toThink")
mode2 = to_think.add_mode("m2")
to_think.set_timing_strategy(mode2, TimingStrategy.IMMEDIATE)
to_think.set_enabling_conditions(mode2, jobs, cpu, 1)
to_think.set_firing_outcome(mode2, jobs, think, 1)

routing = model.init_routing_matrix()
routing.set(jobs, jobs, think, to_cpu, 1.0)
routing.set(jobs, jobs, to_cpu, cpu, 1.0)
routing.set(jobs, jobs, cpu, to_think, 1.0)
routing.set(jobs, jobs, to_think, think, 1.0)
model.link(routing)
think.set_state(population)
cpu.set_state(0)

print(SolverLDES(model, seed=23000, samples=200000).avg_table())

reference = Network("Reference")
delay = Delay(reference, "Think")
queue = Queue(reference, "CPU", SchedStrategy.FCFS)
reference_jobs = ClosedClass(reference, "Jobs", population, delay, 0)
delay.set_service(reference_jobs, Exp(0.5))
queue.set_service(reference_jobs, Exp(1.5))
reference.link(Network.serial_routing(delay, queue))
print(SolverMVA(reference).avg_table())
