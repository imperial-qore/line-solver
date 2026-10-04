"""Tutorial 15: two-class finite-capacity loss network."""

from line_solver import *


model = Network("LossNetwork")
source = Source(model, "Source")
delay = Delay(model, "Delay")
sink = Sink(model, "Sink")

class1 = OpenClass(model, "Class1")
class2 = OpenClass(model, "Class2")
source.set_arrival(class1, Exp(0.3))
source.set_arrival(class2, Exp(0.2))
delay.set_service(class1, Exp(1.0))
delay.set_service(class2, Exp(0.8))

routing = model.init_routing_matrix()
routing.set(class1, class1, source, delay, 1.0)
routing.set(class1, class1, delay, sink, 1.0)
routing.set(class2, class2, source, delay, 1.0)
routing.set(class2, class2, delay, sink, 1.0)
model.link(routing)

region = model.add_region(delay)
region.set_global_max_jobs(5)
region.set_class_max_jobs(class1, 3)
region.set_class_max_jobs(class2, 3)
region.set_drop_rule(class1, True)
region.set_drop_rule(class2, True)

print(SolverNC(model).avg_table())
print(SolverNC(model, method="mci", samples=100000, seed=23000).avg_table())
print(SolverLDES(model, seed=23000, samples=100000).avg_table())
