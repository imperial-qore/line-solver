from line_solver import *
import numpy as np


def json_classcap_droprule():
    """Open 2-class network with per-class capacity, drop rules, and load-dependence."""
    model = Network('ClassCap_DropRule')

    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.PS)
    sink = Sink(model, 'Sink')

    queue.set_number_of_servers(2)
    queue.set_capacity(20)

    class1 = OpenClass(model, 'Class1')
    class2 = OpenClass(model, 'Class2')

    source.set_arrival(class1, Exp(1))
    source.set_arrival(class2, Exp(0.5))
    queue.set_service(class1, Exp(3))
    queue.set_service(class2, Exp(2))

    # Per-class capacity
    queue.set_class_capacity(class1, 8)
    queue.set_class_capacity(class2, 15)

    # Per-class drop rules. WAITQ is not usable here: no solver honours "wait
    # upstream" for an open class at a plain finite capacity (CTMC drops the
    # arrival, JMT blocks at the source and ignores the cap), so the capacity
    # refresh rejects that combination. BAS is the honoured blocking policy.
    queue.set_drop_rule(class1, DropStrategy.DROP)
    queue.set_drop_rule(class2, DropStrategy.BAS)

    # Load-dependent scaling
    queue.set_load_dependence(np.array([1.0, 0.9, 0.8, 0.7]))

    model.link(Network.serial_routing(source, queue, sink))
    return model
