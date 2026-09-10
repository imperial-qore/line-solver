from line_solver import *


def json_signal_classes():
    """Open network with G-network negative signal."""
    model = Network('SignalClasses')

    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    class1 = OpenClass(model, 'Class1')
    signal1 = OpenSignal(model, 'Signal1', SignalType.NEGATIVE)
    signal1.forJobClass(class1)
    signal1.setRemovalPolicy(RemovalPolicy.RANDOM)

    source.set_arrival(class1, Exp(2))
    source.set_arrival(signal1, Exp(0.5))
    queue.set_service(class1, Exp(5))
    queue.set_service(signal1, Immediate())

    P = model.init_routing_matrix()
    P.set(class1, class1, source, queue, 1.0)
    P.set(class1, class1, queue, sink, 1.0)
    P.set(signal1, signal1, source, queue, 1.0)
    P.set(signal1, signal1, queue, sink, 1.0)
    model.link(P)
    return model
