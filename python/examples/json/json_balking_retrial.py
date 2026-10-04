from line_solver import *


def json_balking_retrial():
    """Open 2-class network with balking, retrial, and patience."""
    model = Network('Balking_Retrial')

    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    queue.set_capacity(15)

    class1 = OpenClass(model, 'Class1')
    class2 = OpenClass(model, 'Class2')

    source.set_arrival(class1, Exp(1))
    source.set_arrival(class2, Exp(0.5))
    queue.set_service(class1, Exp(2))
    queue.set_service(class2, Exp(3))

    # Class1: balking based on queue length
    queue.set_balking(class1, BalkingStrategy.QUEUE_LENGTH,
                      [(5, 10, 0.3), (11, float('inf'), 1.0)])

    # Class1: patience (reneging)
    queue.set_patience(class1, Exp(0.1))

    # Class2: retrial with max attempts
    queue.set_retrial(class2, Exp(0.5), 3)

    model.link(Network.serial_routing(source, queue, sink))
    return model
