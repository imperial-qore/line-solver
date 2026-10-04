from line_solver import *


def json_fcr_details():
    """Open 2-class network with a Finite Capacity Region."""
    model = Network('FCR_Details')

    source = Source(model, 'Source')
    queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.PS)
    sink = Sink(model, 'Sink')

    class1 = OpenClass(model, 'Class1')
    class2 = OpenClass(model, 'Class2')

    source.set_arrival(class1, Exp(1))
    source.set_arrival(class2, Exp(0.5))
    queue1.set_service(class1, Exp(3))
    queue1.set_service(class2, Exp(2))
    queue2.set_service(class1, Exp(4))
    queue2.set_service(class2, Exp(3))

    model.link(Network.serial_routing(source, queue1, queue2, sink))

    # Finite Capacity Region
    fcr = Region([queue1, queue2], [class1, class2])
    fcr.set_global_max_jobs(15)
    fcr.set_class_max_jobs(class1, 8)
    fcr.set_class_max_jobs(class2, 10)
    fcr.set_class_weight(class1, 1.0)
    fcr.set_class_weight(class2, 2.0)
    fcr.set_class_size(class1, 1)
    fcr.set_class_size(class2, 3)
    return model
