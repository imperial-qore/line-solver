from line_solver import *


def json_hetero_servers():
    """Closed 2-class network with heterogeneous servers."""
    model = Network('HeteroServers')

    delay = Delay(model, 'Delay')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    queue.set_number_of_servers(3)

    class1 = ClosedClass(model, 'Class1', 3, delay)
    class2 = ClosedClass(model, 'Class2', 2, delay)

    delay.set_service(class1, Exp(1))
    delay.set_service(class2, Exp(2))

    # Default service (required before hetero setup)
    queue.set_service(class1, Exp(3))
    queue.set_service(class2, Exp(2))

    # Heterogeneous server types
    fast = ServerType('Fast', 2, [class1, class2])
    slow = ServerType('Slow', 1, [class1])

    queue.add_server_type(fast)
    queue.add_server_type(slow)

    queue.set_hetero_sched_policy(HeteroSchedPolicy.ORDER)

    queue.set_hetero_service(class1, fast, Exp(5))
    queue.set_hetero_service(class2, fast, Exp(3))
    queue.set_hetero_service(class1, slow, Exp(1))

    P = model.init_routing_matrix()
    P.set(class1, class1, delay, queue, 1.0)
    P.set(class1, class1, queue, delay, 1.0)
    P.set(class2, class2, delay, queue, 1.0)
    P.set(class2, class2, queue, delay, 1.0)
    model.link(P)
    return model
