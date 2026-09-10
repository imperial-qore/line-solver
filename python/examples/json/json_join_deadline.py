from line_solver import *


def json_join_deadline():
    """Closed 1-class fork-join with quorum join and deadline."""
    model = Network('JoinDeadline')

    delay = Delay(model, 'Delay')
    fork = Fork(model, 'Fork')
    queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
    queue3 = Queue(model, 'Queue3', SchedStrategy.FCFS)
    join = Join(model, 'Join', fork)

    class1 = ClosedClass(model, 'Class1', 5, delay)
    class1.setDeadline(5.0)

    delay.set_service(class1, Exp(1))
    queue1.set_service(class1, Exp(3))
    queue2.set_service(class1, Exp(4))
    queue3.set_service(class1, Exp(2))

    join.set_strategy(class1, JoinStrategy.QUORUM)
    join.set_required(class1, 2)

    P = model.init_routing_matrix()
    P.set(class1, class1, delay, fork, 1.0)
    P.set(class1, class1, fork, queue1, 1.0)
    P.set(class1, class1, fork, queue2, 1.0)
    P.set(class1, class1, fork, queue3, 1.0)
    P.set(class1, class1, queue1, join, 1.0)
    P.set(class1, class1, queue2, join, 1.0)
    P.set(class1, class1, queue3, join, 1.0)
    P.set(class1, class1, join, delay, 1.0)
    model.link(P)
    return model
