"""
Mixed Network with Multi-Server PS Queues

This example demonstrates:
- Mixed network with 3 closed jobs and open arrivals
- 5 PS queues with varying numbers of servers (1-5 servers)
- Closed class visits queues 1-4, open class visits queues 1-3 and 5
- Different routing patterns for closed and open classes
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    node = np.empty(5, dtype=object)
    node[0] = Queue(model, 'Queue1', SchedStrategy.PS)
    node[0].set_num_servers(1)
    node[1] = Queue(model, 'Queue2', SchedStrategy.PS)
    node[1].set_num_servers(2)
    node[2] = Queue(model, 'Queue3', SchedStrategy.PS)
    node[2].set_num_servers(3)
    node[3] = Queue(model, 'Queue4', SchedStrategy.PS)  # only closed class
    node[3].set_num_servers(4)
    node[4] = Queue(model, 'Queue5', SchedStrategy.PS)  # only open class
    node[4].set_num_servers(5)

    source = Source(model, 'Source')
    sink = Sink(model, 'Sink')

    jobclass = np.empty(2, dtype=object)
    jobclass[0] = ClosedClass(model, 'ClosedClass', 3, node[0], 0)
    jobclass[1] = OpenClass(model, 'OpenClass', 0)

    for i in range(5):
        node[i].set_service(jobclass[0], Exp(i + 1))
        node[i].set_service(jobclass[1], Exp((i + 1) ** 0.5))

    source.set_arrival(jobclass[1], Exp(0.3))

    # Routing
    P = model.init_routing_matrix()

    # Closed class: serial routing through queues 1-4
    routing_closed = Network.serial_routing([node[0], node[1], node[2], node[3]])
    for i in range(4):
        for j in range(4):
            P.set(jobclass[0], jobclass[0], node[i], node[j], routing_closed[i][j])

    # Open class: source -> queue1 -> queue2 -> queue3 -> queue5 -> sink
    routing_open = Network.serial_routing([source, node[0], node[1], node[2], node[4], sink])
    P.set(jobclass[1], jobclass[1], source, node[0], routing_open[0][1])
    P.set(jobclass[1], jobclass[1], node[0], node[1], routing_open[1][2])
    P.set(jobclass[1], jobclass[1], node[1], node[2], routing_open[2][3])
    P.set(jobclass[1], jobclass[1], node[2], node[4], routing_open[3][4])
    P.set(jobclass[1], jobclass[1], node[4], sink, routing_open[4][5])

    model.link(P)

    # Run multiple solvers
    solver = np.array([], dtype=object)
    solver = np.append(solver, CTMC(model, keep=False, verbose=True, cutoff=3, seed=23000))
    solver = np.append(solver, JMT(model, samples=100000, seed=23000))
    solver = np.append(solver, MVA(model, method='exact'))
    solver = np.append(solver, DES(model, samples=100000, seed=23000))

    avg_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name()}')
        avg_table[s] = solver[s].avg_table()
