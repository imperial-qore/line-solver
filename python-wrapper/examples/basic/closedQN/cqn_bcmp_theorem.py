"""
BCMP Theorem - Product-Form Scheduling Policies

This example demonstrates:
- Comparison of product-form scheduling policies (PS, FCFS, LCFS-PR)
- All three models have identical performance (BCMP theorem)
- 2 closed classes with different service distributions
- Erlang and HyperExp distributions at delay
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    c = 1

    # PS scheduling model
    ps_model = Network('PS scheduling model')

    node = np.empty(2, dtype=object)
    node[0] = Delay(ps_model, 'Delay')
    node[1] = Queue(ps_model, 'Queue1', SchedStrategy.PS)
    node[1].set_num_servers(c)

    jobclass = np.empty(2, dtype=object)
    jobclass[0] = ClosedClass(ps_model, 'Class1', 2, node[0], 0)
    jobclass[1] = ClosedClass(ps_model, 'Class2', 2, node[0], 0)

    node[0].set_service(jobclass[0], Erlang(3, 2))
    node[0].set_service(jobclass[1], HyperExp(0.5, 3.0, 10.0))

    node[1].set_service(jobclass[0], Exp(1))
    node[1].set_service(jobclass[1], Exp(1))

    P = ps_model.init_routing_matrix()
    routing = Network.serial_routing(node)
    P.set(jobclass[0], jobclass[0], node[0], node[1], routing[0][1])
    P.set(jobclass[0], jobclass[0], node[1], node[0], routing[1][0])
    P.set(jobclass[1], jobclass[1], node[0], node[1], routing[0][1])
    P.set(jobclass[1], jobclass[1], node[1], node[0], routing[1][0])
    ps_model.link(P)

    # FCFS scheduling model
    fcfs_model = Network('FCFS scheduling model')

    node = np.empty(2, dtype=object)
    node[0] = Delay(fcfs_model, 'Delay')
    node[1] = Queue(fcfs_model, 'Queue1', SchedStrategy.FCFS)
    node[1].set_num_servers(c)

    jobclass = np.empty(2, dtype=object)
    jobclass[0] = ClosedClass(fcfs_model, 'Class1', 2, node[0], 0)
    jobclass[1] = ClosedClass(fcfs_model, 'Class2', 2, node[0], 0)

    node[0].set_service(jobclass[0], Erlang(3, 2))
    node[0].set_service(jobclass[1], HyperExp(0.5, 3.0, 10.0))

    node[1].set_service(jobclass[0], Exp(1))
    node[1].set_service(jobclass[1], Exp(1))

    P = fcfs_model.init_routing_matrix()
    routing = Network.serial_routing(node)
    P.set(jobclass[0], jobclass[0], node[0], node[1], routing[0][1])
    P.set(jobclass[0], jobclass[0], node[1], node[0], routing[1][0])
    P.set(jobclass[1], jobclass[1], node[0], node[1], routing[0][1])
    P.set(jobclass[1], jobclass[1], node[1], node[0], routing[1][0])
    fcfs_model.link(P)

    # LCFS-PR scheduling model
    lcfspr_model = Network('LCFS-PR scheduling model')

    node = np.empty(2, dtype=object)
    node[0] = Delay(lcfspr_model, 'Delay')
    node[1] = Queue(lcfspr_model, 'Queue1', SchedStrategy.LCFSPR)
    node[1].set_num_servers(c)

    jobclass = np.empty(2, dtype=object)
    jobclass[0] = ClosedClass(lcfspr_model, 'Class1', 2, node[0], 0)
    jobclass[1] = ClosedClass(lcfspr_model, 'Class2', 2, node[0], 0)

    node[0].set_service(jobclass[0], Erlang(3, 2))
    node[0].set_service(jobclass[1], HyperExp(0.5, 3.0, 10.0))

    node[1].set_service(jobclass[0], Exp(1))
    node[1].set_service(jobclass[1], Exp(1))

    P = lcfspr_model.init_routing_matrix()
    routing = Network.serial_routing(node)
    P.set(jobclass[0], jobclass[0], node[0], node[1], routing[0][1])
    P.set(jobclass[0], jobclass[0], node[1], node[0], routing[1][0])
    P.set(jobclass[1], jobclass[1], node[0], node[1], routing[0][1])
    P.set(jobclass[1], jobclass[1], node[1], node[0], routing[1][0])
    lcfspr_model.link(P)

    # Solve all three models
    solver = np.array([], dtype=object)
    solver = np.append(solver, CTMC(ps_model))
    solver = np.append(solver, CTMC(fcfs_model))
    solver = np.append(solver, CTMC(lcfspr_model))

    for s in range(len(solver)):
        print(f'\nMODEL: {solver[s].model.get_name()}')
        avg_table = solver[s].avg_table()
