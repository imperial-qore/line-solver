"""
Multiserver Multiclass Closed Network

This example demonstrates:
- 4 closed classes with different populations
- Multiserver queues (3 servers each)
- Different service distributions (Exp, Erlang)
- Disabled service for specific class-station pairs
- Complex routing between classes
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    node = np.empty(3, dtype=object)
    node[0] = Delay(model, 'Delay')
    node[1] = Queue(model, 'Queue1', SchedStrategy.FCFS)
    node[2] = Queue(model, 'Queue2', SchedStrategy.FCFS)
    node[1].set_num_servers(3)
    node[2].set_num_servers(3)

    jobclass = np.empty(4, dtype=object)
    jobclass[0] = ClosedClass(model, 'Class1', 2, node[0], 0)
    jobclass[1] = ClosedClass(model, 'Class2', 2, node[0], 0)
    jobclass[2] = ClosedClass(model, 'Class3', 2, node[0], 0)
    jobclass[3] = ClosedClass(model, 'Class4', 1, node[0], 0)

    # Set service times
    for r in range(4):
        if r < 3:
            node[0].set_service(jobclass[r], Exp(1) if r < 2 else Exp(10))
        else:
            node[0].set_service(jobclass[r], Exp(1))

    node[1].set_service(jobclass[0], Exp(1))
    node[1].set_service(jobclass[1], Erlang(1, 2))
    node[1].set_service(jobclass[2], Exp(10))
    node[1].set_service(jobclass[3], Exp(1))

    node[2].set_service(jobclass[0], Disabled())
    node[2].set_service(jobclass[1], Disabled())
    node[2].set_service(jobclass[2], Erlang(1, 2))
    node[2].set_service(jobclass[3], Exp(1))

    # Set routing
    M = model.get_number_of_stations()
    K = model.get_number_of_classes()

    my_p = model.init_routing_matrix()

    # Class 1 routing
    pmatrix_1_1 = [[0, 0.5, 0], [1, 0, 0], [1, 0, 0]]
    pmatrix_1_2 = [[0, 0.5, 0], [0, 0, 0], [0, 0, 0]]

    for i in range(3):
        for j in range(3):
            my_p.set(jobclass[0], jobclass[0], node[i], node[j], pmatrix_1_1[i][j])
            my_p.set(jobclass[0], jobclass[1], node[i], node[j], pmatrix_1_2[i][j])

    # Class 2 routing
    pmatrix_2_2 = [[0, 0, 0], [1, 0, 0], [1, 0, 0]]
    pmatrix_2_1 = [[0, 1, 0], [0, 0, 0], [0, 0, 0]]

    for i in range(3):
        for j in range(3):
            my_p.set(jobclass[1], jobclass[1], node[i], node[j], pmatrix_2_2[i][j])
            my_p.set(jobclass[1], jobclass[0], node[i], node[j], pmatrix_2_1[i][j])

    # Class 3 routing
    pmatrix_3_3 = [[0, 0.25, 0.25], [1, 0, 0], [1, 0, 0]]
    pmatrix_3_4 = [[0, 0.5, 0], [0, 0, 0], [0, 0, 0]]

    for i in range(3):
        for j in range(3):
            my_p.set(jobclass[2], jobclass[2], node[i], node[j], pmatrix_3_3[i][j])
            my_p.set(jobclass[2], jobclass[3], node[i], node[j], pmatrix_3_4[i][j])

    # Class 4 routing
    pmatrix_4_4 = [[0, 0, 0], [1, 0, 0], [1, 0, 0]]
    pmatrix_4_3 = [[0, 1, 0], [0, 0, 0], [0, 0, 0]]

    for i in range(3):
        for j in range(3):
            my_p.set(jobclass[3], jobclass[3], node[i], node[j], pmatrix_4_4[i][j])
            my_p.set(jobclass[3], jobclass[2], node[i], node[j], pmatrix_4_3[i][j])

    model.link(my_p)

    # Get state space examples
    space_running = State.from_marginal_and_running(model, node[1], [2, 1, 1, 1], [2, 1, 0, 0])
    print(f'State (Running): {space_running}')

    space_started = State.from_marginal_and_started(model, node[1], [2, 1, 1, 1], [2, 1, 0, 0])
    print(f'State (Started): {space_started}')

    space = State.from_marginal(model, node[1], [2, 1, 1, 1])
    print(f'State (Marginal): {space}')

    # Run solver
    print('\nSOLVER: MVA')
    solver = SolverMVA(model)
    avg_table = solver.getAvgTable()
