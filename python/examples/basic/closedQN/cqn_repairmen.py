"""
Repairman Problem - Classic Closed Queueing Network

This example models a repairman problem:
- 10 machines (jobs) in the system
- Delay station represents operational machines
- Queue station represents repair service
- Machines cycle between operation (delay) and repair (queue)
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    node = np.empty(2, dtype=object)
    node[0] = Delay(model, 'Delay')
    node[1] = Queue(model, 'Queue1', SchedStrategy.FCFS)

    jobclass = ClosedClass(model, 'Class1', 10, node[0], 0)

    node[0].set_service(jobclass, Exp.fit_mean(1.0))  # mean = 1
    node[1].set_service(jobclass, Exp.fit_mean(1.5))  # mean = 1.5

    P = model.init_routing_matrix()
    pmatrix = [[0.7, 0.3], [1.0, 0]]
    for i in range(len(node)):
        for j in range(len(node)):
            P.set(jobclass, jobclass, node[i], node[j], pmatrix[i][j])
    model.link(P)

    # Run multiple solvers
    solver = np.array([], dtype=object)
    solver = np.append(solver, CTMC(model))
    solver = np.append(solver, JMT(model, seed=23000))
    solver = np.append(solver, SSA(model, seed=23000, samples=5000))
    solver = np.append(solver, FLD(model))
    solver = np.append(solver, MVA(model))
    solver = np.append(solver, NC(model, method='exact'))
    solver = np.append(solver, MAM(model))
    solver = np.append(solver, DES(model, seed=23000, samples=5000))

    avg_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name()}')
        avg_table[s] = solver[s].avg_table()
