"""
LCFS Multiclass Closed Queueing Network

This example demonstrates:
- 2-station closed queueing network with LCFS and LCFS-PR scheduling
- 3 classes with 1 job each
- Validation between MVA and CTMC solvers

Reference:
    G. Casale, "A family of multiclass LCFS queueing networks with
    order-dependent product-form solutions", QUESTA 2026.
"""

from line_solver import *
import numpy as np


def cqn_lcfs_multiclass():
    # Service rates: mu(i,r) = service rate at station i for class r
    mu = 1./np.array([[1, 3, 5],    # Station 1 (Queue1)
                      [2, 4, 6]])   # Station 2 (Queue2)
    R = mu.shape[1]  # number of classes

    # Create the network model
    model = Network('LCFS Multiclass Model')

    # Create nodes
    node = np.empty(2, dtype=object)
    node[0] = Queue(model, 'Queue1', SchedStrategy.LCFS)
    node[1] = Queue(model, 'Queue2', SchedStrategy.LCFSPR)

    # Create job classes (one job per class)
    jobclass = np.empty(R, dtype=object)
    for r in range(R):
        jobclass[r] = ClosedClass(model, f'Class{r+1}', 1, node[0], 0)

    # Set service times (exponential distributions)
    for r in range(R):
        node[0].set_service(jobclass[r], Exp(mu[0, r]))
        node[1].set_service(jobclass[r], Exp(mu[1, r]))

    # Set up routing: jobs alternate between the two queues
    P = model.init_routing_matrix()
    for r in range(R):
        P[jobclass[r], jobclass[r]] = np.array([[0, 1], [1, 0]])  # Queue1 -> Queue2 -> Queue1
    model.link(P)

    return model


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = cqn_lcfs_multiclass()

    # Solve with MVA
    print('Solving LCFS+LCFS-PR network with MVA...')
    solver_mva = MVA(model, method='exact')
    print(f'\nSOLVER: {solver_mva.get_name()}')
    avg_table_mva = solver_mva.avg_table()
    print(avg_table_mva)

    # Solve with CTMC for validation
    print('\nSolving with CTMC for validation...')
    solver_ctmc = CTMC(model)
    print(f'\nSOLVER: {solver_ctmc.get_name()}')
    avg_table_ctmc = solver_ctmc.avg_table()
    print(avg_table_ctmc)

    # Compare results
    print('\n=== Comparison ===')
    print('MVA and CTMC results should match for this product-form network.')
