"""
Debug version of lqn_basic.py to trace the routing and visits
"""

from line_solver import *
import numpy as np

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = LayeredNetwork('test_LQN_4')

    P1 = Processor(model, 'P1', 2, SchedStrategy.PS)
    P2 = Processor(model, 'P2', 3, SchedStrategy.PS)

    T1 = Task(model, 'T1', 50, SchedStrategy.REF).on(P1).set_think_time(Exp(1 / 2))
    T2 = Task(model, 'T2', 50, SchedStrategy.FCFS).on(P1).set_think_time(Exp(1 / 3))

    T3 = Task(model, 'T3', 25, SchedStrategy.FCFS).on(P2).set_think_time(Exp(1 / 4))

    E1 = Entry(model, 'E1').on(T1)
    E2 = Entry(model, 'E2').on(T2)
    E3 = Entry(model, 'E3').on(T3)

    A1 = Activity(model, 'AS1', Exp(10)).on(T1).bound_to(E1).synch_call(E2, 1)
    A2 = Activity(model, 'AS2', Exp(20)).on(T2).bound_to(E2).synch_call(E3, 5).replies_to(E2)
    A3 = Activity(model, 'AS3', Exp(50)).on(T3).bound_to(E3).replies_to(E3)

    # Create solver
    solver = LN(model)

    # Print the layers before solving
    print("\n=== LAYER STRUCTURE ===")
    for i, layer in enumerate(solver.ensemble):
        if layer is not None:
            print(f"\nLayer {i}: {layer.name}")
            is_host = i in solver.hostLayerIndices
            print(f"  Type: {'HOST' if is_host else 'TASK'}")

            classes = layer.get_classes()
            nodes = layer.get_nodes()

            client = None
            server = None
            for n in nodes:
                if isinstance(n, Delay):
                    client = n
                elif isinstance(n, Queue):
                    server = n

            print(f"  Nodes: {[n.name for n in nodes]}")
            print(f"  Classes: {[c.name for c in classes]}")

            # Print service times
            print(f"  Client services:")
            for cls in classes:
                if client and hasattr(client, '_service_process') and cls in client._service_process:
                    proc = client._service_process[cls]
                    if proc is not None:
                        if hasattr(proc, 'getMean'):
                            mean = proc.getMean()
                        elif hasattr(proc, 'mean'):
                            mean = proc.mean
                        else:
                            mean = str(proc)
                        print(f"    {cls.name}: {mean}")

            print(f"  Server services:")
            for cls in classes:
                if server and hasattr(server, '_service_process') and cls in server._service_process:
                    proc = server._service_process[cls]
                    if proc is not None:
                        if hasattr(proc, 'getMean'):
                            mean = proc.getMean()
                        elif hasattr(proc, 'mean'):
                            mean = proc.mean
                        else:
                            mean = str(proc)
                        print(f"    {cls.name}: {mean}")

            # Print routing matrix
            sn = layer.get_struct()
            if sn is not None and sn.rt is not None:
                print(f"  Routing matrix shape: {sn.rt.shape}")
                print(f"  nclasses: {sn.nclasses}, nstateful: {sn.nstateful}")

                # Print non-zero routing entries
                rt = sn.rt
                for i in range(rt.shape[0]):
                    for j in range(rt.shape[1]):
                        if abs(rt[i, j]) > 1e-10:
                            from_station = i // sn.nclasses
                            from_class = i % sn.nclasses
                            to_station = j // sn.nclasses
                            to_class = j % sn.nclasses

                            from_class_name = classes[from_class].name if from_class < len(classes) else f"class_{from_class}"
                            to_class_name = classes[to_class].name if to_class < len(classes) else f"class_{to_class}"
                            from_node_name = nodes[from_station].name if from_station < len(nodes) else f"node_{from_station}"
                            to_node_name = nodes[to_station].name if to_station < len(nodes) else f"node_{to_station}"

                            print(f"    P[{from_class_name}, {to_class_name}, {from_node_name}, {to_node_name}] = {rt[i, j]:.4f}")

                # Print visits
                if hasattr(sn, 'visits') and sn.visits is not None:
                    print(f"  Visits:")
                    for c in range(len(sn.visits)):
                        if sn.visits[c] is not None:
                            print(f"    Chain {c}: {sn.visits[c]}")

    # Now run one iteration and check the results
    print("\n=== RUNNING ONE ITERATION ===")

    # Call the iterate method directly
    solver.pre(1)

    # Solve each layer
    layer_results = []
    for e, layer in enumerate(solver.ensemble):
        if layer is not None:
            layer_solver = solver.solvers[e]
            result = layer_solver.get_avg_table()
            layer_results.append(result)

            # Get the result structure
            if hasattr(layer_solver, '_result') and layer_solver._result is not None:
                res = layer_solver._result
                if 'QN' in res:
                    print(f"\nLayer {e} results:")
                    print(f"  QN: {res['QN']}")
                    print(f"  UN: {res['UN']}")
                    print(f"  RN: {res['RN']}")
                    print(f"  TN: {res['TN']}")

                    # Check stability - utilization should be < 1 per server
                    UN = res['UN']
                    nodes = layer.get_nodes()
                    for n_idx, node in enumerate(nodes):
                        if isinstance(node, Queue):
                            nservers = node.get_number_of_servers() if hasattr(node, 'get_number_of_servers') else 1
                            total_util = np.sum(UN[n_idx, :])
                            print(f"  Total utilization at {node.name}: {total_util:.4f} (servers: {nservers})")
                            if total_util > nservers:
                                print(f"    WARNING: Unstable! Utilization exceeds capacity!")
