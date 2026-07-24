"""
Debug version of lqn_basic.py to trace the overflow issue
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

    # Print map contents before iterating
    print("\n=== UPDATE MAPS ===")
    print(f"servt_classes_updmap:\n{solver.servt_classes_updmap}")
    print(f"thinkt_classes_updmap:\n{solver.thinkt_classes_updmap}")
    print(f"call_classes_updmap:\n{solver.call_classes_updmap}")
    print(f"hostLayerIndices: {solver.hostLayerIndices}")
    print(f"taskLayerIndices: {solver.taskLayerIndices}")
    print(f"idxhash: {solver.idxhash}")

    # Override iterate to add debugging
    original_post = solver.post

    def debug_pre(it):
        print(f"\n=== PRE-ITERATION {it} ===")
        if it == 1:
            # Print routing matrix for Layer 0
            layer0 = solver.ensemble[0]
            if layer0:
                sn = layer0.get_struct()
                if sn is not None:
                    print("  Layer 0 routing matrix shape:", sn.rt.shape if sn.rt is not None else "None")

        # Show Layer 0 service times before analysis
        if solver.ensemble and len(solver.ensemble) > 0:
            layer0 = solver.ensemble[0]
            if layer0:
                classes = layer0.get_classes()
                nodes = layer0.get_nodes()
                client = nodes[0] if len(nodes) > 0 else None
                server = nodes[1] if len(nodes) > 1 else None
                if client:
                    print(f"  Layer 0 CLIENT services:")
                    for cls in classes:
                        if hasattr(client, '_service_process') and cls in client._service_process:
                            proc = client._service_process[cls]
                            if proc is not None:
                                if hasattr(proc, 'getMean'):
                                    mean = proc.getMean()
                                elif hasattr(proc, 'mean'):
                                    mean = proc.mean
                                else:
                                    mean = str(proc)
                                print(f"    {cls.name}: {mean:.6e}" if isinstance(mean, float) else f"    {cls.name}: {mean}")
                if server and it >= 3:
                    print(f"  Layer 0 SERVER services (iteration {it}):")
                    for cls in classes:
                        if hasattr(server, '_service_process') and cls in server._service_process:
                            proc = server._service_process[cls]
                            if proc is not None:
                                if hasattr(proc, 'getMean'):
                                    mean = proc.getMean()
                                elif hasattr(proc, 'mean'):
                                    mean = proc.mean
                                else:
                                    mean = str(proc)
                                print(f"    {cls.name}: {mean:.6e}" if isinstance(mean, float) else f"    {cls.name}: {mean}")

                # Check the network struct for think time
                if it == 3:
                    sn = layer0.get_struct()
                    if sn is not None and hasattr(sn, 'rates') and sn.rates is not None:
                        print(f"  Layer 0 sn.rates (iteration {it}):")
                        print(f"    {sn.rates}")
                    if sn is not None and hasattr(sn, 'njobs') and sn.njobs is not None:
                        print(f"  Layer 0 sn.njobs (iteration {it}):")
                        print(f"    {sn.njobs}")

    solver.pre = debug_pre

    def debug_post(it):
        print(f"\n=== ITERATION {it} ===")

        lqn = solver.lqn

        # Show tput values (T1=3, T2=4, T3=5)
        print(f"tput[tidx]: T1={solver.tput[3]:.6e}, T2={solver.tput[4]:.6e}, T3={solver.tput[5]:.6e}")
        print(f"util[tidx]: T1={solver.util[3]:.6e}, T2={solver.util[4]:.6e}, T3={solver.util[5]:.6e}")
        print(f"thinkt[tidx]: T1={solver.thinkt[3]:.6e}, T2={solver.thinkt[4]:.6e}, T3={solver.thinkt[5]:.6e}")
        # Show thinktproc mean values
        for tidx in [3, 4, 5]:
            if solver.thinktproc[tidx] is not None:
                if hasattr(solver.thinktproc[tidx], 'getMean'):
                    mean = solver.thinktproc[tidx].getMean()
                elif hasattr(solver.thinktproc[tidx], 'mean'):
                    mean = solver.thinktproc[tidx].mean
                else:
                    mean = 'unknown'
                print(f"thinktproc[{tidx}] mean: {mean:.6e}" if isinstance(mean, float) else f"thinktproc[{tidx}]: {mean}")

        # Show servt values for activities
        print(f"servt[aidx]: AS1={solver.servt[9]:.6e}, AS2={solver.servt[10]:.6e}, AS3={solver.servt[11]:.6e}")

        # Show servtproc values for activities
        for aidx in [9, 10, 11]:
            if solver.servtproc[aidx] is not None:
                mean = solver.servtproc[aidx].getMean() if hasattr(solver.servtproc[aidx], 'getMean') else float(solver.servtproc[aidx].mean)
                print(f"servtproc[{aidx}] mean: {mean:.6e}")

        # Show servt for entries
        print(f"servt[eidx]: E1={solver.servt[6]:.6e}, E2={solver.servt[7]:.6e}, E3={solver.servt[8]:.6e}")

        # Show callservtproc
        for cidx in [1, 2]:
            if solver.callservtproc[cidx] is not None:
                if hasattr(solver.callservtproc[cidx], 'getMean'):
                    mean = solver.callservtproc[cidx].getMean()
                elif hasattr(solver.callservtproc[cidx], 'mean'):
                    mean = solver.callservtproc[cidx].mean
                else:
                    mean = 'unknown'
                print(f"callservtproc[{cidx}] mean: {mean:.6e}" if isinstance(mean, float) else f"callservtproc[{cidx}]: {mean}")

        # Show callservt values
        print(f"callservt: {solver.callservt[:5]}")
        # Show call_mean values
        for cidx in [1, 2]:
            call_mean = solver._get_call_mean(cidx)
            print(f"call_mean[{cidx}]: {call_mean}")

        # Debug: trace how callservt is computed
        lqn = solver.lqn
        print(f"  Debug callservt computation:")
        for c in range(len(solver.call_classes_updmap)):
            row = solver.call_classes_updmap[c]
            idx = int(row[0])
            cidx = int(row[1])
            nodeidx = int(row[2])
            classidx = int(row[3])
            if nodeidx > 1:
                layer_idx = int(solver.idxhash[idx])
                result = solver.results[-1][layer_idx] if len(solver.results) > 0 and layer_idx < len(solver.results[-1]) else None
                if result is not None and 'RN' in result:
                    RN = result['RN']
                    nodeidx_0 = nodeidx - 1
                    classidx_0 = classidx - 1
                    if nodeidx_0 < RN.shape[0] and classidx_0 < RN.shape[1]:
                        rn_val = RN[nodeidx_0, classidx_0]
                        call_mean = solver._get_call_mean(cidx)
                        computed = rn_val * call_mean
                        print(f"    cidx={cidx}: idx={idx} -> layer_idx={layer_idx}, nodeidx={nodeidx}, classidx={classidx}")
                        print(f"              RN[{nodeidx_0},{classidx_0}]={rn_val:.6e}, call_mean={call_mean}, computed={computed:.6e}")

        # Show njobs
        print(f"njobs[4, :] = {solver.njobs[4, :]}")
        print(f"njobs[5, :] = {solver.njobs[5, :]}")

        # Show layer results
        if solver.results and len(solver.results) > 0:
            for e, result in enumerate(solver.results[-1]):
                if result is not None and 'QN' in result:
                    print(f"  Layer {e}: QN shape={result['QN'].shape}")
                    print(f"  Layer {e}: TN server row={result['TN'][1, :] if result['TN'].shape[0] > 1 else 'N/A'}")
                    if 'UN' in result:
                        print(f"  Layer {e}: UN server row={result['UN'][1, :] if result['UN'].shape[0] > 1 else 'N/A'}")
                    if 'RN' in result:
                        print(f"  Layer {e}: RN server row={result['RN'][1, :] if result['RN'].shape[0] > 1 else 'N/A'}")

        # Call original post FIRST
        result = original_post(it)

        # NOW print updated values
        print(f"\n--- AFTER POST (iteration {it}) ---")
        print(f"callservt: {solver.callservt[:5]}")
        print(f"callservtproc[1] mean: {solver.callservtproc[1].getMean() if solver.callservtproc[1] else 0:.6e}")
        print(f"callservtproc[2] mean: {solver.callservtproc[2].getMean() if solver.callservtproc[2] else 0:.6e}")

        return result

    solver.post = debug_post

    try:
        avg_table = solver.get_avg_table()
        print("\n=== FINAL RESULTS ===")
        print(avg_table)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
