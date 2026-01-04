"""
Code generation for LINE Network models.

This module provides functionality to generate Python code from Network models,
enabling model reproduction and sharing.

Port from:
    - matlab/src/io/QN2MATLAB.m
    - matlab/src/io/LINE2MATLAB.m
    - matlab/src/io/QN2JAVA.m
"""

import sys
import numpy as np
from typing import Optional, Any, TextIO, Union
from io import StringIO

from ..mam import map_mean, map_scv


def qn2python(model: Any, model_name: str = 'my_model',
              file: Optional[Union[str, TextIO]] = None) -> Optional[str]:
    """
    Generate Python code that recreates a Network model.

    Converts a Network model to Python source code that, when executed,
    will create an equivalent model.

    Args:
        model: Network model or NetworkStruct
        model_name: Variable name for the model in generated code
        file: Output file path or file object. If None, returns string.

    Returns:
        Generated Python code as string if file is None, otherwise None.

    References:
        MATLAB: matlab/src/io/QN2MATLAB.m
    """
    # Get network structure
    if hasattr(model, 'getStruct'):
        sn = model.getStruct()
    else:
        sn = model

    # Determine output destination
    close_file = False
    if file is None:
        output = StringIO()
    elif isinstance(file, str):
        output = open(file, 'w')
        close_file = True
    else:
        output = file

    try:
        _generate_python_code(sn, model_name, output)

        if file is None:
            return output.getvalue()
        return None
    finally:
        if close_file:
            output.close()


def _generate_python_code(sn: Any, model_name: str, output: TextIO) -> None:
    """Generate Python code for the network structure."""

    # Header
    output.write("# LINE Network Model - Generated Python Code\n")
    output.write("# This code recreates the network model\n\n")

    # Imports
    output.write("from line_solver import Network, Source, Sink, Queue, Delay\n")
    output.write("from line_solver import OpenClass, ClosedClass\n")
    output.write("from line_solver import Exp, Erlang, HyperExp, APH\n")
    output.write("from line_solver.constants import SchedStrategy\n")
    output.write("import numpy as np\n\n")

    # Initialize model
    output.write(f"# Create network\n")
    output.write(f"{model_name} = Network('{model_name}')\n\n")

    # Get routing matrices
    rt = sn.rt if hasattr(sn, 'rt') else None
    rtnodes = sn.rtnodes if hasattr(sn, 'rtnodes') else None
    has_sink = False
    source_id = None
    PH = sn.proc if hasattr(sn, 'proc') else None

    # Block 1: Nodes
    output.write("# Block 1: Nodes\n")
    output.write("node = {}\n")

    for i in range(sn.nnodes):
        node_name = sn.nodenames[i] if hasattr(sn, 'nodenames') else f'Node{i}'
        node_type = sn.nodetype[i] if hasattr(sn, 'nodetype') else None

        if node_type is not None:
            type_name = node_type.name if hasattr(node_type, 'name') else str(node_type)

            if type_name == 'SOURCE':
                source_id = i
                output.write(f"node[{i}] = Source({model_name}, '{node_name}')\n")
                has_sink = True
            elif type_name == 'DELAY':
                output.write(f"node[{i}] = Delay({model_name}, '{node_name}')\n")
            elif type_name == 'QUEUE':
                ist = sn.nodeToStation[i] if hasattr(sn, 'nodeToStation') else i
                sched = sn.sched[ist] if hasattr(sn, 'sched') else None
                sched_name = sched.name if hasattr(sched, 'name') else 'FCFS'
                output.write(f"node[{i}] = Queue({model_name}, '{node_name}', SchedStrategy.{sched_name})\n")

                # Number of servers
                if hasattr(sn, 'nservers') and sn.nservers[ist] > 1:
                    output.write(f"node[{i}].setNumServers({int(sn.nservers[ist])})\n")
            elif type_name == 'ROUTER':
                output.write(f"node[{i}] = Router({model_name}, '{node_name}')\n")
            elif type_name == 'FORK':
                output.write(f"node[{i}] = Fork({model_name}, '{node_name}')\n")
            elif type_name == 'JOIN':
                # Find associated fork
                if hasattr(sn, 'fj') and sn.fj is not None:
                    fork_idx = np.where(sn.fj[:, i])[0]
                    if len(fork_idx) > 0:
                        output.write(f"node[{i}] = Join({model_name}, '{node_name}', node[{fork_idx[0]}])\n")
                    else:
                        output.write(f"node[{i}] = Join({model_name}, '{node_name}')\n")
                else:
                    output.write(f"node[{i}] = Join({model_name}, '{node_name}')\n")
            elif type_name == 'SINK':
                output.write(f"node[{i}] = Sink({model_name}, '{node_name}')\n")
            elif type_name == 'CLASSSWITCH':
                output.write(f"node[{i}] = Router({model_name}, '{node_name}')  # Class switching embedded in routing\n")

    output.write("\n")

    # Block 2: Classes
    output.write("# Block 2: Classes\n")
    output.write("jobclass = {}\n")

    for k in range(sn.nclasses):
        class_name = sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}'
        njobs = sn.njobs[k] if hasattr(sn, 'njobs') else 0
        priority = sn.classprio[k] if hasattr(sn, 'classprio') else 0

        if njobs > 0 or np.isinf(njobs):
            if np.isinf(njobs):
                output.write(f"jobclass[{k}] = OpenClass({model_name}, '{class_name}', {priority})\n")
            else:
                refstat = sn.refstat[k] if hasattr(sn, 'refstat') else 0
                ref_node = sn.stationToNode[refstat] if hasattr(sn, 'stationToNode') else refstat
                output.write(f"jobclass[{k}] = ClosedClass({model_name}, '{class_name}', {int(njobs)}, node[{ref_node}], {priority})\n")
        else:
            # Artificial class - find first station with non-null rate
            iref = 0
            if PH is not None:
                for ist in range(sn.nstations):
                    if ist < len(PH) and PH[ist] is not None:
                        if k < len(PH[ist]) and PH[ist][k] is not None:
                            if hasattr(PH[ist][k], '__getitem__') and len(PH[ist][k]) > 0:
                                if np.sum(np.abs(PH[ist][k][0])) > 0:
                                    iref = ist
                                    break

            if np.isinf(njobs):
                output.write(f"jobclass[{k}] = OpenClass({model_name}, '{class_name}', {priority})\n")
            else:
                output.write(f"jobclass[{k}] = ClosedClass({model_name}, '{class_name}', {int(njobs)}, node[{iref}], {priority})\n")

    output.write("\n")

    # Block 3: Arrival and Service Processes
    output.write("# Block 3: Arrival and Service Processes\n")

    coarse_tol = 1e-6

    for ist in range(sn.nstations):
        for k in range(sn.nclasses):
            node_type = None
            if hasattr(sn, 'nodetype') and hasattr(sn, 'stationToNode'):
                node_idx = sn.stationToNode[ist]
                node_type = sn.nodetype[node_idx]
                type_name = node_type.name if hasattr(node_type, 'name') else str(node_type)
            else:
                type_name = 'QUEUE'
                node_idx = ist

            # Skip Join nodes
            if type_name == 'JOIN':
                continue

            # Get process parameters
            if PH is not None and ist < len(PH) and PH[ist] is not None:
                if k < len(PH[ist]) and PH[ist][k] is not None:
                    try:
                        scv_ik = map_scv(PH[ist][k])
                        mean_ik = map_mean(PH[ist][k])
                    except Exception:
                        scv_ik = 1.0
                        mean_ik = 1.0
                else:
                    continue
            else:
                continue

            node_name = sn.nodenames[node_idx] if hasattr(sn, 'nodenames') else f'Node{node_idx}'
            class_name = sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}'

            sched = sn.sched[ist] if hasattr(sn, 'sched') else None
            sched_name = sched.name if hasattr(sched, 'name') else 'FCFS'

            if sched_name == 'EXT':
                # Arrival process
                if scv_ik >= 0.5:
                    if abs(scv_ik - 1.0) < coarse_tol:
                        if mean_ik < coarse_tol:
                            output.write(f"node[{node_idx}].setArrival(jobclass[{k}], Immediate())  # ({node_name},{class_name})\n")
                        else:
                            output.write(f"node[{node_idx}].setArrival(jobclass[{k}], Exp.fitMean({mean_ik:.6f}))  # ({node_name},{class_name})\n")
                    else:
                        output.write(f"node[{node_idx}].setArrival(jobclass[{k}], APH.fitMeanAndSCV({mean_ik:.6f}, {scv_ik:.6f}))  # ({node_name},{class_name})\n")
                else:
                    n_phases = max(1, round(1 / scv_ik))
                    if np.isnan(mean_ik):
                        output.write(f"node[{node_idx}].setArrival(jobclass[{k}], Disabled())  # ({node_name},{class_name})\n")
                    else:
                        output.write(f"node[{node_idx}].setArrival(jobclass[{k}], Erlang({n_phases / mean_ik:.6f}, {n_phases}))  # ({node_name},{class_name})\n")
            else:
                # Service process
                if scv_ik >= 0.5:
                    if abs(scv_ik - 1.0) < coarse_tol:
                        if mean_ik < coarse_tol:
                            output.write(f"node[{node_idx}].setService(jobclass[{k}], Immediate())  # ({node_name},{class_name})\n")
                        else:
                            output.write(f"node[{node_idx}].setService(jobclass[{k}], Exp.fitMean({mean_ik:.6f}))  # ({node_name},{class_name})\n")
                    else:
                        output.write(f"node[{node_idx}].setService(jobclass[{k}], APH.fitMeanAndSCV({mean_ik:.6f}, {scv_ik:.6f}))  # ({node_name},{class_name})\n")
                else:
                    n_phases = max(1, round(1 / scv_ik))
                    if np.isnan(mean_ik):
                        output.write(f"node[{node_idx}].setService(jobclass[{k}], Disabled())  # ({node_name},{class_name})\n")
                    else:
                        output.write(f"node[{node_idx}].setService(jobclass[{k}], Erlang({n_phases / mean_ik:.6f}, {n_phases}))  # ({node_name},{class_name})\n")

    output.write("\n")

    # Block 4: Routing
    output.write("# Block 4: Topology (Routing)\n")
    output.write(f"P = {model_name}.initRoutingMatrix()\n")

    if rtnodes is not None:
        for k in range(sn.nclasses):
            for c in range(sn.nclasses):
                for i in range(sn.nnodes):
                    for m in range(sn.nnodes):
                        idx_from = i * sn.nclasses + k
                        idx_to = m * sn.nclasses + c
                        if idx_from < rtnodes.shape[0] and idx_to < rtnodes.shape[1]:
                            prob = rtnodes[idx_from, idx_to]
                            if prob > 0:
                                # Skip if from sink
                                node_type_i = None
                                if hasattr(sn, 'nodetype'):
                                    node_type_i = sn.nodetype[i]
                                    type_name_i = node_type_i.name if hasattr(node_type_i, 'name') else str(node_type_i)
                                    if type_name_i == 'SINK':
                                        continue

                                node_name_i = sn.nodenames[i] if hasattr(sn, 'nodenames') else f'Node{i}'
                                node_name_m = sn.nodenames[m] if hasattr(sn, 'nodenames') else f'Node{m}'
                                class_name_k = sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}'
                                class_name_c = sn.classnames[c] if hasattr(sn, 'classnames') else f'Class{c}'

                                # Fork nodes always have probability 1.0
                                if hasattr(sn, 'nodetype'):
                                    type_name_i = sn.nodetype[i].name if hasattr(sn.nodetype[i], 'name') else str(sn.nodetype[i])
                                    if type_name_i == 'FORK':
                                        prob = 1.0

                                output.write(f"P[{k},{c}][{i},{m}] = {prob}  # ({node_name_i},{class_name_k}) -> ({node_name_m},{class_name_c})\n")

    output.write(f"{model_name}.link(P)\n")


def line2python(model: Any, filename: str) -> None:
    """
    Export a LINE Network model to Python code file.

    Args:
        model: Network model
        filename: Output Python file path

    References:
        MATLAB: matlab/src/io/LINE2MATLAB.m
    """
    model_name = getattr(model, 'name', 'model')
    qn2python(model, model_name, filename)
    print(f"Model exported to: {filename}")


def sn2python(sn: Any, model_name: str = 'model',
              file: Optional[Union[str, TextIO]] = None) -> Optional[str]:
    """
    Generate Python code from NetworkStruct.

    Alias for qn2python that emphasizes NetworkStruct input.

    Args:
        sn: NetworkStruct object
        model_name: Variable name for the model
        file: Output file path or file object

    Returns:
        Generated Python code as string if file is None
    """
    return qn2python(sn, model_name, file)


def qn2java(model: Any, model_name: str = 'myModel',
            file: Optional[Union[str, TextIO]] = None,
            headers: bool = True) -> Optional[str]:
    """
    Generate Java code that recreates a Network model.

    Converts a Network model to Java source code compatible with JLINE.

    Args:
        model: Network model or NetworkStruct
        model_name: Variable name for the model in generated code
        file: Output file path or file object. If None, returns string.
        headers: Whether to include function header/footer

    Returns:
        Generated Java code as string if file is None, otherwise None.

    References:
        MATLAB: matlab/src/io/QN2JAVA.m
    """
    # Get network structure
    if hasattr(model, 'getStruct'):
        sn = model.getStruct()
    else:
        sn = model

    # Determine output destination
    close_file = False
    if file is None:
        output = StringIO()
    elif isinstance(file, str):
        output = open(file, 'w')
        close_file = True
    else:
        output = file

    try:
        _generate_java_code(sn, model_name, output, headers)

        if file is None:
            return output.getvalue()
        return None
    finally:
        if close_file:
            output.close()


def _generate_java_code(sn: Any, model_name: str, output: TextIO, headers: bool) -> None:
    """Generate Java code for the network structure."""

    coarse_tol = 1e-6

    # Get routing matrices and process info
    rt = sn.rt if hasattr(sn, 'rt') else None
    rtnodes = sn.rtnodes if hasattr(sn, 'rtnodes') else None
    has_sink = False
    source_id = None
    PH = sn.proc if hasattr(sn, 'proc') else None

    # Header
    if headers:
        output.write(f"\tpublic static Network ex() {{\n")

    output.write(f'\t\tNetwork model = new Network("{model_name}");\n')
    output.write("\n\t\t// Block 1: nodes\t\t\t\n")

    # Block 1: Write nodes
    for i in range(sn.nnodes):
        node_name = sn.nodenames[i] if hasattr(sn, 'nodenames') else f'Node{i}'
        node_type = sn.nodetype[i] if hasattr(sn, 'nodetype') else None

        if node_type is not None:
            type_name = node_type.name if hasattr(node_type, 'name') else str(node_type)

            if type_name == 'SOURCE':
                source_id = i
                output.write(f'\t\tSource node{i+1} = new Source(model, "{node_name}");\n')
                has_sink = True
            elif type_name == 'DELAY':
                output.write(f'\t\tDelay node{i+1} = new Delay(model, "{node_name}");\n')
            elif type_name == 'QUEUE':
                ist = sn.nodeToStation[i] if hasattr(sn, 'nodeToStation') else i
                sched = sn.sched[ist] if hasattr(sn, 'sched') else None
                sched_name = sched.name if hasattr(sched, 'name') else 'FCFS'
                sched_prop = sched_name.upper()
                output.write(f'\t\tQueue node{i+1} = new Queue(model, "{node_name}", SchedStrategy.{sched_prop});\n')

                # Number of servers
                if hasattr(sn, 'nservers') and sn.nservers[ist] > 1:
                    if np.isinf(sn.nservers[ist]):
                        output.write(f'\t\tnode{i+1}.setNumberOfServers(Integer.MAX_VALUE);\n')
                    else:
                        output.write(f'\t\tnode{i+1}.setNumberOfServers({int(sn.nservers[ist])});\n')
            elif type_name == 'ROUTER':
                output.write(f'\t\tRouter node{i+1} = new Router(model, "{node_name}");\n')
            elif type_name == 'FORK':
                output.write(f'\t\tFork node{i+1} = new Fork(model, "{node_name}");\n')
            elif type_name == 'JOIN':
                # Find associated fork
                if hasattr(sn, 'fj') and sn.fj is not None:
                    fork_idx = np.where(sn.fj[:, i])[0]
                    if len(fork_idx) > 0:
                        output.write(f'\t\tJoin node{i+1} = new Join(model, "{node_name}", node{fork_idx[0]+1});\n')
                    else:
                        output.write(f'\t\tJoin node{i+1} = new Join(model, "{node_name}");\n')
                else:
                    output.write(f'\t\tJoin node{i+1} = new Join(model, "{node_name}");\n')
            elif type_name == 'SINK':
                output.write(f'\t\tSink node{i+1} = new Sink(model, "{node_name}");\n')
            elif type_name == 'CLASSSWITCH':
                output.write(f'\t\tRouter node{i+1} = new Router(model, "{node_name}"); // Dummy node, class switching is embedded in the routing matrix P \n')

    # Block 2: Write classes
    output.write("\n\t\t// Block 2: classes\n")

    for k in range(sn.nclasses):
        class_name = sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}'
        njobs = sn.njobs[k] if hasattr(sn, 'njobs') else 0
        priority = int(sn.classprio[k]) if hasattr(sn, 'classprio') else 0

        if njobs > 0 or np.isinf(njobs):
            if np.isinf(njobs):
                output.write(f'\t\tOpenClass jobclass{k+1} = new OpenClass(model, "{class_name}", {priority});\n')
            else:
                refstat = sn.refstat[k] if hasattr(sn, 'refstat') else 0
                ref_node = sn.stationToNode[refstat] if hasattr(sn, 'stationToNode') else refstat
                output.write(f'\t\tClosedClass jobclass{k+1} = new ClosedClass(model, "{class_name}", {int(njobs)}, node{ref_node+1}, {priority});\n')
        else:
            # Artificial class - find first station with non-null rate
            iref = 0
            if PH is not None:
                for ist in range(sn.nstations):
                    if ist < len(PH) and PH[ist] is not None:
                        if k < len(PH[ist]) and PH[ist][k] is not None:
                            if hasattr(PH[ist][k], '__getitem__') and len(PH[ist][k]) > 0:
                                if np.sum(np.abs(PH[ist][k][0])) > 0:
                                    iref = ist
                                    break

            if np.isinf(njobs):
                output.write(f'\t\tOpenClass jobclass{k+1} = new OpenClass(model, "{class_name}", {priority});\n')
            else:
                output.write(f'\t\tClosedClass jobclass{k+1} = new ClosedClass(model, "{class_name}", {int(njobs)}, node{iref+1}, {priority});\n')

    output.write("\t\t\n")

    # Block 3: Arrival and service processes
    for ist in range(sn.nstations):
        for k in range(sn.nclasses):
            node_type = None
            if hasattr(sn, 'nodetype') and hasattr(sn, 'stationToNode'):
                node_idx = sn.stationToNode[ist]
                node_type = sn.nodetype[node_idx]
                type_name = node_type.name if hasattr(node_type, 'name') else str(node_type)
            else:
                type_name = 'QUEUE'
                node_idx = ist

            # Skip Join nodes
            if type_name == 'JOIN':
                continue

            # Get process parameters
            if PH is not None and ist < len(PH) and PH[ist] is not None:
                if k < len(PH[ist]) and PH[ist][k] is not None:
                    try:
                        scv_ik = map_scv(PH[ist][k])
                        mean_ik = map_mean(PH[ist][k])
                    except Exception:
                        continue
                else:
                    continue
            else:
                continue

            node_name = sn.nodenames[node_idx] if hasattr(sn, 'nodenames') else f'Node{node_idx}'
            class_name = sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}'

            sched = sn.sched[ist] if hasattr(sn, 'sched') else None
            sched_name = sched.name if hasattr(sched, 'name') else 'FCFS'

            # Get schedparam if available
            schedparam = 1.0
            if hasattr(sn, 'schedparam') and sn.schedparam is not None:
                if ist < sn.schedparam.shape[0] and k < sn.schedparam.shape[1]:
                    schedparam = sn.schedparam[ist, k]

            if sched_name == 'EXT':
                # Arrival process
                if scv_ik >= 0.5:
                    if abs(scv_ik - 1.0) < coarse_tol:
                        if mean_ik < coarse_tol:
                            output.write(f'\t\tnode{node_idx+1}.setArrival(jobclass{k+1}, Immediate.getInstance()); // ({node_name},{class_name})\n')
                        else:
                            output.write(f'\t\tnode{node_idx+1}.setArrival(jobclass{k+1}, Exp.fitMean({mean_ik})); // ({node_name},{class_name})\n')
                    else:
                        output.write(f'\t\tnode{node_idx+1}.setArrival(jobclass{k+1}, APH.fitMeanAndSCV({mean_ik},{scv_ik})); // ({node_name},{class_name})\n')
                else:
                    n_phases = max(1, round(1 / scv_ik))
                    if np.isnan(PH[ist][k][0]) if hasattr(PH[ist][k], '__getitem__') else np.isnan(mean_ik):
                        output.write(f'\t\tnode{node_idx+1}.setArrival(jobclass{k+1}, Disabled.getInstance()); // ({node_name},{class_name})\n')
                    else:
                        output.write(f'\t\tnode{node_idx+1}.setArrival(jobclass{k+1}, Erlang({n_phases/mean_ik},{n_phases})); // ({node_name},{class_name})\n')
            else:
                # Service process
                if scv_ik >= 0.5:
                    if abs(scv_ik - 1.0) < coarse_tol:
                        if mean_ik < coarse_tol:
                            output.write(f'\t\tnode{node_idx+1}.setService(jobclass{k+1}, Immediate.getInstance()); // ({node_name},{class_name})\n')
                        else:
                            if schedparam != 1:
                                output.write(f'\t\tnode{node_idx+1}.setService(jobclass{k+1}, Exp.fitMean({mean_ik}), {schedparam}); // ({node_name},{class_name})\n')
                            else:
                                output.write(f'\t\tnode{node_idx+1}.setService(jobclass{k+1}, Exp.fitMean({mean_ik})); // ({node_name},{class_name})\n')
                    else:
                        if schedparam != 1:
                            output.write(f'\t\tnode{node_idx+1}.setService(jobclass{k+1}, APH.fitMeanAndSCV({mean_ik},{scv_ik}), {schedparam}); // ({node_name},{class_name})\n')
                        else:
                            output.write(f'\t\tnode{node_idx+1}.setService(jobclass{k+1}, APH.fitMeanAndSCV({mean_ik},{scv_ik})); // ({node_name},{class_name})\n')
                else:
                    n_phases = max(1, round(1 / scv_ik))
                    if np.isnan(PH[ist][k][0]) if hasattr(PH[ist][k], '__getitem__') else np.isnan(mean_ik):
                        output.write(f'\t\tnode{node_idx+1}.setService(jobclass{k+1}, Disabled.getInstance()); // ({node_name},{class_name})\n')
                    else:
                        if schedparam != 1:
                            output.write(f'\t\tnode{node_idx+1}.setService(jobclass{k+1}, Erlang({n_phases/mean_ik},{n_phases}), {schedparam}); // ({node_name},{class_name})\n')
                        else:
                            output.write(f'\t\tnode{node_idx+1}.setService(jobclass{k+1}, Erlang({n_phases/mean_ik},{n_phases})); // ({node_name},{class_name})\n')

    # Block 4: Topology
    output.write("\n\t\t// Block 3: topology")

    # Handle sink routing for open classes
    if has_sink and source_id is not None and rt is not None:
        for k in range(sn.nclasses):
            if np.isinf(sn.njobs[k]):  # open class
                for ist in range(sn.nstations):
                    # Redirect transitions to source to sink instead
                    idx_from = ist * sn.nclasses + k
                    idx_to_source = source_id * sn.nclasses + k
                    if idx_from < rt.shape[0] and idx_to_source < rt.shape[1]:
                        rt[sn.nstations * sn.nclasses + k, idx_to_source] = rt[idx_from, idx_to_source]
                        rt[idx_from, idx_to_source] = 0

    output.write("\t\n")
    output.write("\t\tRoutingMatrix routingMatrix = model.initRoutingMatrix(); \n")
    output.write("\t\n")

    # Write routing probabilities
    if rtnodes is not None:
        for k in range(sn.nclasses):
            for c in range(sn.nclasses):
                for i in range(sn.nnodes):
                    for m in range(sn.nnodes):
                        idx_from = i * sn.nclasses + k
                        idx_to = m * sn.nclasses + c
                        if idx_from < rtnodes.shape[0] and idx_to < rtnodes.shape[1]:
                            prob = rtnodes[idx_from, idx_to]
                            if prob > 0:
                                # Skip if from sink
                                if hasattr(sn, 'nodetype'):
                                    type_name_i = sn.nodetype[i].name if hasattr(sn.nodetype[i], 'name') else str(sn.nodetype[i])
                                    if type_name_i == 'SINK':
                                        continue
                                    # Fork nodes use fanout value
                                    if type_name_i == 'FORK':
                                        fanout = 1
                                        if hasattr(sn, 'nodeparam') and sn.nodeparam[i] is not None:
                                            if hasattr(sn.nodeparam[i], 'fanOut'):
                                                fanout = sn.nodeparam[i].fanOut
                                        prob = fanout

                                node_name_i = sn.nodenames[i] if hasattr(sn, 'nodenames') else f'Node{i}'
                                node_name_m = sn.nodenames[m] if hasattr(sn, 'nodenames') else f'Node{m}'
                                class_name_k = sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}'
                                class_name_c = sn.classnames[c] if hasattr(sn, 'classnames') else f'Class{c}'

                                output.write(f'\t\troutingMatrix.set(jobclass{k+1}, jobclass{c+1}, node{i+1}, node{m+1}, {prob}); // ({node_name_i},{class_name_k}) -> ({node_name_m},{class_name_c})\n')

    output.write("\n\t\tmodel.link(routingMatrix);\n\n")

    if headers:
        output.write("\t\treturn model;\n")
        output.write("\t}\n")


def lqn2java(model: Any, model_name: str = 'myLayeredModel',
             file: Optional[Union[str, TextIO]] = None) -> Optional[str]:
    """
    Generate Java code that recreates a LayeredNetwork model.

    Converts a LayeredNetwork model to Java source code compatible with JLINE.

    Args:
        model: LayeredNetwork model
        model_name: Variable name for the model in generated code
        file: Output file path or file object. If None, returns string.

    Returns:
        Generated Java code as string if file is None, otherwise None.

    References:
        MATLAB: matlab/src/io/LQN2JAVA.m
    """
    # Get network structure
    if hasattr(model, 'getStruct'):
        sn = model.getStruct()
    else:
        sn = model

    # Determine output destination
    close_file = False
    if file is None:
        output = StringIO()
    elif isinstance(file, str):
        output = open(file, 'w')
        close_file = True
    else:
        output = file

    try:
        _generate_lqn_java_code(sn, model, model_name, output)

        if file is None:
            return output.getvalue()
        return None
    finally:
        if close_file:
            output.close()


def _generate_lqn_java_code(sn: Any, model: Any, model_name: str, output: TextIO) -> None:
    """Generate Java code for a LayeredNetwork structure."""

    # Package and imports
    output.write("package jline.examples;\n\n")
    output.write("import jline.lang.*;\n")
    output.write("import jline.lang.constant.*;\n")
    output.write("import jline.lang.processes.*;\n")
    output.write("import jline.solvers.ln.SolverLN;\n\n")

    clean_name = model_name.replace(' ', '').upper()
    output.write(f"public class TestSolver{clean_name} {{\n\n")
    output.write("\tpublic static void main(String[] args) throws Exception{\n\n")
    output.write(f'\tLayeredNetwork model = new LayeredNetwork("{model_name}");\n')
    output.write("\n")

    # Host processors
    nhosts = sn.nhosts if hasattr(sn, 'nhosts') else 0
    for h in range(nhosts):
        name = sn.names[h] if hasattr(sn, 'names') else f'Host{h}'
        mult = sn.mult[h] if hasattr(sn, 'mult') else 1
        sched = sn.sched[h] if hasattr(sn, 'sched') else None
        sched_str = _get_sched_feature(sched) if sched else 'SchedStrategy.FCFS'

        if np.isinf(mult):
            output.write(f'\tProcessor P{h+1} = new Processor(model, "{name}", Integer.MAX_VALUE, {sched_str});\n')
        else:
            output.write(f'\tProcessor P{h+1} = new Processor(model, "{name}", {int(mult)}, {sched_str});\n')

        # Replication
        repl = sn.repl[h] if hasattr(sn, 'repl') else 1
        if repl != 1:
            output.write(f'\tP{h+1}.setReplication({int(repl)});\n')

    output.write("\n")

    # Tasks
    ntasks = sn.ntasks if hasattr(sn, 'ntasks') else 0
    tshift = sn.tshift if hasattr(sn, 'tshift') else nhosts

    for t in range(ntasks):
        tidx = tshift + t
        name = sn.names[tidx] if hasattr(sn, 'names') else f'Task{t}'
        mult = sn.mult[tidx] if hasattr(sn, 'mult') else 1
        sched = sn.sched[tidx] if hasattr(sn, 'sched') else None
        sched_str = _get_sched_feature(sched) if sched else 'SchedStrategy.FCFS'
        parent = sn.parent[tidx] if hasattr(sn, 'parent') else 1

        if np.isinf(mult):
            output.write(f'\tTask T{t+1} = new Task(model, "{name}", Integer.MAX_VALUE, {sched_str}).on(P{parent});\n')
        else:
            output.write(f'\tTask T{t+1} = new Task(model, "{name}", {int(mult)}, {sched_str}).on(P{parent});\n')

        # Replication
        repl = sn.repl[tidx] if hasattr(sn, 'repl') else 1
        if repl != 1:
            output.write(f'\tT{t+1}.setReplication({int(repl)});\n')

        # Think time
        if hasattr(sn, 'think') and sn.think is not None and tidx < len(sn.think):
            think = sn.think[tidx]
            think_type = sn.think_type[tidx] if hasattr(sn, 'think_type') else None
            if think is not None and think_type is not None:
                think_mean = sn.think_mean[tidx] if hasattr(sn, 'think_mean') else 0
                think_scv = sn.think_scv[tidx] if hasattr(sn, 'think_scv') else 1

                type_name = think_type.name if hasattr(think_type, 'name') else str(think_type)
                if type_name == 'IMMEDIATE':
                    output.write(f'\tT{t+1}.setThinkTime(new Immediate());\n')
                elif type_name == 'EXP':
                    output.write(f'\tT{t+1}.setThinkTime(new Exp({1/think_mean}));\n')
                elif type_name in ['ERLANG', 'HYPEREXP', 'COXIAN', 'APH']:
                    output.write(f'\tT{t+1}.setThinkTime({type_name}.fitMeanAndSCV({think_mean},{think_scv}));\n')

    output.write("\n")

    # Entries
    nentries = sn.nentries if hasattr(sn, 'nentries') else 0
    eshift = sn.eshift if hasattr(sn, 'eshift') else tshift + ntasks

    for e in range(nentries):
        eidx = eshift + e
        name = sn.names[eidx] if hasattr(sn, 'names') else f'Entry{e}'
        parent = sn.parent[eidx] if hasattr(sn, 'parent') else 1
        output.write(f'\tEntry E{e+1} = new Entry(model, "{name}").on(T{parent - tshift});\n')

    output.write("\n")

    # Activities
    nacts = sn.nacts if hasattr(sn, 'nacts') else 0
    ashift = sn.ashift if hasattr(sn, 'ashift') else eshift + nentries

    for a in range(nacts):
        aidx = ashift + a
        name = sn.names[aidx] if hasattr(sn, 'names') else f'Activity{a}'
        parent_tidx = sn.parent[aidx] if hasattr(sn, 'parent') else tshift + 1
        on_task = parent_tidx - tshift

        # Get host demand distribution
        hostdem_type = sn.hostdem_type[aidx] if hasattr(sn, 'hostdem_type') else None
        hostdem_mean = sn.hostdem_mean[aidx] if hasattr(sn, 'hostdem_mean') else 1.0
        hostdem_scv = sn.hostdem_scv[aidx] if hasattr(sn, 'hostdem_scv') else 1.0

        type_name = hostdem_type.name if hasattr(hostdem_type, 'name') else 'EXP'

        if type_name == 'IMMEDIATE':
            output.write(f'\tActivity A{a+1} = new Activity(model, "{name}", new Immediate()).on(T{on_task});')
        elif type_name == 'EXP':
            output.write(f'\tActivity A{a+1} = new Activity(model, "{name}", new Exp({1/hostdem_mean})).on(T{on_task});')
        elif type_name in ['ERLANG', 'HYPEREXP', 'COXIAN', 'APH']:
            output.write(f'\tActivity A{a+1} = new Activity(model, "{name}", {type_name}.fitMeanAndSCV({hostdem_mean},{hostdem_scv})).on(T{on_task});')
        else:
            output.write(f'\tActivity A{a+1} = new Activity(model, "{name}", new Exp({1/hostdem_mean})).on(T{on_task});')

        # Bound to entry
        if hasattr(sn, 'graph') and sn.graph is not None:
            for e in range(nentries):
                eidx = eshift + e
                if eidx < sn.graph.shape[0] and aidx < sn.graph.shape[1]:
                    if sn.graph[eidx, aidx] > 0:
                        output.write(f' A{a+1}.boundTo(E{e+1});')
                        break

        # Calls
        if hasattr(sn, 'callpair') and sn.callpair is not None:
            for c in range(len(sn.callpair)):
                if sn.callpair[c, 0] == aidx:
                    target_entry = sn.callpair[c, 1] - eshift
                    call_mean = sn.callproc_mean[c] if hasattr(sn, 'callproc_mean') else 1
                    call_type = sn.calltype[c] if hasattr(sn, 'calltype') else None
                    call_type_name = call_type.name if hasattr(call_type, 'name') else 'SYNC'
                    if call_type_name == 'SYNC':
                        output.write(f' A{a+1}.synchCall(E{target_entry},{call_mean});')
                    elif call_type_name == 'ASYNC':
                        output.write(f' A{a+1}.asynchCall(E{target_entry},{call_mean});')

        # Replies to
        if hasattr(sn, 'replygraph') and sn.replygraph is not None:
            if a < sn.replygraph.shape[0]:
                for e in range(nentries):
                    if e < sn.replygraph.shape[1] and sn.replygraph[a, e] > 0:
                        output.write(f' A{a+1}.repliesTo(E{e+1});')
                        break

        output.write("\n")

    output.write("\n")

    # Model solution
    output.write("\t// Model solution \n")
    output.write("\tSolverLN solver = new SolverLN(model);\n")
    output.write("\tsolver.getEnsembleAvg();\n")
    output.write("\t}\n}\n")


def _get_sched_feature(sched: Any) -> str:
    """Convert scheduling strategy to Java feature string."""
    if sched is None:
        return 'SchedStrategy.FCFS'
    sched_name = sched.name if hasattr(sched, 'name') else str(sched)
    return f'SchedStrategy.{sched_name.upper()}'


def line2java(model: Any, filename: Optional[str] = None) -> Optional[str]:
    """
    Export a LINE model to Java code file.

    Dispatches to qn2java or lqn2java based on model type.

    Args:
        model: Network or LayeredNetwork model
        filename: Optional output file path. If None, prints to stdout/returns string.

    Returns:
        Generated Java code as string if filename is None

    References:
        MATLAB: matlab/src/io/LINE2JAVA.m
    """
    model_name = getattr(model, 'name', 'model')

    # Check model type
    is_layered = hasattr(model, 'tasks') or (
        hasattr(model, '__class__') and 'Layered' in model.__class__.__name__
    )

    if is_layered:
        return lqn2java(model, model_name, filename)
    else:
        return qn2java(model, model_name, filename)


__all__ = [
    'qn2python',
    'line2python',
    'sn2python',
    'qn2java',
    'lqn2java',
    'line2java',
]
