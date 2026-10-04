"""
Code generation for LINE Network models.

This module provides functionality to generate Python code from Network models,
enabling model reproduction and sharing.

Port from:
    - matlab/src/io/QN2MATLAB.m
    - matlab/src/io/LINE2MATLAB.m
    - matlab/src/io/QN2JAVA.m
    - matlab/src/io/LQN2JAVA.m (lqn2python writes the same layered model as Python)
"""

import sys
import numpy as np
from typing import Optional, Any, TextIO, Union
from io import StringIO

from ..mam import map_mean, map_scv


def _fork_of(sn, join_idx: int, join_name: str) -> int:
    """The node index of the Fork a Join closes, read from sn.fj.

    Raises rather than emitting an unpaired Join: the pairing is a declaration
    carried by the Join, not a derivation from the routing, and a model whose
    Join names no Fork is refused at struct refresh in all four codebases.
    """
    fj = getattr(sn, 'fj', None)
    if fj is not None:
        fork_idx = np.where(np.asarray(fj)[:, join_idx])[0]
        if fork_idx.size > 0:
            return int(fork_idx[0])
    raise ValueError("Join '%s' closes no Fork: the model cannot be written as "
                     "source." % (join_name,))


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
    sn, net = _struct_and_model(model)
    return _emit(lambda _m, name, out: _generate_python_code(sn, name, out, net), model, model_name, file)


def _num(x) -> str:
    """Shortest literal that reads back to the same double, valid in both Python and MATLAB source."""
    x = float(x)
    if np.isinf(x):
        return 'Inf' if x > 0 else '-Inf'
    if np.isnan(x):
        return 'NaN'
    return repr(x)


def _py_num(x) -> str:
    """Python spelling of _num: infinities and NaN go through numpy."""
    s = _num(x)
    return {'Inf': 'np.inf', '-Inf': '-np.inf', 'NaN': 'np.nan'}.get(s, s)


def _count(x, python: bool) -> str:
    """A server count or multiplicity: an integer literal when finite, infinity otherwise."""
    x = float(x)
    if np.isinf(x):
        return 'np.inf' if python else 'Inf'
    return '%d' % int(x) if x == int(x) else repr(x)


def _type_name(value) -> str:
    return value.name if hasattr(value, 'name') else str(value)


def _sched_name(sched) -> str:
    """Name of a scheduling strategy stored either as the enum or as its integer id."""
    if hasattr(sched, 'name'):
        return sched.name
    from ...constants import SchedStrategy
    return SchedStrategy(int(sched)).name


def _qn_process_spec(sn, model, ist: int, k: int):
    """The arrival or service process QN2MATLAB writes for station ist and class k.

    Returns (verb, kind, args) with verb 'setArrival' at an EXT station and 'setService' elsewhere,
    kind one of 'Replayer', 'Trace', 'Immediate', 'Exp', 'APH', 'Erlang', 'Disabled', and args the
    numeric (or file) arguments. Returns None for a Join station, which carries no process. The
    fitting rules are those of matlab/src/io/QN2MATLAB.m: SCV >= 0.5 gives Exp (SCV == 1) or APH,
    a smaller SCV an Erlang of round(1/SCV) phases, and an undefined process Disabled.
    """
    from ...constants import GlobalConstants
    node_idx = int(sn.stationToNode[ist])
    if _type_name(sn.nodetype[node_idx]) == 'JOIN':
        return None
    verb = 'setArrival' if _sched_name(sn.sched[ist]) == 'EXT' else 'setService'
    if model is not None and hasattr(model, 'get_stations'):
        station = model.get_stations()[ist]
        jobclass = model.get_classes()[k]
        getter = station.get_arrival_process if verb == 'setArrival' and hasattr(station, 'get_arrival_process') \
            else getattr(station, 'get_service_process', None)
        proc = getter(jobclass) if getter is not None else None
        pname = type(proc).__name__
        if pname in ('Replayer', 'Trace'):
            return verb, pname, (proc,)
    ph = sn.proc[ist][k] if sn.proc is not None and ist < len(sn.proc) and k < len(sn.proc[ist]) else None
    if ph is None or np.any(np.isnan(np.asarray(ph[0], dtype=float))):
        return verb, 'Disabled', ()
    scv = float(map_scv(ph))
    mean = float(map_mean(ph))
    if scv >= 0.5:
        if scv == 1:
            if mean < GlobalConstants.CoarseTol:
                return verb, 'Immediate', ()
            return verb, 'Exp', (mean,)
        return verb, 'APH', (mean, scv)
    nphases = max(1, int(round(1.0 / scv)))
    return verb, 'Erlang', (nphases / mean, nphases)


def _qn_class_ref_node(sn, k: int) -> int:
    """Node index of the reference station of class k, as QN2MATLAB picks it.

    A class with jobs uses its own reference station. A zero-population class takes its chain's
    reference station (the refstat of the chain's populated closed classes, else its own refstat),
    since refreshRoutingMatrix rejects a chain whose classes name different reference stations;
    failing that, the first station with a non-null rate for it. Port of zeroPopRefNode in
    matlab/src/io/QN2MATLAB.m.
    """
    njobs = sn.njobs[k]
    if njobs > 0:
        return int(sn.stationToNode[int(sn.refstat[k])])
    ist = -1
    chains = np.asarray(getattr(sn, 'chains', np.zeros((0, sn.nclasses))))
    refstat = getattr(sn, 'refstat', None)
    rows = np.where(chains[:, k] > 0)[0] if chains.size > 0 else []
    if len(rows) > 0 and refstat is not None and len(refstat) > k:
        inchain = np.where(chains[rows[0], :] > 0)[0]
        nj = np.asarray(sn.njobs, dtype=float)[inchain]
        cand = [refstat[c] for c in inchain[(nj > 0) & np.isfinite(nj)]] or [refstat[k]]
        c0 = float(cand[0])
        if np.isfinite(c0) and c0 == round(c0) and 0 <= c0 < sn.nstations:
            ist = int(c0)
    if ist < 0:
        for i in range(sn.nstations):
            ph = sn.proc[i][k] if i < len(sn.proc) and k < len(sn.proc[i]) else None
            if ph is not None and np.count_nonzero(np.nan_to_num(np.asarray(ph[0], dtype=float))) > 0:
                ist = i
                break
    if ist < 0:
        raise ValueError("Class '%s' has no jobs, no reference station and no station serving it: "
                         "it cannot be written as source." % sn.classnames[k])
    return int(sn.stationToNode[ist])


def _node_type_text(t) -> str:
    """NodeType.toText spelling of a node type stored as the enum or its integer id."""
    from ...constants import NodeType
    return NodeType.toText(NodeType(int(t)))


def _refuse_node(sn, i: int, who: str):
    """Raise as QN2MATLAB/QN2JAVA do for a node (Cache, Logger, Place, Transition) whose state they do not write."""
    raise ValueError("%s cannot generate code for node %s of type %s."
                     % (who, sn.nodenames[i], _node_type_text(sn.nodetype[i])))


def _qn_node_decls(sn, who: str = 'QN2MATLAB'):
    """(i, kind, name, extra) for every node, in the vocabulary of QN2MATLAB.

    kind is Source, Delay, Queue, Router, Fork, Join or Sink; a ClassSwitch becomes a Router
    whose switching is carried by the routing matrix, as in the MATLAB generator. Any other
    node type is refused by name, as WHO: writing it would rebuild a different model.
    """
    decls = []
    for i in range(sn.nnodes):
        name = sn.nodenames[i]
        t = _type_name(sn.nodetype[i])
        if t == 'SOURCE':
            decls.append((i, 'Source', name, None))
        elif t == 'DELAY':
            decls.append((i, 'Delay', name, None))
        elif t == 'QUEUE':
            ist = int(sn.nodeToStation[i])
            decls.append((i, 'Queue', name, (_sched_name(sn.sched[ist]), sn.nservers[ist])))
        elif t == 'ROUTER':
            decls.append((i, 'Router', name, None))
        elif t == 'FORK':
            decls.append((i, 'Fork', name, None))
        elif t == 'JOIN':
            decls.append((i, 'Join', name, _fork_of(sn, i, name)))
        elif t == 'SINK':
            decls.append((i, 'Sink', name, None))
        elif t == 'CLASSSWITCH':
            decls.append((i, 'ClassSwitch', name, None))
        else:
            _refuse_node(sn, i, who)
    return decls


def _qn_routes(sn):
    """(k, c, i, m, p) for every nonzero entry of sn.rtnodes not leaving a Sink.

    The probability out of a Fork is written as 1.0, since a Fork sends a copy along every
    outgoing link rather than choosing one. A Source row of a closed class is skipped, as in
    QN2MATLAB: a closed job never leaves a Source.
    """
    rtnodes = np.asarray(sn.rtnodes, dtype=float)
    K = sn.nclasses
    routes = []
    for k in range(K):
        for c in range(K):
            for i in range(sn.nnodes):
                ti = _type_name(sn.nodetype[i])
                if ti == 'SINK' or (ti == 'SOURCE' and np.isfinite(sn.njobs[k])):
                    continue
                for m in range(sn.nnodes):
                    p = rtnodes[i * K + k, m * K + c]
                    if p > 0:
                        routes.append((k, c, i, m, 1.0 if ti == 'FORK' else p))
    return routes


def _generate_python_code(sn: Any, model_name: str, output: TextIO, model: Any = None) -> None:
    """Write the Python twin of the script QN2MATLAB writes for the network structure."""
    output.write("# LINE Network Model - Generated Python Code\n")
    output.write("# This code recreates the network model\n\n")
    output.write("from line_solver import *\n")
    output.write("import numpy as np\n\n")
    output.write(f"{model_name} = Network('{model_name}')\n\n")

    output.write("# Block 1: nodes\n")
    output.write("node = {}\n")
    for i, kind, name, extra in _qn_node_decls(sn, 'qn2python'):
        if kind == 'Queue':
            sched, nservers = extra
            output.write(f"node[{i}] = Queue({model_name}, '{name}', SchedStrategy.{sched})\n")
            if nservers > 1:
                output.write(f"node[{i}].setNumberOfServers({_count(nservers, True)})\n")
        elif kind == 'Join':
            output.write(f"node[{i}] = Join({model_name}, '{name}', node[{extra}])\n")
        elif kind == 'ClassSwitch':
            output.write(f"node[{i}] = Router({model_name}, '{name}')  # Class switching is embedded in the routing matrix\n")
        else:
            output.write(f"node[{i}] = {kind}({model_name}, '{name}')\n")

    output.write("\n# Block 2: classes\n")
    output.write("jobclass = {}\n")
    for k in range(sn.nclasses):
        name = sn.classnames[k]
        prio = int(sn.classprio[k])
        if np.isinf(sn.njobs[k]):
            output.write(f"jobclass[{k}] = OpenClass({model_name}, '{name}', {prio})\n")
        else:
            output.write(f"jobclass[{k}] = ClosedClass({model_name}, '{name}', {int(sn.njobs[k])}, "
                         f"node[{_qn_class_ref_node(sn, k)}], {prio})\n")

    output.write("\n# Block 3: arrival and service processes\n")
    for ist in range(sn.nstations):
        for k in range(sn.nclasses):
            spec = _qn_process_spec(sn, model, ist, k)
            if spec is None:
                continue
            verb, kind, args = spec
            i = int(sn.stationToNode[ist])
            if kind == 'Replayer':
                expr = f"Replayer({args[0]._file_path!r})"
            elif kind == 'Trace':
                data = getattr(args[0], '_histogram_data', None)
                data = args[0].trace if data is None else data
                expr = f"Trace(np.array({np.asarray(data).tolist()!r}))"
            elif kind == 'Immediate':
                expr = "Immediate()"
            elif kind == 'Disabled':
                expr = "Disabled()"
            elif kind == 'Exp':
                expr = f"Exp.fitMean({_py_num(args[0])})"
            elif kind == 'APH':
                expr = f"APH.fitMeanAndSCV({_py_num(args[0])}, {_py_num(args[1])})"
            else:
                expr = f"Erlang({_py_num(args[0])}, {args[1]})"
            output.write(f"node[{i}].{verb}(jobclass[{k}], {expr})  # ({sn.nodenames[i]},{sn.classnames[k]})\n")

    output.write("\n# Block 4: topology\n")
    output.write(f"P = {model_name}.initRoutingMatrix()\n")
    for k, c, i, m, p in _qn_routes(sn):
        output.write(f"P.set(jobclass[{k}], jobclass[{c}], node[{i}], node[{m}], {_py_num(p)})  "
                     f"# ({sn.nodenames[i]},{sn.classnames[k]}) -> ({sn.nodenames[m]},{sn.classnames[c]})\n")
    output.write(f"{model_name}.link(P)\n")


def _generate_matlab_code(sn: Any, model_name: str, output: TextIO, model: Any = None) -> None:
    """Port of matlab/src/io/QN2MATLAB.m.

    The script is the one QN2MATLAB writes, with one difference: numbers are written with
    enough digits to read back to the same double, where the MATLAB generator writes %f and
    loses everything past the sixth decimal.
    """
    output.write(f"model = Network('{model_name}');\n")
    output.write("\n%% Block 1: nodes\n")
    for i, kind, name, extra in _qn_node_decls(sn):
        j = i + 1
        if kind == 'Queue':
            sched, nservers = extra
            output.write(f"node{{{j}}} = Queue(model, '{name}', SchedStrategy.{sched});\n")
            if nservers > 1:
                output.write(f"node{{{j}}}.setNumServers({_count(nservers, False)});\n")
        elif kind == 'Delay':
            output.write(f"node{{{j}}} = DelayStation(model, '{name}');\n")
        elif kind == 'Join':
            output.write(f"node{{{j}}} = Join(model, '{name}', node{{{extra + 1}}});\n")
        elif kind == 'ClassSwitch':
            output.write(f"node{{{j}}} = Router(model, '{name}'); % Class switching is embedded in the routing matrix \n")
        else:
            output.write(f"node{{{j}}} = {kind}(model, '{name}');\n")

    output.write("\n%% Block 2: classes\n")
    for k in range(sn.nclasses):
        name = sn.classnames[k]
        prio = int(sn.classprio[k])
        if np.isinf(sn.njobs[k]):
            output.write(f"jobclass{{{k + 1}}} = OpenClass(model, '{name}', {prio});\n")
        else:
            output.write(f"jobclass{{{k + 1}}} = ClosedClass(model, '{name}', {int(sn.njobs[k])}, "
                         f"node{{{_qn_class_ref_node(sn, k) + 1}}}, {prio});\n")
    output.write("\n")

    for ist in range(sn.nstations):
        for k in range(sn.nclasses):
            spec = _qn_process_spec(sn, model, ist, k)
            if spec is None:
                continue
            verb, kind, args = spec
            i = int(sn.stationToNode[ist])
            if kind in ('Replayer', 'Trace'):
                expr = f"{kind}('{args[0]._file_path}')"
            elif kind == 'Immediate':
                expr = "Immediate()"
            elif kind == 'Disabled':
                expr = "Disabled.getInstance()"
            elif kind == 'Exp':
                expr = f"Exp.fitMean({_num(args[0])})"
            elif kind == 'APH':
                expr = f"APH.fitMeanAndSCV({_num(args[0])},{_num(args[1])})"
            else:
                expr = f"Erlang({_num(args[0])},{args[1]})"
            output.write(f"node{{{i + 1}}}.{verb}(jobclass{{{k + 1}}}, {expr}); % ({sn.nodenames[i]},{sn.classnames[k]})\n")

    output.write("\n%% Block 3: topology\n")
    output.write("P = model.initRoutingMatrix(); % initialize routing matrix \n")
    for k, c, i, m, p in _qn_routes(sn):
        output.write(f"P{{{k + 1},{c + 1}}}({i + 1},{m + 1}) = {_num(p)}; "
                     f"% ({sn.nodenames[i]},{sn.classnames[k]}) -> ({sn.nodenames[m]},{sn.classnames[c]})\n")
    output.write("model.link(P);\n")


def _emit(generator, model, model_name, file, *extra) -> Optional[str]:
    """Run a generator into a string, a path or an open file, as the qn2* entry points do."""
    close_file = False
    if file is None:
        output = StringIO()
    elif isinstance(file, str):
        output = open(file, 'w')
        close_file = True
    else:
        output = file
    try:
        generator(model, model_name, output, *extra)
        if file is None:
            return output.getvalue()
        return None
    finally:
        if close_file:
            output.close()


def _struct_and_model(model):
    if hasattr(model, 'getStruct'):
        return model.getStruct(), model
    return model, None


def qn2matlab(model: Any, model_name: str = 'myModel',
              file: Optional[Union[str, TextIO]] = None) -> Optional[str]:
    """
    Generate a MATLAB script that recreates a Network model.

    Args:
        model: Network model or NetworkStruct. With a NetworkStruct the Replayer and
            Trace processes cannot be recognised and are written as their fitted law.
        model_name: Name given to the Network in the script (the variable is `model`)
        file: Output file path or file object. If None, returns string.

    Returns:
        Generated MATLAB code as string if file is None, otherwise None.

    References:
        MATLAB: matlab/src/io/QN2MATLAB.m
    """
    sn, net = _struct_and_model(model)
    return _emit(lambda _m, name, out: _generate_matlab_code(sn, name, out, net), model, model_name, file)


def _py_list(values) -> str:
    return '[' + ', '.join(_py_num(v) for v in np.ravel(np.asarray(values, dtype=float))) + ']'


def _py_matrix(values) -> str:
    rows = np.atleast_2d(np.asarray(values, dtype=float))
    return 'np.array([' + ', '.join(_py_list(r) for r in rows) + '])'


def _py_dist(dist, what: str) -> str:
    """Python constructor that rebuilds DIST exactly, parameter for parameter.

    Raises on a law this generator cannot spell, naming WHAT carries it, rather than writing
    a script whose model differs from the original.
    """
    if isinstance(dist, (int, float, np.integer, np.floating)):
        return _py_num(dist)
    name = type(dist).__name__
    if name == 'Exp':
        return f"Exp({_py_num(dist.getRate())})"
    if name == 'Erlang':
        return f"Erlang({_py_num(dist._phase_rate)}, {int(dist._phases)})"
    if name == 'HyperExp':
        return f"HyperExp({_py_list(dist._probs)}, {_py_list(dist._rates)})"
    if name == 'Det':
        return f"Det({_py_num(dist._value)})"
    if name in ('Immediate', 'Disabled'):
        return f"{name}()"
    if name in ('APH', 'PH'):
        return f"{name}({_py_list(dist.getInitProb())}, {_py_matrix(dist._T)})"
    if name == 'Cox2':
        mu, phi = _coxian_params(dist)
        return f"Cox2({_py_num(mu[0])}, {_py_num(mu[1])}, {_py_num(phi[0])})"
    if name == 'Coxian':
        # the rates as given: -diag(T) is 1/(1/mu), which need not be mu
        mu, phi = _coxian_params(dist)
        return f"Coxian({_py_list(mu)}, {_py_list(phi)})"
    if name == 'MAP':
        return f"MAP({_py_matrix(dist._D0)}, {_py_matrix(dist._D1)})"
    if name == 'ME':
        return f"ME({_py_list(np.ravel(dist._alpha))}, {_py_matrix(dist._A)})"
    if name == 'RAP':
        return f"RAP({_py_matrix(dist._H0)}, {_py_matrix(dist._H1)})"
    if name == 'MMPP2':
        return (f"MMPP2({_py_num(dist._lambda0)}, {_py_num(dist._lambda1)}, "
                f"{_py_num(dist._sigma0)}, {_py_num(dist._sigma1)})")
    if name == 'Normal':
        return f"Normal({_py_num(dist._mean_val)}, {_py_num(dist._std)})"
    if name == 'Binomial':
        return f"Binomial({int(round(dist._n))}, {_py_num(dist._p)})"
    if name == 'DiscreteUniform':
        return f"DiscreteUniform({int(round(dist._a))}, {int(round(dist._b))})"
    if name == 'Zipf':
        return f"Zipf({_py_num(dist._s)}, {int(round(dist._n))})"
    if name == 'DiscreteSampler':
        return f"DiscreteSampler({_py_list(dist._probs)}, {_py_list(dist._values)})"
    if name == 'Replayer':
        return f"Replayer({_replayer_file(dist, 'lqn2python')!r})"
    if name == 'Gamma':
        return f"Gamma({_py_num(dist._shape)}, {_py_num(dist._scale)})"
    if name == 'Weibull':
        return f"Weibull({_py_num(dist._shape)}, {_py_num(dist._scale)})"
    if name == 'Pareto':
        return f"Pareto({_py_num(dist._alpha)}, {_py_num(dist._scale)})"
    if name == 'Lognormal':
        return f"Lognormal({_py_num(dist._mu)}, {_py_num(dist._sigma)})"
    if name == 'Uniform':
        return f"Uniform({_py_num(dist._min)}, {_py_num(dist._max)})"
    if name in ('Geometric', 'Bernoulli'):
        return f"{name}({_py_num(dist._p)})"
    if name == 'Poisson':
        return f"Poisson({_py_num(dist._lambda)})"
    raise ValueError("%s: the %s distribution cannot be written as Python source." % (what, name))


def lqn2python(model: Any, model_name: str = 'myLayeredModel',
               file: Optional[Union[str, TextIO]] = None) -> Optional[str]:
    """
    Generate a native Python script that recreates a LayeredNetwork model.

    The script covers what lqn2java writes (processors, tasks, entries, activities, calls,
    replies, think times, replication and the activity precedences), read from the model
    objects rather than from the struct, so that every distribution is written with its own
    parameters and no fitting takes place. It also carries open arrivals, forwarding, phases,
    activity think times, task priorities, fan-in/fan-out, processor quantum and speed factor,
    and routed call groups. Cache, setup and function tasks, admission constraints and rate
    dependences are refused by name.

    Args:
        model: LayeredNetwork model
        model_name: Name given to the LayeredNetwork in the script (the variable is `model`)
        file: Output file path or file object. If None, returns string.

    Returns:
        Generated Python code as string if file is None, otherwise None.

    References:
        MATLAB: matlab/src/io/LQN2JAVA.m (same walk, Python spelling)
    """
    _require_model(model, 'lqn2python', layered=True)
    return _emit(_generate_lqn_python_code, model, model_name, file)


def _generate_lqn_python_code(model: Any, model_name: str, output: TextIO) -> None:
    """Write the Python script lqn2python describes."""
    procs = list(model.processors)
    tasks = list(model.tasks)
    entries = [e for t in tasks for e in t.entries]
    acts = [a for t in tasks for a in t.activities]
    pvar = {id(p): 'P%d' % (h + 1) for h, p in enumerate(procs)}
    tvar = {id(t): 'T%d' % (i + 1) for i, t in enumerate(tasks)}
    evar = {id(e): 'E%d' % (i + 1) for i, e in enumerate(entries)}
    avar = {id(a): 'A%d' % (i + 1) for i, a in enumerate(acts)}
    ebyname = {e.name: e for e in entries}
    abyname = {a.name: a for a in acts}

    def refuse(obj, what):
        raise ValueError("lqn2python: %s '%s' carries %s, which cannot be written as source yet."
                         % (type(obj).__name__, obj.name, what))

    def check_server(obj):
        if hasattr(obj, 'hasLinearConstraints') and obj.hasLinearConstraints():
            refuse(obj, 'an admission constraint')
        if hasattr(obj, 'hasRateDependence') and obj.hasRateDependence():
            refuse(obj, 'a rate dependence')
        if hasattr(obj, 'hasServerPools') and obj.hasServerPools():
            refuse(obj, 'server pools')

    def act_ref(a):
        a = abyname[a] if isinstance(a, str) else a
        return avar[id(a)]

    output.write("# LINE LayeredNetwork Model - Generated Python Code\n")
    output.write("# This code recreates the layered network model\n\n")
    output.write("from line_solver import *\n")
    output.write("import numpy as np\n\n")
    output.write(f"model = LayeredNetwork('{model_name}')\n\n")

    output.write("# Processors\n")
    for p in procs:
        if type(p).__name__ != 'Processor':
            refuse(p, 'a %s host' % type(p).__name__)
        check_server(p)
        v = pvar[id(p)]
        output.write(f"{v} = Processor(model, '{p.name}', {_count(p.multiplicity, True)}, "
                     f"SchedStrategy.{_sched_name(p.sched_strategy)})\n")
        if p.getReplication() != 1:
            output.write(f"{v}.setReplication({int(p.getReplication())})\n")
        if p.getQuantum() != 0.0:
            output.write(f"{v}.setQuantum({_py_num(p.getQuantum())})\n")
        if p.getSpeedFactor() != 1.0:
            output.write(f"{v}.setSpeedFactor({_py_num(p.getSpeedFactor())})\n")

    output.write("\n# Tasks\n")
    for t in tasks:
        if type(t).__name__ != 'Task':
            refuse(t, 'a %s' % type(t).__name__)
        check_server(t)
        if t.setup_time is not None or t.delay_off_time is not None:
            refuse(t, 'a setup or delay-off time')
        v = tvar[id(t)]
        output.write(f"{v} = Task(model, '{t.name}', {_count(t.multiplicity, True)}, "
                     f"SchedStrategy.{_sched_name(t.sched_strategy)}).on({pvar[id(t.processor)]})\n")
        if t.getReplication() != 1:
            output.write(f"{v}.setReplication({int(t.getReplication())})\n")
        if t.priority != 0:
            output.write(f"{v}.setPriority({int(t.priority)})\n")
        if t.think_time is not None:
            output.write(f"{v}.setThinkTime({_py_dist(t.think_time, 'think time of task ' + t.name)})\n")
        for src, val in t.getFanIn().items():
            output.write(f"{v}.setFanIn({src!r}, {int(val)})\n")
        for dst, val in t.getFanOut().items():
            output.write(f"{v}.setFanOut({dst!r}, {int(val)})\n")

    output.write("\n# Entries\n")
    for e in entries:
        if type(e).__name__ != 'Entry':
            refuse(e, 'an %s' % type(e).__name__)
        v = evar[id(e)]
        output.write(f"{v} = Entry(model, '{e.name}').on({tvar[id(e.task)]})\n")
        if e.getArrival() is not None:
            output.write(f"{v}.setArrival({_py_dist(e.getArrival(), 'arrival of entry ' + e.name)})\n")
    for e in entries:
        for dst, prob in zip(e.getForwardingDests(), e.getForwardingProbs()):
            output.write(f"{evar[id(e)]}.forward({evar[id(ebyname[dst])]}, {_py_num(prob)})\n")

    output.write("\n# Activities\n")
    for a in acts:
        v = avar[id(a)]
        demand = _py_dist(a.host_demand, 'host demand of activity ' + a.name)
        line = f"{v} = Activity(model, '{a.name}', {demand}).on({tvar[id(a.task)]})"
        if a.bound_entry is not None:
            line += f".boundTo({evar[id(a.bound_entry)]})"
        for entry, mean_calls, call_type in a.calls:
            verb = 'synchCall' if _type_name(call_type) == 'SYNC' else 'asynchCall'
            line += f".{verb}({evar[id(entry)]}, {_py_num(mean_calls)})"
        if a.reply_entry is not None:
            line += f".repliesTo({evar[id(a.reply_entry)]})"
        output.write(line + "\n")
        if a.phase != 1:
            output.write(f"{v}.setPhase({int(a.phase)})\n")
        if not (isinstance(a.think_time, (int, float)) and a.think_time == 0):
            output.write(f"{v}.setThinkTime({_py_dist(a.think_time, 'think time of activity ' + a.name)})\n")
        for strategy, group in a.call_groups:
            members = ', '.join(evar[id(e)] for e in group)
            output.write(f"{v}.record_call_group(RoutingStrategy.{_type_name(strategy)}, [{members}])\n")

    precs = [(t, p) for t in tasks for p in t.precedences]
    if precs:
        output.write("\n# Activity precedences\n")
    for t, prec in precs:
        kind = model._get_prec_type(prec)
        pre = [act_ref(x) for x in model._get_prec_pre_activities(prec)]
        post = [act_ref(x) for x in model._get_prec_post_activities(prec)]
        body = [act_ref(x) for x in model._get_prec_activities(prec)]
        kname = _type_name(kind)
        if kname == 'SERIAL':
            expr = f"ActivityPrecedence.Serial({', '.join(body)})"
        elif kname == 'LOOP':
            expr = (f"ActivityPrecedence.Loop({pre[0]}, [{', '.join(body)}], "
                    f"{_py_num(model._get_prec_count(prec))})")
        elif kname == 'PARALLEL' and len(post) > 1:
            expr = f"ActivityPrecedence.AndFork({pre[0]}, [{', '.join(post)}])"
        elif kname == 'PARALLEL':
            quorum = getattr(prec, 'pre_params', None)
            qstr = '' if quorum is None else f", {int(round(float(np.ravel(quorum)[0])))}"
            expr = f"ActivityPrecedence.AndJoin([{', '.join(pre)}], {post[0]}{qstr})"
        elif kname == 'CHOICE' and len(post) > 1:
            probs = model._get_prec_probabilities(prec)
            expr = f"ActivityPrecedence.OrFork({pre[0]}, [{', '.join(post)}], {_py_list(probs)})"
        elif kname == 'CHOICE':
            expr = f"ActivityPrecedence.OrJoin([{', '.join(pre)}], {post[0]})"
        else:
            raise ValueError("lqn2python: a %s precedence of task '%s' cannot be written as source yet."
                             % (kname, t.name))
        output.write(f"{tvar[id(t)]}.addPrecedence({expr})\n")


def _m_scalar(v) -> str:
    """scalar2code of LQN2MATLAB: NaN/Inf literals, an integer below 1e15 as %d, else the first %.15g..%.17g reading back."""
    v = float(v)
    if np.isnan(v):
        return 'NaN'
    if np.isinf(v):
        return 'Inf' if v > 0 else '-Inf'
    if v == round(v) and abs(v) < 1e15:
        return '%d' % int(v)
    for p in (15, 16, 17):
        s = '%.*g' % (p, v)
        if float(s) == v:
            return s
    return repr(v)


def _m_num(v) -> str:
    """num2code of LQN2MATLAB: a scalar, or a row or matrix as [a, b; c, d]; a 1-D array is a row."""
    if isinstance(v, (bool, np.bool_)):
        return 'true' if v else 'false'
    a = np.asarray(v, dtype=float)
    if a.size == 0:
        return '[]'
    if a.ndim == 0:
        return _m_scalar(a)
    if a.ndim > 2:
        raise ValueError('lqn2matlab cannot write an array of more than two dimensions.')
    rows = np.atleast_2d(a)
    if a.size == 1:
        return _m_scalar(a.flat[0])
    return '[' + '; '.join(', '.join(_m_scalar(x) for x in r) for r in rows) + ']'


def _m_col(v) -> str:
    """A vector as a MATLAB column."""
    return _m_num(np.asarray(v, dtype=float).reshape(-1, 1))


def _m_str(s) -> str:
    """MATLAB single-quoted literal."""
    return "'" + str(s).replace("'", "''") + "'"


def _m_names(names) -> str:
    """Cell array of names, from names or named objects."""
    if isinstance(names, str) or hasattr(names, 'name'):
        names = [names]
    return '{' + ', '.join(_m_str(n if isinstance(n, str) else n.name) for n in names) + '}'


def _m_dist(dist, what: str) -> str:
    """dist2code of LQN2MATLAB: the MATLAB constructor that rebuilds DIST with the same parameters.

    A law with no MATLAB spelling here is written as the APH fitted to its mean and SCV, with a
    warning, as LQN2MATLAB does.
    """
    from .logging import line_warning
    if dist is None:
        return '[]'
    if isinstance(dist, (int, float, np.integer, np.floating)):
        return _m_scalar(dist)
    name = type(dist).__name__
    if name in ('Immediate', 'Disabled'):
        return f"{name}()"
    if name == 'Exp':
        return f"Exp({_m_scalar(dist.getRate())})"
    if name == 'Det':
        return f"Det({_m_scalar(dist._value)})"
    if name == 'Erlang':
        return f"Erlang({_m_scalar(dist._phase_rate)}, {int(dist._phases)})"
    if name == 'HyperExp':
        return f"HyperExp({_m_num(np.ravel(dist._probs))}, {_m_num(np.ravel(dist._rates))})"
    if name == 'Gamma':
        return f"Gamma({_m_scalar(dist._shape)}, {_m_scalar(dist._scale)})"
    if name == 'Weibull':
        return f"Weibull({_m_scalar(dist._shape)}, {_m_scalar(dist._scale)})"
    if name == 'Pareto':
        return f"Pareto({_m_scalar(dist._alpha)}, {_m_scalar(dist._scale)})"
    if name == 'Lognormal':
        return f"Lognormal({_m_scalar(dist._mu)}, {_m_scalar(dist._sigma)})"
    if name == 'Uniform':
        return f"Uniform({_m_scalar(dist._min)}, {_m_scalar(dist._max)})"
    if name in ('Geometric', 'Bernoulli'):
        return f"{name}({_m_scalar(dist._p)})"
    if name == 'Poisson':
        return f"Poisson({_m_scalar(dist._lambda)})"
    if name in ('Normal', 'DiscreteUniform', 'Binomial'):
        a, b = {'Normal': ('_mean_val', '_std'), 'DiscreteUniform': ('_a', '_b'), 'Binomial': ('_n', '_p')}[name]
        return f"{name}({_m_scalar(getattr(dist, a))}, {_m_scalar(getattr(dist, b))})"
    if name in ('Coxian', 'Cox2'):
        mu, phi = _coxian_params(dist)
        if name == 'Cox2':
            # MATLAB's Cox2.fitMeanAndSCV returns the three-parameter Coxian(mu1, mu2, phi1)
            return f"Coxian({_m_scalar(mu[0])}, {_m_scalar(mu[1])}, {_m_scalar(phi[0])})"
        return f"Coxian({_m_num(mu)}, {_m_num(phi)})"
    if name in ('APH', 'PH'):
        return f"{name}({_m_num(np.ravel(dist.getInitProb()))}, {_m_num(dist._T)})"
    if name == 'ME':
        return f"ME({_m_num(np.ravel(dist._alpha))}, {_m_num(dist._A)})"
    if name == 'RAP':
        return f"RAP({_m_num(dist._H0)}, {_m_num(dist._H1)})"
    if name == 'MMPP2':
        return (f"MMPP2({_m_scalar(dist._lambda0)}, {_m_scalar(dist._lambda1)}, "
                f"{_m_scalar(dist._sigma0)}, {_m_scalar(dist._sigma1)})")
    if name == 'MAP':
        return f"MAP({_m_num(dist._D0)}, {_m_num(dist._D1)})"
    if name == 'Zipf':
        return f"Zipf({_m_scalar(dist._s)}, {_m_scalar(dist._n)})"
    if name == 'DiscreteSampler':
        return f"DiscreteSampler({_m_num(np.ravel(dist._probs))}, {_m_num(np.ravel(dist._values))})"
    if name in ('Replayer', 'Trace') and getattr(dist, '_file_path', None):
        return f"{name}({_m_str(dist._file_path)})"
    mean, scv = float(dist.getMean()), float(dist.getSCV())
    line_warning('lqn2matlab', 'lqn2matlab writes the %s distribution of %s as an APH fitted to its mean and SCV.',
                 name, what)
    return f"APH.fitMeanAndSCV({_m_scalar(mean)}, {_m_scalar(scv)})"


def lqn2matlab(model: Any, model_name: Optional[str] = None,
               file: Optional[Union[str, TextIO]] = None) -> Optional[str]:
    """
    Generate the MATLAB script that rebuilds a LayeredNetwork, as LQN2MATLAB writes it.

    The script declares, in model order, the processors (with quantum and speed factor,
    replication, admission constraints, load dependence and server pools), the tasks (Task,
    CacheTask, SetupTask, FunctionTask, with think, setup and delay-off times, replication,
    priority, fan-in and fan-out), the entries (Entry, ItemEntry, open arrivals, forwarding),
    the activities (host demand, bound entry, think time, phase, synchronous calls, routed call
    groups, asynchronous calls), the replies and the precedences, all read from the model
    objects so that every value is written with its own parameters. Numbers use the shortest
    spelling that reads back bit for bit. A class or joint dependence is a Python callable,
    which has no MATLAB spelling, and is refused by name.

    Args:
        model: LayeredNetwork model
        model_name: Name given to the LayeredNetwork in the script (the variable is `model`);
            defaults to the model's own name
        file: Output file path or file object. If None, returns string.

    Returns:
        Generated MATLAB code as string if file is None, otherwise None.

    References:
        MATLAB: matlab/src/io/LQN2MATLAB.m
    """
    if model_name is None:
        model_name = model.getName() if hasattr(model, 'getName') else getattr(model, 'name', 'myLayeredModel')
    _require_model(model, 'lqn2matlab', layered=True)
    return _emit(_generate_lqn_matlab_code, model, model_name, file)


def _generate_lqn_matlab_code(model: Any, model_name: str, output: TextIO) -> None:
    """Port of lqn2matlab_generate in matlab/src/io/LQN2MATLAB.m onto the Python layered objects."""
    from .logging import line_warning
    from ...constants import GlobalConstants
    hosts = list(model.processors)
    tasks = list(model.tasks)
    entries = list(model.entries) or [e for t in tasks for e in t.entries]
    acts = list(model.activities) or [a for t in tasks for a in t.activities]
    hidx = {id(h): p + 1 for p, h in enumerate(hosts)}
    tidx = {id(t): i + 1 for i, t in enumerate(tasks)}
    eidx = {id(e): i + 1 for i, e in enumerate(entries)}
    ename = {e.name: i + 1 for i, e in enumerate(entries)}
    abyname = {a.name: a for a in acts}
    w = output.write

    def refuse(obj, what):
        raise ValueError("lqn2matlab: %s '%s' carries %s, which has no MATLAB spelling."
                         % (type(obj).__name__, obj.name, what))

    def entry_ref(e):
        name = e if isinstance(e, str) else e.name
        return 'E{%d}' % ename[name] if name in ename else _m_str(name)

    def act_name(a):
        return a if isinstance(a, str) else a.name

    def time_arg(d):
        # a numeric argument takes the setter's numeric path (Immediate with mean FineTol, else Exp(1/m))
        return _m_scalar(d) if isinstance(d, (int, float, np.integer, np.floating)) else None

    def server_extras(v, el):
        if getattr(el, 'lincon_a', None) is not None and getattr(el, 'lincon_b', None) is not None:
            w(f"{v}.setConstraint({_m_num(el.lincon_a)}, {_m_col(el.lincon_b)});\n")
        for names, coeffs, cap in getattr(el, 'lincon_rows', []):
            w(f"{v}.addConstraint({_m_names(list(names))}, {_m_num(np.ravel(coeffs))}, {_m_scalar(cap)});\n")
        if getattr(el, 'lld_scaling', None) is not None:
            w(f"{v}.setLoadDependence({_m_num(np.ravel(el.lld_scaling))});\n")
        if getattr(el, 'lcd_scaling', None) is not None:
            refuse(el, 'a class dependence, a Python callable')
        if getattr(el, 'ljd_scaling', None) is not None:
            refuse(el, 'a joint dependence, a Python callable')
        for pool in getattr(el, 'server_pools', []):
            w(f"{v}.addServerType(ServerType({_m_str(pool['name'])}, {_m_scalar(pool['count'])}, "
              f"{_m_names(pool['compatible'])}, {_m_scalar(pool['rate'])}));\n")

    def sched(s):
        return 'SchedStrategy.' + _sched_name(s)

    w("% LayeredNetwork generated by LQN2MATLAB\n")
    w(f"model = LayeredNetwork({_m_str(model_name)});\n")
    w("P = {}; T = {}; E = {}; A = {};\n")

    w("\n%% Block 1: processors\n")
    for p, h in enumerate(hosts, 1):
        hv = 'P{%d}' % p
        # Python leaves an unset quantum at 0.0 where MATLAB's Processor defaults to 0.001
        quantum = h.getQuantum() if h.getQuantum() != 0.0 else 0.001
        speed = h.getSpeedFactor()
        if quantum != 0.001 or speed != 1:
            w(f"{hv} = Processor(model, {_m_str(h.name)}, {_m_scalar(h.multiplicity)}, {sched(h.sched_strategy)}, "
              f"{_m_scalar(quantum)}, {_m_scalar(speed)});\n")
        else:
            w(f"{hv} = Processor(model, {_m_str(h.name)}, {_m_scalar(h.multiplicity)}, {sched(h.sched_strategy)});\n")
        if h.getReplication() != 1:
            w(f"{hv}.setReplication({_m_scalar(h.getReplication())});\n")
        server_extras(hv, h)

    w("\n%% Block 2: tasks\n")
    for t, tk in enumerate(tasks, 1):
        tv = 'T{%d}' % t
        cls = type(tk).__name__
        on = 'P{%d}' % hidx[id(tk.processor)]
        if cls == 'CacheTask':
            repl = tk.replacement_strategy
            if not hasattr(repl, 'name'):
                from ...lang.base import ReplacementStrategy
                repl = ReplacementStrategy(int(repl))
            w(f"{tv} = CacheTask(model, {_m_str(tk.name)}, {_m_scalar(tk.total_items)}, {_m_num(tk.cache_capacity)}, "
              f"ReplacementStrategy.{repl.name}, {_m_scalar(tk.multiplicity)}, {sched(tk.sched_strategy)}).on({on});\n")
            if getattr(tk, 'retrieval', False):
                w(f"{tv}.setRetrieval(true);\n")
        elif cls in ('SetupTask', 'FunctionTask'):
            w(f"{tv} = {cls}(model, {_m_str(tk.name)}, {_m_scalar(tk.multiplicity)}, {sched(tk.sched_strategy)}).on({on});\n")
        else:
            if cls != 'Task':
                line_warning('lqn2matlab', 'Task %s is a %s, which lqn2matlab declares as a plain Task.', tk.name, cls)
            w(f"{tv} = Task(model, {_m_str(tk.name)}, {_m_scalar(tk.multiplicity)}, {sched(tk.sched_strategy)}).on({on});\n")
        for setter, d in (('setThinkTime', tk.think_time), ('setSetupTime', tk.setup_time),
                          ('setDelayOffTime', tk.delay_off_time)):
            if d is not None:
                w(f"{tv}.{setter}({time_arg(d) or _m_dist(d, 'task ' + tk.name)});\n")
        if tk.getReplication() != 1:
            w(f"{tv}.setReplication({_m_scalar(tk.getReplication())});\n")
        if tk.priority != 0:
            w(f"{tv}.setPriority({_m_scalar(tk.priority)});\n")
        fanin = tk.getFanIn()
        if len(fanin) > 1:
            refuse(tk, 'a fan-in from %d source tasks, where a MATLAB Task holds one' % len(fanin))
        for src, val in fanin.items():
            w(f"{tv}.setFanIn({_m_str(src)}, {_m_scalar(val)});\n")
        for dst, val in tk.getFanOut().items():
            w(f"{tv}.setFanOut({_m_str(dst)}, {_m_scalar(val)});\n")
        server_extras(tv, tk)

    w("\n%% Block 3: entries\n")
    for e, en in enumerate(entries, 1):
        ev = 'E{%d}' % e
        on = 'T{%d}' % tidx[id(en.task)]
        if type(en).__name__ == 'ItemEntry':
            pop = en.access_prob
            pop = _m_dist(pop, 'entry ' + en.name) if hasattr(pop, 'getMean') or type(pop).__name__ == 'DiscreteSampler' \
                else f"DiscreteSampler({_m_num(np.ravel(np.asarray(pop, dtype=float)))})"
            w(f"{ev} = ItemEntry(model, {_m_str(en.name)}, {_m_scalar(en.total_items)}, {pop}).on({on});\n")
        else:
            w(f"{ev} = Entry(model, {_m_str(en.name)}).on({on});\n")
        if en.getArrival() is not None:
            w(f"{ev}.setArrival({_m_dist(en.getArrival(), 'arrival of entry ' + en.name)});\n")
    for e, en in enumerate(entries, 1):
        for dst, prob in zip(en.getForwardingDests(), en.getForwardingProbs()):
            w(f"E{{{e}}}.forward({entry_ref(dst)}, {_m_scalar(prob)});\n")

    w("\n%% Block 4: activities\n")
    for a, ac in enumerate(acts, 1):
        av = 'A{%d}' % a
        demand = ac.host_demand if ac.host_demand is not None else 0
        dcode = time_arg(demand) or _m_dist(demand, 'host demand of activity ' + ac.name)
        bound = '' if ac.bound_entry is None else f".boundTo({entry_ref(ac.bound_entry)})"
        w(f"{av} = Activity(model, {_m_str(ac.name)}, {dcode}).on(T{{{tidx[id(ac.task)]}}}){bound};\n")
        tt = ac.think_time
        if not (isinstance(tt, (int, float, np.integer, np.floating)) and tt <= GlobalConstants.FineTol) and tt is not None:
            w(f"{av}.setThinkTime({time_arg(tt) or _m_dist(tt, 'think time of activity ' + ac.name)});\n")
        if ac.phase != 1:
            w(f"{av}.setPhase({_m_scalar(ac.phase)});\n")
        for entry, mean_calls, call_type in ac.calls:
            if _type_name(call_type) == 'SYNC':
                w(f"{av}.synchCall({entry_ref(entry)}, {_m_scalar(mean_calls)});\n")
        # the member calls are declared above, so the grouping is recorded without issuing them again
        for strategy, group in ac.call_groups:
            sname = _type_name(strategy)
            if sname not in ('RROBIN', 'JSQ'):
                raise ValueError('Call groups carry RROBIN or JSQ; routing strategy %s cannot be written out.' % sname)
            w(f"{av}.recordCallGroup(RoutingStrategy.{sname}, {_m_names(group)});\n")
        for entry, mean_calls, call_type in ac.calls:
            if _type_name(call_type) != 'SYNC':
                w(f"{av}.asynchCall({entry_ref(entry)}, {_m_scalar(mean_calls)});\n")

    w("\n%% Block 5: replies\n")
    for e, en in enumerate(entries, 1):
        for a, ac in enumerate(acts, 1):
            if ac.reply_entry is not en:
                continue
            if ac.task is not None and _sched_name(ac.task.sched_strategy) == 'REF':
                # repliesTo refuses a reference-task activity, so the reply is recorded directly
                w(f"E{{{e}}}.replyActivity{{end+1}} = {_m_str(ac.name)};\n")
            else:
                w(f"A{{{a}}}.repliesTo(E{{{e}}});\n")

    w("\n%% Block 6: precedences\n")
    for t, tk in enumerate(tasks, 1):
        for prec in tk.precedences:
            for code in _m_prec(model, tk, prec, act_name):
                w(f"T{{{t}}}.addPrecedence({code});\n")


def _m_prec(model, tk, prec, act_name):
    """precCode of LQN2MATLAB: the named ActivityPrecedence factories, activities by name; a Serial chain goes pairwise."""
    kname = _type_name(model._get_prec_type(prec))
    pre = [act_name(x) for x in model._get_prec_pre_activities(prec)]
    post = [act_name(x) for x in model._get_prec_post_activities(prec)]
    body = [act_name(x) for x in model._get_prec_activities(prec)]
    q = _m_str
    if kname == 'SERIAL':
        return [f"ActivityPrecedence.Serial({q(a)}, {q(b)})" for a, b in zip(body[:-1], body[1:])]
    if kname == 'LOOP':
        # Loop(pre, {body..., end}, count): the last activity is the end activity
        return [f"ActivityPrecedence.Loop({q(pre[0])}, {_m_names(body)}, {_m_scalar(model._get_prec_count(prec))})"]
    if kname == 'PARALLEL' and len(pre) == 1 and len(post) > 1:
        return [f"ActivityPrecedence.AndFork({q(pre[0])}, {_m_names(post)})"]
    if kname == 'PARALLEL' and len(post) == 1:
        quorum = getattr(prec, 'pre_params', None)
        qstr = '' if quorum is None else ', ' + _m_num(np.ravel(quorum))
        return [f"ActivityPrecedence.AndJoin({_m_names(pre)}, {q(post[0])}{qstr})"]
    if kname == 'CHOICE' and len(pre) == 1 and len(post) > 1:
        probs = np.ravel(np.asarray(model._get_prec_probabilities(prec), dtype=float))
        return [f"ActivityPrecedence.OrFork({q(pre[0])}, {_m_names(post)}, {_m_num(probs)})"]
    if kname == 'CHOICE' and len(post) == 1:
        return [f"ActivityPrecedence.OrJoin({_m_names(pre)}, {q(post[0])})"]
    if kname == 'CACHE_ACCESS' and len(pre) == 1:
        return [f"ActivityPrecedence.CacheAccess({q(pre[0])}, {_m_names(post)})"]
    raise ValueError("lqn2matlab: a %s precedence of task '%s' cannot be written as source."
                     % (kname, tk.name))


def _is_layered(model: Any) -> bool:
    from line_solver.layered import LayeredNetwork
    return isinstance(model, LayeredNetwork)


def _is_network(model: Any) -> bool:
    from line_solver.lang.network import Network
    return isinstance(model, Network)


def _require_model(model: Any, who: str, layered: bool) -> None:
    """Refuse a struct or any other object where the generator reads the model objects."""
    if (_is_layered(model) if layered else _is_network(model)):
        return
    expected = 'LayeredNetwork' if layered else 'Network'
    hint = ' Pass the model itself, not its getStruct().' if 'Struct' in type(model).__name__ else ''
    raise TypeError('%s supports a %s, got a %s.%s' % (who, expected, type(model).__name__, hint))


def line2python(model: Any, filename: Optional[str] = None) -> Optional[str]:
    """
    Export a LINE model to a Python script.

    Dispatches to qn2python for a Network and to lqn2python for a LayeredNetwork,
    as LINE2MATLAB dispatches to QN2MATLAB and LQN2MATLAB.

    Args:
        model: Network or LayeredNetwork model
        filename: Output Python file path. If None, the code is returned.

    Returns:
        Generated Python code as string if filename is None, otherwise None.

    References:
        MATLAB: matlab/src/io/LINE2MATLAB.m
    """
    _check_line_model(model, 'line2python')
    model_name = model.getName() if hasattr(model, 'getName') else getattr(model, 'name', 'model')
    if _is_layered(model):
        return lqn2python(model, model_name, filename)
    return qn2python(model, model_name, filename)


def _check_line_model(model: Any, who: str) -> None:
    """Refuse anything but a Network or a LayeredNetwork, as LINE2MATLAB and LINE2JAVA do."""
    if not (_is_layered(model) or _is_network(model)):
        raise ValueError('%s supports a Network or a LayeredNetwork, got a %s.' % (who, type(model).__name__))


def line2matlab(model: Any, filename: Optional[str] = None) -> Optional[str]:
    """
    Export a LINE model to a MATLAB script.

    Dispatches to qn2matlab for a Network and to lqn2matlab for a LayeredNetwork,
    and refuses anything else, as LINE2MATLAB does.

    Args:
        model: Network or LayeredNetwork model
        filename: Output .m file path. If None, the code is returned.

    Returns:
        Generated MATLAB code as string if filename is None, otherwise None.

    References:
        MATLAB: matlab/src/io/LINE2MATLAB.m
    """
    _check_line_model(model, 'line2matlab')
    model_name = model.getName() if hasattr(model, 'getName') else getattr(model, 'name', 'model')
    if _is_layered(model):
        return lqn2matlab(model, model_name, filename)
    return qn2matlab(model, model_name, filename)


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


def _ph_disabled(ph, mean) -> bool:
    """QN2JAVA's isnan(PH{ist,k}{1}): a disabled process has a NaN (or missing) representation."""
    if ph is None or not hasattr(ph, '__getitem__') or len(ph) == 0:
        return bool(np.isnan(mean))
    return bool(np.all(np.isnan(np.asarray(ph[0], dtype=float))))


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
                sched_prop = _sched_name(sn.sched[ist])
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
                # The Fork this Join closes, which sn.fj always names; see the
                # python emitter above.
                output.write(f'\t\tJoin node{i+1} = new Join(model, "{node_name}", node{_fork_of(sn, i, node_name) + 1});\n')
            elif type_name == 'SINK':
                output.write(f'\t\tSink node{i+1} = new Sink(model, "{node_name}");\n')
            elif type_name == 'CLASSSWITCH':
                output.write(f'\t\tRouter node{i+1} = new Router(model, "{node_name}"); // Dummy node, class switching is embedded in the routing matrix P \n')
            else:
                _refuse_node(sn, i, 'QN2JAVA')

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
            # zero-population class: the chain's reference station, as refreshRoutingMatrix requires
            if np.isinf(njobs):
                output.write(f'\t\tOpenClass jobclass{k+1} = new OpenClass(model, "{class_name}", {priority});\n')
            else:
                iref = _qn_class_ref_node(sn, k)
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

            sched_name = _sched_name(sn.sched[ist])

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
                    if _ph_disabled(PH[ist][k], mean_ik):
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
                    if _ph_disabled(PH[ist][k], mean_ik):
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
        # as QN2JAVA: an ext column block past the stations takes the open-class returns to the Source
        K = sn.nclasses
        ext = sn.nstations * K
        grown = np.zeros((max(rt.shape[0], ext + K), max(rt.shape[1], ext + K)))
        grown[:rt.shape[0], :rt.shape[1]] = rt
        rt = grown
        rt[ext:ext + K, ext:ext + K] = 0
        src = int(sn.nodeToStation[source_id])
        for k in range(K):
            if np.isinf(sn.njobs[k]):  # open class
                for ist in range(sn.nstations):
                    rt[ist * K + k, ext + k] = rt[ist * K + k, src * K + k]
                    rt[ist * K + k, src * K + k] = 0

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
                                    if type_name_i == 'SOURCE' and np.isfinite(sn.njobs[k]):
                                        continue  # no Source rows for closed classes, as QN2JAVA
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
    Generate a JLINE Java program that recreates a LayeredNetwork model and solves it with SolverLN.

    The program has the layout of MATLAB's LQN2JAVA (class ``TestSolver<NAME>`` in package
    ``jline.examples``, a ``main`` declaring processors, tasks, entries, activities and the
    activity precedences, then ``SolverLN.getEnsembleAvg()``). It is read from the native model
    objects, as lqn2python is, so every law is written with its own parameters (MATLAB's LQN2JAVA
    moment-fits Erlang/HyperExp/Coxian/APH through ``fitMeanAndSCV``) and every number to full
    precision; it covers the same features as lqn2python and refuses the same ones by name.

    Args:
        model: LayeredNetwork model
        model_name: Name given to the LayeredNetwork (and, capitalised, to the class)
        file: Output file path or file object. If None, returns string.

    Returns:
        Generated Java code as string if file is None, otherwise None.

    References:
        MATLAB: matlab/src/io/LQN2JAVA.m
    """
    _require_model(model, 'lqn2java', layered=True)
    return _emit(_generate_lqn_java_code, model, model_name, file)


def _java_double_string(x: float) -> str:
    """Java's Double.toString: the shortest round-trip digits, plain in [1e-3, 1e7) and d.dddE<n> outside it."""
    from decimal import Decimal
    x = float(x)
    if x == 0.0:
        return '-0.0' if np.signbit(x) else '0.0'
    sign = '-' if x < 0 else ''
    t = Decimal(repr(abs(x))).as_tuple()
    digits = ''.join(str(d) for d in t.digits)
    exp = len(digits) + t.exponent - 1
    digits = digits.rstrip('0') or '0'
    if 1e-3 <= abs(x) < 1e7:
        if exp >= 0:
            whole = digits[:exp + 1].ljust(exp + 1, '0')
            frac = digits[exp + 1:] or '0'
            return sign + whole + '.' + frac
        return sign + '0.' + '0' * (-exp - 1) + digits
    return sign + digits[0] + '.' + (digits[1:] or '0') + 'E' + str(exp)


def _java_exact(x) -> str:
    """CodeGenSupport.exact: an integer below 1e15 as 4.0, NaN and infinities by name, anything else as Double.toString."""
    x = float(x)
    if np.isnan(x):
        return 'Double.NaN'
    if np.isinf(x):
        return 'Double.POSITIVE_INFINITY' if x > 0 else 'Double.NEGATIVE_INFINITY'
    if x == np.rint(x) and abs(x) < 1e15:
        return '%d.0' % int(x)
    return _java_double_string(x)


def _java_g(x) -> str:
    """CodeGenSupport.fmtG: MATLAB's %g where it reads back exactly, the exact spelling where %g would round."""
    x = float(x)
    if np.isnan(x) or np.isinf(x):
        return _java_exact(x)
    s = '%g' % x
    return s if float(s) == x else _java_exact(x)


def _java_num(x) -> str:
    """jnum of LQN2JAVA.m: an integer below 2^53 as 4.0, NaN and infinities by name, else the first of %.15g/%.16g/%.17g
    that reads back exactly."""
    x = float(x)
    if np.isnan(x):
        return 'Double.NaN'
    if np.isinf(x):
        return 'Double.POSITIVE_INFINITY' if x > 0 else 'Double.NEGATIVE_INFINITY'
    if x == np.rint(x) and abs(x) < 2.0 ** 53:
        return '%d.0' % int(x)
    for p in (15, 16, 17):
        s = '%.*g' % (p, x)
        if float(s) == x:
            break
    return s


def _java_count(x) -> str:
    """fmtInt: a server count or multiplicity, Integer.MAX_VALUE when infinite."""
    return 'Integer.MAX_VALUE' if np.isinf(float(x)) else '%d' % int(x)


def _java_esc(s) -> str:
    """CodeGenSupport.jstr: the body of a Java string literal."""
    out = str(s).replace('\\', '\\\\').replace('"', '\\"')
    return out.replace('\n', '\\n').replace('\r', '\\r').replace('\t', '\\t')


def _java_str(s) -> str:
    return '"' + _java_esc(s) + '"'


def _java_doubles(values) -> str:
    """jarr of LQN2JAVA.m: a double[] literal."""
    return 'new double[]{' + ', '.join(_java_num(v) for v in np.ravel(np.asarray(values, dtype=float))) + '}'


def _java_matrix(values) -> str:
    """jmat of LQN2JAVA.m: a Matrix from a double[][] literal; a vector is one row."""
    rows = np.atleast_2d(np.asarray(values, dtype=float))
    return 'new Matrix(new double[][]{' + ', '.join(
        '{' + ', '.join(_java_num(v) for v in r) + '}' for r in rows) + '})'


def _java_list(values, fmt) -> str:
    return 'java.util.Arrays.asList(' + ', '.join(fmt(v) for v in values) + ')'


def _coxian_params(dist):
    """(mu, phi) of a Coxian as its constructor took them, phi ending in 1; read off T for one built otherwise."""
    mu = getattr(dist, '_mu_in', None)
    phi = getattr(dist, '_phi_in', None)
    if mu is None or phi is None:
        T = np.asarray(dist._T, dtype=float)
        mu = -np.diag(T)
        phi = -T.sum(axis=1) / mu
        phi[-1] = 1.0
    return np.asarray(mu, dtype=float), np.asarray(phi, dtype=float)


def _replayer_file(dist, who):
    f = getattr(dist, '_file_path', None)
    if f is None or getattr(dist, '_file_path_is_temp', False):
        raise ValueError("%s cannot write a %s built from in-memory samples: it names no trace file."
                         % (who, type(dist).__name__))
    return str(f)


def _java_dist(dist, what: str) -> str:
    """javaDist of LQN2JAVA.m: the JLINE constructor that rebuilds DIST with its own parameters in constructor order.

    Only a law with no spelling here is written as the APH fitted to its mean and SCV, with a warning, as MATLAB does.
    A numeric value is the mean a numeric setter takes (Immediate at or below FineTol, else Exp(1/mean)).
    """
    from .logging import line_warning
    from ...constants import GlobalConstants
    if dist is None:
        return 'new Immediate()'
    if isinstance(dist, (int, float, np.integer, np.floating)):
        m = float(dist)
        return 'new Immediate()' if m <= GlobalConstants.FineTol else 'new Exp(%s)' % _java_num(1.0 / m)
    name = type(dist).__name__
    if name in ('Immediate', 'Disabled'):
        return 'new %s()' % name
    if name == 'Exp':
        return 'new Exp(%s)' % _java_num(dist.getRate())
    if name == 'Det':
        return 'new Det(%s)' % _java_num(dist._value)
    if name in ('Geometric', 'Bernoulli'):
        return 'new %s(%s)' % (name, _java_num(dist._p))
    if name == 'Poisson':
        return 'new Poisson(%s)' % _java_num(dist._lambda)
    if name == 'Erlang':
        return 'new Erlang(%s, %d)' % (_java_num(dist._phase_rate), int(round(dist._phases)))
    if name == 'HyperExp':
        if getattr(dist, '_scalar_form', False):
            return 'new HyperExp(%s, %s, %s)' % (_java_num(dist._probs[0]), _java_num(dist._rates[0]),
                                                 _java_num(dist._rates[1]))
        return 'new HyperExp(%s, %s)' % (_java_doubles(dist._probs), _java_doubles(dist._rates))
    if name == 'Gamma':
        return 'new Gamma(%s, %s)' % (_java_num(dist._shape), _java_num(dist._scale))
    if name == 'Lognormal':
        return 'new Lognormal(%s, %s)' % (_java_num(dist._mu), _java_num(dist._sigma))
    if name == 'Uniform':
        return 'new Uniform(%s, %s)' % (_java_num(dist._min), _java_num(dist._max))
    if name == 'Pareto':
        return 'new Pareto(%s, %s)' % (_java_num(dist._alpha), _java_num(dist._scale))
    if name == 'Normal':
        return 'new Normal(%s, %s)' % (_java_num(dist._mean_val), _java_num(dist._std))
    if name == 'DiscreteUniform':
        return 'new DiscreteUniform(%s, %s)' % (_java_num(dist._a), _java_num(dist._b))
    if name == 'Binomial':
        return 'new Binomial(%d, %s)' % (int(round(dist._n)), _java_num(dist._p))
    if name == 'Weibull':
        return 'new Weibull(%s, %s)' % (_java_num(dist._shape), _java_num(dist._scale))
    if name in ('Coxian', 'Cox2'):
        mu, phi = _coxian_params(dist)
        return 'new Coxian(%s, %s)' % (_java_matrix(mu), _java_matrix(phi))
    if name == 'PH':
        return 'new PH(%s, %s)' % (_java_matrix(dist._alpha), _java_matrix(dist._T))
    if name == 'APH':
        return 'new APH(%s, %s)' % (_java_matrix(dist.getInitProb()), _java_matrix(dist._T))
    if name == 'MAP':
        return 'new MAP(%s, %s)' % (_java_matrix(dist._D0), _java_matrix(dist._D1))
    if name == 'ME':
        return 'new ME(%s, %s)' % (_java_matrix(dist._alpha), _java_matrix(dist._A))
    if name == 'RAP':
        return 'new RAP(%s, %s)' % (_java_matrix(dist._H0), _java_matrix(dist._H1))
    if name == 'MMPP2':
        return 'new MMPP2(%s, %s, %s, %s)' % (_java_num(dist._lambda0), _java_num(dist._lambda1),
                                              _java_num(dist._sigma0), _java_num(dist._sigma1))
    if name == 'Zipf':
        return 'new Zipf(%s, %d)' % (_java_num(dist._s), int(round(dist._n)))
    if name == 'DiscreteSampler':
        return 'new DiscreteSampler(%s, %s)' % (_java_matrix(dist._probs), _java_matrix(dist._values))
    if name in ('Replayer', 'Trace'):
        return 'new %s(%s)' % (name, _java_str(_replayer_file(dist, 'lqn2java')))
    line_warning('lqn2java', 'lqn2java writes the %s distribution of %s as an APH fitted to its mean and SCV.', name, what)
    return 'APH.fitMeanAndSCV(%s, %s)' % (_java_num(dist.getMean()), _java_num(dist.getSCV()))


def _java_class_name(model_name: str) -> str:
    """MATLAB's TestSolver<upper(name without spaces)>, reduced to a valid Java identifier as the JAR does."""
    up = str(model_name).replace(' ', '').upper()
    return 'TestSolver' + ''.join(ch for ch in up if ch.isalnum() or ch in '_$')


def _java_sched(sched) -> str:
    """The constant MATLAB's SchedStrategy.toFeature names: FCFSPRIO is spelled HOL, every other policy by its name."""
    n = _sched_name(sched)
    return 'HOL' if n == 'FCFSPRIO' else n


def _call_count_mean(m) -> float:
    """Mean of the call-count law getStruct builds from a mean m (Immediate, Bernoulli(m) below 1, Geometric(1/m) above),
    which is the value LQN2JAVA.m prints: 1/(1/m) need not be m."""
    from ...constants import GlobalConstants
    m = float(m)
    if np.isnan(m) or m <= GlobalConstants.FineTol:
        return 0.0
    if m < 1:
        return m
    return 1.0 / (1.0 / m)


def _is_ref(task) -> bool:
    return _sched_name(task.sched_strategy) == 'REF'


def _lqn_java_precedence_groups(model, tasks, acts, who):
    """Every declared precedence of TASKS grouped as LQN2JAVA groups them (serial, loop, or-fork, and-fork, or-join,
    and-join), each group in MATLAB's order: by activity index, joins by descending index of the joined activity.

    A serial chain is its consecutive pairs. Returns seven lists of (task number, kind, pre, post, param)."""
    aidx = {id(a): i for i, a in enumerate(acts)}
    abyname = {a.name: a for a in acts}

    def objs(xs):  # a precedence may name its activities or hold them
        return [abyname[x] if isinstance(x, str) else x for x in xs]

    groups = [[] for _ in range(6)]
    for t, tk in enumerate(tasks, 1):
        for prec in tk.precedences:
            kname = _type_name(model._get_prec_type(prec))
            pre = objs(model._get_prec_pre_activities(prec))
            post = objs(model._get_prec_post_activities(prec))
            if kname == 'SERIAL':
                body = objs(model._get_prec_activities(prec))
                for x, y in zip(body[:-1], body[1:]):
                    groups[0].append((aidx[id(x)], (t, 'SERIAL', [x], [y], None)))
            elif kname == 'LOOP':
                body = objs(model._get_prec_activities(prec))
                groups[1].append((aidx[id(pre[0])], (t, 'LOOP', pre, body, model._get_prec_count(prec))))
            elif kname == 'CHOICE' and len(pre) == 1 and len(post) > 1:
                groups[2].append((aidx[id(pre[0])],
                                  (t, 'ORFORK', pre, post, list(model._get_prec_probabilities(prec)))))
            elif kname == 'PARALLEL' and len(pre) == 1 and len(post) > 1:
                groups[3].append((aidx[id(pre[0])], (t, 'ANDFORK', pre, post, None)))
            elif kname == 'CHOICE' and len(post) == 1:
                groups[4].append((-aidx[id(post[0])], (t, 'ORJOIN', pre, post, None)))
            elif kname == 'PARALLEL' and len(post) == 1:
                groups[5].append((-aidx[id(post[0])], (t, 'ANDJOIN', pre, post, getattr(prec, 'pre_params', None))))
            else:
                raise ValueError("%s: a %s precedence of task '%s' cannot be written as source yet."
                                 % (who, kname, tk.name))
    return [[item for _, item in sorted(g, key=lambda kv: kv[0])] for g in groups]


def _generate_lqn_java_code(model: Any, model_name: str, output: TextIO) -> None:
    """Write the Java program lqn2java describes, line for line as LQN2JAVA.m writes it.

    The layout, the order of every declaration and the spelling of every law and number are MATLAB's; where the JAR
    port deliberately departs from MATLAB (call means, loop counts, or-fork probabilities and quorums printed exactly
    where %g would round, precedences read from the declared list), this follows the JAR. Features outside what
    LQN2JAVA.m writes (processor quantum and speed factor, task priority and fan-in/out, open arrivals, forwarding,
    phases, call groups), which the JAR refuses, are written as extra setter lines, so a model without them gives the
    MATLAB text and one with them still rebuilds exactly.
    """
    from ...constants import GlobalConstants
    procs = list(model.processors)
    tasks = list(model.tasks)
    entries = [e for t in tasks for e in t.entries]
    acts = [a for t in tasks for a in t.activities]
    sn = model.getStruct()
    mult = np.ravel(np.asarray(sn.mult, dtype=float))
    pnum = {id(p): h + 1 for h, p in enumerate(procs)}
    tnum = {id(t): i + 1 for i, t in enumerate(tasks)}
    enum = {id(e): i + 1 for i, e in enumerate(entries)}
    ebyname = {e.name: e for e in entries}

    def refuse(obj, what):
        raise ValueError("lqn2java: %s '%s' carries %s, which cannot be written as source yet."
                         % (type(obj).__name__, obj.name, what))

    def check_server(obj):
        if hasattr(obj, 'hasLinearConstraints') and obj.hasLinearConstraints():
            refuse(obj, 'an admission constraint')
        if hasattr(obj, 'hasRateDependence') and obj.hasRateDependence():
            refuse(obj, 'a rate dependence')
        if hasattr(obj, 'hasServerPools') and obj.hasServerPools():
            refuse(obj, 'server pools')

    w = output.write
    w("package jline.examples;\n\n")
    w("import java.util.ArrayList;\n")
    w("import jline.lang.*;\n")
    w("import jline.lang.layered.*;\n")
    w("import jline.lang.constant.*;\n")
    w("import jline.lang.processes.*;\n")
    w("import jline.util.matrix.Matrix;\n")
    w("import jline.solvers.ln.SolverLN;\n\n")
    w("public class %s {\n\n" % _java_class_name(model_name))
    w("\tpublic static void main(String[] args) throws Exception{\n\n")
    w("\tLayeredNetwork model = new LayeredNetwork(%s);\n" % _java_str(model_name))
    w("\n")

    for h, p in enumerate(procs, 1):
        if type(p).__name__ != 'Processor':
            refuse(p, 'a %s host' % type(p).__name__)
        check_server(p)
        w("\tProcessor P%d = new Processor(model, %s, %s, SchedStrategy.%s);\n"
          % (h, _java_str(p.name), _java_count(mult[sn.hshift + h - 1]), _java_sched(p.sched_strategy)))
        if p.getReplication() != 1:
            w("P%d.setReplication(%d);\n" % (h, int(p.getReplication())))
        if p.getQuantum() not in (0.0, 0.001):
            w("\tP%d.setQuantum(%s);\n" % (h, _java_num(p.getQuantum())))
        if p.getSpeedFactor() != 1.0:
            w("\tP%d.setSpeedFactor(%s);\n" % (h, _java_num(p.getSpeedFactor())))
    w("\n")

    for i, t in enumerate(tasks, 1):
        if type(t).__name__ != 'Task':
            refuse(t, 'a %s' % type(t).__name__)
        check_server(t)
        if t.setup_time is not None or t.delay_off_time is not None:
            refuse(t, 'a setup or delay-off time')
        w("\tTask T%d = new Task(model, %s, %s, SchedStrategy.%s).on(P%d);\n"
          % (i, _java_str(t.name), _java_count(mult[sn.tshift + i - 1]), _java_sched(t.sched_strategy),
             pnum[id(t.processor)]))
        if t.getReplication() != 1:
            w("\tT%d.setReplication(%d);\n" % (i, int(t.getReplication())))
        if type(t.think_time).__name__ != 'Disabled':
            w("\tT%d.setThinkTime(%s);\n" % (i, _java_dist(t.think_time, 'think time of task ' + t.name)))
        if t.priority != 0:
            w("\tT%d.setPriority(%d);\n" % (i, int(t.priority)))
        for src, val in t.getFanIn().items():
            w("\tT%d.setFanIn(%s, %d);\n" % (i, _java_str(src), int(val)))
        for dst, val in t.getFanOut().items():
            w("\tT%d.setFanOut(%s, %d);\n" % (i, _java_str(dst), int(val)))
    w("\n")

    for k, e in enumerate(entries, 1):
        if type(e).__name__ != 'Entry':
            refuse(e, 'an %s' % type(e).__name__)
        w("\tEntry E%d = new Entry(model, %s).on(T%d);\n" % (k, _java_str(e.name), tnum[id(e.task)]))
        if e.getArrival() is not None:
            w("\tE%d.setArrival(%s);\n" % (k, _java_dist(e.getArrival(), 'arrival of entry ' + e.name)))
    for k, e in enumerate(entries, 1):
        for dst, prob in zip(e.getForwardingDests(), e.getForwardingProbs()):
            w("\tE%d.forward(%s, %s);\n" % (k, _java_str(dst), _java_num(prob)))
    w("\n")

    for i, a in enumerate(acts, 1):
        w("\tActivity A%d = new Activity(model, %s, %s).on(T%d);"
          % (i, _java_str(a.name), _java_dist(a.host_demand, 'host demand of activity ' + a.name), tnum[id(a.task)]))
        if a.bound_entry is not None:
            w(" A%d.boundTo(E%d);" % (i, enum[id(a.bound_entry)]))
        calls = ''
        for want in ('SYNC', 'ASYNC'):  # getStruct numbers an activity's synchronous calls before its asynchronous ones
            for entry, mean_calls, call_type in a.calls:
                if _type_name(call_type) == want:
                    target = entry if not isinstance(entry, str) else ebyname[entry]
                    calls += '.%s(E%d,%s)' % ('synchCall' if want == 'SYNC' else 'asynchCall', enum[id(target)],
                                              _java_g(_call_count_mean(mean_calls)))
        if calls:
            w(" A%d%s;" % (i, calls))
        # a reference task does not reply, nor is a reply to a reference task's entry written
        if a.reply_entry is not None and not _is_ref(a.task) and not _is_ref(a.reply_entry.task):
            w(" A%d.repliesTo(E%d);" % (i, enum[id(a.reply_entry)]))
        w("\n")
        # the Activity default think time is Immediate, so only a non-default one is written
        th = a.think_time
        if isinstance(th, (int, float, np.integer, np.floating)):
            if float(th) > GlobalConstants.FineTol:
                w("\tA%d.setThinkTime(%s);\n" % (i, _java_dist(th, 'think time of activity ' + a.name)))
        elif th is not None and type(th).__name__ not in ('Immediate', 'Disabled'):
            w("\tA%d.setThinkTime(%s);\n" % (i, _java_dist(th, 'think time of activity ' + a.name)))
        if a.phase != 1:
            w("\tA%d.setPhase(%d);\n" % (i, int(a.phase)))
        for strategy, group in a.call_groups:
            members = _java_list([_java_str(e.name) for e in group], str)
            w("\tA%d.recordCallGroup(RoutingStrategy.%s, %s);\n" % (i, _type_name(strategy), members))
    w("\n")

    aidx = {id(a): n for n, a in enumerate(acts)}

    def by_index(xs):
        return sorted(range(len(xs)), key=lambda j: aidx[id(xs[j])])

    serial, loops, orforks, andforks, orjoins, andjoins = _lqn_java_precedence_groups(model, tasks, acts, 'lqn2java')
    declared = {'precActs': False, 'postActs': False, 'probs': False}

    def declare(var, matlab_space=False):
        if not declared[var]:
            w("\tArrayList<String> %s = new ArrayList<String>();\n" % var)
            declared[var] = True
        else:
            # MATLAB's LQN2JAVA writes the and-fork reassignment with a stray space; kept for identical text
            w(("\t %s = new ArrayList<String>();\n" if matlab_space else "\t%s = new ArrayList<String>();\n") % var)

    for t, _, pre, post, _ in serial:
        w("\tT%d.addPrecedence(ActivityPrecedence.Serial(%s, %s));\n" % (t, _java_str(pre[0].name), _java_str(post[0].name)))
    for t, _, pre, body, count in loops:
        w("\n\t// Loop Activity Precedence \n")
        declare('precActs')
        for b in body:
            w("\tprecActs.add(%s);\n" % _java_str(b.name))
        w("\tT%d.addPrecedence(ActivityPrecedence.Loop(%s, precActs, Matrix.singleton(%s)));\n"
          % (t, _java_str(pre[0].name), _java_g(count)))
    for t, _, pre, post, probs in orforks:
        order = by_index(post)
        w("\n\t// OrFork Activity Precedence \n")
        declare('precActs')
        if not declared['probs']:
            w("\tMatrix probs = new Matrix(1,%d);\n" % len(post))
            declared['probs'] = True
        else:
            w("\tprobs = new Matrix(1,%d);\n" % len(post))
        for j in order:
            w("\tprecActs.add(%s);\n" % _java_str(post[j].name))
        for n, j in enumerate(order):
            w("\tprobs.set(0,%d,%s);\n" % (n, _java_g(probs[j])))
        w("\tT%d.addPrecedence(ActivityPrecedence.OrFork(%s, precActs, probs));\n" % (t, _java_str(pre[0].name)))
    for t, _, pre, post, _ in andforks:
        w("\n\t// AndFork Activity Precedence \n")
        declare('postActs', matlab_space=True)
        for j in by_index(post):
            w("\tpostActs.add(%s);\n" % _java_str(post[j].name))
        w("\tT%d.addPrecedence(ActivityPrecedence.AndFork(%s, postActs));\n" % (t, _java_str(pre[0].name)))
    for t, _, pre, post, _ in orjoins:
        w("\n\t// OrJoin Activity Precedence \n")
        declare('precActs')
        for j in by_index(pre):
            w("\tprecActs.add(%s);\n" % _java_str(pre[j].name))
        w("\tT%d.addPrecedence(ActivityPrecedence.OrJoin(precActs, %s));\n" % (t, _java_str(post[0].name)))
    for t, _, pre, post, quorum in andjoins:
        w("\n\t// AndJoin Activity Precedence \n")
        declare('precActs')
        for j in by_index(pre):
            w("\tprecActs.add(%s);\n" % _java_str(pre[j].name))
        if quorum is None or np.size(quorum) == 0:
            w("\tT%d.addPrecedence(ActivityPrecedence.AndJoin(precActs, %s));\n" % (t, _java_str(post[0].name)))
        else:
            w("\tT%d.addPrecedence(ActivityPrecedence.AndJoin(precActs, %s, Matrix.singleton(%s)));\n"
              % (t, _java_str(post[0].name), _java_g(float(np.ravel(quorum)[0]))))

    w("\n\t// Model solution \n")
    w("\tSolverLN solver = new SolverLN(model);\n")
    w("\tsolver.getEnsembleAvg();\n")
    w("\t}\n}\n")


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
    _check_line_model(model, 'line2java')
    model_name = model.getName() if hasattr(model, 'getName') else getattr(model, 'name', 'model')
    if _is_layered(model):
        return lqn2java(model, model_name, filename)
    else:
        return qn2java(model, model_name, filename)


__all__ = [
    'qn2python',
    'qn2matlab',
    'lqn2python',
    'lqn2matlab',
    'line2python',
    'line2matlab',
    'sn2python',
    'qn2java',
    'lqn2java',
    'line2java',
]
