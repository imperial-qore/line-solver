"""
LINE Model JSON Save/Load

Provides save_model() and load_model() functions for serializing LINE
Network and LayeredNetwork models to/from JSON format conforming to
the line-model.schema.json specification.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import json
import math
import numpy as np
from typing import Any, Dict, Optional, Union


# infinite-multiplicity sentinel is literal Integer.MAX_VALUE, not GlobalConstants.MaxInt; _kb/12-interfaces-and-docs.md linemodel_io.py section.
INF_MULTIPLICITY = 2 ** 31 - 1


def _mult_to_json(mult) -> int:
    """Encode a multiplicity for JSON, mapping non-finite values to the sentinel."""
    if mult is None or not np.isfinite(mult):
        return INF_MULTIPLICITY
    return int(mult)


def _mult_from_json(mult):
    """Decode a JSON multiplicity, mapping the sentinel back to infinity.

    Negative values are also treated as infinite: earlier versions of this
    writer emitted -1 for infinite-server hosts, so files in that format must
    keep loading rather than silently becoming single-server stations.
    """
    if mult is None:
        return 1
    if mult >= INF_MULTIPLICITY or mult < 0:
        return float('inf')
    return mult


# ---------------------------------------------------------------------------
# Distribution serialization
# ---------------------------------------------------------------------------

def _dist_to_json(dist) -> Optional[Dict[str, Any]]:
    """Convert a LINE distribution object to JSON-compatible dict."""
    if dist is None:
        return None

    from ..distributions.continuous import (
        Exp, Det, Erlang, HyperExp, Gamma, Lognormal, Uniform, Immediate, Disabled, Pareto,
        Weibull, Normal, Expolynomial, NHPP, MAPt, PHt
    )
    from ..distributions.discrete import (
        Zipf, DiscreteSampler, Replayer, Geometric, Binomial, Poisson, Bernoulli,
        DiscreteUniform, EmpiricalCDF
    )
    from ..distributions.markovian import (
        PH, APH, Coxian, MAP, MMPP2, ME, RAP, DMAP, MMDP2, MarkedMMPP
    )

    if isinstance(dist, NHPP):
        return {"type": "NHPP", "params": {
            "breakpoints": _to_list(dist.breakpoints),
            "rates": _to_list(dist.rates),
            "cyclic": bool(dist.cyclic)
        }}
    if isinstance(dist, MAPt):
        return {"type": "MAPt", "params": {
            "breakpoints": _to_list(dist.breakpoints),
            "D0": [_matrix_to_list(M) for M in dist.D0],
            "D1": [_matrix_to_list(M) for M in dist.D1],
            "cyclic": bool(dist.cyclic)
        }}
    if isinstance(dist, PHt):
        return {"type": "PHt", "params": {
            "breakpoints": _to_list(dist.breakpoints),
            "alpha": [_to_list(a) for a in dist.alpha],
            "S": [_matrix_to_list(M) for M in dist.S],
            "cyclic": bool(dist.cyclic)
        }}
    if isinstance(dist, Expolynomial):
        return {"type": "Expolynomial", "expolynomial": {
            "density": dist._density, "eft": dist._eft,
            "lft": dist._lft if math.isfinite(dist._lft) else "Inf"
        }}
    if isinstance(dist, Disabled):
        return {"type": "Disabled"}
    if isinstance(dist, Immediate):
        return {"type": "Immediate"}
    if isinstance(dist, Exp):
        return {"type": "Exp", "params": {"lambda": dist._rate}}
    if isinstance(dist, Det):
        return {"type": "Det", "params": {"value": dist._value}}
    if isinstance(dist, Erlang):
        return {"type": "Erlang", "params": {"lambda": dist._phase_rate, "k": int(dist._phases)}}
    if isinstance(dist, HyperExp):
        return {"type": "HyperExp", "params": {
            "p": _to_list(dist._probs), "lambda": _to_list(dist._rates)
        }}
    if isinstance(dist, Gamma):
        # JSON beta is the scale (JAR loader: new Gamma(alpha, beta) with ctor (shape, scale))
        return {"type": "Gamma", "params": {"alpha": dist._shape, "beta": dist._scale}}
    if isinstance(dist, Lognormal):
        return {"type": "Lognormal", "params": {"mu": dist._mu, "sigma": dist._sigma}}
    if isinstance(dist, Uniform):
        return {"type": "Uniform", "params": {"a": dist._min, "b": dist._max}}
    if isinstance(dist, Zipf):
        return {"type": "Zipf", "params": {"s": dist._s, "n": int(dist._n)}}
    if isinstance(dist, Pareto):
        return {"type": "Pareto", "params": {"alpha": dist._alpha, "scale": dist._scale}}
    if isinstance(dist, Weibull):
        return {"type": "Weibull", "params": {"alpha": dist._scale, "beta": dist._shape}}
    if isinstance(dist, Normal):
        return {"type": "Normal", "params": {"mu": dist._mean_val, "sigma": dist._std}}
    if isinstance(dist, Geometric):
        return {"type": "Geometric", "params": {"p": dist._p}}
    if isinstance(dist, Binomial):
        return {"type": "Binomial", "params": {"n": int(dist._n), "p": dist._p}}
    if isinstance(dist, Poisson):
        return {"type": "Poisson", "params": {"lambda": dist._lambda}}
    if isinstance(dist, Bernoulli):
        return {"type": "Bernoulli", "params": {"p": dist._p}}
    if isinstance(dist, DiscreteUniform):
        return {"type": "DiscreteUniform", "params": {"min": dist._a, "max": dist._b}}
    if isinstance(dist, EmpiricalCDF):
        return {"type": "EmpiricalCDF", "params": {
            "x": _to_list(dist._values), "F": _to_list(dist._cdf)
        }}
    if isinstance(dist, DiscreteSampler):
        return {"type": "DiscreteSampler", "params": {
            "p": _to_list(dist._probs), "x": _to_list(dist._values)
        }}
    # Replayer: save file path and mean as fallback, plus APH fit if possible
    if isinstance(dist, Replayer):
        # a Trace built from in-memory data cannot round-trip via its tempfile path; warn and fall back to the APH fit below.
        if getattr(dist, '_file_path_is_temp', False):
            from ..api.io.logging import line_warning
            line_warning('linemodel_io',
                         'Replayer/Trace was built from in-memory samples, so '
                         'its trace lives in the generated tempfile "%s". That '
                         'path will not resolve when this model is reloaded; '
                         'construct the Trace from a persistent file to make it '
                         'round-trip. Saving the fitted APH as a fallback.'
                         % dist._file_path)
        rj = {"type": "Replayer", "params": {"fileName": dist._file_path}}
        try:
            rj["params"]["mean"] = dist.getMean()
        except Exception:
            pass
        try:
            aph = dist.fit_aph()
            if aph is not None and hasattr(aph, '_alpha') and hasattr(aph, '_T'):
                rj["ph"] = {
                    "alpha": _to_list(aph._alpha),
                    "T": _matrix_to_list(aph._T)
                }
        except Exception:
            pass
        return rj
    if isinstance(dist, MMPP2):
        return {"type": "MMPP2", "params": {
            "lambda0": float(dist._lambda0),
            "lambda1": float(dist._lambda1),
            "sigma0": float(dist._sigma0),
            "sigma1": float(dist._sigma1)
        }}
    if isinstance(dist, MMDP2):
        # MMPP2 writes rates (r0/r1), matching MATLAB/JAR, not deterministic times (d0/d1) which made 'MMDP2' mean different processes per codebase.
        return {"type": "MMDP2", "params": {
            "r0": float(dist.r0), "r1": float(dist.r1),
            "sigma0": float(dist.sigma0), "sigma1": float(dist.sigma1)
        }}
    if isinstance(dist, MarkedMMPP):
        # M3PP: full matrix list plus marking-type count K (not inferable from the ctor, which accepts two layouts).
        return {"type": "MarkedMMPP", "params": {
            "D": [_matrix_to_list(m) for m in dist._process],
            "K": int(dist._num_types)
        }}
    if isinstance(dist, DMAP):
        return {"type": "DMAP", "params": {
            "D0": _matrix_to_list(dist._D0),
            "D1": _matrix_to_list(dist._D1)
        }}
    if isinstance(dist, ME):
        return {"type": "ME", "params": {
            "alpha": _to_list(dist._alpha),
            "A": _matrix_to_list(dist._A)
        }}
    if isinstance(dist, RAP):
        return {"type": "RAP", "params": {
            "H0": _matrix_to_list(dist._H0),
            "H1": _matrix_to_list(dist._H1)
        }}
    from ..distributions.markovian import MarkedMAP as _MarkedMAP, BMAP as _BMAP
    # BMAP subclasses MarkedMAP and must be tested first; see _kb/12-interfaces-and-docs.md linemodel_io.py section (isinstance ordering).
    if isinstance(dist, _BMAP):
        return {"type": "BMAP", "params": {
            "D": [_matrix_to_list(Dk) for Dk in dist._process]
        }}
    if isinstance(dist, _MarkedMAP):
        # Marked MAP: {D0, per-mark D1k}; the aggregate is rebuilt on load
        return {"type": "MMAP", "mmap": {
            "D0": _matrix_to_list(dist._process[0]),
            "D1k": [_matrix_to_list(Dk) for Dk in dist._process[1:]]
        }}
    if isinstance(dist, MAP):
        return {"type": "MAP", "map": {
            "D0": _matrix_to_list(dist._D0),
            "D1": _matrix_to_list(dist._D1)
        }}
    # Coxian before APH/PH (subclass check)
    if isinstance(dist, Coxian):
        rates = 1.0 / dist._means  # mu = 1/mean
        phi = np.append(dist._probs, 1.0)  # add implicit last phi=1.0
        return {"type": "Coxian", "params": {
            "mu": _to_list(rates),
            "phi": _to_list(phi)
        }}
    if isinstance(dist, APH):
        return {"type": "APH", "ph": {
            "alpha": _to_list(dist._alpha),
            "T": _matrix_to_list(dist._T)
        }}
    if isinstance(dist, PH):
        return {"type": "PH", "ph": {
            "alpha": _to_list(dist._alpha),
            "T": _matrix_to_list(dist._T)
        }}

    # layered-network synthesized service (bare mean/scv holder) precedes getMean() fallback; _kb/12-interfaces-and-docs.md linemodel_io.py section.
    from ..layered import Distribution as _LayeredDistribution
    if isinstance(dist, _LayeredDistribution):
        mean = float(dist.mean) if dist.mean is not None else 0.0
        scv = float(dist.scv) if dist.scv is not None else 1.0
        if abs(scv - 1.0) < 1e-9:
            rate = (1.0 / mean) if mean > 0 else 0.0
            return {"type": "Exp", "params": {"lambda": rate}}
        return {"type": "Exp", "fit": {"method": "fitMeanAndSCV",
                                       "mean": mean, "scv": scv}}

    # Prior distribution (mixture of alternatives with prior probabilities)
    from ..distributions.continuous import Prior
    if isinstance(dist, Prior):
        alts = []
        for i in range(dist.getNumAlternatives()):
            alt_json = _dist_to_json(dist.getAlternative(i))
            if alt_json is not None:
                alts.append(alt_json)
        return {"type": "Prior", "kind": "discrete", "distributions": alts,
                "probabilities": _to_list(dist.getProbabilities())}

    # unhandled distribution emits type name + first two moments (never collapsed to Exp); _kb/12-interfaces-and-docs.md linemodel_io.py section.
    from ..api.io.logging import line_warning
    real_name = getattr(dist, '_name', None) or type(dist).__name__
    try:
        mean = float(dist.getMean())
    except Exception as err:
        line_warning('linemodel_io',
                     'Distribution "%s" has no JSON representation and its mean '
                     'could not be computed (%s); it is omitted from the saved '
                     'model.' % (real_name, err))
        return None
    try:
        scv = float(dist.getSCV())
    except Exception:
        var = None
        try:
            var = float(dist.getVar())
        except Exception:
            var = None
        if var is None or mean == 0:
            line_warning('linemodel_io',
                         'Distribution "%s" has no JSON representation and no '
                         'computable SCV; saving its mean only, which the reader '
                         'will reconstruct as an exponential.' % real_name)
            return {"type": real_name, "params": {"mean": mean}}
        scv = var / (mean ** 2)
    line_warning('linemodel_io',
                 'Distribution "%s" has no JSON representation; saving its mean '
                 'and SCV only, which the reader will reconstruct as an APH '
                 'matching those two moments.' % real_name)
    return {"type": real_name, "params": {"mean": mean, "scv": scv}}


def _json_to_dist(d: Dict[str, Any]):
    """Convert a JSON dist dict back to a LINE distribution object."""
    from ..distributions.continuous import (
        Exp, Det, Erlang, HyperExp, Gamma, Lognormal, Uniform, Immediate, Disabled, Pareto,
        Weibull, Normal, Expolynomial, NHPP, MAPt, PHt
    )
    from ..distributions.discrete import (
        Zipf, DiscreteSampler, Replayer, Geometric, Binomial, Poisson, Bernoulli,
        DiscreteUniform, EmpiricalCDF
    )
    from ..distributions.markovian import (
        PH, APH, MAP, Coxian, MMPP2, ME, RAP, BMAP, DMAP, MMDP2, MarkedMMPP
    )

    dtype = d.get("type")
    params = d.get("params", {})
    fit = d.get("fit")
    ph = d.get("ph")
    mapspec = d.get("map")

    if dtype == "Expolynomial":
        ep = d.get("expolynomial", {})
        lft = float('inf') if ep.get("lft") == "Inf" else float(ep["lft"])
        return Expolynomial(ep["density"], float(ep["eft"]), lft)

    if dtype == "Disabled":
        return Disabled.getInstance()
    if dtype == "Immediate":
        return Immediate.getInstance()

    # Replayer: try file first, then APH fallback, then Exp fallback
    if dtype == "Replayer":
        file_path = params.get("fileName")
        if file_path:
            import os
            if os.path.isfile(file_path):
                return Replayer(file_path)
        # Fallback to APH if available
        if ph:
            alpha = np.array(ph["alpha"])
            T = np.array(ph["T"])
            return PH(alpha, T)
        # Fallback to Exp with stored mean
        mean = params.get("mean", 1.0)
        return Exp.fit_mean(mean)

    # Direct params
    if params:
        if dtype == "NHPP":
            # Absent 'cyclic' means cyclic, matching the constructor default.
            return NHPP(params["breakpoints"], params["rates"],
                        params.get("cyclic", True))
        if dtype == "MAPt":
            # Absent 'cyclic' means cyclic, matching the constructor default.
            return MAPt(params["breakpoints"], params["D0"], params["D1"],
                        params.get("cyclic", True))
        if dtype == "PHt":
            # Absent 'cyclic' means cyclic, matching the constructor default.
            return PHt(params["breakpoints"], params["alpha"], params["S"],
                       params.get("cyclic", True))
        if dtype == "Exp":
            rate = params.get("lambda", params.get("rate"))
            return Exp(rate)
        if dtype == "Det":
            return Det(params["value"])
        if dtype == "Erlang":
            return Erlang(params["lambda"], int(params["k"]))
        if dtype == "HyperExp":
            p = params["p"]
            lam = params["lambda"]
            if len(p) == 2:
                return HyperExp(p[0], lam[0], lam[1])
            return HyperExp(np.array(p), np.array(lam))
        if dtype == "Gamma":
            return Gamma(params["alpha"], params["beta"])
        if dtype == "Lognormal":
            return Lognormal(params["mu"], params["sigma"])
        if dtype == "Uniform":
            return Uniform(params["a"], params["b"])
        if dtype == "Zipf":
            return Zipf(params["s"], int(params["n"]))
        if dtype == "Pareto":
            return Pareto(params["alpha"], params["scale"])
        if dtype == "Weibull":
            return Weibull(params["beta"], params["alpha"])  # constructor: Weibull(shape, scale)
        if dtype == "Normal":
            return Normal(params["mu"], params["sigma"])
        if dtype == "Geometric":
            return Geometric(params["p"])
        if dtype == "Binomial":
            return Binomial(int(params["n"]), params["p"])
        if dtype == "Poisson":
            return Poisson(params["lambda"])
        if dtype == "Bernoulli":
            return Bernoulli(params["p"])
        if dtype == "DiscreteUniform":
            return DiscreteUniform(int(params["min"]), int(params["max"]))
        if dtype == "DiscreteSampler":
            pv = np.array(params["p"])
            xv = np.array(params["x"])
            return DiscreteSampler(pv, xv)
        if dtype == "Coxian":
            mu = np.array(params["mu"])
            phi = np.array(params["phi"])
            return Coxian(mu, phi)
        if dtype == "EmpiricalCDF":
            return EmpiricalCDF(np.array(params["x"]), np.array(params["F"]))
        if dtype == "ME":
            return ME(np.array(params["alpha"]), np.array(params["A"]))
        if dtype == "RAP":
            return RAP(np.array(params["H0"]), np.array(params["H1"]))
        if dtype == "DMAP":
            return DMAP(np.array(params["D0"]), np.array(params["D1"]))
        if dtype == "BMAP":
            return BMAP([np.array(m) for m in params["D"]])
        if dtype == "MMDP2":
            return MMDP2(params["r0"], params["r1"],
                         params["sigma0"], params["sigma1"])
        if dtype == "MarkedMMPP":
            return MarkedMMPP([np.array(m) for m in params["D"]],
                              int(params["K"]))

    # Prior distribution (mixture of alternatives with prior probabilities)
    if dtype == "Prior":
        from ..distributions.continuous import Prior
        # The continuous form (parameter density plus factory template) is
        # MATLAB/C++ only; the native Prior has no such constructor, and reading
        # its keys as a discrete set would build an EMPTY alternative set
        if d.get("kind", "discrete") == "continuous":
            raise NotImplementedError(
                "a continuous Prior (paramDist plus a factory template) is carried by the MATLAB "
                "and C++ codebases only; the native Prior is discrete. Re-save the model with the "
                "Prior expanded by Prior.discretize, or solve it with lang='matlab'/'cpp'")
        dist_list = d.get("distributions", [])
        prob_list = d.get("probabilities", [])
        alternatives = [_json_to_dist(dd) for dd in dist_list]
        return Prior(alternatives, prob_list)

    # PH/APH representation
    if ph:
        alpha = np.array(ph["alpha"])
        T = np.array(ph["T"])
        if dtype == "APH":
            return APH(alpha, T)
        return PH(alpha, T)

    # MMPP2 (before generic MAP)
    if dtype == "MMPP2" and params:
        return MMPP2(params["lambda0"], params["lambda1"], params["sigma0"], params["sigma1"])

    # Marked MAP: {D0, per-mark D1k} (constructor layout [D0, D11, ..., D1K])
    if dtype == "MMAP" and d.get("mmap"):
        from ..distributions.markovian import MarkedMAP as _MarkedMAP
        mmapspec = d["mmap"]
        D0 = np.array(mmapspec["D0"])
        Dk = [np.array(m) for m in mmapspec["D1k"]]
        return _MarkedMAP([D0] + Dk)

    # MAP representation
    if mapspec:
        D0 = np.array(mapspec["D0"])
        D1 = np.array(mapspec["D1"])
        return MAP(D0, D1)

    # Fit specification
    if fit:
        method = fit["method"]
        if method == "fitMean":
            mean = fit["mean"]
            if dtype == "Exp":
                return Exp.fit_mean(mean)
            if dtype == "Erlang":
                order = fit.get("order", 1)
                return Erlang.fit_mean_and_order(mean, order)
            if dtype == "Det":
                return Det(mean)
            # Default: Exp with given mean
            return Exp.fit_mean(mean)
        if method == "fitMeanAndSCV":
            mean = fit["mean"]
            scv = fit["scv"]
            if dtype == "Erlang":
                return Erlang.fit_mean_and_scv(mean, scv)
            if dtype == "HyperExp":
                return HyperExp.fit_mean_and_scv(mean, scv)
            if dtype == "Gamma":
                return Gamma.fit_mean_and_scv(mean, scv)
            if dtype == "Lognormal":
                return Lognormal.fit_mean_and_scv(mean, scv)
            return Exp.fit_mean(mean)
        if method == "fitMeanAndOrder":
            mean = fit["mean"]
            order = fit["order"]
            if dtype == "Erlang":
                return Erlang.fit_mean_and_order(mean, order)
            return Exp.fit_mean(mean)

    # read side of the writer fallback above: an unrecognised type with mean+scv is reconstructed as a matching APH.
    from ..api.io.logging import line_warning
    if params and "mean" in params:
        mean = float(params["mean"])
        if "scv" in params:
            scv = float(params["scv"])
            line_warning('linemodel_io',
                         'Unrecognised distribution type "%s"; reconstructing an '
                         'APH matching its saved mean (%g) and SCV (%g).'
                         % (dtype, mean, scv))
            return APH.fit_mean_and_scv(mean, scv)
        line_warning('linemodel_io',
                     'Unrecognised distribution type "%s" carrying a mean (%g) '
                     'but no SCV; reconstructing an exponential with that mean.'
                     % (dtype, mean))
        return Exp.fit_mean(mean)
    line_warning('linemodel_io',
                 'Unrecognised distribution type "%s" with no mean; '
                 'reconstructing Exp(1.0). The saved model is not faithful.'
                 % (dtype,))
    return Exp(1.0)


# ---------------------------------------------------------------------------
# Network serialization
# ---------------------------------------------------------------------------

def _network_to_json(model) -> Dict[str, Any]:
    """Convert a Network model to JSON-compatible dict."""
    from ..lang.nodes import Queue, Delay, Source, Sink, Fork, Join, Router, ClassSwitch, Cache, Place, Transition
    from ..lang.classes import (
        OpenClass, ClosedClass, SelfLoopingClass, OpenSignal, ClosedSignal,
        Signal, SignalType, RemovalPolicy
    )
    from ..lang.base import SchedStrategy, NodeType, RoutingStrategy, StatefulNode

    nodes = model.get_nodes()
    classes = model.get_classes()

    # A STATE ON A STRICT SUBSET OF THE STATEFUL NODES DOES NOT TRAVEL, because
    # it is not an initialization: MATLAB reads sn.state through getState, which
    # runs initDefault whenever hasInitState is false, so the rows that are there
    # are dropped rather than combined with default markings for the rest. A
    # document that carried only the named node made the reader mix the two: on
    # Delay->PS with 2 jobs and Q1 alone set to 2, the C++ read (Think=2, Q1=2)
    # and answered getProbSysAggr 0 for a joint state holding 4 of 2 jobs, where
    # MATLAB answers 0.4 for the default marking.
    # A PAS PLACEMENT IS THE EXCEPTION, and it is MATLAB's own: initDefault.m
    # keeps the user row for a pass-and-swap station (its `hasUser` branch) while
    # defaulting every other station, because the ordering is a required input
    # there rather than a default. Dropped from the document, it reaches the
    # reader as the refusal that placement exists to answer.
    fully_initialized = model.has_init_state() if hasattr(model, 'has_init_state') else True

    # AN SPN PLACE IS THE OTHER EXCEPTION, for the same reason as PAS: its
    # INITIAL MARKING is the model, not a default any reader can reconstruct.
    # Dropped from the document, `spn_open_sevenplaces` exported a net whose
    # places all start empty -- the C++ row's JMT reported P1 throughput 1.0111
    # against 2.8376 and lost P5, P6 and P7 from the table entirely.
    def _state_travels(node):
        if fully_initialized:
            return True
        if type(node).__name__ == 'Place':
            return True
        sched = node.get_sched_strategy() if hasattr(node, 'get_sched_strategy') else None
        return sched is not None and str(getattr(sched, 'name', sched)) == 'PAS'

    # Build nodes array
    nodes_json = []
    for node in nodes:
        # Skip implicit ClassSwitch nodes (auto-created by link())
        if isinstance(node, ClassSwitch) and getattr(node, '_auto_added', False):
            continue

        nj = {"name": node.name, "type": _node_type_str(node)}

        # A Logger is defined by the file it writes: without the name it
        # reloads as a tunnel that logs nowhere, so the solver runs and the
        # trace the user asked for is silently absent.
        if _node_type_str(node) == "Logger" and getattr(node, '_file_name', None):
            nj["fileName"] = node._file_name

        # Queue scheduling always recorded, including INF (else load() defaults to FCFS); get_sched_strategy() output normalized to the enum name.
        if isinstance(node, (Queue, Delay)):
            sched = node.get_sched_strategy()
            if sched is not None:
                if not hasattr(sched, 'name'):
                    sched = SchedStrategy(sched)
                nj["scheduling"] = sched.name

        # Servers
        if isinstance(node, Queue) and not isinstance(node, Delay):
            ns = getattr(node, '_number_of_servers', 1)
            if ns is not None and np.isfinite(ns) and int(ns) != 1:
                nj["servers"] = int(ns)

        # Buffer
        if isinstance(node, Queue):
            cap = getattr(node, '_capacity', None)
            if cap is not None and np.isfinite(cap):
                nj["buffer"] = int(cap)

        # per-class buffer capacity serialized for Places too; see _kb/12-interfaces-and-docs.md linemodel_io.py section (finite-capacity Place).
        if isinstance(node, (Queue, Delay, Place)):
            cc_dict = getattr(node, '_class_capacity', {})
            if cc_dict:
                cc_json = {}
                for jc, ccap in cc_dict.items():
                    if np.isfinite(ccap):
                        cc_json[jc.name] = int(ccap)
                if cc_json:
                    nj["classCap"] = cc_json

        # drop rule is a real per-class map; broadcasting one station-wide value discards per-class differences.
        if isinstance(node, (Queue, Delay)):
            from ..lang.base import DropStrategy
            # None is the unset state (Station.__init__), so a station nobody
            # configured writes no dropRule at all -- MATLAB linemodel_save.m
            # gates the same key on `~isempty(node.dropRule)`. An explicitly set
            # rule is written whatever it is: the old `!= "drop"` filter existed
            # because every station used to arrive here preset to DROP, and it
            # silently dropped a user's own setDropRule(DROP).
            dr = getattr(node, '_drop_rule', None)
            dr_json = {}
            for jc in classes:
                rule = dr.get(jc, None) if isinstance(dr, dict) else dr
                if rule is None:
                    continue
                dr_str = _drop_strategy_to_str(rule)
                if dr_str:
                    dr_json[jc.name] = dr_str
            if dr_json:
                nj["dropRule"] = dr_json

        # Load-dependent scaling
        if isinstance(node, (Queue, Delay)):
            lld = getattr(node, '_load_depend_scaling', None)
            if lld is not None:
                # Reduced exactly as Network._refresh_load_dependence reduces it:
                # a (n, R) table written out column by column would not even be
                # a scaling vector on the far side.
                from ..lang.network import lld_scaling_as_1d
                nj["loadDependence"] = {
                    "type": "loadDependent",
                    "scaling": [float(x) for x in lld_scaling_as_1d(lld)]
                }

        # class-dependent scaling materialized over per-class box lattice since callable cannot cross JSON; mirrors MATLAB cd_scaling_table/JAR LineModelIO.
        if isinstance(node, (Queue, Delay)):
            lcd = getattr(node, '_class_depend_scaling', None)
            if lcd is not None:
                maxc = _class_dependence_cutoffs(classes)
                # Declared peak rate scaling per class (Util = T*S/peak);
                # broadcast a scalar to K so the reader restores a vector.
                pk = getattr(node, '_class_depend_scaling_peak', None)
                pk = np.atleast_1d(np.asarray(pk, dtype=float)).ravel() if pk is not None else np.array([])
                if pk.size == 1:
                    pk = np.full(len(classes), pk[0])
                nj["classDependence"] = {
                    "type": "classDependent",
                    "cutoffs": [int(c) for c in maxc],
                    "scaling": _cd_scaling_table(lcd, maxc, len(classes)),
                    "peak": [float(x) for x in pk],
                }

        # joint-dependent (non-product-form eta_i) scaling, materialized over the
        # per-class box lattice just like classDependence; mirrors the JAR
        # jointDependence wire key.
        if isinstance(node, (Queue, Delay)):
            ljd = getattr(node, '_joint_depend_scaling', None)
            if ljd is not None:
                maxc = _class_dependence_cutoffs(classes)
                pk = getattr(node, '_joint_depend_scaling_peak', None)
                pk = np.atleast_1d(np.asarray(pk, dtype=float)).ravel() if pk is not None else np.array([])
                if pk.size == 1:
                    pk = np.full(len(classes), pk[0])
                nj["jointDependence"] = {
                    "type": "jointDependent",
                    "cutoffs": [int(c) for c in maxc],
                    "scaling": _cd_scaling_table(ljd, maxc, len(classes)),
                    "peak": [float(x) for x in pk],
                }

        # Service/arrival distributions
        if isinstance(node, Source):
            svc = {}
            for jc in classes:
                d = node._arrival_process.get(jc)
                if d is not None:
                    dj = _dist_to_json(d)
                    if dj is not None:
                        svc[jc.name] = dj
            if svc:
                nj["service"] = svc
            # batch arrivals (per-class batch-size law) are separate from service, which only spaces the epochs.
            batch = {}
            for jc in classes:
                b = getattr(node, '_arrival_batch', {}).get(jc)
                if b is not None:
                    bj = _dist_to_json(b)
                    if bj is not None:
                        batch[jc.name] = bj
            if batch:
                nj["arrivalBatch"] = batch
            # Marked (MMAP) arrival binding: class names ordered by mark
            marked = getattr(node, '_marked_classes', None)
            if marked:
                nj["markedClasses"] = [c.name for c in marked]
        elif isinstance(node, (Queue, Delay)):
            svc = {}
            for jc in classes:
                d = node._service_process.get(jc)
                if d is not None:
                    dj = _dist_to_json(d)
                    if dj is not None:
                        svc[jc.name] = dj
            if svc:
                nj["service"] = svc

        # queueing place (QPN): serialize embedded-queue scheduling, server count, per-class service and departure discipline so LDES rebuilds a real QPN.
        if isinstance(node, Place) and node.is_queueing():
            psched = getattr(node, '_sched_strategy', None)
            if psched is not None:
                nj["scheduling"] = psched.name if hasattr(psched, 'name') else str(psched).upper()
            pns = getattr(node, '_number_of_servers', 1)
            if pns is not None and np.isfinite(pns) and int(pns) != 1:
                nj["servers"] = int(pns)
            svc = {}
            for jc in classes:
                d = node._service_process.get(jc)
                if d is not None:
                    dj = _dist_to_json(d)
                    if dj is not None:
                        svc[jc.name] = dj
            if svc:
                nj["service"] = svc
            dd = {}
            for jc in classes:
                disc = node._departure_discipline.get(jc)
                if disc is not None:
                    dd[jc.name] = disc.name if hasattr(disc, 'name') else str(disc).upper()
            if dd:
                nj["departureDiscipline"] = dd

        # PAS/OI queue: serialize the swap graph and a macrostate rate table for mu(c); see _kb/12-interfaces-and-docs.md linemodel_io.py section.
        if isinstance(node, Queue) and not isinstance(node, Delay):
            mu_fun = getattr(node, '_svc_rate_fun', None)
            if mu_fun is not None:
                swap = getattr(node, '_swap_graph', None)
                I = len(classes)
                if swap is not None:
                    swap = np.asarray(swap, dtype=float)
                    nj["swapGraph"] = [[float(swap[i, j]) for j in range(I)]
                                       for i in range(I)]

                def _mu_at(nvec):
                    order = []
                    for r in range(I):
                        order.extend([r] * int(nvec[r]))
                    return float(mu_fun(np.array(order, dtype=int)))

                def _counts(total, parts):       # count vectors summing to total
                    if parts == 1:
                        yield (total,)
                        return
                    for first in range(total + 1):
                        for rest in _counts(total - first, parts - 1):
                            yield (first,) + rest

                import itertools as _it

                # Detect per-class SATURATION cutoffs tau_r: the threshold beyond
                # which mu is constant in class r (mu(n)=mu(clamp(n,tau)) for all n).
                # bounded=False means mu keeps growing (e.g. an additive rate
                # sum_r n_r*beta_r), which cannot be reduced to a finite table.
                def _detect_saturation(max_cutoff=256, tol=1e-9):
                    tau = [1] * I
                    changed = True
                    while changed:
                        changed = False
                        for r in range(I):
                            others = [range(tau[j] + 1) if j != r else [tau[r]]
                                      for j in range(I)]
                            sat = True
                            for m in _it.product(*others):
                                mp = list(m); mp[r] += 1
                                if abs(_mu_at(mp) - _mu_at(m)) > tol:
                                    sat = False
                                    break
                            if not sat:
                                tau[r] += 1
                                changed = True
                                if tau[r] > max_cutoff:
                                    return tau, False
                    return tau, True

                tau, bounded = _detect_saturation()
                cap = getattr(node, '_capacity', None)
                finite_cap = cap is not None and np.isfinite(cap)

                rate_tbl = {}
                if not finite_cap:
                    # a finite PAS buffer is required for state-space generation; saturation only compacts the rate table, it cannot bound the queue length.
                    if bounded:
                        raise ValueError(
                            "PAS station '%s' has no finite buffer: LDES state-space "
                            "generation requires one. Its order-independent service "
                            "rate saturates at class counts %s, so call setCap(N) with "
                            "N comfortably above the mean occupancy." %
                            (node.name, tau))
                    raise ValueError(
                        "PAS station '%s' has an unbounded order-independent service "
                        "rate (mu does not saturate) and no finite buffer; call "
                        "setCap(N) on the station." % node.name)

                if bounded:
                    # rate table compacted over the saturation box (clamped counts), keeping it small even for a large buffer.
                    for nvec in _it.product(*[range(t + 1) for t in tau]):
                        if sum(nvec) == 0:
                            continue
                        rate_tbl[",".join(str(x) for x in nvec)] = _mu_at(nvec)
                    nj["oiServiceRate"] = rate_tbl
                    nj["oiCutoffs"] = [int(t) for t in tau]
                else:
                    # Non-saturating but finite buffer: exact simplex table over
                    # 1 <= sum(n) <= cap (states never exceed the buffer).
                    cmax = max(int(cap), 1)
                    for tot in range(1, cmax + 1):
                        for nvec in _counts(tot, I):
                            rate_tbl[",".join(str(x) for x in nvec)] = _mu_at(nvec)
                    nj["oiServiceRate"] = rate_tbl

        # Scheduling params (DPS weights, etc.)
        if isinstance(node, Queue) and hasattr(node, '_sched_param') and node._sched_param:
            sp = {}
            for jc, val in node._sched_param.items():
                if val is not None:
                    sp[jc.name] = float(val)
            if sp:
                nj["schedParams"] = sp

        # ClassSwitch matrix
        if isinstance(node, ClassSwitch):
            csm = node.get_class_switching_matrix()
            if csm is not None:
                cs_dict = {}
                for ri, rc in enumerate(classes):
                    row = {}
                    for ci, cc in enumerate(classes):
                        if ri < csm.shape[0] and ci < csm.shape[1] and csm[ri, ci] != 0:
                            row[cc.name] = float(csm[ri, ci])
                    if row:
                        cs_dict[rc.name] = row
                if cs_dict:
                    nj["classSwitchMatrix"] = cs_dict

        # Join strategy and quorum
        if isinstance(node, Join):
            js_dict = getattr(node, '_join_strategy', {})
            for jc, js in js_dict.items():
                if js is not None and hasattr(js, 'name') and js.name != 'STD':
                    # JoinStrategy.PARTIAL is canonical, so .name is already the
                    # interchange spelling shared with MATLAB.
                    nj["joinStrategy"] = js.name
                    break
            req_dict = getattr(node, '_required', {})
            for jc, req in req_dict.items():
                if req is not None and req > 0:
                    nj["joinQuorum"] = int(req)
                    break

        # Cache config uses flat top-level JAR-compatible keys (native LDES bridge reads this format); _kb/12-interfaces-and-docs.md linemodel_io.py section.
        if isinstance(node, Cache):
            nj["numItems"] = int(node._num_items)
            cap = node._item_level_cap
            if isinstance(cap, np.ndarray):
                cap_list = [int(x) for x in cap]
            elif isinstance(cap, (list, tuple)):
                cap_list = [int(x) for x in cap]
            else:
                cap_list = [int(cap)]
            rs = node._replacement_strategy
            rs_name = rs.name if hasattr(rs, 'name') else str(rs)
            # replacement policy emitted verbatim; CLIMB is rewritten to FIFO unit-capacity-list form at solve time, not at serialization time.
            nj["itemLevelCap"] = cap_list
            nj["replacementStrategy"] = rs_name
            nj["admissionProb"] = float(getattr(node, '_admission_prob', 1.0))
            # per-item storage costs and per-list cost caps (ton21cache Sec. IX)
            isz = getattr(node, '_item_size', None)
            if isz is not None:
                nj["itemSizes"] = [int(x) for x in np.asarray(isz).ravel()]
            ccap = getattr(node, '_cost_cap', None)
            if ccap is not None:
                ccap = np.asarray(ccap).ravel()
                if getattr(node, '_cost_cap_global', False):
                    nj["costCaps"] = int(ccap[0])
                else:
                    nj["costCaps"] = [int(x) for x in ccap]

            # Hit/miss class mappings (full: includes retrieval classes)
            if node._hit_class:
                hit_map = {}
                for in_cls, out_cls in node._hit_class.items():
                    hit_map[in_cls.name] = out_cls.name
                if hit_map:
                    nj["hitClass"] = hit_map
            if node._miss_class:
                miss_map = {}
                for in_cls, out_cls in node._miss_class.items():
                    miss_map[in_cls.name] = out_cls.name
                if miss_map:
                    nj["missClass"] = miss_map

            if getattr(node, "_item_of_class", None):
                item_map = {}
                for jc, item in node._item_of_class.items():
                    if item:
                        item_map[jc.name] = int(item)
                if item_map:
                    nj["itemClass"] = item_map

            # Read popularity distributions (set_read), one per class
            if node._read_process:
                pop_map = {}
                for jc, dist in node._read_process.items():
                    dj = _dist_to_json(dist)
                    if dj is not None:
                        pop_map[jc.name] = dj
                if pop_map:
                    nj["popularity"] = pop_map

            # Access-cost (list-move) structure: per-item graph shared by all
            # classes (set_access_graph), or full per-class accost matrices.
            if getattr(node, '_graph', None) is not None:
                nj["accessGraph"] = [np.asarray(g, dtype=float).tolist()
                                     for g in node._graph]
            elif getattr(node, '_accost', None) is not None:
                nj["accessProb"] = [[np.asarray(m, dtype=float).tolist() for m in row]
                                    for row in node._accost]

            # Initial cache state [class counts | contents | retrieval bitmap]
            cache_state = node.get_state()
            if cache_state is not None and np.asarray(cache_state).size > 0:
                nj["initialState"] = [float(x) for x in np.asarray(cache_state).ravel()]

            # retrieval-system (delayed-hit cache) bookkeeping needed by the simulator to start a fetch on a miss.
            if getattr(node, '_retrieval_system_capacity', 0) > 0:
                by_class = {}
                for jc_idx0, q_indices in node._retrieval_system_queue_indices.items():
                    jc_name = classes[jc_idx0].name
                    by_class.setdefault(jc_name, {})["queues"] = \
                        [nodes[qi].name for qi in q_indices]
                for (item, in_cls), out_cls in node._retrieval_classes.items():
                    items_map = by_class.setdefault(in_cls.name, {}).setdefault("items", {})
                    items_map[str(int(item))] = out_cls.name
                nj["retrievalSystem"] = {
                    "capacity": int(node._retrieval_system_capacity),
                    "byClass": by_class,
                }

        # Setup / delay-off (server vacation). Emitted as a pair keyed by class
        # name, since set_delay_off requires both distributions on reload.
        if isinstance(node, Queue) and not isinstance(node, Delay):
            su_dict = getattr(node, '_setup_time', None)
            doff_dict = getattr(node, '_delay_off_time', None)
            if su_dict and doff_dict:
                su_json = {}
                doff_json = {}
                for jobclass, su_dist in su_dict.items():
                    doff_dist = doff_dict.get(jobclass)
                    if doff_dist is None:
                        continue
                    suj = _dist_to_json(su_dist)
                    doffj = _dist_to_json(doff_dist)
                    if suj is not None and doffj is not None:
                        su_json[jobclass.name] = suj
                        doff_json[jobclass.name] = doffj
                if su_json:
                    nj["setupTime"] = su_json
                    nj["delayOffTime"] = doff_json

        # Server breakdown/repair. The degraded down-server service is written
        # per class, which flattens the class-independent form set_breakdown
        # also accepts: both rebuild the same sn.downServiceRates row.
        if isinstance(node, Queue) and not isinstance(node, Delay) and node.has_breakdown():
            fdist, rdist, dsvc_list = node.get_breakdown()
            fj = _dist_to_json(fdist)
            rj = _dist_to_json(rdist)
            if fj is None or rj is None:
                raise ValueError(
                    'Station "%s" has a breakdown whose failure or repair distribution '
                    'cannot be serialized.' % node.get_name())
            bd = {"failure": fj, "repair": rj}
            down_json = {}
            for r, jobclass in enumerate(classes):
                dsvc = None
                if len(dsvc_list) == 1:
                    dsvc = dsvc_list[0]
                elif len(dsvc_list) > r:
                    dsvc = dsvc_list[r]
                if dsvc is None or type(dsvc).__name__ == 'Disabled':
                    continue
                dj = _dist_to_json(dsvc)
                if dj is not None:
                    down_json[jobclass.name] = dj
            if down_json:
                bd["downService"] = down_json
            nj["breakdown"] = bd

        # polling type/switchover keys written by NAME (python auto()-assigns ids differently from MATLAB/JAR).
        if isinstance(node, Queue) and not isinstance(node, Delay):
            from ..constants import PollingType
            polling_type = getattr(node, '_polling_type', None)
            if polling_type is not None:
                nj["pollingType"] = polling_type.name
                if polling_type == PollingType.KLIMITED:
                    nj["pollingPar"] = int(getattr(node, '_polling_k', 1))

            # under POLLING the switchover key is the departing class alone; otherwise it is the (from,to) pair.
            so_dict = getattr(node, '_switchover', None)
            if so_dict:
                so_list = []
                for key, dist in so_dict.items():
                    dj = _dist_to_json(dist)
                    if dj is None:
                        continue
                    if isinstance(key, tuple):
                        from_cls, to_cls = key
                        so_list.append({
                            "from": from_cls.name,
                            "to": to_cls.name,
                            "distribution": dj
                        })
                    else:
                        so_list.append({
                            "from": key.name,
                            "distribution": dj
                        })
                if so_list:
                    nj["switchoverTimes"] = so_list

        # Heterogeneous server types
        if isinstance(node, Queue) and hasattr(node, 'is_heterogeneous') and node.is_heterogeneous():
            st_arr = []
            for st in node.get_server_types():
                stj = {"name": st.name, "count": st.get_num_of_servers()}
                # Compatible classes
                cc_list = st.get_compatible_classes()
                if cc_list:
                    stj["compatibleClasses"] = [c.name for c in cc_list]
                # Per-class service distributions
                svc = {}
                for jc in classes:
                    dist = node.get_hetero_service(jc, st)
                    if dist is not None:
                        dj = _dist_to_json(dist)
                        if dj is not None:
                            svc[jc.name] = dj
                if svc:
                    stj["service"] = svc
                st_arr.append(stj)
            if st_arr:
                nj["serverTypes"] = st_arr
                hsp = node.get_hetero_sched_policy()
                if hsp is not None:
                    from ..lang.base import HeteroSchedPolicy
                    if hsp != HeteroSchedPolicy.ORDER:
                        nj["heteroSchedPolicy"] = hsp.to_text()

        # Balking
        if isinstance(node, (Queue, Delay)):
            from ..lang.base import BalkingStrategy
            balk_json = {}
            for jc in classes:
                if node.has_balking(jc):
                    strategy, thresholds = node.get_balking(jc)
                    bjc = {"strategy": BalkingStrategy.to_text(strategy)}
                    th_arr = []
                    for th in thresholds:
                        min_j, max_j, prob = th[0], th[1], th[2]
                        tjson = {"minJobs": int(min_j), "probability": float(prob)}
                        if max_j == float('inf') or max_j == float('inf'):
                            tjson["maxJobs"] = -1
                        else:
                            tjson["maxJobs"] = int(max_j)
                        th_arr.append(tjson)
                    bjc["thresholds"] = th_arr
                    balk_json[jc.name] = bjc
            if balk_json:
                nj["balking"] = balk_json

        # Retrial
        if isinstance(node, (Queue, Delay)):
            ret_json = {}
            for jc in classes:
                if node.has_retrial(jc):
                    delay_dist, max_attempts = node.get_retrial(jc)
                    rjc = {
                        "delay": _dist_to_json(delay_dist),
                        "maxAttempts": int(max_attempts)
                    }
                    ret_json[jc.name] = rjc
            if ret_json:
                nj["retrial"] = ret_json

        # Patience
        if isinstance(node, (Queue, Delay)):
            from ..lang.base import ImpatienceType
            pat_json = {}
            for jc in classes:
                if node.has_patience(jc):
                    pat_dist = node.get_patience(jc)
                    pjc = {"distribution": _dist_to_json(pat_dist)}
                    imp_type = node.get_impatience_type(jc)
                    if imp_type is not None:
                        pjc["impatienceType"] = ImpatienceType.to_text(imp_type)
                    pat_json[jc.name] = pjc
            if pat_json:
                nj["patience"] = pat_json

        # Orbit impatience (abandonment from the retrial orbit, distinct from
        # the queue patience above)
        if isinstance(node, (Queue, Delay)):
            orb_json = {}
            for jc in classes:
                if node.has_orbit_impatience(jc):
                    oj = _dist_to_json(node.get_orbit_impatience(jc))
                    if oj is not None:
                        orb_json[jc.name] = oj
            if orb_json:
                nj["orbitImpatience"] = orb_json

        # Batch rejection probability (retrial queues), per class
        if isinstance(node, (Queue, Delay)):
            brp_json = {}
            for jc in classes:
                brp = node.get_batch_reject_probability(jc)
                if brp > 0:
                    brp_json[jc.name] = float(brp)
            if brp_json:
                nj["batchRejectProb"] = brp_json

        # Job parallelism: servers seized at once by a job, per class
        if isinstance(node, Queue):
            par_json = {}
            for jc in classes:
                npar = node.get_server_parallelism(jc)
                if npar > 1:
                    par_json[jc.name] = int(npar)
            if par_json:
                nj["serverParallelism"] = par_json

        # Immediate feedback, per class (node-level; the class-level flag is
        # carried separately on the class object)
        if isinstance(node, Queue):
            imf_json = {}
            for jc in classes:
                if node.has_immediate_feedback(jc):
                    imf_json[jc.name] = True
            if imf_json:
                nj["immediateFeedback"] = imf_json

        # prior probability is meaningless without the matching state space; emit the pair or neither, as the JAR writer does.
        #
        # A trivial [1] prior over one row is NOT emitted here, and does not need
        # to be: `initialState` below already carries that row for every stateful
        # node, and the C++ reader spells it as exactly this pair. What this block
        # is for is a prior over SEVERAL rows, which `initialState` cannot
        # express.
        prior = getattr(node, '_state_prior', None)
        space = getattr(node, '_state_space', None)
        # A NEGATIVE ENTRY IS THE "IGNORE THIS STATION" FLAG of the getProb*
        # family, not a state: the reference sets [-1, ...] on the stations whose
        # marginal it is not asking about, and solver_nc_margaggr reports NaN for
        # them. Read as a state space it would START THE CHAIN there.
        if space is not None and np.asarray(space).size and np.min(np.asarray(space)) < 0:
            space, prior = None, None
        if not _state_travels(node):
            space, prior = None, None
        if prior is not None and np.asarray(prior).size > 1:
            prior = np.asarray(prior, dtype=float).ravel()
            space = None if space is None else np.atleast_2d(np.asarray(space, dtype=float))
            if space is None or space.shape[0] != prior.size:
                from ..api.io.logging import line_warning
                line_warning(
                    "linemodel_save",
                    "Node %s carries a state prior over %d states but a state space of "
                    "%d rows; the prior is not saved."
                    % (node.name, prior.size, 0 if space is None else space.shape[0]))
            else:
                nj["stateSpace"] = [[float(x) for x in row] for row in space]
                nj["statePrior"] = [float(x) for x in prior]

        # Fork tasksPerLink
        if isinstance(node, Fork):
            tpl = getattr(node, '_tasks_per_link', None)
            if tpl is not None and tpl > 1:
                nj["tasksPerLink"] = int(tpl)

        # Join paired fork
        if isinstance(node, Join):
            fork_ref = getattr(node, '_fork', None)
            if fork_ref is not None and hasattr(fork_ref, 'name'):
                nj["forkNode"] = fork_ref.name

        # Transition modes
        if isinstance(node, Transition):
            modes_json = []
            n_modes = node.get_number_of_modes()
            all_nodes = model.get_nodes()
            for mi_idx in range(n_modes):
                mj = {}
                if mi_idx < len(node._mode_names):
                    mj["name"] = node._mode_names[mi_idx]
                else:
                    mj["name"] = f"Mode{mi_idx + 1}"
                # Distribution
                if mi_idx < len(node._distributions) and node._distributions[mi_idx] is not None:
                    dj = _dist_to_json(node._distributions[mi_idx])
                    if dj is not None:
                        mj["distribution"] = dj
                # Timing strategy
                if mi_idx < len(node._timing_strategies):
                    from ..lang.nodes import TimingStrategy
                    if node._timing_strategies[mi_idx] == TimingStrategy.IMMEDIATE:
                        mj["timingStrategy"] = "IMMEDIATE"
                    else:
                        mj["timingStrategy"] = "TIMED"
                # Number of servers
                if mi_idx < len(node._number_of_servers) and node._number_of_servers[mi_idx] > 1:
                    ns_val = node._number_of_servers[mi_idx]
                    if math.isinf(ns_val):
                        mj["numServers"] = "Infinity"
                    else:
                        mj["numServers"] = int(ns_val)
                # Firing priority
                # Omit only when it equals the builder default of 1: an explicit 0 is
                # a legal JMT firing priority and has to survive the round trip (BUG-90).
                if mi_idx < len(node._firing_priorities) and node._firing_priorities[mi_idx] != 1:
                    mj["firingPriority"] = float(node._firing_priorities[mi_idx])
                # Firing weight
                if mi_idx < len(node._firing_weights) and node._firing_weights[mi_idx] != 1.0:
                    mj["firingWeight"] = float(node._firing_weights[mi_idx])
                # Marking-dependent firing-rate multiplier g_mode(marking),
                # materialized over the enabling (place,class) box lattice. Timed
                # modes only; None handle => omitted (unit multiplier).
                frm_list = getattr(node, '_firing_rate_dependence', None)
                is_immediate = (mi_idx < len(node._timing_strategies)
                                and str(getattr(node._timing_strategies[mi_idx], 'name',
                                                node._timing_strategies[mi_idx])).upper() == 'IMMEDIATE')
                if (frm_list is not None and mi_idx < len(frm_list)
                        and frm_list[mi_idx] is not None and not is_immediate):
                    g = frm_list[mi_idx]
                    ec_mat = node._enabling_conditions[mi_idx]
                    slot_idx = []
                    slots = []
                    caps = []
                    for ni in range(ec_mat.shape[0]):
                        for ci in range(ec_mat.shape[1]):
                            if ec_mat[ni, ci] > 0:
                                pcap = getattr(all_nodes[ni], '_capacity', np.inf)
                                cap = int(round(pcap)) if np.isfinite(pcap) else 10
                                slot_idx.append((ni, ci))
                                caps.append(cap)
                                slots.append({"node": all_nodes[ni].name, "class": classes[ci].name})
                    if slot_idx:
                        nnodes_all = len(all_nodes)
                        nclasses_all = len(classes)
                        scaling = {}
                        total = 1
                        for cap in caps:
                            total *= (cap + 1)
                        for li in range(total):
                            rem = li
                            counts = []
                            for cap in caps:
                                counts.append(rem % (cap + 1))
                                rem //= (cap + 1)
                            mm = np.zeros((nnodes_all, nclasses_all))
                            for s, (ni, ci) in enumerate(slot_idx):
                                mm[ni, ci] = counts[s]
                            v = float(g(mm))
                            if not np.isfinite(v):
                                v = 0.0
                            scaling[",".join(str(c) for c in counts)] = v
                        mj["firingRateDependence"] = {
                            "slots": slots,
                            "cutoffs": [int(c) for c in caps],
                            "scaling": scaling,
                        }
                # Enabling conditions
                if mi_idx < len(node._enabling_conditions):
                    ec_mat = node._enabling_conditions[mi_idx]
                    ec_list = []
                    for ni in range(ec_mat.shape[0]):
                        for ci in range(ec_mat.shape[1]):
                            if ec_mat[ni, ci] > 0:
                                ec_list.append({
                                    "node": all_nodes[ni].name,
                                    "class": classes[ci].name,
                                    "count": float(ec_mat[ni, ci])
                                })
                    if ec_list:
                        mj["enablingConditions"] = ec_list
                # Inhibiting conditions
                if mi_idx < len(node._inhibiting_conditions):
                    ic_mat = node._inhibiting_conditions[mi_idx]
                    ic_list = []
                    for ni in range(ic_mat.shape[0]):
                        for ci in range(ic_mat.shape[1]):
                            if np.isfinite(ic_mat[ni, ci]):
                                ic_list.append({
                                    "node": all_nodes[ni].name,
                                    "class": classes[ci].name,
                                    "count": float(ic_mat[ni, ci])
                                })
                    if ic_list:
                        mj["inhibitingConditions"] = ic_list
                # Firing outcomes
                if mi_idx < len(node._firing_outcomes):
                    fo_mat = node._firing_outcomes[mi_idx]
                    fo_list = []
                    for ni in range(fo_mat.shape[0]):
                        for ci in range(fo_mat.shape[1]):
                            if fo_mat[ni, ci] != 0:
                                fo_list.append({
                                    "node": all_nodes[ni].name,
                                    "class": classes[ci].name,
                                    "count": float(fo_mat[ni, ci])
                                })
                    if fo_list:
                        mj["firingOutcomes"] = fo_list
                modes_json.append(mj)
            if modes_json:
                nj["modes"] = modes_json

        # THE INITIAL STATE OF EVERY STATEFUL NODE, not only a Place's token
        # counts and a closed pass-and-swap station's ordered job placement.
        # A reader decides whether the model is initialized by testing EVERY
        # stateful node (`hasInitState`), so a document that names the one node
        # the caller moved and stays silent about the rest reads as
        # uninitialized: the JAR then ran initDefault() and discarded the very
        # state, state space and prior this file carries. init_state_fcfs_nonexp
        # came back over the bridge reporting Prior 1's numbers for all three
        # priors, with nothing in the output saying the priors had been dropped.
        # Emitted only when the node HAS a state, so a model saved before any
        # initDefault() still travels without one and the reader builds its own.
        state = getattr(node, '_state', None)
        # The same "ignore this station" flag is filtered here for the same
        # reason: it is a query argument of getProb*, not a state, and a reader
        # that takes it for one starts the chain in a state the model has not.
        if isinstance(node, StatefulNode) and state is not None and _state_travels(node) \
                and np.asarray(state).size > 0 and np.min(np.asarray(state)) >= 0:
            nj["initialState"] = [float(x) for x in np.atleast_1d(np.asarray(state)).ravel()]

        nodes_json.append(nj)

    # Build classes array
    classes_json = []
    for jc in classes:
        cj = {"name": jc.name}
        if isinstance(jc, Signal) and not isinstance(jc, (OpenSignal, ClosedSignal)):
            # unresolved Signal placeholder must be resolved before the terminal else branch; see _kb/12-interfaces-and-docs.md linemodel_io.py section.
            cj["type"] = "Signal"
            is_open = model.get_index_source_node() >= 0
            cj["openOrClosed"] = "Open" if is_open else "Closed"
            cj["signalType"] = (jc._signal_type.value
                                if hasattr(jc._signal_type, 'value')
                                else str(jc._signal_type))
            if not is_open:
                ref = jc._refstat
                if ref is not None:
                    cj["refNode"] = ref.name
            target = getattr(jc, '_target_job_class', None)
            if target is not None:
                cj["targetClass"] = target.name
            rem_dist = getattr(jc, '_removal_distribution', None)
            if rem_dist is not None:
                cj["removalDistribution"] = _dist_to_json(rem_dist)
            rem_pol = getattr(jc, '_removal_policy', None)
            if rem_pol is not None and rem_pol != RemovalPolicy.RANDOM:
                cj["removalPolicy"] = (rem_pol.value if hasattr(rem_pol, 'value')
                                       else str(rem_pol))
        elif isinstance(jc, OpenSignal):
            cj["type"] = "Signal"
            cj["openOrClosed"] = "Open"
            cj["signalType"] = jc._signal_type.value if hasattr(jc._signal_type, 'value') else str(jc._signal_type)
            target = getattr(jc, '_target_job_class', None)
            if target is not None:
                cj["targetClass"] = target.name
            rem_dist = getattr(jc, '_removal_distribution', None)
            if rem_dist is not None:
                cj["removalDistribution"] = _dist_to_json(rem_dist)
            rem_pol = getattr(jc, '_removal_policy', None)
            if rem_pol is not None and rem_pol != RemovalPolicy.RANDOM:
                cj["removalPolicy"] = rem_pol.value if hasattr(rem_pol, 'value') else str(rem_pol)
        elif isinstance(jc, ClosedSignal):
            cj["type"] = "Signal"
            cj["openOrClosed"] = "Closed"
            cj["signalType"] = jc._signal_type.value if hasattr(jc._signal_type, 'value') else str(jc._signal_type)
            ref = jc._refstat
            if ref is not None:
                cj["refNode"] = ref.name
            target = getattr(jc, '_target_job_class', None)
            if target is not None:
                cj["targetClass"] = target.name
            rem_dist = getattr(jc, '_removal_distribution', None)
            if rem_dist is not None:
                cj["removalDistribution"] = _dist_to_json(rem_dist)
            rem_pol = getattr(jc, '_removal_policy', None)
            if rem_pol is not None and rem_pol != RemovalPolicy.RANDOM:
                cj["removalPolicy"] = rem_pol.value if hasattr(rem_pol, 'value') else str(rem_pol)
        elif isinstance(jc, OpenClass):
            cj["type"] = "Open"
            # An explicitly overridden reference station is real state, so it is
            # carried for open classes too (it defaults to the Source).
            ref = jc._refstat
            if ref is not None:
                cj["refNode"] = ref.name
        elif isinstance(jc, SelfLoopingClass):
            # SelfLoopingClass subclasses ClosedClass and must be tested first; see _kb/12-interfaces-and-docs.md linemodel_io.py section (isinstance ordering).
            cj["type"] = "SelfLooping"
            cj["population"] = int(jc._njobs)
            ref = jc._refstat
            if ref is not None:
                cj["refNode"] = ref.name
        elif isinstance(jc, ClosedClass):
            cj["type"] = "Closed"
            cj["population"] = int(jc._njobs)
            ref = jc._refstat
            if ref is not None:
                cj["refNode"] = ref.name
        else:
            cj["type"] = "Open"

        prio = getattr(jc, '_priority', 0)
        if prio != 0:
            cj["priority"] = int(prio)

        deadline = getattr(jc, '_deadline', float('inf'))
        if np.isfinite(deadline):
            cj["deadline"] = float(deadline)

        # Reference class of its chain (drives the WN residence-time denominator)
        if getattr(jc, '_is_reference_class', False):
            cj["isReferenceClass"] = True

        # Class-level immediate feedback (distinct from the node-level map)
        if getattr(jc, '_immediate_feedback', False):
            cj["immediateFeedback"] = True

        # Class-level patience, distinct from the node-scoped one; a node-scoped
        # patience overrides it.
        pat_dist = getattr(jc, '_patience_distribution', None)
        if pat_dist is not None:
            pat_json = _dist_to_json(pat_dist)
            if pat_json is not None:
                cj["patience"] = pat_json
                pat_type = getattr(jc, '_patience_type', None)
                if pat_type is not None:
                    from ..lang.base import ImpatienceType
                    cj["impatienceType"] = ImpatienceType.to_text(pat_type)

        # Reply signal class (carries sn.syncreply); without it a REPLY signal
        # is inert after a round-trip.
        reply_cls = getattr(jc, '_reply_signal_class', None)
        if reply_cls is not None and hasattr(reply_cls, 'name'):
            cj["replySignalClass"] = reply_cls.name

        # Spawn-on-completion binding (carries sn.classspawn): the class
        # injected at the same station on each completion of this class.
        spawn_cls = getattr(jc, '_spawn_class', None)
        if spawn_cls is not None and hasattr(spawn_cls, 'name'):
            cj["spawnClass"] = spawn_cls.name

        classes_json.append(cj)

    # Build routing matrix
    routing_json = _build_routing_json(model, nodes, classes)

    # Build routing strategies (per-node, per-class).
    #
    # THE SOURCE IS THE REFRESHED STRUCT, NOT THE PER-NODE DECLARATION, exactly
    # as in MATLAB's writer (linemodel_save.m reads sn2.routing(i,r)). RAND is
    # the DEFAULT, so a node that never had setRouting called on it carries an
    # EMPTY `_routing_strategies` and reading that dict wrote no entry at all --
    # the strategy exists only once `_refresh_routing` materializes it into
    # sn.routing. A reader then had nothing but the probability matrix and wrote
    # JMT's EmpiricalStrategy where the model asks for RandomStrategy; the two
    # consume JMT's routing stream differently, so the sample path diverges
    # (sdroute_jsq: one queue 6.6% out with every other cell matching).
    #
    # The name map is MATLAB's `stratNames`, and PROB and DISABLED are absent
    # from it for two different reasons. PROB is already carried by the
    # probability matrix. DISABLED is DERIVED: `_refresh_routing` marks a
    # (node, class) pair the class never visits, every loader skips it on read,
    # and the pairs it lands on include the AUTO-ADDED class-switch nodes the
    # node loop deliberately does not emit -- so writing it produced a
    # strategies entry naming a node absent from "nodes", a dangling reference
    # in every re-entrant and cache model.
    _ROUTING_WIRE_NAMES = (RoutingStrategy.RAND, RoutingStrategy.RROBIN,
                           RoutingStrategy.WRROBIN, RoutingStrategy.JSQ,
                           RoutingStrategy.SQ, RoutingStrategy.FIRING)
    routing_strategies_json = {}
    try:
        sn_routing = np.asarray(model.getStruct().routing)
    except Exception:
        sn_routing = None
    if sn_routing is not None and sn_routing.size:
        wire = dict((int(s.value), s.name) for s in _ROUTING_WIRE_NAMES)
        node_names = [str(n) for n in model.getStruct().nodenames]
        for i, node_name in enumerate(node_names):
            if i >= sn_routing.shape[0]:
                break
            node_strats = {}
            for r, cls in enumerate(classes):
                if r >= sn_routing.shape[1]:
                    break
                name = wire.get(int(sn_routing[i, r]))
                if name is not None:
                    node_strats[cls.name] = name
            if node_strats:
                routing_strategies_json[node_name] = node_strats
    else:
        # No refreshed struct (an unlinked model): fall back to what the nodes
        # themselves declare, which is all there is to report.
        for node in nodes:
            strats = getattr(node, '_routing_strategies', {}) or {}
            node_strats = dict((jc.name, strat.name) for jc, strat in strats.items()
                               if strat in _ROUTING_WIRE_NAMES)
            if node_strats:
                routing_strategies_json[node.name] = node_strats
    # Routing weights (for WRROBIN)
    routing_weights_json = {}
    for node in nodes:
        weights = getattr(node, '_routing_weights', {})
        if weights:
            node_weights = {}
            for jc, dest_weights in weights.items():
                if dest_weights:
                    cls_weights = {}
                    for dest_node, weight in dest_weights.items():
                        dest_name = dest_node.name if hasattr(dest_node, 'name') else str(dest_node)
                        cls_weights[dest_name] = float(weight)
                    if cls_weights:
                        node_weights[jc.name] = cls_weights
            if node_weights:
                routing_weights_json[node.name] = node_weights
    # routing-strategy parameters (SQ k) must be serialized alongside the strategy name; see _kb/12-interfaces-and-docs.md linemodel_io.py section.
    routing_params_json = {}
    for node in nodes:
        strats = getattr(node, '_routing_strategies', {})
        if not strats:
            continue
        node_params = {}
        for jc, strat in strats.items():
            if strat is None:
                continue
            sv = strat.value if hasattr(strat, 'value') else int(strat)
            params = getattr(node, '_routing_params', {}).get(jc, ())
            if not isinstance(params, tuple):
                params = (params,)
            if sv == RoutingStrategy.SQ.value:
                rp = {}
                if len(params) >= 1 and params[0] is not None:
                    rp["d"] = int(params[0])
                if rp:
                    node_params[jc.name] = rp
        if node_params:
            routing_params_json[node.name] = node_params

    # Build finite capacity regions
    fcr_json = []
    regions = getattr(model, '_regions', None) or getattr(model, 'regions', [])
    if regions:
        for region in regions:
            rj = {"name": getattr(region, 'name', 'FCR')}
            region_nodes = getattr(region, '_nodes', None) or getattr(region, 'nodes', [])
            # Stations with per-class details
            stations_json = []
            if region_nodes:
                for rn in region_nodes:
                    sj = {"node": rn.name if hasattr(rn, 'name') else str(rn)}
                    # Per-class classWeight
                    cw_dict = getattr(region, '_class_weight', {})
                    if cw_dict:
                        cw_json = {}
                        for cls, w in cw_dict.items():
                            if w != 1.0:
                                cw_json[cls.name] = float(w)
                        if cw_json:
                            sj["classWeight"] = cw_json
                    # Per-class classSize
                    cs_dict = getattr(region, '_class_size', {})
                    if cs_dict:
                        cs_json = {}
                        for cls, sz in cs_dict.items():
                            if sz != 1:
                                cs_json[cls.name] = int(sz)
                        if cs_json:
                            sj["classSize"] = cs_json
                    stations_json.append(sj)
            rj["stations"] = stations_json
            max_jobs = getattr(region, '_global_max_jobs', None)
            if max_jobs is not None and max_jobs >= 0:
                rj["globalMaxJobs"] = int(max_jobs)
            max_mem = getattr(region, '_global_max_memory', None)
            if max_mem is not None and max_mem >= 0:
                rj["globalMaxMemory"] = int(max_mem)
            # classMaxJobs
            cmj_dict = getattr(region, '_class_max_jobs', {})
            if cmj_dict:
                cmj_json = {}
                for cls, cmj in cmj_dict.items():
                    if cmj >= 0:
                        cmj_json[cls.name] = int(cmj)
                if cmj_json:
                    rj["classMaxJobs"] = cmj_json
            # classMaxMemory folds to an equivalent job cap; see _kb/04-networkstruct.md Python native network.py section (JMT classMemoryConstraint).
            cmm_dict = getattr(region, '_class_max_memory', {})
            if cmm_dict:
                cmm_json = {}
                for cls, cmm in cmm_dict.items():
                    if cmm >= 0:
                        cmm_json[cls.name] = int(cmm)
                if cmm_json:
                    rj["classMaxMemory"] = cmm_json
            # dropRule (Region stores it as _drop_rule, singular)
            dr_dict = getattr(region, '_drop_rule', None) or getattr(region, '_drop_rules', {})
            if dr_dict:
                dr_json = {}
                for cls, dr in dr_dict.items():
                    dr_str = _drop_strategy_to_str(dr) if hasattr(dr, 'name') else str(dr)
                    if dr_str:
                        dr_json[cls.name] = dr_str
                if dr_json:
                    rj["dropRule"] = dr_json
            # Linear admission constraints A*x <= b (Region.set_linear_constraints).
            # Same wire shape as the MATLAB writer: A row-by-row, b as a vector.
            if getattr(region, '_constraint_A', None) is not None \
                    and getattr(region, '_constraint_b', None) is not None:
                A = np.atleast_2d(np.asarray(region._constraint_A, dtype=float))
                b = np.atleast_1d(np.asarray(region._constraint_b, dtype=float))
                rj["constraintA"] = [[float(v) for v in row] for row in A]
                rj["constraintB"] = [float(v) for v in b.ravel()]
            fcr_json.append(rj)

    result = {
        "type": "Network",
        "name": model.name,
        "nodes": nodes_json,
        "classes": classes_json,
        "routing": routing_json
    }
    # NO SIDE TABLE MAY NAME A NODE THE "nodes" ARRAY OMITS. The loop above drops
    # the auto-added ClassSwitch nodes -- link() recreates them on load -- so any
    # per-node table keyed on their names is a dangling reference, and every
    # reader is entitled to trust the key. line-cli builds a stub ClassSwitch for
    # the unknown name and then rejects the model outright ("the class-switch
    # matrix of node 'CS_Fork1_to_Queue1' is not (nclasses x nclasses)"), which is
    # how fj_cs_postfork, fj_cs_prefork and fj_cs_multi_visits lost their whole
    # C++ row. Filtering per TABLE rather than per VALUE is what makes this hold
    # for any marker a future refresh writes onto those nodes, not just DISABLED.
    emitted = set(nj["name"] for nj in nodes_json)
    routing_strategies_json = dict((k, v) for k, v in routing_strategies_json.items()
                                   if k in emitted)
    routing_weights_json = dict((k, v) for k, v in routing_weights_json.items()
                                if k in emitted)
    routing_params_json = dict((k, v) for k, v in routing_params_json.items()
                               if k in emitted)
    if routing_strategies_json:
        result["routingStrategies"] = routing_strategies_json
    if routing_weights_json:
        result["routingWeights"] = routing_weights_json
    if routing_params_json:
        result["routingParams"] = routing_params_json
    # Krzesinski state-dependent routing, carried by NODE NAME so the block is
    # language independent. Branch index 1 is the complement M-V and is written
    # as an empty list, keeping the paper's own numbering.
    # See _kb/16-state-dependent-routing.md
    from ..lang.base import RoutingStrategy as _RS
    for _nd in model._nodes:
        _strats = getattr(_nd, "_routing_strategies", None)
        if not _strats:
            continue
        _done = False
        for _jc, _st in _strats.items():
            _sv = _st.value if hasattr(_st, "value") else int(_st)
            if _sv != int(_RS.SDR):
                continue
            _decl = getattr(_nd, "_routing_params", {}).get(_jc, ())
            if not _decl or not isinstance(_decl[0], dict):
                continue
            _decl = _decl[0]
            _br = [[]]
            for _b in range(1, len(_decl["branch"])):
                _br.append([x.get_name() for x in _decl["branch"][_b]])
            result["stateDepRouting"] = {
                "entry": _nd.get_name(),
                "departure": _decl["departure"].get_name(),
                "class": _jc.get_name(),
                "branches": _br,
                "level": [int(v) for v in _decl["level"]],
                "C": [float(v) for v in _decl["C"]],
                "d": np.asarray(_decl["d"], dtype=float).tolist(),
            }
            _done = True
            break
        if _done:
            break
    if fcr_json:
        result["finiteCapacityRegions"] = fcr_json

    # Global (Whittle) dependence phi(n): materialized over the lattice of the
    # WHOLE network state, unlike the per-station classDependence/jointDependence
    # blocks. Only the (station,class) slots a class can actually occupy carry a
    # coordinate, which is what keeps the lattice finite.
    if getattr(model, 'get_global_dependence', None) is not None \
            and model.get_global_dependence() is not None:
        result["globalDependence"] = _gd_block(model)

    rewards_json = _rewards_to_json(model)
    if rewards_json:
        result["rewards"] = rewards_json

    return result


def _rewards_to_json(model) -> list:
    """Serialize the model's rewards in the declarative {name, type, node, class} form.

    Only rewards created through a Reward.* template carry the structural
    metadata needed to reproduce them. A reward defined from a bare lambda (or
    via Reward.custom) is not reproducible from JSON: warn and omit it rather
    than emit a reward that would be wrong on reload.
    """
    from ..api.io.logging import line_warning
    from ..lang.reward import RewardDescriptor

    rewards = getattr(model, '_rewards', None)
    if not rewards:
        return []

    # rewards emitted in NAME order to byte-match the MATLAB/JAR writers (JAR HashMap iteration order is not insertion order).
    rewards_json = []
    for name in sorted(rewards.keys()):
        fn = rewards[name]
        if not isinstance(fn, RewardDescriptor):
            line_warning('linemodel_io',
                         'Reward "%s" is defined by a bare callable and cannot be serialized to JSON; '
                         'it is omitted from the saved model. Use a Reward.* template '
                         '(Reward.queue_length/utilization/blocking) for a serializable reward.' % name)
            continue
        if fn.kind == 'Custom':
            line_warning('linemodel_io',
                         'Reward "%s" is a custom reward wrapping an arbitrary function and cannot be '
                         'serialized to JSON; it is omitted from the saved model.' % name)
            continue
        if fn.node is None:
            line_warning('linemodel_io',
                         'Reward "%s" of type %s has no associated node and cannot be serialized to JSON; '
                         'it is omitted from the saved model.' % (name, fn.kind))
            continue
        rj = {"name": name, "type": fn.kind, "node": fn.node.name}
        if fn.jobclass is not None:
            rj["class"] = fn.jobclass.name
        rewards_json.append(rj)
    return rewards_json


def _json_to_rewards(model, data, node_map, class_map) -> None:
    """Restore declarative rewards onto a freshly loaded Network."""
    from ..api.io.logging import line_warning
    from ..lang.reward import Reward

    for rw in data.get("rewards", []):
        name = rw.get("name")
        rtype = rw.get("type")
        if name is None or rtype is None:
            line_warning('linemodel_io', 'Ignoring a reward entry without a "name" or "type" field.')
            continue
        node_name = rw.get("node")
        node = node_map.get(node_name) if node_name is not None else None
        if node is None:
            line_warning('linemodel_io',
                         'Reward "%s" refers to node "%s", which is not defined in this model; '
                         'the reward is ignored.' % (name, node_name))
            continue
        jobclass = None
        class_name = rw.get("class")
        if class_name is not None:
            jobclass = class_map.get(class_name)
            if jobclass is None:
                line_warning('linemodel_io',
                             'Reward "%s" refers to class "%s", which is not defined in this model; '
                             'the reward is ignored.' % (name, class_name))
                continue
        if rtype == 'QLen':
            model.set_reward(name, Reward.queue_length(node, jobclass))
        elif rtype == 'Util':
            model.set_reward(name, Reward.utilization(node, jobclass))
        elif rtype == 'Blocking':
            model.set_reward(name, Reward.blocking(node))
        else:
            line_warning('linemodel_io',
                         'Reward "%s" has type "%s", for which no reward template is implemented; '
                         'the reward is ignored.' % (name, rtype))


def _build_routing_json(model, nodes, classes) -> Dict[str, Any]:
    """Build routing JSON from the model's routing matrix."""
    matrix = {}

    # Use sn.rtnodes as primary source - it has computed routing for ALL
    # routing strategies (PROB, RROBIN, WRROBIN, etc.)
    try:
        sn = model.getStruct()
        if sn is not None and hasattr(sn, 'rtnodes') and sn.rtnodes is not None:
            rt = sn.rtnodes
            K = len(classes)
            N = len(nodes)
            # Identify auto-added ClassSwitch nodes to collapse routing through them
            auto_cs_indices = set()
            from ..lang.nodes import ClassSwitch as ClassSwitchNode
            from ..lang.nodes import Source as SourceNode, Sink as SinkNode
            from ..lang.nodes import Cache as CacheNode
            sink_indices = set()
            source_indices = set()
            cache_indices = set()
            for i, nd in enumerate(nodes):
                if isinstance(nd, ClassSwitchNode) and getattr(nd, '_auto_added', False):
                    auto_cs_indices.add(i)
                if isinstance(nd, SinkNode):
                    sink_indices.add(i)
                if isinstance(nd, SourceNode):
                    source_indices.add(i)
                if isinstance(nd, CacheNode):
                    cache_indices.add(i)

            # collect (cache,read,hit/miss) tuples to exclude from routing: Cache reconstructs this internally from its own properties.
            cache_internal_cs = set()
            for ci in cache_indices:
                cnode = nodes[ci]
                for attr in ('_hit_class', '_miss_class'):
                    mapping = getattr(cnode, attr, None)
                    if mapping:
                        for read_cls, hm_cls in mapping.items():
                            r_idx = classes.index(read_cls) if read_cls in classes else -1
                            s_idx = classes.index(hm_cls) if hm_cls in classes else -1
                            if r_idx >= 0 and s_idx >= 0:
                                cache_internal_cs.add((ci, r_idx, s_idx))

            # non-visiting (node,class) pairs require BOTH zero nodevisits AND DISABLED procid; nodevisits alone unreliable on fork-join class-switching models.
            import numpy as np
            from ..constants import ProcessType
            non_visiting = set()  # (node_idx, class_idx)
            if hasattr(sn, 'nodevisits') and sn.nodevisits is not None:
                # Collect all (node, class) with zero visits
                zero_visits = set()
                for i_n in range(N):
                    for r_c in range(K):
                        has_visit = False
                        for c in range(sn.nchains):
                            nv = sn.nodevisits[c]
                            if i_n < nv.shape[0] and r_c < nv.shape[1]:
                                if abs(nv[i_n, r_c]) > 1e-14:
                                    has_visit = True
                                    break
                        if not has_visit:
                            zero_visits.add((i_n, r_c))
                # non-visiting classification mirrors the write-side rule above (nodevisits alone is unreliable on fork-join class-switching models).
                _routed_pairs = set()  # (node_name, class_idx)
                rm = getattr(model, '_routing_matrix', None)
                if rm is not None and hasattr(rm, '_routes'):
                    _rt_routes = getattr(rm, '_original_routes', None) or rm._routes
                    for (cs, _cd), rd in _rt_routes.items():
                        cs_idx = classes.index(cs) if cs in classes else -1
                        cd_idx = classes.index(_cd) if _cd in classes else -1
                        for (ns, _nd), prob in rd.items():
                            if prob > 0:
                                _routed_pairs.add((ns.name, cs_idx))
                for (i_n, r_c) in zero_visits:
                    ist = sn.nodeToStation[i_n] if hasattr(sn, 'nodeToStation') else -1
                    if ist >= 0 and ist < sn.nstations:
                        try:
                            pid = sn.procid[ist, r_c]
                            if pid == ProcessType.DISABLED:
                                non_visiting.add((i_n, r_c))
                        except (IndexError, TypeError):
                            pass
                    elif ist < 0:
                        # Non-station node: only mark non-visiting if
                        # it has no explicit route in the original routing
                        nname = nodes[i_n].name if i_n < len(nodes) else ''
                        if (nname, r_c) not in _routed_pairs:
                            non_visiting.add((i_n, r_c))

            if auto_cs_indices:
                # collapse routing through auto-CS nodes: for i->CS->j, effective prob = sum over CS of rt[i,CS]*rt[CS,j].
                dim = N * K
                rt_arr = np.zeros((dim, dim))
                for a in range(dim):
                    for b in range(dim):
                        rt_arr[a, b] = rt[a, b]
                # Cache-internal class-switch edges zeroed BEFORE the auto-CS collapse; see _kb/12-interfaces-and-docs.md linemodel_io.py section.
                for (ci, r_idx, s_idx) in cache_internal_cs:
                    for j in range(N):
                        rt_arr[ci * K + r_idx, j * K + s_idx] = 0.0
                # Iteratively eliminate auto-CS nodes
                for cs_idx in auto_cs_indices:
                    for r_in in range(K):
                        cs_col = cs_idx * K + r_in
                        for i in range(N):
                            if i in auto_cs_indices:
                                continue
                            for r_src in range(K):
                                src = i * K + r_src
                                p_to_cs = rt_arr[src, cs_col]
                                if p_to_cs < 1e-14:
                                    continue
                                # Route through CS to all successors
                                for j in range(N):
                                    for s_dst in range(K):
                                        dst = j * K + s_dst
                                        p_from_cs = rt_arr[cs_col, dst]
                                        if p_from_cs < 1e-14:
                                            continue
                                        rt_arr[src, dst] += p_to_cs * p_from_cs
                                rt_arr[src, cs_col] = 0.0
                # routing extraction skips auto-CS nodes, Sink->Source loopback (DTMC artifact), and Cache-internal cross-class entries.
                for r in range(K):
                    for s in range(K):
                        from_to = {}
                        for i in range(N):
                            if i in auto_cs_indices:
                                continue
                            # Cache-internal class-switching entries (read->hit/miss) are skipped; reconstructed from Cache properties on load.
                            if (i, r, s) in cache_internal_cs:
                                continue
                            # Skip entries where source class r has zero
                            # visits at node i (DTMC artifacts)
                            if (i, r) in non_visiting:
                                continue
                            for j in range(N):
                                if j in auto_cs_indices:
                                    continue
                                # Skip cross-class Sink->Source entries (internal DTMC artifacts)
                                if r != s and i in sink_indices and j in source_indices:
                                    continue
                                val = rt_arr[i * K + r, j * K + s]
                                if val > 1e-14:
                                    ni = nodes[i].name
                                    nj = nodes[j].name
                                    if ni not in from_to:
                                        from_to[ni] = {}
                                    from_to[ni][nj] = float(val)
                        if from_to:
                            key = f"{classes[r].name},{classes[s].name}"
                            matrix[key] = from_to
            else:
                # explicit ClassSwitch nodes: undo Pcs-applied cross-class routing before writing so JSON stores same-class routing (Pcs already captures switching).
                import numpy as np
                explicit_cs_indices = set()
                explicit_cs_pcs = {}  # cs_idx -> Pcs matrix (K x K)
                for i, nd in enumerate(nodes):
                    if isinstance(nd, ClassSwitchNode) and not getattr(nd, '_auto_added', False):
                        explicit_cs_indices.add(i)
                        explicit_cs_pcs[i] = nd.get_class_switching_matrix()

                # base routing for explicit CS nodes recovered as P_route(CS,j,s) = rtnodes[CS*K+r,j*K+s] / Pcs[r,s].
                cs_base_route = {}  # (cs_idx, s, j) -> P_route
                for cs_idx in explicit_cs_indices:
                    pcs = explicit_cs_pcs[cs_idx]
                    for s in range(K):
                        # Find any source class r with Pcs[r,s] > 0
                        ref_r = None
                        for r in range(K):
                            if pcs is not None and pcs[r, s] > 1e-14:
                                ref_r = r
                                break
                        if ref_r is None:
                            continue
                        for j in range(N):
                            val = rt[cs_idx * K + ref_r, j * K + s]
                            if val > 1e-14:
                                p_route = float(val / pcs[ref_r, s])
                                cs_base_route[(cs_idx, s, j)] = p_route

                for r in range(K):
                    for s in range(K):
                        from_to = {}
                        for i in range(N):
                            if i in explicit_cs_indices:
                                # Use base routing (same-class only)
                                if r == s:
                                    for j in range(N):
                                        val = cs_base_route.get((i, s, j), 0.0)
                                        if val > 1e-14:
                                            ni = nodes[i].name
                                            nj = nodes[j].name
                                            if ni not in from_to:
                                                from_to[ni] = {}
                                            from_to[ni][nj] = val
                                # Skip cross-class entries for explicit CS
                                continue
                            # Skip cross-class entries from Cache nodes (Cache
                            # handles class-switching internally via hitClass/missClass)
                            if r != s and i in cache_indices:
                                continue
                            # with explicit CS nodes, non-CS nodes store only same-class entries (cross-class ones are CS-propagated, would trigger spurious auto-CS insertion).
                            if r != s and explicit_cs_indices:
                                continue
                            # Skip entries where source class r has zero
                            # visits at node i (DTMC artifacts)
                            if (i, r) in non_visiting:
                                continue
                            for j in range(N):
                                # Skip cross-class Sink->Source entries (DTMC artifacts)
                                if r != s and i in sink_indices and j in source_indices:
                                    continue
                                val = rt[i * K + r, j * K + s]
                                if val > 1e-14:
                                    ni = nodes[i].name
                                    nj = nodes[j].name
                                    if ni not in from_to:
                                        from_to[ni] = {}
                                    from_to[ni][nj] = float(val)
                        if from_to:
                            key = f"{classes[r].name},{classes[s].name}"
                            matrix[key] = from_to
    except Exception:
        pass

    # Fallback to _routing_matrix._routes if rtnodes approach didn't work
    if not matrix:
        rm = getattr(model, '_routing_matrix', None)
        if rm is not None and hasattr(rm, '_routes'):
            routes = rm._routes
            original_routes = getattr(rm, '_original_routes', None)
            if original_routes is not None:
                routes = original_routes
            for (cs, cd), route_dict in routes.items():
                key = f"{cs.name},{cd.name}"
                from_to = {}
                for (ns, nd), prob in route_dict.items():
                    if prob > 0:
                        if ns.name not in from_to:
                            from_to[ns.name] = {}
                        from_to[ns.name][nd.name] = prob
                if from_to:
                    matrix[key] = from_to

    return {"type": "matrix", "matrix": matrix}


def _json_to_network(data: Dict[str, Any]):
    """Reconstruct a Network from JSON data."""
    from ..lang.network import Network
    from ..lang.nodes import Queue, Delay, Source, Sink, Fork, Join, Router, ClassSwitch, Cache, Place, Transition
    from ..lang.classes import (
        OpenClass, ClosedClass, SelfLoopingClass, OpenSignal, ClosedSignal,
        SignalType, RemovalPolicy
    )
    from ..lang.routing import RoutingMatrix
    from ..lang.base import SchedStrategy, ReplacementStrategy
    from ..distributions.continuous import Disabled

    model = Network(data.get("name", "model"))

    # Create nodes
    node_map = {}
    for nd in data.get("nodes", []):
        name = nd["name"]
        ntype = nd["type"]
        node = _create_node(model, nd)
        node_map[name] = node

    # Deferred Join fork linking
    for nd in data.get("nodes", []):
        if nd.get("type") == "Join" and "forkNode" in nd:
            join_node = node_map.get(nd["name"])
            fork_node = node_map.get(nd["forkNode"])
            if join_node is not None and fork_node is not None:
                join_node._fork = fork_node

    # Create classes
    class_map = {}
    for cd in data.get("classes", []):
        cname = cd["name"]
        ctype = cd["type"]
        if ctype == "Open":
            prio = cd.get("priority", 0)
            jc = OpenClass(model, cname, prio)
            # An explicitly carried reference station overrides the default
            # (the Source); absent key means "leave the default".
            ref_name = cd.get("refNode")
            if ref_name:
                ref_node = node_map.get(ref_name)
                if ref_node is not None:
                    jc.set_reference_station(ref_node)
        elif ctype == "SelfLooping":
            pop = cd["population"]
            ref_name = cd.get("refNode")
            ref_node = node_map.get(ref_name) if ref_name else None
            prio = cd.get("priority", 0)
            jc = SelfLoopingClass(model, cname, pop, ref_node, prio)
        elif ctype == "Closed":
            pop = cd["population"]
            ref_name = cd.get("refNode")
            ref_node = node_map.get(ref_name) if ref_name else None
            prio = cd.get("priority", 0)
            jc = ClosedClass(model, cname, pop, ref_node, prio)
        elif ctype == "Signal":
            prio = cd.get("priority", 0)
            sig_type_str = cd.get("signalType", "negative")
            sig_type = SignalType(sig_type_str)
            # bare Signal placeholder resolved open/closed as MATLAB resolveSignals: closed absent a Source node, open otherwise; explicit openOrClosed honored.
            open_or_closed = cd.get("openOrClosed")
            if open_or_closed is None:
                has_source = any(isinstance(n, Source) for n in node_map.values())
                open_or_closed = "Open" if has_source else "Closed"
            if open_or_closed == "Closed":
                ref_name = cd.get("refNode")
                ref_node = node_map.get(ref_name) if ref_name else None
                if ref_node is None:
                    # placeholder without a refNode falls back to first non-Source station (zero signal population makes any valid closed reference station admissible).
                    from ..lang.base import Station
                    for n in node_map.values():
                        if isinstance(n, Station) and not isinstance(n, Source):
                            ref_node = n
                            break
                if ref_node is None:
                    raise ValueError(
                        "Reference station not found for closed signal '%s'" % cname)
                jc = ClosedSignal(model, cname, sig_type, ref_node, prio)
            else:
                jc = OpenSignal(model, cname, sig_type, prio)
            # Removal distribution
            rem_dist_json = cd.get("removalDistribution")
            if rem_dist_json is not None:
                rem_dist = _json_to_dist(rem_dist_json)
                if rem_dist is not None:
                    jc.setRemovalDistribution(rem_dist)
            # Removal policy
            rem_pol_str = cd.get("removalPolicy")
            if rem_pol_str is not None:
                jc.setRemovalPolicy(RemovalPolicy(rem_pol_str))
        else:
            jc = OpenClass(model, cname)
        deadline = cd.get("deadline")
        if deadline is not None and np.isfinite(deadline):
            jc._deadline = float(deadline)
        if cd.get("isReferenceClass"):
            jc.set_reference_class(True)
        if cd.get("immediateFeedback"):
            jc.setImmediateFeedback(True)
        # Class-level patience (node-scoped patience overrides it on read too,
        # because it is applied later, in the per-node loop below).
        pat_json = cd.get("patience")
        if pat_json is not None:
            pat_dist = _json_to_dist(pat_json)
            if pat_dist is not None:
                from ..lang.base import ImpatienceType
                imp_str = cd.get("impatienceType")
                imp_type = (ImpatienceType.from_text(imp_str) if imp_str
                            else ImpatienceType.RENEGING)
                jc.set_patience(imp_type, pat_dist)
        class_map[cname] = jc

    # Resolve signal targetClass associations
    for cd in data.get("classes", []):
        if cd.get("type") != "Signal":
            continue
        target_name = cd.get("targetClass")
        if target_name and cd["name"] in class_map and target_name in class_map:
            sig_cls = class_map[cd["name"]]
            target_cls = class_map[target_name]
            sig_cls.forJobClass(target_cls)
            # placeholder refstat resolution prefers the target class's reference station, required by the chain-level consistency check for REPLY signals.
            if isinstance(sig_cls, ClosedSignal) and not cd.get("refNode"):
                target_ref = target_cls.get_reference_station()
                if target_ref is not None:
                    sig_cls.set_reference_station(target_ref)

    # resolve replySignalClass associations once every class exists.
    for cd in data.get("classes", []):
        reply_name = cd.get("replySignalClass")
        if reply_name and cd["name"] in class_map and reply_name in class_map:
            class_map[cd["name"]].set_reply_signal_class(class_map[reply_name])

    # Resolve spawnClass associations (all classes are now created).
    for cd in data.get("classes", []):
        spawn_name = cd.get("spawnClass")
        if spawn_name and cd["name"] in class_map and spawn_name in class_map:
            class_map[cd["name"]].set_spawn_class(class_map[spawn_name])

    # Set service/arrival distributions
    for nd in data.get("nodes", []):
        name = nd["name"]
        node = node_map[name]
        # PAS/OI queues have no per-class service (set_service_rate_function only); representative per-class rates are derived from the rebuilt mu(c).
        is_oi = isinstance(node, Queue) and not isinstance(node, Delay) \
            and node.get_sched_strategy() in (SchedStrategy.PAS, SchedStrategy.OI)
        svc = {} if is_oi else nd.get("service", {})
        for cname, dist_json in svc.items():
            jc = class_map.get(cname)
            if jc is None:
                continue
            dist = _json_to_dist(dist_json)
            if dist is not None:
                if isinstance(node, Source):
                    node.set_arrival(jc, dist)
                elif isinstance(node, (Queue, Delay, Place)):
                    # a queueing Place is a Station, not a Queue, and needs explicit naming or it reloads with no service (pass-through).
                    node.set_service(jc, dist)

        # Batch arrivals, applied after set_arrival because set_arrival_batch
        # validates the batch law independently of the interarrival process.
        batch_json = nd.get("arrivalBatch")
        if batch_json and isinstance(node, Source):
            for cname, bj in batch_json.items():
                if cname not in class_map:
                    continue
                bdist = _json_to_dist(bj)
                if bdist is not None:
                    node.set_arrival_batch(class_map[cname], bdist)

        # marked (MMAP) arrival: rebind the shared MarkedMAP so mark k drives class k, overwriting the per-class copies set above.
        marked_names = nd.get("markedClasses")
        if marked_names and isinstance(node, Source):
            from ..distributions.markovian import MarkedMAP as _MarkedMAP
            marked_classes = [class_map[c] for c in marked_names if c in class_map]
            if marked_classes:
                first = node.get_arrival(marked_classes[0])
                if isinstance(first, _MarkedMAP):
                    node.set_marked_arrival(first, marked_classes)

        # PAS/OI: rebuild mu(c) from rate table + swap graph, not per-class representative fallback; _kb/12-interfaces-and-docs.md linemodel_io.py section.
        if "oiServiceRate" in nd and isinstance(node, Queue) \
                and not isinstance(node, Delay) \
                and node.get_sched_strategy() in (SchedStrategy.PAS, SchedStrategy.OI):
            K = len(class_map)
            # OI queues keep a fixed zero swap graph (set_swap_graph rejects them)
            sg = nd.get("swapGraph")
            if sg is not None and node.get_sched_strategy() == SchedStrategy.PAS:
                node.set_swap_graph(np.asarray(sg, dtype=float))
            rate_tbl = {k: float(v) for k, v in nd["oiServiceRate"].items()}
            max_rate = max(rate_tbl.values()) if rate_tbl else 0.0
            cutoffs = nd.get("oiCutoffs")

            def _mu_fun(order, _tbl=rate_tbl, _cut=cutoffs, _K=K, _max=max_rate):
                order = np.atleast_1d(np.asarray(order, dtype=int)).ravel()
                if order.size == 0:
                    return 0.0
                cnt = [0] * _K
                for cid in order:
                    if 0 <= cid < _K:
                        cnt[cid] += 1
                if _cut is not None:
                    for r in range(_K):
                        if r < len(_cut) and cnt[r] > _cut[r]:
                            cnt[r] = _cut[r]
                return _tbl.get(",".join(str(x) for x in cnt), _max)

            node.set_service(_mu_fun)

        # Per-class buffer capacity
        cc = nd.get("classCap", {})
        if cc and isinstance(node, (Queue, Delay)):
            for cname, cap_val in cc.items():
                jc = class_map.get(cname)
                if jc is not None and int(cap_val) > 0:
                    node.set_class_capacity(jc, int(cap_val))

        # drop rule applied through the per-class setter for every entry; reading only the first entry used to impose one class's rule on the whole station.
        dr = nd.get("dropRule", {})
        if dr and isinstance(node, (Queue, Delay)):
            for cname, dr_str in dr.items():
                jc = class_map.get(cname)
                if jc is not None:
                    node.set_drop_rule(jc, _str_to_drop_strategy(dr_str))

        # Load-dependent scaling
        ld = nd.get("loadDependence")
        if ld and isinstance(node, (Queue, Delay)):
            ld_type = ld.get("type", "loadDependent")
            if ld_type == "loadDependent" and "scaling" in ld:
                node.set_load_dependence(np.array(ld["scaling"]))

        # Class-dependent scaling beta_{i,r}(n): rebuild the callable from the
        # materialized lattice table written by _cd_scaling_table.
        cdep = nd.get("classDependence")
        if cdep and isinstance(node, (Queue, Delay)):
            cdep_type = cdep.get("type", "classDependent")
            if cdep_type == "classDependent" and "scaling" in cdep:
                cd_callable = _cd_table_to_callable(cdep["scaling"], cdep.get("cutoffs"),
                                                    len(class_map))
                cd_peak = cdep.get("peak")
                if cd_peak is None or len(cd_peak) == 0:
                    # Legacy JSON without an explicit peak: derive it from the
                    # handle over the population lattice (cutoffs).
                    from ..api.pfqn.conv import _cd_peak_scaling
                    _cut = cdep.get("cutoffs")
                    _nk = np.asarray(_cut, dtype=int) if _cut else np.ones(len(class_map), dtype=int)
                    cd_peak = _cd_peak_scaling(cd_callable, _nk, len(class_map))
                node.set_class_dependence(cd_callable, cd_peak)

        # Joint-dependent scaling eta_i(n) (non-product-form): rebuild the
        # callable from the materialized lattice table, twin of classDependence.
        jdep = nd.get("jointDependence")
        if jdep and isinstance(node, (Queue, Delay)):
            jdep_type = jdep.get("type", "jointDependent")
            if jdep_type == "jointDependent" and "scaling" in jdep:
                jd_callable = _cd_table_to_callable(jdep["scaling"], jdep.get("cutoffs"),
                                                    len(class_map))
                jd_peak = jdep.get("peak")
                if jd_peak is None or len(jd_peak) == 0:
                    from ..api.pfqn.conv import _cd_peak_scaling
                    _cut = jdep.get("cutoffs")
                    _nk = np.asarray(_cut, dtype=int) if _cut else np.ones(len(class_map), dtype=int)
                    jd_peak = _cd_peak_scaling(jd_callable, _nk, len(class_map))
                node.set_joint_dependence(jd_callable, jd_peak)

        # Join strategy and quorum
        if isinstance(node, Join):
            js_str = nd.get("joinStrategy")
            if js_str is not None:
                from ..lang.base import JoinStrategy
                # Map aliases from other codebases; PARTIAL and QUORUM are the
                # same member here, Quorum is the JAR's spelling.
                _js_map = {'Quorum': 'QUORUM', 'Partial': 'PARTIAL'}
                js_key = _js_map.get(js_str, js_str)
                js = getattr(JoinStrategy, js_key, JoinStrategy.STD)
                for jc in class_map.values():
                    node.set_strategy(jc, js)
            jq = nd.get("joinQuorum")
            if jq is not None:
                for jc in class_map.values():
                    node.set_required(jc, int(jq))

        # Heterogeneous server types
        st_arr = nd.get("serverTypes")
        if st_arr and isinstance(node, Queue):
            from ..lang.constant.server_type import ServerType
            for stData in st_arr:
                st_name = stData["name"]
                st_count = stData["count"]
                st = ServerType(st_name, st_count)
                # Compatible classes
                for cc_name in stData.get("compatibleClasses", []):
                    jc = class_map.get(cc_name)
                    if jc is not None:
                        st.add_compatible(jc)
                node.add_server_type(st)
                # Per-class service distributions
                for cname, dist_json in stData.get("service", {}).items():
                    jc = class_map.get(cname)
                    if jc is not None:
                        dist = _json_to_dist(dist_json)
                        if dist is not None:
                            node.set_hetero_service(jc, st, dist)
            # Scheduling policy
            hsp_str = nd.get("heteroSchedPolicy")
            if hsp_str:
                from ..lang.base import HeteroSchedPolicy
                hsp = HeteroSchedPolicy.from_text(hsp_str)
                node.set_hetero_sched_policy(hsp)

        # Balking
        balk_data = nd.get("balking")
        if balk_data and isinstance(node, (Queue, Delay)):
            from ..lang.base import BalkingStrategy
            for cname, bjc in balk_data.items():
                jc = class_map.get(cname)
                if jc is None:
                    continue
                strategy = BalkingStrategy.from_text(bjc["strategy"])
                thresholds = []
                for td in bjc.get("thresholds", []):
                    min_j = td["minJobs"]
                    max_j = td["maxJobs"]
                    if max_j < 0:
                        max_j = float('inf')
                    prob = td["probability"]
                    thresholds.append((min_j, max_j, prob))
                node.set_balking(jc, strategy, thresholds)

        # Retrial
        ret_data = nd.get("retrial")
        if ret_data and isinstance(node, (Queue, Delay)):
            for cname, rjc in ret_data.items():
                jc = class_map.get(cname)
                if jc is None:
                    continue
                delay_dist = _json_to_dist(rjc["delay"])
                max_attempts = rjc.get("maxAttempts", -1)
                if delay_dist is not None:
                    node.set_retrial(jc, delay_dist, max_attempts)

        # Patience
        pat_data = nd.get("patience")
        if pat_data and isinstance(node, (Queue, Delay)):
            from ..lang.base import ImpatienceType
            for cname, pjc in pat_data.items():
                jc = class_map.get(cname)
                if jc is None:
                    continue
                pat_dist = _json_to_dist(pjc["distribution"])
                imp_type = ImpatienceType.RENEGING
                imp_str = pjc.get("impatienceType")
                if imp_str:
                    imp_type = ImpatienceType.from_text(imp_str)
                if pat_dist is not None:
                    node.set_patience(jc, pat_dist, imp_type)

        # Orbit impatience
        orb_data = nd.get("orbitImpatience")
        if orb_data and isinstance(node, (Queue, Delay)):
            for cname, dist_json in orb_data.items():
                jc = class_map.get(cname)
                if jc is None:
                    continue
                orb_dist = _json_to_dist(dist_json)
                if orb_dist is not None:
                    node.set_orbit_impatience(jc, orb_dist)

        # Batch rejection probability (retrial queues)
        brp_data = nd.get("batchRejectProb")
        if brp_data and isinstance(node, (Queue, Delay)):
            for cname, p in brp_data.items():
                jc = class_map.get(cname)
                if jc is None:
                    continue
                node.set_batch_reject_probability(jc, float(p))

        # Job parallelism: servers seized at once by a job, per class
        par_data = nd.get("serverParallelism")
        if par_data and isinstance(node, Queue):
            for cname, npar in par_data.items():
                jc = class_map.get(cname)
                if jc is None:
                    continue
                node.set_server_parallelism(jc, int(npar))

        # Immediate feedback, per class
        imf_data = nd.get("immediateFeedback")
        if imf_data and isinstance(node, Queue):
            for cname, enabled in imf_data.items():
                jc = class_map.get(cname)
                if jc is not None and enabled:
                    node.set_immediate_feedback(jc)

        # prior probability restored after the state space, since the prior indexes its rows.
        prior = nd.get("statePrior")
        if prior is not None and hasattr(node, 'set_state_prior'):
            space = nd.get("stateSpace")
            if space is not None and hasattr(node, 'set_state_space'):
                node.set_state_space(np.atleast_2d(np.asarray(space, dtype=float)))
            node.set_state_prior(np.asarray(prior, dtype=float))

        # Departure discipline of a queueing Place's depository
        dd_data = nd.get("departureDiscipline")
        if dd_data and isinstance(node, Place):
            from ..constants import DepartureDiscipline
            for cname, disc_str in dd_data.items():
                jc = class_map.get(cname)
                if jc is None:
                    continue
                # The JAR emits "Normal"/"FIFO" and Python "NORMAL"/"FIFO", so
                # match case-insensitively on the name.
                disc = (DepartureDiscipline.FIFO
                        if str(disc_str).upper() == "FIFO"
                        else DepartureDiscipline.NORMAL)
                node.set_departure_discipline(jc, disc)

        # Set scheduling params
        sp = nd.get("schedParams", {})
        if sp and isinstance(node, Queue):
            for cname, val in sp.items():
                jc = class_map.get(cname)
                if jc is not None:
                    node._sched_param[jc] = val

        # Set class switch matrix (dict format: "classSwitchMatrix")
        csm = nd.get("classSwitchMatrix")
        if csm and isinstance(node, ClassSwitch):
            classes_list = model.get_classes()
            K = len(classes_list)
            mat = np.zeros((K, K))
            class_idx = {c.name: i for i, c in enumerate(classes_list)}
            for from_name, to_dict in csm.items():
                ri = class_idx.get(from_name, -1)
                if ri < 0:
                    continue
                for to_name, prob in to_dict.items():
                    ci = class_idx.get(to_name, -1)
                    if ci >= 0:
                        mat[ri, ci] = prob
            node.set_class_switching_matrix(mat)
        # Legacy 2D array format: "csMatrix" (from older JAR saves)
        elif not csm and isinstance(node, ClassSwitch):
            cs_arr = nd.get("csMatrix")
            if cs_arr and isinstance(cs_arr, list):
                mat = np.array(cs_arr, dtype=float)
                node.set_class_switching_matrix(mat)

        # cache hit/miss mappings support both python's nested 'cache' key and the JAR's flat top-level keys.
        if isinstance(node, Cache):
            cc = nd.get("cache", {})

            # PER KEY, NOT PER NODE. MATLAB's linemodel_save writes SOME cache
            # keys nested under 'cache' (hitClass, missClass, popularity) and
            # others flat on the node (retrievalSystem), so choosing one source
            # wholesale drops whatever the other one holds -- silently, since a
            # cache with no retrieval system is a legal model. That is how a
            # MATLAB-exported delayed-hit model arrived here as a PLAIN cache:
            # hit + miss summed to 1, the delayed-hit column was zero, and both
            # lang='python' and lang='cpp' reported it without complaint.
            class _CacheSrc(object):
                def __init__(self, nested, flat):
                    self._nested, self._flat = nested, flat

                def get(self, key, default=None):
                    v = self._nested.get(key)
                    if v is None:
                        v = self._flat.get(key)
                    return default if v is None else v

            cache_src = _CacheSrc(cc, nd)
            # Hit class mapping
            hc = cache_src.get("hitClass", {})
            for in_name, out_name in hc.items():
                in_cls = class_map.get(in_name)
                out_cls = class_map.get(out_name)
                if in_cls is not None and out_cls is not None:
                    node.set_hit_class(in_cls, out_cls)
            # Miss class mapping
            mc = cache_src.get("missClass", {})
            for in_name, out_name in mc.items():
                in_cls = class_map.get(in_name)
                out_cls = class_map.get(out_name)
                if in_cls is not None and out_cls is not None:
                    node.set_miss_class(in_cls, out_cls)
            # itemClass: for a cache network, the item each per-item class reads
            # (Cache.set_item_read_classes). Recorded directly, since popularity and the
            # hit/miss switches come from their own keys; inferring it from a one-hot
            # popularity would be ambiguous against a genuine single-item popularity.
            ic = cache_src.get("itemClass", {})
            for cname, item in ic.items():
                jc = class_map.get(cname)
                if jc is not None:
                    node.set_item_of_class(jc, int(item))
            # Popularity distributions (set_read)
            pop = cache_src.get("popularity", {})
            for cname, dist_json in pop.items():
                jc = class_map.get(cname)
                if jc is not None:
                    pop_dist = _json_to_dist(dist_json)
                    if pop_dist is not None and not isinstance(pop_dist, Disabled):
                        node.set_read(jc, pop_dist)

            # retrieval-system bookkeeping restored; retrieval classes/routing/service reconstruct generically.
            rs = cache_src.get("retrievalSystem")
            if rs:
                node._retrieval_system_capacity = int(rs.get("capacity", 0))
                for jc_name, entry in rs.get("byClass", {}).items():
                    jobin = class_map.get(jc_name)
                    if jobin is None:
                        continue
                    q_idx = []
                    for qn in entry.get("queues", []):
                        qnode = node_map.get(qn)
                        if qnode is not None:
                            q_idx.append(qnode.get_index0())
                    node._retrieval_system_queue_indices[jobin.get_index0()] = q_idx
                    for it_str, rc_name in entry.get("items", {}).items():
                        rc = class_map.get(rc_name)
                        if rc is not None:
                            node.set_retrieval_class(jobin, rc, int(it_str))
                            node._retrieval_class_indices.add(rc.get_index0())

            # Access-cost (list-move) structure
            ag = cache_src.get("accessGraph")
            if ag:
                node.set_access_graph([np.asarray(g, dtype=float) for g in ag])
            else:
                ap = cache_src.get("accessProb")
                if ap:
                    node._accost = [[np.asarray(m, dtype=float) for m in row]
                                    for row in ap]

            # Initial cache state [class counts | contents | retrieval bitmap]
            ist = cache_src.get("initialState")
            if ist:
                node.set_state(np.asarray(ist, dtype=float))

    # Configure Transition modes
    for nd in data.get("nodes", []):
        if nd.get("type") != "Transition" or "modes" not in nd:
            continue
        tnode = node_map[nd["name"]]
        for md in nd["modes"]:
            mode_name = md.get("name", "Mode")
            mode = tnode.add_mode(mode_name)
            # Distribution
            dist_data = md.get("distribution")
            if dist_data:
                dist = _json_to_dist(dist_data)
                if dist is not None:
                    tnode.set_distribution(mode, dist)
            # Timing strategy
            ts = md.get("timingStrategy")
            if ts:
                from ..lang.nodes import TimingStrategy
                if ts == "IMMEDIATE":
                    tnode.set_timing_strategy(mode, TimingStrategy.IMMEDIATE)
                else:
                    tnode.set_timing_strategy(mode, TimingStrategy.TIMED)
            # server count sentinel normalized from wire 'Infinity' onto GlobalConstants.MaxInt (a float infinity would not survive integer-valued exports).
            ns = md.get("numServers")
            if ns is not None:
                from ..constants import GlobalConstants
                if isinstance(ns, str):
                    ns = GlobalConstants.MaxInt if ns.lower() == 'infinity' else float(ns)
                elif math.isinf(ns):
                    ns = GlobalConstants.MaxInt
                if ns > 1:
                    tnode.set_number_of_servers(mode, ns)
            # Firing priority
            fp = md.get("firingPriority")
            if fp is not None:
                tnode.set_firing_priorities(mode, fp)
            # Firing weight
            fw = md.get("firingWeight")
            if fw is not None:
                tnode.set_firing_weights(mode, fw)
            # Enabling conditions
            for ec in md.get("enablingConditions", []):
                ec_node = node_map.get(ec["node"])
                ec_cls = class_map.get(ec["class"])
                if ec_node is not None and ec_cls is not None:
                    tnode.set_enabling_conditions(mode, ec_cls, ec_node, ec["count"])
            # Inhibiting conditions
            for ic in md.get("inhibitingConditions", []):
                ic_node = node_map.get(ic["node"])
                ic_cls = class_map.get(ic["class"])
                if ic_node is not None and ic_cls is not None:
                    tnode.set_inhibiting_conditions(mode, ic_cls, ic_node, ic["count"])
            # Firing outcomes
            for fo in md.get("firingOutcomes", []):
                fo_node = node_map.get(fo["node"])
                fo_cls = class_map.get(fo["class"])
                if fo_node is not None and fo_cls is not None:
                    tnode.set_firing_outcome(mode, fo_cls, fo_node, fo["count"])
            # Marking-dependent firing-rate multiplier (after enabling/timing/
            # distribution are set so the setter guard sees the final mode state).
            frm = md.get("firingRateDependence")
            if frm:
                g = _firingdep_table_to_handle(frm, node_map, class_map)
                if g is not None:
                    tnode.set_firing_rate_dependence(mode, g)

    # Build routing
    routing_data = data.get("routing", {})
    if routing_data.get("type") == "matrix":
        P = RoutingMatrix(model)
        matrix = routing_data.get("matrix", {})
        # Cache-internal cross-class routing entries excluded from load: Cache reconstructs them internally, so applying them would double-count.
        cache_internal_cs = set()
        sink_names = set()
        source_names = set()
        for nname, nobj in node_map.items():
            if isinstance(nobj, Cache):
                for attr in ('_hit_class', '_miss_class'):
                    mapping = getattr(nobj, attr, None)
                    if mapping:
                        for read_cls, hm_cls in mapping.items():
                            cache_internal_cs.add((nname, read_cls.name, hm_cls.name))
            elif isinstance(nobj, Sink):
                sink_names.add(nname)
            elif isinstance(nobj, Source):
                source_names.add(nname)
        for key, from_to in matrix.items():
            parts = key.split(",")
            if len(parts) != 2:
                continue
            cs_name, cd_name = parts[0].strip(), parts[1].strip()
            cs = class_map.get(cs_name)
            cd = class_map.get(cd_name)
            if cs is None or cd is None:
                continue
            for from_name, to_dict in from_to.items():
                ns = node_map.get(from_name)
                if ns is None:
                    continue
                for to_name, prob in to_dict.items():
                    nd = node_map.get(to_name)
                    if nd is None:
                        continue
                    # Skip Cache-internal class-switching self-loop entries
                    if cs_name != cd_name and from_name == to_name and \
                       (from_name, cs_name, cd_name) in cache_internal_cs:
                        continue
                    # Skip cross-class Sink→Source entries (DTMC artifacts
                    # that Python reconstructs internally)
                    if cs_name != cd_name and from_name in sink_names \
                       and to_name in source_names:
                        continue
                    P.set(cs, cd, ns, nd, prob)
        model.link(P)

    # Restore routing strategies (per-node, per-class)
    routing_strats = data.get("routingStrategies", {})
    routing_params = data.get("routingParams", {})
    if routing_strats:
        from ..lang.base import RoutingStrategy
        for node_name, class_strats in routing_strats.items():
            node = node_map.get(node_name)
            if node is None:
                continue
            node_params = routing_params.get(node_name, {})
            for cls_name, strat_name in class_strats.items():
                cls = class_map.get(cls_name)
                if cls is None:
                    continue
                strat = getattr(RoutingStrategy, strat_name, None)
                # PROB is matrix link(P) already applied; RAND is DECLARED not derived, dropping it wrote JMT Empirical where the model asks for Random.
                if strat is None or strat == RoutingStrategy.PROB:
                    continue
                # RoutingStrategy.DISABLED is a derived marker, skipped on load; see _kb/07-cross-language-parity.md RoutingStrategy.DISABLED is derived section.
                if strat == RoutingStrategy.DISABLED:
                    continue
                rp = node_params.get(cls_name)
                if strat == RoutingStrategy.SQ:
                    if not rp or "d" not in rp:
                        from ..api.io.logging import line_warning
                        line_warning('linemodel_io',
                                     'Node "%s" routes class "%s" by SQ '
                                     'but the model carries no routingParams.d; '
                                     'falling back to d=2.'
                                     % (node_name, cls_name))
                        node.set_routing(cls, strat)
                    else:
                        node.set_routing(cls, strat, int(rp["d"]))
                else:
                    node.set_routing(cls, strat)

    # Restore the global (Whittle) dependence phi(n) from the materialized slot
    # lattice. Slots are matched by NAME so a station or class reordering on the
    # writing side cannot silently shift a coordinate.
    gdep = data.get("globalDependence")
    if gdep and gdep.get("type", "globalDependent") == "globalDependent" and "scaling" in gdep:
        gd_callable, gd_peak, gd_cut = _gd_block_to_callable(gdep, node_map, class_map, model)
        if gd_callable is not None:
            model.set_global_dependence(gd_callable, gd_peak, gd_cut)

    # Restore Krzesinski state-dependent routing. Written by NODE NAME, so it is
    # restored after every node exists and after link(P), whose uniform
    # placeholder in the entry row this block supersedes.
    sdr_data = data.get("stateDepRouting")
    if sdr_data:
        _entry = node_map.get(sdr_data["entry"])
        _dep = node_map.get(sdr_data["departure"])
        _cls = class_map.get(sdr_data["class"])
        if _entry is not None and _dep is not None and _cls is not None:
            _branches = [[]]
            for _b in range(1, len(sdr_data["branches"])):
                _branches.append([node_map[n] for n in sdr_data["branches"][_b]])
            _entry.set_state_dep_routing(
                _cls, _dep, _branches,
                [int(v) for v in sdr_data["level"]],
                [float(v) for v in sdr_data["C"]],
                np.asarray(sdr_data["d"], dtype=float))

    # Restore routing weights (for WRROBIN)
    routing_weights = data.get("routingWeights", {})
    if routing_weights:
        from ..lang.base import RoutingStrategy as RS
        for node_name, class_weights in routing_weights.items():
            node = node_map.get(node_name)
            if node is None:
                continue
            for cls_name, dest_weights in class_weights.items():
                cls = class_map.get(cls_name)
                if cls is None:
                    continue
                for dest_name, weight in dest_weights.items():
                    dest = node_map.get(dest_name)
                    if dest is not None:
                        node.set_routing(cls, RS.WRROBIN, dest, weight)

    # Restore setup / delay-off, polling type and switchover times
    for nd_data in data.get("nodes", []):
        node = node_map.get(nd_data["name"])
        if node is None:
            continue

        # Setup / delay-off. The writer emits the two maps together.
        su_json = nd_data.get("setupTime")
        doff_json = nd_data.get("delayOffTime")
        if su_json and doff_json:
            for cname, su_data in su_json.items():
                jobclass = class_map.get(cname)
                if jobclass is None or cname not in doff_json:
                    continue
                su_dist = _json_to_dist(su_data)
                doff_dist = _json_to_dist(doff_json[cname])
                if su_dist is not None and doff_dist is not None:
                    node.set_delay_off(jobclass, su_dist, doff_dist)

        # Server breakdown/repair, with the optional per-class degraded
        # down-server service.
        bd_json = nd_data.get("breakdown")
        if bd_json:
            if "failure" not in bd_json or "repair" not in bd_json:
                raise ValueError('Node "%s": "breakdown" requires both a "failure" and a '
                                 '"repair" distribution.' % nd_data["name"])
            fdist = _json_to_dist(bd_json["failure"])
            rdist = _json_to_dist(bd_json["repair"])
            down_list = []
            down_json = bd_json.get("downService") or {}
            if down_json:
                # Indexed by class, in the order the classes were declared.
                order = [cd["name"] for cd in data.get("classes", [])]
                down_list = [_json_to_dist(down_json[cn]) if cn in down_json else None
                             for cn in order]
            node.set_breakdown(fdist, rdist, down_list if down_list else None)

        # Polling type, restored by name. This must precede the switchover
        # restore: set_polling_type resets every class to Immediate.
        pt_name = nd_data.get("pollingType")
        if pt_name:
            from ..constants import PollingType
            polling_type = PollingType.fromString(pt_name)
            if polling_type is not None:
                if polling_type == PollingType.KLIMITED:
                    node.set_polling_type(polling_type, int(nd_data.get("pollingPar", 1)))
                else:
                    node.set_polling_type(polling_type)

        # Switchover times: entries without "to" carry the per-class polling
        # form, entries with "to" the (from, to) pair form.
        so_list = nd_data.get("switchoverTimes")
        if not so_list:
            continue
        for so in so_list:
            from_cls = class_map.get(so.get("from"))
            dist_data = so.get("distribution")
            if from_cls is None or not dist_data:
                continue
            dist = _json_to_dist(dist_data)
            if dist is None:
                continue
            if so.get("to") is not None:
                to_cls = class_map.get(so.get("to"))
                if to_cls is None:
                    continue
                node.set_switchover(from_cls, to_cls, dist)
            else:
                node.set_switchover(from_cls, dist)

    # Restore finite capacity regions
    fcr_list = data.get("finiteCapacityRegions", [])
    for rj in fcr_list:
        # Support both new "stations" format and old "nodes" format
        region_nodes = []
        if "stations" in rj:
            for sj in rj["stations"]:
                n = node_map.get(sj.get("node"))
                if n is not None:
                    region_nodes.append(n)
        elif "nodes" in rj:
            for nname in rj["nodes"]:
                n = node_map.get(nname)
                if n is not None:
                    region_nodes.append(n)
        max_jobs = rj.get("globalMaxJobs", -1)
        if region_nodes:
            # region-load failures must surface, not be swallowed; add_region takes (name,*nodes), separate job-cap setter (positional max_jobs -> bogus node).
            region = model.addRegion(rj.get("name", "FCR"), *region_nodes)
            if max_jobs is not None and max_jobs > 0:
                region.set_global_max_jobs(int(max_jobs))
            # globalMaxMemory
            gmm = rj.get("globalMaxMemory")
            if gmm is not None:
                region.set_global_max_memory(gmm)
            # classMaxJobs
            for cls_name, cmj in rj.get("classMaxJobs", {}).items():
                cls = class_map.get(cls_name)
                if cls is not None:
                    region.set_class_max_jobs(cls, cmj)
            # classMaxMemory
            for cls_name, cmm in rj.get("classMaxMemory", {}).items():
                cls = class_map.get(cls_name)
                if cls is not None:
                    region.set_class_max_memory(cls, int(cmm))
            # dropRule
            for cls_name, dr_str in rj.get("dropRule", {}).items():
                cls = class_map.get(cls_name)
                if cls is not None:
                    region.set_drop_rule(cls, _str_to_drop_strategy(dr_str))
            # Per-station classWeight and classSize
            if "stations" in rj:
                for sj in rj["stations"]:
                    for cls_name, w in sj.get("classWeight", {}).items():
                        cls = class_map.get(cls_name)
                        if cls is not None:
                            region.set_class_weight(cls, w)
                    for cls_name, sz in sj.get("classSize", {}).items():
                        cls = class_map.get(cls_name)
                        if cls is not None:
                            region.set_class_size(cls, sz)
            # Linear admission constraints A*x <= b
            if rj.get("constraintA") is not None and rj.get("constraintB") is not None:
                region.set_linear_constraints(
                    np.atleast_2d(np.asarray(rj["constraintA"], dtype=float)),
                    np.asarray(rj["constraintB"], dtype=float))

    # per-node initial state restored for EVERY stateful node: a Place's token
    # counts and a closed PAS station's ordered job placement are two readings
    # of the same field, and a station whose state the document names is
    # initialized whether or not it is one of those two. Restoring only those
    # left the model reading as uninitialized, so the first solve rebuilt the
    # default state and the document's own was discarded.
    from ..lang.base import StatefulNode as _StatefulNode
    for nd_data in data.get("nodes", []):
        init_st = nd_data.get("initialState")
        if init_st is None:
            continue
        node = node_map.get(nd_data.get("name"))
        if node is None or not isinstance(node, _StatefulNode):
            continue
        if isinstance(node, Place):
            node.set_state(init_st)
            continue
        setter = getattr(node, 'set_state', None) or getattr(node, 'setState', None)
        if setter is not None:
            state_row = np.atleast_1d(np.asarray(init_st, dtype=float))
            setter(state_row)
            # AND ITS ONE-ROW STATE SPACE AND TRIVIAL PRIOR. The writer omits a
            # [1] prior over one row because `initialState` already carries that
            # row, so restoring the row alone leaves a node whose state, space
            # and prior disagree -- which is not what init_default or any
            # init_from_marginal* builds: all of them set the TRIO together, and
            # a solver that indexes the state space finds it empty. A space or
            # prior the document DID carry was installed earlier and is kept.
            space = node.get_state_space() if hasattr(node, 'get_state_space') else None
            if space is None or np.asarray(space, dtype=object).size == 0:
                if hasattr(node, 'set_state_space'):
                    node.set_state_space(np.atleast_2d(state_row))
                if hasattr(node, 'setStatePrior'):
                    node.setStatePrior(np.array([1.0]))

    _json_to_rewards(model, data, node_map, class_map)

    return model


def _create_node(model, nd: Dict[str, Any]):
    """Create a node from JSON data (before classes are set up)."""
    from ..lang.nodes import Queue, Delay, Source, Sink, Fork, Join, Router, ClassSwitch, Cache, Place, Transition
    from ..lang.base import SchedStrategy, ReplacementStrategy

    name = nd["name"]
    ntype = nd["type"]

    if ntype == "Source":
        return Source(model, name)
    elif ntype == "Sink":
        return Sink(model, name)
    elif ntype == "Delay":
        return Delay(model, name)
    elif ntype == "Queue":
        sched = _str_to_sched_strategy(nd.get("scheduling", "FCFS"))
        node = Queue(model, name, sched)
        servers = _servers_from_json(nd.get("servers", 1))
        if servers > 1:
            node.set_number_of_servers(servers)
        buf = nd.get("buffer")
        if buf is not None:
            node._capacity = buf
        return node
    elif ntype == "Fork":
        node = Fork(model, name)
        tpl = nd.get("tasksPerLink")
        if tpl is not None and tpl > 1:
            node.set_tasks_per_link(tpl)
        return node
    elif ntype == "Join":
        node = Join(model, name)
        # joinStrategy and joinQuorum are applied in post-class phase
        # (set_strategy/set_required require class objects)
        return node
    elif ntype == "Router":
        return Router(model, name)
    elif ntype == "Logger":
        # The wire carries the base file name only; the directory is the
        # model's log path, defaulted here so a round-trip does not fail on a
        # model that saved cleanly.
        from ..lang.nodes import Logger
        if not model.get_log_path():
            import os
            model.set_log_path(os.getcwd())
        return Logger(model, name, nd.get("fileName", "default.csv"))
    elif ntype == "ClassSwitch":
        return ClassSwitch(model, name)
    elif ntype == "Cache":
        # Support both Python format (nested "cache" key) and JAR format
        # (flat top-level with JAR key names)
        cc = nd.get("cache", {})
        # read per KEY, nested first: a nested block that carries only some of
        # the keys (MATLAB nests hitClass/missClass/popularity) must not shadow
        # the flat siblings, or the cache silently reloads at the 10-item,
        # capacity-1, LRU defaults.
        nitems = cc.get("items", nd.get("numItems", nd.get("items", 10)))
        cap = cc.get("capacity", nd.get("itemLevelCap", nd.get("capacity", 1)))
        repl_str = cc.get("replacement",
                          nd.get("replacementStrategy", nd.get("replacement", "LRU")))
        # admissionProb read from the flat key first (MATLAB/JAR emit flat), since nested-vs-flat if/else would miss a flat key alongside a nested block.
        qadm = cc.get("admissionProb", nd.get("admissionProb", 1.0))
        if isinstance(cap, list):
            cap = np.array(cap)
        repl = getattr(ReplacementStrategy, repl_str, ReplacementStrategy.LRU)
        cache = Cache(model, name, nitems, cap, repl)
        if qadm != 1.0:
            cache.set_admission_prob(float(qadm))
        isz = cc.get("itemSizes", nd.get("itemSizes"))
        if isz is not None:
            cache.set_item_sizes(np.asarray(isz, dtype=float).ravel())
        ccap = cc.get("costCaps", nd.get("costCaps"))
        if ccap is not None:
            cache.set_cost_caps(np.asarray(ccap, dtype=float).ravel())
        return cache
    elif ntype == "Place":
        # a queueing Place's embedded scheduling strategy must be restored, or it reloads at the INF default and silently becomes a delay.
        sched_str = nd.get("scheduling")
        if sched_str is None:
            return Place(model, name)
        sched = _str_to_sched_strategy(sched_str)
        node = Place(model, name, sched)
        servers = nd.get("servers")
        if servers is not None:
            node.set_number_of_servers(_servers_from_json(servers))
        return node
    elif ntype == "Transition":
        return Transition(model, name)
    else:
        # Default to Queue FCFS
        return Queue(model, name, SchedStrategy.FCFS)


def _str_to_sched_strategy(s: str):
    """Convert a JSON scheduling name to a SchedStrategy enum member.

    Enum names cross the wire uppercase. FCFSPRIO is the MATLAB/JAR alias for
    HOL (MATLAB gives both the same numeric id), so it is resolved to HOL.
    """
    from ..lang.base import SchedStrategy
    _sched_aliases = {"FCFSPRIO": "HOL"}
    key = str(s).upper()
    key = _sched_aliases.get(key, key)
    sched = getattr(SchedStrategy, key, None)
    if sched is None:
        from ..api.io.logging import line_warning
        line_warning('linemodel_io',
                     'Unrecognised scheduling strategy "%s"; defaulting to '
                     'FCFS.' % (s,))
        return SchedStrategy.FCFS
    return sched


def _servers_from_json(servers):
    """Decode a JSON server count, mapping the "Infinity" form to infinity."""
    if isinstance(servers, str):
        if servers.lower() == "infinity":
            return float('inf')
        return float(servers)
    return servers


def _drop_strategy_to_str(ds) -> Optional[str]:
    """Convert a DropStrategy enum to JSON string."""
    from ..lang.base import DropStrategy
    _map = {
        DropStrategy.DROP: "drop",
        DropStrategy.WaitingQueue: "waitingQueue",
        DropStrategy.BAS: "blockingAfterService",
        DropStrategy.RETRIAL: "retrial",
        DropStrategy.RETRIAL_WITH_LIMIT: "retrialWithLimit",
    }
    return _map.get(ds, None)


def _str_to_drop_strategy(s: str):
    """Convert a JSON drop rule string to DropStrategy enum.

    An unrecognised string resolves to WaitingQueue, matching MATLAB
    (linemodel_load.m str_to_droprule) and the JAR (parseDropStrategy).
    Python used to default to DROP, which turned an unreadable rule into
    silent job loss.
    """
    from ..lang.base import DropStrategy
    _map = {
        "drop": DropStrategy.DROP,
        "waitingQueue": DropStrategy.WaitingQueue,
        "blockingAfterService": DropStrategy.BAS,
        "retrial": DropStrategy.RETRIAL,
        "retrialWithLimit": DropStrategy.RETRIAL_WITH_LIMIT,
    }
    return _map.get(s, DropStrategy.WaitingQueue)


def _node_type_str(node) -> str:
    """Get the schema node type string for a node object."""
    from ..lang.nodes import (Queue, Delay, Source, Sink, Fork, Join, Router, ClassSwitch, Cache,
                              Place, Transition, Logger)
    if isinstance(node, Source):
        return "Source"
    if isinstance(node, Sink):
        return "Sink"
    if isinstance(node, Delay):
        return "Delay"
    if isinstance(node, Cache):
        return "Cache"
    if isinstance(node, Place):
        return "Place"
    if isinstance(node, Transition):
        return "Transition"
    if isinstance(node, Queue):
        return "Queue"
    if isinstance(node, Fork):
        return "Fork"
    if isinstance(node, Join):
        return "Join"
    if isinstance(node, Router):
        return "Router"
    if isinstance(node, ClassSwitch):
        return "ClassSwitch"
    if isinstance(node, Logger):
        return "Logger"
    # No silent default: a node type this writer does not know was saved as an
    # FCFS Queue, which reloads with a station and a service process the model
    # never declared and shifts every station index after it.
    raise TypeError('Node "%s" is of class %s, which has no model.json node type; add it to '
                    '_node_type_str rather than letting it default.'
                    % (node.get_name(), type(node).__name__))


# ---------------------------------------------------------------------------
# LayeredNetwork serialization
# ---------------------------------------------------------------------------

def _lincon_to_json(elem, col_names):
    """
    Admission constraint rows of a Task or Processor on the wire.

    Both declaration forms are normalised to the named form, so the wire is
    order-independent: a positional setConstraint(A, b) matrix is resolved
    against COL_NAMES (the element's entries, or its tasks) at write time.
    COL_NAMES must be in the same declaration order the positional columns
    assume.
    """
    if not hasattr(elem, 'hasLinearConstraints') or not elem.hasLinearConstraints():
        return []
    rows_json = []
    A, b = elem.getLinearConstraints()
    if A is not None and b is not None:
        for r in range(A.shape[0]):
            nz = np.flatnonzero(A[r, :])
            if nz.size == 0:
                continue
            if nz.max() >= len(col_names):
                raise ValueError(f"Admission constraint on {elem.name} references column "
                                 f"{int(nz.max()) + 1} but the element has only "
                                 f"{len(col_names)} operands.")
            rows_json.append({"operands": [col_names[j] for j in nz],
                              "coeffs": [float(A[r, j]) for j in nz],
                              "cap": float(b[r])})
    for names, coeffs, cap in elem.lincon_rows:
        rows_json.append({"operands": list(names),
                          "coeffs": [float(c) for c in coeffs],
                          "cap": float(cap)})
    return rows_json


def _apply_lincon(elem, rows):
    """
    Replay admission constraint rows from the wire onto a Task or Processor.
    Rows name their operands, so no column order is assumed and the referenced
    entries or tasks need not exist yet.
    """
    for row in rows or []:
        ops = row.get("operands")
        if isinstance(ops, str):
            ops = [ops]
        elem.addConstraint(list(ops), row.get("coeffs"), row.get("cap"))


def _layered_to_json(model) -> Dict[str, Any]:
    """Convert a LayeredNetwork to JSON-compatible dict."""
    from ..layered import (
        LayeredNetwork, Processor, Task, Entry, Activity,
        ActivityPrecedence, PrecedenceType, CallType, Distribution
    )
    from ..constants import SchedStrategy

    result = {
        "type": "LayeredNetwork",
        "name": model.name,
    }

    # LQN hosts schema always carries multiplicity/scheduling/quantum/speedFactor explicitly (never omitted), since JAR/native INF-vs-PS defaults differ.
    procs = []
    for p in model.processors:
        pj = {"name": p.name}
        pj["multiplicity"] = _mult_to_json(p.multiplicity)
        sched = p.sched_strategy
        sname = (sched.name if hasattr(sched, 'name') else str(sched)) if sched is not None else "PS"
        pj["scheduling"] = sname
        pj["quantum"] = p.getQuantum()
        pj["speedFactor"] = p.getSpeedFactor()
        repl = p.getReplication()
        if repl > 1:
            pj["replication"] = repl
        # Admission constraints: columns of a host constraint are its tasks
        rows_json = _lincon_to_json(p, [t.name for t in getattr(p, 'tasks', [])])
        if rows_json:
            pj["admissionConstraints"] = rows_json
        procs.append(pj)
    result["hosts"] = procs

    # Tasks. Canonical schema keys the parent processor as "host" and always
    # carries multiplicity/scheduling.
    tasks = []
    for t in model.tasks:
        tj = {"name": t.name}
        # Host (parent processor)
        if t.processor is not None:
            tj["host"] = t.processor.name
        # Multiplicity (infinite-server tasks serialize as INF_MULTIPLICITY)
        tj["multiplicity"] = _mult_to_json(t.multiplicity)
        # Scheduling
        sched = t.sched_strategy
        sname = (sched.name if hasattr(sched, 'name') else str(sched)) if sched is not None else "FCFS"
        tj["scheduling"] = sname
        # think time emitted both as the canonical scalar (thinkTimeMean/SCV) and the richer distribution object (native round-trip).
        if t.think_time is not None:
            tt_mean = _get_dist_mean_safe(t.think_time)
            if tt_mean > 0:
                tj["thinkTime"] = _layered_dist_to_json(t.think_time)
                tj["thinkTimeMean"] = tt_mean
                try:
                    tj["thinkTimeSCV"] = float(t.think_time.getSCV())
                except Exception:
                    tj["thinkTimeSCV"] = 1.0
        # Fan in/out
        if t._fan_in:
            tj["fanIn"] = dict(t._fan_in)
        if t._fan_out:
            tj["fanOut"] = dict(t._fan_out)
        # Replication
        repl = t.getReplication()
        if repl > 1:
            tj["replication"] = repl
        # SetupTask detection
        from ..layered import SetupTask as _SetupTask
        if isinstance(t, _SetupTask) or getattr(t, '_is_setup_task', False):
            tj["taskType"] = "SetupTask"
        # Setup time / delay-off time (on any Task)
        if t.setup_time is not None:
            st_mean = _get_dist_mean_safe(t.setup_time)
            if st_mean > 0:
                tj["setupTime"] = _layered_dist_to_json(t.setup_time)
        if t.delay_off_time is not None:
            dot_mean = _get_dist_mean_safe(t.delay_off_time)
            if dot_mean > 0:
                tj["delayOffTime"] = _layered_dist_to_json(t.delay_off_time)
        # CacheTask detection
        from ..layered import CacheTask as _CacheTask
        if isinstance(t, _CacheTask):
            tj["taskType"] = "CacheTask"
            tj["totalItems"] = t.total_items
            tj["cacheCapacity"] = t.cache_capacity
            rs = t.replacement_strategy
            tj["replacementStrategy"] = rs.name if hasattr(rs, 'name') else str(rs)
        # Admission constraints: columns of a task constraint are its entries
        rows_json = _lincon_to_json(t, [e.name for e in getattr(t, 'entries', [])])
        if rows_json:
            tj["admissionConstraints"] = rows_json
        tasks.append(tj)
    result["tasks"] = tasks

    # Entries
    entries = []
    for e in model.entries:
        ej = {"name": e.name}
        if e.task is not None:
            ej["task"] = e.task.name
        # Entry arrival distribution
        arrival = getattr(e, '_arrival', None)
        if arrival is not None:
            adj = _dist_to_json(arrival)
            if adj is not None:
                ej["arrival"] = adj
        # ItemEntry detection
        from ..layered import ItemEntry as _ItemEntry
        if isinstance(e, _ItemEntry):
            ej["entryType"] = "ItemEntry"
            ej["totalItems"] = e.total_items
            if e.access_prob is not None:
                # Save access probability as a list if it's array-like
                ap = e.access_prob
                if hasattr(ap, 'tolist'):
                    ej["accessProb"] = ap.tolist()
                elif isinstance(ap, (list, tuple)):
                    ej["accessProb"] = list(ap)
                else:
                    # Try to serialize as distribution
                    dj = _dist_to_json(ap)
                    if dj is not None:
                        ej["accessProb"] = dj
        entries.append(ej)
    result["entries"] = entries

    # Activities
    activities = []
    for a in model.activities:
        aj = {"name": a.name}
        if a.task is not None:
            aj["task"] = a.task.name
        # Host demand
        if a.host_demand is not None:
            aj["hostDemand"] = _layered_dist_to_json(a.host_demand)
        # boundToEntry is the canonical LQN JSON key; the native reader also accepts the legacy boundTo.
        if a.bound_entry is not None:
            aj["boundToEntry"] = a.bound_entry.name
        # Replies to entry
        if a.reply_entry is not None:
            aj["repliesTo"] = a.reply_entry.name
        # calls key the destination entry as 'dest' (canonical) with mean call count; the native reader also accepts the legacy 'entry' key.
        synch_calls = []
        asynch_calls = []
        for entry, mean_calls, call_type in a.calls:
            call_obj = {"dest": entry.name, "mean": mean_calls}
            if call_type == CallType.SYNC:
                synch_calls.append(call_obj)
            elif call_type == CallType.ASYNC:
                asynch_calls.append(call_obj)
        if synch_calls:
            aj["synchCalls"] = synch_calls
        if asynch_calls:
            aj["asynchCalls"] = asynch_calls
        activities.append(aj)
    result["activities"] = activities

    # Precedences
    precedences = []
    for t in model.tasks:
        for p in t.precedences:
            pj = {"task": t.name}
            if p.prec_type == PrecedenceType.SERIAL:
                pj["type"] = "Serial"
                pj["activities"] = [a.name for a in p.activities]
            elif p.prec_type == PrecedenceType.PARALLEL:
                if len(p.pre_activities) == 1 and len(p.post_activities) > 1:
                    pj["type"] = "AndFork"
                    pj["activities"] = [p.pre_activities[0].name] + [a.name for a in p.post_activities]
                elif len(p.pre_activities) > 1 and len(p.post_activities) == 1:
                    pj["type"] = "AndJoin"
                    pj["activities"] = [a.name for a in p.pre_activities] + [p.post_activities[0].name]
                else:
                    pj["type"] = "Serial"
                    pj["activities"] = [a.name for a in (p.pre_activities + p.post_activities)]
            elif p.prec_type == PrecedenceType.CHOICE:
                if len(p.pre_activities) == 1 and len(p.post_activities) > 1:
                    pj["type"] = "OrFork"
                    pj["activities"] = [p.pre_activities[0].name] + [a.name for a in p.post_activities]
                    if p.probabilities:
                        pj["probabilities"] = list(p.probabilities)
                elif len(p.pre_activities) > 1 and len(p.post_activities) == 1:
                    pj["type"] = "OrJoin"
                    pj["activities"] = [a.name for a in p.pre_activities] + [p.post_activities[0].name]
                else:
                    pj["type"] = "Serial"
                    pj["activities"] = [a.name for a in (p.pre_activities + p.post_activities)]
            elif p.prec_type == PrecedenceType.LOOP:
                pj["type"] = "Loop"
                pj["activities"] = [a.name for a in p.activities]
                if p.pre_activities:
                    pj["preActivity"] = p.pre_activities[0].name
                if p.count != 1.0:
                    pj["loopCount"] = p.count
            elif p.prec_type == PrecedenceType.CACHE_ACCESS:
                pj["type"] = "CacheAccess"
                pj["activities"] = [p.pre_activities[0].name] + [a.name for a in p.post_activities]
            else:
                continue
            precedences.append(pj)
    if precedences:
        result["precedences"] = precedences

    return result


def _json_to_layered(data: Dict[str, Any]):
    """Reconstruct a LayeredNetwork from JSON data."""
    from ..layered import (
        LayeredNetwork, Processor, Task, Entry, Activity,
        ActivityPrecedence, PrecedenceType, CallType, Distribution as LDist,
        SetupTask, CacheTask, ItemEntry
    )
    from ..constants import SchedStrategy
    from ..lang.base import ReplacementStrategy
    from ..distributions.continuous import Immediate

    model = LayeredNetwork(data.get("name", "model"))

    # Create processors (Python schema: "processors", JAR schema: "hosts")
    proc_map = {}
    for pd in data.get("processors", data.get("hosts", [])):
        name = pd["name"]
        mult = _mult_from_json(pd.get("multiplicity"))
        sched_str = pd.get("scheduling", "INF")
        sched = getattr(SchedStrategy, sched_str, SchedStrategy.INF)
        proc = Processor(model, name, mult, sched)
        q = pd.get("quantum", 0.0)
        if q > 0:
            proc.setQuantum(q)
        sf = pd.get("speedFactor", 1.0)
        if sf != 1.0:
            proc.setSpeedFactor(sf)
        repl = pd.get("replication")
        if repl is not None and repl > 1:
            proc.setReplication(repl)
        _apply_lincon(proc, pd.get("admissionConstraints"))
        proc_map[name] = proc

    # Create tasks
    task_map = {}
    for td in data.get("tasks", []):
        name = td["name"]
        mult = _mult_from_json(td.get("multiplicity"))
        sched_str = td.get("scheduling", "INF")
        sched = getattr(SchedStrategy, sched_str, SchedStrategy.INF)
        proc_name = td.get("processor", td.get("host"))
        proc = proc_map.get(proc_name) if proc_name else None
        task_type = td.get("taskType", "Task")
        if task_type in ("SetupTask", "FunctionTask"):
            # FunctionTask is the legacy name of SetupTask on the wire
            task = SetupTask(model, name, mult, sched)
        elif task_type == "CacheTask":
            total_items = td.get("totalItems", 1)
            cache_cap = td.get("cacheCapacity", 1)
            rs_name = td.get("replacementStrategy", "FIFO")
            rs = getattr(ReplacementStrategy, rs_name, ReplacementStrategy.FIFO)
            task = CacheTask(model, name, total_items, cache_cap, rs, mult)
        else:
            task = Task(model, name, mult, sched)
        if proc is not None:
            task.on(proc)
        # Think time (Python schema: "thinkTime" as dist, JAR schema: "thinkTimeMean"/"thinkTimeSCV")
        tt = td.get("thinkTime")
        if tt is not None:
            task.set_think_time(_json_to_native_or_layered_dist(tt))
        else:
            tt_mean = td.get("thinkTimeMean", 0.0)
            if tt_mean > 0:
                from ..distributions.continuous import Exp
                task.set_think_time(Exp(1.0 / tt_mean))
        # Setup time
        st = td.get("setupTime")
        if st is not None:
            task.set_setup_time(_json_to_native_or_layered_dist(st))
        else:
            st_mean = td.get("setupTimeMean", 0.0)
            if st_mean > 1e-8:
                from ..distributions.continuous import Exp
                task.set_setup_time(Exp(1.0 / st_mean))
        # Delay-off time
        dot = td.get("delayOffTime")
        if dot is not None:
            task.set_delay_off_time(_json_to_native_or_layered_dist(dot))
        else:
            dot_mean = td.get("delayOffTimeMean", 0.0)
            if dot_mean > 1e-8:
                from ..distributions.continuous import Exp
                task.set_delay_off_time(Exp(1.0 / dot_mean))
        # Fan in/out
        fi = td.get("fanIn", {})
        for src, val in fi.items():
            task.setFanIn(src, val)
        fo = td.get("fanOut", {})
        for dst, val in fo.items():
            task.setFanOut(dst, val)
        # Replication
        repl = td.get("replication")
        if repl is not None and repl > 1:
            task.setReplication(repl)
        _apply_lincon(task, td.get("admissionConstraints"))
        task_map[name] = task

    # Create entries
    entry_map = {}
    for ed in data.get("entries", []):
        name = ed["name"]
        task_name = ed.get("task")
        task = task_map.get(task_name) if task_name else None
        entry_type = ed.get("entryType", "Entry")
        if entry_type == "ItemEntry":
            total_items = ed.get("totalItems", 1)
            access_prob = ed.get("accessProb")
            # access_prob may be a list or a distribution dict
            if isinstance(access_prob, dict):
                access_prob = _json_to_dist(access_prob)
            entry = ItemEntry(model, name, total_items, access_prob)
        else:
            entry = Entry(model, name)
        if task is not None:
            entry.on(task)
        # Entry arrival distribution
        arv = ed.get("arrival")
        if arv is not None:
            arv_dist = _json_to_dist(arv)
            if arv_dist is not None:
                entry.setArrival(arv_dist)
        entry_map[name] = entry

    # Create activities
    act_map = {}
    for ad in data.get("activities", []):
        name = ad["name"]
        task_name = ad.get("task")
        task = task_map.get(task_name) if task_name else None
        hd = ad.get("hostDemand")
        host_demand = _json_to_native_or_layered_dist(hd) if hd else Immediate.getInstance()
        act = Activity(model, name, host_demand)
        if task is not None:
            act.on(task)
        # Bound to (Python schema: "boundTo", JAR schema: "boundToEntry")
        bt = ad.get("boundTo", ad.get("boundToEntry"))
        if bt and bt in entry_map:
            act.bound_to(entry_map[bt])
        # Replies to
        rt = ad.get("repliesTo")
        if rt and rt in entry_map:
            act.replies_to(entry_map[rt])
        # Calls (Python schema: "entry", JAR schema: "dest")
        for sc in ad.get("synchCalls", []):
            ename = sc.get("entry", sc.get("dest"))
            mean = sc.get("mean", 1.0)
            if ename and ename in entry_map:
                act.synch_call(entry_map[ename], mean)
        for ac in ad.get("asynchCalls", []):
            ename = ac.get("entry", ac.get("dest"))
            mean = ac.get("mean", 1.0)
            if ename and ename in entry_map:
                act.asynch_call(entry_map[ename], mean)
        act_map[name] = act

    # Create precedences (Python schema: "type"/"activities", JAR schema: "preActs"/"postActs"/"preType"/"postType")
    for pd in data.get("precedences", []):
        task_name = pd.get("task")
        task = task_map.get(task_name)
        if task is None:
            continue

        # Detect schema: JAR uses preActs/postActs, Python uses type/activities
        if "preActs" in pd or "postActs" in pd:
            # JAR schema
            pre_names = pd.get("preActs", [])
            post_names = pd.get("postActs", [])
            pre_type = pd.get("preType", "pre")
            post_type = pd.get("postType", "post")
            pre_acts = [act_map[n] for n in pre_names if n in act_map]
            post_acts = [act_map[n] for n in post_names if n in act_map]
            probs = pd.get("probabilities", [])

            if pre_type == "pre" and post_type == "post":
                # Simple serial: pre→post
                if len(pre_acts) == 1 and len(post_acts) == 1:
                    task.add_precedence(ActivityPrecedence.Serial(pre_acts + post_acts))
            elif pre_type == "pre" and post_type in ("and-fork", "post-AND"):
                if len(pre_acts) >= 1 and len(post_acts) >= 1:
                    task.add_precedence(ActivityPrecedence.AndFork(pre_acts[0], post_acts))
            elif pre_type in ("and-join", "pre-AND") and post_type == "post":
                if len(pre_acts) >= 1 and len(post_acts) >= 1:
                    task.add_precedence(ActivityPrecedence.AndJoin(pre_acts, post_acts[0]))
            elif pre_type == "pre" and post_type in ("or-fork", "post-OR"):
                if len(pre_acts) >= 1 and len(post_acts) >= 1:
                    post_params = pd.get("postParams", probs)
                    if not post_params:
                        n = len(post_acts)
                        post_params = [1.0 / n] * n
                    task.add_precedence(ActivityPrecedence.OrFork(pre_acts[0], post_acts, post_params))
            elif pre_type in ("or-join", "pre-OR") and post_type == "post":
                if len(pre_acts) >= 1 and len(post_acts) >= 1:
                    task.add_precedence(ActivityPrecedence.OrJoin(pre_acts, post_acts[0]))
            elif pre_type == "pre" and post_type in ("loop", "post-LOOP"):
                count = pd.get("loopCount", None)
                if count is None:
                    # JAR format uses postParams for loop count
                    post_params = pd.get("postParams", [])
                    count = post_params[0] if post_params else 1.0
                if len(pre_acts) >= 1 and len(post_acts) >= 1:
                    task.add_precedence(ActivityPrecedence.Loop(pre_acts[0], post_acts, count))
            elif pre_type == "pre" and post_type == "post-CACHE":
                if len(pre_acts) >= 1 and len(post_acts) >= 1:
                    task.add_precedence(ActivityPrecedence.CacheAccess(pre_acts[0], post_acts))
        else:
            # Python schema
            ptype = pd.get("type", "Serial")
            act_names = pd.get("activities", [])
            acts = [act_map[n] for n in act_names if n in act_map]

            if ptype == "Serial" and len(acts) >= 2:
                task.add_precedence(ActivityPrecedence.Serial(acts))
            elif ptype == "AndFork" and len(acts) >= 2:
                task.add_precedence(ActivityPrecedence.AndFork(acts[0], acts[1:]))
            elif ptype == "AndJoin" and len(acts) >= 2:
                task.add_precedence(ActivityPrecedence.AndJoin(acts[:-1], acts[-1]))
            elif ptype == "OrFork" and len(acts) >= 2:
                probs = pd.get("probabilities", [])
                if not probs:
                    n = len(acts) - 1
                    probs = [1.0 / n] * n
                task.add_precedence(ActivityPrecedence.OrFork(acts[0], acts[1:], probs))
            elif ptype == "OrJoin" and len(acts) >= 2:
                task.add_precedence(ActivityPrecedence.OrJoin(acts[:-1], acts[-1]))
            elif ptype == "Loop":
                count = pd.get("loopCount", 1.0)
                pre_name = pd.get("preActivity")
                if pre_name and pre_name in act_map:
                    pre_act = act_map[pre_name]
                    task.add_precedence(ActivityPrecedence.Loop(pre_act, acts, count))
                elif len(acts) >= 2:
                    # Legacy format: first activity is pre, rest is loop body
                    task.add_precedence(ActivityPrecedence.Loop(acts[0], acts[1:], count))
            elif ptype == "CacheAccess" and len(acts) >= 2:
                task.add_precedence(ActivityPrecedence.CacheAccess(acts[0], acts[1:]))

    return model


def _layered_dist_to_json(dist) -> Dict[str, Any]:
    """Convert a layered network distribution (dataclass or native) to JSON."""
    from ..layered import Distribution as LDist
    if isinstance(dist, LDist):
        # Dataclass distribution — use Exp with fitMean
        if dist.scv == 0:
            return {"type": "Det", "params": {"value": dist.mean}}
        return {"type": "Exp", "params": {"lambda": 1.0 / dist.mean if dist.mean > 0 else 1e8}}
    # Try native distribution
    return _dist_to_json(dist) or {"type": "Exp", "params": {"lambda": 1.0}}


def _json_to_layered_dist(d: Dict[str, Any]):
    """Convert JSON dist to a layered Distribution (dataclass)."""
    from ..layered import Distribution as LDist
    # Use native distribution conversion, then get mean
    native = _json_to_dist(d)
    if native is not None:
        try:
            mean = native.getMean()
            scv = native.getSCV()
            return LDist(mean=mean, scv=scv)
        except Exception:
            return LDist(mean=1.0, scv=1.0)
    return LDist(mean=1.0, scv=1.0)


def _json_to_native_or_layered_dist(d: Dict[str, Any]):
    """Convert JSON dist to a native distribution if possible, else a layered Distribution dataclass.
    Native distributions have getMean()/getSCV() methods needed by solvers."""
    native = _json_to_dist(d)
    if native is not None:
        return native
    # Fallback to layered Distribution dataclass
    from ..layered import Distribution as LDist
    return LDist(mean=1.0, scv=1.0)


def _get_dist_mean_safe(dist) -> float:
    """Get mean from any distribution safely."""
    if dist is None:
        return 0.0
    if isinstance(dist, (int, float)):
        return float(dist)
    if hasattr(dist, 'mean') and not callable(dist.mean):
        return dist.mean
    if hasattr(dist, 'getMean'):
        return dist.getMean()
    return 0.0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_GD_MAX_LATTICE = 200000


def _gd_slots(sn, wcut):
    """
    Varying (station, class) coordinates of the global-dependence lattice, with
    their cutoffs. A Source holds no jobs and a class with zero per-class capacity
    at a station never appears there, so those entries are pinned to 0 and carry no
    coordinate. That restriction is lossless: no DEP or PHASE event ever fires at a
    slot the class cannot occupy, so phi is never read there.
    """
    from ..lang.base import NodeType
    M, K = int(sn.nstations), int(sn.nclasses)
    slot_st, slot_cl, cuts = [], [], []
    for i in range(M):
        if sn.nodetype[int(sn.stationToNode[i])] == NodeType.SOURCE:
            continue
        for r in range(K):
            cap = float(sn.classcap[i, r])
            if not cap > 0:
                continue
            nj = float(sn.njobs[r])
            c = int(round(nj)) if math.isfinite(nj) else int(wcut)
            if math.isfinite(cap):
                c = min(c, int(round(cap)))
            slot_st.append(i)
            slot_cl.append(r)
            cuts.append(max(c, 0))
    return slot_st, slot_cl, cuts


def _gd_block(model) -> Dict[str, Any]:
    """Materialize the network-level global (Whittle) dependence onto the wire."""
    sn = model.getStruct()
    M, K = int(sn.nstations), int(sn.nclasses)
    phi = model.get_global_dependence()
    peak = np.asarray(model.get_global_dependence_peak(), dtype=float).reshape(M, K)
    wcut = int(model.get_global_dependence_cutoff())
    slot_st, slot_cl, cuts = _gd_slots(sn, wcut)

    total = 1
    for c in cuts:
        total *= (c + 1)
    if total > _GD_MAX_LATTICE:
        raise ValueError(
            "The global dependence lattice has %d points (%d varying station-class "
            "slots with cutoffs %s), above the wire limit of %d. Lower the cutoff "
            "argument of set_global_dependence, or solve the model natively."
            % (total, len(cuts), cuts, _GD_MAX_LATTICE))

    stations = [model.getNodes()[int(sn.stationToNode[i])].getName() for i in range(M)]
    classes = [jc.getName() for jc in model.getClasses()]
    return {
        "type": "globalDependent",
        "stations": stations,
        "classes": classes,
        "slots": [{"station": stations[slot_st[s]], "class": classes[slot_cl[s]]}
                  for s in range(len(cuts))],
        "cutoffs": [int(c) for c in cuts],
        "cutoff": wcut,
        "scaling": _gd_scaling_table(phi, slot_st, slot_cl, cuts, M, K),
        "peak": [float(x) for x in peak.reshape(-1)],
    }


def _gd_scaling_table(phi, slot_st, slot_cl, cuts, M, K) -> Dict[str, list]:
    """
    Tabulate phi over the slot box lattice. Key: comma-joined 0-based slot counts
    in slot order. Value: the FULL (M, K) scaling flattened row-major, so the
    reader restores the matrix without re-deriving which return form was used.
    """
    tbl = {}
    P = len(cuts)
    shp = [c + 1 for c in cuts]
    total = 1
    for x in shp:
        total *= x
    for i in range(total):
        li, c = i, [0] * P
        for d in range(P):
            c[d] = li % shp[d]
            li //= shp[d]
        n = np.zeros((M, K))
        for d in range(P):
            n[slot_st[d], slot_cl[d]] = c[d]
        v = np.asarray(phi(n), dtype=float)
        if v.size == 1:
            v = np.full((M, K), float(v.reshape(-1)[0]))
        elif v.shape == (M,) or v.shape == (M, 1):
            v = np.repeat(v.reshape(M, 1), K, axis=1)
        v = np.where(np.isfinite(v), v, 0.0)
        key = ",".join(str(x) for x in c) if P else "0"
        tbl[key] = [float(x) for x in v.reshape(-1)]
    return tbl


def _gd_block_to_callable(gdep, node_map, class_map, model):
    """
    Rebuild the global (Whittle) dependence from the slot lattice. Slots carry
    station and class names, resolved through the model's own index spaces. The
    population is clamped to the tabulated cutoffs, which is the same saturation
    the writer's box lattice declares.
    """
    sn = model.getStruct()
    M, K = int(sn.nstations), int(sn.nclasses)
    slot_st, slot_cl = [], []
    for sm in gdep.get("slots", []):
        node = node_map.get(sm["station"])
        cls = class_map.get(sm["class"])
        if node is None or cls is None:
            return None, None, 10
        slot_st.append(int(sn.nodeToStation[node._index]))
        slot_cl.append(int(cls._index))
    cuts = [int(c) for c in gdep.get("cutoffs", [])]
    wcut = int(gdep.get("cutoff", 10) or 10)

    tbl = {}
    for key, vals in gdep["scaling"].items():
        tbl[key] = np.asarray(vals, dtype=float).reshape(M, K)
    peak = gdep.get("peak")
    peak = np.asarray(peak, dtype=float).reshape(M, K) if peak else np.ones((M, K))
    ones = np.ones((M, K))
    P = len(cuts)

    def _phi(n):
        if P == 0:
            return tbl.get("0", ones)
        c = [min(max(int(round(n[slot_st[s], slot_cl[s]])), 0), cuts[s]) for s in range(P)]
        return tbl.get(",".join(str(x) for x in c), ones)

    return _phi, peak, wcut


def _class_dependence_cutoffs(classes) -> list:
    """
    Per-class cutoffs for materializing a class-dependence callable onto a
    bounded lattice. A closed class cannot exceed its population; an open class
    is unbounded, so it gets the same saturation cutoff the OI/PAS rate table
    uses, beyond which beta is taken to be constant. Mirrors the rule in the
    MATLAB writer (linemodel_save) and the JAR (LineModelIO).
    """
    from ..lang.classes import ClosedClass
    maxc = []
    for jc in classes:
        if isinstance(jc, ClosedClass) and math.isfinite(jc.getPopulation()):
            maxc.append(int(round(jc.getPopulation())))
        else:
            maxc.append(10)
    return maxc


def _cd_scaling_table(beta, maxc, K) -> Dict[str, list]:
    """
    Materialize the class-dependence callable beta(n) over the box lattice
    0 <= n[r] <= maxc[r]. Keyed by the comma-joined 0-based per-class counts,
    matching the JAR reader. beta may return a scalar (one scaling shared by
    every class) or a length-K vector; the scalar form is broadcast to K entries
    so the reader is uniform.
    """
    tbl = {}
    total = 1
    for c in maxc:
        total *= (c + 1)
    for li in range(total):
        rem = li
        n = [0] * K
        for r in range(K):
            n[r] = rem % (maxc[r] + 1)
            rem //= (maxc[r] + 1)
        v = beta(np.array(n, dtype=float))
        v = np.atleast_1d(np.asarray(v, dtype=float)).flatten()
        if v.size == 1:
            v = np.repeat(v[0], K)
        v = np.where(np.isfinite(v), v, 0.0)
        tbl[",".join(str(x) for x in n)] = [float(x) for x in v[:K]]
    return tbl


def _firingdep_table_to_handle(frm, node_map, class_map):
    """
    Rebuild a marking-dependent firing-rate dependence callable g(M) from the
    materialized lattice written by the transition-mode serializer. M is the
    node-indexed marking matrix; the enabling (place,class) slots are read off,
    clamped to the cutoffs, and the tabulated scalar multiplier is looked up
    (default 1 outside the tabulated range). Mirrors firingdep_table_to_handle
    in linemodel_load.m and the JAR reader.
    """
    slots = frm.get("slots")
    scaling = frm.get("scaling")
    if not slots or scaling is None:
        return None
    slot_idx = []
    for sm in slots:
        nd = node_map.get(sm["node"])
        jc = class_map.get(sm["class"])
        if nd is None or jc is None:
            return None
        slot_idx.append((nd.get_index0(), jc.get_index0()))
    cutoffs = [int(c) for c in frm.get("cutoffs", [])]
    table = {k: float(v) for k, v in scaling.items()}

    def _g(M):
        Mm = np.atleast_2d(np.asarray(M, dtype=float))
        key = []
        for s, (ni, ci) in enumerate(slot_idx):
            c = int(round(Mm[ni, ci]))
            if c < 0:
                c = 0
            if s < len(cutoffs) and c > cutoffs[s]:
                c = cutoffs[s]
            key.append(c)
        return table.get(",".join(str(x) for x in key), 1.0)

    return _g


def _cd_table_to_callable(tbl: Dict[str, list], cutoffs, K):
    """
    Rebuild a class-dependence callable from a materialized lattice table. The
    population is clamped to the cutoffs, so beta saturates beyond the tabulated
    range exactly as the table intends. A state absent from the table is neutral
    (beta=1), leaving the nominal service rate unscaled.
    """
    cut = [int(c) for c in cutoffs] if cutoffs else None

    def _beta(n):
        nn = np.atleast_1d(np.asarray(n, dtype=float)).flatten()
        key = []
        for r in range(K):
            v = int(round(nn[r])) if r < nn.size else 0
            if v < 0:
                v = 0
            if cut is not None and r < len(cut) and v > cut[r]:
                v = cut[r]
            key.append(v)
        row = tbl.get(",".join(str(x) for x in key))
        if row is None:
            return np.ones(K)
        out = np.ones(K)
        for r in range(min(K, len(row))):
            out[r] = row[r]
        return out

    return _beta


def _to_list(arr) -> list:
    """Convert numpy array or list to plain Python list."""
    if isinstance(arr, np.ndarray):
        return arr.tolist()
    if isinstance(arr, (list, tuple)):
        return list(arr)
    return [arr]


def _matrix_to_list(mat) -> list:
    """Convert 2D numpy array to nested list."""
    if isinstance(mat, np.ndarray):
        return mat.tolist()
    return mat


class _JSONEncoder(json.JSONEncoder):
    """Custom JSON encoder that handles numpy types."""
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.bool_):
            return bool(obj)
        return super().default(obj)


def _wire_nonfinite(doc):
    """Rewrite every non-finite number into the wire's own spelling.

    JSON HAS NO INFINITY LITERAL. `json.dump` writes a bare `Infinity` /
    `NaN` anyway, which is a Python extension no strict parser takes:
    `linemodel_save` (the reference) writes an infinite scalar as the STRING
    "Infinity" / "-Infinity" and a NaN as `null`, for every numeric field
    rather than a named few, and the C++ reader `num_from_json` decodes
    exactly that pair. So a model whose Source carries its declared state --
    an open class holds an infinite population, so `initialState` and
    `stateSpace` carry the sentinel the moment the model is initialized --
    left here as a document nlohmann refuses at PARSE time. That is how
    fcr_oqndrop[CROSS:Matlab->Python] failed: the native `common/ldes` engine
    rejected the file, SolverLDES fell back to `ldes.jar` without saying so,
    and the two sides of the row were then compared across two engines and
    differed by 1.3e-2 with nothing in the output naming the cause.
    """
    if isinstance(doc, dict):
        return dict((k, _wire_nonfinite(v)) for k, v in doc.items())
    if isinstance(doc, (list, tuple)):
        return [_wire_nonfinite(v) for v in doc]
    if isinstance(doc, np.ndarray):
        return _wire_nonfinite(doc.tolist())
    if isinstance(doc, (float, np.floating)) and not isinstance(doc, bool):
        v = float(doc)
        if math.isnan(v):
            return None
        if math.isinf(v):
            return "Infinity" if v > 0 else "-Infinity"
    return doc


# ---------------------------------------------------------------------------
# Workflow serialization
# ---------------------------------------------------------------------------

def _workflow_to_json(model) -> Dict[str, Any]:
    """Convert a Workflow model to a JSON-compatible dict."""
    from ..distributions.continuous import Exp

    result = {
        "type": "Workflow",
        "name": model.name,
    }

    # Activities
    acts_json = []
    for act in model.getActivities():
        act_obj = {"name": act.name}
        if act._distribution is not None:
            act_obj["hostDemand"] = _dist_to_json(act._distribution)
        else:
            # Activity was created from a float mean — serialize as Exp
            act_obj["hostDemand"] = _dist_to_json(Exp(1.0 / act._host_demand_mean))
        acts_json.append(act_obj)
    result["activities"] = acts_json

    # Precedences
    precs_json = []
    for prec in model.getPrecedences():
        prec_obj = {
            "preActs": list(prec.pre_acts),
            "postActs": list(prec.post_acts),
            "preType": prec.pre_type,
            "postType": prec.post_type,
        }
        if prec.pre_params is not None and len(prec.pre_params) > 0:
            prec_obj["preParams"] = _to_list(prec.pre_params)
        if prec.post_params is not None and len(prec.post_params) > 0:
            prec_obj["postParams"] = _to_list(prec.post_params)
        precs_json.append(prec_obj)
    result["precedences"] = precs_json

    return result


def _json_to_workflow(data: Dict[str, Any]):
    """Convert a JSON dict back to a Workflow model."""
    from ..lang.workflow import Workflow, ActivityPrecedence

    name = data.get("name", "workflow")
    wf = Workflow(name)

    # Activities
    if "activities" in data:
        for act_data in data["activities"]:
            act_name = act_data["name"]
            if "hostDemand" in act_data:
                host_demand = _json_to_dist(act_data["hostDemand"])
                wf.addActivity(act_name, host_demand)
            else:
                wf.addActivity(act_name, 1.0)

    # Precedences
    if "precedences" in data:
        for prec_data in data["precedences"]:
            pre_acts = list(prec_data["preActs"])
            post_acts = list(prec_data["postActs"])
            pre_type = prec_data.get("preType", "pre")
            post_type = prec_data.get("postType", "post")
            pre_params = None
            if "preParams" in prec_data:
                pre_params = np.atleast_1d(np.array(prec_data["preParams"], dtype=float))
            post_params = None
            if "postParams" in prec_data:
                post_params = np.atleast_1d(np.array(prec_data["postParams"], dtype=float))
            wf.addPrecedence(ActivityPrecedence(
                pre_acts=pre_acts,
                post_acts=post_acts,
                pre_type=pre_type,
                post_type=post_type,
                pre_params=pre_params,
                post_params=post_params
            ))

    return wf


# ---------------------------------------------------------------------------
# Environment serialization
# ---------------------------------------------------------------------------

def _environment_to_json(model) -> Dict[str, Any]:
    """Convert an Environment model to a JSON-compatible dict."""
    result = {
        "type": "Environment",
        "name": model.name,
        "numStages": model.num_stages,
    }

    # Stages
    stages_json = []
    for i in range(model.num_stages):
        stage_name = model._stage_names[i] if i < len(model._stage_names) else None
        if stage_name is None:
            continue
        stage_obj = {"name": stage_name}
        # The stage TYPE, which every reader already looks for and no writer
        # emitted: an Environment round-tripped through JSON came back with its
        # stage types blanked, so getStageTable and any consumer keying off
        # UP/DOWN read a different environment from the one that was saved.
        stage_type = model._stage_types[i] if i < len(model._stage_types) else ''
        if stage_type:
            stage_obj["type"] = str(stage_type)
        stage_model = model.get_model(i)
        if stage_model is not None:
            stage_obj["model"] = _network_to_json(stage_model)
        stages_json.append(stage_obj)
    result["stages"] = stages_json

    # Transitions
    trans_json = []
    for i in range(model.num_stages):
        for j in range(model.num_stages):
            if model.env[i][j] is not None:
                trans_obj = {
                    "from": i,
                    "to": j,
                    "distribution": _dist_to_json(model.env[i][j]),
                }
                trans_json.append(trans_obj)
    result["transitions"] = trans_json

    # node-failure declarative record also carries queue-length reset policies (callables), not recoverable from stage/transition structure alone.
    node_failures_json = _node_failures_to_json(model)
    if node_failures_json:
        result["nodeFailures"] = node_failures_json

    return result


def _node_failures_to_json(model) -> list:
    """Serialize the environment's node breakdown/repair descriptors."""
    from ..api.io.logging import line_warning

    node_failures = getattr(model, '_node_failures', None)
    if not node_failures:
        return []

    nf_json = []
    for nf in node_failures:
        node_name = nf['node']
        breakdown = _dist_to_json(nf['breakdown'])
        if breakdown is None:
            line_warning('linemodel_io',
                         'Node failure on "%s" has a breakdown distribution that cannot be serialized; '
                         'the nodeFailures entry is omitted.' % node_name)
            continue
        down_service = _dist_to_json(nf['downService'])
        if down_service is None:
            line_warning('linemodel_io',
                         'Node failure on "%s" has a down-service distribution that cannot be serialized; '
                         'the nodeFailures entry is omitted.' % node_name)
            continue
        nj = {"node": node_name, "breakdownRate": breakdown}
        if nf['repair'] is not None:
            repair = _dist_to_json(nf['repair'])
            if repair is None:
                line_warning('linemodel_io',
                             'Node failure on "%s" has a repair distribution that cannot be serialized; '
                             'the nodeFailures entry is omitted.' % node_name)
                continue
            nj["repairRate"] = repair
        nj["downService"] = down_service
        if nf['breakdownResetPolicy'] == 'custom':
            line_warning('linemodel_io',
                         'Node failure on "%s" uses a custom breakdown reset function, which cannot be '
                         'serialized to JSON; the saved model falls back to the "keep" policy on reload.'
                         % node_name)
        else:
            nj["breakdownResetPolicy"] = nf['breakdownResetPolicy']
        if nf['repairResetPolicy']:
            if nf['repairResetPolicy'] == 'custom':
                line_warning('linemodel_io',
                             'Node failure on "%s" uses a custom repair reset function, which cannot be '
                             'serialized to JSON; the saved model falls back to the "keep" policy on reload.'
                             % node_name)
            else:
                nj["repairResetPolicy"] = nf['repairResetPolicy']
        nf_json.append(nj)
    return nf_json


def _node_failure_fields(nf: Dict[str, Any]):
    """Decode the distributions and reset policies of one "nodeFailures" entry.

    Note that breakdownRate/repairRate carry full distributions, not scalar rates.
    """
    node_name = nf.get("node")
    if not nf.get("breakdownRate"):
        raise ValueError('Node failure on "%s" is missing the required "breakdownRate" field.' % node_name)
    if not nf.get("downService"):
        raise ValueError('Node failure on "%s" is missing the required "downService" field.' % node_name)
    breakdown_dist = _json_to_dist(nf["breakdownRate"])
    down_service_dist = _json_to_dist(nf["downService"])
    repair_dist = _json_to_dist(nf["repairRate"]) if nf.get("repairRate") else None
    reset_b = nf.get("breakdownResetPolicy") or 'keep'
    reset_r = nf.get("repairResetPolicy") or 'keep'
    return breakdown_dist, repair_dist, down_service_dist, reset_b, reset_r


def _json_to_environment(data: Dict[str, Any]):
    """Convert a JSON dict back to an Environment model."""
    from ..environment import Environment

    name = data.get("name", "env")
    num_stages = data.get("numStages", 0)
    env = Environment(name, num_stages)

    # nodeFailures: expand model to UP/DOWN stages if absent, else restore callable reset policies; _kb/12-interfaces-and-docs.md linemodel_io.py section.
    nf_arr = data.get("nodeFailures", [])
    stages = data.get("stages", [])
    declared_names = [sd.get("name", "Stage" + str(i)) for i, sd in enumerate(stages)]

    macro_mode = bool(nf_arr)
    for nf in nf_arr:
        if "node" not in nf:
            raise ValueError('A "nodeFailures" entry is missing the required "node" field.')
        if ('DOWN_%s' % nf["node"]) in declared_names:
            macro_mode = False
            break

    if macro_mode:
        if len(stages) != 1:
            raise ValueError('"nodeFailures" expands the base model into the UP and DOWN_<node> stages, so '
                             '"stages" must declare exactly one stage, holding the base (UP) model.')
        if data.get("transitions"):
            raise ValueError('"nodeFailures" implies the breakdown and repair transitions; "transitions" '
                             'must not be declared alongside it.')
        if not stages[0].get("model"):
            raise ValueError('"nodeFailures" requires the base stage to carry a "model".')
        base_model = _json_to_network(stages[0]["model"])
        for nf in nf_arr:
            breakdown, repair, down_service, reset_b, reset_r = _node_failure_fields(nf)
            if repair is None:
                env.add_node_breakdown(base_model, nf["node"], breakdown, down_service, reset_b)
            else:
                env.add_node_failure_repair(base_model, nf["node"], breakdown, repair, down_service,
                                            reset_b, reset_r)
    else:
        # Stages
        for i, stage_data in enumerate(stages):
            stage_name = declared_names[i]
            stage_type = stage_data.get("type", "")
            stage_model = None
            if "model" in stage_data:
                stage_model = _json_to_network(stage_data["model"])
            if stage_model is not None:
                env.add_stage(i, stage_name, stage_type, stage_model)

        # Transitions
        if "transitions" in data:
            for trans_data in data["transitions"]:
                from_idx = trans_data["from"]
                to_idx = trans_data["to"]
                dist = _json_to_dist(trans_data["distribution"])
                env.add_transition(from_idx, to_idx, dist)

        # Re-attach the node-failure descriptors and their reset policies to the
        # stages just built, so that the environment serializes back identically.
        for nf in nf_arr:
            breakdown, repair, down_service, reset_b, reset_r = _node_failure_fields(nf)
            env.register_node_failure(nf["node"], breakdown, repair, down_service, reset_b, reset_r)

    env.init()
    return env


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def save_model(model, filename: str) -> None:
    """
    Save a LINE model to a JSON file.

    Supports Network, LayeredNetwork, Workflow, and Environment models.
    The output conforms to the line-model.schema.json specification.

    Args:
        model: A Network, LayeredNetwork, Workflow, or Environment instance.
        filename: Output file path (should end in .json).

    Example:
        >>> model = Network('M/M/1')
        >>> # ... define model ...
        >>> save_model(model, 'mymodel.json')
    """
    from ..lang.network import Network
    from ..layered import LayeredNetwork
    from ..lang.workflow import Workflow
    from ..environment import Environment

    if isinstance(model, LayeredNetwork):
        model_json = _layered_to_json(model)
    elif isinstance(model, Network):
        model_json = _network_to_json(model)
    elif isinstance(model, Workflow):
        model_json = _workflow_to_json(model)
    elif isinstance(model, Environment):
        model_json = _environment_to_json(model)
    else:
        raise TypeError(f"Unsupported model type: {type(model)}")

    doc = {
        "format": "line-model",
        "version": "1.0",
        "model": model_json
    }

    with open(filename, 'w') as f:
        # allow_nan=False so a non-finite that escapes the rewrite is reported
        # here rather than written as a literal no other reader accepts.
        json.dump(_wire_nonfinite(doc), f, indent=2, cls=_JSONEncoder,
                  allow_nan=False)


def load_model(filename: str):
    """
    Load a LINE model from a JSON file.

    Returns a Network, LayeredNetwork, Workflow, or Environment
    depending on the model type in the JSON file.

    Args:
        filename: Path to a .json file conforming to line-model.schema.json.

    Returns:
        Network, LayeredNetwork, Workflow, or Environment instance.

    Example:
        >>> model = load_model('mymodel.json')
        >>> solver = SolverMVA(model)
        >>> print(solver.avg_table())
    """
    with open(filename, 'r') as f:
        doc = json.load(f)

    model_data = doc.get("model", {})
    mtype = model_data.get("type")

    if mtype == "Network":
        return _json_to_network(model_data)
    elif mtype == "LayeredNetwork":
        return _json_to_layered(model_data)
    elif mtype == "Workflow":
        return _json_to_workflow(model_data)
    elif mtype == "Environment":
        return _json_to_environment(model_data)
    else:
        raise ValueError(f"Unsupported model type: {mtype}")
