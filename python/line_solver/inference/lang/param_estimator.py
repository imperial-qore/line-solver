import numpy as np
from scipy.interpolate import interp1d

from line_solver import MetricType, Network, Delay, Queue, Source, Sink, Exp
from line_solver import ClosedClass, OpenClass, SchedStrategy, SolverMVA

from line_solver.inference.lang.sampled_metric import SampledMetric


class ParamEstimator:
    """Service demand estimator for queueing network models.

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """

    def __init__(self, model, options=None):
        if options is None:
            options = ParamEstimator.default_options()
        self.model = model
        self.options = options
        n_nodes = model.getNumberOfNodes()
        n_classes = model.getNumberOfClasses()
        # samples[i][r] is a list of SampledMetric objects
        self.samples = [[[] for _ in range(n_classes)] for _ in range(n_nodes)]
        # samplesAggr[i] is a list of aggregate SampledMetric objects
        self.samples_aggr = [[] for _ in range(n_nodes)]

    def add_samples(self, sample_data):
        i = self.model.get_node_index(sample_data.node) - 1  # 0-based
        if sample_data.is_aggregate():
            self.samples_aggr[i].append(sample_data)
        else:
            r = self.model.get_class_index(sample_data.jobclass) - 1  # 0-based
            self.samples[i][r].append(sample_data)

    def get_data(self):
        return self.samples

    def get_data_aggr(self):
        return self.samples_aggr

    def get_arvr(self, node, jobclass):
        i = self.model.get_node_index(node) - 1
        r = self.model.get_class_index(jobclass) - 1
        for sm in self.samples[i][r]:
            if sm.type == MetricType.ArvR:
                return sm
        return None

    def get_util(self, node, jobclass):
        i = self.model.get_node_index(node) - 1
        r = self.model.get_class_index(jobclass) - 1
        for sm in self.samples[i][r]:
            if sm.type == MetricType.Util:
                return sm
        return None

    def get_respt(self, node, jobclass):
        i = self.model.get_node_index(node) - 1
        r = self.model.get_class_index(jobclass) - 1
        for sm in self.samples[i][r]:
            if sm.type == MetricType.RespT:
                return sm
        return None

    def get_aggr_util(self, node):
        i = self.model.get_node_index(node) - 1
        for sm in self.samples_aggr[i]:
            if sm.type == MetricType.Util:
                return sm
        return None

    def get_qlen(self, node, jobclass, ev=None):
        i = self.model.get_node_index(node) - 1
        r = self.model.get_class_index(jobclass) - 1
        node_data = self.samples[i][r]
        if ev is None:
            result = [sm for sm in node_data if sm.type == MetricType.QLen]
            if len(result) == 1:
                return result[0]
            return result if result else None
        else:
            for sm in node_data:
                if (sm.type == MetricType.QLen and sm.cond is not None
                        and sm.cond.node == ev.node
                        and sm.cond.jobclass == ev.jobclass
                        and sm.cond.event == ev.event):
                    return sm
            return None

    def get_aggr_qlen(self, node, ev=None):
        i = self.model.get_node_index(node) - 1
        node_data = self.samples_aggr[i]
        if ev is None:
            for sm in node_data:
                if sm.type == MetricType.QLen:
                    return sm
            return None
        else:
            for sm in node_data:
                if (sm.type == MetricType.QLen and sm.cond is not None
                        and sm.cond.node == ev.node
                        and getattr(sm.cond, 'jobclass', getattr(sm.cond, 'class', None)) == getattr(ev, 'jobclass', getattr(ev, 'class', None))
                        and sm.cond.event == ev.event):
                    return sm
            return None

    def get_tput(self, node, jobclass):
        i = self.model.get_node_index(node) - 1
        r = self.model.get_class_index(jobclass) - 1
        for sm in self.samples[i][r]:
            if sm.type == MetricType.Tput:
                return sm
        return None

    def auto_method(self):
        """Automatically select the best estimation method based on available data."""
        has_arvr = False
        has_respt = False
        has_util = False
        has_qlen = False
        has_tput = False
        has_trace = False
        has_aggr_util = False
        has_aggr_qlen = False

        for i in range(len(self.samples)):
            for r in range(len(self.samples[i])):
                for sm in self.samples[i][r]:
                    if sm.type == MetricType.ArvR:
                        has_arvr = True
                    if sm.type == MetricType.RespT:
                        has_respt = True
                    if sm.type == MetricType.Util:
                        has_util = True
                    if sm.type == MetricType.QLen:
                        has_qlen = True
                    if sm.type == MetricType.Tput:
                        has_tput = True
                    if sm.is_trace():
                        has_trace = True

        for i in range(len(self.samples_aggr)):
            for sm in self.samples_aggr[i]:
                if sm.type == MetricType.Util:
                    has_aggr_util = True
                if sm.type == MetricType.QLen:
                    has_aggr_qlen = True

        if has_trace and has_respt and has_aggr_qlen:
            method = 'erps'
        elif has_trace and has_respt and has_arvr:
            method = 'mlps'
        elif has_arvr and has_respt and (has_util or has_aggr_util):
            method = 'ubo'
        elif has_arvr and (has_util or has_aggr_util):
            method = 'ubr'
        elif has_qlen:
            method = 'qmle'
        else:
            raise ValueError(
                'Insufficient data to automatically select an estimation method. '
                'Please set options["method"] manually.')

        self.options['method'] = method
        return method

    def interpolate(self):
        """Interpolate all data across all available timestamps."""
        tunion = np.array([], dtype=float)

        for i in range(len(self.samples)):
            for r in range(len(self.samples[i])):
                for sm in self.samples[i][r]:
                    if sm is not None:
                        tunion = np.union1d(tunion, sm.t)

        for i in range(len(self.samples_aggr)):
            for sm in self.samples_aggr[i]:
                if sm is not None:
                    tunion = np.union1d(tunion, sm.t)

        if len(tunion) == 0:
            return

        for i in range(len(self.samples)):
            for r in range(len(self.samples[i])):
                for sm in self.samples[i][r]:
                    if sm is not None and len(sm.t) > 1:
                        f = interp1d(sm.t, sm.data, kind='cubic',
                                     fill_value='extrapolate')
                        sm.data = f(tunion)
                        sm.t = tunion.copy()

        for i in range(len(self.samples_aggr)):
            for sm in self.samples_aggr[i]:
                if sm is not None and len(sm.t) > 1:
                    f = interp1d(sm.t, sm.data, kind='cubic',
                                 fill_value='extrapolate')
                    sm.data = f(tunion)
                    sm.t = tunion.copy()

    def estimate_at(self, nodes):
        """Dispatch to appropriate estimator and update model parameters."""
        if not isinstance(nodes, (list, tuple)):
            nodes = [nodes]

        method = self.options['method']
        if method == 'ubr':
            est_val = self._estimator_ubr(nodes)
        elif method == 'ubo':
            est_val = self._estimator_ubo(nodes)
        elif method == 'erps':
            est_val = self._estimator_erps(nodes)
        elif method == 'ekf':
            est_val = self._estimator_ekf(nodes)
        elif method == 'mcmc':
            est_val = self._estimator_mcmc(nodes)
        elif method == 'mle':
            est_val = self._estimator_mle(nodes)
        elif method == 'mlps':
            est_val = self._estimator_mlps(nodes)
        elif method == 'fmlps':
            est_val = self._estimator_fmlps(nodes)
        elif method == 'qmle':
            est_val = self._estimator_qmle(nodes)
        elif method == 'gibbs':
            est_val = self._estimator_gibbs(nodes)
        elif method == 'vi':
            est_val = self._estimator_variational(nodes)
        else:
            raise ValueError(f'Unknown inference method: {method}.')

        # Update model parameters
        sn = (self.model.refreshStruct(), self.model.getStruct())[1]
        classes = self.model.getClasses()
        est_val = self._as_estimate_matrix(est_val, len(nodes), int(sn.nclasses))
        for n_idx in range(len(nodes)):
            nd = nodes[n_idx]
            if isinstance(nd, (Source, Sink)):
                continue
            for r in range(sn.nclasses):
                val = est_val[n_idx, r] if est_val.ndim == 2 else est_val[r]
                if not np.isnan(val) and val > 0:
                    nd.setService(classes[r], Exp.fitMean(val))
        self.model.reset()

        return est_val

    def _as_estimate_matrix(self, est_val, n_nodes, n_classes):
        """Put an estimator's output into the one shape callers may rely on.

        The estimators disagree on what they hand back -- `ubo` returns a
        (nodes x classes) matrix, `ekf` and `mle` a per-class vector, and some
        report a flat vector over EVERY node of the model -- and every caller
        then has to guess. MATLAB and the JAR both settle this by returning a
        (nodes x classes) matrix from `estimateAt`, so that is the shape imposed
        here, with one concession to Python: a single requested node yields the
        per-class vector directly, which is what `estVal[r]` in the examples
        reads and what MATLAB's own `estVal(r)` linear indexing amounts to.
        """
        arr = np.asarray(est_val, dtype=float).ravel()
        if arr.size == n_nodes * n_classes:
            rows = n_nodes
        elif n_classes > 0 and arr.size % n_classes == 0:
            # an estimator that reports every node of the model, not only the
            # requested ones; keep all of its rows
            rows = arr.size // n_classes
        else:
            raise ValueError(
                'estimator returned %d values, which is not a whole number of '
                'per-class rows (%d classes)' % (arr.size, n_classes))
        matrix = arr.reshape(rows, n_classes)
        return matrix[0] if rows == 1 else matrix

    def _build_closed_equivalent_for_ps(self, node):
        """Build a closed equivalent model for open/mixed networks."""
        sn = (self.model.refreshStruct(), self.model.getStruct())[1]
        R = sn.nclasses
        classes = self.model.getClasses()

        Nopen = self.options.get('openPopulation', 100)

        N = np.zeros(R)
        Z = np.zeros(R)

        all_nodes = self.model.getNodes()
        for nd in all_nodes:
            if isinstance(nd, Delay):
                for r in range(R):
                    if sn.njobs[r] < np.inf:
                        svc = nd.getService(classes[r])
                        Z[r] += svc.getMean()
            elif isinstance(nd, Source):
                for r in range(R):
                    if sn.njobs[r] == np.inf:
                        svc = nd.getService(classes[r])
                        lambda_r = 1.0 / svc.getMean()
                        Z[r] = Nopen / lambda_r

        for r in range(R):
            if sn.njobs[r] < np.inf:
                N[r] = sn.njobs[r]
            else:
                N[r] = Nopen

        delay_rate = 1.0 / Z

        eq_model = Network('closed_equiv')
        eq_delay = Delay(eq_model, 'Think')
        eq_queue = Queue(eq_model, 'Queue1', SchedStrategy.PS)
        eq_queue.setNumberOfServers(node.getNumberOfServers())

        eq_classes = []
        for r in range(R):
            eq_class = ClosedClass(eq_model, f'Class{r + 1}', int(N[r]), eq_delay, 0)
            eq_classes.append(eq_class)
            eq_delay.setService(eq_class, Exp(delay_rate[r]))
            queue_svc = node.getService(classes[r])
            eq_queue.setService(eq_class, queue_svc)

        P = eq_model.initRoutingMatrix()
        for r in range(R):
            P.set(eq_classes[r], eq_classes[r], eq_delay, eq_queue, 1.0)
            P.set(eq_classes[r], eq_classes[r], eq_queue, eq_delay, 1.0)
        eq_model.link(P)

        return eq_model, eq_queue

    # ---- Estimator methods ----

    def _estimator_ubr(self, nodes):
        from line_solver.inference.api._estimators import estimator_ubr
        return estimator_ubr(self, nodes)

    def _estimator_ubo(self, nodes):
        from line_solver.inference.api._estimators import estimator_ubo
        return estimator_ubo(self, nodes)

    def _estimator_erps(self, nodes):
        from line_solver.inference.api._estimators import estimator_erps
        return estimator_erps(self, nodes)

    def _estimator_ekf(self, nodes):
        from line_solver.inference.api._estimators import estimator_ekf
        return estimator_ekf(self, nodes)

    def _estimator_mcmc(self, nodes):
        from line_solver.inference.api._estimators import estimator_mcmc
        return estimator_mcmc(self, nodes)

    def _estimator_mle(self, nodes):
        from line_solver.inference.api._estimators import estimator_mle
        return estimator_mle(self, nodes)

    def _estimator_mlps(self, nodes):
        from line_solver.inference.api._estimators import estimator_mlps
        return estimator_mlps(self, nodes)

    def _estimator_fmlps(self, nodes):
        from line_solver.inference.api._estimators import estimator_fmlps
        return estimator_fmlps(self, nodes)

    def _estimator_qmle(self, nodes):
        from line_solver.inference.api._estimators import estimator_qmle
        return estimator_qmle(self, nodes)

    def _estimator_gibbs(self, nodes):
        from line_solver.inference.api._estimators import estimator_gibbs
        return estimator_gibbs(self, nodes)

    def _estimator_variational(self, nodes):
        from line_solver.inference.api._estimators import estimator_variational
        return estimator_variational(self, nodes)

    @staticmethod
    def default_options():
        return {
            'verbose': 1,
            'method': 'ubr',
            'variant': 'default',
            'iter_max': 1000,
            'tol': 1e-3,
            'solver': SolverMVA,
            'openPopulation': 100,
        }

    @staticmethod
    def get_required_metrics(method):
        descs = {
            'ubr': 'ArvR (per-class) + Util (per-class or aggregate)',
            'ubo': 'ArvR (per-class) + RespT (per-class) + Util (aggregate)',
            'erps': 'RespT (per-class) + QLen (aggregate, conditional on class arrivals). PS stations only.',
            'ekf': 'RespT (per-class) + Util (aggregate). Sequential/recursive estimation.',
            'mcmc': 'QLen (aggregate). Gibbs sampling with MCMC. Open/mixed via closed equivalence.',
            'mle': 'ArvR (per-class) + RespT (per-class) + Util (aggregate)',
            'mlps': 'ArvR (per-class, trace) + RespT (per-class, trace). PS stations only.',
            'fmlps': 'ArvR (per-class, trace) + RespT (per-class, trace). PS stations only.',
            'qmle': 'QLen (per-class). Open/mixed via closed equivalence.',
            'gibbs': 'ArvR (per-class, trace) + RespT (per-class, trace) + Tput (per-class). Gibbs sampling.',
            'vi': 'QLen (per-class, timeseries) at every station. Variational inference over transition counts; noisy readings, Gamma posteriors.',
        }
        return descs.get(method, f'Unknown method: {method}')

    # MATLAB-compatible aliases
    addSamples = add_samples
    getData = get_data
    getDataAggr = get_data_aggr
    getArvR = get_arvr
    getUtil = get_util
    getRespT = get_respt
    getAggrUtil = get_aggr_util
    getQLen = get_qlen
    getAggrQLen = get_aggr_qlen
    getTput = get_tput
    autoMethod = auto_method
    estimateAt = estimate_at
    defaultOptions = default_options
    getRequiredMetrics = get_required_metrics
    buildClosedEquivalentForPS = _build_closed_equivalent_for_ps
