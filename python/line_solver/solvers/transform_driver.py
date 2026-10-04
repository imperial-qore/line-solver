"""Solver-agnostic driver of a model TRANSFORMATION, the sibling of
:mod:`fork_join_driver`.

Where the fork-join driver drives ONE transformation (MMT/HT), this drives
whichever one ``options.config['transform']`` names: the strategy rewrites the
model into subproblems, each subproblem is solved by a REAL solver, and the
strategy maps the metrics back onto the original classes and stations::

    expand -> [ for e in subproblems: solve(e); couple(e) ] -> converged -> lift

Mirrors MATLAB ``@NetworkSolver/transformSolve.m`` and
``matlab/src/solvers/TR/``.

THE INNER SOLVER IS THE OUTER SOLVER. Each subproblem is solved by
``type(self)(submodel, ...)``, so a transformation written once serves every
solver rather than the one it was first written for. ``_run_chain_aggregation``
used to hard-wire ``SolverCTMC(chain_model, ...)``, which is what kept a
transform with nothing CTMC-specific in it out of reach of MVA, NC and FLD.

SINGLE PASS BY DEFAULT: unless ``expand`` sets ``ctx['iterated']`` the loop runs
one sweep and lifts. Coupling is Gauss-Seidel by construction, since ``couple``
runs inside the sweep rather than after it.
"""

import time

#: Sweep cap when the caller's options carry none; the pfqn_bklc default.
_DEFAULT_ITER_MAX = 1000


class TransformResult:
    """What a transformed solve answers with, in ORIGINAL coordinates."""

    def __init__(self, Q, U, R, T, C, X, lG, runtime, iters, method):
        self.Q = Q
        self.U = U
        self.R = R
        self.T = T
        self.C = C
        self.X = X
        self.lG = lG
        self.runtime = runtime
        self.iter = iters
        self.method = method


def transform_method(method):
    """Normalise a transform method name onto its canonical name and strategy.

    An unrecognised method name is an error rather than a silent 'none': a mistyped
    transform that quietly solved the untransformed model would report a
    plausible number for the wrong problem.
    """
    if method is None or method == '':
        return 'none', None
    method_name = str(method).lower()
    if method_name in ('none', 'off'):
        return 'none', None
    if method_name in ('chains', 'chain', 'chainaggr', 'chain_aggregation'):
        return 'chains', _chains_strategy
    if method_name in ('lc', 'loadconceal', 'load_concealment', 'thinning'):
        return 'lc', _lc_strategy
    raise ValueError(
        "options.config['transform']=%r is not a known model transformation. "
        "Use 'none', 'chains' or 'lc'." % (method,))


def _chains_strategy(phase, solver, *args):
    """Chain aggregation: collapse every chain onto a single class.

    ``ModelAdapter.aggregate_chains`` builds the collapsed model and
    ``sn_deaggregate_chain_results`` maps its metrics back through alpha, the
    per-station share of the chain's visits each class carries. EXACT on a
    product-form model, an approximation otherwise, because one aggregate
    service law, fitted to the alpha-weighted first two moments, replaces the
    per-class ones.

    Single pass: one solve of the aggregate determines the answer, so
    ``ctx['iterated']`` stays False and no couple/converged phase exists.
    """
    if phase == 'expand':
        sn, = args
        from ..api.io.model_adapter import ModelAdapter
        # With one class per chain the aggregation is the IDENTITY, and
        # aggregate_chains then returns no deaggregation tables at all.
        # Refusing by name beats failing later on a None field.
        if int(sn.nchains) >= int(sn.nclasses):
            raise ValueError(
                'chain aggregation needs more classes than chains: this model has %d classes '
                'in %d chains, so the transform is the identity.'
                % (int(sn.nclasses), int(sn.nchains)))
        chain_model, alpha, deagg = ModelAdapter.aggregate_chains(solver.model)
        ctx = {'sn_orig': sn, 'alpha': alpha, 'deagg': deagg, 'iterated': False}
        return [chain_model], ctx
    if phase == 'lift':
        ctx, res = args
        from ..api.sn.deaggregate import sn_deaggregate_chain_results
        r = res[0]
        deagg = ctx['deagg']
        # ST is left None so the per-class service times are recovered from
        # sn.rates. X is the per-chain SYSTEM throughput, which is why the
        # driver collects it on its own channel rather than from getAvg.
        return sn_deaggregate_chain_results(
            ctx['sn_orig'], deagg.L_chain, None, deagg.ST_chain, deagg.V_chain,
            ctx['alpha'], r['Q'], r['U'], r['R'], r['T'], None, r['X'])
    raise ValueError('Unknown chain-transformation phase: %s' % phase)


class TransformSolveMixin:
    """Run a model transformation with this solver as its own inner solver."""

    def _transform_inner_config(self, depth):
        """Options the subproblem is solved under.

        The inner solve must not re-enter the driver, so the method name is cleared
        and the depth advanced. A KERNEL selection inside the inner solver is
        not a transform and is correctly not cut.
        """
        cfg = dict(self.options.config or {})
        cfg['transform'] = 'none'
        cfg['transform_depth'] = depth + 1
        return cfg

    def _transform_requested(self):
        """The canonical method name asked for, or None when no transform is wanted."""
        method_name = (self.options.config or {}).get('transform')
        if not method_name or str(method_name).lower() == 'none':
            return None
        return transform_method(method_name)[0]

    def _transform_publish(self, tr, method):
        """Store a transformed solve's answer in THIS solver's result container.

        The default is the dict shape SolverMVA and SolverNC use. SolverCTMC
        overrides it because it keeps a dataclass-like `_QRFResult` instead;
        python has no single result container, which is the same obstacle
        ForkJoinDriverMixin meets with `_fj_publish`.
        """
        import numpy as np

        from ..api.sn.getters import sn_get_arvr_from_tput

        sn = self._sn if getattr(self, '_sn', None) is not None else self._get_network_struct()
        AN = np.atleast_2d(np.asarray(sn_get_arvr_from_tput(sn, tr.T, None), dtype=float))
        self._result = {
            'QN': tr.Q, 'UN': tr.U, 'RN': tr.R, 'TN': tr.T, 'AN': AN,
            'XN': tr.X, 'WN': tr.R, 'CN': tr.C, 'lG': tr.lG,
            'runtime': tr.runtime, 'lastiter': tr.iter, 'iter': tr.iter,
            'method': method,
        }
        return self._result

    def maybe_transform(self):
        """Run a requested transformation, publish it, and report that we did.

        This is the seam every solver calls at the top of its runAnalyzer, the
        counterpart of MATLAB's runAnalyzerPreamble branch. Returns False when
        no transform was asked for, so the caller proceeds normally.
        """
        method_name = self._transform_requested()
        if method_name is None:
            return False
        sn = self._sn if getattr(self, '_sn', None) is not None else self._get_network_struct()
        tr = self.transform_solve(sn)
        self._transform_publish(tr, '%s/%s' % (self.options.method, tr.method))
        return True

    def transform_solve(self, sn):
        t0 = time.time()
        cfg = self.options.config or {}
        method_name, strategy = transform_method(cfg.get('transform'))
        if method_name == 'none':
            raise ValueError("transform_solve called with transform='none'.")
        depth = int(cfg.get('transform_depth') or 0)
        if depth > 0:
            raise ValueError(
                "a model transformation (%r) cannot be nested inside another one; "
                "options.config['transform_depth'] is %d." % (method_name, depth))

        submodels, ctx = strategy('expand', self, sn)
        iterated = bool(ctx.get('iterated', False))
        inner_cfg = self._transform_inner_config(depth)

        res = [None] * len(submodels)
        iters = 0
        # Not every solver's options carry iter_max (SolverCTMCOptions does not),
        # and defaulting an ITERATED strategy to one sweep would silently return
        # an unconverged answer that looks converged. The fallback is the same
        # 1000 the pfqn_bklc kernel uses.
        iter_max = int(getattr(self.options, 'iter_max', None) or _DEFAULT_ITER_MAX)
        for it in range(1, max(1, iter_max) + 1):
            iters = it
            for e, submodel in enumerate(submodels):
                # A FRESH solver per sweep, deliberately: getAvg caches on the
                # solver, so an iterated strategy re-solving a mutated model
                # through one handle would read the previous sweep and converge
                # to the wrong point.
                inner = type(self)(submodel, config=inner_cfg,
                                   method=self.options.method,
                                   verbose=self.options.verbose)
                res[e] = _collect_inner(inner)
                if iterated:
                    submodels, ctx = strategy('couple', self, ctx, submodels, res, e)
            if not iterated:
                break
            done, ctx = strategy('converged', self, ctx, res, it)
            if done:
                break

        out = strategy('lift', self, ctx, res)
        return TransformResult(out.Q, out.U, out.R, out.T, out.C, out.X,
                               res[0]['lG'], time.time() - t0, iters, method_name)


def _collect_inner(inner):
    """The inner solve answers on FOUR channels, not one.

    ``getAvg`` returns (Q, U, R, T, A, W) whose sixth entry is the RESIDENCE
    time, so the system throughput and the normalizing constant are separate
    asks.
    """
    import numpy as np

    avg = inner.getAvg()
    r = {'Q': avg[0], 'U': avg[1], 'R': avg[2], 'T': avg[3]}
    r['X'] = np.atleast_1d(np.asarray(inner.getAvgSysTput(), dtype=float)).ravel()
    result = getattr(inner, '_result', None)
    r['lG'] = getattr(result, 'lG', None)
    r['method'] = getattr(result, 'method', '')
    return r


class _LcIdentityResult:
    """The chain tables re-indexed onto classes, when the two coincide."""

    def __init__(self, Q, U, R, T, C, X):
        self.Q = Q
        self.U = U
        self.R = R
        self.T = T
        self.C = C
        self.X = X


def _lc_seed(L, N, Z):
    """Step 1 of Algorithm 2, reproducing pfqn_bklc's own seed exactly.

    The saddle point utilizations of Birman-Kogan Corollary 1, with the same
    fallbacks and the same capacity clamp, so the transformation starts the
    iteration from the same point as the kernel.
    """
    import numpy as np

    from ..api.pfqn.bk import pfqn_bk

    R = L.shape[1]
    try:
        X = np.asarray(pfqn_bk(L, N, Z)[2], dtype=float).ravel()
    except Exception:
        X = np.zeros(R)
    if X.size != R:
        X = np.zeros(R)
    X = np.where(np.isfinite(X) & (X >= 0), X, 0.0)
    for r in range(R):
        if X[r] == 0 and N[r] > 0:
            denom = Z[r] + float(L[:, r].sum())
            X[r] = N[r] / denom if denom > 0 else 0.0
    # A chain cannot draw more than the capacity of its own slowest station.
    for r in range(R):
        cap = float(L[:, r].max()) if L.shape[0] else 0.0
        if cap > 0:
            X[r] = min(X[r], 1.0 / cap)
    return X


def _lc_conceal_all(ctx, submodels):
    """Re-conceal every subproblem against the current throughput vector."""
    import numpy as np

    from ..distributions import dist_scale_rate

    L, X = ctx['L'], ctx['X']
    for l, m in enumerate(submodels):
        A = 1.0 - (L.dot(X) - L[:, l] * X[l])
        A = np.maximum(A, ctx['fine_tol'])
        jobclass = m._classes[0]
        for i in range(ctx['M']):
            if ctx['is_delay'][i]:
                continue
            # dist_scale_rate multiplies the RATE, so a factor of A_i divides
            # the mean by A_i: exactly the concealed demand L(i,l)/A_i.
            m._stations[i].setService(jobclass,
                                     dist_scale_rate(ctx['base_service'][l][i], float(A[i])))
        m.refreshStruct()
    return submodels


def _lc_chain_tput(r):
    import numpy as np

    x = np.atleast_1d(np.asarray(r['X'], dtype=float)).ravel()
    if x.size == 0 or not np.isfinite(x[0]) or x[0] < 0:
        return 0.0
    return float(x[0])


def _lc_strategy(phase, solver, *args):
    """Load concealment (Birman-Kogan Algorithm 2) as a transformation.

    Chain aggregation, then a Gauss-Seidel sweep in which chain l is solved on
    its own against the residual capacity the others leave it,
    A_i = 1 - sum_{k!=l} L(i,k) X_k, so it sees the concealed demand
    L(i,l)/A_i.

    The per-chain subproblem is a real single-class Network, so the inner solve
    is the CALLER'S own solver. It is NOT a more accurate `lc` than the
    `pfqn_bklc` kernel: the kernel sees only L, while the chain aggregation
    refits the chain service law to two moments.

    The tolerance is FIXED at 1e-10 and is deliberately not options.iter_tol: a
    looser one stops the sweep at a different iteration in each codebase.
    """
    import numpy as np

    if phase == 'expand':
        sn, = args
        from ..api.io.model_adapter import ModelAdapter
        from ..api.sn.demands import sn_get_demands_chain
        from ..constants import GlobalConstants
        from ..lang.base import SchedStrategy

        chain_model, alpha, deagg = ModelAdapter.aggregate_chains(solver.model)
        sn_chain = chain_model.getStruct()
        M = int(sn_chain.nstations)
        dem = sn_get_demands_chain(sn_chain)
        Lc = np.atleast_2d(np.asarray(dem.Lchain, dtype=float))
        # The concealment is over CHAINS, and the aggregated model is supposed
        # to carry one class per chain, so the demand matrix's column count and
        # the class count must agree.
        R = int(Lc.shape[1])
        if int(sn_chain.nclasses) != R:
            raise ValueError(
                'the chain-aggregated model has %d classes but %d chains: load concealment '
                'needs one class per chain.' % (int(sn_chain.nclasses), R))
        Nc = np.asarray(dem.Nchain, dtype=float).ravel()
        is_delay = np.asarray(
            [int(sn_chain.sched[i]) == int(SchedStrategy.INF) for i in range(M)], dtype=bool)
        # The concealment slows QUEUEING stations only: a delay holds no queue,
        # so zeroing its rows makes A come out as exactly 1 there.
        L = Lc.copy()
        L[is_delay, :] = 0.0
        Z = Lc[is_delay, :].sum(axis=0) if is_delay.any() else np.zeros(R)

        # One single-class model per chain, built ONCE and re-concealed in place.
        orig_classes = list(chain_model._classes)
        submodels, base_service = [], []
        for l in range(R):
            m = chain_model
            for k in range(R):
                if k != l:
                    m = ModelAdapter.remove_class(m, orig_classes[k])
            submodels.append(m)
            # Only the queueing stations are ever concealed, so only their
            # service law is retained.
            base_service.append([None if is_delay[i]
                                 else m._stations[i].getService(m._classes[0])
                                 for i in range(M)])

        # ONE CLASS PER CHAIN is the ordinary case for load concealment, and
        # there the chain aggregation is the IDENTITY: aggregate_chains returns
        # no deaggregation tables, because the chain answer already IS the class
        # answer. The lift then just re-indexes chains onto their classes.
        ctx = {
            'sn_orig': sn, 'alpha': alpha, 'deagg': deagg,
            'L': L, 'N': Nc, 'Z': Z, 'M': M, 'R': R,
            'identity': int(sn.nchains) >= int(sn.nclasses), 'Korig': int(sn.nclasses),
            'is_delay': is_delay, 'base_service': base_service,
            'X': _lc_seed(L, Nc, Z), 'Xold': -np.ones(R),
            'tol': 1e-10, 'fine_tol': GlobalConstants.FineTol,
            'iterated': True,
        }
        return _lc_conceal_all(ctx, submodels), ctx

    if phase == 'couple':
        # GAUSS-SEIDEL: chain e's throughput is published the moment it is known,
        # so chain e+1 of this same sweep already sees it.
        ctx, submodels, res, e = args
        ctx['X'][e] = _lc_chain_tput(res[e])
        return _lc_conceal_all(ctx, submodels), ctx

    if phase == 'converged':
        ctx, res, it = args
        done = bool(np.max(np.abs(ctx['X'] - ctx['Xold']))
                    <= ctx['tol'] * max(1.0, float(np.max(np.abs(ctx['X'])))))
        ctx['Xold'] = ctx['X'].copy()
        return done, ctx

    if phase == 'lift':
        ctx, res = args
        from ..api.sn.deaggregate import sn_deaggregate_chain_results
        M, R = ctx['M'], ctx['R']

        def col(a):
            a = np.atleast_2d(np.asarray(a, dtype=float))
            out = np.zeros(M)
            n = min(M, a.shape[0])
            out[:n] = a[:n, 0]
            return out

        Q = np.column_stack([col(res[l]['Q']) for l in range(R)])
        U = np.column_stack([col(res[l]['U']) for l in range(R)])
        Rr = np.column_stack([col(res[l]['R']) for l in range(R)])
        T = np.column_stack([col(res[l]['T']) for l in range(R)])
        X = np.array([_lc_chain_tput(res[l]) for l in range(R)], dtype=float)

        if ctx['identity']:
            K = ctx['Korig']
            Qo = np.zeros((M, K)); Uo = np.zeros((M, K))
            Ro = np.zeros((M, K)); To = np.zeros((M, K)); Xo = np.zeros(K)
            for c in range(R):
                k = int(np.atleast_1d(ctx['sn_orig'].inchain[c]).ravel()[0])
                Qo[:, k] = Q[:, c]; Uo[:, k] = U[:, c]
                Ro[:, k] = Rr[:, c]; To[:, k] = T[:, c]; Xo[k] = X[c]
            return _LcIdentityResult(Qo, Uo, Ro, To, Ro.sum(axis=0), Xo)

        deagg = ctx['deagg']
        return sn_deaggregate_chain_results(
            ctx['sn_orig'], deagg.L_chain, None, deagg.ST_chain, deagg.V_chain,
            ctx['alpha'], Q, U, Rr, T, None, X)

    raise ValueError('Unknown load-concealment phase: %s' % phase)
