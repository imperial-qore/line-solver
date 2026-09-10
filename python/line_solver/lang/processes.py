"""
Markov process classes (pure Python).

Native ports of MATLAB's lang/processes MarkovProcess and MarkovChain:
container classes for a CTMC/DTMC with optional state space labels,
supporting steady-state solution (numeric and sympy-symbolic), time
reversal, uniformization, random generation, trajectory sampling,
estimation from sampled system trajectories, and plotting.
"""

import numpy as np

from ..api.mc.ctmc import (
    ctmc_makeinfgen,
    ctmc_solve,
    ctmc_solve_reducible,
    ctmc_timereverse,
    issym,
)
from ..api.mc.dtmc import (
    dtmc_makestochastic,
    dtmc_solve,
    dtmc_solve_reducible,
    dtmc_timereverse,
)
from ..api.io.logging import line_error
from ..lib.kpctoolbox.mc import (
    ctmc_rand,
    ctmc_randomization,
    ctmc_simulate,
    dtmc_rand,
)


def _sym_makestochastic(P):
    """Row-normalize a sympy matrix to a stochastic matrix."""
    import sympy
    M = P.copy() if isinstance(P, sympy.MatrixBase) else sympy.Matrix(P.tolist())
    n = M.rows
    for i in range(n):
        row_sum = sum(M[i, j] for j in range(M.cols))
        if row_sum != 0:
            for j in range(M.cols):
                M[i, j] = M[i, j] / row_sum
    return M


def _joint_sample_state(sa):
    """Joint time-aligned state trajectory and event times from a sampled
    system trajectory (SampleResult, equivalent object, or dict).

    Accepts .state either as the joint state matrix or as a list of per-node
    time-aligned matrices (MATLAB sa.state cells), which are concatenated
    column-wise.
    """
    if isinstance(sa, dict):
        state = sa['state']
        t = sa['t']
    else:
        state = sa.state
        t = sa.t
    if isinstance(state, (list, tuple)):
        state = np.hstack([np.atleast_2d(np.asarray(s, dtype=float))
                           for s in state])
    else:
        state = np.atleast_2d(np.asarray(state, dtype=float))
    t = np.asarray(t, dtype=float).flatten()
    return state, t


def _plot_chain(matrix, state_space, fmt):
    """Plot a labeled transition digraph with matplotlib (circular layout)."""
    import matplotlib.pyplot as plt
    A = np.asarray(matrix, dtype=np.float64)
    n = A.shape[0]
    if state_space is not None and np.size(state_space) > 0:
        space = np.atleast_2d(np.asarray(state_space))
        node_lbl = [','.join(str(int(v)) if float(v).is_integer() else str(v)
                             for v in space[s, :]) for s in range(n)]
    else:
        node_lbl = [str(s + 1) for s in range(n)]
    theta = 2 * np.pi * np.arange(n) / max(n, 1)
    xs, ys = np.cos(theta), np.sin(theta)
    fig, ax = plt.subplots()
    ax.scatter(xs, ys, s=600, facecolors='white', edgecolors='black', zorder=2)
    for s in range(n):
        ax.annotate(node_lbl[s], (xs[s], ys[s]), ha='center', va='center',
                    zorder=3)
    for i in range(n):
        for j in range(n):
            if i != j and A[i, j] != 0:
                dx, dy = xs[j] - xs[i], ys[j] - ys[i]
                ax.annotate('', xytext=(xs[i] + 0.12 * dx, ys[i] + 0.12 * dy),
                            xy=(xs[j] - 0.12 * dx, ys[j] - 0.12 * dy),
                            arrowprops=dict(arrowstyle='-|>', color='black',
                                            connectionstyle='arc3,rad=0.1'),
                            zorder=1)
                ax.annotate(fmt % A[i, j],
                            (xs[i] + 0.4 * dx + 0.08 * dy,
                             ys[i] + 0.4 * dy - 0.08 * dx),
                            ha='center', va='center', fontsize=8, zorder=3)
    ax.set_aspect('equal')
    ax.axis('off')
    plt.show()


class MarkovProcess:
    """A continuous-time Markov chain (CTMC).

    Args:
        infgen: Infinitesimal generator matrix (numpy array or sympy matrix);
            the diagonal is corrected via ctmc_makeinfgen
        isfinite: True if the state space is finite (default True)
        state_space: Optional state space matrix, one row per state
    """

    def __init__(self, infgen, isfinite=True, state_space=None):
        self.name = 'MarkovProcess'
        self.infGen = ctmc_makeinfgen(infgen)
        self.isfinite = isfinite
        self.stateSpace = state_space

    def toMarkovChain(self, q=None):
        """Uniformized DTMC: P = Q/q + I with default q = max|Q| + rand."""
        if issym(self.infGen):
            if q is None:
                line_error('toMarkovChain',
                           'A uniformization constant q must be provided for '
                           'symbolic CTMCs.')
            import sympy
            P = self.infGen / q + sympy.eye(self.infGen.rows)
        else:
            P, q = ctmc_randomization(np.asarray(self.infGen, dtype=np.float64), q)
        A = MarkovChain(P)
        A.setStateSpace(self.stateSpace)
        return A

    def toDTMC(self, q=None):
        """Alias for toMarkovChain for backwards compatibility."""
        return self.toMarkovChain(q)

    def toTimeReversed(self):
        """Time-reversed CTMC: Q*_ij = pi_j Q_ji / pi_i."""
        if issym(self.infGen):
            import sympy
            pi = ctmc_solve(self.infGen)
            n = self.infGen.rows
            Qrev = sympy.zeros(n, n)
            for i in range(n):
                for j in range(n):
                    Qrev[i, j] = pi[j] * self.infGen[j, i] / pi[i]
            return MarkovProcess(Qrev)
        return MarkovProcess(ctmc_timereverse(self.infGen))

    def setStateSpace(self, state_space):
        self.stateSpace = state_space

    def getGenerator(self):
        """Return the infinitesimal generator matrix."""
        return self.infGen

    def getProbState(self, state):
        """Probability of a single state by Cramer's rule.

        Returns (pi_i, num, den) with pi_i = det(Q_i)/det(Q), where Q has its
        first column replaced by ones (normalization) and Q_i additionally has
        row i replaced by the first unit row.
        """
        from ..api.pfqn.utils import matchrow
        i = matchrow(np.atleast_2d(np.asarray(self.stateSpace)),
                     np.asarray(state).flatten())
        if i == 0:
            line_error('getProbState', 'State not found in the state space.')
        i = i - 1  # matchrow is 1-based
        if issym(self.infGen):
            import sympy
            Q = self.infGen.copy()
            n = Q.rows
            for r in range(n):
                Q[r, 0] = 1
            Q_i = Q.copy()
            for c in range(n):
                Q_i[i, c] = 0
            Q_i[i, 0] = 1
            num = Q_i.det()
            den = Q.det()
            pi_i = sympy.simplify(num / den)
        else:
            Q = np.asarray(self.infGen, dtype=np.float64).copy()
            Q[:, 0] = 1
            Q_i = Q.copy()
            Q_i[i, :] = 0
            Q_i[i, 0] = 1
            num = np.linalg.det(Q_i)
            den = np.linalg.det(Q)
            pi_i = num / den
        return pi_i, num, den

    def solve(self):
        """Steady-state distribution (exact rational for symbolic input)."""
        if issym(self.infGen):
            return ctmc_solve(self.infGen)
        return ctmc_solve_reducible(self.infGen)

    def transient(self, pi0=None, t=1.0, method='unif'):
        """Distribution at time t from pi0 (uniform if None).

        method: 'unif' for Jensen uniformization (default), 'foxglynn' for the
        Fox-Glynn weights, which avoid evaluating the Poisson terms directly.
        """
        Q = np.asarray(self.infGen, dtype=np.float64)
        pi0 = self._init_or_uniform(pi0)
        if str(method).lower() == 'foxglynn':
            from ..api.mc.foxglynn import ctmc_foxglynn
            return ctmc_foxglynn(pi0, Q, t)
        from ..lib.kpctoolbox.mc import ctmc_uniformization
        pi_t, _ = ctmc_uniformization(pi0, Q, t)
        return pi_t

    def solveRelative(self, refstate=0):
        """Equilibrium distribution relative to refstate, i.e. with p(refstate)=1.

        Unnormalized by construction, so it is defined even where the
        normalizing constant is not; refstate is 0-based here and 1-based in
        MATLAB, as elsewhere between the two codebases.
        """
        from ..lib.kpctoolbox.mc import ctmc_relsolve
        return ctmc_relsolve(np.asarray(self.infGen, dtype=np.float64), refstate)

    def aggregate(self, MS, method='courtois', param=None):
        """Aggregation-disaggregation over the macrostate partition MS.

        MS is a list of lists of 0-based state indices, one per macrostate.
        method: 'courtois' (nearly-completely-decomposable approximation, param
        is the randomization rate q), 'kms' (Koury-McAllister-Stewart) or
        'takahashi' (param is the iteration count, default 10), or 'multi'
        (param is the second-level partition MSS).

        Returns (p, eps, epsMAX): the approximate stationary vector, the
        nearly-complete-decomposability index of the partition, and the largest
        index for which the approximation is meant to hold.
        """
        from ..api.mc.aggregation import (ctmc_courtois, ctmc_kms, ctmc_multi,
                                          ctmc_takahashi)
        Q = np.asarray(self.infGen, dtype=np.float64)
        method = str(method).lower()
        if method == 'courtois':
            res = ctmc_courtois(Q, MS, param)
        elif method == 'kms':
            res = ctmc_kms(Q, MS, 10 if param is None else int(param))
        elif method == 'takahashi':
            res = ctmc_takahashi(Q, MS, 10 if param is None else int(param))
        elif method == 'multi':
            if param is None:
                line_error('aggregate', "The 'multi' method requires the second-level partition MSS.")
            res = ctmc_multi(Q, MS, param)
        else:
            line_error('aggregate', "Unknown aggregation method '%s'." % method)
        return res.p, res.eps, res.epsMAX

    # Alias, under the name the JAR must use since 'transient' is a Java keyword
    transientProb = transient

    def timeAverage(self, pi0=None, t=1.0):
        """Time-averaged distribution over [0,t] and its endpoint."""
        from ..api.mc.ctmc import ctmc_timeaverage
        Q = np.asarray(self.infGen, dtype=np.float64)
        out = ctmc_timeaverage(self._init_or_uniform(pi0), Q, t)
        return out[0], out[1]

    def sens(self, dQ):
        """Sensitivity of the stationary distribution to a scalar parameter."""
        from ..api.mc.ctmc import ctmc_sens
        return ctmc_sens(np.asarray(self.infGen, dtype=np.float64),
                         np.asarray(dQ, dtype=np.float64), self.solve())

    def stochComp(self, I=None):
        """Stochastic complement of the states I, a generator on that subset.

        Returns the complement itself; stochCompFull additionally returns the
        four blocks of the partitioned generator and the return-path term T.
        """
        from ..api.mc.ctmc import ctmc_stochcomp
        return ctmc_stochcomp(np.asarray(self.infGen, dtype=np.float64), I)['S']

    def stochCompFull(self, I=None):
        """Stochastic complement of the states I with its blocks.

        Returns a dict with the complement 'S', the blocks 'Q11', 'Q12', 'Q21',
        'Q22' of the generator partitioned by I and its complement, and the
        return-path term 'T' = Q12*inv(-Q22)*Q21, so that S = Q11 + T. Twin of
        the MATLAB [S,Q11,Q12,Q21,Q22,T] = ctmc.stochCompFull(I) and of the JAR
        stochCompFull, which return the same six matrices.
        """
        from ..api.mc.ctmc import ctmc_stochcomp
        return ctmc_stochcomp(np.asarray(self.infGen, dtype=np.float64), I)

    def isFeasible(self):
        """True when the generator is a valid one."""
        from ..api.mc.ctmc import ctmc_isfeasible
        return bool(ctmc_isfeasible(np.asarray(self.infGen, dtype=np.float64)))

    def toEmbedded(self):
        """Embedded jump chain, i.e. the DTMC of the states visited at
        transition epochs.

        Unlike toDTMC (uniformization) it does not preserve the stationary
        distribution, since it drops the holding times; an absorbing state
        stays absorbing.
        """
        Q = np.asarray(self.infGen, dtype=np.float64)
        n = Q.shape[0]
        exit_rate = -np.diag(Q)
        P = Q - np.diag(np.diag(Q))
        for i in range(n):
            if exit_rate[i] > 0:
                P[i, :] = P[i, :] / exit_rate[i]
            else:
                P[i, i] = 1.0
        A = MarkovChain(P)
        A.setStateSpace(self.stateSpace)
        return A

    def _init_or_uniform(self, pi0):
        """The given initial distribution as a row vector, or the uniform one."""
        n = np.asarray(self.infGen).shape[0]
        if pi0 is None or np.asarray(pi0).size != n:
            return np.ones(n) / n
        return np.asarray(pi0, dtype=np.float64).flatten()

    def sample(self, n=1, seed=None):
        """Sample a trajectory of n transitions.

        The initial state is drawn from a random distribution, as in MATLAB's
        ctmc_simulate default. Returns (sojourn_times, states), as in MATLAB.
        """
        if issym(self.infGen):
            line_error('sample',
                       'MarkovProcess.sample does not support symbolic CTMCs.')
        Q = np.asarray(self.infGen, dtype=np.float64)
        return ctmc_simulate(Q, None, n, seed=seed)

    def plot(self):
        """Plot the transition digraph with rate-labeled edges."""
        if issym(self.infGen):
            line_error('plot',
                       'CTMC.plot does not support symbolic CTMCs.')
        Q0 = np.asarray(self.infGen, dtype=np.float64).copy()
        np.fill_diagonal(Q0, 0.0)
        _plot_chain(Q0, self.stateSpace, '%.4f')

    @staticmethod
    def rand(n_states):
        """Create a random CTMC with n_states states."""
        return MarkovProcess(ctmc_rand(n_states))

    @staticmethod
    def fromSampleSysAggr(sa):
        """Estimate a CTMC from a sampled system trajectory.

        Transition counts between consecutive distinct joint states are
        normalized into a DTMC and divided by the empirical mean holding
        times, as in MATLAB MarkovProcess.fromSampleSysAggr.
        """
        state, t = _joint_sample_state(sa)
        state_space, state_hash = np.unique(state, axis=0, return_inverse=True)
        m = state_space.shape[0]
        counts = np.zeros((m, m))
        hold_time = np.zeros(m)
        for i in range(1, len(state_hash)):
            counts[state_hash[i - 1], state_hash[i]] += 1
            hold_time[state_hash[i - 1]] += t[i] - t[i - 1]
        with np.errstate(divide='ignore', invalid='ignore'):
            hold_time = hold_time / counts.sum(axis=1)
            inf_gen = ctmc_makeinfgen(
                dtmc_makestochastic(counts) / hold_time[:, None])
        return MarkovProcess(inf_gen, True, state_space)


class MarkovChain:
    """A discrete-time Markov chain (DTMC).

    Args:
        trans_mat: Transition matrix (numpy array or sympy matrix); rows are
            normalized via dtmc_makestochastic
        isfinite: True if the state space is finite (default True)
    """

    def __init__(self, trans_mat, isfinite=True):
        self.name = 'MarkovChain'
        if issym(trans_mat):
            self.transMat = _sym_makestochastic(trans_mat)
        else:
            self.transMat = dtmc_makestochastic(trans_mat)
        self.stateSpace = None
        self.isfinite = isfinite

    def toMarkovProcess(self):
        """Embedded CTMC: Q = P - I."""
        if issym(self.transMat):
            import sympy
            Q = self.transMat - sympy.eye(self.transMat.rows)
        else:
            Q = np.asarray(self.transMat, dtype=np.float64) - \
                np.eye(self.transMat.shape[0])
        A = MarkovProcess(Q)
        A.setStateSpace(self.stateSpace)
        return A

    def toCTMC(self):
        """Alias for toMarkovProcess for backwards compatibility."""
        return self.toMarkovProcess()

    def toTimeReversed(self):
        """Time-reversed DTMC: P*_ij = pi_j P_ji / pi_i."""
        if issym(self.transMat):
            import sympy
            n = self.transMat.rows
            # exact stationary vector via the embedded symbolic CTMC
            pi = ctmc_solve(self.transMat - sympy.eye(n))
            Prev = sympy.zeros(n, n)
            for i in range(n):
                for j in range(n):
                    Prev[i, j] = pi[j] * self.transMat[j, i] / pi[i]
            return MarkovChain(Prev)
        return MarkovChain(dtmc_timereverse(self.transMat))

    def getTransMat(self):
        """Return the transition matrix."""
        return self.transMat

    def solve(self):
        """Stationary distribution of the DTMC. Twin of MarkovProcess.solve."""
        if issym(self.transMat):
            return dtmc_solve(self.transMat)
        return dtmc_solve_reducible(self.transMat)

    def transient(self, pi0=None, steps=1):
        """Distribution at each step 0,...,steps from pi0 (uniform if None)."""
        from ..api.mc.dtmc import dtmc_transient
        P = np.asarray(self.transMat, dtype=np.float64)
        if pi0 is None or np.asarray(pi0).size != P.shape[0]:
            pi0 = np.ones(P.shape[0]) / P.shape[0]
        return dtmc_transient(P, np.asarray(pi0, dtype=np.float64).flatten(), int(steps))

    # Alias, under the name the JAR must use since 'transient' is a Java keyword
    transientProb = transient

    def hittingTime(self, target_states):
        """Mean number of steps to reach any state in target_states."""
        from ..api.mc.dtmc import dtmc_hitting_time
        return dtmc_hitting_time(np.asarray(self.transMat, dtype=np.float64), target_states)

    def stochComp(self, keep_states):
        """Stochastic complement of the kept states, a DTMC on that subset."""
        from ..api.mc.dtmc import dtmc_stochcomp
        return dtmc_stochcomp(np.asarray(self.transMat, dtype=np.float64), keep_states)

    def stochCompFull(self, keep_states):
        """Stochastic complement of the kept states with its blocks.

        Returns a dict with the complement 'S' and the blocks 'P11', 'P12',
        'P21', 'P22' of the transition matrix partitioned by the kept and the
        eliminated states. Twin of the MATLAB
        [S,P11,P12,P21,P22] = dtmc.stochCompFull(I) and of the JAR
        stochCompFull.
        """
        from ..api.mc.dtmc import dtmc_stochcomp_full
        return dtmc_stochcomp_full(np.asarray(self.transMat, dtype=np.float64), keep_states)

    def transientUnif(self, pi0=None, t=1.0):
        """Distribution at time t of the DTMC seen through uniformization.

        The chain is read as the randomized image of a CTMC, so t is continuous
        here, unlike the step count taken by transient. Wraps dtmc_uniformization.
        """
        from ..lib.kpctoolbox.mc import dtmc_uniformization
        P = np.asarray(self.transMat, dtype=np.float64)
        if pi0 is None or np.asarray(pi0).size != P.shape[0]:
            pi0 = np.ones(P.shape[0]) / P.shape[0]
        pi_t, _ = dtmc_uniformization(np.asarray(pi0, dtype=np.float64).flatten(), P, t)
        return pi_t

    def isFeasible(self):
        """True when the transition matrix is stochastic."""
        from ..api.mc.dtmc import dtmc_isfeasible
        return bool(dtmc_isfeasible(np.asarray(self.transMat, dtype=np.float64)))

    def sample(self, n=1, seed=None):
        """Simulate n steps of the DTMC from a random initial state.

        Returns the sequence of visited state indices (0-based).
        """
        if issym(self.transMat):
            line_error('sample',
                       'MarkovChain.sample does not support symbolic chains.')
        if seed is not None:
            np.random.seed(seed)
        from ..lib.kpctoolbox.mc import dtmc_simulate
        P = np.asarray(self.transMat, dtype=np.float64)
        pi0 = np.random.rand(P.shape[0])
        pi0 = pi0 / pi0.sum()
        return dtmc_simulate(P, pi0, n)

    def setStateSpace(self, state_space):
        self.stateSpace = state_space

    def plot(self):
        """Plot the transition digraph with probability-labeled edges."""
        if issym(self.transMat):
            line_error('plot',
                       'DTMC.plot does not support symbolic chains.')
        _plot_chain(self.transMat, self.stateSpace, '%.2f')

    @staticmethod
    def rand(n_states):
        """Create a random DTMC with n_states states."""
        return MarkovChain(dtmc_rand(n_states))

    @staticmethod
    def fromSampleSysAggr(sa):
        """Estimate a DTMC from a sampled system trajectory.

        Transition counts between consecutive joint states are normalized
        into a stochastic matrix, as in MATLAB MarkovChain.fromSampleSysAggr.
        """
        state, _ = _joint_sample_state(sa)
        state_space, state_hash = np.unique(state, axis=0, return_inverse=True)
        m = state_space.shape[0]
        counts = np.zeros((m, m))
        for i in range(1, len(state_hash)):
            counts[state_hash[i - 1], state_hash[i]] += 1
        dtmc_obj = MarkovChain(dtmc_makestochastic(counts), True)
        dtmc_obj.setStateSpace(state_space)
        return dtmc_obj
