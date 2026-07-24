
import numpy as np
from scipy.special import gammaln, logsumexp

def nchoosekln(n, k):
    if k < 0 or k > n:
        raise ValueError("Invalid k")
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)

def factln(n):
    if n < 0:
        raise ValueError("Invalid n")
    return gammaln(n + 1)

def exp_stable(x):
    if x < -700:
        return 0.0
    if x > 700:
        return np.inf
    return np.exp(x)

def log_sum_exp(x):
    if len(x) == 0:
        return -np.inf
    return logsumexp(x)

class ProbabilityResults:
    def __init__(self, logNormConstAggr=np.nan, sys_aggr_prob=None, computation_method="binomial_approx"):
        self.logNormConstAggr = logNormConstAggr
        self.sys_aggr_prob = sys_aggr_prob
        self.computation_method = computation_method
        self.marginal_probs = {}
        self.aggr_probs = {}
        self._aggr_cache = {}

    def has_marginal(self, ist, r):
        return (ist, r) in self.marginal_probs

    def get_marginal(self, ist, r):
        return self.marginal_probs.get((ist, r))

    def set_marginal(self, ist, r, probs):
        self.marginal_probs[(ist, r)] = probs

    def clear(self):
        self.logNormConstAggr = np.nan
        self.sys_aggr_prob = None
        self.marginal_probs = {}
        self._aggr_cache = {}

class SolverResults:
    def __init__(self, Q=None, U=None, R=None, T=None, X=None, runtime=0.0, method="", prob=None):
        self.Q = Q
        self.U = U
        self.R = R
        self.T = T
        self.X = X
        self.runtime = runtime
        self.method = method
        self.prob = prob if prob else ProbabilityResults()

    def is_complete(self):
        return self.Q is not None and self.U is not None and self.R is not None and self.T is not None and self.X is not None

    def clear(self):
        self.Q = None
        self.U = None
        self.runtime = 0.0
        self.method = ""

    def get_queue_length(self, i, r):
        return self.Q[i-1, r-1] if self.Q is not None else 0

    def get_utilization(self, i, r):
        return self.U[i-1, r-1] if self.U is not None else 0

    def get_response_time(self, i, r):
        return self.R[i-1, r-1] if self.R is not None else 0

    def get_throughput(self, i, r):
        return self.T[i-1, r-1] if self.T is not None else 0

    def get_sys_throughput(self, r):
        return self.X[r-1] if self.X is not None else 0

def schmidt_binomial_prob_aggr(Q, N, nir, i):
    """
    Compute aggregate probability using Schmidt 1997 binomial approximation.

    Computes P(n_i) = product over classes r of P(n_ir | N_r)
    where P(n_ir | N_r) = C(N_r, n_ir) * p_ir^n_ir * (1-p_ir)^(N_r - n_ir)
    and p_ir = Q[i,r] / N[r]

    Args:
        Q: Mean queue length matrix [stations x classes]
        N: Population vector [classes]
        nir: State vector at station i [classes] - number of jobs per class
        i: Station index

    Returns:
        (log_prob, prob): Log probability and probability
    """
    Q = np.atleast_2d(Q)
    N = np.atleast_1d(N)
    nir = np.atleast_1d(nir)

    num_classes = len(N)
    log_prob = 0.0

    for r in range(num_classes):
        N_r = int(N[r])
        n_r = int(nir[r])

        # Check validity
        if n_r < 0 or n_r > N_r:
            return -np.inf, 0.0

        # Compute probability p = Q[i,r] / N[r]
        if N_r == 0:
            if n_r == 0:
                continue  # P(0 | 0) = 1, log(1) = 0
            else:
                return -np.inf, 0.0

        p = Q[i, r] / N_r if N_r > 0 else 0.0

        # Clamp p to [0, 1]
        p = max(0.0, min(1.0, p))

        # Compute log binomial probability
        # log(C(N,n) * p^n * (1-p)^(N-n))
        log_binom = nchoosekln(N_r, n_r)

        if n_r > 0:
            if p <= 0:
                return -np.inf, 0.0
            log_binom += n_r * np.log(p)

        if n_r < N_r:
            if p >= 1:
                return -np.inf, 0.0
            log_binom += (N_r - n_r) * np.log(1 - p)

        log_prob += log_binom

    prob = exp_stable(log_prob)
    return log_prob, prob


def schmidt_binomial_prob_marg(Q, N, i, r):
    """
    Compute marginal distribution of class r at station i using binomial approximation.

    Args:
        Q: Mean queue length matrix [stations x classes]
        N: Population vector [classes]
        i: Station index
        r: Class index

    Returns:
        (states, probs): Arrays of possible states [0, 1, ..., N_r] and their probabilities
    """
    Q = np.atleast_2d(Q)
    N = np.atleast_1d(N)

    N_r = int(N[r])

    if N_r == 0:
        return np.array([0]), np.array([1.0])

    # Compute probability p = Q[i,r] / N[r]
    p = Q[i, r] / N_r if N_r > 0 else 0.0
    p = max(0.0, min(1.0, p))

    states = np.arange(N_r + 1)
    probs = np.zeros(N_r + 1)

    for n in range(N_r + 1):
        # log(C(N,n) * p^n * (1-p)^(N-n))
        log_binom = nchoosekln(N_r, n)

        if n > 0:
            if p <= 0:
                probs[n] = 0.0
                continue
            log_binom += n * np.log(p)

        if n < N_r:
            if p >= 1:
                probs[n] = 0.0
                continue
            log_binom += (N_r - n) * np.log(1 - p)

        probs[n] = exp_stable(log_binom)

    # Normalize to ensure sum = 1 (handles numerical issues)
    total = np.sum(probs)
    if total > 0:
        probs = probs / total

    return states, probs


def schmidt_binomial_prob_sys(Q, N, state):
    """
    Compute joint system state probability using binomial approximation.

    Assumes independence across stations (product form approximation).

    Args:
        Q: Mean queue length matrix [stations x classes]
        N: Population vector [classes]
        state: State matrix [stations x classes] - number of jobs at each station/class

    Returns:
        (log_prob, prob): Log probability and probability
    """
    Q = np.atleast_2d(Q)
    N = np.atleast_1d(N)
    state = np.atleast_2d(state)

    num_stations, num_classes = Q.shape

    # Validate state dimensions
    if state.shape[1] != num_classes:
        raise ValueError(f"State has {state.shape[1]} classes but Q has {num_classes}")

    log_prob = 0.0

    # For each station, compute the binomial probability
    for i in range(state.shape[0]):
        for r in range(num_classes):
            N_r = int(N[r])
            n_r = int(state[i, r])

            # Check validity
            if n_r < 0 or n_r > N_r:
                return -np.inf, 0.0

            if N_r == 0:
                if n_r == 0:
                    continue
                else:
                    return -np.inf, 0.0

            # Compute probability p = Q[i,r] / N[r]
            p = Q[i, r] / N_r if N_r > 0 else 0.0
            p = max(0.0, min(1.0, p))

            # Compute log binomial probability
            log_binom = nchoosekln(N_r, n_r)

            if n_r > 0:
                if p <= 0:
                    return -np.inf, 0.0
                log_binom += n_r * np.log(p)

            if n_r < N_r:
                if p >= 1:
                    return -np.inf, 0.0
                log_binom += (N_r - n_r) * np.log(1 - p)

            log_prob += log_binom

    prob = exp_stable(log_prob)
    return log_prob, prob
