"""Exact reference for the renv_fourstages_repairmen random-environment model.

SolverENV is an approximation, and until 2026-07-31 the only thing its results
were checked against was the OTHER codebase's approximation: the parity harness
carried a 5% tolerance on this model and nobody had an exact number.

The model is a closed 2-station network, N = 30: a Delay whose exponential rate
is modulated by a 4-state environment (rates 4, 3, 2, 1 across stages) and a PS
Queue with constant rate 1 in every stage. Because the queue's rate never varies
and N is large, the queue is saturated throughout.

THE ENVIRONMENT IS NOT EXPONENTIAL. The example declares each transition as
`APH.fit_mean_and_scv(1/rate, 0.5)`, and scv = 1/2 with mean 2 is Erlang-2 with
unit phase rates. Stage e races one such clock per outgoing transition, all
reset on entry, so the holding time is a minimum of Erlang-2 variables and the
process is semi-Markov, not a CTMC over the four stages. Treating it as
exponential is the easy mistake here and it moves the answer by 3%: it gives
0.419934 instead of the 0.432291 below.

Expanding the Erlang phases restores a CTMC: 14 environment states (2 for the
one-clock stage, 4 for each two-clock stage) x 31 job counts = 434 states, which
this test solves directly. Two independent checks confirm the construction:

  * the stage marginals come out as [4, 4, 2, 3]/13, matching MATLAB's
    `Environment.probEnv` to six digits;
  * that same vector follows from the semi-Markov formula pi ~ pj / lambda with
    the embedded chain pj = [5, 8, 4, 6]/23 and lambda = [0.5, 0.8, 0.8, 0.8],
    where 1/0.8 = 1.25 = E[min of two iid Erlang-2(1)] exactly.

DO NOT use a MAP on the Delay as an oracle for this model. Encoding the
modulation as D0 = Q_env - diag(mu), D1 = diag(mu) looks equivalent and is not:
a MAP restarts each service from the MAP's ARRIVAL-stationary phase vector,
biased toward the fast phases because completions are more frequent there,
rather than from the time-stationary one. LDES simulates that encoding
faithfully and returns the encoding's own answer, which is not this model's.
LDES cannot express a shared environment.

Measured against the exact values below on 2026-07-31:
  MATLAB SolverENV   Delay QLen 0.444640 (+2.9%), Tput 0.977960 (exact is 1)
  python  SolverENV  Delay QLen 0.487179 (+12.7%), Tput 1.000000
Neither is asserted here: both are approximations, and pinning either would bake
in its error. MATLAB is much the closer of the two on queue length but is the
only one wrong on throughput, where a closed tandem must carry equal flow at
both stations. Improve the blends against THESE numbers, not against each other.
"""
import numpy as np

N = 30
MU = np.array([4.0, 3.0, 2.0, 1.0])   # Delay service rate per stage
NU = 1.0                              # PS queue rate, the same in every stage
TARGETS = {0: [1], 1: [2, 3], 2: [0, 3], 3: [0, 1]}   # nonzero envRates entries
PHASE_RATE = 1.0                      # Erlang-2 with mean 2

EXACT_DELAY_QLEN = 0.432291
EXACT_TPUT = 1.0
EXACT_STAGE_PROB = np.array([4.0, 4.0, 2.0, 3.0]) / 13.0


def _env_states():
    """(stage, per-clock Erlang phase); every clock resets when a stage begins."""
    states = []
    for e in range(4):
        for ph in np.ndindex(*([2] * len(TARGETS[e]))):
            states.append((e, tuple(ph)))
    return states


def _build_generator():
    env = _env_states()
    eidx = {s: i for i, s in enumerate(env)}
    E = len(env)
    entry = lambda e: eidx[(e, tuple([0] * len(TARGETS[e])))]
    size = (N + 1) * E
    idx = lambda n, s: n * E + s
    gen = np.zeros((size, size))
    for n in range(N + 1):
        for s, (e, ph) in enumerate(env):
            i = idx(n, s)
            if n > 0:
                gen[i, idx(n - 1, s)] += n * MU[e]    # INF servers at the Delay
            if N - n > 0:
                gen[i, idx(n + 1, s)] += NU           # single PS server, rate 1
            for k, tgt in enumerate(TARGETS[e]):
                if ph[k] == 0:                        # Erlang phase 1 -> 2
                    nph = list(ph)
                    nph[k] = 1
                    gen[i, idx(n, eidx[(e, tuple(nph))])] += PHASE_RATE
                else:                                 # clock fires, stage changes
                    gen[i, idx(n, entry(tgt))] += PHASE_RATE
    np.fill_diagonal(gen, -gen.sum(1))
    return gen, env, idx


def _stationary(gen):
    n = gen.shape[0]
    rhs = np.zeros(n + 1)
    rhs[-1] = 1.0
    return np.linalg.lstsq(np.vstack([gen.T, np.ones(n)]), rhs, rcond=None)[0]


def test_joint_chain_mean_queue_length_and_throughput():
    gen, env, idx = _build_generator()
    p = _stationary(gen)
    qlen = p @ np.repeat(np.arange(N + 1), len(env))
    busy = sum(p[idx(n, s)] for n in range(N) for s in range(len(env)))
    assert abs(qlen - EXACT_DELAY_QLEN) < 1e-6
    assert abs(NU * busy - EXACT_TPUT) < 1e-6


def test_stage_marginals_match_probenv():
    """The blend weights MATLAB computes are correct; the error is elsewhere."""
    gen, env, idx = _build_generator()
    p = _stationary(gen)
    pe = np.zeros(4)
    for s, (e, _) in enumerate(env):
        pe[e] += sum(p[idx(n, s)] for n in range(N + 1))
    np.testing.assert_allclose(pe, EXACT_STAGE_PROB, rtol=1e-6)


def test_semi_markov_derivation_of_stage_probabilities():
    """pi ~ pj / lambda, with lambda set by E[min of two iid Erlang-2(1)] = 5/4."""
    P = np.array([[0, 1, 0, 0], [0, 0, .5, .5], [.5, 0, 0, .5], [.5, .5, 0, 0]])
    rhs = np.zeros(5)
    rhs[-1] = 1.0
    pj = np.linalg.lstsq(np.vstack([P.T - np.eye(4), np.ones(4)]), rhs,
                         rcond=None)[0]
    lam = np.array([0.5, 0.8, 0.8, 0.8])
    pi = pj / lam
    pi /= pi.sum()
    np.testing.assert_allclose(pi, EXACT_STAGE_PROB, rtol=1e-6)
