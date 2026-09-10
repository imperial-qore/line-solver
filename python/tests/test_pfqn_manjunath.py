"""
Validation of the Manjunath-Sikdar transform for product-form queueing networks
(pfqn_manjunath).

THE TWO ORACLES DO NOT COME OUT OF THE IMPLEMENTATION. With no extra rows the
transform computes the ordinary closed-network normalizing constant, for which
pfqn_ca's convolution recursion is an independent exact algorithm sharing no
code; the two are algebraic identities for the same sum, so agreement to 1e-13 is
the correct expectation and not a tolerance chosen to pass. With extra rows
pfqn_ca has nothing to say, and the oracle becomes `bcmp_enum` below, which sums
the BCMP product form over the enumerated state space and applies each row by
direct comparison -- the very enumeration the transform exists to avoid, so a
coefficient-domain defect cannot hide behind a shared traversal.

The constrained expectations are additionally pinned to the MATLAB reference
implementation's output, reproduced to 12 digits: the transform is exact, so
there is no tolerance to hide behind and any drift is a defect rather than noise.

Ported from matlab/src/api/pfqn/test_pfqn_manjunath.m.
"""

import itertools
import math

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_ca, pfqn_manjunath


def bcmp_enum(L, N, Z, A, b, sense):
    """Direct sum of the BCMP product form over the enumerated closed state
    space, keeping the states that satisfy every extra row."""
    M, R = L.shape
    Mz = Z.shape[0]
    S = M + Mz
    allocs = [[c for c in itertools.product(range(N[r] + 1), repeat=S) if sum(c) == N[r]]
              for r in range(R)]
    g = 0.0
    for combo in itertools.product(*allocs):
        n = np.array(combo).T                      # (S, R)
        ok = True
        for j in range(len(b)):
            v = float(A[j, :] @ n.ravel(order='F'))
            if sense[j] == 'E':
                ok &= abs(v - b[j]) < 1e-9
            elif sense[j] == 'L':
                ok &= v <= b[j] + 1e-9
            else:
                ok &= v > b[j] + 1e-9
        if not ok:
            continue
        t = 1.0
        for i in range(M):
            t *= math.factorial(int(n[i, :].sum()))
            for r in range(R):
                t *= L[i, r] ** n[i, r] / math.factorial(int(n[i, r]))
        for k in range(Mz):
            for r in range(R):
                t *= Z[k, r] ** n[M + k, r] / math.factorial(int(n[M + k, r]))
        g += t
    return g


UNCONSTRAINED = [
    (np.array([[1.0], [2.0]]), [4], None),
    (np.array([[1.0, 2.0], [3.0, 1.0]]), [2, 3], None),
    (np.array([[1.0, 2.0], [3.0, 1.0]]), [2, 3], np.array([[0.5, 1.5]])),
    (np.array([[0.4, 0.2], [0.9, 0.7], [0.1, 1.1]]), [3, 2], np.array([[1.0, 2.0]])),
    (None, [2, 1], np.array([[1.0, 3.0]])),
    (np.array([[1.0, 2.0], [3.0, 1.0]]), [0, 0], np.array([[1.0, 1.0]])),
    (np.array([[5.0, 1.0], [1.0, 6.0]]), [6, 5], np.array([[2.0, 3.0]])),
]


@pytest.mark.parametrize("idx", range(len(UNCONSTRAINED)))
def test_unconstrained_matches_pfqn_ca(idx):
    """With no extra rows the transform IS the closed-network constant."""
    L, N, Z = UNCONSTRAINED[idx]
    R = len(N)
    Lc = np.zeros((0, R)) if L is None else L
    Zc = np.zeros((1, R)) if Z is None else Z
    _, lG_ca = pfqn_ca(Lc, np.array(N, dtype=float), Zc.sum(axis=0))
    _, lG_mj, _ = pfqn_manjunath(L, N, Z)
    assert abs(lG_ca - lG_mj) / max(1.0, abs(lG_ca)) < 1e-13


# The shared constrained fixture: three queueing stations and one delay, two
# classes, so the occupancy is 4-by-2 and every extra row has 8 coefficients.
L_C = np.array([[1.0, 2.0], [3.0, 1.0], [0.5, 0.5]])
N_C = [3, 2]
Z_C = np.array([[1.0, 2.0]])
S_C = 4

_a1 = np.zeros(S_C * 2); _a1[0] = 1; _a1[0 + S_C] = 1          # jobs at queue 1
_a2 = np.zeros(S_C * 2); _a2[0] = 2; _a2[1] = 1                # weighted budget
_a2[0 + S_C] = 1; _a2[1 + S_C] = 3
_a3 = np.zeros(S_C * 2); _a3[3] = 1                            # class 1 at the delay

CONSTRAINED = [
    (_a1[None, :], [2], 'L', 2214.58333333333),
    (_a1[None, :], [2], 'E', 526.541666666667),
    (_a1[None, :], [1], 'G', 1000.29166666667),
    (_a2[None, :], [6], 'L', 1611.83333333333),
    (np.vstack([_a1, _a2]), [2, 6], 'LL', 1326.58333333333),
    (np.vstack([_a1, _a2]), [2, 6], 'EL', 353.041666666667),
    (np.vstack([_a1, _a2]), [1, 6], 'GL', 638.291666666667),
    (np.vstack([_a1, _a2]), [1, 5], 'GG', 559.25),
    (np.vstack([_a1, _a3]), [2, 1], 'LL', 2159.6875),
]


@pytest.mark.parametrize("idx", range(len(CONSTRAINED)))
def test_constrained_matches_enumeration_and_matlab(idx):
    A, b, sense, matlab = CONSTRAINED[idx]
    g_enum = bcmp_enum(L_C, N_C, Z_C, A, np.array(b, dtype=float), sense)
    G, lG, _ = pfqn_manjunath(L_C, N_C, Z_C, A, b, sense)
    assert abs(g_enum - G) / abs(g_enum) < 1e-12
    assert abs(matlab - G) / abs(matlab) < 1e-12
    assert abs(math.exp(lG) - G) / G < 1e-12


def test_population_row_is_redundant():
    """Every admissible state holds sum(N) jobs in total, so declaring that as an
    extra row must change nothing -- a direct check that an equality row is
    discharged by picking a coefficient rather than by summing one."""
    G_ref, _, _ = pfqn_manjunath(L_C, N_C, Z_C)
    total = sum(N_C)
    ones = np.ones((1, S_C * 2))
    G_eq, _, _ = pfqn_manjunath(L_C, N_C, Z_C, ones, [total], 'E')
    G_le, _, _ = pfqn_manjunath(L_C, N_C, Z_C, ones, [total], 'L')
    G_gt, _, _ = pfqn_manjunath(L_C, N_C, Z_C, ones, [total], 'G')
    assert abs(G_eq - G_ref) / G_ref < 1e-13
    assert abs(G_le - G_ref) / G_ref < 1e-13
    assert G_gt == 0.0


def test_trivial_rows_are_decided_not_carried():
    G_ref, _, _ = pfqn_manjunath(L_C, N_C, Z_C)
    z = np.zeros((1, S_C * 2))
    assert abs(pfqn_manjunath(L_C, N_C, Z_C, z, [0], 'E')[0] - G_ref) / G_ref < 1e-13
    assert pfqn_manjunath(L_C, N_C, Z_C, z, [3], 'E')[0] == 0.0
    assert pfqn_manjunath(L_C, N_C, Z_C, z, [-1], 'L')[0] == 0.0
    assert abs(pfqn_manjunath(L_C, N_C, Z_C, _a1[None, :], [-1], 'G')[0]
               - G_ref) / G_ref < 1e-13


def test_peak_is_the_class_lattice_when_unconstrained():
    """With no extra rows the only live axes are the classes, so the realised
    cost must be exactly the lattice pfqn_ca walks."""
    _, _, peak = pfqn_manjunath(L_C, N_C, Z_C)
    assert peak == (N_C[0] + 1) * (N_C[1] + 1)


@pytest.mark.parametrize("A,b,sense", [
    (np.concatenate([[0.5], np.zeros(S_C * 2 - 1)])[None, :], [1], 'L'),   # fractional A
    (np.concatenate([[-1.0], np.zeros(S_C * 2 - 1)])[None, :], [1], 'L'),  # negative A
    (_a1[None, :], [1.5], 'L'),                                            # fractional b
    (_a1[None, :], [1], 'X'),                                              # bad sense
])
def test_refusals(A, b, sense):
    with pytest.raises(ValueError):
        pfqn_manjunath(L_C, N_C, Z_C, A, b, sense)


def test_negative_population_is_empty_not_an_error():
    G, lG, _ = pfqn_manjunath(L_C, [-1, 2], Z_C)
    assert G == 0.0
    assert lG == -np.inf


# ---------------------------------------------------------------------------
# The per-class decomposition (stats=True)
# ---------------------------------------------------------------------------
# Reference instance: PS queue (demands 1, 2) + delay (think 2, 4), N = [4 4]
# both starting at the delay, with 2*n1 + 3*n2 <= 10 on the queue occupancy.
# The expectations are the stationary law of an INDEPENDENTLY built exact CTMC
# under HOLD truncation (a refused admission is a deleted transition), which
# agrees with the truncated product form to 8.2e-17 -- so these are not the
# routine restating itself.
L_S = np.array([[1.0, 2.0]])
Z_S = np.array([[2.0, 4.0]])
N_S = [4, 4]
A_S = np.zeros((1, 4)); A_S[0, 0] = 2; A_S[0, 2] = 3
B_S = [10]

STATS_REF = {
    'Q':       [1.88612099644128, 1.50177935943061],
    'X':       [0.544483985765125, 0.224199288256228],
    'U':       [0.544483985765125, 0.448398576512456],
    'think':   [1.08896797153025, 0.896797153024911],
    'blocked': [1.02491103202847, 1.60142348754448],
    'delay':   [2.11387900355872, 2.49822064056939],
}


@pytest.mark.parametrize("field", sorted(STATS_REF))
def test_stats_match_the_exact_hold_chain(field):
    _, _, _, st = pfqn_manjunath(L_S, N_S, Z_S, A_S, B_S, 'L', stats=True)
    assert np.allclose(getattr(st, field), STATS_REF[field], rtol=1e-12, atol=0)


def test_stats_conserve_the_population():
    """A blocked job never leaves the delay, so nothing escapes the accounting:
    queue + thinking + held must be exactly N, class by class."""
    _, _, _, st = pfqn_manjunath(L_S, N_S, Z_S, A_S, B_S, 'L', stats=True)
    assert np.allclose(st.Q + st.think + st.blocked, N_S, rtol=0, atol=1e-12)
    assert np.allclose(st.delay, st.think + st.blocked, rtol=0, atol=1e-12)


def test_unconstrained_throughput_reduces_to_the_classical_ratio():
    """With no rows the shifted right-hand side is the original one, so the
    throughput identity collapses to the textbook X_r = G(N-e_r)/G(N)."""
    _, _, _, st = pfqn_manjunath(L_S, N_S, Z_S, stats=True)
    _, lG = pfqn_ca(L_S, np.array(N_S, dtype=float), Z_S.sum(axis=0))
    for r in range(2):
        Nr = list(N_S); Nr[r] -= 1
        _, lGr = pfqn_ca(L_S, np.array(Nr, dtype=float), Z_S.sum(axis=0))
        assert abs(st.X[r] - math.exp(lGr - lG)) < 1e-13
    # and nothing can be held when nothing is constrained
    assert np.allclose(st.blocked, 0.0, atol=1e-12)


def test_stats_refuse_two_delay_stations():
    with pytest.raises(ValueError, match="one delay station"):
        pfqn_manjunath(L_S, N_S, np.vstack([Z_S, Z_S]), np.zeros((1, 6)), B_S, 'L', stats=True)


def test_stats_refuse_two_queueing_stations():
    # The delay->q1->q2->delay cycle makes the chain irreversible, so Kelly
    # truncation no longer holds and the product form is not the stationary law.
    with pytest.raises(ValueError, match="one queueing station"):
        pfqn_manjunath(np.vstack([L_S, L_S]), N_S, Z_S, np.zeros((1, 6)), B_S, 'L', stats=True)


def test_stats_refuse_a_delay_inside_the_region():
    A_bad = np.zeros((1, 4)); A_bad[0, 1] = 1        # column 1 is (delay, class 1)
    with pytest.raises(ValueError, match="OUTSIDE"):
        pfqn_manjunath(L_S, N_S, Z_S, A_bad, B_S, 'L', stats=True)


def test_stats_are_opt_in():
    """The three-tuple contract is unchanged for callers that never ask."""
    assert len(pfqn_manjunath(L_S, N_S, Z_S, A_S, B_S, 'L')) == 3
    assert len(pfqn_manjunath(L_S, N_S, Z_S, A_S, B_S, 'L', stats=True)) == 4
