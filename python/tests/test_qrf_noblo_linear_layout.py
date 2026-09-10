"""Regression test for the QRF no-blocking linear-constraint layout fix.

The decision vector consumed by sub_qrfvar uses a COMPACT per-station-phase
layout (phase index ranges over 1..K[i], not 1..max(K)). An earlier revision of
qrf_noblo_mmi_linear._deltap2/_deltae assumed a max(K)-padded layout; the two
agree only when all K are equal, so for mixed-phase models the linear
constraint matrices indexed the wrong columns and silently disagreed with the
callback constraints (sub_qrfcon_noblo). This test pins the fix by asserting
the linear builder reproduces the callback constraints exactly on a mixed-phase
model (Erlang-2 station beside an exponential station).
"""

import numpy as np

from line_solver.api.mapqn.qrf_noblo_common import (
    sub_qrfcon_noblo, build_q_ld)
from line_solver.api.mapqn.qrf_noblo_mmi_linear import (
    _build_linear_constraints, _compact_num_vars, _deltap2)


def _model_C():
    # station 0: Erlang-2 (mean 1), station 1: exp(1), cyclic routing, N=2
    M, N, MR = 2, 2, 1
    K = np.array([2, 1])
    mu = np.zeros((M, 2, 2)); v = np.zeros((M, 2, 2))
    mu[0, 1, 0] = 2.0; v[0, 0, 1] = 2.0; mu[1, 0, 0] = 1.0
    rt = np.array([[0.0, 1.0], [1.0, 0.0]])
    BB = np.zeros((1, M)); F = np.full(M, N, dtype=int)
    # 5D q: the linear builder mirrors MATLAB qrf_noblo_mmi_linear.m, which
    # always carries the emitting station's population index (and alpha) in q.
    q = build_q_ld(M, K, mu, v, rt, N)
    return M, N, MR, K, BB, F, q


def test_deltap2_is_sequential_bijection_mixed_K():
    M, N, MR, K, _, _, _ = _model_C()
    ctr = 0
    for j in range(M):
        for nj in range(N + 1):
            for k in range(K[j]):
                for i in range(M):
                    for ni in range(N + 1):
                        for h in range(K[i]):
                            for m in range(MR):
                                assert _deltap2(j, nj, k, i, ni, h, m,
                                                M, N, K, MR) == ctr
                                ctr += 1
    # p2 count == (N+1)^2 * MR * SK^2 with SK = sum K = 3 -> 9*9 = 81
    assert ctr == (N + 1) ** 2 * MR * int(K.sum()) ** 2 == 81


def test_linear_builder_matches_callback_mixed_K():
    M, N, MR, K, BB, F, q = _model_C()
    nv = _compact_num_vars(K, N, MR)
    Aeq, beq, Aub, bub = _build_linear_constraints(q, M, MR, BB, F, N, K, nv)
    Aeq = Aeq.toarray(); beq = np.asarray(beq, float)
    has_ub = Aub.shape[0] > 0
    if has_ub:
        Aub = Aub.toarray(); bub = np.asarray(bub, float)
    rng = np.random.RandomState(1)
    for _ in range(25):
        x = rng.rand(nv)
        c, ceq = sub_qrfcon_noblo(x, q, M, MR, BB, F, N, K)
        # builder eq: Aeq x - beq must equal callback ceq
        assert np.max(np.abs((Aeq @ x - beq) - ceq)) < 1e-9
        if has_ub:
            # callback c <= 0  ==  Aub x - bub
            assert np.max(np.abs((Aub @ x - bub) - c)) < 1e-9
