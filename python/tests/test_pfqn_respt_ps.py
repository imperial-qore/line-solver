"""Tests for the processor-sharing sojourn-time moments of Mitra and Morrison
(1983): qsys_mm1_ps for the open station and pfqn_respt_ps_moments for the
closed terminal-driven system.

Mirrors jar/src/test/java/jline/api/PfqnResptPsTest.java and
line-test.git/test/testsAPI/test_pfqn_respt_ps.m one for one.

The mathematics is anchored to results that do not come from that paper: the
single-class M/M/1-PS second moment of Coffman, Muntz and Trotter (1970), and
the trivial one-job closed system whose sojourn time is the service time itself.
The paper's own two routes are then cross-checked against each other, the exact
route being the reference and the asymptotic one required to converge to it as
the expansion parameter grows. The literal values are the MATLAB and JAR
outputs, so a divergence in any codebase shows up here.
"""

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_respt_ps_moments
from line_solver.api.qsys import qsys_mm1_ps


def test_open_single_class_cmt():
    """The single-class M/M/1-PS second moment is Coffman, Muntz and Trotter's
    4/(mu^2 (1-rho)^2 (2-rho))."""
    mu = 2.0
    for rho in [0.2, 0.5, 0.8]:
        W, W2, alpha = qsys_mm1_ps([rho * mu], [mu])
        assert W[0] == pytest.approx(1.0 / (mu * (1 - rho)), rel=1e-12)
        cmt = 4.0 / (mu ** 2 * (1 - rho) ** 2 * (2 - rho))
        assert W2[0] == pytest.approx(cmt, rel=1e-12)
        assert alpha == pytest.approx(1 - rho, rel=1e-12)


def test_open_multiclass_parity():
    """Cross-codebase values: MATLAB qsys_mm1_ps and JAR Qsys_mm1_ps agree."""
    W, W2, alpha = qsys_mm1_ps([0.3, 0.4], [1.0, 3.0])
    assert W[0] == pytest.approx(1.764705882353, rel=1e-10)
    assert W[1] == pytest.approx(0.588235294118, rel=1e-10)
    assert W2[0] == pytest.approx(7.750865051903, rel=1e-10)
    assert W2[1] == pytest.approx(0.927201263144, rel=1e-10)
    assert alpha == pytest.approx(1 - 0.3 - 0.4 / 3.0, rel=1e-12)


def test_open_unstable_raises():
    """An unstable station has no moments; a negative variance would be worse
    than an error."""
    with pytest.raises(ValueError):
        qsys_mm1_ps([1.0, 1.0], [1.0, 2.0])


def test_closed_single_job():
    """A lone job is never delayed, so its sojourn time is its exponential
    service time: W = S and E[W^2] = 2 S^2."""
    S = 0.5
    W, W2, out = pfqn_respt_ps_moments([S], [1], [5.0], 'exact')
    assert W[0] == pytest.approx(S, rel=1e-12)
    assert W2[0] == pytest.approx(2 * S ** 2, rel=1e-12)


def test_closed_exact_parity():
    W, W2, out = pfqn_respt_ps_moments([1.0, 0.5], [4, 2], [50.0, 100.0], 'exact')
    assert W[0] == pytest.approx(1.072358918582, rel=1e-9)
    assert W[1] == pytest.approx(0.544470487795, rel=1e-9)
    assert W2[0] == pytest.approx(2.372013893895, rel=1e-9)
    assert W2[1] == pytest.approx(0.623810622242, rel=1e-9)
    assert out.method == ['exact', 'exact']
    assert list(out.nstates) == [12, 10]


def test_closed_asymptotic_parity():
    W, W2, out = pfqn_respt_ps_moments([1.0, 0.5], [4, 2], [50.0, 100.0],
                                       'asymptotic')
    assert W[0] == pytest.approx(1.072160744545, rel=1e-9)
    assert W2[0] == pytest.approx(2.370214800752, rel=1e-9)
    assert W2[1] == pytest.approx(0.623075482738, rel=1e-9)
    assert out.c0[0] == pytest.approx(2.392420108971, rel=1e-9)
    assert out.c1[0] == pytest.approx(-4.441061643908, rel=1e-9)
    assert out.expansionParam == pytest.approx(200.0, rel=1e-12)


def test_asymptotic_converges():
    """Two terms of an expansion in 1/Nexp must get better as Nexp grows, and
    the exact route is what they must approach."""
    q = [1.0, 2.0]
    K = [3, 2]
    prev = np.inf
    for scale in [1, 2, 4, 8]:
        Z = [50.0 * scale, 100.0 * scale]
        S = [1 / q[0], 1 / q[1]]
        N = [K[0] + 1, K[1]]
        _, W2e, _ = pfqn_respt_ps_moments(S, N, Z, 'exact')
        _, W2a, _ = pfqn_respt_ps_moments(S, N, Z, 'asymptotic')
        err = abs(W2a[0] / W2e[0] - 1)
        assert err < prev
        prev = err
    assert prev < 1e-5


def test_closed_leading_term_is_open():
    """With lambda_j = K_j/Z_j held fixed while K_j grows, the closed system
    tends to the open one; c0 is already exactly the open second moment."""
    q = [1.0, 3.0]
    lam = [0.3, 0.4]
    _, W2o, _ = qsys_mm1_ps(lam, q)
    for K in [10, 100, 1000]:
        S = [1 / q[0], 1 / q[1]]
        N = [K + 1, K]
        Z = [K / lam[0], K / lam[1]]
        _, _, out = pfqn_respt_ps_moments(S, N, Z, 'asymptotic')
        assert out.c0[0] == pytest.approx(W2o[0], rel=1e-10)


def test_auto_routes():
    """auto is exact on a small state space and asymptotic on a large one."""
    _, _, small = pfqn_respt_ps_moments([1.0], [5], [50.0])
    assert small.method[0] == 'exact'
    _, _, big = pfqn_respt_ps_moments([1.0], [100000], [5.0e6])
    assert big.method[0] == 'asymptotic'


def test_empty_class_is_nan():
    """An unpopulated class has no sojourn time, which is a blank and not a
    zero."""
    W, W2, out = pfqn_respt_ps_moments([1.0, 0.5], [3, 0], [50.0, 100.0])
    assert np.isnan(W[1])
    assert not np.isnan(W[0])
    assert out.method[1] == 'none'


def test_matches_independent_simulation():
    """A sample-path check that owes nothing to the paper: simulate the closed
    terminal-driven system directly and compare the first two sojourn-time
    moments. The tolerance is the sampling error of the run length below."""
    rng = np.random.default_rng(20260724)
    p = [0.02, 0.01]
    q = [1.0, 2.0]
    pop = [4, 2]
    R = 2
    in_cpu = [[] for _ in range(R)]
    n_think = list(pop)
    t = 0.0
    s1 = np.zeros(R)
    s2 = np.zeros(R)
    cnt = np.zeros(R)
    horizon = 4.0e5
    while t < horizon:
        n = sum(len(x) for x in in_cpu)
        rates = [p[j] * n_think[j] for j in range(R)]
        rates += [q[j] * len(in_cpu[j]) / n if n > 0 else 0.0 for j in range(R)]
        tot = sum(rates)
        t += rng.exponential(1.0 / tot)
        u = rng.random() * tot
        acc = 0.0
        for k, rt in enumerate(rates):
            acc += rt
            if u <= acc:
                break
        if k < R:
            n_think[k] -= 1
            in_cpu[k].append(t)
        else:
            j = k - R
            i = rng.integers(len(in_cpu[j]))    # exchangeable within a class
            a = in_cpu[j].pop(i)
            n_think[j] += 1
            if t > 0.1 * horizon:
                w = t - a
                s1[j] += w
                s2[j] += w * w
                cnt[j] += 1
    W, W2, _ = pfqn_respt_ps_moments([1 / q[0], 1 / q[1]], pop,
                                     [1 / p[0], 1 / p[1]], 'exact')
    for j in range(R):
        assert s1[j] / cnt[j] == pytest.approx(W[j], rel=0.03)
        assert s2[j] / cnt[j] == pytest.approx(W2[j], rel=0.05)
