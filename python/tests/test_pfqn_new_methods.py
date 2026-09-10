"""Numerical validation of pfqn_rgf, pfqn_dnc, pfqn_nintmva, pfqn_tay and pfqn_hst
against their source papers and against the exact routes already in the API."""

import numpy as np
import pytest

from line_solver.api.pfqn import (pfqn_ca, pfqn_dnc, pfqn_hst, pfqn_mci, pfqn_mva,
                                  pfqn_nintmva, pfqn_rgf, pfqn_rgfmc,
                                  pfqn_tay)

DOWDY = np.array([1 / 20.3, 0.6114 / 10.1, 0.3886 / 1.2])


def _lg_ca(L, n, z):
    return pfqn_ca(np.asarray(L).reshape(-1, 1), np.array([float(n)]), np.array([float(z)]))[1]


def test_rgf_matches_convolution():
    cases = [[0.5, 0.3, 0.2], [1.0, 1.0, 1.0, 0.4], [0.05, 0.05, 0.05, 0.05, 0.9]]
    for L in cases:
        for n in range(1, 21):
            for z in (0.0, 3.0, 12.5):
                lg_ca = _lg_ca(L, n, z)
                lg_rgf = pfqn_rgf(L, n, z)[1]
                assert lg_rgf == pytest.approx(lg_ca, rel=1e-9, abs=1e-9)


def test_rgf_reproduces_coury_harrison_tables():
    # Coury-Harrison (1997) Sec. 4: IS terminals of load 1, three groups of m devices
    # with per-device loads 0.001/0.002/0.003, one CPU of load 0.004. The m=20, k=300
    # entry is printed 227.52 in both tables, a dropped-zero typo for 227.052.
    p = [0.001, 0.002, 0.003]
    pub = {2: [19.67, 97.83, 192.44, 249.67],
           10: [18.75, 92.33, 178.49, 242.84],
           20: [17.71, 86.444, 164.93, 227.052]}
    for m, expected in pub.items():
        L = [p[0]] * m + [p[1]] * m + [p[2]] * m + [0.004]
        for c, k in enumerate([20, 100, 200, 300]):
            T = np.exp(pfqn_rgf(L, k - 1, 1.0)[1] - pfqn_rgf(L, k, 1.0)[1])
            assert T == pytest.approx(expected[c], rel=5e-4)


def test_rgfmc_matches_convolution_with_think_times():
    # Multiclass RGF: iterated residues (Harrison-Coury 2002 Thm 1) with the
    # delay carried by the Bertozzi-McKenna truncation, which neither RGF paper
    # has.  Exact wherever the alternating residue sum keeps its significance.
    cases = [(np.array([[0.9, 0.4], [0.6, 0.7], [0.3, 1.1]]), [3, 4]),
             (np.array([[0.5, 0.5], [0.2, 0.9]]), [2, 5]),
             (np.array([[0.2, 0.4, 1.2], [0.7, 1.3, 0.3], [1.4, 0.2, 0.7]]),
              [2, 3, 2])]
    for L, N in cases:
        R = L.shape[1]
        for Z in ([0.0] * R, [2.0] + [0.0] * (R - 1), [1.3] * R):
            N = np.asarray(N, dtype=float)
            Z = np.asarray(Z, dtype=float)
            lg_ca = pfqn_ca(L, N, Z)[1]
            assert pfqn_rgfmc(L, N, Z)[1] == pytest.approx(lg_ca, abs=1e-7)


def test_rgfmc_is_insensitive_to_the_eliminated_population():
    # Harrison-Lee sec. 4: a class removed by residues enters only as a pole
    # ORDER, so its population is nearly free -- but only while it carries no
    # think time, since the truncation reintroduces exactly that dependence.
    L = np.array([[0.9, 0.4], [0.6, 0.7], [0.3, 1.1], [0.5, 0.2]])
    for n2 in (20, 200, 2000):
        N = np.array([6.0, float(n2)])
        Z = np.array([1.5, 0.0])
        assert pfqn_rgfmc(L, N, Z)[1] == pytest.approx(pfqn_ca(L, N, Z)[1],
                                                       abs=1e-7)


def test_rgfmc_is_exact_or_refuses_never_wrong():
    # The residue elimination is exact in exact arithmetic, so the only failure
    # mode is cancellation.  Over a random sweep every model must either come
    # back at the convolution's value or raise; a wrong number is the one
    # outcome the cancellation guard exists to make impossible.
    rng = np.random.default_rng(11)
    answered = refused = 0
    for _ in range(120):
        M = int(rng.integers(2, 5))
        R = int(rng.integers(2, 4))
        L = rng.uniform(0.1, 1.5, (M, R))
        N = rng.integers(1, 6, R).astype(float)
        Z = rng.choice([0.0, 0.0, 1.3, 2.5], R)
        try:
            lg = pfqn_rgfmc(L, N, Z)[1]
        except ValueError:
            refused += 1
            continue
        answered += 1
        assert lg == pytest.approx(pfqn_ca(L, N, Z)[1], abs=1e-7)
    assert answered > refused          # the guard must not swallow the method


def test_rgfmc_refuses_rather_than_returning_a_wrong_constant():
    # Near-coincident loads over the eliminated class make the alternating sum
    # cancel; the routine must refuse, not answer.  Same model at a wide
    # maxcancel answers and is wrong, which is what the guard exists to stop.
    L = np.array([[0.81733524, 1.24343101, 0.86870538],
                  [1.4732791, 0.38631325, 0.87522251]])
    N = np.array([2.0, 3.0, 2.0])
    Z = np.zeros(3)
    exact = pfqn_ca(L, N, Z)[1]
    with pytest.raises(ValueError):
        pfqn_rgfmc(L, N, Z, maxcancel=1.0)
    loose = pfqn_rgfmc(L, N, Z, maxcancel=1e9)[1]
    assert abs(loose - exact) > 1e-9


def test_rgfmc_single_class_delegates_to_rgf():
    L = np.array([[0.5], [0.3], [0.9]])
    assert pfqn_rgfmc(L, [4.0], [2.0])[1] == pytest.approx(
        pfqn_rgf([0.5, 0.3, 0.9], 4, 2.0)[1], rel=1e-12)


def test_dnc_is_exact_at_integer_populations():
    for L in (DOWDY, np.array([0.5, 0.5, 0.5, 0.2, 0.2, 0.9])):
        for n in range(1, 13):
            assert pfqn_dnc(L, n)[2] == pytest.approx(_lg_ca(L, n, 0.0), rel=1e-8, abs=1e-8)


def test_dnc_and_nintmva_reproduce_dowdy_gordon_figure5():
    # Dowdy-Gordon (1984) Sec. 5 at DMPavg = 2.5: DNC 2.999, aMVA 2.991.
    assert pfqn_dnc(DOWDY, 2.5)[0] == pytest.approx(2.999, abs=1e-3)
    assert pfqn_nintmva(DOWDY, 2.5, 0.0)[0] == pytest.approx(2.991, abs=1e-3)
    xs = [pfqn_dnc(DOWDY, n)[0] for n in np.arange(2.0, 3.001, 0.25)]
    assert all(b > a for a, b in zip(xs, xs[1:]))


def test_nintmva_matches_exact_mva_at_integer_populations():
    for z in (0.0, 5.0):
        for n in range(1, 13):
            res = pfqn_mva(DOWDY.reshape(-1, 1), np.array([float(n)]), np.array([z]))
            X, Q = np.ravel(res[0])[0], np.ravel(res[2])   # (XN, CN, QN, UN, RN, TN, AN)
            Xi, Qi, _, _ = pfqn_nintmva(DOWDY, n, z)
            assert Xi == pytest.approx(X, rel=1e-12)
            assert np.allclose(Qi, Q, atol=1e-12)


def test_tay_reproduces_survey_example4():
    # Tay's Example 4 (CN82 Fig. 2 network) at N = (8,1); the survey's Tay row is
    # X1 = 0.654, X2 = 0.336, L11 = 5.234, L12 = 2.766, L22 = 1.000.
    X, Q, U, R, it, Qarr = pfqn_tay(np.array([[1.0, 1.0]]), np.array([8.0, 1.0]),
                                    np.array([8.0, 0.0]))
    assert X[0] == pytest.approx(0.654, abs=1e-3)
    assert X[1] == pytest.approx(0.336, abs=1e-3)
    assert X[0] * 8 == pytest.approx(5.234, abs=1e-3)
    assert Q[0, 0] == pytest.approx(2.766, abs=1e-3)
    assert Q[0, 1] == pytest.approx(1.000, abs=1e-3)
    # Table B holds the arrival-instant queue lengths, not re-solves at N - e_r.
    assert Qarr[0, 0, 0] == pytest.approx(2.23, abs=5e-3)
    assert Qarr[0, 0, 1] == pytest.approx(1.98, abs=5e-3)
    assert Qarr[0, 1, 0] == pytest.approx(1.00, abs=5e-3)


def test_tay_handles_empty_classes():
    L = np.array([[1.0, 1.0], [2.0, 1.0]])
    X, Q, U, R, it, _ = pfqn_tay(L, np.array([3.0, 0.0]), np.array([0.0, 0.0]))
    assert X[1] == 0.0
    assert X[0] > 0
    assert np.allclose(Q[:, 1], 0.0)


def test_tay_beats_bard_schweitzer_on_random_models():
    from line_solver.api.pfqn import pfqn_bs
    rng = np.random.default_rng(3)
    worst_tay = worst_bs = 0.0
    for _ in range(40):
        M, R = int(rng.integers(2, 5)), int(rng.integers(2, 4))
        D = 0.2 + rng.random((M, R))
        N = rng.integers(1, 7, R).astype(float)
        Z = rng.random(R) * 3 if rng.random() < 0.5 else np.zeros(R)
        Xe = np.ravel(pfqn_mva(D, N, Z)[0])
        Xt = np.ravel(pfqn_tay(D, N, Z)[0])
        Xb = np.ravel(pfqn_bs(D, N, Z)[0])
        assert np.all(np.isfinite(Xt)) and np.all(Xt > 0)
        worst_tay = max(worst_tay, float(np.max(np.abs(Xt - Xe) / Xe)))
        worst_bs = max(worst_bs, float(np.max(np.abs(Xb - Xe) / Xe)))
    assert worst_tay < worst_bs


def test_hst_reproduces_suri_figure1():
    # Suri (1983) Figure 1: M = 4, N = 7, workloads (1,1,1,2), station 4 analysed.
    L = np.array([1.0, 1.0, 1.0, 2.0])
    s = pfqn_hst(L, 7, 0.0, 3)
    assert s['X'] == pytest.approx(0.48, abs=5e-3)
    assert s['U'] == pytest.approx(0.96, abs=5e-3)
    assert np.allclose(s['p'], [0.037, 0.058, 0.087, 0.124, 0.166, 0.198, 0.198, 0.132],
                       atol=1e-3)
    assert np.allclose(np.abs(s['c']), [0.023, 0.055, 0.097, 0.145, 0.186, 0.193, 0.132],
                       atol=1e-3)
    assert s['total'] == pytest.approx(0.831, abs=1e-3)
    assert s['worst'] == pytest.approx(0.102, abs=1.5e-3)
    # Lemma 3.1: the total equals Q_i(N) - Q_i(N-1).
    assert s['total'] == pytest.approx(s['Q'] - pfqn_hst(L, 6, 0.0, 3)['Q'], abs=1e-10)
    assert pfqn_hst(L, 7, 0.0, 0)['total'] == pytest.approx(0.057, abs=1.5e-3)


def test_mci_lhs_variant_reduces_variance():
    # The Latin-hypercube sampler is the variance-reduced MCI variant; the antithetic
    # one is not (see the pfqn_mci docstring), so only 'lhsmci' is asserted on.
    D = np.array([[0.4, 0.2], [0.3, 0.5], [0.1, 0.6]])
    N = np.array([6.0, 4.0])
    Z = np.array([1.0, 2.0])
    lGex = pfqn_ca(D, N, Z)[1]
    sd = {}
    for variant in ('imci', 'amci', 'lhsmci'):
        errs = []
        for s in range(30):
            np.random.seed(5000 + s)
            errs.append(pfqn_mci(D, N, Z, 20000, variant)[1] - lGex)
        assert abs(float(np.mean(errs))) < 5e-3        # all three stay unbiased
        sd[variant] = float(np.std(errs, ddof=1))
    assert sd['lhsmci'] < sd['imci']
