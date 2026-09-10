"""
Tests for the Ott-Yashkov sojourn time transform of the M/G/1-PS queue.

The references are exact where one exists: the conditional mean is x/(1-rho) for
any service law, the M/M/1-PS unconditional moments are the Coffman-Muntz-Trotter
values that ``qsys_mm1_ps`` returns, and the k = 0 atom is (1-rho)*exp(-lam*x).
The conditional CDF values were additionally checked against 40000 replications of
an exact event-driven processor-sharing simulation, and the unconditional CDF
against SolverLDES.
"""

import numpy as np
import pytest

from line_solver.api.qsys import qsys_mg1_ps, qsys_mm1_ps

LAM = 0.7
MU = 1.0
EXP_ALPHA = [1.0]
EXP_T = [[-1.0]]


def test_conditional_mean_is_insensitive():
    x = [0.25, 1.0, 4.0]
    r = qsys_mg1_ps(LAM, EXP_ALPHA, EXP_T, x=x)
    assert np.allclose(r["meanCond"], np.array(x) / (1 - r["rho"]), atol=1e-12)
    assert abs(r["meanUncond"] - 1.0 / (MU * (1 - LAM / MU))) < 1e-12


def test_mm1ps_moments_match_coffman_muntz_trotter():
    r = qsys_mg1_ps(LAM, EXP_ALPHA, EXP_T)
    rho = LAM / MU
    m2exact = 4.0 / (MU ** 2 * (1 - rho) ** 2 * (2 - rho))
    assert abs(r["m2Uncond"] - m2exact) < 1e-6 * m2exact
    W, W2, _ = qsys_mm1_ps([LAM], [MU])
    assert abs(r["meanUncond"] - W[0]) < 1e-12
    assert abs(r["m2Uncond"] - W2[0]) < 1e-6 * W2[0]


def test_transform_at_the_origin_and_at_large_arguments():
    r = qsys_mg1_ps(LAM, EXP_ALPHA, EXP_T)
    assert abs(np.real(r["lstCond"](0.0, 2.0)) - 1.0) < 1e-12
    atom = (1 - r["rho"]) * np.exp(-LAM * 2.0)
    assert abs(np.real(r["lstExcess"](1e8, 2.0)) - atom) < 1e-7


def test_conditional_distribution_against_simulation():
    t = [1.2, 1.5, 1.8, 2.4, 3.0, 4.0, 6.0, 10.0, 16.0]
    sim = [0.207350, 0.282000, 0.342050, 0.481000, 0.584625,
           0.715825, 0.866050, 0.969875, 0.996150]
    r = qsys_mg1_ps(LAM, EXP_ALPHA, EXP_T, x=[1.0], t=t)
    cdf = r["cdfCond"][0]
    assert np.allclose(cdf, sim, atol=5e-3)
    assert np.all(np.diff(cdf) >= 0)
    assert abs(r["atomCond"][0] - (1 - r["rho"]) * np.exp(-LAM)) < 1e-12


def test_no_mass_below_the_service_requirement():
    r = qsys_mg1_ps(LAM, EXP_ALPHA, EXP_T, x=[2.0], t=[0.5, 1.9, 2.0])
    assert r["cdfCond"][0, 0] == 0.0
    assert r["cdfCond"][0, 1] == 0.0
    assert abs(r["cdfCond"][0, 2] - r["atomCond"][0]) < 1e-12


def test_density_is_undefined_on_the_lattice():
    r = qsys_mg1_ps(LAM, EXP_ALPHA, EXP_T, x=[1.0], t=[1.5, 2.0, 3.0])
    assert np.isfinite(r["pdfCond"][0, 0])
    assert np.isnan(r["pdfCond"][0, 1])
    assert np.isnan(r["pdfCond"][0, 2])


def test_phase_type_agrees_with_the_transform_handle():
    mu, lam = 4.0, 1.2
    x, s = [0.25, 1.0, 3.0], [0.2, 2.0]
    ph = qsys_mg1_ps(lam, [1.0, 0.0], [[-mu, mu], [0.0, -mu]], x=x, s=s)
    hd = qsys_mg1_ps(lam, lambda tau: (mu / (mu + tau)) ** 2, 2.0 / mu, x=x, s=s,
                     pdf=lambda y: mu ** 2 * y * np.exp(-mu * y))
    assert np.allclose(ph["lstCondVal"], hd["lstCondVal"], atol=1e-8)
    assert np.allclose(ph["lstUncondVal"], hd["lstUncondVal"], atol=1e-8)
    assert abs(ph["m2Uncond"] - hd["m2Uncond"]) < 1e-4 * ph["m2Uncond"]


def test_unconditional_distribution_is_a_proper_law():
    mu, lam = 4.0, 1.2
    t = [0.1, 0.5, 1.0, 2.0, 5.0, 20.0]
    r = qsys_mg1_ps(lam, [1.0, 0.0], [[-mu, mu], [0.0, -mu]], t=t)
    cdf = r["cdfUncond"]
    assert np.all(np.diff(cdf) >= 0)
    assert 0 < cdf[0] and cdf[-1] > 0.999
    assert np.all(r["pdfUncond"] >= 0)
    # the mean recovered from the density matches the exact m1/(1-rho)
    assert abs(r["meanUncond"] - 0.5 / (1 - lam * 0.5)) < 1e-12


def test_repeated_pole_falls_back_to_the_companion_matrix():
    """At s = -(sqrt(mu)-sqrt(lam))**2 the two poles of P coincide, so the residue
    formula is unusable and the companion-matrix branch must take over. There the
    confluent partial fraction is exact: Ahat/(tau-r)**2 inverts to
    (a1 + (a1*r + a0)*x)*exp(r*x)."""
    rho = LAM / MU
    sstar = -(np.sqrt(MU) - np.sqrt(LAM)) ** 2
    assert np.ptp(np.roots([1.0, MU - sstar - LAM, -sstar * MU])) == 0.0, "the poles coincide"
    x = np.array([0.25, 1.0, 3.0, 10.0])
    r = qsys_mg1_ps(LAM, EXP_ALPHA, EXP_T, x=list(x), s=[sstar])
    a1 = 1 - rho
    a0 = (1 - rho) * (MU - LAM) + sstar * rho
    rr = -(MU - sstar - LAM) / 2
    exact = (1 - rho) / ((a1 + (a1 * rr + a0) * x) * np.exp(rr * x))
    assert np.allclose(r["lstCondVal"][:, 0], exact, rtol=1e-12)
    # and the residue branch either side of the confluence brackets it
    lo = np.real(r["lstCond"](sstar + 1e-6, 1.0))
    hi = np.real(r["lstCond"](sstar - 1e-6, 1.0))
    assert lo < r["lstCondVal"][1, 0] < hi


def test_malformed_inputs_are_rejected():
    with pytest.raises(ValueError):
        qsys_mg1_ps(2.5, EXP_ALPHA, EXP_T)
    with pytest.raises(ValueError):
        qsys_mg1_ps(LAM, [0.5], EXP_T)
    with pytest.raises(ValueError):
        qsys_mg1_ps(LAM, EXP_ALPHA, EXP_T, nterms=40)
    with pytest.raises(ValueError):
        qsys_mg1_ps(LAM, lambda tau: 1.0 / (1.0 + tau), 1.0, t=[1.0])
