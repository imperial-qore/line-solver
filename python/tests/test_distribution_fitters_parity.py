"""Parity guard for the distribution-class moment fitters.

Every value here was produced by the MATLAB reference
(`matlab/src/lang/processes/*.m`) and agrees with native Python to machine
precision; the JAR agrees where it implements the same method. The algebra
behind these fitters is proved in `sage/proofs/distribution_fitters.py`.

Two things this catches that a round trip cannot:
- a fitter that goes MISSING in one codebase (several of these did not exist in
  native Python until 2026-07-22);
- an approximation whose constant drifts, e.g. the Justus exponent in the
  Weibull fit, where the realized SCV is deliberately NOT the requested one.
"""
import warnings

import pytest

from line_solver.distributions.continuous import (Erlang, Exp, Gamma, HyperExp,
                                                  Lognormal, Normal, Pareto,
                                                  Weibull)
from line_solver.distributions.markovian import APH, PH, Cox2, Coxian, MMPP2

warnings.simplefilter('ignore')


def _mean_scv(d):
    m = d.getMean()
    return m, d.getVar() / (m * m)


# (label, factory, expected mean, expected SCV) -- MATLAB reference values
CASES = [
    ('Weibull 0.25', lambda: Weibull.fit_mean_and_scv(1, 0.25), 1.0, 0.24547114143809),
    ('Weibull 4', lambda: Weibull.fit_mean_and_scv(1, 4), 1.0, 5.92941859548084),
    ('Gamma', lambda: Gamma.fit_mean_and_scv(2, 0.5), 2.0, 0.5),
    ('Lognormal', lambda: Lognormal.fit_mean_and_scv(2, 0.5), 2.0, 0.5),
    ('Pareto', lambda: Pareto.fit_mean_and_scv(2, 0.5), 2.0, 0.5),
    ('Normal', lambda: Normal.fit_mean_and_scv(2, 0.25), 2.0, 0.25),
    ('Exp', lambda: Exp.fit_mean_and_scv(2, 1), 2.0, 1.0),
    ('Erlang scv', lambda: Erlang.fit_mean_and_scv(1, 0.25), 1.0, 0.25),
    ('Erlang fit', lambda: Erlang.fit(1, 0.25, 0), 1.0, 0.25),
    ('HyperExp', lambda: HyperExp.fit_mean_and_scv(1, 4), 1.0, 4.0),
    ('HyperExp balanced', lambda: HyperExp.fit_mean_and_scv_balanced(1, 4), 1.0, 4.0),
    ('HyperExp mean', lambda: HyperExp.fit_mean(2), 2.0, 1.0),
    ('HyperExp rate', lambda: HyperExp.fit_rate(0.5), 2.0, 1.0),
    ('Cox2 scv>1', lambda: Cox2.fit_mean_and_scv(1, 4), 1.0, 4.0),
    ('Cox2 scv<1', lambda: Cox2.fit_mean_and_scv(1, 0.7), 1.0, 0.7),
    ('Cox2 central', lambda: Cox2.fit_central(1.9, 3.89, 2.06170806220559), 1.9, 1.07756232687094),
    ('Coxian central', lambda: Coxian.fit_central(1, 0.99, 1.999), 1.0, 0.99),
    ('Coxian central fallback', lambda: Coxian.fit_central(2, 2.0, 3.0), 2.0, 0.5),
    ('APH scv', lambda: APH.fit_mean_and_scv(1, 4), 1.0, 4.0),
    ('APH fit', lambda: APH.fit(1, 4, 3), 1.0, 4.0),
    ('PH scv', lambda: PH.fit_mean_and_scv(1, 4), 1.0, 4.0),
    ('MMPP2 raw+decay', lambda: MMPP2.fitRawMomentsAndACFDecay(1, 5, 45, 0.3), 1.0, 4.0),
    ('MMPP2 raw+lag1', lambda: MMPP2.fitRawMomentsAndACFLag1(1, 5, 45, 0.1125), 1.0, 4.0),
    ('MMPP2 raw+idc', lambda: MMPP2.fitRawMomentsAndIDC(1, 5, 45, 6), 1.0, 4.0),
    ('MMPP2 central+decay', lambda: MMPP2.fitCentralAndACFDecay(1, 4, 3.5, 0.3), 1.0, 4.0),
]


@pytest.mark.parametrize("label,factory,mean,scv", CASES)
def test_matches_the_matlab_reference(label, factory, mean, scv):
    m, s = _mean_scv(factory())
    assert m == pytest.approx(mean, rel=1e-9)
    assert s == pytest.approx(scv, rel=1e-9)


def test_weibull_uses_the_justus_shape():
    # k = CV^(-1.086) with CV = sqrt(scv), and the scale makes the MEAN exact.
    # The SCV is deliberately approximate; pin the exponent through a realized
    # value so a change of constant cannot pass unnoticed.
    import numpy as np
    from scipy.special import gamma as gamma_fn
    w = Weibull.fit_mean_and_scv(1.0, 0.5)
    k = np.sqrt(0.5) ** (-1.086)
    assert w.shape == pytest.approx(k, rel=1e-12)
    assert w.scale == pytest.approx(1.0 / gamma_fn(1 + 1.0 / k), rel=1e-12)
    assert w.getMean() == pytest.approx(1.0, rel=1e-12)


def test_mmpp2_factories_realize_the_requested_autocorrelation():
    from line_solver.api.mam.map_analysis import map_gamma2
    import numpy as np
    d = MMPP2.fitRawMomentsAndACFDecay(1, 5, 45, 0.3)
    D0, D1 = (np.asarray(x, float) for x in d.getRepresentation()[:2])
    assert map_gamma2(D0, D1) == pytest.approx(0.3, rel=1e-8)
    # lag-1 form: rho1 = gamma2*(1-1/SCV)/2 = 0.3*(1-1/4)/2 = 0.1125
    d = MMPP2.fitRawMomentsAndACFLag1(1, 5, 45, 0.1125)
    D0, D1 = (np.asarray(x, float) for x in d.getRepresentation()[:2])
    assert map_gamma2(D0, D1) == pytest.approx(0.3, rel=1e-8)


def test_coxian_fit_central_matches_the_third_moment():
    # Coxian.fitCentral delegates to the exact three-moment Cox2 fit in MATLAB
    # and the JAR, and only falls back to two moments when that solution misses
    # the SCV by more than 1%. Native Python used to discard the skew outright,
    # which silently returned a different arrival process (Cox/M/1 QLen 0.99502
    # instead of 0.99517). Mean and SCV alone cannot catch that, so pin the skew.
    cx = Coxian.fit_central(1, 0.99, 1.999)
    assert cx.getSkew() == pytest.approx(1.999, rel=1e-9)
    ref = Cox2.fit_central(1, 0.99, 1.999)
    assert sorted(cx.means) == pytest.approx(sorted(ref.means), rel=1e-12)
    assert list(cx.probs) == pytest.approx(list(ref.probs), rel=1e-12)
    # Second argument is the VARIANCE, as in MATLAB/JAR fitCentral(MEAN,VAR,SKEW).
    scaled = Coxian.fit_central(2, 4 * 0.99, 1.999)
    assert scaled.getVar() / 4.0 == pytest.approx(0.99, rel=1e-9)
