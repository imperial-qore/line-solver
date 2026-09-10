"""Gerasimov's residue closed form for the normalizing constant, generalized to R classes.

A. I. Gerasimov, "On Normalizing Constants in Multiclass Queueing Networks",
Operations Research 43(4):704-711, 1995.

The paper's own worked example (Figure 1, Table I, and the two closed forms in the
Appendix) pins the R = 2 algorithm rather than this port. One superscript is lost in
the scan of the second Appendix form -- its leading term prints as x21/(x21-x11)^3
where the series it sums requires x21^(N1+4)/(x21-x11)^3, which is what is asserted
here; with that reading it agrees with pfqn_ca to 1e-14, and with the printed reading
it is off by a factor of order one.

Everything else is checked against pfqn_ca, which is exact and independent.
"""

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_ca, pfqn_gerasimov

Y = 0.05202          # Table I: x12 = x32
X11 = 0.06627        # Table I: x11


def relerr(a, b):
    return abs(a - b) / abs(b)


def test_appendix_form_i_equal_loads():
    """Appendix (i), x11 = x21 = c: G = (y^2/6) c^N1 (N1^3+9N1^2+26N1+18)."""
    c = 0.08
    for n1 in range(1, 7):
        L = np.array([[c, Y], [c, 0.0], [0.0, Y]])
        G, _ = pfqn_gerasimov(L, [n1, 2])
        Gf = (Y ** 2 / 6) * c ** n1 * (n1 ** 3 + 9 * n1 ** 2 + 26 * n1 + 18)
        assert relerr(G, Gf) < 1e-12


def test_appendix_form_ii_distinct_loads():
    """Appendix (ii), x11 != x21, with the leading superscript restored."""
    x21 = 5.0
    for n1 in range(1, 7):
        L = np.array([[X11, Y], [x21, 0.0], [0.0, Y]])
        G, _ = pfqn_gerasimov(L, [n1, 2])
        Gf = (Y ** 2 / X11) * (
            x21 ** (n1 + 4) / (x21 - X11) ** 3
            + (n1 ** 2 + 7 * n1 + 12) * X11 ** (n1 + 2) / (2 * (X11 - x21))
            - (n1 + 4) * X11 ** (n1 + 3) / (X11 - x21) ** 2
            + X11 ** (n1 + 4) / (x21 - X11) ** 3
            - x21 ** (n1 + 1))
        assert relerr(G, Gf) < 1e-7


CASES = [
    # (label, L, N, Z)
    ("fig1", [[X11, Y], [5.0, 0.0], [0.0, Y]], [3, 2], [0, 0]),
    ("coincident poles x11=x21", [[0.08, Y], [0.08, 0.0], [0.0, Y]], [4, 2], [0, 0]),
    ("tied class-2 demands", [[1, 2], [3, 2], [5, 7]], [3, 2], [0, 0]),
    ("identical station rows", [[1, 2], [1, 2], [5, 7]], [3, 3], [0, 0]),
    ("two stations unvisited by class 2", [[1, 0], [3, 0], [5, 7]], [4, 2], [0, 0]),
    ("think time", [[1, 2], [3, 4]], [3, 2], [0.5, 1.5]),
    ("think time one class", [[1, 2], [3, 4], [2, 1]], [2, 3], [1.0, 0.0]),
    ("single class", [[1], [2], [3]], [6], [0]),
    ("single class with delay", [[1], [2], [3]], [6], [2]),
    ("empty class", [[1, 2], [3, 4]], [4, 0], [0, 0]),
    ("three classes", [[1, 2, 3], [3, 4, 1], [2, 1, 2]], [2, 2, 2], [0, 0, 0]),
    ("three classes with delay", [[1, 2, 3], [3, 4, 1], [2, 1, 2]], [3, 2, 1], [0.7, 0, 0.3]),
    ("four classes", [[1, 2, 3, 1], [3, 4, 1, 2], [2, 1, 2, 3]], [2, 1, 2, 1], [0, 0, 0, 0]),
    ("four classes with delay", [[1, 2, 3, 1], [3, 4, 1, 2], [2, 1, 2, 3]], [2, 1, 2, 1],
     [0.3, 0.2, 0, 0.1]),
]


@pytest.mark.parametrize("label,L,N,Z", CASES, ids=[c[0] for c in CASES])
def test_matches_convolution(label, L, N, Z):
    """The degeneracies the paper's hypotheses exclude are ordinary cases here."""
    L = np.array(L, dtype=float)
    N = np.array(N, dtype=float)
    Z = np.array(Z, dtype=float)
    Gg, _ = pfqn_gerasimov(L, N, Z)
    Gc, _ = pfqn_ca(L, N, Z)
    assert relerr(Gg, Gc) < 1e-12


def test_eliminated_population_is_free():
    """A population removed by residues enters only as a pole ORDER, so it costs
    nothing: the answer must stay exact as it grows by three orders of magnitude."""
    L = np.array([[1.0, 2], [3, 1], [2, 4], [0.5, 0.7]])
    for n2 in (10, 100, 2000, 20000):
        _, lGg = pfqn_gerasimov(L, [6, n2])
        _, lGc = pfqn_ca(L, np.array([6.0, n2]), np.zeros(2))
        assert relerr(lGg, lGc) < 1e-12


def test_maxterms_refuses_rather_than_truncates():
    """A truncated residue sum is not a bound on G, it is a wrong number."""
    L = np.round(np.random.RandomState(0).rand(6, 3) * 9 + 1) / 3
    with pytest.raises(ValueError, match="maxterms"):
        pfqn_gerasimov(L, [4, 8, 8], maxterms=50)


def test_rejects_noninteger_population():
    with pytest.raises(ValueError, match="integer"):
        pfqn_gerasimov(np.array([[1.0, 2.0]]), [1.5, 2])
