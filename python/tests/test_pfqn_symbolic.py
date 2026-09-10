"""
Symbolic scaling functions in the gld family.

pfqn_gldsingle and pfqn_gld never COMPARE a rate, only divide by one, so they
can take a rate matrix whose entries are sympy expressions and need no
threshold declared anywhere. pfqn_lldsingle and pfqn_lld cannot: their whole
saving rests on locating the index past which a row is constant, and that
comparison has no truth value on a symbol, so they hand symbolic input back to
the gld family.

Every assertion here pins the SUBSTITUTED value against the numeric routine, so
a wrong expression cannot pass by being merely well-formed.
"""

import numpy as np
import pytest

sp = pytest.importorskip("sympy")

from line_solver.api.pfqn import pfqn_gld, pfqn_gldsingle, pfqn_lld, pfqn_lldsingle


def test_symbolic_rates_single_class():
    m1, m2, c = sp.symbols("m1 m2 c", positive=True)
    L = np.array([[2.0], [3.0]], dtype=object)
    mu = np.array([[m1, m2, c], [1, 1, 1]], dtype=object)
    G = sp.simplify(pfqn_gldsingle(L, np.array([3.0]), mu).G)
    assert G.free_symbols == {m1, m2, c}
    got = float(G.subs({m1: 1, m2: 1, c: 1}))
    ref = float(np.exp(pfqn_gldsingle(np.array([[2.0], [3.0]]), np.array([3.0]),
                                      np.ones((2, 3))).lG))
    assert got == pytest.approx(ref, rel=1e-12)


def test_multiserver_written_symbolically_needs_no_threshold():
    """
    The point of the gld family: a multiserver row given as sym(1), sym(2),
    sym(2), sym(2) is handled without anyone identifying where it settles.
    """
    mu = np.array([[sp.Integer(1), sp.Integer(2), sp.Integer(2), sp.Integer(2)],
                   [1, 1, 1, 1]], dtype=object)
    got = float(sp.simplify(pfqn_gldsingle(np.array([[2.0], [3.0]], dtype=object),
                                           np.array([4.0]), mu).G))
    ref = float(np.exp(pfqn_gldsingle(np.array([[2.0], [3.0]]), np.array([4.0]),
                                      np.array([np.minimum(np.arange(1, 5), 2),
                                                np.ones(4)])).lG))
    assert got == pytest.approx(ref, rel=1e-12)


def test_symbolic_demands_and_rates():
    a, b, m1, m2 = sp.symbols("a b m1 m2", positive=True)
    L = np.array([[a], [b]], dtype=object)
    mu = np.array([[m1, m1], [m2, m2]], dtype=object)
    G = sp.simplify(pfqn_gldsingle(L, np.array([2.0]), mu).G)
    assert G.free_symbols == {a, b, m1, m2}
    got = float(G.subs({a: 2, b: 3, m1: 1, m2: 1}))
    ref = float(np.exp(pfqn_gldsingle(np.array([[2.0], [3.0]]), np.array([2.0]),
                                      np.ones((2, 2))).lG))
    assert got == pytest.approx(ref, rel=1e-12)


def test_gld_multiclass_symbolic_rates():
    m1, m2 = sp.symbols("m1 m2", positive=True)
    L = np.array([[2.0, 1.0], [1.0, 3.0]], dtype=object)
    mu = np.array([[m1, m1], [m2, m2]], dtype=object)
    G = sp.simplify(pfqn_gld(L, np.array([1.0, 1.0]), mu).G)
    assert G.free_symbols == {m1, m2}
    got = float(G.subs({m1: 1, m2: 1}))
    ref = pfqn_gld(np.array([[2.0, 1.0], [1.0, 3.0]]), np.array([1.0, 1.0]),
                   np.ones((2, 2))).G
    assert got == pytest.approx(ref, rel=1e-12)


def test_gld_single_station_multinomial_is_exact():
    """
    The M == 1 closed form must build the multinomial as a ratio of factorials.
    Routing it through exp(factln(...)) would leave a float carried as a
    rational approximation, e.g. exp(2473854946935173/2251799813685248) in
    place of the integer 3, which no simplification recovers.
    """
    a, b = sp.symbols("a b", positive=True)
    G = sp.simplify(pfqn_gld(np.array([[a, b]], dtype=object),
                             np.array([1.0, 2.0]),
                             np.ones((1, 3), dtype=object)).G)
    assert G == 3 * a * b ** 2
    assert not G.atoms(sp.exp)


def test_lld_family_delegates_symbolic_to_gld():
    """pfqn_lld cannot find a threshold on a symbol, so it must hand over."""
    m1, m2 = sp.symbols("m1 m2", positive=True)
    L = np.array([[2.0, 1.0], [1.0, 3.0]], dtype=object)
    mu = np.array([[m1, m1], [m2, m2]], dtype=object)
    N = np.array([1.0, 1.0])
    assert sp.simplify(pfqn_lld(L, N, mu).G - pfqn_gld(L, N, mu).G) == 0


def test_numeric_path_untouched_by_the_symbolic_guard():
    """The is_sym guard is False on numbers: nothing may move."""
    rng = np.random.default_rng(3)
    for M in (1, 2, 4):
        for Ntot in (1, 5, 12):
            L = 1 + 9 * rng.random((M, 1))
            mu = np.minimum(np.arange(1, Ntot + 1), 3) * np.ones((M, 1))
            a = pfqn_gldsingle(L, np.array([float(Ntot)]), mu).lG
            b = pfqn_lldsingle(L, np.array([float(Ntot)]), mu).lG
            assert a == b
