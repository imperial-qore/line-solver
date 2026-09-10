"""Conformance of the ported RODAS against the unmodified Hairer-Wanner Fortran.

THE ORACLE IS THE ORIGINAL, NOT A CLOSED FORM. rodas.py is a transliteration,
so the property that matters is not that it solves these problems well but that
it solves them IDENTICALLY to the Fortran it came from: a translation defect
shows up as a last-digit drift long before it shows up as a wrong answer. The
expected values below were produced by

    gfortran -O2 -std=legacy rodas.f dc_decsol.f decsol.f

on the sources named in the module docstring, printed at 17 significant digits
so the decimal round-trips through a double. They are the SAME constants
cpp/tests/test_rodas.cpp carries, so the two ports are pinned to one oracle
rather than to each other.

THE STEP COUNTERS ARE CHECKED TOO, and are the more sensitive test: nfcn, njac
and nstep diverge on any change to the step-size controller or the error
estimate, while the final value can absorb a small perturbation and still look
right.

BIT EQUALITY IS ASSERTED FOR THE FOUR ROBERTSON CASES, because a tolerance
there would let a real translation drift pass unnoticed, which is the entire
failure mode this file exists to catch. They are short enough that the step-size
controller never diverges, and that was verified on both cluster images.

THE TWO VAN DER POL CASES CANNOT BE, and the reason is libm rather than this
port. The controller calls pow(err, 1/4) at the FAC1 and FACGUS lines, and
glibc's pow is not correctly rounded, so it is not portable across glibc
releases: feeding both images the same 400k inputs, pow(x, 0.25) hashes
differently under glibc 2.31 (20.04) and 2.35 (22.04) while sqrt(sqrt(x))
hashes identically. One ulp there re-times an accepted step, and from there the
step sequence, the work counters and the last digits of y all depend on the
host. Measured across both images:

  analytic  J: counters identical, y agrees to ~2e-13 relative
  numerical J: nfcn 32755 vs 32663, nstep 5467 vs 5449, y to ~1.1e-9

Both solves are run at rtol 1e-8, so agreement is asserted at 1e-8. Asserting
tighter than the accuracy actually asked of the integrator would be pinning the
host's libm rather than this port -- which is what the bit exact form of these
two cases was doing, and it is why they failed on every 20.04 worker.
Do NOT tighten these back without removing the pow call first: replacing
pow(x,1/4) by sqrt(sqrt(x)) is bit portable (sqrt is correctly rounded by
IEEE 754) and was measured to give identical results on both images, but it
moves the answer off the vendored Fortran's, so it is a separate decision.
This mirrors cpp/tests/test_rodas.cpp, which reached the same split.

The cases reach three different IJOB paths in DECOMR/SLVROD -- identity, banded
and full mass matrix -- crossed with analytic and numerical Jacobians, because
those are separate code paths and a transliteration error can hit one and not
the others. The singular mass matrix (index-1 DAE) is the case LINE needs; the
rest guard it against regressions elsewhere in the file.
"""

import numpy as np
import pytest

from line_solver.solvers.solver_fld.ode.rodas import rodas


# --- Robertson's problem as an index-1 DAE:  M y' = f(y),  M = diag(1,1,0).
#     The third equation is the algebraic constraint y1+y2+y3 = 1.
def _frob(x, y, f):
    f[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2]
    f[1] = 0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * y[1] * y[1]
    f[2] = y[0] + y[1] + y[2] - 1.0


def _jrob(x, y, dfy):
    dfy[0, 0] = -0.04
    dfy[0, 1] = 1.0e4 * y[2]
    dfy[0, 2] = 1.0e4 * y[1]
    dfy[1, 0] = 0.04
    dfy[1, 1] = -1.0e4 * y[2] - 6.0e7 * y[1]
    dfy[1, 2] = -1.0e4 * y[1]
    dfy[2, 0] = 1.0
    dfy[2, 1] = 1.0
    dfy[2, 2] = 1.0


def _mdiag3(am):
    """Banded storage, MLMAS=MUMAS=0, i.e. the diagonal -- reaches IJOB=3."""
    am[0, 0] = 1.0
    am[0, 1] = 1.0
    am[0, 2] = 0.0


def _mfull3(am):
    """The same mass matrix stored full -- reaches IJOB=5, and must agree."""
    am[...] = 0.0
    am[0, 0] = 1.0
    am[1, 1] = 1.0


# --- Van der Pol at eps=1e-6: stiff, M = identity, reaches IJOB=1.
def _fvdp(x, y, f):
    eps = 1.0e-6
    f[0] = y[1]
    f[1] = ((1.0 - y[0] * y[0]) * y[1] - y[0]) / eps


def _jvdp(x, y, dfy):
    eps = 1.0e-6
    dfy[0, 0] = 0.0
    dfy[0, 1] = 1.0
    dfy[1, 0] = (-2.0 * y[0] * y[1] - 1.0) / eps
    dfy[1, 1] = (1.0 - y[0] * y[0]) / eps


def _run_rob(ijac, mass_full, mlmas, mumas, rt, at):
    y = np.array([1.0, 0.0, 0.0])
    return rodas(3, _frob, 0.0, y, 0.4, h=1.0e-6, rtol=rt, atol=at, itol=0,
                 jac=_jrob, ijac=ijac, mljac=3, mujac=3, ifcn=0,
                 mas=(_mfull3 if mass_full else _mdiag3), imas=1,
                 mlmas=mlmas, mumas=mumas, iout=0)


def _run_vdp(ijac, rt, at):
    y = np.array([2.0, -0.6])
    return rodas(2, _fvdp, 0.0, y, 2.0, h=1.0e-6, rtol=rt, atol=at, itol=0,
                 jac=_jvdp, ijac=ijac, mljac=2, mujac=2, ifcn=0,
                 imas=0, mlmas=0, mumas=0, iout=0)


def test_dae_banded_mass_numerical_jacobian():
    r = _run_rob(0, False, 0, 0, 1.0e-8, 1.0e-10)
    assert r.idid == 1
    assert r.nfcn == 281
    assert r.njac == 46
    assert r.nstep == 47
    assert r.naccpt == 46
    assert r.y[0] == 0.98517211385739623
    assert r.y[1] == 0.33863953610620744e-4
    assert r.y[2] == 0.14794022188993190e-1


def test_dae_banded_mass_analytic_jacobian():
    r = _run_rob(1, False, 0, 0, 1.0e-8, 1.0e-10)
    assert r.idid == 1
    assert r.nfcn == 281
    assert r.njac == 46
    assert r.nstep == 47
    assert r.y[0] == 0.98517211385689285
    assert r.y[1] == 0.33863953595594975e-4
    assert r.y[2] == 0.14794022189511580e-1


def test_dae_tight_tolerance():
    r = _run_rob(1, False, 0, 0, 1.0e-11, 1.0e-13)
    assert r.idid == 1
    assert r.nfcn == 1510
    assert r.njac == 250
    assert r.nstep == 252
    assert r.y[0] == 0.98517211386097781
    assert r.y[1] == 0.33863953787836347e-4
    assert r.y[2] == 0.14794022185234504e-1


def test_full_mass_matrix_matches_banded():
    """Same problem, different storage, therefore a different IJOB path through
    DECOMR and SLVROD. Agreement here is what says the two paths were ported
    correctly, without needing an external reference for either."""
    full = _run_rob(1, True, 3, 3, 1.0e-8, 1.0e-10)
    band = _run_rob(1, False, 0, 0, 1.0e-8, 1.0e-10)
    assert full.idid == 1
    assert full.nfcn == band.nfcn
    assert full.nstep == band.nstep
    for i in range(3):
        assert full.y[i] == band.y[i]
    assert full.y[0] == 0.98517211385689285


def _check_count(got, want):
    """Work counters move with the step sequence, so they are bounded rather
    than fixed. The measured spread is under 0.4%; 5% still catches a real
    regression in the work a solve costs."""
    assert 0.95 * want <= got <= 1.05 * want


def test_stiff_ode_identity_mass_analytic_jacobian():
    r = _run_vdp(1, 1.0e-8, 1.0e-10)
    assert r.idid == 1
    _check_count(r.nfcn, 32446)
    _check_count(r.njac, 5401)
    _check_count(r.nstep, 5409)
    assert r.y[0] == pytest.approx(0.17061674639889246e1, rel=1e-8)
    assert r.y[1] == pytest.approx(-0.89280998821562574, rel=1e-8)


def test_stiff_ode_identity_mass_numerical_jacobian():
    r = _run_vdp(0, 1.0e-8, 1.0e-10)
    assert r.idid == 1
    _check_count(r.nfcn, 32755)
    _check_count(r.njac, 5420)
    _check_count(r.nstep, 5467)
    assert r.y[0] == pytest.approx(0.17061674637556288e1, rel=1e-8)
    assert r.y[1] == pytest.approx(-0.89280998834260528, rel=1e-8)


def test_dense_output_reproduces_the_endpoint_and_the_constraint():
    """CONTRO, the interpolant the fluid DAE reads its output grid off.

    Two things are asserted because two things can break independently: the
    interpolant must reproduce the integrator's own value where the two are
    evaluated at the same point (the right endpoint of the last accepted step),
    and it must satisfy the ALGEBRAIC row everywhere in between -- y1+y2+y3 = 1
    is an equation of the DAE, not a quantity that drifts, so an interpolant
    that only got the differential rows right would fail here and nowhere else.
    """
    grid = [0.05 * k for k in range(1, 9)]
    seen = []

    def solout(nr, xold, x, y, dense):
        while len(seen) < len(grid) and grid[len(seen)] <= x + 1e-13:
            tq = grid[len(seen)]
            # nr <= 1 is the call rodas makes BEFORE any step, where cont holds
            # no coefficients yet and the state at that point is y itself.
            if nr <= 1:
                seen.append((tq, np.array(y)))
            else:
                seen.append((tq, np.array([dense.value(i, tq)
                                           for i in range(3)])))
        return 0

    y = np.array([1.0, 0.0, 0.0])
    r = rodas(3, _frob, 0.0, y, 0.4, h=1.0e-6, rtol=1.0e-8, atol=1.0e-10,
              itol=0, jac=_jrob, ijac=1, mljac=3, mujac=3, ifcn=0,
              mas=_mdiag3, imas=1, mlmas=0, mumas=0,
              solout=solout, iout=1)
    assert r.idid == 1
    assert len(seen) == len(grid)
    # the algebraic row, at every interpolated point
    for tq, xs in seen:
        assert abs(xs[0] + xs[1] + xs[2] - 1.0) < 1e-9, "constraint at t=%g" % tq
    # dense output at the horizon agrees with the integrator's own y there
    tail = []

    def solout_end(nr, xold, x, y, dense):
        if abs(x - 0.4) < 1e-15 and nr > 1:
            tail.append(np.array([dense.value(i, 0.4) for i in range(3)]))
        return 0

    y2 = np.array([1.0, 0.0, 0.0])
    r2 = rodas(3, _frob, 0.0, y2, 0.4, h=1.0e-6, rtol=1.0e-8, atol=1.0e-10,
               itol=0, jac=_jrob, ijac=1, mljac=3, mujac=3, ifcn=0,
               mas=_mdiag3, imas=1, mlmas=0, mumas=0,
               solout=solout_end, iout=1)
    assert r2.idid == 1
    assert tail, "rodas never landed on the horizon"
    np.testing.assert_allclose(tail[-1], r2.y, rtol=1e-12, atol=1e-14)


def test_iout_does_not_change_the_trajectory():
    """Turning dense output on must cost steps, not accuracy: `cont` is written
    from the stage values that already exist, so the step sequence is the same
    one. A port that fed CONT back into the step would fail here."""
    quiet = _run_rob(1, False, 0, 0, 1.0e-8, 1.0e-10)

    y = np.array([1.0, 0.0, 0.0])
    loud = rodas(3, _frob, 0.0, y, 0.4, h=1.0e-6, rtol=1.0e-8, atol=1.0e-10,
                 itol=0, jac=_jrob, ijac=1, mljac=3, mujac=3, ifcn=0,
                 mas=_mdiag3, imas=1, mlmas=0, mumas=0,
                 solout=lambda *a: 0, iout=1)
    assert loud.nstep == quiet.nstep
    assert loud.nfcn == quiet.nfcn
    for i in range(3):
        assert loud.y[i] == quiet.y[i]


def test_solout_can_stop_the_integration():
    def solout(nr, xold, x, y, dense):
        return -1 if x > 0.1 else 0

    y = np.array([1.0, 0.0, 0.0])
    r = rodas(3, _frob, 0.0, y, 0.4, h=1.0e-6, rtol=1.0e-8, atol=1.0e-10,
              itol=0, jac=_jrob, ijac=1, mljac=3, mujac=3, ifcn=0,
              mas=_mdiag3, imas=1, mlmas=0, mumas=0, solout=solout, iout=1)
    assert r.idid == 2
    assert r.x < 0.4


def test_tolerances_too_small_is_refused_by_name():
    y = np.array([1.0, 0.0, 0.0])
    with pytest.raises(Exception):
        rodas(3, _frob, 0.0, y, 0.4, h=1.0e-6, rtol=1.0e-18, atol=1.0e-20,
              itol=0, jac=_jrob, ijac=1, mljac=3, mujac=3, ifcn=0,
              mas=_mdiag3, imas=1, mlmas=0, mumas=0, iout=0)
