"""
Grid-independent sojourn-weighted average of a cache mean-field transient.

Python twin of MATLAB ``cache_sojourn_ode.m``, JAR ``CacheSojourn`` and C++
``cache_sojourn.h``. For a drift dx/dt = f(x) and a phase-type clock with row
phase vector phi, dphi/dt = phi A, density g(t) = phi(t) c, it returns

    xbar = int_{t0}^{t1} x g dt / int_{t0}^{t1} g dt,   wtot = int_{t0}^{t1} g dt,

by augmenting the state with phi, int g x and int g, so the integrator itself
carries the integral. The value is accurate to the ODE tolerance and does NOT
depend on any output grid: a Riemann sum over the grid did, and that is what made
the adaptive (MATLAB) and fixed-grid (JAR/Python/C++) ENV cache mean fields
differ in the third digit.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm


def cache_sojourn_clock(D0, pie, lam, t0=0.0):
    """Holding-time clock of a stage whose holding time is the MAP with sub-generator
    D0 and initial vector pie, in a drift time unit that is ``lam`` times real time:
    A = D0/lam, c = -D0 1/lam, phi0 = pie expm(A t0), so int g = F(t1/lam) - F(t0/lam).
    """
    D0 = np.asarray(D0, dtype=float)
    A = D0 / lam
    c = -D0 @ np.ones(D0.shape[0]) / lam
    phi0 = np.asarray(pie, dtype=float).ravel()
    if t0 != 0.0:
        phi0 = phi0 @ expm(A * t0)
    return {'A': A, 'phi0': phi0, 'c': c}


def cache_sojourn_ode(drift, t0, t1, x0, clock, rtol=1e-8, atol=1e-10):
    """Integrate dx/dt = drift(x) over [t0, t1] from x0 together with the clock and
    return dict(xbar=..., wtot=...); xbar is None when wtot is not positive."""
    x0 = np.asarray(x0, dtype=float).ravel()
    n = x0.size
    A = np.asarray(clock['A'], dtype=float)
    c = np.asarray(clock['c'], dtype=float).ravel()
    phi0 = np.asarray(clock['phi0'], dtype=float).ravel()
    nph = A.shape[0]
    At = A.T

    def rhs(_t, z):
        x = z[:n]
        phi = z[n:n + nph]
        g = float(phi @ c)
        return np.concatenate([drift(x), At @ phi, g * x, [g]])

    z0 = np.concatenate([x0, phi0, np.zeros(n), [0.0]])
    sol = solve_ivp(rhs, (t0, t1), z0, method='LSODA', rtol=rtol, atol=atol)
    z1 = sol.y[:, -1]
    wtot = float(z1[-1])
    xbar = z1[n + nph:n + nph + n] / wtot if wtot > 0 else None
    return {'xbar': xbar, 'wtot': wtot}
