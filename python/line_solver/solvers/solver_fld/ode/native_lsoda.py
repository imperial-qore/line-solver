"""The in-tree LSODA (`line_solver.lib.lsoda`) behind scipy's OdeSolver protocol.

`lib/lsoda.py` is a transcription of liblsoda checked against the C reference
vectors in `python/tests/test_lsoda.py`. This module is what lets the fluid path
drive it: an `OdeSolver` subclass, so it can be handed to `solve_ivp` as
``method=NativeLSODA`` or stepped by hand the way
`solver_fld/methods/closing.py` steps its integrator to reproduce
``odeset('NonNegative')``.

NOTHING SELECTS IT BY DEFAULT. The fluid path keeps scipy's compiled LSODA and
BDF; this is reached only through ``options.odesolver``, which is the Python
twin of assigning `options.odesolvers.accurateStiffOdeSolver = @lsoda_*` in
MATLAB.

Two classes are exported because the stiff slot must not be answered by the
Adams half:

  NativeLSODA       - the auto-switcher, as the library ships it
  NativeLSODAStiff  - starts on BDF and stays there (`force_stiff`), the same
                      pin as LSODA.setForceStiff(true) in the JAR

The dense output is LSODA's own `intdy`, i.e. the Nordsieck history evaluated by
Horner, and it SNAPSHOTS that history at construction so the object stays valid
after the solver has stepped on.
"""

import numpy as np
from scipy.integrate import OdeSolver, DenseOutput

from ....lib.lsoda import LSODAStepper, LSODAError


class LsodaDenseOutput(DenseOutput):
    """`intdy` at k = 0 on a snapshot of the Nordsieck history."""

    def __init__(self, t_old, t, stepper):
        super().__init__(t_old, t)
        c = stepper.c
        self.nq = c.nq
        self.tn = c.tn
        self.h = c.h
        self.neq = stepper.neq
        # rows 1..nq+1 are the history; row 0 is unused, as in the C
        self.yh = [np.array(c.yh[j][1:stepper.neq + 1], copy=True)
                   for j in range(0, c.nq + 2)]

    def _horner(self, t):
        s = (t - self.tn) / self.h
        dky = self.yh[self.nq + 1].copy()
        for j in range(self.nq - 1, -1, -1):
            dky = self.yh[j + 1] + s * dky
        return dky

    def _call_impl(self, t):
        t = np.asarray(t, dtype=float)
        if t.ndim == 0:
            return self._horner(float(t))
        out = np.empty((self.neq, t.size))
        for k in range(t.size):
            out[:, k] = self._horner(float(t.flat[k]))
        return out


class NativeLSODA(OdeSolver):
    """LSODA from `line_solver.lib.lsoda`, stepped through scipy's protocol.

    Each `step()` is one call of the C driver with itask = 5: one internal step,
    never past `t_bound`. That is the mode scipy's own LSODA wrapper uses, and
    it is what makes the solver land exactly on the horizon rather than
    overshooting it and interpolating back.

    Parameters beyond the OdeSolver contract:

    max_steps : int
        The step budget of a single `advance` call before LSODA gives up. It is
        a give-up threshold, not a switch: nothing follows it here.
    force_stiff : bool
        Class attribute on `NativeLSODAStiff`; see the module docstring.
    """

    force_stiff = False

    def __init__(self, fun, t0, y0, t_bound, first_step=None, min_step=0.0,
                 max_step=np.inf, rtol=1e-3, atol=1e-6, max_steps=100000,
                 vectorized=False, **extraneous):
        super().__init__(fun, t0, y0, t_bound, vectorized, support_complex=False)
        if max_step <= 0 or min_step < 0:
            raise ValueError('max_step must be positive and min_step nonnegative')
        hmax = 0.0 if not np.isfinite(max_step) else float(max_step)
        h0 = 0.0
        if first_step is not None:
            if first_step <= 0:
                raise ValueError('first_step must be positive')
            h0 = self.direction * float(first_step)

        def _rhs(t, yv, ydot, _data):
            ydot[:] = self.fun(t, yv)

        self._stepper = LSODAStepper(
            _rhs, self.y, self.t, rtol=rtol, atol=atol,
            max_steps=int(max_steps), hmax=hmax, hmin=float(min_step), h0=h0,
            tcrit=float(t_bound), force_stiff=self.force_stiff)

    @property
    def h_abs(self):
        """The step size just taken, which `closing.py` carries into a restart."""
        return abs(self._stepper.h)

    def _step_impl(self):
        try:
            state = self._stepper.advance(self.t_bound, itask=5)
        except LSODAError as err:
            return False, str(err)
        if state <= 0:
            return False, self._stepper.message
        self.t = self._stepper.t
        self.y = self._stepper.y
        # njev is scipy's Jacobian counter; nfev is maintained by the wrapped fun
        self.njev = self._stepper.nje
        return True, None

    def _dense_output_impl(self):
        return LsodaDenseOutput(self.t_old, self.t, self._stepper)


class NativeLSODAStiff(NativeLSODA):
    """`NativeLSODA` with the BDF half pinned, for the stiff slot.

    WHERE THE PIN IS NOT SAFE, measured rather than assumed: pinning BDF means
    correcting with a finite-difference Jacobian from the very first step, and
    LSODA sizes that difference as ``r = max(sqrt(eps)*|y_j|, r0/ewt_j)``. A
    component that is EXACTLY zero at t0 under an atol far below the scale of
    the drift takes the second branch with r ~ 1e-19, so its Jacobian column is
    rounding noise and the corrector converges to nonsense. Robertson started at
    y = (1, 0, 0) with atol = (1e-6, 1e-10, 1e-6) does exactly that here: the
    solve runs away to y1 = -1.5e7 and still reports success, and it flips on a
    last-bit change to the right-hand side. The auto-switcher escapes it by
    starting on Adams, which needs no Jacobian, and switching once the solution
    is smooth; starting the pinned solver from t = 0.4 instead is also fine.

    That is not the shape of a fluid drift, and the pin is measured sound on
    HIRES, van der Pol at mu = 1000, the Oregonator and the fluid models
    (`tests/test_lsoda_fluid.py`). Note also that Robertson at rtol = 1e-4 is
    knife-edge for this family generally: the compiled C itself runs away at a
    uniform atol = 1e-6, with the pinned solver the one that stays correct.
    """

    force_stiff = True
