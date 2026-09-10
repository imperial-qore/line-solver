"""
Pure Python port of the LSODA ODE solver from liblsoda (C version).

LSODA solves initial value problems for stiff or nonstiff systems of
first order ODEs, dy/dt = f(t,y), with automatic method switching
between Adams (nonstiff) and BDF (stiff) methods.

Original authors:
  Linda R. Petzold and Alan C. Hindmarsh, Lawrence Livermore National Laboratory.

C version by Hon Wah Tam (Wolfram Research) and Yu Feng (Carnegie Mellon).

The MIT License
Copyright (c) 2011 McWilliam Cosmology Center, Carnegie Mellon University.

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Python port for the LINE solver project.
"""

import math
import numpy as np

# Constants
ETA = 2.2204460492503131e-16
SQRTETA = 1.4901161193847656e-08
CCMAX = 0.3
MAXCOR = 3
MSBP = 20
MXNCF = 10
RATIO = 5.0

# Stability region limits for Adams method (orders 0..12)
sm1 = [0., 0.5, 0.575, 0.55, 0.45, 0.35, 0.25, 0.2, 0.15, 0.1, 0.075, 0.05, 0.025]

# Method switch constants (precomputed from tesco/elco)
_cm1 = [
    0.0, 2.0, 5.999999999999999, 4.0,
    1.5777777777777777, 0.44444444444444453, 0.09721974866042459,
    0.01755731922398589, 0.002677042079828849, 0.00035127324808082747,
    4.013484869498161e-05, 4.011251133543618e-06, 3.541120754086548e-07,
]

_cm2 = [
    0.0, 2.0, 1.5, 0.6666666666666667,
    0.20833333333333337, 0.04999999999999999, 0.09721974866042459,
    0.01755731922398589, 0.002677042079828849, 0.00035127324808082747,
    4.013484869498161e-05, 4.011251133543618e-06, 3.541120754086548e-07,
]


class LSODAError(Exception):
    """Exception raised for LSODA solver errors."""
    pass


# ---- BLAS-like helper functions (1-based indexing on numpy arrays) ----
# These operate on arrays that are allocated with size n+1, where index 0
# is unused, matching the original C/Fortran convention.

def _vmnorm(n, v, w):
    """Weighted max-norm: max(|v[i]| * w[i]) for i=1..n."""
    vm = 0.0
    for i in range(1, n + 1):
        vm = max(vm, abs(v[i]) * w[i])
    return vm


def _fnorm(n, a, w):
    """Weighted matrix norm consistent with vmnorm.
    fnorm = max_i (w[i] * sum_j |a[i][j]| / w[j])
    """
    an = 0.0
    for i in range(1, n + 1):
        s = 0.0
        for j in range(1, n + 1):
            s += abs(a[i][j]) / w[j]
        an = max(an, s * w[i])
    return an


def _idamax(n, dx, incx):
    """Index of element with max absolute value (1-based)."""
    if n <= 0:
        return 0
    if n == 1 or incx <= 0:
        return 1
    if incx != 1:
        dmax = abs(dx[1])
        xindex = 1
        ii = 2
        for i in range(1 + incx, n * incx + 1, incx):
            xmag = abs(dx[i])
            if xmag > dmax:
                xindex = ii
                dmax = xmag
            ii += 1
        return xindex
    # incx == 1
    dmax = abs(dx[1])
    xindex = 1
    for i in range(2, n + 1):
        xmag = abs(dx[i])
        if xmag > dmax:
            xindex = i
            dmax = xmag
    return xindex


def _dscal(n, da, dx, incx):
    """Scale vector: dx = da * dx."""
    if n <= 0:
        return
    if incx != 1:
        for i in range(1, n * incx + 1, incx):
            dx[i] = da * dx[i]
        return
    for i in range(1, n + 1):
        dx[i] = da * dx[i]


def _daxpy(n, da, dx, incx, dy, incy):
    """dy = da * dx + dy."""
    if n <= 0 or da == 0.0:
        return
    if incx != incy or incx < 1:
        ix = 1
        iy = 1
        if incx < 0:
            ix = (-n + 1) * incx + 1
        if incy < 0:
            iy = (-n + 1) * incy + 1
        for _ in range(n):
            dy[iy] += da * dx[ix]
            ix += incx
            iy += incy
        return
    if incx == 1:
        for i in range(1, n + 1):
            dy[i] += da * dx[i]
        return
    for i in range(1, n * incx + 1, incx):
        dy[i] = da * dx[i] + dy[i]


def _ddot(n, dx, incx, dy, incy):
    """Inner product dx . dy."""
    dotprod = 0.0
    if n <= 0:
        return dotprod
    if incx != incy or incx < 1:
        ix = 1
        iy = 1
        if incx < 0:
            ix = (-n + 1) * incx + 1
        if incy < 0:
            iy = (-n + 1) * incy + 1
        for _ in range(n):
            dotprod += dx[ix] * dy[iy]
            ix += incx
            iy += incy
        return dotprod
    if incx == 1:
        for i in range(1, n + 1):
            dotprod += dx[i] * dy[i]
        return dotprod
    for i in range(1, n * incx + 1, incx):
        dotprod += dx[i] * dy[i]
    return dotprod


def _dgefa(a, n, ipvt):
    """Gaussian elimination with partial pivoting (LU factorization).
    Returns info (0 = success, k = singular at row k).

    Note: The C version uses pointer arithmetic (a[k]+k-1, a[k]+k) to shift
    array start positions. Here we use explicit index offsets instead.
    """
    info = 0
    for k in range(1, n):
        # Find pivot: idamax on a[k][k..n] (n-k+1 elements starting at k)
        dmax = abs(a[k][k])
        j = k
        for ii in range(k + 1, n + 1):
            xmag = abs(a[k][ii])
            if xmag > dmax:
                j = ii
                dmax = xmag
        ipvt[k] = j
        if a[k][j] == 0.0:
            info = k
            continue
        if j != k:
            a[k][j], a[k][k] = a[k][k], a[k][j]
        t = -1.0 / a[k][k]
        # dscal on a[k][k+1..n]
        for ii in range(k + 1, n + 1):
            a[k][ii] *= t
        # Column elimination
        for i in range(k + 1, n + 1):
            t = a[i][j]
            if j != k:
                a[i][j] = a[i][k]
                a[i][k] = t
            # daxpy: a[i][k+1..n] += t * a[k][k+1..n]
            for ii in range(k + 1, n + 1):
                a[i][ii] += t * a[k][ii]
    ipvt[n] = n
    if a[n][n] == 0.0:
        info = n
    return info


def _dgesl(a, n, ipvt, b, job):
    """Solve a*x=b (job=0) or a^T*x=b (job!=0) using LU factors from dgefa.

    Note: The C version uses pointer arithmetic (a[k]+k, b+k) to shift
    array start positions. Here we use explicit index offsets instead.
    """
    if job == 0:
        # Forward solve: L * y = b
        for k in range(1, n + 1):
            # ddot(k-1, a[k][1..k-1], b[1..k-1])
            t = 0.0
            for ii in range(1, k):
                t += a[k][ii] * b[ii]
            b[k] = (b[k] - t) / a[k][k]
        # Back solve: U * x = y
        for k in range(n - 1, 0, -1):
            # ddot(n-k, a[k][k+1..n], b[k+1..n])
            t = 0.0
            for ii in range(k + 1, n + 1):
                t += a[k][ii] * b[ii]
            b[k] += t
            j = ipvt[k]
            if j != k:
                b[j], b[k] = b[k], b[j]
    else:
        # Forward: solve U^T * y = b
        for k in range(1, n):
            j = ipvt[k]
            t = b[j]
            if j != k:
                b[j] = b[k]
                b[k] = t
            # daxpy(n-k, t, a[k][k+1..n], b[k+1..n])
            for ii in range(k + 1, n + 1):
                b[ii] += t * a[k][ii]
        # Back: solve L^T * x = y
        for k in range(n, 0, -1):
            b[k] /= a[k][k]
            t = -b[k]
            # daxpy(k-1, t, a[k][1..k-1], b[1..k-1])
            for ii in range(1, k):
                b[ii] += t * a[k][ii]


# ---- Internal LSODA routines ----

class _LSODACommon:
    """Internal solver state, equivalent to lsoda_common_t."""
    __slots__ = [
        'yh', 'wm', 'ewt', 'savf', 'acor', 'ipvt',
        'h', 'hu', 'rc', 'tn', 'tsw', 'pdnorm',
        'crate', 'el', 'elco', 'tesco',
        'hold', 'rmax', 'pdest', 'pdlast',
        'ialth', 'ipup', 'nslp', 'icount', 'irflag',
        'imxer', 'illin', 'nhnil', 'nslast',
        'jcur', 'meth', 'mused', 'nq', 'nst',
        'ncf', 'nfe', 'nje', 'nqu', 'miter',
        'hmin_val', 'told_corr',
    ]

    def __init__(self):
        for attr in self.__slots__:
            setattr(self, attr, 0 if attr not in ('h', 'hu', 'rc', 'tn', 'tsw',
                                                    'pdnorm', 'crate', 'hold',
                                                    'rmax', 'pdest', 'pdlast') else 0.0)
        self.el = None
        self.elco = None
        self.tesco = None
        self.yh = None
        self.wm = None
        self.ewt = None
        self.savf = None
        self.acor = None
        self.ipvt = None


class LSODAOptions:
    """Options for the LSODA solver."""
    __slots__ = [
        'ixpr', 'mxstep', 'mxhnil', 'mxordn', 'mxords',
        'tcrit', 'h0', 'hmax', 'hmin', 'hmxi', 'itask',
        'rtol', 'atol', 'force_stiff',
    ]

    def __init__(self):
        self.ixpr = 0
        self.mxstep = 0
        self.mxhnil = 0
        self.mxordn = 0
        self.mxords = 0
        self.tcrit = 0.0
        self.h0 = 0.0
        self.hmax = 0.0
        self.hmin = 0.0
        self.hmxi = 0.0
        self.itask = 0
        self.rtol = None
        self.atol = None
        # Start on BDF and stay there, never switching to Adams. The Adams half
        # loses its stability bound at a fixed point: the corrector converges on
        # the `del <= 100*pnorm*ETA` branch before `pdest` is ever formed, so
        # `scaleh`'s pdh guard never binds and h runs up to hmax. Same pin as
        # LSODA.setForceStiff(true) in the JAR and lsoda_matlab's forceStiff.
        self.force_stiff = False


def _alloc_mem(neq, mxordn, mxords):
    """Allocate internal arrays (1-based indexing: index 0 is unused)."""
    lenyh = 1 + max(mxordn, mxords)
    c = _LSODACommon()
    # yh[0..lenyh][0..neq] -- history array
    c.yh = [np.zeros(neq + 1) for _ in range(lenyh + 1)]
    # wm[0..neq][0..neq] -- work matrix
    c.wm = [np.zeros(neq + 1) for _ in range(neq + 1)]
    # 1-D work arrays
    c.ewt = np.zeros(neq + 1)
    c.savf = np.zeros(neq + 1)
    c.acor = np.zeros(neq + 1)
    c.ipvt = np.zeros(neq + 1, dtype=np.int64)
    # Coefficient arrays
    c.el = np.zeros(14)
    c.elco = [np.zeros(14) for _ in range(14)]
    c.tesco = [np.zeros(4) for _ in range(14)]
    return c


def _cfode(c, meth_val):
    """Compute method coefficients for Adams (meth=1) or BDF (meth=2)."""
    if meth_val == 1:
        c.elco[1][1] = 1.0
        c.elco[1][2] = 1.0
        c.tesco[1][1] = 0.0
        c.tesco[1][2] = 2.0
        c.tesco[2][1] = 1.0
        c.tesco[12][3] = 0.0
        pc = np.zeros(14)
        pc[1] = 1.0
        rqfac = 1.0
        for nq in range(2, 13):
            rq1fac = rqfac
            rqfac /= float(nq)
            nqm1 = nq - 1
            fnqm1 = float(nqm1)
            nqp1 = nq + 1
            pc[nq] = 0.0
            for i in range(nq, 1, -1):
                pc[i] = pc[i - 1] + fnqm1 * pc[i]
            pc[1] = fnqm1 * pc[1]
            pint = pc[1]
            xpin = pc[1] / 2.0
            tsign = 1.0
            for i in range(2, nq + 1):
                tsign = -tsign
                pint += tsign * pc[i] / float(i)
                xpin += tsign * pc[i] / float(i + 1)
            c.elco[nq][1] = pint * rq1fac
            c.elco[nq][2] = 1.0
            for i in range(2, nq + 1):
                c.elco[nq][i + 1] = rq1fac * pc[i] / float(i)
            agamq = rqfac * xpin
            ragq = 1.0 / agamq
            c.tesco[nq][2] = ragq
            if nq < 12:
                c.tesco[nqp1][1] = ragq * rqfac / float(nqp1)
            c.tesco[nqm1][3] = ragq
    else:
        # BDF method (meth == 2)
        pc = np.zeros(14)
        pc[1] = 1.0
        rq1fac = 1.0
        for nq in range(1, 6):
            fnq = float(nq)
            nqp1 = nq + 1
            pc[nqp1] = 0.0
            for i in range(nq + 1, 1, -1):
                pc[i] = pc[i - 1] + fnq * pc[i]
            pc[1] *= fnq
            for i in range(1, nqp1 + 1):
                c.elco[nq][i] = pc[i] / pc[2]
            c.elco[nq][2] = 1.0
            c.tesco[nq][1] = rq1fac
            c.tesco[nq][2] = float(nqp1) / c.elco[nq][1]
            c.tesco[nq][3] = float(nq + 2) / c.elco[nq][1]
            rq1fac /= fnq


def _resetcoeff(c):
    """Reset el vector when order nq changes."""
    el0 = c.el[1]
    for i in range(1, c.nq + 2):
        c.el[i] = c.elco[c.nq][i]
    c.rc = c.rc * c.el[1] / el0


def _scaleh(c, rh, neq, opt):
    """Scale step size h by factor rh and rescale yh array."""
    rh = min(rh, c.rmax)
    rh = rh / max(1.0, abs(c.h) * opt.hmxi * rh)
    # Stability region check for Adams method
    if c.meth == 1:
        c.irflag = 0
        pdh = max(abs(c.h) * c.pdlast, 0.000001)
        if rh * pdh * 1.00001 >= sm1[c.nq]:
            rh = sm1[c.nq] / pdh
            c.irflag = 1
    r = 1.0
    for j in range(2, c.nq + 2):
        r *= rh
        for i in range(1, neq + 1):
            c.yh[j][i] *= r
    c.h *= rh
    c.rc *= rh
    c.ialth = c.nq + 1


def _intdy(c, t, k, dky, neq):
    """Interpolate k-th derivative at time t."""
    if k < 0 or k > c.nq:
        return -1
    tp = c.tn - c.hu - 100.0 * ETA * (c.tn + c.hu)
    if (t - tp) * (t - c.tn) > 0.0:
        return -2
    s = (t - c.tn) / c.h
    ic = 1
    for jj in range(c.nq + 1 - k, c.nq + 1):
        ic *= jj
    co = float(ic)
    for i in range(1, neq + 1):
        dky[i] = co * c.yh[c.nq + 1][i]
    for j in range(c.nq - 1, k - 1, -1):
        jp1 = j + 1
        ic = 1
        for jj in range(jp1 - k, j + 1):
            ic *= jj
        co = float(ic)
        for i in range(1, neq + 1):
            dky[i] = co * c.yh[jp1][i] + s * dky[i]
    if k != 0:
        r = c.h ** float(-k)
        for i in range(1, neq + 1):
            dky[i] *= r
    return 0


def _solsy(c, y, neq):
    """Solve the linear system using LU factors."""
    _dgesl(c.wm, neq, c.ipvt, y, 0)
    return 1


def _prja(c, y, neq, f, data, opt=None):
    """Compute and process Jacobian P = I - h*el[1]*J by finite differences.

    UNDER force_stiff THE INCREMENT CARRIES MATLAB NUMJAC'S FLOOR. LSODA sizes
    the difference as max(sqrt(eps)*|y_j|, r0/ewt_j), and a component sitting at
    EXACTLY zero under a tight atol takes the second branch with r ~ 1e-19: the
    column is then rounding noise divided by that, the Newton matrix is
    meaningless, and the corrector converges to nonsense. The auto-switcher never
    meets this because it starts on Adams and only takes a Jacobian once the
    solution is smooth; pinning BDF takes one at t0, so it needs the rule ode15s
    uses, sqrt(eps)*max(|y_j|, atol_j/rtol_j), which is exactly why ode15s has
    never had this failure. Applied ONLY under the pin, so the auto-switcher and
    the C reference vectors it reproduces are untouched.
    """
    c.nje += 1
    hl0 = c.h * c.el[1]
    fac = _vmnorm(neq, c.savf, c.ewt)
    r0 = 1000.0 * abs(c.h) * ETA * float(neq) * fac
    if r0 == 0.0:
        r0 = 1.0
    pin = opt is not None and opt.force_stiff
    for j in range(1, neq + 1):
        yj = y[j]
        r = max(SQRTETA * abs(yj), r0 / c.ewt[j])
        if pin and opt.rtol[j] > 0.0:
            r = max(r, SQRTETA * max(abs(yj), opt.atol[j] / opt.rtol[j]))
        y[j] += r
        fac_val = -hl0 / r
        f(c.tn, y, c.acor, data)
        for i in range(1, neq + 1):
            c.wm[i][j] = (c.acor[i] - c.savf[i]) * fac_val
        y[j] = yj
    c.nfe += neq
    # Compute norm of Jacobian
    c.pdnorm = _fnorm(neq, c.wm, c.ewt) / abs(hl0)
    # Add identity matrix
    for i in range(1, neq + 1):
        c.wm[i][i] += 1.0
    # LU decomposition
    ier = _dgefa(c.wm, neq, c.ipvt)
    if ier != 0:
        return 0
    return 1


def _corfailure(c, told, neq):
    """Handle corrector convergence failure."""
    c.ncf += 1
    c.rmax = 2.0
    c.tn = told
    for j in range(c.nq, 0, -1):
        for i1 in range(j, c.nq + 1):
            for i in range(1, neq + 1):
                c.yh[i1][i] -= c.yh[i1 + 1][i]
    hmin = c.hmin_val
    if abs(c.h) <= hmin * 1.00001 or c.ncf == MXNCF:
        return 2
    c.ipup = c.miter
    return 1


def _correction(c, y, pnorm, neq, f, data, opt=None):
    """Corrector iteration. Returns (corflag, del_val, delp, m).
    corflag: 0=converged, 1=reduce h & redo, 2=failure.
    """
    m = 0
    rate = 0.0
    del_val = 0.0
    delp = 0.0
    for i in range(1, neq + 1):
        y[i] = c.yh[1][i]
    f(c.tn, y, c.savf, data)
    c.nfe += 1

    while True:
        if m == 0:
            if c.ipup > 0:
                ierpj = _prja(c, y, neq, f, data, opt)
                c.jcur = 1
                c.ipup = 0
                c.rc = 1.0
                c.nslp = c.nst
                c.crate = 0.7
                if not ierpj:
                    return _corfailure(c, c.told_corr, neq), del_val, delp, m
            for i in range(1, neq + 1):
                c.acor[i] = 0.0

        if c.miter == 0:
            # Functional iteration
            for i in range(1, neq + 1):
                c.savf[i] = c.h * c.savf[i] - c.yh[2][i]
                y[i] = c.savf[i] - c.acor[i]
            del_val = _vmnorm(neq, y, c.ewt)
            for i in range(1, neq + 1):
                y[i] = c.yh[1][i] + c.el[1] * c.savf[i]
                c.acor[i] = c.savf[i]
        else:
            # Chord method
            for i in range(1, neq + 1):
                y[i] = c.h * c.savf[i] - (c.yh[2][i] + c.acor[i])
            _solsy(c, y, neq)
            del_val = _vmnorm(neq, y, c.ewt)
            for i in range(1, neq + 1):
                c.acor[i] += y[i]
                y[i] = c.yh[1][i] + c.el[1] * c.acor[i]

        # Test for convergence
        if del_val <= 100.0 * pnorm * ETA:
            break
        if m != 0 or c.meth != 1:
            if m != 0:
                rm = 1024.0
                if del_val <= 1024.0 * delp:
                    rm = del_val / delp
                rate = max(rate, rm)
                c.crate = max(0.2 * c.crate, rm)
            conit = 0.5 / float(c.nq + 2)
            dcon = del_val * min(1.0, 1.5 * c.crate) / (c.tesco[c.nq][2] * conit)
            if dcon <= 1.0:
                c.pdest = max(c.pdest, rate / abs(c.h * c.el[1]))
                if c.pdest != 0.0:
                    c.pdlast = c.pdest
                break

        m += 1
        if m == MAXCOR or (m >= 2 and del_val > 2.0 * delp):
            if c.miter == 0 or c.jcur == 1:
                return _corfailure(c, c.told_corr, neq), del_val, delp, m
            c.ipup = c.miter
            m = 0
            rate = 0.0
            del_val = 0.0
            for i in range(1, neq + 1):
                y[i] = c.yh[1][i]
            f(c.tn, y, c.savf, data)
            c.nfe += 1
        else:
            delp = del_val
            f(c.tn, y, c.savf, data)
            c.nfe += 1

    return 0, del_val, delp, m


def _methodswitch(c, dsm, pnorm, neq, opt):
    """Consider switching between Adams and BDF methods. Returns rh or None."""
    if opt.force_stiff:
        return None
    mxordn = opt.mxordn
    mxords = opt.mxords

    if c.meth == 1:
        # Currently Adams, consider switching to BDF
        if c.nq > 5:
            return None
        if dsm <= 100.0 * pnorm * ETA or c.pdest == 0.0:
            if c.irflag == 0:
                return None
            rh2 = 2.0
            nqm2 = min(c.nq, mxords)
        else:
            exsm = 1.0 / float(c.nq + 1)
            rh1 = 1.0 / (1.2 * dsm ** exsm + 0.0000012)
            rh1it = 2.0 * rh1
            pdh = c.pdlast * abs(c.h)
            if pdh * rh1 > 0.00001:
                rh1it = sm1[c.nq] / pdh
            rh1 = min(rh1, rh1it)
            if c.nq > mxords:
                nqm2 = mxords
                lm2 = mxords + 1
                exm2 = 1.0 / float(lm2)
                lm2p1 = lm2 + 1
                dm2 = _vmnorm(neq, c.yh[lm2p1], c.ewt) / _cm2[mxords]
                rh2 = 1.0 / (1.2 * dm2 ** exm2 + 0.0000012)
            else:
                dm2 = dsm * (_cm1[c.nq] / _cm2[c.nq])
                rh2 = 1.0 / (1.2 * dm2 ** exsm + 0.0000012)
                nqm2 = c.nq
            if rh2 < RATIO * rh1:
                return None
        # Switch to BDF
        c.icount = 20
        c.meth = 2
        c.miter = 2
        c.pdlast = 0.0
        c.nq = nqm2
        return rh2
    else:
        # Currently BDF, consider switching to Adams
        exsm = 1.0 / float(c.nq + 1)
        if mxordn < c.nq:
            nqm1 = mxordn
            lm1 = mxordn + 1
            exm1 = 1.0 / float(lm1)
            lm1p1 = lm1 + 1
            dm1 = _vmnorm(neq, c.yh[lm1p1], c.ewt) / _cm1[mxordn]
            rh1 = 1.0 / (1.2 * dm1 ** exm1 + 0.0000012)
        else:
            dm1 = dsm * (_cm2[c.nq] / _cm1[c.nq])
            rh1 = 1.0 / (1.2 * dm1 ** exsm + 0.0000012)
            nqm1 = c.nq
            exm1 = exsm
        rh1it = 2.0 * rh1
        pdh = c.pdnorm * abs(c.h)
        if pdh * rh1 > 0.00001:
            rh1it = sm1[nqm1] / pdh
        rh1 = min(rh1, rh1it)
        rh2 = 1.0 / (1.2 * dsm ** exsm + 0.0000012)
        if rh1 * RATIO < 5.0 * rh2:
            return None
        alpha = max(0.001, rh1)
        dm1 *= alpha ** exm1
        if dm1 <= 1000.0 * ETA * pnorm:
            return None
        # Switch to Adams
        c.icount = 20
        c.meth = 1
        c.miter = 0
        c.pdlast = 0.0
        c.nq = nqm1
        return rh1


def _orderswitch(c, rhup, dsm, kflag, maxord, neq, opt):
    """Order selection. Returns (orderflag, rh).
    orderflag: 0=no change, 1=change h only, 2=change h and nq.
    """
    exsm = 1.0 / float(c.nq + 1)
    rhsm = 1.0 / (1.2 * dsm ** exsm + 0.0000012)
    rhdn = 0.0
    if c.nq != 1:
        ddn = _vmnorm(neq, c.yh[c.nq + 1], c.ewt) / c.tesco[c.nq][1]
        exdn = 1.0 / float(c.nq)
        rhdn = 1.0 / (1.3 * ddn ** exdn + 0.0000013)

    # Stability region limits for Adams
    if c.meth == 1:
        pdh = max(abs(c.h) * c.pdlast, 0.000001)
        if c.nq + 1 < maxord + 1:
            rhup = min(rhup, sm1[c.nq + 1] / pdh)
        rhsm = min(rhsm, sm1[c.nq] / pdh)
        if c.nq > 1:
            rhdn = min(rhdn, sm1[c.nq - 1] / pdh)
        c.pdest = 0.0

    if rhsm >= rhup:
        if rhsm >= rhdn:
            newq = c.nq
            rh = rhsm
        else:
            newq = c.nq - 1
            rh = rhdn
            if kflag < 0 and rh > 1.0:
                rh = 1.0
    else:
        if rhup <= rhdn:
            newq = c.nq - 1
            rh = rhdn
            if kflag < 0 and rh > 1.0:
                rh = 1.0
        else:
            rh = rhup
            if rh >= 1.1:
                r = c.el[c.nq + 1] / float(c.nq + 1)
                c.nq += 1
                for i in range(1, neq + 1):
                    c.yh[c.nq + 1][i] = c.acor[i] * r
                return 2, rh
            else:
                c.ialth = 3
                return 0, rh

    # Stability bypass test for Adams
    if c.meth == 1:
        if rh * pdh * 1.00001 < sm1[newq]:
            if kflag == 0 and rh < 1.1:
                c.ialth = 3
                return 0, rh
    else:
        if kflag == 0 and rh < 1.1:
            c.ialth = 3
            return 0, rh

    if kflag <= -2:
        rh = min(rh, 0.2)

    if newq == c.nq:
        return 1, rh
    c.nq = newq
    return 2, rh


def _stoda(c, y, jstart, neq, f, data, opt):
    """One step of the integration. Returns kflag."""
    kflag = 0
    told = c.tn
    c.ncf = 0
    c.hmin_val = opt.hmin
    hmin = opt.hmin
    mxords = opt.mxords
    mxordn = opt.mxordn

    maxord = mxordn
    if c.meth == 2:
        maxord = mxords

    if jstart == 0:
        c.nq = 1
        c.ialth = 2
        c.rmax = 10000.0
        c.rc = 0.0
        c.crate = 0.7
        c.hold = c.h
        c.nslp = 0
        c.ipup = c.miter
        c.el[1] = 1.0
        c.icount = 20
        c.irflag = 0
        c.pdest = 0.0
        c.pdlast = 0.0
        # cfode(1) in the C, which assumes meth = 1 at the start; under
        # force_stiff the start is meth = 2 and the tables must match it
        _cfode(c, c.meth)
        _resetcoeff(c)

    if jstart == -1:
        c.ipup = c.miter
        if c.ialth == 1:
            c.ialth = 2
        if c.meth != c.mused:
            _cfode(c, c.meth)
            c.ialth = c.nq + 1
            _resetcoeff(c)
        if c.h != c.hold:
            rh = c.h / c.hold
            c.h = c.hold
            _scaleh(c, rh, neq, opt)

    if jstart == -2:
        if c.h != c.hold:
            rh = c.h / c.hold
            c.h = c.hold
            _scaleh(c, rh, neq, opt)

    dsm = 0.0
    while True:
        # Prediction loop
        c.jcur = 0
        while True:
            if abs(c.rc - 1.0) > CCMAX:
                c.ipup = c.miter
            if c.nst >= c.nslp + MSBP:
                c.ipup = c.miter
            c.tn += c.h
            # Pascal triangle multiplication on yh
            for j in range(c.nq, 0, -1):
                for i1 in range(j, c.nq + 1):
                    for i in range(1, neq + 1):
                        c.yh[i1][i] += c.yh[i1 + 1][i]
            pnorm = _vmnorm(neq, c.yh[1], c.ewt)
            c.told_corr = told
            corflag, del_val, delp, m = _correction(c, y, pnorm, neq, f, data, opt)
            if corflag == 0:
                break
            if corflag == 1:
                rh = max(0.25, hmin / abs(c.h))
                _scaleh(c, rh, neq, opt)
                continue
            if corflag == 2:
                kflag = -2
                c.hold = c.h
                return kflag

        # Local error test
        if m == 0:
            dsm = del_val / c.tesco[c.nq][2]
        if m > 0:
            dsm = _vmnorm(neq, c.acor, c.ewt) / c.tesco[c.nq][2]

        if dsm <= 1.0:
            # Successful step
            kflag = 0
            c.nst += 1
            c.hu = c.h
            c.nqu = c.nq
            c.mused = c.meth
            for j in range(1, c.nq + 2):
                r = c.el[j]
                for i in range(1, neq + 1):
                    c.yh[j][i] += r * c.acor[i]
            c.icount -= 1
            if c.icount < 0:
                rh_sw = _methodswitch(c, dsm, pnorm, neq, opt)
                if rh_sw is not None and c.meth != c.mused:
                    rh_sw = max(rh_sw, hmin / abs(c.h))
                    _scaleh(c, rh_sw, neq, opt)
                    c.rmax = 10.0
                    # endstoda
                    r = 1.0 / c.tesco[c.nqu][2]
                    for i in range(1, neq + 1):
                        c.acor[i] *= r
                    c.hold = c.h
                    break

            # No method switch -- usual step/order selection
            c.ialth -= 1
            if c.ialth == 0:
                rhup = 0.0
                if c.nq + 1 != maxord + 1:
                    for i in range(1, neq + 1):
                        c.savf[i] = c.acor[i] - c.yh[maxord + 1][i]
                    dup = _vmnorm(neq, c.savf, c.ewt) / c.tesco[c.nq][3]
                    exup = 1.0 / float(c.nq + 2)
                    rhup = 1.0 / (1.4 * dup ** exup + 0.0000014)
                orderflag, rh = _orderswitch(c, rhup, dsm, kflag, maxord, neq, opt)
                if orderflag == 0:
                    # endstoda
                    r = 1.0 / c.tesco[c.nqu][2]
                    for i in range(1, neq + 1):
                        c.acor[i] *= r
                    c.hold = c.h
                    break
                if orderflag == 1:
                    rh = max(rh, hmin / abs(c.h))
                    _scaleh(c, rh, neq, opt)
                    c.rmax = 10.0
                    # endstoda
                    r = 1.0 / c.tesco[c.nqu][2]
                    for i in range(1, neq + 1):
                        c.acor[i] *= r
                    c.hold = c.h
                    break
                if orderflag == 2:
                    _resetcoeff(c)
                    rh = max(rh, hmin / abs(c.h))
                    _scaleh(c, rh, neq, opt)
                    c.rmax = 10.0
                    # endstoda
                    r = 1.0 / c.tesco[c.nqu][2]
                    for i in range(1, neq + 1):
                        c.acor[i] *= r
                    c.hold = c.h
                    break

            if c.ialth > 1 or c.nq + 1 == maxord + 1:
                # endstoda
                r = 1.0 / c.tesco[c.nqu][2]
                for i in range(1, neq + 1):
                    c.acor[i] *= r
                c.hold = c.h
                break
            for i in range(1, neq + 1):
                c.yh[maxord + 1][i] = c.acor[i]
            # endstoda
            r = 1.0 / c.tesco[c.nqu][2]
            for i in range(1, neq + 1):
                c.acor[i] *= r
            c.hold = c.h
            break
        else:
            # Error test failed
            kflag -= 1
            c.tn = told
            for j in range(c.nq, 0, -1):
                for i1 in range(j, c.nq + 1):
                    for i in range(1, neq + 1):
                        c.yh[i1][i] -= c.yh[i1 + 1][i]
            c.rmax = 2.0
            if abs(c.h) <= hmin * 1.00001:
                kflag = -1
                c.hold = c.h
                return kflag
            if kflag > -3:
                orderflag, rh = _orderswitch(c, 0.0, dsm, kflag, maxord, neq, opt)
                if orderflag == 1 or orderflag == 0:
                    if orderflag == 0:
                        rh = min(rh, 0.2)
                    rh = max(rh, hmin / abs(c.h))
                    _scaleh(c, rh, neq, opt)
                if orderflag == 2:
                    _resetcoeff(c)
                    rh = max(rh, hmin / abs(c.h))
                    _scaleh(c, rh, neq, opt)
                continue
            else:
                # 3 or more failures
                if kflag == -10:
                    kflag = -1
                    c.hold = c.h
                    return kflag
                else:
                    rh = 0.1
                    rh = max(hmin / abs(c.h), rh)
                    c.h *= rh
                    for i in range(1, neq + 1):
                        y[i] = c.yh[1][i]
                    f(c.tn, y, c.savf, data)
                    c.nfe += 1
                    for i in range(1, neq + 1):
                        c.yh[2][i] = c.h * c.savf[i]
                    c.ipup = c.miter
                    c.ialth = 5
                    if c.nq == 1:
                        continue
                    c.nq = 1
                    _resetcoeff(c)
                    continue

    return kflag


# ---- Public interface ----

def _tol_array(tol, neq):
    """Tolerances as a 1-based array of length neq+1, index 0 unused."""
    if np.isscalar(tol):
        return np.full(neq + 1, float(tol))
    out = np.zeros(neq + 1)
    out[1:] = np.asarray(tol, dtype=float)
    return out


class LSODAStepper:
    """The C driver `lsoda()` as a stateful object: one call advances the solve.

    `lsoda()` below drives it with itask=1, integrating to each requested output
    time, which is the mode the C reference benchmarks are stated for. A
    STEPPING interface needs itask=5 instead -- one internal step, never past
    `tcrit` -- which is the mode scipy's own LSODA wrapper uses and the one
    `solver_fld` drives through `ode/native_lsoda.py`.

    Parameters mirror `lsoda()`, plus:

    tcrit : float
        The instant itask=4/5 must not step past.
    force_stiff : bool
        Start on BDF and never switch to Adams; see `LSODAOptions.force_stiff`.

    The public state is `t`, `y`, `state` (2 after a successful call, negative
    after a soft failure) and `message`. A soft failure leaves `y` at the last
    accepted step; an illegal input or a failure at the very first step raises
    `LSODAError`, which is what the C calls a hard failure.
    """

    def __init__(self, f, y0, t0, rtol=1e-6, atol=1e-6, max_steps=500,
                 mxordn=12, mxords=5, hmax=0.0, hmin=0.0, h0=0.0, tcrit=0.0,
                 force_stiff=False, data=None):
        y0 = np.asarray(y0, dtype=float).ravel()
        neq = len(y0)
        if neq < 1:
            raise LSODAError("neq = %d is less than 1" % neq)
        if hmax < 0.0 or hmin < 0.0:
            raise LSODAError("hmax and hmin must be nonnegative")
        self.neq = neq
        self.data = data
        self._user_f = f
        self._ydot_buf = np.zeros(neq)

        opt = LSODAOptions()
        opt.rtol = _tol_array(rtol, neq)
        opt.atol = _tol_array(atol, neq)
        if np.any(opt.rtol[1:] < 0.0):
            raise LSODAError("rtol = %g is less than 0." % np.min(opt.rtol[1:]))
        if np.any(opt.atol[1:] < 0.0):
            raise LSODAError("atol = %g is less than 0." % np.min(opt.atol[1:]))
        opt.mxstep = max_steps if max_steps > 0 else 500
        opt.mxhnil = 10
        opt.mxordn = min(mxordn if mxordn > 0 else 100, 12)
        opt.mxords = min(mxords if mxords > 0 else 100, 5)
        opt.h0 = h0
        opt.hmax = hmax
        opt.hmin = hmin
        opt.hmxi = 1.0 / hmax if hmax > 0 else 0.0
        opt.tcrit = tcrit
        opt.itask = 1
        opt.ixpr = 0
        opt.force_stiff = bool(force_stiff)
        self.opt = opt

        self.c = _alloc_mem(neq, opt.mxordn, opt.mxords)
        self.state = 1
        self.t = float(t0)
        self.message = ''
        # JSTART lives in the FORTRAN common block, i.e. it PERSISTS across
        # calls: dstoda leaves it at 1 and the driver overrides it with -1 when
        # a method switch has to be completed on the next step. liblsoda made it
        # a local and rebuilds it as 1 on every continuation call, which drops
        # that -1. With itask=1 over sparse output times the switch is normally
        # completed inside the same call and nothing shows; in STEPPING mode
        # (itask=2/5) every step returns, so the switch was never completed, the
        # elco tables stayed on the old method and Robertson ran away to
        # y1 = -1.9e7 in 8e6 function evaluations. Persisting it is the
        # reference's own rule, not a repair on top of it.
        self.jstart = 0
        self._y = np.zeros(neq + 1)
        self._y[1:neq + 1] = y0

    # -- read-only views on the internal state --

    @property
    def y(self):
        return self._y[1:self.neq + 1].copy()

    @property
    def h(self):
        return self.c.h

    @property
    def nfe(self):
        return self.c.nfe

    @property
    def nje(self):
        return self.c.nje

    @property
    def nst(self):
        return self.c.nst

    @property
    def meth(self):
        return self.c.meth

    def interpolate(self, t, k=0):
        """The k-th derivative at t from the Nordsieck history, i.e. `intdy`."""
        dky = np.zeros(self.neq + 1)
        iflag = _intdy(self.c, t, k, dky, self.neq)
        if iflag != 0:
            raise LSODAError("intdy refused t = %g (iflag = %d)" % (t, iflag))
        return dky[1:self.neq + 1].copy()

    def _f(self, t, y1, ydot1, data):
        self._user_f(t, y1[1:self.neq + 1], self._ydot_buf, data)
        ydot1[1:self.neq + 1] = self._ydot_buf[:]

    def advance(self, tout, itask=1):
        """One call of the C driver, integrating from t towards tout.

        Returns the state: 2 on success, negative on a soft failure (the reason
        is in `message`).
        """
        c = self.c
        opt = self.opt
        neq = self.neq
        y = self._y
        opt.itask = itask
        h0 = opt.h0
        ihit = False

        if self.state == 1 and (tout - self.t) * h0 < 0.0:
            raise LSODAError("tout = %g behind t = %g, the integration direction "
                             "is given by %g" % (tout, self.t, h0))
        if self.state == 3:
            self.jstart = -1

        if self.state == 1:
            c.meth = 2 if opt.force_stiff else 1
            if opt.force_stiff:
                c.miter = 2
            c.tn = self.t
            c.tsw = self.t
            if itask == 4 or itask == 5:
                if (opt.tcrit - tout) * (tout - self.t) < 0.0:
                    raise LSODAError("itask = 4 or 5 and tcrit behind tout")
                if h0 != 0.0 and (self.t + h0 - opt.tcrit) * h0 > 0.0:
                    h0 = opt.tcrit - self.t
            self.jstart = 0
            c.nq = 1
            self._f(self.t, y, c.yh[2], self.data)
            c.nfe = 1
            for i in range(1, neq + 1):
                c.yh[1][i] = y[i]
            for i in range(1, neq + 1):
                c.ewt[i] = opt.rtol[i] * abs(y[i]) + opt.atol[i]
                c.ewt[i] = 1.0 / c.ewt[i]
                if c.ewt[i] <= 0.0:
                    raise LSODAError("ewt[%d] = %g <= 0" % (i, c.ewt[i]))
            if h0 == 0.0:
                tdist = abs(tout - self.t)
                w0 = max(abs(self.t), abs(tout))
                if tdist < 2.0 * ETA * w0:
                    raise LSODAError("tout too close to t to start integration")
                tol_val = 0.0
                for i in range(1, neq + 1):
                    tol_val = max(tol_val, opt.rtol[i])
                if tol_val <= 0.0:
                    for i in range(1, neq + 1):
                        ayi = abs(y[i])
                        if ayi != 0.0:
                            tol_val = max(tol_val, opt.atol[i] / ayi)
                tol_val = max(tol_val, 100.0 * ETA)
                tol_val = min(tol_val, 0.001)
                sum_val = _vmnorm(neq, c.yh[2], c.ewt)
                sum_val = 1.0 / (tol_val * w0 * w0) + tol_val * sum_val * sum_val
                h0 = 1.0 / math.sqrt(sum_val)
                h0 = min(h0, tdist)
                h0 *= (1.0 if tout - self.t >= 0.0 else -1.0)
            rh = abs(h0) * opt.hmxi
            if rh > 1.0:
                h0 /= rh
            c.h = h0
            for i in range(1, neq + 1):
                c.yh[2][i] *= h0

        if self.state == 2 or self.state == 3:
            c.nslast = c.nst
            if itask == 1:
                if (c.tn - tout) * c.h >= 0.0:
                    return self._intdy_return(tout, itask)
            elif itask == 3:
                tp = c.tn - c.hu * (1.0 + 100.0 * ETA)
                if (tp - tout) * c.h > 0.0:
                    raise LSODAError("itask = %d and tout behind tcur - hu" % itask)
                if (c.tn - tout) * c.h >= 0.0:
                    return self._success_return(itask, False)
            elif itask == 4 or itask == 5:
                # case 4 falls through into case 5 in the C driver
                if itask == 4:
                    if (c.tn - opt.tcrit) * c.h > 0.0:
                        raise LSODAError("itask = 4 or 5 and tcrit behind tcur")
                    if (opt.tcrit - tout) * c.h < 0.0:
                        raise LSODAError("itask = 4 or 5 and tcrit behind tout")
                    if (c.tn - tout) * c.h >= 0.0:
                        return self._intdy_return(tout, itask)
                else:
                    if (c.tn - opt.tcrit) * c.h > 0.0:
                        raise LSODAError("itask = 4 or 5 and tcrit behind tcur")
                hmx = abs(c.tn) + abs(c.h)
                ihit = abs(c.tn - opt.tcrit) <= (100.0 * ETA * hmx)
                if ihit:
                    self.t = opt.tcrit
                    return self._success_return(itask, ihit)
                tnext = c.tn + c.h * (1.0 + 4.0 * ETA)
                if (tnext - opt.tcrit) * c.h > 0.0:
                    c.h = (opt.tcrit - c.tn) * (1.0 - 4.0 * ETA)
                    if self.state == 2:
                        self.jstart = -2
            elif itask != 2:
                raise LSODAError("illegal itask = %d" % itask)

        while True:
            if self.state != 1 or c.nst != 0:
                if (c.nst - c.nslast) >= opt.mxstep:
                    return self._soft_failure(
                        -1, "%d steps taken before reaching tout" % opt.mxstep)
                bad = 0
                for i in range(1, neq + 1):
                    c.ewt[i] = opt.rtol[i] * abs(c.yh[1][i]) + opt.atol[i]
                    c.ewt[i] = 1.0 / c.ewt[i]
                    if c.ewt[i] <= 0.0 and bad == 0:
                        bad = i
                if bad:
                    return self._soft_failure(
                        -6, "ewt[%d] = %g <= 0." % (bad, c.ewt[bad]))
            tolsf = ETA * _vmnorm(neq, c.yh[1], c.ewt)
            if tolsf > 0.01:
                tolsf *= 200.0
                if c.nst == 0:
                    raise LSODAError("at start of problem, too much accuracy requested "
                                     "for precision of machine, suggested scaling "
                                     "factor = %g" % tolsf)
                return self._soft_failure(
                    -2, "at t = %g, too much accuracy requested for precision of "
                        "machine, suggested scaling factor = %g" % (self.t, tolsf))
            if c.tn + c.h == c.tn:
                c.nhnil += 1
                if c.nhnil <= opt.mxhnil:
                    import warnings
                    warnings.warn("lsoda: internal t=%g and h=%g are such that t+h=t"
                                  % (c.tn, c.h))

            kflag = _stoda(c, y, self.jstart, neq, self._f, self.data, opt)

            if kflag == 0:
                self.jstart = 1
                if c.meth != c.mused:
                    c.tsw = c.tn
                    self.jstart = -1
                if itask == 1:
                    if (c.tn - tout) * c.h < 0.0:
                        continue
                    return self._intdy_return(tout, itask)
                if itask == 2:
                    return self._success_return(itask, ihit)
                if itask == 3:
                    if (c.tn - tout) * c.h >= 0.0:
                        return self._success_return(itask, ihit)
                    continue
                if itask == 4:
                    if (c.tn - tout) * c.h >= 0.0:
                        return self._intdy_return(tout, itask)
                    hmx = abs(c.tn) + abs(c.h)
                    ihit = abs(c.tn - opt.tcrit) <= (100.0 * ETA * hmx)
                    if ihit:
                        return self._success_return(itask, ihit)
                    tnext = c.tn + c.h * (1.0 + 4.0 * ETA)
                    if (tnext - opt.tcrit) * c.h <= 0.0:
                        continue
                    c.h = (opt.tcrit - c.tn) * (1.0 - 4.0 * ETA)
                    self.jstart = -2
                    continue
                if itask == 5:
                    hmx = abs(c.tn) + abs(c.h)
                    ihit = abs(c.tn - opt.tcrit) <= (100.0 * ETA * hmx)
                    return self._success_return(itask, ihit)

            if kflag == -1 or kflag == -2:
                big = 0.0
                c.imxer = 1
                for i in range(1, neq + 1):
                    size = abs(c.acor[i]) * c.ewt[i]
                    if big < size:
                        big = size
                        c.imxer = i
                if kflag == -1:
                    return self._soft_failure(
                        -4, "at t = %g and step size h = %g, the error test failed "
                            "repeatedly or with abs(h) = hmin" % (c.tn, c.h))
                return self._soft_failure(
                    -5, "at t = %g and step size h = %g, the corrector convergence "
                        "failed repeatedly or with abs(h) = hmin" % (c.tn, c.h))

    # -- the three exits of the C driver --

    def _success_return(self, itask, ihit):
        c = self.c
        for i in range(1, self.neq + 1):
            self._y[i] = c.yh[1][i]
        self.t = c.tn
        if (itask == 4 or itask == 5) and ihit:
            self.t = self.opt.tcrit
        self.state = 2
        return self.state

    def _intdy_return(self, tout, itask):
        c = self.c
        iflag = _intdy(c, tout, 0, self._y, self.neq)
        if iflag != 0:
            import warnings
            warnings.warn("lsoda: trouble from intdy, itask = %d, tout = %g"
                          % (itask, tout))
            for i in range(1, self.neq + 1):
                self._y[i] = c.yh[1][i]
        self.t = tout
        self.state = 2
        return self.state

    def _soft_failure(self, code, message):
        c = self.c
        for i in range(1, self.neq + 1):
            self._y[i] = c.yh[1][i]
        self.t = c.tn
        self.state = code
        self.message = message
        return self.state


def lsoda(f, y0, t_span, t_eval=None, rtol=1e-6, atol=1e-6, max_steps=500,
          h0=0.0, hmax=0.0, hmin=0.0, mxordn=12, mxords=5, data=None,
          dense_output=False, force_stiff=False):
    """Solve an ODE system using the LSODA method.

    Parameters
    ----------
    f : callable
        Right-hand side function f(t, y, ydot, data).
        Must write derivatives into ydot (1-based array, indices 1..neq).
        When called from this interface, f receives 0-based arrays and
        the wrapper handles conversion.
    y0 : array_like
        Initial state vector (0-based, length neq).
    t_span : tuple (t0, tf)
        Integration interval.
    t_eval : array_like, optional
        Times at which to store the solution. If None, only returns
        the final state.
    rtol : float or array_like
        Relative tolerance(s).
    atol : float or array_like
        Absolute tolerance(s).
    max_steps : int
        Maximum number of steps between output points.
    h0 : float
        Initial step size (0 = automatic).
    hmax : float
        Maximum step size (0 = no limit).
    hmin : float
        Minimum step size.
    mxordn : int
        Maximum Adams order (default 12).
    mxords : int
        Maximum BDF order (default 5).
    data : object
        Extra data passed to f.
    dense_output : bool
        If True, return interpolation-capable solution (not yet implemented).
    force_stiff : bool
        Start on BDF and never switch to Adams; see `LSODAOptions.force_stiff`
        and the caveat on `solver_fld.ode.native_lsoda.NativeLSODAStiff`.

    Returns
    -------
    result : LSODAResult
        Object with attributes:
        - t : array of times
        - y : array of states (shape: len(t) x neq)
        - nfe : number of function evaluations
        - nje : number of Jacobian evaluations
        - nst : number of steps
        - success : bool
        - message : str
    """
    y0 = np.asarray(y0, dtype=float)
    t0, tf = float(t_span[0]), float(t_span[1])

    stepper = LSODAStepper(f, y0, t0, rtol=rtol, atol=atol, max_steps=max_steps,
                           mxordn=mxordn, mxords=mxords, hmax=hmax, hmin=hmin,
                           h0=h0, force_stiff=force_stiff, data=data)

    if t_eval is not None:
        t_out_list = np.asarray(t_eval, dtype=float)
    else:
        t_out_list = np.array([tf])

    t_results = [t0]
    y_results = [y0.copy()]
    success = True
    message = "Integration successful"

    for tout in t_out_list:
        state = stepper.advance(float(tout), itask=1)
        t_results.append(float(tout))
        y_results.append(stepper.y)
        if state <= 0:
            success = False
            message = stepper.message
            break

    return LSODAResult(
        t=np.array(t_results),
        y=np.array(y_results),
        nfe=stepper.nfe,
        nje=stepper.nje,
        nst=stepper.nst,
        success=success,
        message=message,
    )


class LSODAResult:
    """Result of an LSODA integration."""

    def __init__(self, t, y, nfe, nje, nst, success, message):
        self.t = t
        self.y = y
        self.nfe = nfe
        self.nje = nje
        self.nst = nst
        self.success = success
        self.message = message

    def __repr__(self):
        return (f"LSODAResult(success={self.success}, nst={self.nst}, "
                f"nfe={self.nfe}, nje={self.nje}, message='{self.message}')")


def lsoda_odeint(f, y0, t_eval, rtol=1e-6, atol=1e-6, max_steps=500,
                 h0=0.0, hmax=0.0, hmin=0.0, data=None):
    """Convenience wrapper matching scipy.integrate.odeint style.

    Parameters
    ----------
    f : callable
        f(t, y) -> dydt. Takes and returns 0-based arrays.
    y0 : array_like
        Initial conditions (0-based).
    t_eval : array_like
        Output times (first element is t0).
    rtol, atol : float or array_like
        Tolerances.
    max_steps : int
        Maximum steps between output points.
    h0, hmax, hmin : float
        Step size controls.
    data : object
        Extra data passed to f.

    Returns
    -------
    result : LSODAResult
        Result with t and y arrays.
    """
    t_eval = np.asarray(t_eval, dtype=float)
    y0 = np.asarray(y0, dtype=float)
    neq = len(y0)
    _ydot = np.zeros(neq)

    def _f_bridge(t_val, y_arr, ydot_arr, udata):
        result = f(t_val, y_arr)
        ydot_arr[:] = result

    return lsoda(_f_bridge, y0, (t_eval[0], t_eval[-1]),
                 t_eval=t_eval[1:], rtol=rtol, atol=atol,
                 max_steps=max_steps, h0=h0, hmax=hmax, hmin=hmin,
                 data=data)
