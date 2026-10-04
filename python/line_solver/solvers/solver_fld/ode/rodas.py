"""
RODAS -- Rosenbrock method of order (3)4 for stiff and differential-algebraic
systems  M y' = f(x,y),  including a SINGULAR mass matrix (index-1 DAE).

PORTED THIRD-PARTY CODE -- do not edit to fix a LINE bug; fix the caller.

Source: E. Hairer and G. Wanner, rodas.f / dc_decsol.f / decsol.f, version of
October 28, 1996, as published with "Solving Ordinary Differential Equations
II. Stiff and Differential-Algebraic Problems", Springer Series in
Computational Mathematics 14.

WHY THIS EXISTS IN PYTHON. scipy has no mass matrix at all: `solve_ivp` solves
y' = f, and `BDF`/`Radau` take no M. An index-1 DAE with a SINGULAR M -- which
is what the fluid `dae` method integrates, one algebraic row per closed chain --
therefore has no scipy route, and the native Python implementation may not
depend on a Fortran/C extension. This is the same solver `cpp/third_party/
rodas.hpp` vendors as an f2c translation and `matlab/lib/thirdparty/rodas`
carries, so all four codebases run one step sequence.

THE TRANSLATION IS INDEX-FOR-INDEX. Every array here is allocated one longer
than it needs to be and addressed from 1, exactly as the Fortran does, and the
loop bounds are the Fortran's own. That is deliberate: a 0-based rewrite of
2000 lines of Fortran is where an off-by-one hides, and it would be invisible
until it changed an answer nobody has a reference for. Index 0 of every array
is unused.

WHERE NUMPY REPLACES A LOOP IT IS ELEMENTWISE, NEVER A REDUCTION. A vectorised
`a[k+1:n+1, j] += a[k+1:n+1, k] * t` performs the same multiply and the same
add on each element in the same order as the Fortran DO loop, so it is bit
identical. A vectorised SUM is not: numpy pairwise-sums, Fortran accumulates
sequentially, and the two differ in the last bits. Every reduction in this file
is therefore still a Python loop -- the error norm, the inner sums of DECOMR
and SLVROD -- and that is not an oversight to be optimised away.

Reachable IJOB values are 1..5 (and 11..15 when M1 > 0, the second-order form);
RODAS itself never selects 6 or 7, which belong to RADAU5's Hessenberg option
and are absent here rather than transliterated dead.

NOT REENTRANT ACROSS A SINGLE `RodasDense`. The dense-output state travels with
the result object rather than in a module global, so concurrent integrations
are fine; it is the caller's `cont` array that must not be shared.
"""

import math

import numpy as np

__all__ = ["rodas", "RodasResult", "RodasDense", "RodasError"]


class RodasError(RuntimeError):
    """Raised for an input RODAS itself rejects before integrating."""


# IDID values, as the Fortran sets them.
IDID_SUCCESS = 1
IDID_SOLOUT_STOP = 2
IDID_NMAX = -2
IDID_STEP_TOO_SMALL = -3
IDID_SINGULAR = -4


class _Linal(object):
    """The /LINAL/ COMMON block: band offsets shared by DECOMR and SLVROD."""

    __slots__ = ("mle", "mue", "mbjac", "mbb", "mdiag", "mdiff", "mbdiag")

    def __init__(self):
        self.mle = 0
        self.mue = 0
        self.mbjac = 0
        self.mbb = 0
        self.mdiag = 0
        self.mdiff = 0
        self.mbdiag = 0


class RodasDense(object):
    """
    The /CONROS/ COMMON block, and the third-order interpolant it feeds.

    A Rosenbrock method chooses its step from the local error, so its accepted
    points are wherever the stiffness put them and never the ones a caller asked
    for. This is the interpolant RODAS carries for exactly that, valid over the
    step just accepted, so an arbitrary output grid costs no extra step.
    """

    __slots__ = ("cont", "xold", "h", "n")

    def __init__(self, cont, n):
        self.cont = cont
        self.xold = 0.0
        self.h = 0.0
        self.n = n

    def value(self, i, x):
        """CONTRO: component `i` (0-BASED, as a Python caller expects) at `x`."""
        cont = self.cont
        n = self.n
        ii = i + 1
        s = (x - self.xold) / self.h
        return cont[ii] * (1 - s) + s * (
            cont[ii + n] + (1 - s) * (cont[ii + 2 * n] + s * cont[ii + 3 * n])
        )


class RodasResult(object):
    """What the integration ended at, and what it cost."""

    __slots__ = ("x", "y", "h", "idid", "nfcn", "njac", "nstep",
                 "naccpt", "nrejct", "ndec", "nsol")

    def __init__(self):
        self.x = 0.0
        self.y = None
        self.h = 0.0
        self.idid = 0
        self.nfcn = 0
        self.njac = 0
        self.nstep = 0
        self.naccpt = 0
        self.nrejct = 0
        self.ndec = 0
        self.nsol = 0

    @property
    def success(self):
        return self.idid > 0

    def __repr__(self):
        return ("RodasResult(idid=%d, x=%r, nfcn=%d, njac=%d, nstep=%d, "
                "naccpt=%d, nrejct=%d)" % (self.idid, self.x, self.nfcn,
                                           self.njac, self.nstep,
                                           self.naccpt, self.nrejct))


# ---------------------------------------------------------------------------
# decsol.f: the real factorisations RODAS uses. DECC/SOLC/DECHC/SOLBC belong to
# RADAU5's complex pair and DECH/SOLH to its Hessenberg option; RODAS reaches
# none of them, so they are absent rather than carried as dead code.
# ---------------------------------------------------------------------------

def _dec(n, a, ip):
    """
    DEC: LU by Gaussian elimination with partial pivoting, Moler's ACM 423.

    `a` is (ndim+1, n+1), 1-based. On return the strict upper triangle plus the
    diagonal is U and the strict lower triangle holds the NEGATED multipliers,
    which is what SOL expects. ip[n] carries (-1)**(interchanges), or 0 when the
    matrix was found singular.

    Returns ier: 0 if nonsingular, else the stage k at which the pivot vanished.
    """
    ier = 0
    ip[n] = 1
    if n != 1:
        nm1 = n - 1
        for k in range(1, nm1 + 1):
            kp1 = k + 1
            # The Fortran scans i = k+1..n replacing on STRICTLY greater, so it
            # keeps the first index attaining the maximum, k included. argmax
            # keeps the first too.
            m = k + int(np.argmax(np.abs(a[k:n + 1, k])))
            ip[k] = m
            t = a[m, k]
            if m != k:
                ip[n] = -ip[n]
                a[m, k] = a[k, k]
                a[k, k] = t
            if t == 0.0:
                ip[n] = 0
                return k
            t = 1.0 / t
            a[kp1:n + 1, k] = -a[kp1:n + 1, k] * t
            for j in range(kp1, n + 1):
                t = a[m, j]
                a[m, j] = a[k, j]
                a[k, j] = t
                if t != 0.0:
                    a[kp1:n + 1, j] += a[kp1:n + 1, k] * t
    if a[n, n] == 0.0:
        ip[n] = 0
        return n
    return ier


def _sol(n, a, b, ip):
    """SOL: forward/back substitution against the factors DEC left in `a`."""
    if n != 1:
        nm1 = n - 1
        for k in range(1, nm1 + 1):
            kp1 = k + 1
            m = ip[k]
            t = b[m]
            b[m] = b[k]
            b[k] = t
            b[kp1:n + 1] += a[kp1:n + 1, k] * t
        for kb in range(1, nm1 + 1):
            km1 = n - kb
            k = km1 + 1
            b[k] /= a[k, k]
            t = -b[k]
            b[1:km1 + 1] += a[1:km1 + 1, k] * t
    b[1] /= a[1, 1]


def _decb(n, a, ml, mu, ip):
    """
    DECB: the banded counterpart of DEC. The matrix arrives in LINPACK band
    storage, its diagonals in rows ML+1 .. 2*ML+MU+1 of `a`.
    """
    ier = 0
    ip[n] = 1
    md = ml + mu + 1
    md1 = md + 1
    ju = 0
    if ml != 0 and n != 1:
        if n >= mu + 2:
            for j in range(mu + 2, n + 1):
                a[1:ml + 1, j] = 0.0
        nm1 = n - 1
        for k in range(1, nm1 + 1):
            kp1 = k + 1
            m = md
            mdl = min(ml, n - k) + md
            for i in range(md1, mdl + 1):
                if abs(a[i, k]) > abs(a[m, k]):
                    m = i
            ip[k] = m + k - md
            t = a[m, k]
            if m != md:
                ip[n] = -ip[n]
                a[m, k] = a[md, k]
                a[md, k] = t
            if t == 0.0:
                ip[n] = 0
                return k
            t = 1.0 / t
            a[md1:mdl + 1, k] = -a[md1:mdl + 1, k] * t
            ju = min(max(ju, mu + ip[k]), n)
            mm = md
            if ju >= kp1:
                for j in range(kp1, ju + 1):
                    m -= 1
                    mm -= 1
                    t = a[m, j]
                    if m != mm:
                        a[m, j] = a[mm, j]
                        a[mm, j] = t
                    if t != 0.0:
                        jk = j - k
                        # Left as the Fortran's loop: the vectorised form needs
                        # md1-jk > 0, which holds by the band structure, but a
                        # negative start index would wrap in numpy rather than
                        # raise, so the guard is worth more than the speed.
                        for i in range(md1, mdl + 1):
                            a[i - jk, j] += a[i, k] * t
    if a[md, n] == 0.0:
        ip[n] = 0
        return n
    return ier


def _solb(n, a, ml, mu, b, ip):
    """SOLB: forward/back substitution against the band factors of DECB."""
    md = ml + mu + 1
    md1 = md + 1
    mdm = md - 1
    nm1 = n - 1
    if ml != 0:
        if n == 1:
            b[1] /= a[md, 1]
            return
        for k in range(1, nm1 + 1):
            m = ip[k]
            t = b[m]
            b[m] = b[k]
            b[k] = t
            mdl = min(ml, n - k) + md
            b[md1 + k - md:mdl + 1 + k - md] += a[md1:mdl + 1, k] * t
    for kb in range(1, nm1 + 1):
        k = n + 1 - kb
        b[k] /= a[md, k]
        t = -b[k]
        kmd = md - k
        lm = max(1, kmd + 1)
        if lm <= mdm:
            b[lm - kmd:mdm + 1 - kmd] += a[lm:mdm + 1, k] * t
    b[1] /= a[md, 1]


# ---------------------------------------------------------------------------
# rodas.f: the method coefficients
# ---------------------------------------------------------------------------

def _rocoe(meth):
    """
    ROCOE: the tableau. METH 1 is the default; 2 is Hairer-Wanner's second set
    and 3 is Steinebach's (1993) order-4 set for linear parabolic problems.
    The BET2P..BET4P of the original are computed and unused, exactly as there.
    """
    c = {}
    if meth == 1:
        c.update(c2=.386, c3=.21, c4=.63,
                 d1=.25, d2=-.1043, d3=.1035, d4=-.03620000000000023,
                 a21=1.544, a31=.9466785280815826, a32=.2557011698983284,
                 a41=3.314825187068521, a42=2.896124015972201,
                 a43=.9986419139977817, a51=1.221224509226641,
                 a52=6.019134481288629, a53=12.53708332932087,
                 a54=-.687886036105895,
                 c21=-5.6688, c31=-2.430093356833875, c32=-.2063599157091915,
                 c41=-.1073529058151375, c42=-9.594562251023355,
                 c43=-20.47028614809616, c51=7.496443313967647,
                 c52=-10.24680431464352, c53=-33.99990352819905,
                 c54=11.7089089320616, c61=8.083246795921522,
                 c62=-7.981132988064893, c63=-31.52159432874371,
                 c64=16.31930543123136, c65=-6.058818238834054,
                 gamma=.25,
                 d21=10.12623508344586, d22=-7.487995877610167,
                 d23=-34.80091861555747, d24=-7.992771707568823,
                 d25=1.025137723295662, d31=-.6762803392801253,
                 d32=6.087714651680015, d33=16.43084320892478,
                 d34=24.76722511418386, d35=-6.594389125716872)
    elif meth == 2:
        c.update(c2=.3507221, c3=.2557041, c4=.681779,
                 d1=.25, d2=-.06902209999999998, d3=-9.671999999999459e-4,
                 d4=-.08797900000000025,
                 a21=1.4028884, a31=.6581212688557198, a32=-1.320936088384301,
                 a41=7.131197445744498, a42=16.02964143958207,
                 a43=-5.561572550509766, a51=22.73885722420363,
                 a52=67.38147284535289, a53=-31.2187749303856,
                 a54=.7285641833203814,
                 c21=-5.1043536, c31=-2.899967805418783, c32=4.040399359702244,
                 c41=-32.64449927841361, c42=-99.35311008728094,
                 c43=49.99119122405989, c51=-76.46023087151691,
                 c52=-278.5942120829058, c53=153.9294840910643,
                 c54=10.97101866258358, c61=-76.29701586804983,
                 c62=-294.2795630511232, c63=162.0029695867566,
                 c64=23.6516690309527, c65=-7.652977706771382,
                 gamma=.25,
                 d21=-38.71940424117216, d22=-135.8025833007622,
                 d23=64.51068857505875, d24=-4.192663174613162,
                 d25=-2.53193205033506, d31=-14.99268484949843,
                 d32=-76.30242396627033, d33=58.65928432851416,
                 d34=16.61359034616402, d35=-.6758691794084156)
    elif meth == 3:
        gamma = .25
        c.update(gamma=gamma, c2=gamma * 3., c3=.21, c4=.63,
                 d1=.25, d2=-.5, d3=-.023504, d4=-.0362,
                 a21=3., a31=1.831036793486759, a32=.4955183967433795,
                 a41=2.304376582692669, a42=-.05249275245743001,
                 a43=-1.176798761832782, a51=-7.170454962423024,
                 a52=-4.741636671481785, a53=-16.31002631330971,
                 a54=-1.062004044111401,
                 c21=-12., c31=-8.791795173947035, c32=-2.207865586973518,
                 c41=10.81793056857153, c42=6.780270611428266,
                 c43=19.5348594464241, c51=34.19095006749676,
                 c52=15.49671153725963, c53=54.7476087596413,
                 c54=14.16005392148534, c61=34.62605830930532,
                 c62=15.30084976114473, c63=56.99955578662667,
                 c64=18.40807009793095, c65=-5.714285714285717,
                 d21=25.09876703708589, d22=11.62013104361867,
                 d23=28.49148307714626, d24=-5.664021568594133, d25=0.,
                 d31=1.638054557396973, d32=-.7373619806678748,
                 d33=8.47791821923899, d34=15.9925314877952,
                 d35=-1.882352941176471)
    else:
        raise RodasError("rodas: CURIOUS INPUT IWORK(2)=%d" % meth)
    return c


# ---------------------------------------------------------------------------
# dc_decsol.f: build and factor E = fac1*M - J, and solve against it
# ---------------------------------------------------------------------------

def _decomr(n, fjac, fmas, mlmas, mumas, m1, m2, nm1, fac1, e1, ip1, ijob,
            lin):
    """DECOMR: assemble E1 for this IJOB and factor it. Returns ier."""
    if ijob in (1, 11):
        if ijob == 1:
            for j in range(1, n + 1):
                e1[1:n + 1, j] = -fjac[1:n + 1, j]
                e1[j, j] += fac1
            return _dec(n, e1, ip1)
        for j in range(1, nm1 + 1):
            jm1 = j + m1
            e1[1:nm1 + 1, j] = -fjac[1:nm1 + 1, jm1]
            e1[j, j] += fac1
        return _decomr_l45(fjac, m1, m2, nm1, fac1, e1, ip1)

    if ijob in (2, 12):
        if ijob == 2:
            for j in range(1, n + 1):
                e1[1 + lin.mle:lin.mbjac + lin.mle + 1, j] = \
                    -fjac[1:lin.mbjac + 1, j]
                e1[lin.mdiag, j] += fac1
            return _decb(n, e1, lin.mle, lin.mue, ip1)
        for j in range(1, nm1 + 1):
            jm1 = j + m1
            e1[1 + lin.mle:lin.mbjac + lin.mle + 1, j] = \
                -fjac[1:lin.mbjac + 1, jm1]
            e1[lin.mdiag, j] += fac1
        return _decomr_l46(fjac, m1, m2, nm1, fac1, e1, ip1, lin)

    if ijob in (3, 13):
        nn = n if ijob == 3 else nm1
        for j in range(1, nn + 1):
            jm1 = j if ijob == 3 else j + m1
            e1[1:nn + 1, j] = -fjac[1:nn + 1, jm1]
            for i in range(max(1, j - mumas), min(nn, j + mlmas) + 1):
                e1[i, j] += fac1 * fmas[i - j + lin.mbdiag, j]
        if ijob == 3:
            return _dec(n, e1, ip1)
        return _decomr_l45(fjac, m1, m2, nm1, fac1, e1, ip1)

    if ijob in (4, 14):
        nn = n if ijob == 4 else nm1
        for j in range(1, nn + 1):
            jm1 = j if ijob == 4 else j + m1
            e1[1 + lin.mle:lin.mbjac + lin.mle + 1, j] = \
                -fjac[1:lin.mbjac + 1, jm1]
            for i in range(1, lin.mbb + 1):
                ib = i + lin.mdiff
                e1[ib, j] += fac1 * fmas[i, j]
        if ijob == 4:
            return _decb(n, e1, lin.mle, lin.mue, ip1)
        return _decomr_l46(fjac, m1, m2, nm1, fac1, e1, ip1, lin)

    if ijob in (5, 15):
        nn = n if ijob == 5 else nm1
        for j in range(1, nn + 1):
            jm1 = j if ijob == 5 else j + m1
            e1[1:nn + 1, j] = fmas[1:nn + 1, j] * fac1 - fjac[1:nn + 1, jm1]
        if ijob == 5:
            return _dec(n, e1, ip1)
        return _decomr_l45(fjac, m1, m2, nm1, fac1, e1, ip1)

    # 6 is "THIS OPTION IS NOT PROVIDED" upstream; 7..10 belong to RADAU5.
    return 0


def _decomr_l45(fjac, m1, m2, nm1, fac1, e1, ip1):
    """DECOMR label 45: fold the second-order block into E1, then factor full."""
    mm = m1 // m2
    for j in range(1, m2 + 1):
        for i in range(1, nm1 + 1):
            s = 0.0
            for k in range(0, mm):
                s = (s + fjac[i, j + k * m2]) / fac1
            e1[i, j] -= s
    return _dec(nm1, e1, ip1)


def _decomr_l46(fjac, m1, m2, nm1, fac1, e1, ip1, lin):
    """DECOMR label 46: the same fold, banded."""
    mm = m1 // m2
    for j in range(1, m2 + 1):
        for i in range(1, lin.mbjac + 1):
            s = 0.0
            for k in range(0, mm):
                s = (s + fjac[i, j + k * m2]) / fac1
            e1[i + lin.mle, j] -= s
    return _decb(nm1, e1, lin.mle, lin.mue, ip1)


def _slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, m2, nm1, fac1,
            e, ip, dy, ak, fx, ynew, hd, ijob, stage1, lin):
    """SLVROD: one stage. `ak` receives the solve; `dy`, `fx`, `ynew` are read."""
    if hd == 0.0:
        ak[1:n + 1] = dy[1:n + 1]
    else:
        ak[1:n + 1] = dy[1:n + 1] + hd * fx[1:n + 1]

    if ijob in (1, 11):
        if stage1:
            ak[1:n + 1] += ynew[1:n + 1]
        if ijob == 1:
            _sol(n, e, ak, ip)
            return
        _slvrod_l48(fjac, m1, m2, nm1, fac1, e, ip, ak)
        return

    if ijob in (2, 12):
        if stage1:
            ak[1:n + 1] += ynew[1:n + 1]
        if ijob == 2:
            _solb(n, e, lin.mle, lin.mue, ak, ip)
            return
        _slvrod_l45(fjac, mljac, mujac, m1, m2, nm1, fac1, e, ip, ak, lin)
        return

    if ijob == 3:
        if stage1:
            for i in range(1, n + 1):
                s = 0.0
                for j in range(max(1, i - mlmas), min(n, i + mumas) + 1):
                    s += fmas[i - j + lin.mbdiag, j] * ynew[j]
                ak[i] += s
        _sol(n, e, ak, ip)
        return

    if ijob in (13, 14):
        if stage1:
            ak[1:m1 + 1] += ynew[1:m1 + 1]
            for i in range(1, nm1 + 1):
                s = 0.0
                for j in range(max(1, i - mlmas), min(nm1, i + mumas) + 1):
                    s += fmas[i - j + lin.mbdiag, j] * ynew[j + m1]
                ak[i + m1] += s
        if ijob == 14:
            _slvrod_l45(fjac, mljac, mujac, m1, m2, nm1, fac1, e, ip, ak, lin)
        else:
            _slvrod_l48(fjac, m1, m2, nm1, fac1, e, ip, ak)
        return

    if ijob == 4:
        if stage1:
            for i in range(1, n + 1):
                s = 0.0
                for j in range(max(1, i - mlmas), min(n, i + mumas) + 1):
                    s += fmas[i - j + lin.mbdiag, j] * ynew[j]
                ak[i] += s
        _solb(n, e, lin.mle, lin.mue, ak, ip)
        return

    if ijob == 5:
        if stage1:
            for i in range(1, n + 1):
                s = 0.0
                for j in range(1, n + 1):
                    s += fmas[i, j] * ynew[j]
                ak[i] += s
        _sol(n, e, ak, ip)
        return

    if ijob == 15:
        if stage1:
            ak[1:m1 + 1] += ynew[1:m1 + 1]
            for i in range(1, nm1 + 1):
                s = 0.0
                for j in range(1, nm1 + 1):
                    s += fmas[i, j] * ynew[j + m1]
                ak[i + m1] += s
        _slvrod_l48(fjac, m1, m2, nm1, fac1, e, ip, ak)
        return

    if ijob == 6:
        # "THIS OPTION IS NOT PROVIDED" upstream, and it solves only under
        # stage1 there. Kept identical rather than tidied.
        if stage1:
            for i in range(1, n + 1):
                s = 0.0
                for j in range(1, n + 1):
                    s += fmas[i, j] * ynew[j]
                ak[i] += s
            _solb(n, e, lin.mle, lin.mue, ak, ip)
        return


def _slvrod_l48(fjac, m1, m2, nm1, fac1, e, ip, ak):
    """SLVROD label 48: the second-order elimination, full Jacobian."""
    mm = m1 // m2
    for j in range(1, m2 + 1):
        s = 0.0
        for k in range(mm - 1, -1, -1):
            jkm = j + k * m2
            s = (ak[jkm] + s) / fac1
            ak[1 + m1:nm1 + m1 + 1] += fjac[1:nm1 + 1, jkm] * s
    _sol(nm1, e, ak[m1:], ip)
    for i in range(m1, 0, -1):
        ak[i] = (ak[i] + ak[m2 + i]) / fac1


def _slvrod_l45(fjac, mljac, mujac, m1, m2, nm1, fac1, e, ip, ak, lin):
    """SLVROD label 45: the second-order elimination, banded Jacobian."""
    mm = m1 // m2
    for j in range(1, m2 + 1):
        s = 0.0
        for k in range(mm - 1, -1, -1):
            jkm = j + k * m2
            s = (ak[jkm] + s) / fac1
            for i in range(max(1, j - mujac), min(nm1, j + mljac) + 1):
                ak[i + m1] += fjac[i + mujac + 1 - j, jkm] * s
    _solb(nm1, e, lin.mle, lin.mue, ak[m1:], ip)
    for i in range(m1, 0, -1):
        ak[i] = (ak[i] + ak[m2 + i]) / fac1


# ---------------------------------------------------------------------------
# rodas.f: the driver and the core integrator
# ---------------------------------------------------------------------------

def rodas(n, fcn, x, y, xend, h=1e-6, rtol=1e-6, atol=1e-6, itol=0,
          jac=None, ijac=0, mljac=None, mujac=0,
          dfx=None, idfx=0, ifcn=0,
          mas=None, imas=0, mlmas=0, mumas=0,
          solout=None, iout=0,
          nmax=100000, meth=1, pred=True, uround=1e-16, hmax=None,
          fac1=None, fac2=None, safe=0.9, m1=0, m2=0):
    """
    Integrate M y' = f(x,y) from `x` to `xend`, RODAS order (3)4.

    The signature is the Fortran's, with the workspace arithmetic replaced by
    ordinary allocation and the SWITCHES kept: a caller that knows rodas.f knows
    this. Callbacks take and fill 0-BASED numpy views, which is the only place
    this deviates from the original's conventions and is the boundary a Python
    caller actually touches:

        fcn(x, y, f)                fill f (length n)
        jac(x, y, dfy)              fill dfy (ldjac x n), only if ijac == 1
        dfx(x, y, fx)               fill fx (length n), only if idfx == 1
        mas(am)                     fill am (ldmas x nm1), only if imas == 1
        solout(nr, xold, x, y, dense) -> int, negative to stop

    `mljac` defaults to n, i.e. a full Jacobian. `hmax` defaults to xend - x,
    `fac1`/`fac2` to the Fortran's 5 and 1/6 written as their reciprocals.

    Returns a RodasResult; `y` is modified in place as the Fortran modifies it.
    """
    if mljac is None:
        mljac = n
    if nmax <= 0:
        raise RodasError("rodas: WRONG INPUT IWORK(1)=%d" % nmax)
    if meth <= 0 or meth >= 4:
        raise RodasError("rodas: CURIOUS INPUT IWORK(2)=%d" % meth)

    nm1 = n - m1
    if m1 == 0:
        m2 = n
    if m2 == 0:
        m2 = m1
    if m1 < 0 or m2 < 0 or m1 + m2 > n:
        raise RodasError("rodas: CURIOUS INPUT FOR IWORK(9,10)=%d %d" % (m1, m2))
    if uround < 1e-16 or uround >= 1.0:
        raise RodasError("rodas: COEFFICIENTS HAVE 16 DIGITS, UROUND=%r" % uround)
    if hmax is None:
        hmax = xend - x
    fac1 = 5.0 if fac1 is None else 1.0 / fac1
    fac2 = .16666666666666666 if fac2 is None else 1.0 / fac2
    if fac1 < 1.0 or fac2 > 1.0:
        raise RodasError("rodas: CURIOUS INPUT WORK(3,4)")
    if safe <= .001 or safe >= 1.0:
        raise RodasError("rodas: CURIOUS INPUT FOR WORK(5)=%r" % safe)

    rtolv = np.zeros(n + 1)
    atolv = np.zeros(n + 1)
    if itol == 0:
        rtolv[1] = float(rtol)
        atolv[1] = float(atol)
        if atolv[1] <= 0.0 or rtolv[1] <= uround * 10.0:
            raise RodasError("rodas: TOLERANCES ARE TOO SMALL")
    else:
        rtolv[1:n + 1] = np.asarray(rtol, dtype=float).reshape(n)
        atolv[1:n + 1] = np.asarray(atol, dtype=float).reshape(n)
        for i in range(1, n + 1):
            if atolv[i] <= 0.0 or rtolv[i] <= uround * 10.0:
                raise RodasError("rodas: TOLERANCES(%d) ARE TOO SMALL" % i)

    autnms = (ifcn == 0)
    implct = (imas != 0)
    jband = (mljac < nm1)

    if jband:
        ldjac = mljac + mujac + 1
        lde = mljac + ldjac
    else:
        mljac = nm1
        mujac = nm1
        ldjac = nm1
        lde = nm1

    if implct:
        if mlmas != nm1:
            ldmas = mlmas + mumas + 1
            ijob = 4 if jband else 3
        else:
            ldmas = nm1
            ijob = 5
        if mlmas > mljac or mumas > mujac:
            raise RodasError(
                'rodas: BANDWITH OF "MAS" NOT LARGER THAN BANDWITH OF "JAC"')
    else:
        ldmas = 0
        ijob = 2 if jband else 1
    ldmas2 = max(1, ldmas)

    res = _roscor(n, fcn, x, y, xend, hmax, h, rtolv, atolv, itol,
                  jac, ijac, mljac, mujac, dfx, idfx, mas,
                  mlmas, mumas, solout, iout, nmax, uround, meth, ijob,
                  fac1, fac2, safe, autnms, implct, jband, pred,
                  ldjac, lde, ldmas2, m1, m2, nm1)
    return res


def _roscor(n, fcn, x, y0, xend, hmax, h, rtol, atol, itol,
            jac, ijac, mljac, mujac, dfx, idfx, mas,
            mlmas, mumas, solout, iout, nmax, uround, meth, ijob,
            fac1, fac2, safe, autnms, implct, jband, pred,
            ldjac, lde, ldmas, m1, m2, nm1):
    """ROSCOR: the core integrator. `y0` is the caller's array, modified."""
    res = RodasResult()

    y = np.zeros(n + 1)
    y[1:n + 1] = np.asarray(y0, dtype=float).reshape(n)

    ynew = np.zeros(n + 1)
    dy1 = np.zeros(n + 1)
    dy = np.zeros(n + 1)
    ak1 = np.zeros(n + 1)
    ak2 = np.zeros(n + 1)
    ak3 = np.zeros(n + 1)
    ak4 = np.zeros(n + 1)
    ak5 = np.zeros(n + 1)
    ak6 = np.zeros(n + 1)
    fx = np.zeros(n + 1)
    cont = np.zeros(4 * n + 1)
    fjac = np.zeros((ldjac + 1, n + 1))
    e = np.zeros((lde + 1, nm1 + 1))
    fmas = np.zeros((ldmas + 1, nm1 + 1))
    ip = np.zeros(nm1 + 1, dtype=np.int64)

    lin = _Linal()
    dense = RodasDense(cont, n)
    nn2 = 2 * n
    nn3 = 3 * n
    lrc = 4 * n

    yview = y[1:]
    ynewview = ynew[1:]   # the stages evaluate FCN at YNEW, not at Y

    if implct:
        mas(fmas[1:, 1:])

    c = _rocoe(meth)
    a21 = c["a21"]; a31 = c["a31"]; a32 = c["a32"]
    a41 = c["a41"]; a42 = c["a42"]; a43 = c["a43"]
    a51 = c["a51"]; a52 = c["a52"]; a53 = c["a53"]; a54 = c["a54"]
    c21 = c["c21"]; c31 = c["c31"]; c32 = c["c32"]
    c41 = c["c41"]; c42 = c["c42"]; c43 = c["c43"]
    c51 = c["c51"]; c52 = c["c52"]; c53 = c["c53"]; c54 = c["c54"]
    c61 = c["c61"]; c62 = c["c62"]; c63 = c["c63"]; c64 = c["c64"]; c65 = c["c65"]
    gamma = c["gamma"]
    c2 = c["c2"]; c3 = c["c3"]; c4 = c["c4"]
    d1 = c["d1"]; d2 = c["d2"]; d3 = c["d3"]; d4 = c["d4"]
    d21 = c["d21"]; d22 = c["d22"]; d23 = c["d23"]; d24 = c["d24"]; d25 = c["d25"]
    d31 = c["d31"]; d32 = c["d32"]; d33 = c["d33"]; d34 = c["d34"]; d35 = c["d35"]

    if m1 > 0:
        ijob += 10

    posneg = math.copysign(1.0, xend - x)
    hmaxn = min(abs(hmax), abs(xend - x))
    if abs(h) <= uround * 10.0:
        h = 1e-6
    h = min(abs(h), hmaxn)
    h = math.copysign(h, posneg)
    reject = False
    last = False
    nsing = 0
    hd1 = hd2 = hd3 = hd4 = 0.0
    hacc = 0.0
    erracc = 0.0
    hopt = h

    lin.mbdiag = mumas + 1
    if jband:
        lin.mle = mljac
        lin.mue = mujac
        lin.mbjac = mljac + mujac + 1
        lin.mbb = mlmas + mumas + 1
        lin.mdiag = lin.mle + lin.mue + 1
        lin.mdiff = lin.mle + lin.mue - mumas

    nfcn = njac = nstep = naccpt = nrejct = ndec = nsol = 0

    def _finish(idid):
        res.x = x
        res.y = y[1:n + 1].copy()
        res.h = h
        res.idid = idid
        res.nfcn = nfcn
        res.njac = njac
        res.nstep = nstep
        res.naccpt = naccpt
        res.nrejct = nrejct
        res.ndec = ndec
        res.nsol = nsol
        if isinstance(y0, np.ndarray):
            y0[...] = y[1:n + 1]
        return res

    if iout != 0:
        dense.xold = x
        dense.h = h
        irtrn = solout(naccpt + 1, dense.xold, x, yview, dense)
        if irtrn is not None and irtrn < 0:
            return _finish(IDID_SOLOUT_STOP)

    # --- BASIC INTEGRATION STEP (label 1) ---
    while True:
        if nstep > nmax:
            return _finish(IDID_NMAX)
        if abs(h) * .1 <= abs(x) * uround:
            return _finish(IDID_STEP_TOO_SMALL)
        if last:
            h = hopt
            return _finish(IDID_SUCCESS)
        hopt = h
        if (x + h * 1.0001 - xend) * posneg >= 0.0:
            h = xend - x
            last = True

        # --- COMPUTATION OF THE JACOBIAN ---
        fcn(x, yview, dy1[1:])
        nfcn += 1
        njac += 1
        if ijac == 0:
            if jband:
                mujacp = mujac + 1
                md = min(lin.mbjac, n)
                for mm in range(1, m1 // m2 + 2):
                    for k in range(1, md + 1):
                        j = k + (mm - 1) * m2
                        while True:
                            ak2[j] = y[j]
                            ak3[j] = math.sqrt(uround * max(1e-5, abs(y[j])))
                            y[j] += ak3[j]
                            j += md
                            if j > mm * m2:
                                break
                        fcn(x, yview, ak1[1:])
                        j = k + (mm - 1) * m2
                        j1 = k
                        lbeg = max(1, j1 - mujac) + m1
                        while True:
                            lend = min(m2, j1 + mljac) + m1
                            y[j] = ak2[j]
                            mujacj = mujacp - j1 - m1
                            for l in range(lbeg, lend + 1):
                                fjac[l + mujacj, j] = (ak1[l] - dy1[l]) / ak3[j]
                            j += md
                            j1 += md
                            lbeg = lend + 1
                            if j > mm * m2:
                                break
            else:
                for i in range(1, n + 1):
                    ysafe = y[i]
                    delt = math.sqrt(uround * max(1e-5, abs(ysafe)))
                    y[i] = ysafe + delt
                    fcn(x, yview, ak1[1:])
                    fjac[1:n + 1 - m1, i] = \
                        (ak1[m1 + 1:n + 1] - dy1[m1 + 1:n + 1]) / delt
                    y[i] = ysafe
        else:
            jac(x, yview, fjac[1:, 1:])

        if not autnms:
            if idfx == 0:
                delt = math.sqrt(uround * max(1e-5, abs(x)))
                xdelt = x + delt
                fcn(xdelt, yview, ak1[1:])
                fx[1:n + 1] = (ak1[1:n + 1] - dy1[1:n + 1]) / delt
            else:
                dfx(x, yview, fx[1:])

        # --- COMPUTE THE STAGES (label 2) ---
        while True:
            fac = 1.0 / (h * gamma)
            ier = _decomr(n, fjac, fmas, mlmas, mumas, m1, m2, nm1, fac,
                          e, ip, ijob, lin)
            if ier != 0:
                # --- SINGULAR MATRIX (label 80) ---
                nsing += 1
                if nsing >= 5:
                    return _finish(IDID_SINGULAR)
                h *= .5
                reject = True
                last = False
                continue
            ndec += 1

            hc21 = c21 / h; hc31 = c31 / h; hc32 = c32 / h
            hc41 = c41 / h; hc42 = c42 / h; hc43 = c43 / h
            hc51 = c51 / h; hc52 = c52 / h; hc53 = c53 / h; hc54 = c54 / h
            hc61 = c61 / h; hc62 = c62 / h; hc63 = c63 / h
            hc64 = c64 / h; hc65 = c65 / h
            if not autnms:
                hd1 = h * d1
                hd2 = h * d2
                hd3 = h * d3
                hd4 = h * d4

            # --- THE STAGES ---
            _slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, m2, nm1,
                    fac, e, ip, dy1, ak1, fx, ynew, hd1, ijob, False, lin)
            ynew[1:n + 1] = y[1:n + 1] + a21 * ak1[1:n + 1]
            fcn(x + c2 * h, ynewview, dy[1:])
            ynew[1:n + 1] = hc21 * ak1[1:n + 1]
            _slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, m2, nm1,
                    fac, e, ip, dy, ak2, fx, ynew, hd2, ijob, True, lin)

            ynew[1:n + 1] = y[1:n + 1] + a31 * ak1[1:n + 1] + a32 * ak2[1:n + 1]
            fcn(x + c3 * h, ynewview, dy[1:])
            ynew[1:n + 1] = hc31 * ak1[1:n + 1] + hc32 * ak2[1:n + 1]
            _slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, m2, nm1,
                    fac, e, ip, dy, ak3, fx, ynew, hd3, ijob, True, lin)

            ynew[1:n + 1] = (y[1:n + 1] + a41 * ak1[1:n + 1]
                             + a42 * ak2[1:n + 1] + a43 * ak3[1:n + 1])
            fcn(x + c4 * h, ynewview, dy[1:])
            ynew[1:n + 1] = (hc41 * ak1[1:n + 1] + hc42 * ak2[1:n + 1]
                             + hc43 * ak3[1:n + 1])
            _slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, m2, nm1,
                    fac, e, ip, dy, ak4, fx, ynew, hd4, ijob, True, lin)

            ynew[1:n + 1] = (y[1:n + 1] + a51 * ak1[1:n + 1]
                             + a52 * ak2[1:n + 1] + a53 * ak3[1:n + 1]
                             + a54 * ak4[1:n + 1])
            fcn(x + h, ynewview, dy[1:])
            ak6[1:n + 1] = (hc52 * ak2[1:n + 1] + hc54 * ak4[1:n + 1]
                            + hc51 * ak1[1:n + 1] + hc53 * ak3[1:n + 1])
            _slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, m2, nm1,
                    fac, e, ip, dy, ak5, fx, ak6, 0.0, ijob, True, lin)

            # ------------ EMBEDDED SOLUTION ---------------
            ynew[1:n + 1] += ak5[1:n + 1]
            fcn(x + h, ynewview, dy[1:])
            cont[1:n + 1] = (hc61 * ak1[1:n + 1] + hc62 * ak2[1:n + 1]
                             + hc65 * ak5[1:n + 1] + hc64 * ak4[1:n + 1]
                             + hc63 * ak3[1:n + 1])
            _slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, m2, nm1,
                    fac, e, ip, dy, ak6, fx, cont, 0.0, ijob, True, lin)

            # ------------ NEW SOLUTION ---------------
            ynew[1:n + 1] += ak6[1:n + 1]
            nsol += 6
            nfcn += 5

            # ------------ DENSE OUTPUT ----------
            if iout != 0:
                cont[1:n + 1] = y[1:n + 1]
                cont[1 + nn2:n + 1 + nn2] = (
                    d21 * ak1[1:n + 1] + d22 * ak2[1:n + 1]
                    + d23 * ak3[1:n + 1] + d24 * ak4[1:n + 1]
                    + d25 * ak5[1:n + 1])
                cont[1 + nn3:n + 1 + nn3] = (
                    d31 * ak1[1:n + 1] + d32 * ak2[1:n + 1]
                    + d33 * ak3[1:n + 1] + d34 * ak4[1:n + 1]
                    + d35 * ak5[1:n + 1])

            # --- ERROR ESTIMATION ---
            nstep += 1
            # A SEQUENTIAL SUM, as the Fortran accumulates it. np.sum would
            # pairwise-sum and differ in the last bits, which is enough to move
            # a step-size decision and desynchronise the whole trajectory.
            err = 0.0
            for i in range(1, n + 1):
                if itol == 0:
                    sk = atol[1] + rtol[1] * max(abs(y[i]), abs(ynew[i]))
                else:
                    sk = atol[i] + rtol[i] * max(abs(y[i]), abs(ynew[i]))
                q = ak6[i] / sk
                err += q * q
            err = math.sqrt(err / n)

            # --- COMPUTATION OF HNEW, .2 <= HNEW/H <= 6 ---
            fac = max(fac2, min(fac1, err ** 0.25 / safe))
            hnew = h / fac

            if err <= 1.0:
                # --- STEP IS ACCEPTED ---
                naccpt += 1
                if pred:
                    if naccpt > 1:
                        facgus = (hacc / h) * (err * err / erracc) ** 0.25 / safe
                        facgus = max(fac2, min(fac1, facgus))
                        fac = max(fac, facgus)
                        hnew = h / fac
                    hacc = h
                    erracc = max(.01, err)
                y[1:n + 1] = ynew[1:n + 1]
                dense.xold = x
                x += h
                if iout != 0:
                    cont[n + 1:2 * n + 1] = y[1:n + 1]
                    dense.h = h
                    irtrn = solout(naccpt + 1, dense.xold, x, yview, dense)
                    if irtrn is not None and irtrn < 0:
                        return _finish(IDID_SOLOUT_STOP)
                if abs(hnew) > hmaxn:
                    hnew = posneg * hmaxn
                if reject:
                    hnew = posneg * min(abs(hnew), abs(h))
                reject = False
                h = hnew
                break  # back to label 1
            else:
                # --- STEP IS REJECTED ---
                reject = True
                last = False
                h = hnew
                if naccpt >= 1:
                    nrejct += 1
                continue  # back to label 2
