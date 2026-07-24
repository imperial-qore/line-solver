"""
Laplace-domain transient analysis of single-class open queues via the transient
QBD method plus numerical inverse Laplace transform.

Supports single-server MAP/MAP/1 (infinite buffer) and MAP/MAP/1/N (finite
buffer), reading arrival and service uniformly as (D0, D1) MAPs from sn.proc;
this subsumes M/M/1, M/PH/1 and correlated-MAP arrival/service that the
libQBD/expm fast path (ldqbd_transient) cannot represent exactly.

Python-native port of the MATLAB transient-QBD solver (no JVM). The transient
transforms mam_transient2_open / mam_transient2 are transcribed from the MATLAB
reference using 1-based-padded lists (index 0 unused) to mirror the block
indexing exactly.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np

from line_solver.lib.thirdparty.iltcme import matlab_ilt, matlab_ilt_matrix
from line_solver.api.mam.map_analysis import map_lambda, map_pie, map_prob


# ---------------------------------------------------------------------------
# Small linear-algebra helpers
# ---------------------------------------------------------------------------

def _mrdivide(X, Y):
    """MATLAB X / Y = X * inv(Y), solved as (Y' \\ X')'."""
    return np.linalg.solve(Y.conj().T, X.conj().T).conj().T if np.iscomplexobj(Y) or np.iscomplexobj(X) \
        else np.linalg.solve(Y.T, X.T).T


def _mldivide(A, B):
    """MATLAB A \\ B = inv(A) * B."""
    return np.linalg.solve(A, B)


def _mxpow(A, k):
    """Matrix power returning the identity when k == 0."""
    if k == 0:
        return np.eye(A.shape[0], dtype=A.dtype)
    return np.linalg.matrix_power(A, k)


def _krons(A, B):
    """Kronecker sum A (+) B = kron(A, I) + kron(I, B)."""
    na = A.shape[0]
    nb = B.shape[0]
    return np.kron(A, np.eye(nb)) + np.kron(np.eye(na), B)


def _eye_like(A):
    return np.eye(A.shape[0], dtype=complex)


# ---------------------------------------------------------------------------
# Fundamental matrices (cyclic reduction), complex-capable
# ---------------------------------------------------------------------------

def qbd_fundmat(B, L, F, matrices='G', precision=1e-14, max_num_it=50):
    """Cyclic-reduction fundamental matrices G (and optionally R) of a QBD.

    Operates on the raw level blocks (B, L, F) and supports complex L (as in a
    Laplace-domain analysis where L is shifted by -s*I). Returns a tuple in the
    order of the character codes in *matrices* (subset of 'G', 'R').
    """
    m = L.shape[0]
    II = np.eye(m, dtype=complex)
    lamb = np.max(-np.real(np.diag(L)))
    Bm = B / lamb
    Lm = L / lamb + II
    Fm = F / lamb

    BF = _mldivide(II - Lm, II)
    BB = BF @ Fm
    BF = BF @ Bm
    G = BF.copy()
    PI = BB.copy()
    check = 1.0
    numit = 0
    while check > precision and numit < max_num_it:
        Lstar = BF @ BB + BB @ BF
        Bstar = BB @ BB
        Fstar = BF @ BF
        BB = _mldivide(II - Lstar, II)
        BF = BB @ Fstar
        BB = BB @ Bstar
        G = G + PI @ BF
        PI = PI @ BB
        check = min(np.linalg.norm(BB, np.inf), np.linalg.norm(BF, np.inf))
        numit += 1

    R = None
    outs = []
    for c in matrices:
        if c == 'G':
            outs.append(G)
        elif c == 'R':
            if R is None:
                R = Fm @ _mldivide(II - (Lm + Fm @ G), II)
            outs.append(R)
        else:
            raise ValueError("unknown matrix code '%s' in '%s'" % (c, matrices))
    return tuple(outs)


# ---------------------------------------------------------------------------
# Transient transform: open (infinite) QBD
# ---------------------------------------------------------------------------

def mam_transient2_open(B, L, F, Lv, T, n, m, s):
    """Laplace-domain transient V(s,n,m) for an open (infinite) QBD.

    Block lists B, L, F, Lv and threshold list T are 1-based (index 0 unused),
    with K = len(T) - 1 regimes; the last regime repeats to infinity.
    """
    K = len(T) - 1  # number of regimes (1-based lists padded at index 0)

    Gs = [None] * (K + 1)
    Rs = [None] * (K + 1)
    Ghs = [None] * (K + 1)
    Rhs = [None] * (K + 1)
    for k in range(1, K + 1):
        if k < K and T[k + 1] - T[k] == 1:
            continue
        Ik = _eye_like(L[k])
        Gs[k], Rs[k] = qbd_fundmat(B[k], L[k] - s * Ik, F[k], 'GR')
        Ghs[k], Rhs[k] = qbd_fundmat(F[k], L[k] - s * Ik, B[k], 'GR')

    SvHn = [None] * (K + 1)
    SvH0 = [None] * (K + 1)
    SvHhn = [None] * (K + 1)
    SvHh0 = [None] * (K + 1)
    for k in range(1, K):
        NN = Lv[k].shape[0]
        if T[k + 1] - T[k] > 1:
            d = T[k + 1] - T[k]
            Ik = np.eye(NN, dtype=complex)
            num = np.block([[_mxpow(Ghs[k], d - 1), Gs[k]], [Ghs[k], _mxpow(Gs[k], d - 1)]])
            den = np.block([[Ik, _mxpow(Gs[k], d)], [_mxpow(Ghs[k], d), Ik]])
            SH = _mrdivide(num, den)
            SvHn[k] = SH[0:NN, 0:NN]
            SvH0[k] = SH[0:NN, NN:2 * NN]
            SvHhn[k] = SH[NN:2 * NN, 0:NN]
            SvHh0[k] = SH[NN:2 * NN, NN:2 * NN]
        else:
            NN1 = Lv[k + 1].shape[0]
            SvH0[k] = np.zeros((NN1, NN), dtype=complex)
            SvHh0[k] = np.eye(NN, dtype=complex)
            SvHhn[k] = np.zeros((NN, NN1), dtype=complex)
            SvHn[k] = np.eye(NN1, dtype=complex)

    SY = [None] * (K + 1)
    SY[K] = Gs[K]
    for k in range(K - 1, 0, -1):
        NNk1 = Lv[k + 1].shape[0]
        SY[k] = SvH0[k] + SvHn[k] @ _mldivide(
            s * np.eye(NNk1) - Lv[k + 1] - F[k + 1] @ SY[k + 1] - B[k] @ SvHhn[k],
            B[k] @ SvHh0[k])

    NN1 = Lv[1].shape[0]
    SYh = [None] * (K + 1)
    if K >= 1:
        SYh[1] = SvHhn[1] + SvHh0[1] @ _mldivide(
            s * np.eye(NN1) - Lv[1] - F[1] @ SvH0[1], F[1] @ SvHn[1])
    for k in range(2, K):
        NNk = Lv[k].shape[0]
        SYh[k] = SvHhn[k] + SvHh0[k] @ _mldivide(
            s * np.eye(NNk) - Lv[k] - B[k - 1] @ SYh[k - 1] - F[k] @ SvH0[k], F[k] @ SvHn[k])

    SV = [[None] * (K + 1) for _ in range(K + 1)]
    for l in range(0, K):
        if l == 0:
            SV[1][1] = _mldivide(s * np.eye(NN1) - Lv[1] - F[1] @ SY[1], np.eye(NN1))
        else:
            NNl1 = Lv[l + 1].shape[0]
            SV[l + 1][l + 1] = _mldivide(
                s * np.eye(NNl1) - Lv[l + 1] - F[l + 1] @ SY[l + 1] - B[l] @ SYh[l], np.eye(NNl1))
        for k in range(l + 1, K):
            NNk1 = Lv[k + 1].shape[0]
            SV[k + 1][l + 1] = _mldivide(
                s * np.eye(NNk1) - Lv[k + 1] - F[k + 1] @ SY[k + 1] - B[k] @ SvHhn[k],
                B[k] @ SvHh0[k] @ SV[k][l + 1])
        for k in range(l - 1, 0, -1):
            NNk1 = Lv[k + 1].shape[0]
            SV[k + 1][l + 1] = _mldivide(
                s * np.eye(NNk1) - Lv[k + 1] - F[k + 1] @ SvH0[k + 1] - B[k] @ SYh[k],
                F[k + 1] @ SvHn[k + 1] @ SV[k + 2][l + 1])
        if l > 0:
            SV[1][l + 1] = _mldivide(
                s * np.eye(NN1) - Lv[1] - F[1] @ SvH0[1], F[1] @ SvHn[1] @ SV[2][l + 1])

    kn = _regime_of(T, n, K)
    km = _regime_of(T, m, K)

    NN = Lv[kn].shape[0]
    II = np.eye(NN, dtype=complex)

    Vu = None
    Lu = None
    threshold_n = (T[kn] == n)
    if threshold_n:
        if T[km] == m:
            return SV[kn][km]
        Vl = SV[kn][km]
        Ll = T[km]
        if km < K:
            Vu = SV[kn][km + 1]
            Lu = T[km + 1]
    else:
        if kn < K:
            d1 = T[kn + 1] - n
            num = np.block([[_mxpow(Ghs[kn], d1 - 1), Gs[kn]], [Ghs[kn], _mxpow(Gs[kn], d1 - 1)]])
            den = np.block([[II, _mxpow(Gs[kn], d1)], [_mxpow(Ghs[kn], d1), II]])
            Tmp = _mrdivide(num, den)
            HTnn = Tmp[0:NN, 0:NN]
            HTn0 = Tmp[0:NN, NN:2 * NN]
            HhTnn = Tmp[NN:2 * NN, 0:NN]
            HhTn0 = Tmp[NN:2 * NN, NN:2 * NN]

            d2 = n - T[kn]
            num = np.block([[_mxpow(Ghs[kn], d2 - 1), Gs[kn]], [Ghs[kn], _mxpow(Gs[kn], d2 - 1)]])
            den = np.block([[II, _mxpow(Gs[kn], d2)], [_mxpow(Ghs[kn], d2), II]])
            Tmp = _mrdivide(num, den)
            HnTn = Tmp[0:NN, 0:NN]
            HnT0 = Tmp[0:NN, NN:2 * NN]
            HhnTn = Tmp[NN:2 * NN, 0:NN]
            HhnT0 = Tmp[NN:2 * NN, NN:2 * NN]

            NNkn1 = Lv[kn + 1].shape[0]
            Yn = HTn0 + HTnn @ _mldivide(
                s * np.eye(NNkn1) - Lv[kn + 1] - F[kn + 1] @ SY[kn + 1] - B[kn] @ HhTnn,
                B[kn] @ HhTn0)
        else:
            d2 = n - T[kn]
            num = np.block([[_mxpow(Ghs[kn], d2 - 1), Gs[kn]], [Ghs[kn], _mxpow(Gs[kn], d2 - 1)]])
            den = np.block([[II, _mxpow(Gs[kn], d2)], [_mxpow(Ghs[kn], d2), II]])
            Tmp = _mrdivide(num, den)
            HnTn = Tmp[0:NN, 0:NN]
            HnT0 = Tmp[0:NN, NN:2 * NN]
            HhnTn = Tmp[NN:2 * NN, 0:NN]
            HhnT0 = Tmp[NN:2 * NN, NN:2 * NN]
            HTnn = HTn0 = HhTnn = HhTn0 = None
            Yn = Gs[kn]

        if kn == 1:
            Yhn = HhnTn + HhnT0 @ _mldivide(
                s * np.eye(NN1) - Lv[1] - F[1] @ HnT0, F[1] @ HnTn)
        else:
            NNkn = Lv[kn].shape[0]
            Yhn = HhnTn + HhnT0 @ _mldivide(
                s * np.eye(NNkn) - Lv[kn] - B[kn - 1] @ SYh[kn - 1] - F[kn] @ HnT0, F[kn] @ HnTn)

        Mkn = s * np.eye(L[kn].shape[0]) - L[kn]
        if T[km] < n:
            Vnl = _mldivide(Mkn - B[kn] @ HhnTn - F[kn] @ Yn, B[kn] @ HhnT0 @ SV[kn][km])
        else:
            Vnl = _mldivide(Mkn - F[kn] @ HTn0 - B[kn] @ Yhn, F[kn] @ HTnn @ SV[kn + 1][km])
        if m == T[km]:
            return Vnl

        Vnu = None
        if km < K:
            if T[km + 1] < n:
                Vnu = _mldivide(Mkn - B[kn] @ HhnTn - F[kn] @ Yn, B[kn] @ HhnT0 @ SV[kn][km + 1])
            else:
                Vnu = _mldivide(Mkn - F[kn] @ HTn0 - B[kn] @ Yhn, F[kn] @ HTnn @ SV[kn + 1][km + 1])
        Vnn = _mldivide(s * II - L[kn] - B[kn] @ Yhn - F[kn] @ Yn, II)
        if n == m:
            return Vnn
        if km == K and n < m:
            if kn < K:
                return Vnl @ _mxpow(Rs[km], m - T[km])
            return Vnn @ _mxpow(Rs[km], m - n)
        if kn != km:
            Vu, Vl, Lu, Ll = Vnu, Vnl, T[km + 1], T[km]
        elif n <= m:
            Vu, Vl, Lu, Ll = Vnu, Vnn, T[km + 1], n
        else:
            Vu, Vl, Lu, Ll = Vnn, Vnl, n, T[km]

    if km == K and n < m:
        if kn < K:
            return Vl @ _mxpow(Rs[km], m - T[km])
        return SV[kn][kn] @ _mxpow(Rs[km], m - n)

    NN = Rs[km].shape[0]
    II = np.eye(NN, dtype=complex)
    Zden = np.block([[II, _mxpow(Rs[km], Lu - Ll)], [_mxpow(Rhs[km], Lu - Ll), II]])
    Znum = np.block([[_mxpow(Rs[km], m - Ll)], [_mxpow(Rhs[km], Lu - m)]])
    Z = _mldivide(Zden, Znum)
    return Vl @ Z[0:NN, :] + Vu @ Z[NN:2 * NN, :]


# ---------------------------------------------------------------------------
# Transient transform: closed (finite) QBD
# ---------------------------------------------------------------------------

def mam_transient2(B, L, F, Lv, T, n, m, s):
    """Laplace-domain transient V(s,n,m) for a finite QBD.

    Block lists are 1-based (index 0 unused); K = len(T) - 2 regimes with
    Lv[K+1] the top boundary level (T is 1-based of length K+2).
    """
    K = len(T) - 2  # T is 1-based with entries T[1..K+1]

    Gs = [None] * (K + 1)
    Rs = [None] * (K + 1)
    Ghs = [None] * (K + 1)
    Rhs = [None] * (K + 1)
    for k in range(1, K + 1):
        if T[k + 1] - T[k] == 1:
            continue
        Ik = _eye_like(L[k])
        Gs[k], Rs[k] = qbd_fundmat(B[k], L[k] - s * Ik, F[k], 'GR')
        Ghs[k], Rhs[k] = qbd_fundmat(F[k], L[k] - s * Ik, B[k], 'GR')

    SvHn = [None] * (K + 2)
    SvH0 = [None] * (K + 2)
    SvHhn = [None] * (K + 2)
    SvHh0 = [None] * (K + 2)
    for k in range(1, K + 1):
        NN = Lv[k].shape[0]
        if T[k + 1] - T[k] > 1:
            d = T[k + 1] - T[k]
            Ik = np.eye(NN, dtype=complex)
            num = np.block([[_mxpow(Ghs[k], d - 1), Gs[k]], [Ghs[k], _mxpow(Gs[k], d - 1)]])
            den = np.block([[Ik, _mxpow(Gs[k], d)], [_mxpow(Ghs[k], d), Ik]])
            SH = _mrdivide(num, den)
            SvHn[k] = SH[0:NN, 0:NN]
            SvH0[k] = SH[0:NN, NN:2 * NN]
            SvHhn[k] = SH[NN:2 * NN, 0:NN]
            SvHh0[k] = SH[NN:2 * NN, NN:2 * NN]
        else:
            NN1 = Lv[k + 1].shape[0]
            SvH0[k] = np.zeros((NN1, NN), dtype=complex)
            SvHh0[k] = np.eye(NN, dtype=complex)
            SvHhn[k] = np.zeros((NN, NN1), dtype=complex)
            SvHn[k] = np.eye(NN1, dtype=complex)

    NNK = Lv[K + 1].shape[0]
    SY = [None] * (K + 1)
    SY[K] = SvH0[K] + SvHn[K] @ _mldivide(
        s * np.eye(NNK) - Lv[K + 1] - B[K] @ SvHhn[K], B[K] @ SvHh0[K])
    for k in range(K - 1, 0, -1):
        NNk1 = Lv[k + 1].shape[0]
        SY[k] = SvH0[k] + SvHn[k] @ _mldivide(
            s * np.eye(NNk1) - Lv[k + 1] - F[k + 1] @ SY[k + 1] - B[k] @ SvHhn[k],
            B[k] @ SvHh0[k])

    NN1 = Lv[1].shape[0]
    SYh = [None] * (K + 1)
    SYh[1] = SvHhn[1] + SvHh0[1] @ _mldivide(
        s * np.eye(NN1) - Lv[1] - F[1] @ SvH0[1], F[1] @ SvHn[1])
    for k in range(2, K + 1):
        NNk = Lv[k].shape[0]
        SYh[k] = SvHhn[k] + SvHh0[k] @ _mldivide(
            s * np.eye(NNk) - Lv[k] - B[k - 1] @ SYh[k - 1] - F[k] @ SvH0[k], F[k] @ SvHn[k])

    SV = [[None] * (K + 2) for _ in range(K + 2)]
    for l in range(0, K + 1):
        if l == 0:
            SV[1][1] = _mldivide(s * np.eye(NN1) - Lv[1] - F[1] @ SY[1], np.eye(NN1))
        elif l == K:
            SV[K + 1][K + 1] = _mldivide(s * np.eye(NNK) - Lv[K + 1] - B[K] @ SYh[K], np.eye(NNK))
        else:
            NNl1 = Lv[l + 1].shape[0]
            SV[l + 1][l + 1] = _mldivide(
                s * np.eye(NNl1) - Lv[l + 1] - F[l + 1] @ SY[l + 1] - B[l] @ SYh[l], np.eye(NNl1))
        for k in range(l + 1, K):
            NNk1 = Lv[k + 1].shape[0]
            SV[k + 1][l + 1] = _mldivide(
                s * np.eye(NNk1) - Lv[k + 1] - F[k + 1] @ SY[k + 1] - B[k] @ SvHhn[k],
                B[k] @ SvHh0[k] @ SV[k][l + 1])
        if l < K:
            SV[K + 1][l + 1] = _mldivide(
                s * np.eye(NNK) - Lv[K + 1] - B[K] @ SvHhn[K],
                B[K] @ SvHh0[K] @ SV[K][l + 1])
        for k in range(l - 1, 0, -1):
            NNk1 = Lv[k + 1].shape[0]
            SV[k + 1][l + 1] = _mldivide(
                s * np.eye(NNk1) - Lv[k + 1] - F[k + 1] @ SvH0[k + 1] - B[k] @ SYh[k],
                F[k + 1] @ SvHn[k + 1] @ SV[k + 2][l + 1])
        if l > 0:
            SV[1][l + 1] = _mldivide(
                s * np.eye(NN1) - Lv[1] - F[1] @ SvH0[1], F[1] @ SvHn[1] @ SV[2][l + 1])

    kn = _regime_of(T, n, K + 1)
    km = _regime_of(T, m, K + 1)

    NN = Lv[kn].shape[0]
    II = np.eye(NN, dtype=complex)

    if T[kn] == n:
        if T[km] == m:
            return SV[kn][km]
        Vu = SV[kn][km + 1]
        Vl = SV[kn][km]
        Lu = T[km + 1]
        Ll = T[km]
    else:
        d1 = T[kn + 1] - n
        num = np.block([[_mxpow(Ghs[kn], d1 - 1), Gs[kn]], [Ghs[kn], _mxpow(Gs[kn], d1 - 1)]])
        den = np.block([[II, _mxpow(Gs[kn], d1)], [_mxpow(Ghs[kn], d1), II]])
        Tmp = _mrdivide(num, den)
        HTnn = Tmp[0:NN, 0:NN]
        HTn0 = Tmp[0:NN, NN:2 * NN]
        HhTnn = Tmp[NN:2 * NN, 0:NN]
        HhTn0 = Tmp[NN:2 * NN, NN:2 * NN]

        d2 = n - T[kn]
        num = np.block([[_mxpow(Ghs[kn], d2 - 1), Gs[kn]], [Ghs[kn], _mxpow(Gs[kn], d2 - 1)]])
        den = np.block([[II, _mxpow(Gs[kn], d2)], [_mxpow(Ghs[kn], d2), II]])
        Tmp = _mrdivide(num, den)
        HnTn = Tmp[0:NN, 0:NN]
        HnT0 = Tmp[0:NN, NN:2 * NN]
        HhnTn = Tmp[NN:2 * NN, 0:NN]
        HhnT0 = Tmp[NN:2 * NN, NN:2 * NN]

        if kn == K:
            Yn = HTn0 + HTnn @ _mldivide(
                s * np.eye(NNK) - Lv[K + 1] - B[K] @ HhTnn, B[K] @ HhTn0)
        else:
            NNkn1 = Lv[kn + 1].shape[0]
            Yn = HTn0 + HTnn @ _mldivide(
                s * np.eye(NNkn1) - Lv[kn + 1] - F[kn + 1] @ SY[kn + 1] - B[kn] @ HhTnn,
                B[kn] @ HhTn0)
        if kn == 1:
            Yhn = HhnTn + HhnT0 @ _mldivide(
                s * np.eye(NN1) - Lv[1] - F[1] @ HnT0, F[1] @ HnTn)
        else:
            NNkn = Lv[kn].shape[0]
            Yhn = HhnTn + HhnT0 @ _mldivide(
                s * np.eye(NNkn) - Lv[kn] - B[kn - 1] @ SYh[kn - 1] - F[kn] @ HnT0, F[kn] @ HnTn)

        Mkn = s * np.eye(L[kn].shape[0]) - L[kn]
        if T[km] < n:
            Vnl = _mldivide(Mkn - B[kn] @ HhnTn - F[kn] @ Yn, B[kn] @ HhnT0 @ SV[kn][km])
        else:
            Vnl = _mldivide(Mkn - F[kn] @ HTn0 - B[kn] @ Yhn, F[kn] @ HTnn @ SV[kn + 1][km])
        if m == T[km]:
            return Vnl
        if T[km + 1] < n:
            Vnu = _mldivide(Mkn - B[kn] @ HhnTn - F[kn] @ Yn, B[kn] @ HhnT0 @ SV[kn][km + 1])
        else:
            Vnu = _mldivide(Mkn - F[kn] @ HTn0 - B[kn] @ Yhn, F[kn] @ HTnn @ SV[kn + 1][km + 1])
        Vnn = _mldivide(s * II - L[kn] - B[kn] @ Yhn - F[kn] @ Yn, II)
        if n == m:
            return Vnn
        if kn != km:
            Vu, Vl, Lu, Ll = Vnu, Vnl, T[km + 1], T[km]
        elif n <= m:
            Vu, Vl, Lu, Ll = Vnu, Vnn, T[km + 1], n
        else:
            Vu, Vl, Lu, Ll = Vnn, Vnl, n, T[km]

    NN = Rs[km].shape[0]
    II = np.eye(NN, dtype=complex)
    Zden = np.block([[II, _mxpow(Rs[km], Lu - Ll)], [_mxpow(Rhs[km], Lu - Ll), II]])
    Znum = np.block([[_mxpow(Rs[km], m - Ll)], [_mxpow(Rhs[km], Lu - m)]])
    Z = _mldivide(Zden, Znum)
    return Vl @ Z[0:NN, :] + Vu @ Z[NN:2 * NN, :]


def _regime_of(T, level, kmax):
    """Return the 1-based regime index k such that T[k] <= level < T[k+1],
    clamped to kmax when level is at or beyond the last threshold."""
    for k in range(1, len(T)):
        if T[k] > level:
            return k - 1
    return kmax


# ---------------------------------------------------------------------------
# sn -> blocks bridge + metric reduction
# ---------------------------------------------------------------------------

def _proc_ph_entry(sn, idx):
    """Extract the per-station process entry from the native sn.proc, which may
    be a list, a dict keyed by idx, or a dict keyed by (idx, class)."""
    proc = sn.proc
    if isinstance(proc, list):
        return proc[idx][0] if isinstance(proc[idx], list) else proc[idx]
    if isinstance(proc, dict):
        return proc.get((idx, 0), proc.get(idx, None))
    return proc[idx]


def _as_map_block(x):
    """Coerce a proc element to a 2D float matrix, or None when it is not a
    matrix (scalar distribution parameters, parameter dicts, empty)."""
    if x is None or isinstance(x, dict) or np.isscalar(x):
        return None
    arr = np.asarray(x)
    if arr.dtype == object or arr.ndim == 0 or arr.size == 0:
        return None
    return np.atleast_2d(arr.astype(float))


def _ph_pair_to_map(a0, M0):
    """Return (D0, D1) for a proc entry stored in the native [alpha, T] phase-type
    form, or None when the pair is not that form.

    Unlike MATLAB, which always stores sn.proc{i,r} as the MAP pair {D0, D1},
    the native Python struct records PH/APH/Coxian/Cox2 as [alpha, T]: the
    initial phase probability row and the subgenerator. The two forms are
    distinguished by alpha being a nonnegative row of length T.shape[0] summing
    to one. The equivalent renewal MAP is D0 = T, D1 = t0 * alpha with the exit
    rates t0 = -T * 1 re-initialized by alpha."""
    a0 = np.asarray(a0, dtype=float).ravel()
    M0 = np.atleast_2d(np.asarray(M0, dtype=float))
    if M0.shape[0] != M0.shape[1] or a0.size != M0.shape[0]:
        return None
    if a0.min() < -1e-12 or abs(a0.sum() - 1.0) > 1e-6:
        return None
    t0 = -M0 @ np.ones(M0.shape[0])
    return M0, np.atleast_2d(np.outer(t0, a0))


def _read_map(entry):
    """Return (D0, D1) from a native proc entry in any supported format:
    {'rate': mu}, [D0, D1], [alpha, T], {0: D0, 1: D1}, or a bare D0. Entries
    that are distribution-parameter forms rather than MAP blocks (e.g. Erlang
    [{'k','mu'}], Gamma [shape, scale]) yield (None, None) so callers decline
    instead of crashing."""
    if isinstance(entry, dict) and 'rate' in entry:
        r = float(entry['rate'])
        return np.array([[-r]]), np.array([[r]])
    if isinstance(entry, (list, tuple)):
        D0 = _as_map_block(entry[0]) if len(entry) > 0 else None
        if D0 is None:
            return None, None
        D1 = _as_map_block(entry[1]) if len(entry) > 1 else None
        if D1 is not None:
            ph = _ph_pair_to_map(entry[0], entry[1])
            if ph is not None:
                return ph
        return D0, D1
    if isinstance(entry, dict) and 0 in entry:
        D0 = _as_map_block(entry[0])
        if D0 is None:
            return None, None
        D1 = _as_map_block(entry.get(1)) if 1 in entry else None
        return D0, D1
    D0 = _as_map_block(entry)
    return (D0, None) if D0 is not None else (None, None)


def transient_qbd_applicable(sn):
    """True when transient analysis should use the Laplace transient QBD solver
    (single-server open queue with non-Poisson arrival or correlated-MAP
    service) rather than the libQBD/expm fast path."""
    if sn.nclasses != 1:
        return False
    njobs = np.asarray(sn.njobs).flatten()
    if not np.isinf(njobs[0]):
        return False

    source_idx, queue_idx = _find_source_queue(sn)
    if source_idx is None or queue_idx is None:
        return False
    nservers = np.asarray(sn.nservers).flatten()
    if int(nservers[queue_idx]) != 1:
        return False

    Da0, Da1 = _read_map(_proc_ph_entry(sn, source_idx))
    Ds0, Ds1 = _read_map(_proc_ph_entry(sn, queue_idx))
    if Da1 is None or Ds1 is None:
        return False   # cannot form MAP blocks; leave to the fast path
    arrival_is_poisson = (Da0.shape[0] == 1)
    service_is_renewal = _is_renewal_map(Ds0, Ds1)
    return (not arrival_is_poisson) or (not service_is_renewal)


def _is_renewal_map(D0, D1):
    ns = D0.shape[0]
    if ns == 1:
        return True
    t_exit = D1 @ np.ones(ns)
    alpha = map_pie(D0, D1)
    return np.linalg.norm(D1 - np.outer(t_exit, alpha), 'fro') < 1e-9 * max(1.0, np.linalg.norm(D1, 'fro'))


def _find_source_queue(sn):
    sched = sn.sched if isinstance(sn.sched, dict) else {}
    source_idx = None
    queue_idx = None
    for i in range(sn.nstations):
        sval = sched.get(i, None) if isinstance(sched, dict) else None
        sname = getattr(sval, 'name', None) if sval is not None else None
        if sname is None:
            arr = np.asarray(sn.sched).flatten() if not isinstance(sn.sched, dict) else None
            if arr is not None:
                sname = None
                sval = arr[i]
        if sname == 'EXT' or sval == 16:
            source_idx = i
        elif sname == 'FCFS' or sval == 0:
            queue_idx = i
    return source_idx, queue_idx


def solver_mam_transient_qbd(sn, options):
    """Transient analysis of a single-class open queue via Laplace transient QBD.

    Returns (Qt, Ut, Tt): each an M x K list-of-lists with the queue station's
    entry a (nTimePoints, 2) array of [metric, time]."""
    M = sn.nstations
    K = sn.nclasses
    if K != 1:
        raise ValueError("Transient QBD method requires a single-class model.")
    njobs = np.asarray(sn.njobs).flatten()
    if not np.isinf(njobs[0]):
        raise ValueError("Transient QBD method requires an open model.")

    source_idx, queue_idx = _find_source_queue(sn)
    if source_idx is None or queue_idx is None:
        raise ValueError("Transient QBD method requires exactly one Source and one FCFS Queue.")
    nservers = np.asarray(sn.nservers).flatten()
    if int(nservers[queue_idx]) != 1:
        raise ValueError("Transient QBD (Laplace) method supports single-server queues only.")

    Da0, Da1 = _read_map(_proc_ph_entry(sn, source_idx))
    Ds0, Ds1 = _read_map(_proc_ph_entry(sn, queue_idx))
    na = Da0.shape[0]
    ns = Ds0.shape[0]
    Ina = np.eye(na)
    Ins = np.eye(ns)

    Lrep = _krons(Da0, Ds0)
    Frep = np.kron(Da1, Ins)
    Brep = np.kron(Ina, Ds1)
    Lv0 = np.kron(Da0, Ins)
    B0 = np.kron(Ina, Ds1)
    F0 = np.kron(Da1, Ins)

    pi_arr = map_prob(Da0, Da1)
    pi_svc = map_prob(Ds0, Ds1)
    pi0 = np.kron(pi_arr, pi_svc).reshape(1, -1)

    s_exit = Ds1 @ np.ones(ns)
    w_dep = np.kron(np.ones(na), s_exit).reshape(-1, 1)
    w_one = np.ones((na * ns, 1))

    cap = np.asarray(sn.cap).flatten()
    buf_cap = cap[queue_idx]
    is_finite = not np.isinf(buf_cap)

    t_start, t_end = options.timespan[0], options.timespan[1]
    dur = t_end - t_start
    n_time = int(min(101, max(11, round(dur * 10))))
    times = np.linspace(t_start, t_end, n_time)
    # see _kb/06-solver-catalog.md (MAM: "MAP/MAP/1 exact fast-path, LDQBD,
    # transient QBD, traffic superposition") for the t=0 ILT singularity
    pos_mask = times > 0
    tpos = times[pos_mask]
    max_fn_evals = 100
    if hasattr(options, 'iter_max') and options.iter_max and options.iter_max > 1:
        max_fn_evals = int(min(1000, max(11, round(options.iter_max))))

    if is_finite:
        Ncap = int(buf_cap)
        Lv_top = np.kron(Da0 + Da1, Ins) + np.kron(Ina, Ds0)
        # 1-based padded lists; K=1 regime, T=[_,0,Ncap], Lv[_,Lv0,Lv_top]
        Bc = [None, Brep]
        Lc = [None, Lrep]
        Fc = [None, Frep]
        Lvc = [None, Lv0, Lv_top]
        Tc = [None, 0, Ncap]

        def en_fun(s):
            acc = 0.0
            for mm in range(1, Ncap + 1):
                V = mam_transient2(Bc, Lc, Fc, Lvc, Tc, 0, mm, s)
                acc = acc + mm * (pi0 @ V @ np.ones((V.shape[1], 1)))
            return acc

        def dep_fun(s):
            acc = 0.0
            for mm in range(1, Ncap + 1):
                V = mam_transient2(Bc, Lc, Fc, Lvc, Tc, 0, mm, s)
                acc = acc + pi0 @ V @ w_dep
            return acc

        def p0_fun(s):
            V = mam_transient2(Bc, Lc, Fc, Lvc, Tc, 0, 0, s)
            return pi0 @ V @ w_one

        EN = np.zeros(n_time); DEP = np.zeros(n_time); P0 = np.ones(n_time)
        EN[pos_mask] = matlab_ilt(lambda s: complex(np.asarray(en_fun(s)).ravel()[0]), tpos, max_fn_evals)
        DEP[pos_mask] = matlab_ilt(lambda s: complex(np.asarray(dep_fun(s)).ravel()[0]), tpos, max_fn_evals)
        P0[pos_mask] = matlab_ilt(lambda s: complex(np.asarray(p0_fun(s)).ravel()[0]), tpos, max_fn_evals)
    else:
        # Open: K=2 regimes, T=[_,0,1]; regime 1 = level 0, regime 2 repeats.
        Bo = [None, B0, Brep]
        Lo = [None, None, Lrep]
        Fo = [None, F0, Frep]
        Lvo = [None, Lv0, Lrep]
        To = [None, 0, 1]
        Krg = 2

        def _R_of(s):
            IK = np.eye(Lo[Krg].shape[0], dtype=complex)
            _, R = qbd_fundmat(Bo[Krg], Lo[Krg] - s * IK, Fo[Krg], 'GR')
            return R, IK

        # Homogeneous open case: V(s,0,m) = V(s,0,1) * R^(m-1), so the level sums
        # have exact closed forms (no truncation):
        #   E[N](s) = pi0 * V(s,0,1) * (I-R)^-2 * 1
        #   Tput(s) = pi0 * V(s,0,1) * (I-R)^-1 * wDep
        def en_fun(s):
            R, IK = _R_of(s)
            V1 = mam_transient2_open(Bo, Lo, Fo, Lvo, To, 0, 1, s)
            ImR = IK - R
            ones_col = np.ones((IK.shape[0], 1))
            return pi0 @ V1 @ _mldivide(ImR, _mldivide(ImR, ones_col))

        def dep_fun(s):
            R, IK = _R_of(s)
            V1 = mam_transient2_open(Bo, Lo, Fo, Lvo, To, 0, 1, s)
            return pi0 @ V1 @ _mldivide(IK - R, w_dep)

        def p0_fun(s):
            V = mam_transient2_open(Bo, Lo, Fo, Lvo, To, 0, 0, s)
            return pi0 @ V @ w_one

        EN = np.zeros(n_time); DEP = np.zeros(n_time); P0 = np.ones(n_time)
        EN[pos_mask] = matlab_ilt(lambda s: complex(np.asarray(en_fun(s)).ravel()[0]), tpos, max_fn_evals)
        DEP[pos_mask] = matlab_ilt(lambda s: complex(np.asarray(dep_fun(s)).ravel()[0]), tpos, max_fn_evals)
        P0[pos_mask] = matlab_ilt(lambda s: complex(np.asarray(p0_fun(s)).ravel()[0]), tpos, max_fn_evals)

    U = 1.0 - P0

    Qt = [[None] * K for _ in range(M)]
    Ut = [[None] * K for _ in range(M)]
    Tt = [[None] * K for _ in range(M)]
    Qt[queue_idx][0] = np.column_stack([EN, times])
    Ut[queue_idx][0] = np.column_stack([U, times])
    Tt[queue_idx][0] = np.column_stack([DEP, times])
    return Qt, Ut, Tt
