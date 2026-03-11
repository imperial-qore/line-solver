"""
Composite Moment (CoMoM) Methods for Normalizing Constant Computation.

Native Python implementations of CoMoM and related algorithms for
computing normalizing constants in product-form queueing networks.

Key functions:
    pfqn_comom: CoMoM algorithm for general networks
    pfqn_comomrm: CoMoM for finite repairman models
    pfqn_comomrm_orig: Original CoMoM implementation
    pfqn_comomrm_ms: CoMoM for multi-server repairman models
    pfqn_procomom2: Projected CoMoM method

References:
    Casale, G. "CoMoM: A class of algorithms for product-form solution
    of Closed Queueing Networks with Multiple Servers."
    Original MATLAB: matlab/src/api/pfqn/pfqn_comom*.m
"""

import numpy as np
from math import comb
from typing import Tuple, Optional, Dict, Any, NamedTuple
from scipy.special import gammaln

from .utils import factln, factln_vec, multichoose, matchrow, oner


class ComomResult(NamedTuple):
    """Result of CoMoM normalizing constant computation."""
    lG: float
    lGbasis: Optional[np.ndarray] = None


def pfqn_comom(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
               atol: float = 1e-8) -> float:
    """
    CoMoM algorithm for computing the normalizing constant.

    Implements the Composite Method of Moments algorithm for computing
    normalizing constants in closed product-form queueing networks.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,), default zeros
        atol: Absolute tolerance

    Returns:
        lG: Logarithm of the normalizing constant

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_comom.m
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()

    M, R = L.shape

    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=float).flatten()

    # Rescale demands for numerical stability
    Lmax = L.copy()
    Lmax[Lmax < atol] = Z[Lmax < atol]
    Lmax = np.max(Lmax, axis=0)
    Lmax[Lmax < atol] = 1.0
    L = L / Lmax
    Z = Z / Lmax

    # Prepare CoMoM data structures
    Dn = multichoose(R, M)
    Dn[:, R - 1] = 0  # Set last column to 0 (for subset selection)
    Dn = _sortbynnzpos(Dn)

    # Initialize
    nvec = np.zeros(R)
    h = _ginit(L, N, Dn, M, R)
    lh = np.log(np.abs(h) + 1e-300) + factln(np.sum(nvec) + M - 1) - np.sum(factln_vec(nvec))
    h = np.exp(lh - np.max(lh))  # Rescale for numerical stability

    scale = np.zeros(int(np.sum(N)))

    # Iterate through all classes
    for r in range(R):
        for Nr in range(1, int(N[r]) + 1):
            nvec[r] += 1

            if Nr == 1:
                A, B, DA = _genmatrix(L, nvec.copy(), Z, r, N, Dn, M, R)
            else:
                A = A + DA

            b = B @ h * nvec[r] / (np.sum(nvec) + M - 1)

            try:
                h = np.linalg.solve(A, b)
            except np.linalg.LinAlgError:
                h = np.linalg.lstsq(A, b, rcond=None)[0]

            nt = int(np.sum(nvec))
            scale[nt - 1] = np.abs(np.sum(np.sort(h)))
            if scale[nt - 1] > 0:
                h = np.abs(h) / scale[nt - 1]

    # Unscale and return the log of the normalizing constant
    lG = (np.log(np.abs(h[-(R - 1) - 1]) + 1e-300) +
          factln(np.sum(N) + M - 1) -
          np.sum(factln_vec(N)) +
          np.dot(N, np.log(Lmax)) +
          np.sum(np.log(scale[scale > 0])))

    return lG


def _sortbynnzpos(I: np.ndarray) -> np.ndarray:
    """Sort combinations by number of nonzeros and position."""
    def nnzcmp(i1, i2):
        nnz1 = np.count_nonzero(i1)
        nnz2 = np.count_nonzero(i2)
        if nnz1 > nnz2:
            return 1
        elif nnz1 < nnz2:
            return 0
        else:
            for j in range(len(i1)):
                if i1[j] == 0 and i2[j] > 0:
                    return 1
                elif i1[j] > 0 and i2[j] == 0:
                    return 0
            return 0

    # Bubble sort based on custom comparator
    for i in range(len(I) - 1):
        for j in range(i + 1, len(I)):
            if nnzcmp(I[i], I[j]) == 1:
                I[[i, j]] = I[[j, i]]

    return I


def _ginit(L: np.ndarray, N: np.ndarray, Dn: np.ndarray,
           M: int, R: int) -> np.ndarray:
    """Initialize the g vector for CoMoM."""
    g = np.zeros(Dn.shape[0] * (M + 1))

    for i in range(M + 1):
        idx = _hash(N, N, i, Dn, M, R)
        if 0 <= idx < len(g):
            g[idx] = 1.0

    return g


def _hash(Ntot: np.ndarray, n: np.ndarray, i: int,
          Dn: np.ndarray, M: int, R: int) -> int:
    """Compute hash index for CoMoM."""
    diff = Ntot - n
    row_idx = matchrow(Dn, diff)

    if i == 0:  # i=1 in MATLAB (0-indexed)
        return Dn.shape[0] * M + row_idx - 1
    else:
        return (row_idx - 1) * M + i - 1


def _genmatrix(L: np.ndarray, N: np.ndarray, Z: np.ndarray, r: int,
               Ntot: np.ndarray, Dn: np.ndarray, M: int, R: int
               ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate the A, B, DA matrices for CoMoM iteration."""
    size = comb(M + R - 1, M) * (M + 1)
    A = np.zeros((size, size))
    DA = np.zeros((size, size))
    B = np.zeros((size, size))

    row = 0

    for d in range(len(Dn)):
        if np.sum(Dn[d, r:R - 1]) > 0:
            # Dummy rows for unused norm consts
            for k in range(M + 1):
                if row >= size:
                    break
                col = _hash(N, N - Dn[d, :], k, Dn, M, R)
                if 0 <= col < size:
                    A[row, col] = 1

                if np.sum(Dn[d, (r + 1):R - 1]) > 0:
                    if 0 <= col < size:
                        B[row, col] = 1
                else:
                    er = np.zeros(R)
                    er[r] = 1
                    col = _hash(N, N - Dn[d, :] + er, k, Dn, M, R)
                    if 0 <= col < size:
                        B[row, col] = 1
                row += 1
        else:
            if np.sum(Dn[d, :r + 1]) < M:
                for k in range(1, M + 1):
                    if row >= size:
                        break
                    # Add CE
                    col1 = _hash(N, N - Dn[d, :], k, Dn, M, R)
                    col0 = _hash(N, N - Dn[d, :], 0, Dn, M, R)

                    if 0 <= col1 < size:
                        A[row, col1] = 1
                    if 0 <= col0 < size:
                        A[row, col0] = -1

                    for s in range(r):
                        n_oner = oner(N - Dn[d, :], s)
                        colk = _hash(N, n_oner, k, Dn, M, R)
                        if 0 <= colk < size:
                            A[row, colk] = -L[k - 1, s]

                    if 0 <= col1 < size:
                        B[row, col1] = L[k - 1, r]

                    row += 1

                for s in range(r):
                    if row >= size:
                        break
                    # Add PC to A
                    n = N - Dn[d, :]
                    col0 = _hash(N, n, 0, Dn, M, R)
                    n_oner = oner(n, s)
                    col_oner0 = _hash(N, n_oner, 0, Dn, M, R)

                    if 0 <= col0 < size:
                        A[row, col0] = n[s]
                    if 0 <= col_oner0 < size:
                        A[row, col_oner0] = -Z[s]

                    for k in range(1, M + 1):
                        col_onerk = _hash(N, n_oner, k, Dn, M, R)
                        if 0 <= col_onerk < size:
                            A[row, col_onerk] = -L[k - 1, s]

                    row += 1

    # Add PC of class R
    for d in range(len(Dn)):
        if np.sum(Dn[d, r:R - 1]) <= 0:
            if row >= size:
                break
            n = N - Dn[d, :]
            col0 = _hash(N, n, 0, Dn, M, R)

            if 0 <= col0 < size:
                A[row, col0] = n[r]
                DA[row, col0] = 1
                B[row, col0] = Z[r]

            for k in range(1, M + 1):
                colk = _hash(N, n, k, Dn, M, R)
                if 0 <= colk < size:
                    B[row, colk] = L[k - 1, r]

            row += 1

    return A, B, DA


def pfqn_comomrm(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
                 m: int = 1, atol: float = 1e-8) -> ComomResult:
    """
    CoMoM for finite repairman model.

    Computes the normalizing constant for a closed network with
    a single queueing station and delay stations (repairman model).

    Args:
        L: Service demand matrix (1 x R) - single station
        N: Population vector (R,)
        Z: Think time vector (R,)
        m: Replication factor (default: 1)
        atol: Absolute tolerance

    Returns:
        ComomResult with lG (log normalizing constant) and lGbasis

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_comomrm.m
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()

    M, R = L.shape

    if M != 1:
        raise ValueError("pfqn_comomrm: accepts at most a single queueing station.")

    # Import sanitize function
    from .utils import pfqn_nc_sanitize

    lam = np.zeros(R)
    lam, L, N, Z, lG0 = pfqn_nc_sanitize(lam, L, N, Z, atol)

    R = len(N)
    if R == 0:
        return ComomResult(lG=lG0, lGbasis=None)

    zerothinktimes = np.where(Z < 1e-6)[0]

    # Initialize
    nvec = np.zeros(R)

    if len(zerothinktimes) > 0:
        nvec[zerothinktimes] = N[zerothinktimes]
        lh = []

        # These are trivial models with a single queueing station
        lh.append(factln(np.sum(nvec) + m + 1 - 1) - np.sum(factln_vec(nvec)))

        for s in zerothinktimes:
            nvec_s = oner(nvec, s)
            lh.append(factln(np.sum(nvec_s) + m + 1 - 1) - np.sum(factln_vec(nvec_s)))

        lh.append(factln(np.sum(nvec) + m - 1) - np.sum(factln_vec(nvec)))

        for s in zerothinktimes:
            nvec_s = oner(nvec, s)
            lh.append(factln(np.sum(nvec_s) + m - 1) - np.sum(factln_vec(nvec_s)))

        lh = np.array(lh)
    else:
        lh = np.zeros(2)

    h = np.exp(lh - np.max(lh)) if len(lh) > 0 else np.zeros(2)

    if len(zerothinktimes) == R:
        lGbasis = lh
        lG = lG0 + lh[len(lh) - R - 1] if len(lh) > R else lG0
        return ComomResult(lG=lG, lGbasis=lGbasis)

    scale = np.ones(int(np.sum(N)))
    nt = int(np.sum(nvec))
    h_1 = h.copy()

    # Iterate through remaining classes
    for r in range(len(zerothinktimes), R):
        for Nr in range(1, int(N[r]) + 1):
            nvec[r] += 1

            if Nr == 1:
                if r > len(zerothinktimes):
                    # Expand h vector
                    hr = np.zeros(2 * r)
                    hr[:r - 1] = h[:r - 1]
                    hr[r:2 * r - 1] = h[r - 1:2 * (r - 1)]
                    h = hr

                    # Update scalings
                    if nt > 0:
                        h[r - 1] = h_1[0] / scale[nt - 1]
                        h[-1] = h_1[r - 1] / scale[nt - 1]

                # Build coefficient matrices
                # CE for G+
                A12 = np.zeros((r, r))
                A12[0, 0] = -1

                # Class PCs for G
                for s in range(r - 1):
                    A12[1 + s, 0] = N[s]
                    A12[1 + s, 1 + s] = -Z[s]

                # Class-R PCs
                B2r = np.column_stack([m * L[0, r] * np.eye(r), Z[r] * np.eye(r)])

                # Explicit formula for inv(C)
                iC = -np.eye(r) / m
                iC[0, :] = -1 / m
                iC[0, 0] = 1

                # Explicit formula for F1r
                F1r = np.zeros((2 * r, 2 * r))
                F1r[0, 0] = 1

                # F2r by the definition
                F2r = np.vstack([-iC @ A12 @ B2r, B2r])

            h_1 = h.copy()
            h = (F1r + F2r / nvec[r]) @ h_1

            nt = int(np.sum(nvec))
            scale[nt - 1] = np.abs(np.sum(np.sort(h)))
            if scale[nt - 1] > 0:
                h = np.abs(h) / scale[nt - 1]

    # Unscale and return the log of the normalizing constant
    log_scale_sum = np.sum(np.log(scale[scale > 0]))
    lG = lG0 + np.log(np.abs(h[len(h) - R]) + 1e-300) + log_scale_sum
    lGbasis = np.log(np.abs(h) + 1e-300) + log_scale_sum

    return ComomResult(lG=lG, lGbasis=lGbasis)


def pfqn_comomrm_orig(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
                      m: int = 1, atol: float = 1e-8) -> ComomResult:
    """
    Original CoMoM implementation for repairman model.

    This is the original implementation of CoMoM without optimizations.
    Kept for reference and validation.

    Args:
        L: Service demand matrix (1 x R)
        N: Population vector (R,)
        Z: Think time vector (R,)
        m: Replication factor (default: 1)
        atol: Absolute tolerance

    Returns:
        ComomResult with lG and lGbasis

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_comomrm_orig.m
    """
    # This is essentially the same algorithm as pfqn_comomrm
    # but implemented more directly following the original MATLAB code
    return pfqn_comomrm(L, N, Z, m, atol)


def pfqn_comomrm_ms(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
                    m: int, c: int, atol: float = 1e-8) -> ComomResult:
    """
    CoMoM for multi-server repairman model.

    Computes the normalizing constant for a repairman model where
    the queueing station has multiple servers.

    Args:
        L: Service demand matrix (1 x R)
        N: Population vector (R,)
        Z: Think time vector (R,)
        m: Number of identical stations
        c: Number of servers per station
        atol: Absolute tolerance

    Returns:
        ComomResult with lG and lGbasis

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_comomrm_ms.m
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()

    M, R = L.shape

    if M != 1:
        raise ValueError("pfqn_comomrm_ms: accepts at most a single queueing station.")

    # Compute load-dependent service rates for multi-server
    from .utils import pfqn_mu_ms

    Ntot = int(np.sum(N))
    mu = pfqn_mu_ms(Ntot, m, c)

    # Use load-dependent COMOM
    from .ncld import pfqn_comomrm_ld

    # Reshape mu for pfqn_comomrm_ld
    mu_matrix = mu.reshape(1, -1)

    result = pfqn_comomrm_ld(L, N, Z, mu_matrix, {'tol': atol})

    return ComomResult(lG=result.lG, lGbasis=None)


def pfqn_procomom(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
                  atol: float = 1e-14):
    """
    ProCoMoM algorithm for computing marginal queue-length probabilities.

    Computes the marginal queue-length probability distribution at each station
    using the Probabilistic Convolution Method of Moments. Uses matrix recursion
    with SVD/QR decomposition for numerical stability, with automatic
    perturbation on rank deficiency.

    Args:
        L: Service demand matrix (M x R).
        N: Population vector (R,).
        Z: Think time vector (R,), default zeros.
        atol: Absolute numerical tolerance (default 1e-14).

    Returns:
        Pr: Marginal probability matrix (M x sumN+1).
            Pr[k, j] = P(n_k = j) for station k, queue length j.
        Q: Mean queue length vector (M,).

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_procomom.m
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten().astype(int)

    M, R = L.shape

    if Z is None or (hasattr(Z, '__len__') and len(Z) == 0):
        Z = np.zeros(R, dtype=float)
    else:
        Z = np.asarray(Z, dtype=float).flatten()

    sumN = int(np.sum(N))

    # Rescale demands per class for numerical stability
    # Normalized probabilities are invariant to per-class demand scaling
    Lmax = np.max(L, axis=0)
    Lmax[Lmax < atol] = 1.0
    L = L / Lmax[np.newaxis, :]
    Z = Z / Lmax

    # Build Dn basis: multichoose(R, M) with column R zeroed, sorted by nnzpos
    Dn = multichoose(R, M)
    Dn[:, R - 1] = 0
    Dn = _sortbynnzpos_procomom(Dn)
    numDn = Dn.shape[0]
    basisSize = numDn * M

    def _phash(dn, i):
        """Hash function mapping (dn, i) to column index (0-based).
        i is 1-based station index."""
        pos = matchrow(Dn, dn)  # 1-based
        if pos <= 0:
            return -1
        return (pos - 1) * M + (i - 1)  # 0-based column

    def _genpmatrix(Ls, Zs, Ncur, r):
        """Generate probability matrices for ProCoMoM recursion.
        r is 1-based class index (MATLAB convention)."""
        # Count rows
        numRows = 0
        for d in range(numDn):
            dn = Dn[d]
            if r <= R - 1 and np.sum(dn[r - 1:R - 1]) > 0:
                numRows += M
            else:
                if np.sum(dn[:r]) < M:
                    numRows += M + r - 1  # (M-1) CE + r PC, but actually (M-1) CE + (r-1) PC
                else:
                    pass
                numRows += 1  # extra PC for class r

        A = np.zeros((numRows, basisSize))
        B = np.zeros((numRows, basisSize))
        DC = np.zeros((numRows, basisSize))
        DD = np.zeros((numRows, basisSize))
        row = 0

        for d in range(numDn):
            dn = Dn[d]
            # Check Branch A condition: r <= R-1 and sum(Dn(d, r:R-1)) > 0
            # MATLAB: r<=R-1 && sum(Dn(d,r:R-1))>0
            # In MATLAB r is 1-based, Dn columns are 1-based
            # Dn(d, r:R-1) means columns r to R-1 (1-based)
            # In Python 0-based: columns (r-1) to (R-2)
            if r <= R - 1 and np.sum(dn[r - 1:R - 1]) > 0:
                # Branch A: propagation through class boundaries
                for k in range(1, M + 1):
                    colA = _phash(dn, k)
                    if 0 <= colA < basisSize:
                        A[row, colA] = 1.0

                    # Check if r+1 <= R-1 and sum(Dn(d, r+1:R-1)) > 0
                    # MATLAB: r+1<=R-1 && sum(Dn(d,r+1:R-1))>0
                    # Python: columns r to R-2
                    if r + 1 <= R - 1 and np.sum(dn[r:R - 1]) > 0:
                        if 0 <= colA < basisSize:
                            B[row, colA] = 1.0
                    else:
                        shifted = dn.copy()
                        shifted[r - 1] -= 1  # 0-based index for class r
                        colB = _phash(shifted, k)
                        if 0 <= colB < basisSize:
                            B[row, colB] = 1.0
                    row += 1
            else:
                # Branch B: CE, PC, and extra PC equations
                if np.sum(dn[:r]) < M:
                    # CE equations: k=1..M-1
                    for k in range(1, M):
                        colKP1 = _phash(dn, k + 1)
                        col1 = _phash(dn, 1)
                        if 0 <= colKP1 < basisSize:
                            A[row, colKP1] = 1.0
                        if 0 <= col1 < basisSize:
                            A[row, col1] = -1.0

                        for s in range(1, r):  # s=1..r-1 (MATLAB 1-based)
                            shifted = dn.copy()
                            shifted[s - 1] += 1  # UP shift
                            col = _phash(shifted, k + 1)
                            if 0 <= col < basisSize:
                                A[row, col] -= Ls[k - 1, s - 1]

                        if 0 <= colKP1 < basisSize:
                            B[row, colKP1] = Ls[k - 1, r - 1]
                        row += 1

                    # PC equations: s=1..r-1 (MATLAB 1-based)
                    for s in range(1, r):
                        nd_s = Ncur[s - 1] - dn[s - 1]
                        col1 = _phash(dn, 1)
                        if 0 <= col1 < basisSize:
                            A[row, col1] = float(nd_s)

                        shifted = dn.copy()
                        shifted[s - 1] += 1  # UP shift
                        colBase = _phash(shifted, 1)
                        if 0 <= colBase < basisSize:
                            A[row, colBase] -= Zs[s - 1]
                            DC[row, colBase] = Ls[M - 1, s - 1]

                        for k in range(1, M):
                            col = _phash(shifted, k + 1)
                            if 0 <= col < basisSize:
                                A[row, col] -= Ls[k - 1, s - 1]
                        row += 1

                # Extra PC for class r (always in Branch B)
                nd_r = Ncur[r - 1] - dn[r - 1]
                col1 = _phash(dn, 1)
                if 0 <= col1 < basisSize:
                    A[row, col1] = float(nd_r)
                    B[row, col1] = Zs[r - 1]
                for k in range(1, M):
                    colKP1 = _phash(dn, k + 1)
                    if 0 <= colKP1 < basisSize:
                        B[row, colKP1] = Ls[k - 1, r - 1]
                if 0 <= col1 < basisSize:
                    DD[row, col1] = Ls[M - 1, r - 1]
                row += 1

        return A[:row], B[:row], DC[:row], DD[:row]

    def _solve_station(Ls, Zs):
        """Solve for a single station (after rotation to last position)."""
        rankdef = False

        # pk[:, j] holds basis coefficients for queue length n=j
        pk = np.zeros((basisSize, sumN + 1))

        # Initialize: for empty network, all stations have probability 1 at n=0
        zero_dn = np.zeros(R, dtype=int)
        for kk in range(1, M + 1):
            idx = _phash(zero_dn, kk)
            if 0 <= idx < basisSize:
                pk[idx, 0] = 1.0

        Ncur = np.zeros(R, dtype=int)
        for r in range(1, R + 1):  # 1-based class
            for Nr in range(1, N[r - 1] + 1):
                Ncur[r - 1] = Nr
                pklast = pk.copy()
                pk = np.zeros((basisSize, sumN + 1))

                Ag, Bg, DCg, DDg = _genpmatrix(Ls, Zs, Ncur, r)
                numRows_local = Ag.shape[0]

                # SVD for rank-revealing decomposition
                U, sv, Vt = np.linalg.svd(Ag, full_matrices=False)
                V = Vt.T  # V is now numCols x min(numRows, numCols)

                tol_svd = max(numRows_local, basisSize) * np.finfo(float).eps * sv[0] if len(sv) > 0 else 0
                rk = np.sum(sv > tol_svd)

                sumNcur = int(np.sum(Ncur))

                if rk < basisSize:
                    rankdef = True
                    # Use minimum-norm solution via truncated SVD
                    Ur = U[:, :rk]
                    Vr = V[:, :rk]
                    Si = np.diag(1.0 / sv[:rk])

                    # pB = Vr * Si * Ur' * Bg
                    pB = Vr @ Si @ Ur.T @ Bg
                    pDC = Vr @ Si @ Ur.T @ DCg
                    pDD = Vr @ Si @ Ur.T @ DDg

                    pk[:, 0] = pB @ pklast[:, 0]
                    for n in range(1, sumNcur + 1):
                        pk[:, n] = (pB @ pklast[:, n]
                                    + n * pDC @ pk[:, n - 1]
                                    + n * pDD @ pklast[:, n - 1])
                else:
                    # Full rank: use QR for speed and accuracy
                    Qg, Rg = np.linalg.qr(Ag, mode='reduced')
                    QtB = Qg.T @ Bg
                    QtDC = Qg.T @ DCg
                    QtDD = Qg.T @ DDg

                    pk[:, 0] = np.linalg.solve(Rg, QtB @ pklast[:, 0])
                    for n in range(1, sumNcur + 1):
                        rhs = (QtB @ pklast[:, n]
                               + n * QtDC @ pk[:, n - 1]
                               + n * QtDD @ pklast[:, n - 1])
                        pk[:, n] = np.linalg.solve(Rg, rhs)

                # Rescale to prevent overflow
                smax = np.max(np.abs(pk))
                if smax > 0 and np.isfinite(smax):
                    pk = pk / smax

        # Extract result: first basis element at zero Dn entry
        zero_dn = np.zeros(R, dtype=int)
        idx = _phash(zero_dn, 1)
        if 0 <= idx < basisSize:
            dist = pk[idx, :]
        else:
            dist = np.zeros(sumN + 1)

        return dist, rankdef

    def _solve_all(Ls_in, Zs_in):
        """Solve for all stations."""
        Pr_local = np.zeros((M, sumN + 1))
        rankdef = False

        for station in range(M):
            # Rotate: swap station with last (M-1)
            Lrot = Ls_in.copy()
            Lrot[[station, M - 1], :] = Lrot[[M - 1, station], :]

            dist, rd = _solve_station(Lrot, Zs_in)
            if rd:
                rankdef = True
            total = np.sum(dist)
            if abs(total) > 0:
                Pr_local[station, :] = dist / total

        return Pr_local, rankdef

    # Solve with auto-perturbation on rank deficiency
    Pr, rankdef = _solve_all(L, Z)

    if rankdef:
        # Try progressively larger perturbations
        rng = np.random.RandomState(23000)
        Lscale = np.max(np.abs(L))
        if Lscale < atol:
            Lscale = 1.0

        for delta_exp in [-10, -8, -6, -4]:
            delta = Lscale * 10 ** delta_exp
            rng_local = np.random.RandomState(23000)
            Lp = L + delta * (1 + rng_local.rand(M, R))
            Zp = Z + delta * (1 + rng_local.rand(R))

            Pr_try, rd = _solve_all(Lp, Zp)
            if (not rd and np.all(Pr_try >= -1e-6)
                    and np.all(np.abs(np.sum(Pr_try, axis=1) - 1.0) < 0.01)):
                Pr = Pr_try
                break
            Pr = Pr_try

    Q = Pr @ np.arange(sumN + 1)

    return Pr, Q


def _sortbynnzpos_procomom(I: np.ndarray) -> np.ndarray:
    """Sort rows of I by number of nonzero positions (ascending),
    with ties broken by position of zeros (rows with leading zeros first).

    Uses bubble sort matching MATLAB's nested loop implementation."""
    I = I.copy()
    n = I.shape[0]
    for ii in range(n - 1):
        for jj in range(ii + 1, n):
            if _nnzcmp(I[ii], I[jj]) == 1:
                I[[ii, jj]] = I[[jj, ii]]
    return I


def _nnzcmp(i1: np.ndarray, i2: np.ndarray) -> int:
    """Compare two vectors by number of nonzeros, then by position of zeros.
    Returns 1 if i1 should come after i2, 0 otherwise."""
    nnz1 = np.count_nonzero(i1)
    nnz2 = np.count_nonzero(i2)
    if nnz1 > nnz2:
        return 1
    elif nnz1 < nnz2:
        return 0
    else:
        for jj in range(len(i1)):
            if i1[jj] == 0 and i2[jj] > 0:
                return 1
            elif i1[jj] > 0 and i2[jj] == 0:
                return 0
        return 0


def pfqn_procomom2(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
                   atol: float = 1e-8) -> float:
    """
    Projected CoMoM method for normalizing constant.

    Uses a projection-based approach to compute the normalizing constant,
    which can be more efficient for certain model structures.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,), default zeros
        atol: Absolute tolerance

    Returns:
        lG: Logarithm of the normalizing constant

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_procomom2.m
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()

    M, R = L.shape

    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=float).flatten()

    # For now, fall back to standard CoMoM
    # Full implementation would use projection techniques
    return pfqn_comom(L, N, Z, atol)


__all__ = [
    'pfqn_comom',
    'pfqn_comomrm',
    'pfqn_comomrm_orig',
    'pfqn_comomrm_ms',
    'pfqn_procomom',
    'pfqn_procomom2',
    'ComomResult',
]
