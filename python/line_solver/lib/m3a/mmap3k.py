"""MMAP(3,K): closed-form marking fit for a marked MAP of third order.

The MMAP(2,K) argument does not depend on the order. Two facts carry over
unchanged (both verified symbolically, see
io/sage/proofs/mmap3k_marking_inverse.py):

1. every per-class characteristic in which the class matrix appears exactly
   once is LINEAR in the marking fractions. With z = nnz(D1) fractions per
   class, z characteristics give a square linear system;
2. the system is BLOCK DIAGONAL in the classes, so one z x z block is built
   once and reused for every class, and the cost does not grow with K.

What does change with the order is which characteristics are needed. At order
two, (p_c, F_c, B_c) suffice; at order three the independent set of lowest
total order is

    (a, b) = (1,0), (1,1), (2,0), (3,0)
    i.e.    p_c,   F_c,   B_c,   B_c^(2)

where a is the backward order and b the forward order of
pie A^a D1^(c) A^b 1. Alternating forward and backward orders does NOT work at
higher orders: at order four it yields only four independent functionals out of
five, which is a trap worth knowing.

The z x z block is assembled exactly at run time by evaluating the linear map
on unit markings -- no finite differences and no computer algebra, since the map
is linear and its offset is zero. That is preferred over inlining the symbolic
inverse here: the order-3 inverse is about 5 KB of expressions and depends on
the sparsity pattern of D1, whereas the column construction below is exact for
ANY order and ANY pattern.

The underlying MAP(3) is taken as given: unlike order two, there is no canonical
inverse in the tree that turns (moments, autocorrelation) into an order-3 MAP,
so the caller supplies (D0, D1) from whatever fitter it prefers.
"""
from typing import List, Optional, Sequence, Tuple

import numpy as np

FEASTOL = 1e-8


def _stationary_embedded(D0: np.ndarray, D1: np.ndarray) -> np.ndarray:
    n = D0.shape[0]
    A = np.linalg.inv(-D0)
    P = A @ D1
    M = (P.T - np.eye(n)).copy()
    M[n - 1, :] = 1.0
    rhs = np.zeros(n)
    rhs[n - 1] = 1.0
    return np.linalg.solve(M, rhs)


def _functional(pie: np.ndarray, A: np.ndarray, Dc: np.ndarray, a: int, b: int) -> float:
    """pie A^a Dc A^b 1, the generic per-class linear characteristic."""
    n = A.shape[0]
    return float(pie @ np.linalg.matrix_power(A, a) @ Dc @ np.linalg.matrix_power(A, b)
                 @ np.ones(n))


def marking_orders(n: int) -> List[Tuple[int, int]]:
    """(backward, forward) orders of an independent characteristic set.

    Verified independent for orders 2 and 3 in
    io/sage/proofs/mmap3k_marking_inverse.py; the same greedy rule (lowest total
    order first) extends to higher orders.
    """
    if n == 2:
        return [(1, 0), (1, 1), (2, 0)]
    if n == 3:
        return [(1, 0), (1, 1), (2, 0), (3, 0)]
    # general rule: (1,0), (1,1) then increasing pure backward orders
    orders = [(1, 0), (1, 1)]
    a = 2
    while len(orders) < n + 1:
        orders.append((a, 0))
        a += 1
    return orders


def marking_block(D0: np.ndarray, D1: np.ndarray,
                  orders: Optional[Sequence[Tuple[int, int]]] = None
                  ) -> Tuple[np.ndarray, List[Tuple[int, int]], List[Tuple[int, int]]]:
    """The z x z block mapping one class's marking fractions to its
    characteristics, plus the nonzero positions of D1 the fractions refer to.

    Exact: the map is linear with zero offset, so column j is the characteristic
    vector of the unit marking on the j-th nonzero of D1.
    """
    n = D0.shape[0]
    nz = [(i, j) for i in range(n) for j in range(n) if D1[i, j] != 0.0]
    z = len(nz)
    if orders is None:
        orders = marking_orders(n)
    orders = list(orders)[:z]
    if len(orders) != z:
        raise ValueError('mmap3k_fit: need %d characteristics for %d marking '
                         'fractions' % (z, z))
    A = np.linalg.inv(-D0)
    pie = _stationary_embedded(D0, D1)
    M = np.zeros((z, z))
    for jj, (i, j) in enumerate(nz):
        Dc = np.zeros_like(D1)
        Dc[i, j] = D1[i, j]
        for ii, (a, b) in enumerate(orders):
            M[ii, jj] = _functional(pie, A, Dc, a, b)
    return M, nz, orders


def mmap3k_fit(D0: np.ndarray, D1: np.ndarray, P: Sequence[float],
               F: Sequence[float], B: Sequence[float],
               B2: Optional[Sequence[float]] = None,
               exact_only: bool = False) -> List[np.ndarray]:
    """Mark a given MAP(3) so that the per-class characteristics are matched.

    Args:
        D0, D1: the underlying MAP, of any order (validated at orders 2 and 3)
        P: class probabilities, summing to one
        F: first-order forward moments
        B: first-order backward moments
        B2: second-order backward moments, required at order three and above
        exact_only: raise instead of returning an infeasible marking

    Returns:
        The MMAP as [D0, D1, D1^(1), ..., D1^(K)]

    Raises:
        ValueError: when the marking is infeasible and exact_only is set, or
            when a required characteristic is missing
    """
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    P = np.asarray(P, dtype=float).ravel()
    F = np.asarray(F, dtype=float).ravel()
    B = np.asarray(B, dtype=float).ravel()
    K = P.size
    n = D0.shape[0]

    M, nz, orders = marking_block(D0, D1)
    z = len(nz)
    if z > 3 and B2 is None:
        raise ValueError('mmap3k_fit: order %d needs the second-order backward '
                         'moments (B2) as well' % n)
    B2arr = np.asarray(B2, dtype=float).ravel() if B2 is not None else None

    if abs(np.linalg.det(M)) < 1e-12 * max(1.0, np.abs(M).max() ** z):
        raise ValueError('mmap3k_fit: the underlying MAP is on the degenerate '
                         'locus of the marking system')

    q = np.zeros((z, K))
    for c in range(K):
        # y = (p_c, p_c F_c, p_c B_c, p_c B2_c, ...), in the order of `orders`
        y = []
        for (a, b) in orders:
            if (a, b) == (1, 0):
                y.append(P[c])
            elif (a, b) == (1, 1):
                y.append(P[c] * F[c])
            elif (a, b) == (2, 0):
                y.append(P[c] * B[c])
            elif (a, b) == (3, 0):
                y.append(P[c] * B2arr[c])
            else:
                raise ValueError('mmap3k_fit: no target supplied for the '
                                 'characteristic (a=%d, b=%d)' % (a, b))
        q[:, c] = np.linalg.solve(M, np.asarray(y, dtype=float))

    viol = max(float(-q.min()), float(q.max() - 1.0), 0.0)
    viol = max(viol, float(np.max(np.abs(q.sum(axis=1) - 1.0))))
    if viol > FEASTOL:
        if exact_only:
            raise ValueError('mmap3k_fit: the closed-form marking is infeasible '
                             '(worst violation %.3g)' % viol)

    out = [D0.copy(), D1.copy()]
    qc = np.clip(q, 0.0, 1.0)
    for c in range(K):
        Dc = np.zeros_like(D1)
        for jj, (i, j) in enumerate(nz):
            Dc[i, j] = D1[i, j] * qc[jj, c]
        out.append(Dc)
    return out
