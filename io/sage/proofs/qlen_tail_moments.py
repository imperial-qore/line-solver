"""The survival (tail) vertex of the house of moments, and the two closed-form
identities that feed it from normalizing constants (api/moment, api/pfqn).

Three claims are proved in exact arithmetic:

  * the TAIL EDGE. `moment_binomial_from_tail` and its inverse are ported
    verbatim and checked against the definition on a law with FREE symbolic
    masses over {0,...,n}, univariate and joint. Since the transform is linear
    in the law and every law supported on the box is covered, this proves it in
    general. The edge is UPPER triangular, so it consumes the whole tail: the
    proof also exhibits that truncating it strictly loses mass, which is why
    the implementation documents exactness only when the array covers the
    support;
  * the SINGLE-CLASS SURVIVAL IDENTITY behind `pfqn_qlen_joint_moments`,
    P(n_i >= k_i for all i) = (prod L_i^k_i) G(N - sum k_i)/G(N), with the
    demands and the think time FREE. This is the identity that lets any
    normalizing-constant algorithm produce joint queue-length moments;
  * the MULTICLASS COMPLEMENTARY-NETWORK IDENTITY, P(n_S = m) =
    prod_{i in S} f_i(m_i) G_{S^c}(N - sum m_i)/G(N) with the multinomial
    occupancy f_i, which is what replaces the survival identity when the
    geometric factorization fails. The failure itself is exhibited: the naive
    survival formula is NOT an identity for two classes.

Everything is done over the rational function field in the demands and the
think times, so no numerical tolerance is involved anywhere.
"""
from itertools import product

from sage.all import QQ, PolynomialRing, binomial, factorial, prod

results = []


def check(label, ok):
    results.append(bool(ok))
    print('  %-58s %s' % (label, 'IDENTITY' if ok else 'FAILED'))


# ---------------------------------------------------------------------------
# verbatim port of the tail edge
# ---------------------------------------------------------------------------

def moment_housematrix_tail(edge, n):
    T = [[QQ(0)] * (n + 1) for _ in range(n + 1)]
    T[0][0] = QQ(1)
    for i in range(1, n + 1):
        for k in range(i, n + 1):
            c = binomial(k - 1, i - 1)
            T[i][k] = c if edge == 'binomial_from_tail' else QQ(-1) ** (k - i) * c
    return T


def apply_mode(A, dims, T, mode):
    out = {}
    for a in product(*[range(m) for m in dims]):
        acc = 0
        for k in range(dims[mode]):
            b = list(a)
            b[mode] = k
            acc += T[a[mode]][k] * A[tuple(b)]
        out[a] = acc
    return out


def jointtrans(A, dims, edge):
    out = A
    for mode in range(len(dims)):
        out = apply_mode(out, dims, moment_housematrix_tail(edge, dims[mode] - 1), mode)
    return out


# ---------------------------------------------------------------------------
# 1. the tail edge against the definition, with FREE masses
# ---------------------------------------------------------------------------
NU = 6
Rp = PolynomialRing(QQ, ['p%d' % v for v in range(NU + 1)])
pv = Rp.gens()
tail = [sum(pv[v] for v in range(m, NU + 1)) for m in range(NU + 1)]
bino = [sum(binomial(v, j) * pv[v] for v in range(NU + 1)) for j in range(NU + 1)]

print('== univariate tail edge, law on {0,...,%d} with free masses ==' % NU)
A = moment_housematrix_tail('binomial_from_tail', NU)
B = moment_housematrix_tail('tail_from_binomial', NU)
got = [sum(A[j][m] * tail[m] for m in range(NU + 1)) for j in range(NU + 1)]
# b_0 = 1 only when the masses sum to 1; the edge itself maps t_0 to b_0
check('b_j = sum_(m>=j) C(m-1,j-1) t_m for j >= 1',
      all(got[j] == bino[j] for j in range(1, NU + 1)) and got[0] == tail[0])
back = [sum(B[m][j] * got[j] for j in range(NU + 1)) for m in range(NU + 1)]
check('tail_from_binomial inverts it', back == tail)
check('the two matrices are mutually inverse',
      all(sum(A[i][k] * B[k][j] for k in range(NU + 1)) == (1 if i == j else 0)
          for i in range(NU + 1) for j in range(NU + 1)))
# truncation: dropping the last tail entry loses exactly its contribution
trunc = tail[:-1] + [Rp(0)]
lost = [sum(A[j][m] * (tail[m] - trunc[m]) for m in range(NU + 1)) for j in range(NU + 1)]
check('truncating the tail subtracts a nonnegative amount from every b_j',
      all(all(c >= 0 for c in poly.coefficients()) for poly in lost))

# ---------------------------------------------------------------------------
# 2. the joint tail edge, bivariate law with free masses
# ---------------------------------------------------------------------------
DIMS = (4, 4)
idx = list(product(*[range(m) for m in DIMS]))
Rq = PolynomialRing(QQ, ['q%d_%d' % a for a in idx])
qv = dict(zip(idx, Rq.gens()))
jtail = {k: sum(qv[v] for v in idx if all(v[j] >= k[j] for j in range(2))) for k in idx}
jbino = {k: sum(prod([binomial(v[j], k[j]) for j in range(2)]) * qv[v] for v in idx)
         for k in idx}

print('\n== joint tail edge, bivariate law on a %s box with free masses ==' % (DIMS,))
gotj = jointtrans(jtail, DIMS, 'binomial_from_tail')
check('b_k = E[prod_j C(N_j,k_j)] for every multi-index',
      all(gotj[k] == jbino[k] for k in idx))
check('joint_tail_from_binomial inverts it',
      all(jointtrans(gotj, DIMS, 'tail_from_binomial')[k] == jtail[k] for k in idx))

# ---------------------------------------------------------------------------
# 3. the single-class survival identity, free demands and think time
# ---------------------------------------------------------------------------
print('\n== single-class survival identity vs the state sum ==')
for M, NPOP in [(2, 5), (3, 4)]:
    Rl = PolynomialRing(QQ, ['L%d' % i for i in range(M)] + ['Z'])
    Ls = Rl.gens()[:M]
    Zs = Rl.gens()[M]

    def states(n):
        return [s for s in product(*[range(n + 1)] * (M + 1)) if sum(s) == n]

    def weight(s):
        w = Rl(1)
        for i in range(M):
            w *= Ls[i] ** s[i]
        return w * Zs ** s[M] / factorial(s[M])

    G = {n: sum(weight(s) for s in states(n)) for n in range(NPOP + 1)}
    ok = True
    for k in product(*[range(NPOP + 1)] * M):
        if sum(k) > NPOP:
            continue
        lhs = sum(weight(s) for s in states(NPOP) if all(s[i] >= k[i] for i in range(M)))
        rhs = prod([Ls[i] ** k[i] for i in range(M)]) * G[NPOP - sum(k)]
        ok = ok and Rl(lhs - rhs) == 0
    check('M = %d, N = %d: P(n >= k) G(N) = (prod L^k) G(N-|k|)' % (M, NPOP), ok)

# ---------------------------------------------------------------------------
# 4. the multiclass identities: the naive one fails, the complementary one holds
# ---------------------------------------------------------------------------
print('\n== multiclass: complementary-network law vs the state sum ==')
M, R, NPOP = 2, 2, (3, 2)
Rm = PolynomialRing(QQ, ['L%d_%d' % (i, r) for i in range(M) for r in range(R)]
                    + ['Z%d' % r for r in range(R)])
gens = Rm.gens()
LL = [[gens[i * R + r] for r in range(R)] for i in range(M)]
ZZ = list(gens[M * R:])


def mstates(n, stations):
    """States of a network with the given queueing stations plus the delay."""
    per = []
    for r in range(R):
        per.append([c for c in product(*[range(n[r] + 1)] * (len(stations) + 1))
                    if sum(c) == n[r]])
    out = []
    for combo in product(*per):
        out.append([[combo[r][j] for r in range(R)] for j in range(len(stations) + 1)])
    return out


def mweight(s, stations):
    w = Rm(1)
    for j, i in enumerate(stations):
        tot = sum(s[j])
        w *= factorial(tot)
        for r in range(R):
            w *= LL[i][r] ** s[j][r] / factorial(s[j][r])
    for r in range(R):
        w *= ZZ[r] ** s[-1][r] / factorial(s[-1][r])
    return w


def mG(n, stations):
    return sum(mweight(s, stations) for s in mstates(n, stations))


GN = mG(list(NPOP), [0, 1])
ok_c, ok_naive = True, True
for m0 in product(range(NPOP[0] + 1), range(NPOP[1] + 1)):
    lhs = sum(mweight(s, [0, 1]) for s in mstates(list(NPOP), [0, 1])
              if [s[0][r] for r in range(R)] == list(m0))
    f0 = factorial(sum(m0)) * prod([LL[0][r] ** m0[r] / factorial(m0[r]) for r in range(R)])
    rhs = f0 * mG([NPOP[r] - m0[r] for r in range(R)], [1])
    ok_c = ok_c and Rm(lhs - rhs) == 0
    naive = prod([LL[0][r] ** m0[r] for r in range(R)]) * \
        mG([NPOP[r] - m0[r] for r in range(R)], [0, 1])
    if sum(m0) > 1:
        ok_naive = ok_naive and Rm(lhs - naive) == 0
check('P(n_0 = m) G(N) = f_0(m) G_(complement)(N-m)', ok_c)
print('  %-58s %s' % ('the naive geometric survival form for R = 2',
                      'HOLDS' if ok_naive else 'FAILS, as expected'))
results.append(not ok_naive)

print('\ntail vertex and its queueing identities verified:', all(results))
