"""
PAS_CYCLIC_NORMCONST_BRUTEFORCE  Brute-force normalizing constant of a closed
cyclic network of two pass-and-swap (PAS) queues.

  Topology:   --> Queue1 --> Queue2 -->   (cyclic, R closed classes)

A pass-and-swap queue (Dorsman & Gardner 2024, Queueing Syst. 107:205-256) is an
order-independent (OI) station: its total service rate is a function mu_i(c) of
the ordered class sequence c at the station, and the stationary law is
product-form and INVARIANT to the swap graph. Hence the swap graph is irrelevant
to the normalizing constant and is not needed here.

For a closed network of OI/PAS stations the unnormalized stationary weight of a
global ordered state (c1,c2) factorizes as

      w(c1,c2) = Phi_1(c1) * Phi_2(c2),
      Phi_i(c) = prod_{j=1..n} e_{i,c_j} / mu_i(c_1,...,c_j),

where e_{i,r} is the per-class visit ratio (relative arrival rate at station i)
and mu_i(prefix) is the total service rate of the ordered prefix. The
normalizing constant is the exact sum over the whole state space,

      G = sum_{(c1,c2)} Phi_1(c1) * Phi_2(c2)
        = sum_{m<=K} S_1(m) * S_2(K-m),

obtained by (i) splitting the population vector K into the per-class counts m at
station 1 and K-m at station 2, and (ii) summing the OI balance function over
every ordered arrangement (multiset permutation) of each station's count vector.
S_i(n) is that per-station sum.

Returns the scalar normalizing constant G.
"""

import itertools


def oi_sum(counts, e, mu):
    """Sum of the OI balance function over every ordered arrangement (distinct
    multiset permutation) of the per-class multiset `counts`.

    Identical-class jobs yield a single ordered state, so the recursion branches
    on the next CLASS placed, generating each distinct sequence exactly once.
    """
    return _oi_dfs(list(counts), [], 1.0, e, mu)


def _oi_dfs(counts, prefix, w, e, mu):
    if all(c == 0 for c in counts):
        return w                       # completed an ordered state
    total = 0.0
    for r in range(len(counts)):
        if counts[r] > 0:
            c2 = list(counts)
            c2[r] -= 1
            p2 = prefix + [r]
            total += _oi_dfs(c2, p2, w * e[r] / mu(p2), e, mu)
    return total


def pas_rate(c, beta, k):
    """Total OI service rate of the ordered prefix c (0-based class ids) for a
    k-server PAS station: the base rates of the first min(n,k) positions."""
    if not c:
        return 0.0
    return sum(beta[j] for j in c[:min(len(c), k)])


def enum_splits(K):
    """Every per-class count vector m with 0 <= m_r <= K_r."""
    return [list(m) for m in itertools.product(*[range(k + 1) for k in K])]


def norm_const(K, e1, e2, mu1, mu2):
    """Normalizing constant for an arbitrary population K (used for X = G/G)."""
    G = 0.0
    for m in enum_splits(K):
        G += oi_sum(m, e1, mu1) * oi_sum([K[r] - m[r] for r in range(len(K))],
                                         e2, mu2)
    return G


def selftest_single_class():
    """Single class, two single-server PAS stations reduces to a closed
    two-queue Gordon-Newell network whose normalizing constant is the geometric
    sum G = sum_{n=0..N} (1/mu1)^n (1/mu2)^(N-n). Verify the brute-force
    enumeration reproduces it."""
    N = 5
    mu1, mu2 = 3.0, 2.0
    e1 = [1.0]
    e2 = [1.0]
    muf1 = lambda c: pas_rate(c, [mu1], 1)
    muf2 = lambda c: pas_rate(c, [mu2], 1)
    Gbf = norm_const([N], e1, e2, muf1, muf2)
    Gcf = sum((1 / mu1) ** n * (1 / mu2) ** (N - n) for n in range(N + 1))
    assert abs(Gbf - Gcf) <= 1e-12 * Gcf, \
        'self-test failed: brute-force %.15g vs closed-form %.15g' % (Gbf, Gcf)
    print('Self-test (single-class Gordon-Newell): PASS (G=%.12g)' % Gbf)


def pas_cyclic_normconst_bruteforce():
    # ---- model parameters ---------------------------------------------
    K = [2, 2]                  # per-class closed populations (R classes)
    R = len(K)

    # Visit ratios. Cyclic two-queue with no class switching => every class
    # visits each station once per cycle, so e1 = e2 = 1 (any common scaling
    # of e cancels in normalized metrics).
    e1 = [1.0] * R
    e2 = [1.0] * R

    # Total OI service-rate functions mu_i(c) of the ordered prefix c (list of
    # 0-based class ids). ORDER-INDEPENDENCE requires the TOTAL rate to depend
    # only on the customer multiset; a single-server class-dependent head rate
    # would violate this and is NOT a valid PAS/OI queue. Two valid choices:
    #   station 1: M/M/k OI queue, class-independent  -> mu1(c) = min(n,k1)*s1
    #   station 2: infinite-server, class-dependent   -> mu2(c) = sum_j beta2(c_j)
    # For station 2 the prefix sums depend on order, so the ordered-state
    # enumeration is genuinely needed (it does not collapse to class counts).
    s1, k1 = 1.0, 2             # station 1: 2-server OI (M/M/2)
    beta2 = [1.5, 1.0]          # station 2: per-class infinite-server rates
    mu1 = lambda c: min(len(c), k1) * s1
    mu2 = lambda c: sum(beta2[j] for j in c)

    # ---- brute-force enumeration ---------------------------------------
    G = 0.0
    for m in enum_splits(K):
        S1 = oi_sum(m, e1, mu1)                              # orderings at station 1
        S2 = oi_sum([K[r] - m[r] for r in range(R)], e2, mu2)  # at station 2
        G += S1 * S2

    print('Closed cyclic PAS network: R=%d classes, population K=[%s]'
          % (R, ' '.join(str(k) for k in K)))
    print('Brute-force normalizing constant G = %.15g' % G)

    # ---- per-class throughputs from the normalizing constant -----------
    # For product-form closed networks the chain throughput satisfies
    # X_r = e_r * G(K - 1_r) / G(K). With e_r = 1 this is just the ratio of
    # normalizing constants at populations K-1_r and K.
    X = [0.0] * R
    for r in range(R):
        if K[r] > 0:
            Km = list(K)
            Km[r] -= 1
            X[r] = norm_const(Km, e1, e2, mu1, mu2) / G
    print('Per-class throughput X = [%s]' % ' '.join('%.6g' % x for x in X))

    # ---- self-test: single-class reduction to Gordon-Newell ------------
    selftest_single_class()
    return G


if __name__ == '__main__':
    pas_cyclic_normconst_bruteforce()
