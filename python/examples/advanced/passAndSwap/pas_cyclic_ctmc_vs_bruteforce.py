"""
Validate the CTMC (and LDES) solvers on a CLOSED CYCLIC network of two
pass-and-swap (PAS / order-independent) queues against the exact
product-form brute-force normalizing constant.

    Topology:   --> Queue1 --> Queue2 -->   (cyclic, R closed classes)

PAS queues are order-independent, so the network is product-form. The
brute force enumerates the whole ORDERED state space to obtain the
normalizing constant G, hence exact per-class throughputs and mean queue
lengths, which are compared against CTMC (exact) and LDES
(simulation, agreement only within Monte-Carlo noise).

Both stations must be valid order-independent (OI) queues: the TOTAL
service rate mu_i(c) must depend only on the customer multiset.
  station 1: M/M/k OI queue, class-independent  -> mu1(c) = min(n,k1)*s1
  station 2: infinite-server, class-dependent   -> mu2(c) = sum_j beta2[c_j]
For station 2 the prefix sums depend on order, so the ordered-state
enumeration is genuinely needed (it does not collapse to class counts).

Class indices passed to mu(.) and the swap graph are 0-based (Pythonic).
"""

import itertools
import numpy as np
from line_solver import *

# ---- shared model parameters ------------------------------------------
K = [2, 2]                 # per-class closed populations
s1, k1 = 1.0, 2            # station 1: 2-server OI (M/M/2)
beta2 = [1.5, 1.0]         # station 2: per-class infinite-server rates
R = len(K)


def mu1(c):
    c = np.asarray(c)
    return float(min(len(c), k1) * s1)


def mu2(c):
    c = np.asarray(c, dtype=int)
    return float(sum(beta2[ci] for ci in c))


# =======================================================================
# (1) brute-force product-form reference
# =======================================================================
def oi_dfs(counts, prefix, w, e, mu):
    if all(x == 0 for x in counts):
        return w
    s = 0.0
    for r in range(len(counts)):
        if counts[r] > 0:
            counts[r] -= 1
            prefix.append(r)
            s += oi_dfs(counts, prefix, w * e[r] / mu(prefix), e, mu)
            prefix.pop()
            counts[r] += 1
    return s


def oi_sum(counts, e, mu):
    # sum of the OI balance function over every distinct ordered arrangement
    return oi_dfs(list(counts), [], 1.0, e, mu)


def enum_splits(Kvec):
    return list(itertools.product(*[range(k + 1) for k in Kvec]))


def norm_const(Kvec, e1, e2, f1, f2):
    G = 0.0
    for m in enum_splits(Kvec):
        comp = [Kvec[r] - m[r] for r in range(len(Kvec))]
        G += oi_sum(m, e1, f1) * oi_sum(comp, e2, f2)
    return G


def bruteforce(Kvec, f1, f2):
    Rn = len(Kvec)
    e1 = [1.0] * Rn
    e2 = [1.0] * Rn
    G = norm_const(Kvec, e1, e2, f1, f2)
    X = np.zeros(Rn)
    for r in range(Rn):                       # X_r = e_r * G(K-1_r)/G(K), e_r=1
        Km = list(Kvec)
        Km[r] -= 1
        X[r] = norm_const(Km, e1, e2, f1, f2) / G
    Q = np.zeros((2, Rn))                      # Q(1,r) = (1/G) sum_m m_r S1 S2
    for m in enum_splits(Kvec):
        comp = [Kvec[r] - m[r] for r in range(Rn)]
        w = oi_sum(m, e1, f1) * oi_sum(comp, e2, f2)
        for r in range(Rn):
            Q[0, r] += m[r] * w
    Q[0, :] /= G
    Q[1, :] = np.array(Kvec) - Q[0, :]
    return G, X, Q


Gbf, Xbf, Qbf = bruteforce(K, mu1, mu2)
print("Brute-force normalizing constant G = %.12g\n" % Gbf)

# =======================================================================
# (2) LINE closed PAS network
# =======================================================================
def build_model():
    model = Network('PAScyclic')
    q1 = Queue(model, 'PASQueue1', SchedStrategy.PAS)
    q2 = Queue(model, 'PASQueue2', SchedStrategy.PAS)
    classes = [ClosedClass(model, 'Class%d' % (r + 1), K[r], q1) for r in range(R)]
    q1.setService(mu1)
    q2.setService(mu2)
    q1.setSwapGraph(np.zeros((R, R)))
    q1.setNumberOfServers(k1)
    q1.setCap(sum(K))
    q2.setSwapGraph(np.zeros((R, R)))
    q2.setNumberOfServers(sum(K))
    q2.setCap(sum(K))
    model.addLink(q1, q2)
    model.addLink(q2, q1)
    for r in range(R):
        q1.setProbRouting(classes[r], q2, 1.0)
        q2.setProbRouting(classes[r], q1, 1.0)
    return model


def extract(tab):
    # pull X (chain throughput, from q1) and Q[station,class] from an AvgTable
    names = ['PASQueue1', 'PASQueue2']
    col = 'Node' if 'Node' in tab.columns else 'Station'
    X = np.zeros(R)
    Q = np.zeros((2, R))
    for i, nm in enumerate(names):
        for r in range(R):
            row = tab[(tab[col] == nm) & (tab['JobClass'] == 'Class%d' % (r + 1))]
            Q[i, r] = float(row['QLen'].iloc[0])
            if i == 0:
                X[r] = float(row['Tput'].iloc[0])
    return X, Q


# =======================================================================
# (3) CTMC (exact) comparison
# =======================================================================
ctmc = CTMC(build_model(), cutoff=sum(K))
Tab_c = ctmc.getAvgTable()
print(Tab_c)
Xc, Qc = extract(Tab_c)               # one visit/station -> chain throughput

print("\n--- CTMC vs brute-force ---")
print("%-8s %-8s %14s %14s %10s" % ('metric', 'class', 'brute-force', 'CTMC', 'abs.err'))
for r in range(R):
    print("Tput     Class%-2d %14.10f %14.10f %10.2e" % (r + 1, Xbf[r], Xc[r], abs(Xbf[r] - Xc[r])))
for i in range(2):
    for r in range(R):
        print("QLen q%-1d   Class%-2d %14.10f %14.10f %10.2e" %
              (i + 1, r + 1, Qbf[i, r], Qc[i, r], abs(Qbf[i, r] - Qc[i, r])))

tol = 1e-9
errX = np.max(np.abs(Xbf - Xc))
errQ = np.max(np.abs(Qbf - Qc))
print("\nCTMC:  max|dX| = %.3e   max|dQ| = %.3e" % (errX, errQ))
assert errX <= tol and errQ <= tol, "CTMC vs brute-force mismatch"
print("PASS: CTMC matches brute-force product form within %.1e." % tol)

# =======================================================================
# (4) LDES (simulation) comparison
# =======================================================================
ldes = LDES(build_model(), samples=200000, seed=23000)
Xl, Ql = extract(ldes.getAvgTable())

print("\n--- LDES vs brute-force (simulation, 2e5 samples) ---")
for r in range(R):
    print("Tput     Class%-2d %14.6f %14.6f %8.2f%%" %
          (r + 1, Xbf[r], Xl[r], 100 * abs(Xbf[r] - Xl[r]) / Xbf[r]))
for i in range(2):
    for r in range(R):
        print("QLen q%-1d   Class%-2d %14.6f %14.6f %8.2f%%" %
              (i + 1, r + 1, Qbf[i, r], Ql[i, r], 100 * abs(Qbf[i, r] - Ql[i, r]) / Qbf[i, r]))

tol_sim = 0.02
relX = np.max(np.abs(Xbf - Xl) / Xbf)
relQ = np.max(np.abs(Qbf - Ql) / Qbf)
print("\nLDES:  max rel|dX| = %.2f%%   max rel|dQ| = %.2f%%" % (100 * relX, 100 * relQ))
assert relX <= tol_sim and relQ <= tol_sim, "LDES vs brute-force exceeds 2%"
print("PASS: LDES matches brute-force within %.0f%% (simulation noise)." % (100 * tol_sim))
