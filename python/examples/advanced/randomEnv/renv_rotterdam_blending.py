"""
renv_rotterdam_blending - Exact reproduction of the fyp26 SOQN blending result

Reproduces the 24-environment exponential blending accuracy result of the fyp26
Rotterdam container-terminal study (Dhingra 2018), using the same
level-dependent QBD blocks and the cyclic resolvent blending that the ENV
state-vector analyzer uses internally.

Per-environment model (PoissonSOQNSpec): a semi-open SOQN with N tokens, a
Poisson arrival rate lambda_h (hour h of the day), and two flow-equivalent
server (FES) subnetworks in tandem (upstream S1, downstream S2), each a
load-dependent station whose throughput curve mu1(n)/mu2(n) is calibrated by
exact MVA on the underlying closed subnetwork. The LD-QBD state is (n,k): level
n = jobs upstream of S2 (S1 + external backlog), phase k = S2 occupancy.

The 24 hourly environments are visited cyclically with exponential sojourns
(mean = hour duration). Each visit applies the resolvent s*pi*(sI-Q)^{-1}
(exit = time-average, memoryless), entry vectors are chained to an L1 fixed
point, and the blend weights stages by their time fraction. Mean external wait W
and external queue length Qex are compared to the Dhingra simulation (Table 13).
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from line_solver import (ClosedClass, Delay, Exp, GlobalConstants, MVA, Network,
                         Queue, SchedStrategy, VerboseLevel)


def fes_curve(N, p, which):
    """Load-dependent FES throughput curve mu(l), l = 0..N, by exact MVA on the
    closed subnetwork at each population l (Norton flow-equivalent)."""
    mu = np.zeros(N + 1)
    for l in range(1, N + 1):
        if which == 'up':
            m = Network('S1')
            eg = Queue(m, 'EntryGates', SchedStrategy.FCFS)
            eg.setNumberOfServers(p['numEntryServers'])
            tv = Delay(m, 'TravelToStack')
            st = [Queue(m, 'Stack%d' % (i + 1), SchedStrategy.FCFS)
                  for i in range(p['numStacks'])]
            for q in st:
                q.setNumberOfServers(1)
            cls = ClosedClass(m, 'Trucks', l, eg)
            eg.setService(cls, Exp(1.0 / p['entryServiceTime']))
            tv.setService(cls, Exp(1.0 / p['travelToStackTime']))
            for q in st:
                q.setService(cls, Exp(1.0 / p['stackServiceTime']))
            ns = 2 + p['numStacks']
            R = np.zeros((ns, ns))
            R[0, 1] = 1.0                                   # EntryGates -> TravelToStack
            R[1, 2:ns] = 1.0 / p['numStacks']               # TravelToStack -> Stack_i
            R[2:ns, 0] = 1.0                                # Stack_i -> EntryGates
            P = m.initRoutingMatrix()
            P.set(cls, cls, R)
            m.link(P)
            T = MVA(m, 'exact', verbose=False).getAvg()[3]
            mu[l] = float(np.sum(T[0, :]))                  # throughput at EntryGates
        else:
            m = Network('S2')
            tv = Delay(m, 'TravelToExit')
            xg = Queue(m, 'ExitGates', SchedStrategy.FCFS)
            xg.setNumberOfServers(p['numExitServers'])
            cls = ClosedClass(m, 'Trucks', l, tv)
            tv.setService(cls, Exp(1.0 / p['travelToExitTime']))
            xg.setService(cls, Exp(1.0 / p['exitServiceTime']))
            P = m.initRoutingMatrix()
            P.set(cls, cls, np.array([[0.0, 1.0], [1.0, 0.0]]))
            m.link(P)
            T = MVA(m, 'exact', verbose=False).getAvg()[3]
            mu[l] = float(np.sum(T[1, :]))                  # throughput at ExitGates
    return mu


def soqn_generator(N, lam, mu1, mu2, tail_factor):
    """Flat generator of the Poisson SOQN LD-QBD. State (n,k), flat index n*M+k.
    Level n = 0..Mtr (S1 + backlog), phase k = 0..N (S2 occupancy)."""
    M = N + 1
    Mtr = N + tail_factor * N
    dim = (Mtr + 1) * M
    rows, cols, vals = [], [], []
    for n in range(Mtr + 1):
        at_top = (n == Mtr)
        lam_eff = 0.0 if at_top else lam   # no arrivals leave the truncated chain at top
        for k in range(N + 1):
            row = n * M + k
            m1 = mu1[min(n, N - k)]
            m2 = mu2[k]
            rows.append(row); cols.append(row); vals.append(-(lam_eff + m1 + m2))
            if k >= 1:                     # S2 completion within level: (n,k)->(n,k-1)
                rows.append(row); cols.append(n * M + (k - 1)); vals.append(m2)
            if n < Mtr:                    # arrival up: (n,k)->(n+1,k)
                rows.append(row); cols.append((n + 1) * M + k); vals.append(lam)
            if n >= 1 and k <= N - 1:      # S1 completion down: (n,k)->(n-1,k+1)
                rows.append(row); cols.append((n - 1) * M + (k + 1)); vals.append(m1)
    return sparse.csc_matrix((vals, (rows, cols)), shape=(dim, dim))


def _resolvent(pi_row, s, Q, I):
    """s * pi * (sI - Q)^{-1}, i.e. the time-average of an exponentially
    terminated sojourn started at pi. Solved on the transpose, since spsolve
    takes a column system."""
    return s * spsolve((s * I - Q).T.tocsc(), np.asarray(pi_row).ravel())


def blend_soqn(N, lam, dur_min, frac_w, mu1, mu2, tail_factor):
    """Cyclic exponential blending over the 24 hourly SOQN environments, via the
    resolvent s*pi*(sI-Q)^{-1}, to an L1 fixed point; then blend and extract
    W, Qex."""
    K = len(lam)
    M = N + 1
    Mtr = N + tail_factor * N
    dim = (Mtr + 1) * M
    Q = [soqn_generator(N, lam[h], mu1, mu2, tail_factor) for h in range(K)]

    pi0 = np.zeros(dim)
    pi0[0] = 1.0                                   # empty system
    pi_enter = [pi0.copy() for _ in range(K)]

    I = sparse.eye(dim, format='csc')
    for _ in range(200):
        prev = [v.copy() for v in pi_enter]
        for h in range(K):
            s = 1.0 / dur_min[h]
            pi_enter[(h + 1) % K] = _resolvent(pi_enter[h], s, Q[h], I)
        l1 = max(float(np.sum(np.abs(pi_enter[h] - prev[h]))) for h in range(K))
        if l1 < 1e-10:
            break

    # Blend (exp sojourn: exit = time-average), weighted by time fraction
    pi_avg = np.zeros(dim)
    for h in range(K):
        s = 1.0 / dur_min[h]
        pi_avg += frac_w[h] * _resolvent(pi_enter[h], s, Q[h], I)
    pi_avg = np.clip(pi_avg, 0.0, None)
    pi_avg /= pi_avg.sum()

    # Metrics (extractMetricsPoissonSOQN)
    QLex = QL1 = QL2 = tput = 0.0
    for n in range(Mtr + 1):
        for k in range(N + 1):
            p = pi_avg[n * M + k]
            if p == 0.0:
                continue
            QLex += max(0, n + k - N) * p
            QL1 += min(n, N - k) * p
            QL2 += k * p
            if k >= 1:
                tput += mu2[k] * p
    return (QLex + QL1 + QL2) / tput, QLex     # W (Little), Qex


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.SILENT)

    # Dhingra 24-hour schedule and simulation ground truth (Table 13)
    daily_rate_hr = np.array([6, 30, 40, 62, 76, 79, 119, 164, 152, 130, 79, 70,
                              57, 57, 113, 130, 162, 202, 148, 118, 92, 62, 36, 8],
                             dtype=float)
    daily_dur_hr = np.array([0.51, 0.84, 0.76, 0.98, 0.67, 2.33, 1.80, 0.62, 0.36,
                             1.30, 1.02, 0.98, 0.86, 1.56, 1.30, 1.57, 0.34, 0.22,
                             0.47, 1.33, 1.56, 0.93, 0.71, 1.11], dtype=float)
    lam = daily_rate_hr * 0.5 / 60.0        # arrivals per minute (buildLambdasMin)
    dur_min = daily_dur_hr * 60.0           # hour durations in minutes
    frac_w = dur_min / dur_min.sum()        # time-fraction blend weights f_i

    sim_N = np.arange(24, 35)
    sim_W = np.array([154.0661, 118.3936, 99.962, 88.2379, 80.9272, 75.0887,
                      70.745, 67.4326, 64.9719, 62.7913, 61.2534])
    sim_Qex = np.array([89.8508, 63.3970, 49.3650, 40.8783, 35.3546, 30.6543,
                        27.4111, 24.7218, 22.3622, 20.3893, 18.9252])

    tail_factor = 15                        # M_trunc = N*(1+tailFactor) (ACC_TAIL)
    params = {'numEntryServers': 6, 'numStacks': 29, 'entryServiceTime': 6.0,
              'travelToStackTime': 5.6, 'stackServiceTime': 6.0,
              'numExitServers': 6, 'travelToExitTime': 5.6, 'exitServiceTime': 6.0}

    print('Rotterdam SOQN exponential blending vs Dhingra simulation (tailFactor=%d)'
          % tail_factor)
    print('%4s | %10s %10s %7s | %10s %10s %7s'
          % ('N', 'W_blend', 'W_sim', 'err%', 'Qex_bl', 'Qex_sim', 'err%'))
    for N in (24, 28, 34):
        mu1 = fes_curve(N, params, 'up')
        mu2 = fes_curve(N, params, 'down')
        W_bl, Qex_bl = blend_soqn(N, lam, dur_min, frac_w, mu1, mu2, tail_factor)
        j = int(np.flatnonzero(sim_N == N)[0])
        print('%4d | %10.4f %10.4f %6.2f | %10.4f %10.4f %6.2f'
              % (N, W_bl, sim_W[j], 100 * abs(W_bl - sim_W[j]) / sim_W[j],
                 Qex_bl, sim_Qex[j], 100 * abs(Qex_bl - sim_Qex[j]) / sim_Qex[j]))
