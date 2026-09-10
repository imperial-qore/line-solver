/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAPQN_MAPQN_QR_BOUNDS_BAS_H
#define LINE_API_MAPQN_MAPQN_QR_BOUNDS_BAS_H

/**
 * Quadratic-reduction bound on the utilization of one queue of a closed MAP
 * queueing network with a FINITE-CAPACITY station under blocking-after-service.
 *
 * Templated port of matlab/lib/qrf/qrf_bas.m (ground truth), whose own origin is
 * the AMPL model qrboundsbas_skel.mod.
 *
 * WHAT BAS MEANS HERE. Station f has capacity F(f). A job completing at station
 * j and routed to a full f cannot move, so j is BLOCKED: it holds the completed
 * job and serves nothing until f frees a slot. The set of currently blocked
 * stations, and the order in which they blocked, is the blocking configuration
 * m in 0..MR-1. Every probability variable carries it, which is what separates
 * this model from mapqn_bnd_qr_ld: the state is (pairwise occupancy, blocking
 * configuration), not occupancy alone.
 *
 * THE BLOCKING TABLES ARE SUPPLIED, NOT DERIVED. BB, MM, ZZ, ZM and MM1 are
 * caller data here exactly as they are in the reference, where
 * example_bas_small.m writes them out literally. Deriving them from (M, f)
 * would be a different routine and is deliberately not attempted.
 *
 * FAMILY-BY-FAMILY TRAPS, all verified against the reference and recorded in
 * _kb/03-api-layer.md. Read them before editing any family below:
 *  - THM30 and THM3 do NOT share a loop nest. THM30's RHS sums over hj with a
 *    coefficient independent of hj (a pure multiplicity); THM3's sums over hi
 *    with a coefficient that depends on it. Copying one into the other scales
 *    the row by K and still solves.
 *  - THM3f is pinned to configuration 0, not looped over m: below capacity
 *    there is no blocking, so only the unblocked configuration contributes.
 *  - THM3L couples TWO configurations through mp = MM1(m,j), and its RHS
 *    variable is the diagonal in mp, not in m.
 *  - THM4 is accumulated as (sum nt p2 - N sum p2) and then NEGATED, giving
 *    N P(j at (nj,kj), i nonempty) <= sum_t E[n_t ...].
 *  - SYMMETRY skips pairs both of whose members are already pinned to zero, so
 *    it must run AFTER the ZERO pass. The emission order below is load bearing.
 *
 * UPPER BOUNDS ARE INFINITE, unlike every other bound in this family. The
 * reference initializes ub = inf and only the ZERO families pin anything; the
 * variables are bounded above through ONE instead. lp::LpModel's default is
 * exactly lb = 0 with a free upper bound, so nothing is set here beyond the
 * ZERO pins.
 *
 * COST. MR * B^2 + sum_i K(i) columns with B = (N+1) sum_i K(i), so the model
 * is MR times the load-dependent one. The tableau here is dense, which confines
 * the port to small instances; the reference's own paper instance is ~6e4
 * columns and needs a sparse revised simplex.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mapqn/mapqn_params.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/lp_highs.h"
#include "line/util/simplex.h"

namespace line {
namespace mapqn {

/**
 * Parameters of the BAS bound, mirroring the reference's `params`.
 *
 * Queues, phases and blocking configurations are 0-based here; the reference is
 * 1-based in all three, so `f`, the entries of `MM` and the entries of `MM1`
 * all drop by one relative to a MATLAB script. `MM1` keeps the reference's
 * "absent" marker as a negative entry rather than 0.
 */
template <class T>
struct QrBasParams {
    int M = 0;                    ///< number of queues
    int N = 0;                    ///< total population
    int f = 0;                    ///< index of the finite-capacity queue, 0-based
    std::vector<int> F;           ///< (M) capacity of each queue
    std::vector<int> K;           ///< (M) number of phases of each queue
    std::vector<Matrix<T>> mu;    ///< mu[i] is K(i) x K(i), completion rates
    std::vector<Matrix<T>> v;     ///< v[i] is K(i) x K(i), background rates
    Matrix<T> r;                  ///< (M x M) routing probabilities
    int MR = 0;                   ///< number of blocking configurations
    std::vector<std::vector<int> > BB;   ///< (MR x M) 1 if queue i is blocked in m
    std::vector<std::vector<int> > MM;   ///< (MR x 2) blocking order, 0-based queue indices
    std::vector<int> ZZ;          ///< (MR) number of blocked queues in m
    int ZM = 0;                   ///< maximum blocking depth
    std::vector<std::vector<int> > MM1;  ///< (MR x M) extended order; negative means absent

    void validate() const {
        if (M <= 0) throw InputError("qrf_bas: M must be positive");
        if (N < 1) throw InputError("qrf_bas: N must be at least 1");
        if (f < 0 || f >= M) throw InputError("qrf_bas: f out of range");
        if (static_cast<int>(F.size()) != M) throw InputError("qrf_bas: F has the wrong length");
        if (static_cast<int>(K.size()) != M) throw InputError("qrf_bas: K has the wrong length");
        if (static_cast<int>(mu.size()) != M || static_cast<int>(v.size()) != M)
            throw InputError("qrf_bas: mu and v must have one entry per queue");
        for (int i = 0; i < M; ++i) {
            if (K[i] <= 0) throw InputError("qrf_bas: every queue needs at least one phase");
            if (F[i] < 0 || F[i] > N) throw InputError("qrf_bas: F(i) must lie in 0..N");
            const std::size_t k = static_cast<std::size_t>(K[i]);
            if (mu[i].rows() != k || mu[i].cols() != k)
                throw InputError("qrf_bas: mu{i} must be K(i) x K(i)");
            if (v[i].rows() != k || v[i].cols() != k)
                throw InputError("qrf_bas: v{i} must be K(i) x K(i)");
        }
        if (r.rows() != static_cast<std::size_t>(M) || r.cols() != static_cast<std::size_t>(M))
            throw InputError("qrf_bas: r must be M x M");
        if (MR < 1) throw InputError("qrf_bas: MR must be at least 1");
        if (static_cast<int>(BB.size()) != MR || static_cast<int>(MM.size()) != MR ||
            static_cast<int>(ZZ.size()) != MR || static_cast<int>(MM1.size()) != MR)
            throw InputError("qrf_bas: the blocking tables must have MR rows");
        for (int m = 0; m < MR; ++m) {
            if (static_cast<int>(BB[m].size()) != M || static_cast<int>(MM1[m].size()) != M)
                throw InputError("qrf_bas: BB and MM1 must have M columns");
            if (MM[m].size() < 2) throw InputError("qrf_bas: MM must have two columns");
        }
        // ZM IS max(ZZ), and a larger one couples THM3I to a depth that has no
        // configuration, emptying the polytope instead of failing cleanly. Zero
        // is legal and means no configuration blocks anything, which is what
        // `qrf_bas.m` runs when its `for z = 0:(ZM-1)` is empty.
        int zmax = 0;
        for (std::size_t m = 0; m < ZZ.size(); ++m) zmax = std::max(zmax, ZZ[m]);
        if (ZM != zmax)
            throw InputError("qrf_bas: ZM must equal max(ZZ) (got ZM = " + std::to_string(ZM) +
                             ", max(ZZ) = " + std::to_string(zmax) + ")");
    }
};

/** Result of a BAS bound solve. */
template <class T>
struct QrBasResult {
    bool ok = false;      ///< the LP reached an optimal vertex
    std::string status;   ///< textual LP status
    T objective = T();    ///< the bound on the utilization of the target queue
    std::vector<T> U;     ///< (M) utilization of each queue at the optimal vertex
    /**
     * (M) P(n_i >= 1) over EVERY configuration, blocked ones included. This is
     * what `U` used to hold; at a BAS station it counts a blocked server as
     * busy, so it is occupancy and NOT utilization. Kept because it is the
     * quantity the QRF papers report.
     */
    std::vector<T> occupancy;
    Matrix<T> e;          ///< M x max(K) effective per-phase utilizations
    std::vector<T> x;     ///< full solution vector, indexed by QrBasIndex
    std::size_t num_vars = 0;
    std::size_t num_rows = 0;
    std::size_t iterations = 0;
};

/**
 * Variable layout: p2(j,nj,kj,i,ni,hi,m) then e(i,ki).
 *
 * The pairwise half-index is the one the rest of the mapqn family uses, so the
 * blocking configuration is the FASTEST varying subscript here. The reference
 * nests it differently (j, nj, kj, i, m, ni, hi); order does not affect the
 * polytope, only how the two assemblies diff.
 */
struct QrBasIndex {
    int M = 0, N = 0, MR = 0;
    std::vector<int> K;
    std::vector<int> cumK;
    std::size_t sumK = 0, block = 0, off_e = 0, total = 0;

    QrBasIndex() {}
    QrBasIndex(int m, int n, const std::vector<int>& k, int mr) : M(m), N(n), MR(mr), K(k) {
        cumK.assign(static_cast<std::size_t>(M) + 1, 0);
        for (int i = 0; i < M; ++i) cumK[i + 1] = cumK[i] + K[i];
        sumK = static_cast<std::size_t>(cumK[M]);
        block = static_cast<std::size_t>(N + 1) * sumK;
        off_e = block * block * static_cast<std::size_t>(MR);
        total = off_e + sumK;
    }

    std::size_t half(int i, int ni, int h) const {
        return static_cast<std::size_t>(N + 1) * static_cast<std::size_t>(cumK[i]) +
               static_cast<std::size_t>(ni) * static_cast<std::size_t>(K[i]) +
               static_cast<std::size_t>(h);
    }
    std::size_t p2(int j, int nj, int kj, int i, int ni, int hi, int m) const {
        return (half(j, nj, kj) * block + half(i, ni, hi)) * static_cast<std::size_t>(MR) +
               static_cast<std::size_t>(m);
    }
    std::size_t e(int i, int ki) const {
        return off_e + static_cast<std::size_t>(cumK[i]) + static_cast<std::size_t>(ki);
    }
    std::size_t num_vars() const { return total; }
};

namespace detail {

/** q(i,j,k,h), the same non-load-dependent rate the rest of the family uses. */
template <class T>
T bas_rate(const QrBasParams<T>& p, int i, int j, int k, int h) {
    const std::size_t ki = static_cast<std::size_t>(k), hi = static_cast<std::size_t>(h);
    const std::size_t ii = static_cast<std::size_t>(i), ji = static_cast<std::size_t>(j);
    if (j != i) return T(p.r(ii, ji) * p.mu[i](ki, hi));
    return T(p.v[i](ki, hi) + p.r(ii, ii) * p.mu[i](ki, hi));
}

/**
 * ZERO1..ZERO8: the states that carry no mass, pinned as ub = 0.
 *
 * Returns the indicator, because SYMMETRY skips pairs both of whose members are
 * pinned and therefore depends on this pass having run.
 */
template <class T>
std::vector<char> bas_zero(const QrBasParams<T>& p, const QrBasIndex& x, lp::LpModel<T>& m) {
    std::vector<char> zero(x.num_vars(), 0);
    for (int j = 0; j < p.M; ++j) {
        for (int nj = 0; nj <= p.N; ++nj) {
            for (int kj = 0; kj < p.K[j]; ++kj) {
                for (int i = 0; i < p.M; ++i) {
                    for (int ni = 0; ni <= p.N; ++ni) {
                        for (int hi = 0; hi < p.K[i]; ++hi) {
                            for (int mm = 0; mm < p.MR; ++mm) {
                                bool z = false;
                                if (i == j && nj == ni && hi != kj) z = true;          // ZERO1
                                if (i == j && nj != ni) z = true;                      // ZERO2
                                if (i != j && nj + ni > p.N) z = true;                 // ZERO3
                                if (nj > p.F[j]) z = true;                             // ZERO6
                                if (mm >= 1 && p.BB[mm][j] == 1 && nj == 0) z = true;  // ZERO5
                                if (mm >= 1 && p.BB[mm][j] == 1 && i != j && i != p.f &&
                                    ni + nj + p.F[p.f] > p.N)
                                    z = true;  // ZERO7
                                if (j == p.f && nj >= 1 && nj <= p.F[p.f] - 1 && mm >= 1)
                                    z = true;  // ZERO8
                                if (z) {
                                    const std::size_t idx = x.p2(j, nj, kj, i, ni, hi, mm);
                                    m.fix(idx, T());
                                    zero[idx] = 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // ZERO4: for m >= 1 and j != f, the finite queue below capacity carries no
    // mass in a blocking configuration.
    for (int j = 0; j < p.M; ++j) {
        if (j == p.f) continue;
        for (int nj = 0; nj <= p.N; ++nj) {
            for (int kj = 0; kj < p.K[j]; ++kj) {
                for (int mm = 1; mm < p.MR; ++mm) {
                    for (int nf = 0; nf <= p.F[p.f] - 1; ++nf) {
                        for (int hf = 0; hf < p.K[p.f]; ++hf) {
                            const std::size_t idx = x.p2(j, nj, kj, p.f, nf, hf, mm);
                            m.fix(idx, T());
                            zero[idx] = 1;
                        }
                    }
                }
            }
        }
    }
    return zero;
}

/** ONE: the diagonal marginal of each queue normalizes to one. */
template <class T>
void bas_one(const QrBasParams<T>& p, const QrBasIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int nj = 0; nj <= p.N; ++nj)
            for (int kj = 0; kj < p.K[j]; ++kj)
                for (int mm = 0; mm < p.MR; ++mm) m.row_add_int(x.p2(j, nj, kj, j, nj, kj, mm), 1);
        m.emit_eq_int(1);
    }
}

/** SYMMETRY: p2 is symmetric under swapping its halves, within a configuration. */
template <class T>
void bas_symmetry(const QrBasParams<T>& p, const QrBasIndex& x, lp::LpModel<T>& m,
                  const std::vector<char>& zero) {
    for (int j = 0; j < p.M; ++j) {
        const int njmax = p.N < p.F[j] ? p.N : p.F[j];
        for (int nj = 0; nj <= njmax; ++nj) {
            for (int kj = 0; kj < p.K[j]; ++kj) {
                for (int i = j + 1; i < p.M; ++i) {
                    const int nimax = p.N < p.F[i] ? p.N : p.F[i];
                    for (int ni = 0; ni <= nimax; ++ni) {
                        if (nj + ni > p.N) continue;
                        for (int hi = 0; hi < p.K[i]; ++hi) {
                            for (int mm = 0; mm < p.MR; ++mm) {
                                const std::size_t a = x.p2(j, nj, kj, i, ni, hi, mm);
                                const std::size_t b = x.p2(i, ni, hi, j, nj, kj, mm);
                                if (zero[a] && zero[b]) continue;
                                if (a == b) continue;
                                m.row_add_int(a, 1);
                                m.row_add_int(b, -1);
                                m.emit_eq_int(0);
                            }
                        }
                    }
                }
            }
        }
    }
}

/** MARGINALS: the pairwise law agrees with its own marginal at every level. */
template <class T>
void bas_marginals(const QrBasParams<T>& p, const QrBasIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int kj = 0; kj < p.K[j]; ++kj) {
            const int njmax = p.N < p.F[j] ? p.N : p.F[j];
            for (int nj = 0; nj <= njmax; ++nj) {
                for (int i = 0; i < p.M; ++i) {
                    if (i == j) continue;
                    for (int mm = 0; mm < p.MR; ++mm) {
                        m.row_add_int(x.p2(j, nj, kj, j, nj, kj, mm), 1);
                        const int cap = (p.N - nj) < p.F[i] ? (p.N - nj) : p.F[i];
                        for (int ni = 0; ni <= cap; ++ni)
                            for (int hi = 0; hi < p.K[i]; ++hi)
                                m.row_add_int(x.p2(j, nj, kj, i, ni, hi, mm), -1);
                        m.emit_eq_int(0);
                    }
                }
            }
        }
    }
}

/**
 * UEFF: e(i,ki) is the probability that queue i is busy in phase ki AND NOT
 * blocked. It is the effective utilization, which is what the objective and
 * THM1 read, and it is why a blocked station contributes nothing to throughput.
 */
template <class T>
void bas_ueff(const QrBasParams<T>& p, const QrBasIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int ki = 0; ki < p.K[i]; ++ki) {
            m.row_add_int(x.e(i, ki), -1);
            for (int j = 0; j < p.M; ++j) {
                const int njmax = p.N < p.F[j] ? p.N : p.F[j];
                for (int nj = 0; nj <= njmax; ++nj) {
                    for (int kj = 0; kj < p.K[j]; ++kj) {
                        for (int mm = 0; mm < p.MR; ++mm) {
                            if (p.BB[mm][i] != 0) continue;
                            const int nimax = p.N < p.F[i] ? p.N : p.F[i];
                            for (int ni = 1; ni <= nimax; ++ni)
                                m.row_add_int(x.p2(j, nj, kj, i, ni, ki, mm), 1);
                        }
                    }
                }
            }
            m.emit_eq_int(0);
        }
    }
}

/** THM1: phase balance on the effective utilizations. */
template <class T>
void bas_thm1(const QrBasParams<T>& p, const QrBasIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int ki = 0; ki < p.K[i]; ++ki) {
            for (int j = 0; j < p.M; ++j)
                for (int hi = 0; hi < p.K[i]; ++hi)
                    if (j != i || hi != ki) m.row_add(x.e(i, ki), bas_rate(p, i, j, ki, hi));
            for (int j = 0; j < p.M; ++j)
                for (int hi = 0; hi < p.K[i]; ++hi)
                    if (j != i || hi != ki)
                        m.row_add(x.e(i, hi), T(-bas_rate(p, i, j, hi, ki)));
            m.emit_eq_int(0);
        }
    }
}

/** THM2: the queue-length theorem conditioned on (j,nj,kj,m). */
template <class T>
void bas_thm2(const QrBasParams<T>& p, const QrBasIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int kj = 0; kj < p.K[j]; ++kj) {
            for (int nj = 0; nj <= p.F[j]; ++nj) {
                for (int mm = 0; mm < p.MR; ++mm) {
                    m.row_add_int(x.p2(j, nj, kj, j, nj, kj, mm), -p.N);
                    for (int i = 0; i < p.M; ++i)
                        for (int ni = 1; ni <= p.F[i]; ++ni)
                            for (int ki = 0; ki < p.K[i]; ++ki)
                                m.row_add_int(x.p2(j, nj, kj, i, ni, ki, mm), ni);
                    m.emit_eq_int(0);
                }
            }
        }
    }
}

/** COR1: the second moment of the population. */
template <class T>
void bas_cor1(const QrBasParams<T>& p, const QrBasIndex& x, lp::LpModel<T>& m) {
    for (int mm = 0; mm < p.MR; ++mm)
        for (int i = 0; i < p.M; ++i)
            for (int j = 0; j < p.M; ++j)
                for (int nj = 1; nj <= p.F[j]; ++nj)
                    for (int ni = 1; ni <= p.F[i]; ++ni)
                        for (int ki = 0; ki < p.K[i]; ++ki)
                            for (int kj = 0; kj < p.K[j]; ++kj)
                                m.row_add_int(x.p2(j, nj, kj, i, ni, ki, mm), ni * nj);
    m.emit_eq_int(p.N * p.N);
}

/**
 * THM30: level-crossing balance at an empty station i != f, per arrival phase.
 *
 * NOTE the RHS multiplicity: the first two RHS blocks sum over the OTHER
 * queue's phase (hj) while the coefficient q(i,j,ki,ui) does not depend on it.
 * That is deliberate and differs from THM3 below, which sums over i's own phase
 * with a coefficient that does depend on it.
 */
template <class T>
void bas_thm30(const QrBasParams<T>& p, const QrBasIndex& x, lp::LpModel<T>& m) {
    const int f = p.f;
    for (int i = 0; i < p.M; ++i) {
        if (i == f) continue;
        for (int ui = 0; ui < p.K[i]; ++ui) {
            for (int j = 0; j < p.M; ++j) {
                if (j == i || j == f) continue;
                for (int nj = 1; nj <= p.F[j]; ++nj)
                    for (int kj = 0; kj < p.K[j]; ++kj)
                        for (int hj = 0; hj < p.K[j]; ++hj)
                            for (int mm = 0; mm < p.MR; ++mm)
                                if (p.BB[mm][j] == 0)
                                    m.row_add(x.p2(j, nj, kj, i, 0, ui, mm),
                                              bas_rate(p, j, i, kj, hj));
            }
            for (int nj = 1; nj <= p.F[f]; ++nj)
                for (int kj = 0; kj < p.K[f]; ++kj)
                    for (int hj = 0; hj < p.K[f]; ++hj)
                        for (int mm = 0; mm < p.MR; ++mm)
                            if (p.MM[mm][0] != i)
                                m.row_add(x.p2(f, nj, kj, i, 0, ui, mm),
                                          bas_rate(p, f, i, kj, hj));

            for (int j = 0; j < p.M; ++j) {
                if (j == i || j == f) continue;
                for (int nj = 0; nj <= p.F[j]; ++nj)
                    for (int ki = 0; ki < p.K[i]; ++ki)
                        for (int hj = 0; hj < p.K[j]; ++hj)
                            for (int mm = 0; mm < p.MR; ++mm)
                                if (p.BB[mm][i] == 0)
                                    m.row_add(x.p2(j, nj, hj, i, 1, ki, mm),
                                              T(-bas_rate(p, i, j, ki, ui)));
            }
            for (int nj = 0; nj <= p.F[f] - 1; ++nj)
                for (int ki = 0; ki < p.K[i]; ++ki)
                    for (int hj = 0; hj < p.K[f]; ++hj)
                        for (int mm = 0; mm < p.MR; ++mm)
                            if (p.BB[mm][i] == 0)
                                m.row_add(x.p2(f, nj, hj, i, 1, ki, mm),
                                          T(-bas_rate(p, i, f, ki, ui)));
            for (int mm = 0; mm < p.MR; ++mm) {
                if (!(p.BB[mm][i] == 1 && p.MM[mm][0] == i)) continue;
                for (int kf = 0; kf < p.K[f]; ++kf)
                    for (int pf = 0; pf < p.K[f]; ++pf)
                        for (int w = 0; w < p.M; ++w)
                            if (w != f && w != i)
                                m.row_add(x.p2(f, p.F[f], kf, i, 1, ui, mm),
                                          T(-bas_rate(p, f, w, kf, pf)));
            }
            m.emit_eq_int(0);
        }
    }
}

/**
 * THM3: level-crossing balance between ni and ni+1 at a station i != f.
 *
 * NOTE the RHS multiplicity differs from THM30: here the sum runs over i's own
 * phase hi and the coefficient q(i,j,ki,hi) depends on it.
 */
template <class T>
void bas_thm3(const QrBasParams<T>& p, const QrBasIndex& x, lp::LpModel<T>& m) {
    const int f = p.f;
    for (int i = 0; i < p.M; ++i) {
        if (i == f) continue;
        for (int ni = 1; ni <= p.F[i] - 1; ++ni) {
            for (int j = 0; j < p.M; ++j) {
                if (j == i || j == f) continue;
                for (int nj = 1; nj <= p.F[j]; ++nj)
                    for (int kj = 0; kj < p.K[j]; ++kj)
                        for (int hj = 0; hj < p.K[j]; ++hj)
                            for (int ui = 0; ui < p.K[i]; ++ui)
                                for (int mm = 0; mm < p.MR; ++mm)
                                    if (p.BB[mm][j] == 0)
                                        m.row_add(x.p2(j, nj, kj, i, ni, ui, mm),
                                                  bas_rate(p, j, i, kj, hj));
            }
            for (int nj = 1; nj <= p.F[f]; ++nj)
                for (int kj = 0; kj < p.K[f]; ++kj)
                    for (int hj = 0; hj < p.K[f]; ++hj)
                        for (int ui = 0; ui < p.K[i]; ++ui)
                            for (int mm = 0; mm < p.MR; ++mm)
                                if (p.MM[mm][0] != i)
                                    m.row_add(x.p2(f, nj, kj, i, ni, ui, mm),
                                              bas_rate(p, f, i, kj, hj));

            for (int j = 0; j < p.M; ++j) {
                if (j == i || j == f) continue;
                for (int nj = 0; nj <= p.F[j]; ++nj)
                    for (int ki = 0; ki < p.K[i]; ++ki)
                        for (int hi = 0; hi < p.K[i]; ++hi)
                            for (int uj = 0; uj < p.K[j]; ++uj)
                                for (int mm = 0; mm < p.MR; ++mm)
                                    if (p.BB[mm][i] == 0)
                                        m.row_add(x.p2(j, nj, uj, i, ni + 1, ki, mm),
                                                  T(-bas_rate(p, i, j, ki, hi)));
            }
            for (int nj = 0; nj <= p.F[f] - 1; ++nj)
                for (int ki = 0; ki < p.K[i]; ++ki)
                    for (int hi = 0; hi < p.K[i]; ++hi)
                        for (int uj = 0; uj < p.K[f]; ++uj)
                            for (int mm = 0; mm < p.MR; ++mm)
                                if (p.BB[mm][i] == 0)
                                    m.row_add(x.p2(f, nj, uj, i, ni + 1, ki, mm),
                                              T(-bas_rate(p, i, f, ki, hi)));
            for (int mm = 0; mm < p.MR; ++mm) {
                if (!(p.BB[mm][i] == 1 && p.MM[mm][0] == i)) continue;
                for (int ki = 0; ki < p.K[i]; ++ki)
                    for (int kf = 0; kf < p.K[f]; ++kf)
                        for (int pf = 0; pf < p.K[f]; ++pf)
                            for (int w = 0; w < p.M; ++w)
                                if (w != f && w != i)
                                    m.row_add(x.p2(f, p.F[f], kf, i, ni + 1, ki, mm),
                                              T(-bas_rate(p, f, w, kf, pf)));
            }
            m.emit_eq_int(0);
        }
    }
}

/**
 * THM3f: level-crossing balance at the finite queue itself.
 *
 * The RHS is pinned to configuration 0: with f below capacity nothing is
 * blocked, so only the unblocked configuration can supply the departure.
 */
template <class T>
void bas_thm3f(const QrBasParams<T>& p, const QrBasIndex& x, lp::LpModel<T>& m) {
    const int f = p.f;
    for (int ni = 0; ni <= p.F[f] - 1; ++ni) {
        for (int j = 0; j < p.M; ++j) {
            if (j == f) continue;
            for (int nj = 1; nj <= p.F[j]; ++nj)
                for (int kj = 0; kj < p.K[j]; ++kj)
                    for (int hj = 0; hj < p.K[j]; ++hj)
                        for (int uf = 0; uf < p.K[f]; ++uf)
                            for (int mm = 0; mm < p.MR; ++mm)
                                if (p.BB[mm][j] == 0)
                                    m.row_add(x.p2(j, nj, kj, f, ni, uf, mm),
                                              bas_rate(p, j, f, kj, hj));
        }
        for (int j = 0; j < p.M; ++j) {
            if (j == f) continue;
            for (int nj = 0; nj <= p.F[j]; ++nj)
                for (int kf = 0; kf < p.K[f]; ++kf)
                    for (int hf = 0; hf < p.K[f]; ++hf)
                        for (int uj = 0; uj < p.K[j]; ++uj)
                            m.row_add(x.p2(j, nj, uj, f, ni + 1, kf, 0),
                                      T(-bas_rate(p, f, j, kf, hf)));
        }
        m.emit_eq_int(0);
    }
}

/** THM3I: balance across blocking depth z, at the finite queue's capacity. */
template <class T>
void bas_thm3i(const QrBasParams<T>& p, const QrBasIndex& x, lp::LpModel<T>& m) {
    const int f = p.f;
    for (int z = 0; z <= p.ZM - 1; ++z) {
        for (int j = 0; j < p.M; ++j) {
            if (j == f) continue;
            for (int nj = 1; nj <= p.F[j]; ++nj)
                for (int kj = 0; kj < p.K[j]; ++kj)
                    for (int hj = 0; hj < p.K[j]; ++hj)
                        for (int uf = 0; uf < p.K[f]; ++uf)
                            for (int mm = 0; mm < p.MR; ++mm)
                                if (p.BB[mm][j] == 0 && p.ZZ[mm] == z)
                                    m.row_add(x.p2(j, nj, kj, f, p.F[f], uf, mm),
                                              bas_rate(p, j, f, kj, hj));
        }
        for (int j = 0; j < p.M; ++j) {
            if (j == f) continue;
            for (int nj = 0; nj <= p.F[j]; ++nj)
                for (int kf = 0; kf < p.K[f]; ++kf)
                    for (int hf = 0; hf < p.K[f]; ++hf)
                        for (int uj = 0; uj < p.K[j]; ++uj)
                            for (int mm = 0; mm < p.MR; ++mm)
                                if (p.ZZ[mm] == z + 1)
                                    m.row_add(x.p2(j, nj, uj, f, p.F[f], kf, mm),
                                              T(-bas_rate(p, f, j, kf, hf)));
        }
        m.emit_eq_int(0);
    }
}

/**
 * THM3L: the maximum-blocking-depth closure.
 *
 * The only family that couples two blocking configurations: mp = MM1(m,j) is a
 * second configuration index and the RHS variable is the DIAGONAL in mp, not
 * in m. Fires only for the m at depth ZM - 1.
 */
template <class T>
void bas_thm3l(const QrBasParams<T>& p, const QrBasIndex& x, lp::LpModel<T>& m) {
    const int f = p.f;
    for (int mm = 0; mm < p.MR; ++mm) {
        if (p.ZZ[mm] != p.ZM - 1) continue;
        for (int j = 0; j < p.M; ++j) {
            if (j == f || p.BB[mm][j] != 0 || p.MM1[mm][j] < 0) continue;
            for (int nj = 1; nj <= p.F[j]; ++nj)
                for (int kj = 0; kj < p.K[j]; ++kj)
                    for (int hj = 0; hj < p.K[j]; ++hj)
                        for (int uf = 0; uf < p.K[f]; ++uf)
                            m.row_add(x.p2(j, nj, kj, f, p.F[f], uf, mm),
                                      bas_rate(p, j, f, kj, hj));
        }
        for (int j = 0; j < p.M; ++j) {
            if (j == f || p.BB[mm][j] != 0 || p.MM1[mm][j] < 0) continue;
            const int mp = p.MM1[mm][j];
            for (int kf = 0; kf < p.K[f]; ++kf)
                for (int uf = 0; uf < p.K[f]; ++uf)
                    for (int w = 0; w < p.M; ++w)
                        if (w != f)
                            m.row_add(x.p2(f, p.F[f], kf, f, p.F[f], kf, mp),
                                      T(-bas_rate(p, f, w, kf, uf)));
        }
        m.emit_eq_int(0);
    }
}

/**
 * THM4: the QMIN inequality.
 *
 * The reference accumulates (sum nt p2 - N sum p2) and emits the NEGATED row,
 * so the constraint is N P(j at (nj,kj), i nonempty) <= sum_t E[n_t ...].
 * Emitted directly in that sense here.
 */
template <class T>
void bas_thm4(const QrBasParams<T>& p, const QrBasIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int kj = 0; kj < p.K[j]; ++kj) {
            for (int i = 0; i < p.M; ++i) {
                for (int mm = 0; mm < p.MR; ++mm) {
                    for (int t = 0; t < p.M; ++t)
                        for (int ht = 0; ht < p.K[t]; ++ht)
                            for (int nj = 0; nj <= p.F[j]; ++nj)
                                for (int nt = 1; nt <= p.F[t]; ++nt)
                                    m.row_add_int(x.p2(j, nj, kj, t, nt, ht, mm), -nt);
                    for (int hi = 0; hi < p.K[i]; ++hi)
                        for (int nj = 0; nj <= p.F[j]; ++nj)
                            for (int ni = 1; ni <= p.F[i]; ++ni)
                                m.row_add_int(x.p2(j, nj, kj, i, ni, hi, mm), p.N);
                    m.emit_le_int(0);
                }
            }
        }
    }
}

}  // namespace detail

/**
 * Bound the utilization of one queue over the BAS polytope.
 *
 * The utilization is sum over m, k and n >= 1 of p2(i,n,k,i,n,k,m), i.e. the
 * probability that queue i holds at least one job, in ANY blocking
 * configuration. That is the raw utilization; the per-phase EFFECTIVE
 * utilization, which excludes the blocked configurations, is returned in `e`.
 *
 * @param p               network and blocking parameters, all 0-based
 * @param objective_queue queue index, 0..M-1
 * @param sense           Max for an upper bound, Min for a lower bound
 */
template <class T>
QrBasResult<T> mapqn_qr_bounds_bas(const QrBasParams<T>& p, int objective_queue,
                                   MapqnSense sense = MapqnSense::Min) {
    p.validate();
    if (objective_queue < 0 || objective_queue >= p.M)
        throw InputError("qrf_bas: objective_queue out of range");

    const QrBasIndex x(p.M, p.N, p.K, p.MR);
    lp::LpModel<T> m(x.num_vars());
    // No upper bounds are set: the reference leaves ub = inf and bounds the
    // variables through ONE. LpModel's default is exactly lb = 0, ub free.

    // Families in the reference's emission order. ZERO must precede SYMMETRY,
    // which skips pairs both of whose members it pinned.
    const std::vector<char> zero = detail::bas_zero(p, x, m);
    detail::bas_one(p, x, m);
    detail::bas_symmetry(p, x, m, zero);
    detail::bas_marginals(p, x, m);
    detail::bas_ueff(p, x, m);
    detail::bas_thm1(p, x, m);
    detail::bas_thm2(p, x, m);
    detail::bas_cor1(p, x, m);
    detail::bas_thm30(p, x, m);
    detail::bas_thm3(p, x, m);
    detail::bas_thm3f(p, x, m);
    detail::bas_thm3i(p, x, m);
    detail::bas_thm3l(p, x, m);
    detail::bas_thm4(p, x, m);

    // UTILIZATION of the target queue, not occupancy. The cost used to sum the
    // diagonal p2 over ALL configurations, i.e. P(n_i >= 1) with the BLOCKED
    // ones included; a blocked BAS server holds a job it has already finished
    // and does no work, so that is occupancy. On cqn_bas_blocking it gave
    // U = 1 against an exact utilization of 0.590164, and the error propagated
    // into the derived throughput through U = X*V*s.
    //
    // e is already the right quantity -- see bas_ueff, which restricts to the
    // configurations where i is NOT blocked -- and optimising it rather than
    // reading it out afterwards is what keeps the answer a BOUND.
    //
    // The 1/M matches THIS port's UEFF, which emits one row per (i,ki) with j
    // summed INSIDE, leaving e scaled by M. The python and JAR ports emit one
    // row per (j,i,ki) instead and carry no such factor: the scale belongs to
    // the formulation, not to the definition of e.
    const T one = num_traits<T>::from_int(1);
    const T inv_m = one / num_traits<T>::from_int(p.M);
    for (int ki = 0; ki < p.K[objective_queue]; ++ki)
        m.set_cost(x.e(objective_queue, ki), inv_m);
    m.set_maximize(sense == MapqnSense::Max);

    // lp_solve, not simplex_solve: this model outgrows the dense tableau
    // quickly, so a double instantiation hands wide models to HiGHS while
    // Rational stays on the exact path at compile time.
    const lp::LpSolution<T> sol = lp::lp_solve(m);

    QrBasResult<T> out;
    out.status = lp::lp_status_name(sol.status);
    out.ok = sol.ok();
    out.objective = sol.objective;
    out.x = sol.x;
    out.num_vars = m.num_vars();
    out.num_rows = m.num_rows();
    out.iterations = sol.iterations;
    if (!out.ok) return out;

    // U(i) at the optimal vertex, read from the SAME quantity the cost vector
    // optimises so the objective queue's entry is a genuine bound. `qrf_bas.m`
    // reads all M out of ONE solve, which is what the analyzer needs to rebuild
    // the table; the other stations are incidental values at that vertex.
    const T zero_u = num_traits<T>::from_int(0);
    out.U.assign(static_cast<std::size_t>(p.M), zero_u);
    for (int i = 0; i < p.M; ++i)
        for (int ki = 0; ki < p.K[i]; ++ki)
            out.U[static_cast<std::size_t>(i)] += sol.x[x.e(i, ki)] * inv_m;

    out.occupancy.assign(static_cast<std::size_t>(p.M), zero_u);
    for (int i = 0; i < p.M; ++i)
        for (int mm = 0; mm < p.MR; ++mm)
            for (int ki = 0; ki < p.K[i]; ++ki)
                for (int ni = 1; ni <= p.F[i]; ++ni)
                    out.occupancy[static_cast<std::size_t>(i)] +=
                        sol.x[x.p2(i, ni, ki, i, ni, ki, mm)];

    std::size_t maxK = 0;
    for (int i = 0; i < p.M; ++i)
        if (static_cast<std::size_t>(p.K[i]) > maxK) maxK = static_cast<std::size_t>(p.K[i]);
    out.e = Matrix<T>(static_cast<std::size_t>(p.M), maxK);
    for (int i = 0; i < p.M; ++i)
        for (int ki = 0; ki < p.K[i]; ++ki)
            out.e(static_cast<std::size_t>(i), static_cast<std::size_t>(ki)) = sol.x[x.e(i, ki)];
    return out;
}

}  // namespace mapqn
}  // namespace line

#endif  // LINE_API_MAPQN_MAPQN_QR_BOUNDS_BAS_H
