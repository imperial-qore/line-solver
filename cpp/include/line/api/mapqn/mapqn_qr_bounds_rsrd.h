/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAPQN_MAPQN_QR_BOUNDS_RSRD_H
#define LINE_API_MAPQN_MAPQN_QR_BOUNDS_RSRD_H

/**
 * Quadratic-reduction bound on the utilization of one queue of a closed MAP
 * queueing network under RS-RD blocking (repetitive service, random
 * destination).
 *
 * Templated port of matlab/lib/qrf/qrf_rsrd.m (ground truth), whose own origin
 * is the AMPL model qrboundsrsrd_skel.mod.
 *
 * WHAT RS-RD MEANS HERE. A job completing at i and routed to a FULL destination
 * is not held: service is repeated at i and a fresh destination drawn. Nothing
 * is blocked in the BAS sense, so unlike mapqn_qr_bounds_bas there is no
 * blocking-configuration index and the model is far narrower. What the blocking
 * costs instead is EFFECTIVE service: Ueff(i,k,n) subtracts the mass whose
 * chosen destination was full, and pb(i) accumulates the difference.
 *
 * ITS HALF-INDEX IS NOT THE FAMILY'S. Populations run 0..F(i), not 0..N, so the
 * pairwise half-index is built on F(i)+1 and this model cannot share a layout
 * with any other mapqn bound. Every population loop below is capped by F, and
 * that is deliberate, not an optimization.
 *
 * THE ONE TRAP, and it is invisible in the output. THM1 is AGGREGATED over
 * nj >= 1: ONE row per (j,kj), not one per (j,kj,nj). The per-nj form is
 * strictly stronger, solves cleanly, and over-tightens the polytope -- on the
 * paper's M = 5, N = 20 instance it returns U1min = 0.92508 against the
 * published 0.87058. THM1c is the separate nj = 0 row and is NOT part of that
 * aggregation. The reference carries this warning in its own comment; it is
 * repeated here because a porter reading only the loop nest would not see it.
 *
 * LOAD DEPENDENCE IS LIVE. Unlike qrf_bas, q here carries alpha(i,n), so the
 * rate depends on the population at the departing queue. The reference indexes
 * it as q{i,j}(k,h,n+1), a 1-based population axis whose first slot is
 * population 0; rsrd_rate below takes the population directly.
 *
 * ARITHMETIC. Assembly is +, -, * on the model data and lp::simplex_solve uses
 * Bland's rule with no tolerance, so at T = line::Rational the returned bound
 * is the EXACT optimum of the exact polytope.
 *
 * COST. B^2 + 2 sum_i K(i) F(i) + M columns with B = sum_i (F(i)+1) K(i). No
 * MR factor, so it is much smaller than the BAS model at equal size, but the
 * tableau is still dense: see the ceiling recorded in _kb/03-api-layer.md.
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

/** Parameters of the RS-RD bound, mirroring the reference's `params`. */
template <class T>
struct QrRsrdParams {
    int M = 0;                    ///< number of queues
    int N = 0;                    ///< total population
    std::vector<int> F;           ///< (M) capacity of each queue
    std::vector<int> K;           ///< (M) number of phases of each queue
    std::vector<Matrix<T>> mu;    ///< mu[i] is K(i) x K(i), completion rates
    std::vector<Matrix<T>> v;     ///< v[i] is K(i) x K(i), background rates
    Matrix<T> r;                  ///< (M x M) routing probabilities
    std::vector<std::vector<T> > alpha;  ///< optional (M) x (N+1) load scalings; empty means 1

    void validate() const {
        if (M <= 0) throw InputError("qrf_rsrd: M must be positive");
        if (N < 1) throw InputError("qrf_rsrd: N must be at least 1");
        if (static_cast<int>(F.size()) != M) throw InputError("qrf_rsrd: F has the wrong length");
        if (static_cast<int>(K.size()) != M) throw InputError("qrf_rsrd: K has the wrong length");
        if (static_cast<int>(mu.size()) != M || static_cast<int>(v.size()) != M)
            throw InputError("qrf_rsrd: mu and v must have one entry per queue");
        for (int i = 0; i < M; ++i) {
            if (K[i] <= 0) throw InputError("qrf_rsrd: every queue needs at least one phase");
            if (F[i] < 1 || F[i] > N) throw InputError("qrf_rsrd: F(i) must lie in 1..N");
            const std::size_t k = static_cast<std::size_t>(K[i]);
            if (mu[i].rows() != k || mu[i].cols() != k)
                throw InputError("qrf_rsrd: mu{i} must be K(i) x K(i)");
            if (v[i].rows() != k || v[i].cols() != k)
                throw InputError("qrf_rsrd: v{i} must be K(i) x K(i)");
        }
        if (r.rows() != static_cast<std::size_t>(M) || r.cols() != static_cast<std::size_t>(M))
            throw InputError("qrf_rsrd: r must be M x M");
        if (!alpha.empty() && static_cast<int>(alpha.size()) != M)
            throw InputError("qrf_rsrd: alpha must have M rows when given");
    }
};

/** Result of an RS-RD bound solve. */
template <class T>
struct QrRsrdResult {
    bool ok = false;
    std::string status;
    T objective = T();    ///< the bound on the utilization of the target queue
    std::vector<T> U;     ///< (M) utilization of each queue
    std::vector<T> Ueff;  ///< (M) effective utilization of each queue
    std::vector<T> pb;    ///< (M) blocking probability of each queue
    std::vector<T> x;
    std::size_t num_vars = 0;
    std::size_t num_rows = 0;
    std::size_t iterations = 0;
};

/**
 * Variable layout: p2(j,nj,kj,i,ni,hi), then U(i,k,n) and Ueff(i,k,n) for
 * n >= 1, then pb(i).
 *
 * The half-index runs over (queue, population 0..F, phase), so its stride is
 * F(i)+1 rather than the N+1 the rest of the mapqn family uses.
 */
struct QrRsrdIndex {
    int M = 0, N = 0;
    std::vector<int> K, F, base, cumU;
    std::size_t B = 0, off_U = 0, off_Ueff = 0, off_pb = 0, total = 0;

    QrRsrdIndex() {}
    QrRsrdIndex(int m, int n, const std::vector<int>& k, const std::vector<int>& f)
        : M(m), N(n), K(k), F(f) {
        base.assign(static_cast<std::size_t>(M) + 1, 0);
        cumU.assign(static_cast<std::size_t>(M) + 1, 0);
        for (int i = 0; i < M; ++i) {
            base[i + 1] = base[i] + (F[i] + 1) * K[i];
            cumU[i + 1] = cumU[i] + K[i] * F[i];
        }
        B = static_cast<std::size_t>(base[M]);
        off_U = B * B;
        off_Ueff = off_U + static_cast<std::size_t>(cumU[M]);
        off_pb = off_Ueff + static_cast<std::size_t>(cumU[M]);
        total = off_pb + static_cast<std::size_t>(M);
    }

    std::size_t half(int i, int ni, int h) const {
        return static_cast<std::size_t>(base[i]) +
               static_cast<std::size_t>(ni) * static_cast<std::size_t>(K[i]) +
               static_cast<std::size_t>(h);
    }
    std::size_t p2(int j, int nj, int kj, int i, int ni, int hi) const {
        return half(j, nj, kj) * B + half(i, ni, hi);
    }
    /** n is 1-based here: U is only defined for a busy queue. */
    std::size_t U(int i, int k, int n) const {
        return off_U + static_cast<std::size_t>(cumU[i]) +
               static_cast<std::size_t>(k) * static_cast<std::size_t>(F[i]) +
               static_cast<std::size_t>(n - 1);
    }
    std::size_t Ueff(int i, int k, int n) const {
        return off_Ueff + static_cast<std::size_t>(cumU[i]) +
               static_cast<std::size_t>(k) * static_cast<std::size_t>(F[i]) +
               static_cast<std::size_t>(n - 1);
    }
    std::size_t pb(int i) const { return off_pb + static_cast<std::size_t>(i); }
    std::size_t num_vars() const { return total; }
};

namespace detail {

/** q(i,j,k,h,n): the load-dependent rate at population n of queue i. */
template <class T>
T rsrd_rate(const QrRsrdParams<T>& p, int i, int j, int k, int h, int n) {
    const std::size_t ki = static_cast<std::size_t>(k), hi = static_cast<std::size_t>(h);
    const std::size_t ii = static_cast<std::size_t>(i), ji = static_cast<std::size_t>(j);
    T a = num_traits<T>::from_int(1);
    if (!p.alpha.empty() && n >= 0 &&
        static_cast<std::size_t>(n) < p.alpha[static_cast<std::size_t>(i)].size())
        a = p.alpha[static_cast<std::size_t>(i)][static_cast<std::size_t>(n)];
    if (j != i) return T(a * p.r(ii, ji) * p.mu[i](ki, hi));
    return T(a * (p.v[i](ki, hi) + p.r(ii, ii) * p.mu[i](ki, hi)));
}

/** Bounds: every variable lies in [0,1]. */
template <class T>
void rsrd_bounds(const QrRsrdParams<T>& p, const QrRsrdIndex& x, lp::LpModel<T>& m) {
    const T one = num_traits<T>::from_int(1);
    for (std::size_t j = 0; j < x.num_vars(); ++j) m.set_bounds(j, T(), one);
}

/** ONE: the diagonal marginal of each queue normalizes to one. */
template <class T>
void rsrd_one(const QrRsrdParams<T>& p, const QrRsrdIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int nj = 0; nj <= p.F[j]; ++nj)
            for (int kj = 0; kj < p.K[j]; ++kj) m.row_add_int(x.p2(j, nj, kj, j, nj, kj), 1);
        m.emit_eq_int(1);
    }
}

/**
 * ZERO1/2/3/6/7: the states that carry no mass, pinned as ub = 0.
 *
 * ZERO6 and ZERO7 are the capacity-aware members: a pair (or a single queue)
 * holding so few jobs that the REMAINING capacity of the network cannot absorb
 * the rest of the population is impossible.
 */
template <class T>
void rsrd_zero(const QrRsrdParams<T>& p, const QrRsrdIndex& x, lp::LpModel<T>& m) {
    int totF = 0;
    for (int i = 0; i < p.M; ++i) totF += p.F[i];
    for (int j = 0; j < p.M; ++j) {
        for (int nj = 0; nj <= p.F[j]; ++nj) {
            for (int kj = 0; kj < p.K[j]; ++kj) {
                for (int i = 0; i < p.M; ++i) {
                    for (int ni = 0; ni <= p.F[i]; ++ni) {
                        for (int hi = 0; hi < p.K[i]; ++hi) {
                            bool z = false;
                            if (i == j && nj == ni && hi != kj) z = true;  // ZERO1
                            if (i == j && nj != ni) z = true;              // ZERO2
                            if (i != j && nj + ni > p.N) z = true;         // ZERO3
                            if (i != j && p.N - nj - ni > totF - p.F[i] - p.F[j])
                                z = true;  // ZERO6
                            if (z) m.fix(x.p2(j, nj, kj, i, ni, hi), T());
                        }
                    }
                }
            }
        }
        for (int nj = 0; nj <= p.F[j]; ++nj)  // ZERO7
            for (int kj = 0; kj < p.K[j]; ++kj)
                if (p.N - nj > totF - p.F[j]) m.fix(x.p2(j, nj, kj, j, nj, kj), T());
    }
}

/** SYMMETRY: p2 is symmetric under swapping its halves. */
template <class T>
void rsrd_symmetry(const QrRsrdParams<T>& p, const QrRsrdIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j)
        for (int nj = 0; nj <= p.F[j]; ++nj)
            for (int kj = 0; kj < p.K[j]; ++kj)
                for (int i = j + 1; i < p.M; ++i)
                    for (int ni = 0; ni <= p.F[i]; ++ni)
                        for (int hi = 0; hi < p.K[i]; ++hi) {
                            const std::size_t a = x.p2(j, nj, kj, i, ni, hi);
                            const std::size_t b = x.p2(i, ni, hi, j, nj, kj);
                            if (a == b) continue;
                            m.row_add_int(a, 1);
                            m.row_add_int(b, -1);
                            m.emit_eq_int(0);
                        }
}

/** MARGINALS: the pairwise law agrees with its own marginal at every level. */
template <class T>
void rsrd_marginals(const QrRsrdParams<T>& p, const QrRsrdIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j)
        for (int kj = 0; kj < p.K[j]; ++kj)
            for (int nj = 0; nj <= p.F[j]; ++nj)
                for (int i = 0; i < p.M; ++i) {
                    if (i == j) continue;
                    m.row_add_int(x.p2(j, nj, kj, j, nj, kj), 1);
                    for (int ni = 0; ni <= p.F[i]; ++ni)
                        for (int hi = 0; hi < p.K[i]; ++hi)
                            m.row_add_int(x.p2(j, nj, kj, i, ni, hi), -1);
                    m.emit_eq_int(0);
                }
}

/** UCLASSIC: U(i,k,n) is the diagonal of p2, the classical utilization. */
template <class T>
void rsrd_uclassic(const QrRsrdParams<T>& p, const QrRsrdIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i)
        for (int ki = 0; ki < p.K[i]; ++ki)
            for (int ni = 1; ni <= p.F[i]; ++ni) {
                m.row_add_int(x.U(i, ki, ni), 1);
                m.row_add_int(x.p2(i, ni, ki, i, ni, ki), -1);
                m.emit_eq_int(0);
            }
}

/**
 * UEFFS: the effective utilization removes the mass whose chosen destination is
 * full, weighted by the routing probability of choosing it.
 */
template <class T>
void rsrd_ueffs(const QrRsrdParams<T>& p, const QrRsrdIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i)
        for (int ki = 0; ki < p.K[i]; ++ki)
            for (int ni = 1; ni <= p.F[i]; ++ni) {
                m.row_add_int(x.Ueff(i, ki, ni), 1);
                m.row_add_int(x.p2(i, ni, ki, i, ni, ki), -1);
                for (int j = 0; j < p.M; ++j) {
                    if (j == i) continue;
                    if (!(p.r(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) > T()))
                        continue;
                    for (int hj = 0; hj < p.K[j]; ++hj)
                        m.row_add(x.p2(i, ni, ki, j, p.F[j], hj),
                                  p.r(static_cast<std::size_t>(i), static_cast<std::size_t>(j)));
                }
                m.emit_eq_int(0);
            }
}

/** PBLOCK: pb(i) is the total mass lost to full destinations. */
template <class T>
void rsrd_pblock(const QrRsrdParams<T>& p, const QrRsrdIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        m.row_add_int(x.pb(i), 1);
        for (int ki = 0; ki < p.K[i]; ++ki)
            for (int ni = 1; ni <= p.F[i]; ++ni) {
                m.row_add_int(x.U(i, ki, ni), -1);
                m.row_add_int(x.Ueff(i, ki, ni), 1);
            }
        m.emit_eq_int(0);
    }
}

/** PBB: pb(i) cannot exceed the probability that some destination is full. */
template <class T>
void rsrd_pbb(const QrRsrdParams<T>& p, const QrRsrdIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        m.row_add_int(x.pb(i), 1);
        for (int j = 0; j < p.M; ++j) {
            if (j == i) continue;
            if (!(p.r(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) > T())) continue;
            for (int hj = 0; hj < p.K[j]; ++hj)
                m.row_add_int(x.p2(j, p.F[j], hj, j, p.F[j], hj), -1);
        }
        m.emit_le_int(0);
    }
}

/**
 * THM2: phase balance on the EFFECTIVE utilizations for the routing-away terms
 * and on the classical diagonal for the self-routing ones. A repeated service
 * does not change the phase balance, which is why the two levels appear in one
 * row.
 */
template <class T>
void rsrd_thm2(const QrRsrdParams<T>& p, const QrRsrdIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int ki = 0; ki < p.K[i]; ++ki) {
            for (int ni = 1; ni <= p.F[i]; ++ni) {
                for (int j = 0; j < p.M; ++j) {
                    if (j == i) continue;
                    for (int hi = 0; hi < p.K[i]; ++hi) {
                        if (hi == ki) continue;
                        m.row_add(x.Ueff(i, ki, ni), rsrd_rate(p, i, j, ki, hi, ni));
                    }
                }
                for (int hi = 0; hi < p.K[i]; ++hi) {
                    if (hi == ki) continue;
                    m.row_add(x.p2(i, ni, ki, i, ni, ki), rsrd_rate(p, i, i, ki, hi, ni));
                }
            }
            for (int ni = 1; ni <= p.F[i]; ++ni) {
                for (int j = 0; j < p.M; ++j) {
                    if (j == i) continue;
                    for (int hi = 0; hi < p.K[i]; ++hi) {
                        if (hi == ki) continue;
                        m.row_add(x.Ueff(i, hi, ni), T(-rsrd_rate(p, i, j, hi, ki, ni)));
                    }
                }
                for (int hi = 0; hi < p.K[i]; ++hi) {
                    if (hi == ki) continue;
                    m.row_add(x.p2(i, ni, hi, i, ni, hi), T(-rsrd_rate(p, i, i, hi, ki, ni)));
                }
            }
            m.emit_eq_int(0);
        }
    }
}

/**
 * THM1: the queue-length theorem, AGGREGATED over nj >= 1.
 *
 * ONE row per (j,kj). Emitting it per-nj is strictly stronger and silently
 * over-tightens the polytope; see the header.
 */
template <class T>
void rsrd_thm1(const QrRsrdParams<T>& p, const QrRsrdIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int kj = 0; kj < p.K[j]; ++kj) {
            for (int nj = 1; nj <= p.F[j]; ++nj) {
                m.row_add_int(x.p2(j, nj, kj, j, nj, kj), -p.N);
                for (int i = 0; i < p.M; ++i)
                    for (int ni = 1; ni <= p.F[i]; ++ni)
                        for (int hi = 0; hi < p.K[i]; ++hi)
                            m.row_add_int(x.p2(j, nj, kj, i, ni, hi), ni);
            }
            m.emit_eq_int(0);
        }
    }
}

/** THM1c: the same theorem at nj = 0, where only i != j can hold jobs. */
template <class T>
void rsrd_thm1c(const QrRsrdParams<T>& p, const QrRsrdIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int kj = 0; kj < p.K[j]; ++kj) {
            m.row_add_int(x.p2(j, 0, kj, j, 0, kj), -p.N);
            for (int i = 0; i < p.M; ++i) {
                if (i == j) continue;
                for (int ni = 1; ni <= p.F[i]; ++ni)
                    for (int hi = 0; hi < p.K[i]; ++hi)
                        m.row_add_int(x.p2(j, 0, kj, i, ni, hi), ni);
            }
            m.emit_eq_int(0);
        }
    }
}

/** THM3a: level-crossing balance between ni and ni+1, for 1 <= ni <= F(i)-1. */
template <class T>
void rsrd_thm3a(const QrRsrdParams<T>& p, const QrRsrdIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int ni = 1; ni <= p.F[i] - 1; ++ni) {
            for (int j = 0; j < p.M; ++j) {
                if (j == i) continue;
                for (int kj = 0; kj < p.K[j]; ++kj)
                    for (int hj = 0; hj < p.K[j]; ++hj)
                        for (int ui = 0; ui < p.K[i]; ++ui)
                            for (int nj = 1; nj <= p.F[j]; ++nj)
                                m.row_add(x.p2(j, nj, kj, i, ni, ui),
                                          rsrd_rate(p, j, i, kj, hj, nj));
            }
            for (int j = 0; j < p.M; ++j) {
                if (j == i) continue;
                for (int ki = 0; ki < p.K[i]; ++ki)
                    for (int hi = 0; hi < p.K[i]; ++hi)
                        for (int uj = 0; uj < p.K[j]; ++uj)
                            for (int nj = 0; nj <= p.F[j] - 1; ++nj)
                                m.row_add(x.p2(i, ni + 1, ki, j, nj, uj),
                                          T(-rsrd_rate(p, i, j, ki, hi, ni + 1)));
            }
            m.emit_eq_int(0);
        }
    }
}

/** THM3b: the same balance at ni = 0, per arrival phase. */
template <class T>
void rsrd_thm3b(const QrRsrdParams<T>& p, const QrRsrdIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int ui = 0; ui < p.K[i]; ++ui) {
            for (int j = 0; j < p.M; ++j) {
                if (j == i) continue;
                for (int kj = 0; kj < p.K[j]; ++kj)
                    for (int hj = 0; hj < p.K[j]; ++hj)
                        for (int nj = 1; nj <= p.F[j]; ++nj)
                            m.row_add(x.p2(j, nj, kj, i, 0, ui), rsrd_rate(p, j, i, kj, hj, nj));
            }
            for (int j = 0; j < p.M; ++j) {
                if (j == i) continue;
                for (int ki = 0; ki < p.K[i]; ++ki)
                    for (int nj = 0; nj <= p.F[j] - 1; ++nj)
                        for (int hj = 0; hj < p.K[j]; ++hj)
                            m.row_add(x.p2(i, 1, ki, j, nj, hj),
                                      T(-rsrd_rate(p, i, j, ki, ui, 1)));
            }
            m.emit_eq_int(0);
        }
    }
}

/** QBAL: queue balance, six terms, three on each side. */
template <class T>
void rsrd_qbal(const QrRsrdParams<T>& p, const QrRsrdIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int ki = 0; ki < p.K[i]; ++ki) {
            for (int hi = 0; hi < p.K[i]; ++hi) {  // LHS 1
                if (hi == ki) continue;
                for (int j = 0; j < p.M; ++j) {
                    if (j == i) continue;
                    for (int ni = 1; ni <= p.F[i]; ++ni)
                        for (int uj = 0; uj < p.K[j]; ++uj)
                            for (int nj = 0; nj <= p.F[j] - 1; ++nj)
                                m.row_add(x.p2(i, ni, ki, j, nj, uj),
                                          T(rsrd_rate(p, i, j, ki, hi, ni) *
                                            num_traits<T>::from_int(ni)));
                }
            }
            for (int hi = 0; hi < p.K[i]; ++hi) {  // LHS 2
                if (hi == ki) continue;
                for (int ni = 1; ni <= p.F[i]; ++ni)
                    m.row_add(x.p2(i, ni, ki, i, ni, ki),
                              T(rsrd_rate(p, i, i, ki, hi, ni) * num_traits<T>::from_int(ni)));
            }
            for (int j = 0; j < p.M; ++j) {  // LHS 3
                if (j == i) continue;
                for (int hi = 0; hi < p.K[i]; ++hi)
                    for (int ni = 1; ni <= p.F[i]; ++ni)
                        for (int uj = 0; uj < p.K[j]; ++uj) {
                            const int cap =
                                (p.F[j] - 1) < (p.N - ni) ? (p.F[j] - 1) : (p.N - ni);
                            for (int nj = 0; nj <= cap; ++nj)
                                m.row_add(x.p2(i, ni, hi, j, nj, uj),
                                          rsrd_rate(p, i, j, hi, ki, ni));
                        }
            }
            for (int j = 0; j < p.M; ++j) {  // RHS 1
                if (j == i) continue;
                for (int hj = 0; hj < p.K[j]; ++hj)
                    for (int ni = 0; ni <= p.F[i] - 1; ++ni)
                        for (int uj = 0; uj < p.K[j]; ++uj)
                            for (int nj = 1; nj <= p.F[j]; ++nj)
                                m.row_add(x.p2(i, ni, ki, j, nj, hj),
                                          T(-rsrd_rate(p, j, i, hj, uj, nj)));
            }
            for (int hi = 0; hi < p.K[i]; ++hi) {  // RHS 2
                if (hi == ki) continue;
                for (int ni = 1; ni <= p.F[i]; ++ni)
                    m.row_add(x.p2(i, ni, hi, i, ni, hi),
                              T(-(rsrd_rate(p, i, i, hi, ki, ni) * num_traits<T>::from_int(ni))));
            }
            for (int hi = 0; hi < p.K[i]; ++hi) {  // RHS 3
                if (hi == ki) continue;
                for (int j = 0; j < p.M; ++j) {
                    if (j == i) continue;
                    for (int ni = 1; ni <= p.F[i]; ++ni)
                        for (int uj = 0; uj < p.K[j]; ++uj)
                            for (int nj = 0; nj <= p.F[j] - 1; ++nj)
                                m.row_add(x.p2(i, ni, hi, j, nj, uj),
                                          T(-(rsrd_rate(p, i, j, hi, ki, ni) *
                                              num_traits<T>::from_int(ni))));
                }
            }
            m.emit_eq_int(0);
        }
    }
}

/**
 * THM4: the QMIN inequality. The reference accumulates
 * (sum nt p2 - N sum p2) and emits the NEGATED row, so the constraint is
 * N P(j at (nj,kj), i nonempty) <= sum_t E[n_t ...].
 */
template <class T>
void rsrd_thm4(const QrRsrdParams<T>& p, const QrRsrdIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j)
        for (int kj = 0; kj < p.K[j]; ++kj)
            for (int i = 0; i < p.M; ++i) {
                for (int t = 0; t < p.M; ++t)
                    for (int ht = 0; ht < p.K[t]; ++ht)
                        for (int nj = 0; nj <= p.F[j]; ++nj)
                            for (int nt = 1; nt <= p.F[t]; ++nt)
                                m.row_add_int(x.p2(j, nj, kj, t, nt, ht), -nt);
                for (int hi = 0; hi < p.K[i]; ++hi)
                    for (int nj = 0; nj <= p.F[j]; ++nj)
                        for (int ni = 1; ni <= p.F[i]; ++ni)
                            m.row_add_int(x.p2(j, nj, kj, i, ni, hi), p.N);
                m.emit_le_int(0);
            }
}

}  // namespace detail

/**
 * Bound the utilization of one queue over the RS-RD polytope.
 *
 * @param p               network parameters, all 0-based
 * @param objective_queue queue index, 0..M-1
 * @param sense           Max for an upper bound, Min for a lower bound
 */
template <class T>
QrRsrdResult<T> mapqn_qr_bounds_rsrd(const QrRsrdParams<T>& p, int objective_queue,
                                     MapqnSense sense = MapqnSense::Min) {
    p.validate();
    if (objective_queue < 0 || objective_queue >= p.M)
        throw InputError("qrf_rsrd: objective_queue out of range");

    const QrRsrdIndex x(p.M, p.N, p.K, p.F);
    lp::LpModel<T> m(x.num_vars());
    detail::rsrd_bounds(p, x, m);

    detail::rsrd_one(p, x, m);
    detail::rsrd_zero(p, x, m);
    detail::rsrd_symmetry(p, x, m);
    detail::rsrd_marginals(p, x, m);
    detail::rsrd_uclassic(p, x, m);
    detail::rsrd_ueffs(p, x, m);
    detail::rsrd_pblock(p, x, m);
    detail::rsrd_pbb(p, x, m);
    detail::rsrd_thm2(p, x, m);
    detail::rsrd_thm1(p, x, m);
    detail::rsrd_thm1c(p, x, m);
    detail::rsrd_thm3a(p, x, m);
    detail::rsrd_thm3b(p, x, m);
    detail::rsrd_qbal(p, x, m);
    detail::rsrd_thm4(p, x, m);

    const T one = num_traits<T>::from_int(1);
    for (int ki = 0; ki < p.K[objective_queue]; ++ki)
        for (int ni = 1; ni <= p.F[objective_queue]; ++ni)
            m.set_cost(x.p2(objective_queue, ni, ki, objective_queue, ni, ki), one);
    m.set_maximize(sense == MapqnSense::Max);

    // lp_solve, not simplex_solve: this model outgrows the dense tableau
    // quickly, so a double instantiation hands wide models to HiGHS while
    // Rational stays on the exact path at compile time.
    const lp::LpSolution<T> sol = lp::lp_solve(m);

    QrRsrdResult<T> out;
    out.status = lp::lp_status_name(sol.status);
    out.ok = sol.ok();
    out.objective = sol.objective;
    out.x = sol.x;
    out.num_vars = m.num_vars();
    out.num_rows = m.num_rows();
    out.iterations = sol.iterations;
    if (!out.ok) return out;

    out.U.assign(static_cast<std::size_t>(p.M), T());
    out.Ueff.assign(static_cast<std::size_t>(p.M), T());
    out.pb.assign(static_cast<std::size_t>(p.M), T());
    for (int i = 0; i < p.M; ++i) {
        for (int ki = 0; ki < p.K[i]; ++ki) {
            for (int ni = 1; ni <= p.F[i]; ++ni) {
                out.U[static_cast<std::size_t>(i)] += sol.x[x.U(i, ki, ni)];
                out.Ueff[static_cast<std::size_t>(i)] += sol.x[x.Ueff(i, ki, ni)];
            }
        }
        out.pb[static_cast<std::size_t>(i)] = sol.x[x.pb(i)];
    }
    return out;
}

}  // namespace mapqn
}  // namespace line

#endif  // LINE_API_MAPQN_MAPQN_QR_BOUNDS_RSRD_H
