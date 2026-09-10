/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAPQN_MAPQN_P1_COMMON_H
#define LINE_API_MAPQN_MAPQN_P1_COMMON_H

/**
 * The variable layout and the constraint families shared by the two reductions
 * that carry singly-indexed probabilities: the linear reduction
 * (mapqn_bnd_lr.m) and the general quadratic reduction (mapqn_bnd_qr.m).
 *
 * The two references duplicate these fifteen families verbatim -- ZER1..ZER4,
 * CEQU, ONE1, UTLB, UTLC, QLEN, CLEN, ONE, POPC, SRVB, UUB1, QUB1 are the same
 * loops with the same coefficients in both files -- and the QR file then adds
 * the p2 level on top. Emitting each family from ONE function here is the same
 * discipline mapqn_qr_common.h applies to the load-dependent pair, and for the
 * same reason: a correction to a family must not be able to land in one entry
 * point and miss the other.
 *
 * INDEX LAYOUT. Blocks in the order the references register them, minus the
 * three inert ones (see below). Populations are 0..N as written; queues and
 * phases are 0-based, so every reference loop `for i = 1:M` becomes
 * `for i = 0; i < M; ++i` and every subscript except a population drops by one.
 * `half(i,ni,h)` is shared by the inner half of p1 and p1c and by both halves
 * of p2, which is what makes the QR projection families PI21/PI22 plain index
 * arithmetic and PI23 a swap.
 *
 * INERT VARIABLES OMITTED. Both references register UP(j,k,i,h), QP(j,k,i,h)
 * and I_var(j,k,i), give them upper bounds, and then reference them in no
 * constraint and in no objective (verified in both files: the three index
 * arrays appear only at their own registration and bounding). A variable with
 * no row and no cost cannot move the optimum, so they are not allocated.
 */

#include <cstddef>
#include <vector>

#include "line/api/mapqn/mapqn_params.h"
#include "line/num/number.h"
#include "line/util/simplex.h"

namespace line {
namespace mapqn {

/**
 * Variable layout of the p1-level models.
 *
 * `with_p2` selects the quadratic reduction's extra block; the linear reduction
 * carries no joint variables at all, and allocating them would square the model
 * for nothing.
 */
struct MapqnP1Index {
    int M = 0, N = 0;
    std::vector<int> K;
    std::vector<int> cumK;   ///< cumK[i] = sum_{i' < i} K(i')
    std::size_t sumK = 0;    ///< sum_i K(i)
    std::size_t block = 0;   ///< (N+1) * sumK, the pairwise half-index range
    bool has_p2 = false;
    std::size_t off_U = 0, off_IT = 0, off_Q = 0, off_C = 0;
    std::size_t off_p1 = 0, off_p1c = 0, off_p2 = 0, total = 0;

    MapqnP1Index() {}
    MapqnP1Index(int m, int n, const std::vector<int>& k, bool with_p2)
        : M(m), N(n), K(k), has_p2(with_p2) {
        cumK.assign(static_cast<std::size_t>(M) + 1, 0);
        for (int i = 0; i < M; ++i) cumK[i + 1] = cumK[i] + K[i];
        sumK = static_cast<std::size_t>(cumK[M]);
        block = static_cast<std::size_t>(N + 1) * sumK;
        off_U = 0;
        off_IT = off_U + sumK;
        off_Q = off_IT + sumK;
        off_C = off_Q + sumK;
        off_p1 = off_C + sumK * static_cast<std::size_t>(M);
        off_p1c = off_p1 + sumK * block;
        off_p2 = off_p1c + sumK * block;
        total = has_p2 ? off_p2 + block * block : off_p2;
    }

    /** Flat (queue, phase) index, range sumK. */
    std::size_t pk(int i, int k) const {
        return static_cast<std::size_t>(cumK[i]) + static_cast<std::size_t>(k);
    }
    /** Flat (queue, population, phase) half-index, range block. */
    std::size_t half(int i, int ni, int h) const {
        return static_cast<std::size_t>(N + 1) * static_cast<std::size_t>(cumK[i]) +
               static_cast<std::size_t>(ni) * static_cast<std::size_t>(K[i]) +
               static_cast<std::size_t>(h);
    }

    std::size_t U(int i, int k) const { return off_U + pk(i, k); }
    std::size_t IT(int i, int k) const { return off_IT + pk(i, k); }
    std::size_t Q(int i, int k) const { return off_Q + pk(i, k); }
    std::size_t C(int j, int k, int i) const {
        return off_C + pk(j, k) * static_cast<std::size_t>(M) + static_cast<std::size_t>(i);
    }
    std::size_t p1(int j, int k, int i, int ni, int h) const {
        return off_p1 + pk(j, k) * block + half(i, ni, h);
    }
    std::size_t p1c(int j, int k, int i, int ni, int h) const {
        return off_p1c + pk(j, k) * block + half(i, ni, h);
    }
    std::size_t p2(int j, int nj, int k, int i, int ni, int h) const {
        return off_p2 + half(j, nj, k) * block + half(i, ni, h);
    }

    std::size_t num_vars() const { return total; }
};

namespace detail {

/**
 * q(i,j,k,h): rate at which queue i in phase k moves to phase h while routing a
 * job to queue j.
 *
 * Port of the q{i,j}(k,h) cell built at the head of both references. It is
 * deliberately NOT mapqn_q from mapqn_params.h: neither model has load
 * dependence, so there is no population argument and no alpha factor, and
 * mapqn_q additionally returns 0 at n == 0, which has no counterpart here. The
 * i == j branch adds the background rate v, because a self-routing completion
 * and a phase change without completion are indistinguishable in the marginal.
 */
template <class T>
T qr_rate(const MapqnParams<T>& p, int i, int j, int k, int h) {
    const std::size_t ki = static_cast<std::size_t>(k), hi = static_cast<std::size_t>(h);
    const std::size_t ii = static_cast<std::size_t>(i), ji = static_cast<std::size_t>(j);
    if (j != i) return T(p.r(ii, ji) * p.mu[i](ki, hi));
    return T(p.v[i](ki, hi) + p.r(ii, ii) * p.mu[i](ki, hi));
}

/** Default bounds on the p1-level variables, before ZER pins the zeros. */
template <class T>
void p1_bounds(const MapqnParams<T>& p, const MapqnP1Index& x, lp::LpModel<T>& m) {
    const T one = num_traits<T>::from_int(1);
    const T nn = num_traits<T>::from_int(p.N);
    for (int i = 0; i < p.M; ++i) {
        for (int k = 0; k < p.K[i]; ++k) {
            m.set_bounds(x.U(i, k), T(), one);
            m.set_bounds(x.IT(i, k), T(), one);
            m.set_bounds(x.Q(i, k), T(), nn);
            for (int t = 0; t < p.M; ++t) m.set_bounds(x.C(i, k, t), T(), nn);
            for (int t = 0; t < p.M; ++t) {
                for (int nt = 0; nt <= p.N; ++nt) {
                    for (int h = 0; h < p.K[t]; ++h) {
                        m.set_bounds(x.p1(i, k, t, nt, h), T(), one);
                        m.set_bounds(x.p1c(i, k, t, nt, h), T(), one);
                    }
                }
            }
        }
    }
}

/** ZER1: p1(j,k,j,0,k) = 0. An occupied station cannot hold zero jobs. */
template <class T>
void p1_zer1(const MapqnParams<T>& p, const MapqnP1Index& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j)
        for (int k = 0; k < p.K[j]; ++k) m.fix(x.p1(j, k, j, 0, k), T());
}

/** ZER2: p1(j,k,j,nj,h) = 0 for h != k. One station, one phase. */
template <class T>
void p1_zer2(const MapqnParams<T>& p, const MapqnP1Index& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j)
        for (int k = 0; k < p.K[j]; ++k)
            for (int nj = 0; nj <= p.N; ++nj)
                for (int h = 0; h < p.K[j]; ++h)
                    if (h != k) m.fix(x.p1(j, k, j, nj, h), T());
}

/** ZER3: p1(j,k,i,N,h) = 0 for i != j. A busy j leaves i short of N. */
template <class T>
void p1_zer3(const MapqnParams<T>& p, const MapqnP1Index& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j)
        for (int k = 0; k < p.K[j]; ++k)
            for (int i = 0; i < p.M; ++i)
                if (i != j)
                    for (int h = 0; h < p.K[i]; ++h) m.fix(x.p1(j, k, i, p.N, h), T());
}

/** ZER4: p1c(j,k,j,nj,h) = 0 for nj >= 1. The complement is the idle branch. */
template <class T>
void p1_zer4(const MapqnParams<T>& p, const MapqnP1Index& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j)
        for (int k = 0; k < p.K[j]; ++k)
            for (int nj = 1; nj <= p.N; ++nj)
                for (int h = 0; h < p.K[j]; ++h) m.fix(x.p1c(j, k, j, nj, h), T());
}

/** CEQU: C(j,k,j) = Q(j,k). The self-conditioned length is the mean length. */
template <class T>
void p1_cequ(const MapqnParams<T>& p, const MapqnP1Index& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int k = 0; k < p.K[j]; ++k) {
            m.row_add_int(x.C(j, k, j), 1);
            m.row_add_int(x.Q(j, k), -1);
            m.emit_eq_int(0);
        }
    }
}

/** ONE1: sum over k,h,ni of (p1 + p1c) = 1 for each ordered pair (j,i). */
template <class T>
void p1_one1(const MapqnParams<T>& p, const MapqnP1Index& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int i = 0; i < p.M; ++i) {
            for (int k = 0; k < p.K[j]; ++k) {
                for (int h = 0; h < p.K[i]; ++h) {
                    for (int ni = 0; ni <= p.N; ++ni) {
                        m.row_add_int(x.p1(j, k, i, ni, h), 1);
                        m.row_add_int(x.p1c(j, k, i, ni, h), 1);
                    }
                }
            }
            m.emit_eq_int(1);
        }
    }
}

/** UTLB: U(i,k) = sum over nt,h of p1(i,k,t,nt,h), one row per witness t. */
template <class T>
void p1_utlb(const MapqnParams<T>& p, const MapqnP1Index& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int k = 0; k < p.K[i]; ++k) {
            for (int t = 0; t < p.M; ++t) {
                m.row_add_int(x.U(i, k), 1);
                for (int nt = 0; nt <= p.N; ++nt)
                    for (int h = 0; h < p.K[t]; ++h) m.row_add_int(x.p1(i, k, t, nt, h), -1);
                m.emit_eq_int(0);
            }
        }
    }
}

/** UTLC: IT(i,k) = sum over nt,h of p1c(i,k,t,nt,h), one row per witness t. */
template <class T>
void p1_utlc(const MapqnParams<T>& p, const MapqnP1Index& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int k = 0; k < p.K[i]; ++k) {
            for (int t = 0; t < p.M; ++t) {
                m.row_add_int(x.IT(i, k), 1);
                for (int nt = 0; nt <= p.N; ++nt)
                    for (int h = 0; h < p.K[t]; ++h) m.row_add_int(x.p1c(i, k, t, nt, h), -1);
                m.emit_eq_int(0);
            }
        }
    }
}

/** QLEN: Q(i,k) = sum over ni of ni * p1(i,k,i,ni,k). */
template <class T>
void p1_qlen(const MapqnParams<T>& p, const MapqnP1Index& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int k = 0; k < p.K[i]; ++k) {
            m.row_add_int(x.Q(i, k), 1);
            for (int ni = 0; ni <= p.N; ++ni) m.row_add_int(x.p1(i, k, i, ni, k), -ni);
            m.emit_eq_int(0);
        }
    }
}

/**
 * CLEN: C(j,k,i) = sum over ni,h of ni * p1(j,k,i,ni,h).
 *
 * Without it CEQU defines C only on the diagonal and every upper bound on C for
 * i != j is vacuous.
 */
template <class T>
void p1_clen(const MapqnParams<T>& p, const MapqnP1Index& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int k = 0; k < p.K[j]; ++k) {
            for (int i = 0; i < p.M; ++i) {
                m.row_add_int(x.C(j, k, i), 1);
                for (int ni = 0; ni <= p.N; ++ni)
                    for (int h = 0; h < p.K[i]; ++h) m.row_add_int(x.p1(j, k, i, ni, h), -ni);
                m.emit_eq_int(0);
            }
        }
    }
}

/** ONE: each station is busy in some phase or idle, with probability one. */
template <class T>
void p1_one(const MapqnParams<T>& p, const MapqnP1Index& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int k = 0; k < p.K[j]; ++k) {
            m.row_add_int(x.U(j, k), 1);
            m.row_add_int(x.IT(j, k), 1);
        }
        m.emit_eq_int(1);
    }
}

/** POPC: the mean lengths carry the whole closed population. */
template <class T>
void p1_popc(const MapqnParams<T>& p, const MapqnP1Index& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i)
        for (int k = 0; k < p.K[i]; ++k) m.row_add_int(x.Q(i, k), 1);
    m.emit_eq_int(p.N);
}

/**
 * SRVB: phase balance of the marginal, sum{j,h} q(i,j,k,h) U(i,k) equals
 * sum{j,h} q(i,j,h,k) U(i,h) (AMPL THM1).
 *
 * This is the only family that reads the transition rates into the U variables;
 * without it the phase split of each utilization is free. A single-phase
 * station gives an identically zero row and the references skip it rather than
 * emitting one, so keep the skip and the row counts match.
 */
template <class T>
void p1_srvb(const MapqnParams<T>& p, const MapqnP1Index& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        if (p.K[i] < 2) continue;
        for (int k = 0; k < p.K[i]; ++k) {
            for (int j = 0; j < p.M; ++j) {
                for (int h = 0; h < p.K[i]; ++h) {
                    m.row_add(x.U(i, k), qr_rate(p, i, j, k, h));
                    m.row_add(x.U(i, h), T(-qr_rate(p, i, j, h, k)));
                }
            }
            m.emit_eq_int(0);
        }
    }
}

/** UUB1: a station is busy in at most one phase at a time. */
template <class T>
void p1_uub1(const MapqnParams<T>& p, const MapqnP1Index& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int k = 0; k < p.K[i]; ++k) m.row_add_int(x.U(i, k), 1);
        m.emit_le_int(1);
    }
}

/** QUB1: Q(j,k) <= N U(j,k). An idle phase holds no jobs. */
template <class T>
void p1_qub1(const MapqnParams<T>& p, const MapqnP1Index& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int k = 0; k < p.K[j]; ++k) {
            m.row_add_int(x.Q(j, k), 1);
            m.row_add_int(x.U(j, k), -p.N);
            m.emit_le_int(0);
        }
    }
}

/** Shared argument validation for both p1-level entry points. */
template <class T>
void p1_check_objective(const MapqnParams<T>& p, int objective_queue, int objective_phase) {
    if (p.N < 1) throw InputError("mapqn: N must be at least 1");
    if (objective_queue < 0 || objective_queue >= p.M)
        throw InputError("mapqn: objective_queue out of range");
    if (objective_phase < 0 || objective_phase >= p.K[objective_queue])
        throw InputError("mapqn: objective_phase out of range");
}

}  // namespace detail

}  // namespace mapqn
}  // namespace line

#endif  // LINE_API_MAPQN_MAPQN_P1_COMMON_H
