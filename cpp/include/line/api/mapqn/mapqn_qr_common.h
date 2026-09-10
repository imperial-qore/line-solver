/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAPQN_MAPQN_QR_COMMON_H
#define LINE_API_MAPQN_MAPQN_QR_COMMON_H

/**
 * The constraint families of the quadratic-reduction (QR) polytope.
 *
 * mapqn_bnd_qr_ld and mapqn_bnd_qr_delay are the same polytope up to one
 * family: the delay model adds XZ, the think-time balance between the delay
 * station M and the throughput at queue 1. MATLAB keeps two nearly identical
 * 800-line files (mapqn_bnd_qr_ld.m and mapqn_bnd_qr_delay.m differ only in
 * XZ, the order in which PC2 is emitted, and the shape of the returned
 * marginals); the JAR duplicates them again. Here each family is emitted by
 * one function and each entry point calls the families it needs, so a fix to
 * a family cannot land in one entry point and miss the other.
 *
 * Emitting the families BY NAME is deliberate. The recurring failure mode in
 * this domain is an inventory gap: a file whose every subscript is internally
 * consistent but which omits a whole family, producing a vacuous [0,1] bound
 * that reads as "loose but valid" (see the mapqn notes in _kb/03-api-layer.md).
 * Naming the families makes the inventory diffable against the reference.
 *
 * Sign convention: every family is assembled as a single accumulated row and
 * emitted as `row = 0` or `row >= 0`, with the reference's right-hand terms
 * carried across with a negative sign rather than dropped. QBAL in particular
 * is LHS1 + LHS2 = RHS1 + RHS2; dropping the RHS terms would force-zero the
 * variables they carry and silently tighten the polytope.
 *
 * Variable bounds. The reference sets lb = 0, ub = 1 on every variable and
 * ub = 0 on the states the ZERO families exclude. Those are passed to the
 * solver as bounds, not as rows: line::lp::LpModel takes explicit per-variable
 * bounds, substitutes out any variable with lb == ub (which is every ZERO
 * state, the bulk of the model) and materializes a row only for a finite upper
 * bound that is not also a lower bound. The JAR has to add ~2 nVars explicit
 * rows instead, because Apache Commons SimplexSolver does not box variables
 * and the maximization is otherwise unbounded.
 */

#include <cstddef>
#include <vector>

#include "line/api/mapqn/mapqn_params.h"
#include "line/num/number.h"
#include "line/util/simplex.h"

namespace line {
namespace mapqn {

/**
 * ZERO1/2/3: states that carry no probability mass, imposed as ub = 0.
 *   ZERO1  i == j, nj == ni, h != k   (one queue cannot be in two phases)
 *   ZERO2  i == j, nj != ni           (one queue cannot hold two populations)
 *   ZERO3  i != j, nj + ni > N        (more jobs than the network holds)
 * Returns the indicator so SYMMETRY can skip pairs that are both zeroed, as
 * the reference does.
 */
template <class T>
std::vector<char> qr_zero_bounds(const MapqnParams<T>& p, const P2Index& idx, lp::LpModel<T>& m) {
    const int M = p.M, N = p.N;
    std::vector<char> is_zero(idx.num_vars(), 0);
    for (int j = 0; j < M; ++j)
        for (int nj = 0; nj <= N; ++nj)
            for (int kj = 0; kj < p.K[j]; ++kj)
                for (int i = 0; i < M; ++i)
                    for (int ni = 0; ni <= N; ++ni)
                        for (int hi = 0; hi < p.K[i]; ++hi) {
                            const bool z = (i == j && nj == ni && hi != kj) || (i == j && nj != ni) ||
                                           (i != j && nj + ni > N);
                            if (z) {
                                const std::size_t v = idx(j, nj, kj, i, ni, hi);
                                is_zero[v] = 1;
                                m.set_upper(v, T());
                            }
                        }
    return is_zero;
}

/** ONE: sum over (nj,k) of p2(j,nj,k,j,nj,k) = 1, per queue j. */
template <class T>
void qr_one(const MapqnParams<T>& p, const P2Index& idx, lp::LpModel<T>& m) {
    const T one = num_traits<T>::from_int(1);
    for (int j = 0; j < p.M; ++j) {
        for (int nj = 0; nj <= p.N; ++nj)
            for (int kj = 0; kj < p.K[j]; ++kj) m.row_add(idx(j, nj, kj, j, nj, kj), one);
        m.emit_eq(one);
    }
}

/** SYMMETRY: p2(i,ni,h,j,nj,k) = p2(j,nj,k,i,ni,h), emitted once per pair. */
template <class T>
void qr_symmetry(const MapqnParams<T>& p, const P2Index& idx, lp::LpModel<T>& m,
                 const std::vector<char>& is_zero) {
    const int M = p.M, N = p.N;
    const T one = num_traits<T>::from_int(1);
    const T mone = num_traits<T>::from_int(-1);
    for (int j = 0; j < M; ++j)
        for (int nj = 0; nj <= N; ++nj)
            for (int kj = 0; kj < p.K[j]; ++kj)
                for (int i = j + 1; i < M; ++i)  // i > j only, one ordering
                    for (int ni = 0; ni <= N; ++ni) {
                        if (i != j && nj + ni > N) continue;
                        for (int hi = 0; hi < p.K[i]; ++hi) {
                            const std::size_t a = idx(j, nj, kj, i, ni, hi);
                            const std::size_t b = idx(i, ni, hi, j, nj, kj);
                            if (is_zero[a] && is_zero[b]) continue;
                            if (a == b) continue;
                            m.row_add(a, one);
                            m.row_add(b, mone);
                            m.emit_eq(T());
                        }
                    }
}

/**
 * MARGINALS: p2(j,nj,k,j,nj,k) = sum over (ni <= N-nj, h) of p2(j,nj,k,i,ni,h)
 * for every i != j. The diagonal entry is the marginal of queue j, so the
 * joint over the pair (j,i) must sum back to it.
 */
template <class T>
void qr_marginals(const MapqnParams<T>& p, const P2Index& idx, lp::LpModel<T>& m) {
    const int M = p.M, N = p.N;
    const T one = num_traits<T>::from_int(1);
    const T mone = num_traits<T>::from_int(-1);
    for (int j = 0; j < M; ++j)
        for (int kj = 0; kj < p.K[j]; ++kj)
            for (int nj = 0; nj <= N; ++nj)
                for (int i = 0; i < M; ++i) {
                    if (i == j) continue;
                    m.row_add(idx(j, nj, kj, j, nj, kj), one);
                    for (int ni = 0; ni <= N - nj; ++ni)
                        for (int hi = 0; hi < p.K[i]; ++hi) m.row_add(idx(j, nj, kj, i, ni, hi), mone);
                    m.emit_eq(T());
                }
}

/**
 * THM1 (Little's law in probability form): for each (j,k),
 *   sum_{i,nj>=1,ni>=1,h} ni p2(j,nj,k,i,ni,h) = N sum_{nj>=1} p2(j,nj,k,j,nj,k).
 */
template <class T>
void qr_thm1(const MapqnParams<T>& p, const P2Index& idx, lp::LpModel<T>& m) {
    const int M = p.M, N = p.N;
    for (int j = 0; j < M; ++j)
        for (int kj = 0; kj < p.K[j]; ++kj) {
            for (int i = 0; i < M; ++i)
                for (int nj = 1; nj <= N; ++nj)
                    for (int ni = 1; ni <= N; ++ni)
                        for (int hi = 0; hi < p.K[i]; ++hi)
                            m.row_add_int(idx(j, nj, kj, i, ni, hi), ni);
            for (int nj = 1; nj <= N; ++nj) m.row_add_int(idx(j, nj, kj, j, nj, kj), -N);
            m.emit_eq(T());
        }
}

/** THM1c: the nj = 0 companion of THM1, conditioning on queue j being empty. */
template <class T>
void qr_thm1c(const MapqnParams<T>& p, const P2Index& idx, lp::LpModel<T>& m) {
    const int M = p.M, N = p.N;
    for (int j = 0; j < M; ++j)
        for (int kj = 0; kj < p.K[j]; ++kj) {
            for (int i = 0; i < M; ++i)
                for (int ni = 1; ni <= N; ++ni)
                    for (int hi = 0; hi < p.K[i]; ++hi) m.row_add_int(idx(j, 0, kj, i, ni, hi), ni);
            m.row_add_int(idx(j, 0, kj, j, 0, kj), -N);
            m.emit_eq(T());
        }
}

/** PC2 (second moment): sum_{i,j,ni>=1,nj>=1,h,k} nj ni p2(j,nj,k,i,ni,h) = N^2. */
template <class T>
void qr_pc2(const MapqnParams<T>& p, const P2Index& idx, lp::LpModel<T>& m) {
    const int M = p.M, N = p.N;
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j)
            for (int ni = 1; ni <= N; ++ni)
                for (int nj = 1; nj <= N; ++nj)
                    for (int hi = 0; hi < p.K[i]; ++hi)
                        for (int kj = 0; kj < p.K[j]; ++kj)
                            m.row_add_int(idx(j, nj, kj, i, ni, hi), static_cast<long>(nj) * ni);
    m.emit_eq(num_traits<T>::from_int(static_cast<long>(N) * N));
}

/**
 * XZ (delay model only): the think-time balance
 *   sum_{ni>=1,k} ni p2(M,ni,k,M,ni,k) = (Z/D1) sum_{k,nj>=1} p2(1,nj,k,1,nj,k),
 * i.e. the mean population at the delay station equals Z times the throughput
 * of queue 1, whose service demand is D1.
 */
template <class T>
void qr_xz(const MapqnParams<T>& p, const P2Index& idx, lp::LpModel<T>& m) {
    const int M = p.M, N = p.N;
    const int last = M - 1;
    for (int ni = 1; ni <= N; ++ni)
        for (int kM = 0; kM < p.K[last]; ++kM) m.row_add_int(idx(last, ni, kM, last, ni, kM), ni);
    if (p.D1 == T()) throw InputError("mapqn_bnd_qr_delay: D1 must be nonzero");
    const T ratio = p.Z / p.D1;
    const T mratio = -ratio;
    for (int kj = 0; kj < p.K[0]; ++kj)
        for (int nj = 1; nj <= N; ++nj) m.row_add(idx(0, nj, kj, 0, nj, kj), mratio);
    m.emit_eq(T());
}

/**
 * THM2 (phase balance): for each (i,k) the total rate out of phase k at queue
 * i equals the total rate into it, summed over populations ni >= 1.
 */
template <class T>
void qr_thm2(const MapqnParams<T>& p, const P2Index& idx, lp::LpModel<T>& m) {
    const int M = p.M, N = p.N;
    for (int i = 0; i < M; ++i)
        for (int ki = 0; ki < p.K[i]; ++ki) {
            for (int j = 0; j < M; ++j)
                for (int hi = 0; hi < p.K[i]; ++hi) {
                    if (hi == ki && j == i) continue;
                    for (int ni = 1; ni <= N; ++ni) {
                        const T q_out = mapqn_q(p, i, j, ki, hi, ni);
                        const T q_in = mapqn_q(p, i, j, hi, ki, ni);
                        m.row_add(idx(i, ni, ki, i, ni, ki), q_out);
                        const T mq_in = -q_in;
                        m.row_add(idx(i, ni, hi, i, ni, hi), mq_in);
                    }
                }
            m.emit_eq(T());
        }
}

/**
 * THM3a (population flow balance, 1 <= ni <= N-1): the rate at which queue i
 * is entered while holding ni jobs equals the rate at which it is left while
 * holding ni+1.
 */
template <class T>
void qr_thm3a(const MapqnParams<T>& p, const P2Index& idx, lp::LpModel<T>& m) {
    const int M = p.M, N = p.N;
    for (int i = 0; i < M; ++i)
        for (int ni = 1; ni <= N - 1; ++ni) {
            for (int j = 0; j < M; ++j) {
                if (j == i) continue;
                for (int kj = 0; kj < p.K[j]; ++kj)
                    for (int hj = 0; hj < p.K[j]; ++hj)
                        for (int u = 0; u < p.K[i]; ++u)
                            for (int nj = 1; nj <= N - ni; ++nj)
                                m.row_add(idx(j, nj, kj, i, ni, u), mapqn_q(p, j, i, kj, hj, nj));
            }
            for (int j = 0; j < M; ++j) {
                if (j == i) continue;
                for (int ki = 0; ki < p.K[i]; ++ki)
                    for (int hi = 0; hi < p.K[i]; ++hi) {
                        const T qv = mapqn_q(p, i, j, ki, hi, ni + 1);
                        const T mqv = -qv;
                        m.row_add(idx(i, ni + 1, ki, i, ni + 1, ki), mqv);
                    }
            }
            m.emit_eq(T());
        }
}

/** THM3b: the ni = 0 boundary case of THM3a, resolved per arrival phase u. */
template <class T>
void qr_thm3b(const MapqnParams<T>& p, const P2Index& idx, lp::LpModel<T>& m) {
    const int M = p.M, N = p.N;
    for (int i = 0; i < M; ++i)
        for (int u = 0; u < p.K[i]; ++u) {
            for (int j = 0; j < M; ++j) {
                if (j == i) continue;
                for (int kj = 0; kj < p.K[j]; ++kj)
                    for (int hj = 0; hj < p.K[j]; ++hj)
                        for (int nj = 1; nj <= N; ++nj)
                            m.row_add(idx(j, nj, kj, i, 0, u), mapqn_q(p, j, i, kj, hj, nj));
            }
            for (int j = 0; j < M; ++j) {
                if (j == i) continue;
                for (int ki = 0; ki < p.K[i]; ++ki) {
                    const T qv = mapqn_q(p, i, j, ki, u, 1);
                    const T mqv = -qv;
                    m.row_add(idx(i, 1, ki, i, 1, ki), mqv);
                }
            }
            m.emit_eq(T());
        }
}

/**
 * QBAL (queue balance): LHS1 + LHS2 = RHS1 + RHS2 for each (i,k).
 *
 * This is the family the reference warns about. The right-hand side carries
 * two distinct blocks -- the arrival flow into queue i (RHS1, itself in two
 * pieces, the ni = 0 term and the ni >= 1 term) and the population-weighted
 * phase inflow (RHS2). Dropping either would leave the variables they touch
 * with no other constraint mentioning them and force them to zero.
 */
template <class T>
void qr_qbal(const MapqnParams<T>& p, const P2Index& idx, lp::LpModel<T>& m) {
    const int M = p.M, N = p.N;
    for (int i = 0; i < M; ++i)
        for (int ki = 0; ki < p.K[i]; ++ki) {
            // LHS1: sum_{h!=k, j, ni>=1} q(i,j,k,h,ni) * ni * p2(i,ni,k,i,ni,k)
            for (int hi = 0; hi < p.K[i]; ++hi) {
                if (hi == ki) continue;
                for (int j = 0; j < M; ++j)
                    for (int ni = 1; ni <= N; ++ni) {
                        const T qv = mapqn_q(p, i, j, ki, hi, ni);
                        const T w = qv * num_traits<T>::from_int(ni);
                        m.row_add(idx(i, ni, ki, i, ni, ki), w);
                    }
            }
            // LHS2: sum_{j!=i, h, ni>=1} q(i,j,h,k,ni) * p2(i,ni,h,i,ni,h)
            for (int j = 0; j < M; ++j) {
                if (j == i) continue;
                for (int hi = 0; hi < p.K[i]; ++hi)
                    for (int ni = 1; ni <= N; ++ni)
                        m.row_add(idx(i, ni, hi, i, ni, hi), mapqn_q(p, i, j, hi, ki, ni));
            }
            // -RHS1a: arrivals finding queue i empty
            for (int j = 0; j < M; ++j) {
                if (j == i) continue;
                for (int u = 0; u < p.K[j]; ++u)
                    for (int w = 0; w < p.K[j]; ++w)
                        for (int nj = 1; nj <= N; ++nj) {
                            const T qv = mapqn_q(p, j, i, u, w, nj);
                            const T mqv = -qv;
                            m.row_add(idx(j, nj, u, i, 0, ki), mqv);
                        }
            }
            // -RHS1b: arrivals finding queue i busy
            for (int j = 0; j < M; ++j) {
                if (j == i) continue;
                for (int u = 0; u < p.K[j]; ++u)
                    for (int w = 0; w < p.K[j]; ++w)
                        for (int nj = 1; nj <= N; ++nj) {
                            const T qv = mapqn_q(p, j, i, u, w, nj);
                            const T mqv = -qv;
                            for (int ni = 1; ni <= N; ++ni) m.row_add(idx(i, ni, ki, j, nj, u), mqv);
                        }
            }
            // -RHS2: population-weighted phase inflow
            for (int hi = 0; hi < p.K[i]; ++hi) {
                if (hi == ki) continue;
                for (int j = 0; j < M; ++j)
                    for (int ni = 1; ni <= N; ++ni) {
                        const T qv = mapqn_q(p, i, j, hi, ki, ni);
                        const T w = qv * num_traits<T>::from_int(ni);
                        const T mw = -w;
                        m.row_add(idx(i, ni, hi, i, ni, hi), mw);
                    }
            }
            m.emit_eq(T());
        }
}

/**
 * COR1a: the order-1 correlation cut, for each (i, kstar, ni = 0..N-2).
 * Blocks A..H follow the reference letter for letter; A and B are the arrival
 * terms, C..H the departure and phase-change terms at populations ni+1 and
 * ni+2.
 */
template <class T>
void qr_cor1a(const MapqnParams<T>& p, const P2Index& idx, lp::LpModel<T>& m) {
    const int M = p.M, N = p.N;
    for (int i = 0; i < M; ++i)
        for (int kstar = 0; kstar < p.K[i]; ++kstar)
            for (int nic = 0; nic <= N - 2; ++nic) {
                // A
                for (int j = 0; j < M; ++j) {
                    if (j == i) continue;
                    for (int kj = 0; kj < p.K[j]; ++kj)
                        for (int hj = 0; hj < p.K[j]; ++hj)
                            for (int u = 0; u < p.K[i]; ++u) {
                                if (u == kstar) continue;
                                for (int nj = 1; nj <= N - nic; ++nj)
                                    m.row_add(idx(j, nj, kj, i, nic, u), mapqn_q(p, j, i, kj, hj, nj));
                            }
                }
                // B
                for (int j = 0; j < M; ++j) {
                    if (j == i) continue;
                    for (int kj = 0; kj < p.K[j]; ++kj)
                        for (int hj = 0; hj < p.K[j]; ++hj)
                            for (int nj = 1; nj <= N - nic; ++nj)
                                m.row_add(idx(j, nj, kj, i, nic + 1, kstar),
                                          mapqn_q(p, j, i, kj, hj, nj));
                }
                // C
                for (int k2 = 0; k2 < p.K[i]; ++k2) {
                    if (k2 == kstar) continue;
                    m.row_add(idx(i, nic + 1, kstar, i, nic + 1, kstar),
                              mapqn_q(p, i, i, kstar, k2, nic + 1));
                }
                // -D
                for (int j = 0; j < M; ++j) {
                    if (j == i) continue;
                    for (int k2 = 0; k2 < p.K[i]; ++k2) {
                        if (k2 == kstar) continue;
                        const T qv = mapqn_q(p, i, j, k2, k2, nic + 1);
                        const T mqv = -qv;
                        m.row_add(idx(i, nic + 1, k2, i, nic + 1, k2), mqv);
                    }
                }
                // -E
                for (int j = 0; j < M; ++j) {
                    if (j == i) continue;
                    for (int k2 = 0; k2 < p.K[i]; ++k2) {
                        if (k2 == kstar) continue;
                        for (int h2 = 0; h2 < p.K[i]; ++h2) {
                            if (h2 == k2) continue;
                            const T qv = mapqn_q(p, i, j, k2, h2, nic + 1);
                            const T mqv = -qv;
                            m.row_add(idx(i, nic + 1, k2, i, nic + 1, k2), mqv);
                        }
                    }
                }
                // -F
                for (int j = 0; j < M; ++j) {
                    if (j == i) continue;
                    for (int k2 = 0; k2 < p.K[i]; ++k2) {
                        if (k2 == kstar) continue;
                        const T qv = mapqn_q(p, i, j, k2, kstar, nic + 2);
                        const T mqv = -qv;
                        m.row_add(idx(i, nic + 2, k2, i, nic + 2, k2), mqv);
                    }
                }
                // -G
                for (int j = 0; j < M; ++j) {
                    if (j == i) continue;
                    const T qv = mapqn_q(p, i, j, kstar, kstar, nic + 2);
                    const T mqv = -qv;
                    m.row_add(idx(i, nic + 2, kstar, i, nic + 2, kstar), mqv);
                }
                // -H
                for (int k2 = 0; k2 < p.K[i]; ++k2) {
                    if (k2 == kstar) continue;
                    const T qv = mapqn_q(p, i, i, k2, kstar, nic + 1);
                    const T mqv = -qv;
                    m.row_add(idx(i, nic + 1, k2, i, nic + 1, k2), mqv);
                }
                m.emit_eq(T());
            }
}

/** COR1b: the ni = N-1 boundary of COR1a (blocks A', C', D', E', H'). */
template <class T>
void qr_cor1b(const MapqnParams<T>& p, const P2Index& idx, lp::LpModel<T>& m) {
    const int M = p.M, N = p.N;
    for (int i = 0; i < M; ++i)
        for (int kstar = 0; kstar < p.K[i]; ++kstar) {
            // A'
            for (int j = 0; j < M; ++j) {
                if (j == i) continue;
                for (int kj = 0; kj < p.K[j]; ++kj)
                    for (int hj = 0; hj < p.K[j]; ++hj)
                        for (int u = 0; u < p.K[i]; ++u) {
                            if (u == kstar) continue;
                            m.row_add(idx(j, 1, kj, i, N - 1, u), mapqn_q(p, j, i, kj, hj, 1));
                        }
            }
            // C'
            for (int k2 = 0; k2 < p.K[i]; ++k2) {
                if (k2 == kstar) continue;
                m.row_add(idx(i, N, kstar, i, N, kstar), mapqn_q(p, i, i, kstar, k2, N));
            }
            // -D'
            for (int j = 0; j < M; ++j) {
                if (j == i) continue;
                for (int k2 = 0; k2 < p.K[i]; ++k2) {
                    if (k2 == kstar) continue;
                    const T qv = mapqn_q(p, i, j, k2, k2, N);
                    const T mqv = -qv;
                    m.row_add(idx(i, N, k2, i, N, k2), mqv);
                }
            }
            // -E'
            for (int j = 0; j < M; ++j) {
                if (j == i) continue;
                for (int k2 = 0; k2 < p.K[i]; ++k2) {
                    if (k2 == kstar) continue;
                    for (int h2 = 0; h2 < p.K[i]; ++h2) {
                        if (h2 == k2) continue;
                        const T qv = mapqn_q(p, i, j, k2, h2, N);
                        const T mqv = -qv;
                        m.row_add(idx(i, N, k2, i, N, k2), mqv);
                    }
                }
            }
            // -H'
            for (int k2 = 0; k2 < p.K[i]; ++k2) {
                if (k2 == kstar) continue;
                const T qv = mapqn_q(p, i, i, k2, kstar, N);
                const T mqv = -qv;
                m.row_add(idx(i, N, k2, i, N, k2), mqv);
            }
            m.emit_eq(T());
        }
}

/**
 * THM4 (QMIN): for each (j,k,i),
 *   sum_{t,h,nj,nt} nt p2(j,nj,k,t,nt,h) >= N sum_{h,nj,ni} p2(j,nj,k,i,ni,h),
 * the only inequality family. The reference passes it to linprog as
 * -row * x <= 0; here it is emitted directly as row >= 0.
 */
template <class T>
void qr_thm4(const MapqnParams<T>& p, const P2Index& idx, lp::LpModel<T>& m) {
    const int M = p.M, N = p.N;
    for (int j = 0; j < M; ++j)
        for (int kj = 0; kj < p.K[j]; ++kj)
            for (int i = 0; i < M; ++i) {
                for (int t = 0; t < M; ++t)
                    for (int ht = 0; ht < p.K[t]; ++ht)
                        for (int nj = 0; nj <= N; ++nj)
                            for (int nt = 0; nt <= N; ++nt)
                                m.row_add_int(idx(j, nj, kj, t, nt, ht), nt);
                for (int hi = 0; hi < p.K[i]; ++hi)
                    for (int nj = 0; nj <= N; ++nj)
                        for (int ni = 0; ni <= N; ++ni)
                            m.row_add_int(idx(j, nj, kj, i, ni, hi), -N);
                m.emit_ge(T());
            }
}

namespace detail {

/** Common tail: solve, then unpack the objective and the diagonal marginals. */
template <class T>
MapqnQrResult<T> qr_finish(const MapqnParams<T>& p, const P2Index& idx, lp::LpModel<T>& m,
                           int objective_queue, int objective_phase, int objective_n,
                           MapqnSense sense) {
    m.set_cost(idx(objective_queue, objective_n, objective_phase, objective_queue, objective_n,
                   objective_phase),
               num_traits<T>::from_int(1));
    m.set_maximize(sense == MapqnSense::Max);

    MapqnQrResult<T> res;
    res.num_vars = m.num_vars();
    res.num_rows = m.num_rows();
    const lp::LpSolution<T> s = lp::simplex_solve(m);
    res.status = lp::lp_status_name(s.status);
    res.iterations = s.iterations;
    res.ok = s.ok();
    if (!res.ok) return res;
    res.objective = s.objective;
    res.x = s.x;
    res.p2marginals.resize(static_cast<std::size_t>(p.M));
    for (int j = 0; j < p.M; ++j) {
        res.p2marginals[j] = Matrix<T>(static_cast<std::size_t>(p.N + 1),
                                       static_cast<std::size_t>(p.K[j]), T());
        for (int nj = 0; nj <= p.N; ++nj)
            for (int kj = 0; kj < p.K[j]; ++kj)
                res.p2marginals[j](static_cast<std::size_t>(nj), static_cast<std::size_t>(kj)) =
                    s.x[idx(j, nj, kj, j, nj, kj)];
    }
    return res;
}

/** Shared argument validation for both entry points. */
template <class T>
void qr_check_objective(const MapqnParams<T>& p, int objective_queue, int objective_phase,
                        int objective_n) {
    if (objective_queue < 0 || objective_queue >= p.M)
        throw InputError("mapqn: objective_queue out of range");
    if (objective_phase < 0 || objective_phase >= p.K[objective_queue])
        throw InputError("mapqn: objective_phase out of range");
    if (objective_n < 0 || objective_n > p.N) throw InputError("mapqn: objective_n out of range");
}

}  // namespace detail

}  // namespace mapqn
}  // namespace line

#endif  // LINE_API_MAPQN_MAPQN_QR_COMMON_H
