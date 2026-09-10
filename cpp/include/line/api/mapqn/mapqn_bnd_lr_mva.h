/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAPQN_MAPQN_BND_LR_MVA_H
#define LINE_API_MAPQN_MAPQN_BND_LR_MVA_H

/**
 * MVA-shaped linear-reduction bound for a closed network of M - 1 exponential
 * queues and ONE MAP queue.
 *
 * Templated port of matlab/lib/qrf/mapqn_bnd_lr_mva.m (ground truth),
 * cross-checked against python/line_solver/api/mapqn/bnd_lr_mva.py.
 *
 * A DIFFERENT MODEL, NOT A VARIANT. mapqn_bnd_lr relaxes onto the marginal
 * probabilities p1(j,k,i,ni,h) and reads utilizations out of them. This one
 * never introduces a probability variable at all: it works directly on the
 * mean-value quantities UN(i,k), QN(i,k) and the conditional length B(j,k,i),
 * in the manner of an MVA recursion turned into a relaxation. The consequence
 * is a model of 2 M K + M^2 K columns instead of one quadratic in (N+1) sum_i
 * K(i), so it is by far the cheapest bound in the family and the only one that
 * stays small as the population grows. Nothing here is shared with
 * mapqn_p1_common.h, and it deliberately does not use MapqnParams: the network
 * shape is different (see LrMvaParams).
 *
 * NETWORK SHAPE. Queues 1..M-1 are exponential with scalar rates muM(i), and
 * queue M is the MAP, whose K levels are the phase process. That asymmetry is
 * carried entirely by q(i,j,k,h) below and by nothing else, which is why the
 * families read uniformly over i even though the stations are not alike.
 *
 * THE FAILURE MODE OF THIS FAMILY IS A VACUOUS BOUND, NOT A CRASH. FLOW, UBAL,
 * QBAL, MCC and MCC2 are the families that read the rates; the rest are
 * structural. Drop the rate-bearing ones and UN = 0 stays feasible while every
 * subscript remains internally consistent. Each family is emitted by a function
 * named after the reference's own numbered section, so the inventory is
 * diffable; do not inline them.
 *
 * ARITHMETIC. Assembly is +, -, * on the model data and lp::simplex_solve uses
 * Bland's rule with no tolerance, so at T = line::Rational the returned bound
 * is the EXACT optimum of the exact polytope. The MATLAB reference defaults to
 * linprog's 'interior-point-legacy' here, which on these badly scaled instances
 * is markedly more accurate than plain 'interior-point'.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mapqn/mapqn_params.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/simplex.h"

namespace line {
namespace mapqn {

/**
 * Parameters of the MVA-shaped LR bound, mirroring the reference's `params`.
 *
 * Queues are 0-based here, so the MAP queue is index M - 1 and muM carries the
 * M - 1 exponential rates of queues 0..M-2. K is a scalar: it is the level
 * count of the single MAP, not a per-queue vector.
 */
template <class T>
struct LrMvaParams {
    int M = 0;         ///< number of queues, the last of which is the MAP
    int N = 0;         ///< closed population
    int K = 0;         ///< number of levels of the MAP queue
    std::vector<T> muM;  ///< (M-1) service rates of the exponential queues
    Matrix<T> muMAP;   ///< (K x K) completion rates of the MAP queue
    Matrix<T> v;       ///< (K x K) level-change rates of the MAP queue
    Matrix<T> r;       ///< (M x M) routing probabilities

    void validate() const {
        if (M <= 1) throw InputError("mapqn_bnd_lr_mva: M must be at least 2");
        if (N < 1) throw InputError("mapqn_bnd_lr_mva: N must be at least 1");
        if (K < 1) throw InputError("mapqn_bnd_lr_mva: K must be at least 1");
        if (static_cast<int>(muM.size()) != M - 1)
            throw InputError("mapqn_bnd_lr_mva: muM must have M-1 entries");
        const std::size_t k = static_cast<std::size_t>(K);
        if (muMAP.rows() != k || muMAP.cols() != k)
            throw InputError("mapqn_bnd_lr_mva: muMAP must be K x K");
        if (v.rows() != k || v.cols() != k)
            throw InputError("mapqn_bnd_lr_mva: v must be K x K");
        if (r.rows() != static_cast<std::size_t>(M) || r.cols() != static_cast<std::size_t>(M))
            throw InputError("mapqn_bnd_lr_mva: r must be M x M");
    }
};

/**
 * The variable family the objective is taken over.
 *
 * The paper states its bounds on the AGGREGATE over levels: U_i(N) = sum_k
 * U_i^k(N) is the utilization of station i, while U_i^k alone is its
 * utilization while the MAP sits in phase k. Pass objective_level = -1 for
 * that aggregate; optimizing the K terms separately and adding them is also a
 * bound but a strictly looser one, since the phases cannot all peak at once.
 */
enum class MapqnObjectiveVar { UN, QN };

/** Result of an MVA-shaped LR bound solve. */
template <class T>
struct MapqnBndLrMvaResult {
    bool ok = false;      ///< the LP reached an optimal vertex
    std::string status;   ///< textual LP status
    T objective = T();    ///< the bound on UN(objective_queue, objective_level)
    Matrix<T> UN;         ///< M x K utilizations
    Matrix<T> QN;         ///< M x K queue lengths
    std::vector<T> x;     ///< full solution vector, indexed by LrMvaIndex
    std::size_t num_vars = 0;
    std::size_t num_rows = 0;
    std::size_t iterations = 0;
};

/** Variable layout: UN, then QN, then B, in the reference's declaration order. */
struct LrMvaIndex {
    int M = 0, K = 0;
    std::size_t off_UN = 0, off_QN = 0, off_B = 0, total = 0;

    LrMvaIndex() {}
    LrMvaIndex(int m, int k) : M(m), K(k) {
        const std::size_t mk = static_cast<std::size_t>(M) * static_cast<std::size_t>(K);
        off_UN = 0;
        off_QN = mk;
        off_B = 2 * mk;
        total = off_B + mk * static_cast<std::size_t>(M);
    }
    std::size_t UN(int i, int k) const {
        return off_UN + static_cast<std::size_t>(i) * static_cast<std::size_t>(K) +
               static_cast<std::size_t>(k);
    }
    std::size_t QN(int i, int k) const {
        return off_QN + static_cast<std::size_t>(i) * static_cast<std::size_t>(K) +
               static_cast<std::size_t>(k);
    }
    std::size_t B(int j, int k, int i) const {
        return off_B +
               (static_cast<std::size_t>(j) * static_cast<std::size_t>(K) +
                static_cast<std::size_t>(k)) *
                   static_cast<std::size_t>(M) +
               static_cast<std::size_t>(i);
    }
    std::size_t num_vars() const { return total; }
};

namespace detail {

/**
 * q(i,j,k,h): rate at which queue i at level k routes a job to queue j and the
 * level becomes h.
 *
 * Port of the nested q of mapqn_bnd_lr_mva.m. The two branches are the whole of
 * the network's asymmetry:
 *   i < M-1  an exponential queue, which cannot change the level, so the rate
 *            is r(i,j) mu(i) on the diagonal k == h and zero off it;
 *   i == M-1 the MAP. Routing away (j < M-1) carries muMAP(k,h). Routing to
 *            itself carries v(k,h) + r(M-1,M-1) muMAP(k,h) and is zero at
 *            k == h, since a self-loop that changes nothing is not an event.
 */
template <class T>
T lr_mva_rate(const LrMvaParams<T>& p, int i, int j, int k, int h) {
    const std::size_t ki = static_cast<std::size_t>(k), hi = static_cast<std::size_t>(h);
    const int last = p.M - 1;
    if (i < last) {
        if (k != h) return T();
        return T(p.r(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) *
                 p.muM[static_cast<std::size_t>(i)]);
    }
    if (j < last)
        return T(p.r(static_cast<std::size_t>(last), static_cast<std::size_t>(j)) *
                 p.muMAP(ki, hi));
    if (k == h) return T();
    return T(p.v(ki, hi) + p.r(static_cast<std::size_t>(last), static_cast<std::size_t>(last)) *
                               p.muMAP(ki, hi));
}

/** Variable bounds: 0 <= UN <= 1, 0 <= QN <= N, 0 <= B <= N. */
template <class T>
void lr_mva_bounds(const LrMvaParams<T>& p, const LrMvaIndex& x, lp::LpModel<T>& m) {
    const T one = num_traits<T>::from_int(1);
    const T nn = num_traits<T>::from_int(p.N);
    for (int i = 0; i < p.M; ++i) {
        for (int k = 0; k < p.K; ++k) {
            m.set_bounds(x.UN(i, k), T(), one);
            m.set_bounds(x.QN(i, k), T(), nn);
            for (int j = 0; j < p.M; ++j) m.set_bounds(x.B(i, k, j), T(), nn);
        }
    }
}

/** 1. QNB: B(j,k,i) <= QN(i,k). A conditional length cannot exceed the mean. */
template <class T>
void lr_mva_qnb(const LrMvaParams<T>& p, const LrMvaIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int k = 0; k < p.K; ++k) {
            for (int j = 0; j < p.M; ++j) {
                m.row_add_int(x.QN(i, k), -1);
                m.row_add_int(x.B(j, k, i), 1);
                m.emit_le_int(0);
            }
        }
    }
}

/** 2. UMAX: a station is busy at one level at a time. */
template <class T>
void lr_mva_umax(const LrMvaParams<T>& p, const LrMvaIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int k = 0; k < p.K; ++k) m.row_add_int(x.UN(i, k), 1);
        m.emit_le_int(1);
    }
}

/** 3. POPCONSTR: the mean lengths carry the whole closed population. */
template <class T>
void lr_mva_popconstr(const LrMvaParams<T>& p, const LrMvaIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i)
        for (int k = 0; k < p.K; ++k) m.row_add_int(x.QN(i, k), 1);
    m.emit_eq_int(p.N);
}

/** 4. FLOW: rate into station i equals rate out of it. */
template <class T>
void lr_mva_flow(const LrMvaParams<T>& p, const LrMvaIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int k = 0; k < p.K; ++k) {
            for (int mm = 0; mm < p.K; ++mm) {
                for (int w = 0; w < p.M; ++w) {
                    m.row_add(x.UN(w, k), lr_mva_rate(p, w, i, k, mm));
                    m.row_add(x.UN(i, mm), T(-lr_mva_rate(p, i, w, mm, k)));
                }
            }
        }
        m.emit_eq_int(0);
    }
}

/** 5. UBAL: level balance of the MAP queue's utilization. */
template <class T>
void lr_mva_ubal(const LrMvaParams<T>& p, const LrMvaIndex& x, lp::LpModel<T>& m) {
    const int last = p.M - 1;
    for (int k = 0; k < p.K; ++k) {
        for (int h = 0; h < p.K; ++h) {
            if (h == k) continue;
            for (int w = 0; w < p.M; ++w) {
                m.row_add(x.UN(last, k), lr_mva_rate(p, last, w, k, h));
                m.row_add(x.UN(last, h), T(-lr_mva_rate(p, last, w, h, k)));
            }
        }
        m.emit_eq_int(0);
    }
}

/** 6. QBAL: level balance of the MAP queue's mean length. */
template <class T>
void lr_mva_qbal(const LrMvaParams<T>& p, const LrMvaIndex& x, lp::LpModel<T>& m) {
    const int last = p.M - 1;
    for (int k = 0; k < p.K; ++k) {
        for (int h = 0; h < p.K; ++h) {
            if (h == k) continue;
            for (int w = 0; w < p.M; ++w) m.row_add(x.QN(last, k), lr_mva_rate(p, last, w, k, h));
        }
        for (int mm = 0; mm < p.K; ++mm)
            for (int j = 0; j < last; ++j)
                m.row_add(x.UN(last, mm), lr_mva_rate(p, last, j, mm, k));
        for (int j = 0; j < last; ++j)
            m.row_add(x.UN(j, k), T(-lr_mva_rate(p, j, last, k, k)));
        for (int h = 0; h < p.K; ++h) {
            if (h == k) continue;
            for (int w = 0; w < p.M; ++w)
                m.row_add(x.QN(last, h), T(-lr_mva_rate(p, last, w, h, k)));
        }
        m.emit_eq_int(0);
    }
}

/** 7. MCC: the marginal-consistency cut, in its (N+1)-weighted form. */
template <class T>
void lr_mva_mcc(const LrMvaParams<T>& p, const LrMvaIndex& x, lp::LpModel<T>& m) {
    const T np1 = num_traits<T>::from_int(p.N + 1);
    for (int i = 0; i < p.M; ++i) {
        for (int k = 0; k < p.K; ++k) {
            for (int mm = 0; mm < p.K; ++mm) {
                for (int w = 0; w < p.M; ++w) {
                    if (w == i) continue;
                    m.row_add(x.QN(i, k), lr_mva_rate(p, i, w, k, mm));
                }
                for (int j = 0; j < p.M; ++j) {
                    if (j == i) continue;
                    const T q = lr_mva_rate(p, j, i, k, mm);
                    m.row_add(x.QN(j, k), q);
                    for (int wp = 0; wp < p.M; ++wp) {
                        if (wp == i || wp == j) continue;
                        m.row_add(x.B(j, k, wp), q);
                    }
                    m.row_add(x.UN(j, k), T(-(np1 * q)));
                }
            }
        }
        m.emit_eq_int(0);
    }
}

/** 8. MCC2: the same cut with the conditional length taken at i itself. */
template <class T>
void lr_mva_mcc2(const LrMvaParams<T>& p, const LrMvaIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int k = 0; k < p.K; ++k) {
            for (int mm = 0; mm < p.K; ++mm) {
                for (int w = 0; w < p.M; ++w) {
                    if (w == i) continue;
                    m.row_add(x.QN(i, k), lr_mva_rate(p, i, w, k, mm));
                }
                for (int j = 0; j < p.M; ++j) {
                    if (j == i) continue;
                    const T q = lr_mva_rate(p, j, i, k, mm);
                    m.row_add(x.B(j, k, i), T(-q));
                    m.row_add(x.UN(j, k), T(-q));
                }
            }
        }
        m.emit_eq_int(0);
    }
}

/** 9. QMAX: QN(w,k) <= N UN(w,k). An idle level holds no jobs. */
template <class T>
void lr_mva_qmax(const LrMvaParams<T>& p, const LrMvaIndex& x, lp::LpModel<T>& m) {
    for (int w = 0; w < p.M; ++w) {
        for (int k = 0; k < p.K; ++k) {
            m.row_add_int(x.QN(w, k), 1);
            m.row_add_int(x.UN(w, k), -p.N);
            m.emit_le_int(0);
        }
    }
}

/** 10. QMIN: N UN(j,k) <= sum_w QN(w,k), for every witness j. */
template <class T>
void lr_mva_qmin(const LrMvaParams<T>& p, const LrMvaIndex& x, lp::LpModel<T>& m) {
    for (int k = 0; k < p.K; ++k) {
        for (int j = 0; j < p.M; ++j) {
            for (int w = 0; w < p.M; ++w) m.row_add_int(x.QN(w, k), -1);
            m.row_add_int(x.UN(j, k), p.N);
            m.emit_le_int(0);
        }
    }
}

}  // namespace detail

/**
 * Bound UN or QN at (objective_queue, objective_level) over the MVA-shaped LR
 * polytope.
 *
 * @param p               network parameters; queues and levels are 0-based, and
 *                        the MAP queue is index M - 1
 * @param objective_queue queue index, 0..M-1
 * @param objective_level level index, 0..K-1, or -1 for the SUM over levels
 * @param sense           Max for an upper bound, Min for a lower bound
 * @param objective_var   UN (default) or QN, the variable family optimized over
 */
template <class T>
MapqnBndLrMvaResult<T> mapqn_bnd_lr_mva(const LrMvaParams<T>& p, int objective_queue,
                                        int objective_level,
                                        MapqnSense sense = MapqnSense::Max,
                                        MapqnObjectiveVar objective_var = MapqnObjectiveVar::UN) {
    p.validate();
    if (objective_queue < 0 || objective_queue >= p.M)
        throw InputError("mapqn_bnd_lr_mva: objective_queue out of range");
    if (objective_level < -1 || objective_level >= p.K)
        throw InputError("mapqn_bnd_lr_mva: objective_level out of range");

    const LrMvaIndex x(p.M, p.K);
    lp::LpModel<T> m(x.num_vars());
    detail::lr_mva_bounds(p, x, m);

    // Families, in the reference's numbered order. Named one per function so
    // the inventory is diffable against mapqn_bnd_lr_mva.m.
    detail::lr_mva_qnb(p, x, m);
    detail::lr_mva_umax(p, x, m);
    detail::lr_mva_popconstr(p, x, m);
    detail::lr_mva_flow(p, x, m);
    detail::lr_mva_ubal(p, x, m);
    detail::lr_mva_qbal(p, x, m);
    detail::lr_mva_mcc(p, x, m);
    detail::lr_mva_mcc2(p, x, m);
    detail::lr_mva_qmax(p, x, m);
    detail::lr_mva_qmin(p, x, m);

    // One level, or their sum when objective_level is -1.
    const int first_level = (objective_level < 0) ? 0 : objective_level;
    const int last_level = (objective_level < 0) ? p.K - 1 : objective_level;
    for (int k = first_level; k <= last_level; ++k) {
        m.set_cost(objective_var == MapqnObjectiveVar::UN ? x.UN(objective_queue, k)
                                                          : x.QN(objective_queue, k),
                   num_traits<T>::from_int(1));
    }
    m.set_maximize(sense == MapqnSense::Max);

    const lp::LpSolution<T> sol = lp::simplex_solve(m);

    MapqnBndLrMvaResult<T> out;
    out.status = lp::lp_status_name(sol.status);
    out.ok = sol.ok();
    out.objective = sol.objective;
    out.x = sol.x;
    out.num_vars = m.num_vars();
    out.num_rows = m.num_rows();
    out.iterations = sol.iterations;
    if (!out.ok) return out;

    out.UN = Matrix<T>(static_cast<std::size_t>(p.M), static_cast<std::size_t>(p.K));
    out.QN = Matrix<T>(static_cast<std::size_t>(p.M), static_cast<std::size_t>(p.K));
    for (int i = 0; i < p.M; ++i) {
        for (int k = 0; k < p.K; ++k) {
            const std::size_t ii = static_cast<std::size_t>(i), kk = static_cast<std::size_t>(k);
            out.UN(ii, kk) = sol.x[x.UN(i, k)];
            out.QN(ii, kk) = sol.x[x.QN(i, k)];
        }
    }
    return out;
}

}  // namespace mapqn
}  // namespace line

#endif  // LINE_API_MAPQN_MAPQN_BND_LR_MVA_H
