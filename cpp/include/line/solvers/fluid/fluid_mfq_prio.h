/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_MFQ_PRIO_H
#define LINE_SOLVERS_FLUID_FLUID_MFQ_PRIO_H

/**
 * The priority branch of the `mfq` method: a port of `solver_mfq_prio.m`.
 *
 * WHEN IT IS REACHED. `solver_fluid_analyzer.m` sends a single-queue open model
 * to `solver_mfq` when every class has the same priority and HERE when they do
 * not. The model is the same one `mfq` solves -- Source -> Queue -> Sink -- but
 * the server now drains the highest-priority fluid first, preemptively, so the
 * per-class levels no longer follow from one fluid-fluid queue.
 *
 * WHAT IS SOLVED. Each class arrives as a MAP. The per-class background chains
 * are SUPERPOSED into one joint chain by a Kronecker sum, and the joint chain
 * modulates a K x N rate matrix, one row per class. That pair (Qjoint, Rjoint)
 * with the constant drain rate d is exactly the input of the fluid priority
 * queue of G. Horvath, "Efficient analysis of the MMAP[K]/PH[K]/1 priority
 * queue", EJOR 246(1):128-139, 2015, already ported as `mam::mfq_prio_queue`.
 * The fluid level of a class is read as its queue length and the fluid sojourn
 * as its response time, the same reading `solver_mfq.m` makes.
 *
 * THE ROW ORDER MATTERS AND IS INVERTED. LINE's `classprio` is a RANK: the
 * SMALLER the value the HIGHER the priority. FluidPrioQueue takes the opposite
 * convention, the LAST row being the highest priority. The reference therefore
 * sorts the open classes by `classprio` DESCENDING before building the rows,
 * and this port keeps that sort stable so ties preserve class order.
 *
 * WHEN IT DECLINES. The reference falls back to the matrix fluid method, with a
 * warning, whenever the fluid priority model degenerates: a non-MAP arrival,
 * class-dependent or non-positive service, or a joint chain with a single state
 * (unmodulated arrivals, where the fluid level is identically zero and the
 * priority structure carries no information). This port reports the same
 * decision through `MfqPrioResult::fallback` and lets the caller re-dispatch,
 * rather than returning zeros that look like an answer.
 *
 * A FLUID LEVEL IS NOT A CUSTOMER COUNT. The level is zero unless the ARRIVAL
 * RATE EXCEEDS d in some background state; a class whose peak rate stays below
 * the drain rate reports QN = 0 exactly, which is correct for the fluid model
 * and is not a failure of the solve.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <string>
#include <vector>

#include "line/api/mam/mfq_prio_queue.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/fluid/fluid_mfq.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace fluid {

/** Per-class metrics of the priority queue, indexed by class. */
struct MfqPrioResult {
    std::vector<double> QN, RN, TN, UN;
    bool fallback = false;    ///< true when the reference would run the matrix method instead
    std::string reason;       ///< why, for the caller's warning
};

namespace detail {

/** Kronecker sum term: kron over j of (j == i ? A_i : I_{n_j}), in j order. */
inline Matrix<double> prio_kron_term(const std::vector<Matrix<double>>& A,
                                     const std::vector<std::size_t>& n, std::size_t i) {
    Matrix<double> term(1, 1, 1.0);
    for (std::size_t j = 0; j < A.size(); ++j) {
        if (j == i)
            term = mam::kron(term, A[i]);
        else
            term = mam::kron(term, line::eye<double>(n[j]));
    }
    return term;
}

}  // namespace detail

/**
 * Solve the single priority fluid queue of `sn`.
 *
 * @param top  the single-queue topology, as `mfq_is_single_queue` found it
 * @param tol  `options.tol`, passed to the Riccati iterations as their precision
 * @param sn the refreshed network struct
 */
template <class T>
MfqPrioResult fluid_mfq_prio(const qn::NetworkStruct<T>& sn, const MfqTopology& top, double tol) {
    const std::size_t K = sn.nclasses;
    if (!top.ok)
        throw UnsupportedError(
            "fluid mfq: the priority fluid queue needs a model of open classes flowing "
            "Source -> Queue -> Sink and nothing else");
    MfqPrioResult out;
    out.QN.assign(K, 0.0);
    out.RN.assign(K, 0.0);
    out.TN.assign(K, 0.0);
    out.UN.assign(K, 0.0);

    // Open classes, ordered so that the LAST row is the highest priority.
    std::vector<std::size_t> ord = top.open_classes;
    std::stable_sort(ord.begin(), ord.end(), [&sn](std::size_t a, std::size_t b) {
        return sn.classes[a].prio > sn.classes[b].prio;
    });
    const std::size_t Kc = ord.size();
    if (Kc == 0) {
        out.fallback = true;
        out.reason = "invalid service rates";
        return out;
    }

    // The drain rate must be one number: a priority fluid queue has ONE server.
    std::vector<double> mu(Kc, 0.0);
    for (std::size_t i = 0; i < Kc; ++i)
        mu[i] = num_traits<T>::to_double(sn.rates(top.queue, ord[i]));
    for (std::size_t i = 0; i < Kc; ++i)
        if (!std::isfinite(mu[i]) || mu[i] <= 0.0) {
            out.fallback = true;
            out.reason = "invalid service rates";
            return out;
        }
    const double fine_tol = 1e-8;  // GlobalConstants.FineTol
    for (std::size_t i = 0; i < Kc; ++i)
        if (std::fabs(mu[i] - mu[0]) > fine_tol * std::max(1.0, mu[0])) {
            out.fallback = true;
            out.reason = "class-dependent service";
            return out;
        }
    const double d = mu[0];

    // Per-class arrival fluid: Q_k = D0 + D1 and the rate vector R_k = D1 e.
    std::vector<Matrix<double>> Qk(Kc), Rk(Kc);
    std::vector<std::size_t> Nk(Kc, 0);
    std::vector<double> lambda(Kc, 0.0);
    for (std::size_t i = 0; i < Kc; ++i) {
        const lang::Distrib<T>& arr = sn.service[top.source][ord[i]];
        const std::size_t n = arr.D0.rows();
        if (n == 0 || arr.D1.rows() != n) {
            out.fallback = true;
            out.reason = "non-MAP arrival";
            return out;
        }
        Qk[i] = Matrix<double>(n, n, 0.0);
        Rk[i] = Matrix<double>(n, 1, 0.0);
        for (std::size_t a = 0; a < n; ++a)
            for (std::size_t b = 0; b < n; ++b) {
                Qk[i](a, b) = num_traits<T>::to_double(arr.D0(a, b)) +
                              num_traits<T>::to_double(arr.D1(a, b));
                Rk[i](a, 0) += num_traits<T>::to_double(arr.D1(a, b));
            }
        Nk[i] = n;
        lambda[i] = num_traits<T>::to_double(sn.rates(top.source, ord[i]));
    }

    std::size_t Njoint = 1;
    for (std::size_t i = 0; i < Kc; ++i) Njoint *= Nk[i];
    if (Njoint < 2) {
        out.fallback = true;
        out.reason = "non-modulated (exponential) arrivals";
        return out;
    }

    // The joint background chain, and the per-class rate in each of its states.
    Matrix<double> Qjoint(Njoint, Njoint, 0.0);
    for (std::size_t i = 0; i < Kc; ++i) {
        const Matrix<double> term = detail::prio_kron_term(Qk, Nk, i);
        for (std::size_t a = 0; a < Njoint; ++a)
            for (std::size_t b = 0; b < Njoint; ++b) Qjoint(a, b) += term(a, b);
    }
    Matrix<double> Rjoint(Kc, Njoint, 0.0);
    std::vector<Matrix<double>> ones_col(Kc);
    for (std::size_t j = 0; j < Kc; ++j) ones_col[j] = Matrix<double>(Nk[j], 1, 1.0);
    for (std::size_t i = 0; i < Kc; ++i) {
        Matrix<double> v(1, 1, 1.0);
        for (std::size_t j = 0; j < Kc; ++j) v = mam::kron(v, (j == i) ? Rk[i] : ones_col[j]);
        for (std::size_t a = 0; a < Njoint; ++a) Rjoint(i, a) = v(a, 0);
    }

    mam::FluidPrioOptions po;
    po.prec = (tol > 0.0) ? tol : 1e-14;
    po.classes.clear();
    for (std::size_t i = 1; i <= Kc; ++i) po.classes.push_back(i);
    mam::FluidPrioResult fl, st;
    try {
        po.flMoms = 1;
        po.stMoms = 0;
        fl = mam::mfq_prio_queue(Qjoint, Rjoint, d, po);
        po.flMoms = 0;
        po.stMoms = 1;
        st = mam::mfq_prio_queue(Qjoint, Rjoint, d, po);
    } catch (const std::exception& e) {
        out.fallback = true;
        out.reason = e.what();
        return out;
    }

    for (std::size_t i = 0; i < Kc; ++i) {
        const std::size_t k = ord[i];
        out.QN[k] = (i < fl.flMoms.size() && !fl.flMoms[i].empty()) ? fl.flMoms[i][0] : 0.0;
        out.RN[k] = (i < st.stMoms.size() && !st.stMoms[i].empty()) ? st.stMoms[i][0] : 0.0;
        out.TN[k] = lambda[i];
        out.UN[k] = std::min(1.0, lambda[i] / d);
    }
    return out;
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_MFQ_PRIO_H
