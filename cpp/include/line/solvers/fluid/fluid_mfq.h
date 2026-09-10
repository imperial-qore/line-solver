/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_MFQ_H
#define LINE_SOLVERS_FLUID_FLUID_MFQ_H

/**
 * The `mfq` method: a port of `solver_mfq.m` and the single-queue gate
 * `fluid_is_single_queue.m`.
 *
 * WHAT MAKES THIS DIFFERENT FROM EVERY OTHER FLUID METHOD HERE. The others
 * approximate a network by following the mean drift of its queues. This one is
 * EXACT, and only works on one queue: a Source feeding a single station whose
 * arrival and service processes are Markov-modulated fluids. It solves the
 * fluid queue analytically rather than integrating anything, which is why it
 * carries a topology gate instead of a tolerance -- there is nothing to
 * converge.
 *
 * WHAT IT SOLVES. Both processes are given as (D0, D1) pairs and converted to
 * BuTools fluid form: the background generator is Q = D0 + D1 and the fluid
 * rate in each phase is the row sum of D1. `mfq_fluflu_sojourn` (already ported
 * from BuTools' FluFluQueue) then returns a matrix-exponential representation
 * of the sojourn time, whose mean is the response time. Queue length follows by
 * Little's law, which is exact in steady state, and throughput is the arrival
 * rate the queue is stable under.
 *
 * THE M/M/1 SHORT CIRCUIT is the reference's. With a single-phase arrival and a
 * single-phase service there is no modulation left, the fluid machinery is
 * degenerate, and the reference falls back to the textbook formulas
 * L = rho/(1-rho), W = 1/(mu - lambda). Reproduced, including the instability
 * check: at rho >= 1 the queue has no stationary distribution and the reference
 * reports infinities rather than a finite wrong answer.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/mfq_fluflu_sojourn.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace fluid {

/** What the single-queue gate found, when it matches. */
struct MfqTopology {
    bool ok = false;
    std::size_t source = 0;  ///< 0-based station index
    std::size_t queue = 0;
    std::size_t cls = 0;                      ///< the first open class, which `mfq` analyzes
    std::vector<std::size_t> open_classes;    ///< every open class, in class order
};

/**
 * Port of `fluid_is_single_queue.m`: the model must be one open class flowing
 * Source -> Queue -> Sink and nothing else.
 */
template <class T>
MfqTopology mfq_is_single_queue(const qn::NetworkStruct<T>& sn) {
    MfqTopology t;
    const std::size_t M = sn.nstations, K = sn.nclasses;
    std::size_t nsrc = 0, nq = 0, nsink = 0;
    for (std::size_t i = 0; i < M; ++i) {
        switch (sn.stations[i].nodetype) {
            case qn::NodeType::Source: ++nsrc; t.source = i; break;
            case qn::NodeType::Sink: ++nsink; break;
            case qn::NodeType::Queue: ++nq; t.queue = i; break;
            default: return t;  // a Delay, Cache or anything else disqualifies
        }
    }
    // A Sink is a node but not always a station in this port, so it is
    // counted when present and not required.
    (void)nsink;
    if (nsrc != 1 || nq != 1) return t;
    // The reference requires AT LEAST one open class, not exactly one: a
    // multiclass single queue is where the priority branch takes over.
    for (std::size_t r = 0; r < K; ++r) {
        if (std::isfinite(sn.classes[r].population)) return t;  // closed class
        t.open_classes.push_back(r);
    }
    if (t.open_classes.empty()) return t;
    t.cls = t.open_classes[0];
    t.ok = true;
    return t;
}

/** The metrics `mfq` reports for its single queue. */
struct MfqResult {
    double QN = 0.0, RN = 0.0, TN = 0.0, UN = 0.0;
    bool unstable = false;
};

namespace detail {

/** Port of `fluid_dist2butools`: Q = D0 + D1, R = diag(row sums of D1). */
template <class T>
void mfq_dist2butools(const lang::Distrib<T>& d, Matrix<double>& Q, Matrix<double>& R) {
    const std::size_t n = d.D0.rows();
    Q = Matrix<double>(n, n, 0.0);
    R = Matrix<double>(n, n, 0.0);
    for (std::size_t a = 0; a < n; ++a) {
        double row = 0.0;
        for (std::size_t b = 0; b < n; ++b) {
            Q(a, b) = num_traits<T>::to_double(d.D0(a, b)) + num_traits<T>::to_double(d.D1(a, b));
            row += num_traits<T>::to_double(d.D1(a, b));
        }
        R(a, a) = row;
    }
}

/**
 * Mean of a matrix-exponential representation: -alpha A^-1 e.
 *
 * The sojourn representation is (alpha, A) with density alpha exp(A t) (-A e),
 * so the first moment is -alpha A^-1 e; solving against A is what avoids
 * forming the inverse.
 */
inline double mfq_me_mean(const mam::MeRepresentation<double>& me) {
    const std::size_t n = me.A.rows();
    if (n == 0) return 0.0;
    std::vector<double> e(n, 1.0);
    // Solve A y = e, then the mean is -alpha . y.
    Matrix<double> A = me.A;
    const std::vector<std::size_t> piv = lu_factor(A);
    std::vector<double> y = e;
    lu_solve(A, piv, y);
    double s = 0.0;
    for (std::size_t i = 0; i < n && i < me.alpha.size(); ++i) s += me.alpha[i] * y[i];
    return -s;
}

}  // namespace detail

/** Solve the single fluid queue of `sn`. */
template <class T>
MfqResult fluid_mfq(const qn::NetworkStruct<T>& sn, const MfqTopology& top, double tol) {
    if (!top.ok)
        throw UnsupportedError(
            "fluid mfq: the method solves a single Markov-modulated fluid queue and needs a model "
            "of exactly one open class flowing Source -> Queue -> Sink");
    const lang::Distrib<T>& arr = sn.service[top.source][top.cls];
    const lang::Distrib<T>& svc = sn.service[top.queue][top.cls];

    mam::Map<double> am, sm;
    const std::size_t na = arr.D0.rows(), ns = svc.D0.rows();
    am.D0 = Matrix<double>(na, na, 0.0);
    am.D1 = Matrix<double>(na, na, 0.0);
    sm.D0 = Matrix<double>(ns, ns, 0.0);
    sm.D1 = Matrix<double>(ns, ns, 0.0);
    for (std::size_t a = 0; a < na; ++a)
        for (std::size_t b = 0; b < na; ++b) {
            am.D0(a, b) = num_traits<T>::to_double(arr.D0(a, b));
            am.D1(a, b) = num_traits<T>::to_double(arr.D1(a, b));
        }
    for (std::size_t a = 0; a < ns; ++a)
        for (std::size_t b = 0; b < ns; ++b) {
            sm.D0(a, b) = num_traits<T>::to_double(svc.D0(a, b));
            sm.D1(a, b) = num_traits<T>::to_double(svc.D1(a, b));
        }
    const double lambda = 1.0 / mam::map_mean(am);
    const double mu = 1.0 / mam::map_mean(sm);

    MfqResult out;
    out.TN = lambda;

    Matrix<double> Qin, Rin, Qout, Rout;
    detail::mfq_dist2butools(arr, Qin, Rin);
    detail::mfq_dist2butools(svc, Qout, Rout);

    // The degenerate case the reference short-circuits: no modulation left.
    const bool simple_exp = (na == 1 && ns == 1);
    if (simple_exp) {
        const double rho = lambda / mu;
        if (rho >= 1.0) {
            out.unstable = true;
            out.QN = std::numeric_limits<double>::infinity();
            out.RN = std::numeric_limits<double>::infinity();
            out.UN = 1.0;
            return out;
        }
        out.QN = rho / (1.0 - rho);
        out.RN = 1.0 / (mu - lambda);
        out.UN = rho;
        return out;
    }

    // The genuinely Markov-modulated case: solve the fluid-fluid queue.
    //
    // A NOTE ON WHAT THIS MEASURES, because it is easy to check against the
    // wrong formula. This is a FLUID queue: work arrives as a continuous
    // stream at a Markov-modulated rate and drains at another, so the sojourn
    // is the delay of a DROP OF FLUID, not a customer's waiting time in the
    // corresponding M/G/1 system. On M/E2/1 with lambda 0.5 and mean service
    // 0.5 the two differ by more than tenfold (0.0417 against 0.625), and the
    // fluid figure is the right one here -- MATLAB's mfq reports exactly
    // 0.0417 on that model.
    const mam::MeRepresentation<double> me =
        mam::mfq_fluflu_sojourn(Qin, Rin, Qout, Rout, /*srv0stop=*/true, /*transToPH=*/false,
                                std::max(tol, 1e-14));
    out.RN = detail::mfq_me_mean(me);
    out.QN = lambda * out.RN;  // Little's law, exact in steady state
    // The reference reports the fluid LEVEL as the utilization on this branch,
    // not lambda/mu: the server is busy exactly while there is fluid to drain.
    out.UN = out.QN;
    (void)mu;
    return out;
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_MFQ_H
