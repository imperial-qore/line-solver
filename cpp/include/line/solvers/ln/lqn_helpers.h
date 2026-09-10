/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_LN_LQN_HELPERS_H
#define LINE_SOLVERS_LN_LQN_HELPERS_H

/**
 * Standalone LQN routines that SolverLN needs but does not contain.
 *
 * Port of matlab/src/solvers/LN/lqn_fwd_rendezvous.m and
 * matlab/src/solvers/LN/lqn_overtake_markov.m.
 *
 * lqn_act_thinktime.m is DELIBERATELY ABSENT. Its whole content -- add the
 * activity think time to the activity's service and residence time, skipping
 * the unset (NaN in MATLAB, `disabled` here) case -- is already inline in the
 * servt_map loop of solver_ln.h. A second copy of a two-line rule that two
 * files would have to keep agreeing on is worse than no helper at all.
 *
 * BOTH ARE WIRED IN. SolverLN::construct calls lqn_fwd_rendezvous before it
 * reads anything else out of the struct, as @@SolverLN/SolverLN.m:162-164 does,
 * and lqn_overtake_markov is reached from update_metrics through the input
 * mapping of lqn_analyzers.h's lqn_overtake_prob_markov. Neither construct is
 * refused by name any longer.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/lqn/lqn_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ln {

using lang::CallType;
using lqn::LqnStruct;

// ---------------------------------------------------------------------------
// lqn_fwd_rendezvous
// ---------------------------------------------------------------------------

/**
 * Replace every forwarding chain reachable from a synchronous call by
 * caller-side pseudo rendezvous calls to the forwarding targets.
 *
 * Port of LQNS Phase::addForwardingRendezvous (phase.cc). A forwarded call
 * blocks the original caller until the LAST task in the chain replies, so the
 * caller's blocking time spans the chain, not just the entry it named. Rather
 * than teach every downstream stage what a FWD arc means, the chain is
 * flattened here into ordinary SYNC arcs from the original calling activity to
 * each entry on the chain, each with mean equal to the original call mean times
 * the product of the forwarding probabilities on the path to it. Layer
 * construction, think times, populations and the interlock analysis then see
 * plain rendezvous arcs and need no forwarding case at all. This is also what
 * LQNS does: interlock.cc drops FWD arcs outright ("Drop forward -- keep rnv")
 * and accounts for forwarding only through these pseudo arcs.
 *
 * The FWD calls survive in the struct but must not contribute blocking anywhere
 * after this point, or the chain is charged twice.
 *
 * Asynchronous calls into a forwarding chain are left alone: LQNS breaks the
 * backward search at a send-no-reply, since nobody is blocked waiting for it.
 *
 * The MATLAB version also rebuilds the Geometric process descriptor of each
 * rewritten call. There is nothing to port: LqnStruct carries only
 * callproc_mean, which is the only field SolverLN reads.
 */
template <class T>
void lqn_fwd_rendezvous(LqnStruct<T>& lqn) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    bool any_fwd = false;
    for (std::size_t c = 1; c <= lqn.ncalls; ++c)
        if (lqn.calltype[c] == CallType::FWD) any_fwd = true;
    if (!any_fwd) return;

    // Frozen: the pseudo calls appended below are themselves SYNC, and
    // rewriting them again would raise the chain to a power of its own
    // probabilities.
    const std::size_t ncalls0 = lqn.ncalls;

    for (std::size_t cidx = 1; cidx <= ncalls0; ++cidx) {
        if (lqn.calltype[cidx] != CallType::SYNC) continue;
        const std::size_t aidx = lqn.callpair_src[cidx];
        const std::size_t tidx = lqn.parent[aidx];
        const T base_mean = lqn.callproc_mean[cidx];
        if (base_mean <= zero) continue;

        // breadth-first walk of the forwarding chain out of the sync target,
        // carrying the probability of the path that reached each entry
        std::vector<std::size_t> frontier{lqn.callpair_dst[cidx]};
        std::vector<T> probs{one};
        std::vector<std::size_t> visited;
        auto seen = [](const std::vector<std::size_t>& v, std::size_t x) {
            for (std::size_t e : v)
                if (e == x) return true;
            return false;
        };

        while (!frontier.empty()) {
            const std::size_t eidx = frontier.front();
            const T p_path = probs.front();
            frontier.erase(frontier.begin());
            probs.erase(probs.begin());
            if (seen(visited, eidx)) continue;
            visited.push_back(eidx);

            for (std::size_t fcidx = 1; fcidx <= ncalls0; ++fcidx) {
                if (lqn.calltype[fcidx] != CallType::FWD) continue;
                if (lqn.callpair_src[fcidx] != eidx) continue;
                const T fprob = lqn.callproc_mean[fcidx];
                const std::size_t tgt = lqn.callpair_dst[fcidx];
                const T pseudo_mean = T(base_mean * p_path * fprob);

                // A chain that comes back to the caller's own task carries no
                // pseudo arc: the caller would appear as its own client.
                if (pseudo_mean > zero && lqn.parent[tgt] != tidx) {
                    // see _kb/06-solver-catalog.md (LN section) for rationale
                    std::size_t mrow = 0;
                    for (std::size_t scan = 1; scan <= lqn.ncalls; ++scan)
                        if (lqn.calltype[scan] == CallType::SYNC &&
                            lqn.callpair_src[scan] == aidx && lqn.callpair_dst[scan] == tgt) {
                            mrow = scan;
                            break;
                        }
                    if (mrow > 0) {
                        // an arc already exists, so the forwarded work is extra
                        // visits on it rather than a second parallel class
                        lqn.callproc_mean[mrow] = T(lqn.callproc_mean[mrow] + pseudo_mean);
                    } else {
                        const std::size_t ncall = lqn.ncalls + 1;
                        lqn.ncalls = ncall;
                        const std::size_t target_tidx = lqn.parent[tgt];
                        lqn.calltype.push_back(CallType::SYNC);
                        lqn.callpair_src.push_back(aidx);
                        lqn.callpair_dst.push_back(tgt);
                        lqn.callproc_mean.push_back(pseudo_mean);
                        lqn.callnames.push_back(lqn.names[aidx] + "=>" + lqn.names[tgt]);
                        lqn.callhashnames.push_back(lqn.hashnames[aidx] + "=>" +
                                                    lqn.hashnames[tgt]);
                        lqn.callsof[aidx].push_back(ncall);
                        lqn.iscaller.set(tidx, target_tidx);
                        lqn.iscaller.set(aidx, target_tidx);
                        lqn.iscaller.set(tidx, tgt);
                        lqn.iscaller.set(aidx, tgt);
                        lqn.issynccaller.set(tidx, target_tidx);
                        lqn.issynccaller.set(aidx, target_tidx);
                        lqn.issynccaller.set(tidx, tgt);
                        lqn.issynccaller.set(aidx, tgt);
                        lqn.graph.set(aidx, tgt, one);
                        lqn.taskgraph.set(tidx, target_tidx, one);
                    }
                }

                if (!seen(visited, tgt) && !seen(frontier, tgt)) {
                    frontier.push_back(tgt);
                    probs.push_back(T(p_path * fprob));
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// lqn_overtake_markov
// ---------------------------------------------------------------------------

namespace detail {

/** True when x is a finite value of T; only a floating T can fail this. */
template <class T>
bool ln_finite(const T& x) {
    return std::isfinite(num_traits<T>::to_double(x));
}

/**
 * Transition rates of one client phase in the LQNS slice chain (slice.cc
 * setRates). The four outputs are the probabilities of, respectively, staying
 * in the client phase (a), completing it and revisiting the server (b),
 * catching the server in the tested phase (c), and being held at another server
 * first (d).
 */
template <class T>
void ln_set_rates(const T& xj, const T& prA, const T& nSlices, const T& service, const T& y_ij,
                  const T& y_ik, const T& t_k, T& a, T& b, T& c, T& d) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    const T y_sum = T(y_ij + y_ik + one);
    const T slice = nSlices != zero ? T(service / nSlices) : zero;

    T q0, q1, q3, q5;
    T temp = T(xj + t_k);
    if (!ln_finite(temp)) {
        q0 = one;
        q3 = zero;
    } else if (temp != zero) {
        q0 = T(xj / temp);
        q3 = T(t_k / temp);
    } else {
        q0 = zero;
        q3 = one;
    }
    temp = T(xj + slice);
    if (!ln_finite(temp)) {
        q1 = zero;
        q5 = one;
    } else if (temp != zero) {
        q1 = T(xj / temp);
        q5 = T(slice / temp);
    } else {
        q1 = one;
        q5 = zero;
    }
    const T q2 = T(y_ik / y_sum), q4 = T(y_ij / y_sum), q6 = T(prA / y_sum);

    a = T(q5 + q1 * q2 * q3);
    b = T(q1 * q6);
    c = T(q1 * q4);
    d = T(q0 * q1 * q2);
}

}  // namespace detail

/**
 * Overtaking probability from the LQNS phased-server Markov chain.
 *
 * Layer-1 port of LQNS V6 (slice.cc setRates and prOvertakingStates,
 * overtake.cc computeOvertaking) for the single-conditioning case, where the
 * calling entry and the conditioning entry coincide. Overtaking is the event
 * that a client's next request reaches the server while the server is still
 * running the second phase of the PREVIOUS request from the same client: the
 * early reply released the client, so two of its requests can be in flight at
 * once and the later one can pass the earlier one. `lqns -t overtaking` prints
 * the same quantity; on its 31-overtaking model, phase 2 of the server gives
 * 0.5.
 *
 * The chain is over client phases 0..maxPhaseA, phase 0 being the client's
 * think slice. Row p of `clientPhases` is
 *   [nSlices, service, y_ij, y_ik, t_k]
 * with nSlices = 1 + the rendezvous calls made in phase p (the slice count the
 * phase is chopped into), service = total host residence of the phase, y_ij =
 * calls to the tested server's task, y_ik = calls to any other task, and t_k =
 * mean time spent at those other tasks. Row 0's service is the think time.
 * `xj` is the server's residence time in the phase being tested, `prVisit` the
 * client entry's visit probability, and y_aj[0] the client's total calls to the
 * server with y_aj[i] the calls made in client phase i.
 *
 * Every rate is a ratio of times (xj/(xj+slice) and so on), so the result is
 * invariant to a common rescaling of all the time inputs, as a probability must
 * be.
 *
 * Reference: Franks and Woodside, "Effectiveness of early replies in
 * client-server systems", Perf. Eval. 36 (1999).
 */
template <class T>
T lqn_overtake_markov(const Matrix<T>& clientPhases, const T& prVisit, const T& xj,
                      const std::vector<T>& y_aj) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    if (clientPhases.rows() < 1 || clientPhases.cols() < 5)
        throw InputError(
            "lqn_overtake_markov: clientPhases must have five columns "
            "[nSlices service y_ij y_ik t_k] and one row per client phase 0..maxPhaseA");
    const std::size_t nStates = clientPhases.rows();
    const std::size_t maxPhaseA = nStates - 1;
    if (y_aj.size() < nStates)
        throw InputError(
            "lqn_overtake_markov: y_aj must hold the total call count and one entry per client "
            "phase");

    // setRates for each client phase. prVisit applies only to the last phase,
    // which is the one that ends the client cycle and issues the next request.
    std::vector<T> a(nStates, zero), b(nStates, zero), c(nStates, zero), d(nStates, zero);
    for (std::size_t p = 0; p < nStates; ++p) {
        const T prA = (p == maxPhaseA) ? prVisit : one;
        detail::ln_set_rates(xj, prA, clientPhases(p, 0), clientPhases(p, 1), clientPhases(p, 2),
                             clientPhases(p, 3), clientPhases(p, 4), a[p], b[p], c[p], d[p]);
    }

    // prOvertakingStates. prod_of_b folds the self-loop at a phase (probability
    // d of a detour to another server) into the phase's forward probability;
    // the denominator closes the cycle over all phases, so `product` is the
    // stationary weight of arriving in phase r having started the cycle in
    // phase i. Both denominators are structurally positive for non-negative
    // inputs, since d = q0*q1*q2 < 1 whenever y_ij + 1 > 0.
    // `next` is the second plane of the reference's PrOT array. Nothing reads it
    // in the single-conditioning case solved here; it is kept because it is half
    // of the state pair the chain produces and dropping it would leave the port
    // silently unable to answer the multi-conditioning case.
    auto prod_of_b = [&](std::size_t k) { return T(b[k] / (one - d[k])); };
    Matrix<T> over(nStates, nStates, zero), next(nStates, nStates, zero);
    for (std::size_t i0 = 0; i0 < nStates; ++i0) {
        T temp = one;
        for (std::size_t r0 = 0; r0 < nStates; ++r0)
            if (r0 != i0) temp = T(temp * prod_of_b(r0));
        T product = T(one / (one - (b[i0] * temp + d[i0])));
        std::size_t r0 = i0;
        for (;;) {
            over(i0, r0) = T(c[i0] * product);
            next(i0, r0) = T(a[i0] * product);
            r0 = (r0 == 0) ? maxPhaseA : r0 - 1;
            product = T(product * prod_of_b(r0));
            if (r0 == i0) break;
        }
    }

    // computeOvertaking with entA == entC. The client's next request leaves
    // from its last phase, so nextProb puts all mass there; the y_aj ratio
    // conditions on the request having been issued in phase i.
    std::vector<T> nextProb(nStates, zero);
    if (maxPhaseA >= 1) nextProb[maxPhaseA] = one;

    T prOt = zero;
    for (std::size_t i = 1; i <= maxPhaseA; ++i) {
        if (clientPhases(i, 2) == zero) continue;  // phase i never calls this server
        T acc = zero;
        for (std::size_t r = 1; r <= maxPhaseA; ++r) acc = T(acc + nextProb[r] * over(i, r));
        // y_aj[i] cannot vanish here: it counts the same calls as clientPhases(i,2)
        prOt = T(prOt + acc * (y_aj[0] / y_aj[i]));
    }
    return prOt;
}

}  // namespace ln
}  // namespace line

#endif  // LINE_SOLVERS_LN_LQN_HELPERS_H
