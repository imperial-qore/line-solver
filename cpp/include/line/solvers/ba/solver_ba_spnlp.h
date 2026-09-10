#pragma once
/**
 * @file solver_ba_spnlp.h
 * @brief Linear-programming bounds on the mean marking and the throughputs of a
 *        stochastic timed Petri net.
 *
 * Port of `matlab/src/solvers/BA/solver_ba_spnlp_analyzer.m`. The polytope and
 * the LP are `spn::spn_lpbnd`; this analyzer maps the LINE model onto them
 * and reads one side of the bracket back per place.
 *
 * METHOD NAMES. Four, in two families: `spnlp.upper` and `spnlp.lower` are the
 * Markovian LP and need exponential firing times; `spnlp.op.upper` and
 * `spnlp.op.lower` drop the second-moment, covariance and Little's-law families
 * and the whole E[X_p e_t] block with them, which is what removes the
 * exponential requirement and admits any phase-type law. The operational pair
 * is much looser, and is the reference's own "without Markovian assumption"
 * column.
 *
 * BOUND CONVENTION. Q(i,0) is the reported side of the bracket on the mean
 * number of tokens in place i. Tp(i,0) is the same side of the bracket on the
 * token throughput of that place, and R follows by Little's law from the two.
 * U(i,0) = Q(i,0) DELIBERATELY: a Place is an INF station and LINE reports
 * U = Q at an infinite server, which is what SolverCTMC and
 * `solver_nc_spn_analyzer` both do on the same net. The reference's place
 * utilization 1 - P(m = 0) is a different quantity and is not this column.
 *
 * SINGLE CLASS ONLY, and column 0 is the only one written, matching
 * `solver_nc_spn.h`: `spn_lpbnd` puts one level per place because
 * `NetworkStruct::transparam` carries no class dimension, and refuses a
 * coloured net rather than collapse it. The MATLAB, JAR and python twins carry
 * the class axis and do not have this restriction.
 *
 * A TRANSITION GETS NO ROW. It is a stateful node and not a station, so it has
 * no station index; the mode throughputs and enabling probabilities the LP also
 * brackets stay inside `spn_lpbnd`'s return value, the same way `spn_metrics`
 * keeps mode_tput and mode_util off the table.
 *
 * HOW TIGHT. The reference's own Table 2 measures it on a four-server
 * production line: the upper side lands 2% to 11% above simulation and the
 * lower side 30% to 40% below it, both comfortably inside the operational
 * bounds it also reports. Expect a usable upper bound and a weak lower one.
 *
 * Reference: Z. Liu (1998). Performance analysis of stochastic timed Petri nets
 * using linear programming approach. IEEE Transactions on Software Engineering
 * 24(11), 1014-1030.
 */

#include <cstddef>
#include <string>

#include "line/api/spn/spn_lpbnd.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/ba/solver_ba_analyzer.h"
#include "line/util/error.h"

namespace line {
namespace ba {

/** Whether a resolved method name belongs to the Petri-net LP family. */
inline bool is_spnlp_method(const std::string& method) {
    return method.rfind("spnlp", 0) == 0;
}

/**
 * Moment-relaxation LP bounds for a stochastic Petri net.
 *
 * @param L   the refreshed struct of a net holding Places and Transitions
 * @param opt the method; every other gate belongs to `spn_lpbnd`
 */
template <class T>
BaSolution<T> solver_ba_spnlp_analyzer(const qn::NetworkStruct<T>& L, const BaOptions& opt) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = L.nstations, K = L.nclasses;

    BaSolution<T> s;
    s.Q = Matrix<T>(M, K, zero);
    s.U = Matrix<T>(M, K, zero);
    s.R = Matrix<T>(M, K, zero);
    s.Tp = Matrix<T>(M, K, zero);
    s.C.assign(K, zero);
    s.X.assign(K, zero);
    s.iter = 1;

    bool markovian;
    bool upper;
    if (opt.method == "spnlp.upper") {
        markovian = true;
        upper = true;
    } else if (opt.method == "spnlp.lower") {
        markovian = true;
        upper = false;
    } else if (opt.method == "spnlp.op.upper") {
        markovian = false;
        upper = true;
    } else if (opt.method == "spnlp.op.lower") {
        markovian = false;
        upper = false;
    } else {
        throw UnsupportedError("solver_ba_spnlp_analyzer: unknown SPN bound method '" +
                               opt.method +
                               "'. Valid: spnlp.upper, spnlp.lower, spnlp.op.upper, "
                               "spnlp.op.lower");
    }

    // ---- model gates ----
    // Every other check belongs to `spn_lpbnd`, which refuses by name on the
    // mode it cannot represent. What must be decided here is only whether this
    // is a Petri net at all, and whether the places carry an embedded queue the
    // relaxation has no variable for.
    bool hasTransition = false;
    for (std::size_t i = 0; i < L.nodes.size() && !hasTransition; ++i)
        hasTransition = L.nodes[i].nodetype == lang::NodeType::Transition;
    if (!hasTransition)
        throw UnsupportedError("solver_ba_spnlp_analyzer: method '" + opt.method +
                               "' bounds a stochastic Petri net; this model has no Transition "
                               "node. Use the queueing-network bound families");
    for (std::size_t i = 0; i < L.nodes.size(); ++i) {
        if (L.nodes[i].nodetype != lang::NodeType::Place) continue;
        const std::size_t ist = L.nodes[i].station;
        if (ist >= 1 && ist <= M && L.stations[ist - 1].sched != SchedStrategy::INF)
            throw UnsupportedError("solver_ba_spnlp_analyzer: method '" + opt.method +
                                   "' does not support queueing places: place " + L.nodes[i].name +
                                   " serves under a non-INF discipline, and the relaxation "
                                   "carries one variable per place marking with no notion of an "
                                   "embedded queue");
    }

    spn::SpnLpOptions lpopt;
    lpopt.markovian = markovian;
    const spn::SpnLpBounds bnd = spn::spn_lpbnd(L, lpopt);

    for (std::size_t pp = 0; pp < bnd.places.size(); ++pp) {
        const std::size_t ist = L.nodes[bnd.places[pp] - 1].station;
        if (ist < 1 || ist > M) continue;
        const double q = upper ? bnd.tokens_hi[pp] : bnd.tokens_lo[pp];
        const double t = upper ? bnd.place_tput_hi[pp] : bnd.place_tput_lo[pp];
        s.Q(ist - 1, 0) = num_traits<T>::from_double(q);
        s.U(ist - 1, 0) = num_traits<T>::from_double(q);
        s.Tp(ist - 1, 0) = num_traits<T>::from_double(t);
        if (t > 0) s.R(ist - 1, 0) = num_traits<T>::from_double(q / t);
    }

    if (!L.classes.empty()) {
        const std::size_t ref = L.classes[0].refstat;
        if (ref >= 1 && ref <= M) s.X[0] = s.Tp(ref - 1, 0);
        T nk = zero;
        for (std::size_t i = 0; i < M; ++i) nk += s.Q(i, 0);
        if (num_traits<T>::to_double(s.X[0]) > 0 && num_traits<T>::to_double(nk) > 0)
            s.C[0] = nk / s.X[0];
    }
    return s;
}

}  // namespace ba
}  // namespace line
