/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_BA_SOLVER_BA_ANALYZER_H
#define LINE_SOLVERS_BA_SOLVER_BA_ANALYZER_H

/**
 * Port of `matlab/src/solvers/BA/solver_ba_analyzer.m`, the bound-analysis
 * handler behind SolverBA.
 *
 * Each method returns ONE side of a bracket, optimistic or pessimistic, never a
 * point estimate. The bound itself is a scalar on the chain throughput; the
 * per-station [Q,U,R,T,C] around it is reconstructed by the utilization law and
 * by the same optimistic/pessimistic residence convention the ABA bound uses,
 * which is why every family agrees on the shape of the answer even where the
 * literature leaves the per-station response time undefined.
 *
 * WHAT A BOUND NEEDS is demands (visits x mean service time) and populations.
 * No service-time distribution enters, so the featset is deliberately narrow:
 * closed, product-form-parameterized models. The single-class families reject a
 * multiclass model by name rather than answer for one class.
 *
 * ARITHMETIC. The asymptotic, balanced, proportional and box families are sums,
 * products, integer powers, minima and divisions, so they are exact in rational
 * arithmetic. Three paths are not and are gated: the geometric family (gb) and
 * the successively-improving bounds (sib) solve a quadratic, and the sb lower
 * bound takes an (N-1)-st root.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "line/api/mapqn/mapqn_bnd_lr_mva.h"
#include "line/api/mapqn/mapqn_bnd_lr_pf.h"
#include "line/api/pfqn/pfqn_cbh.h"
#include "line/api/pfqn/pfqn_harel_bounds.h"
#include "line/api/pfqn/pfqn_ldbcmp.h"
#include "line/api/pfqn/pfqn_looping.h"
#include "line/api/pfqn/pfqn_mcub.h"
#include "line/api/pfqn/pfqn_mwrbb.h"
#include "line/api/pfqn/pfqn_pbh.h"
#include "line/api/pfqn/pfqn_qzgblow.h"
#include "line/api/pfqn/pfqn_qzgbup.h"
#include "line/api/pfqn/pfqn_sib.h"
#include "line/api/pfqn/pfqn_scb.h"
#include "line/api/pfqn/pfqn_ssd.h"
#include "line/api/pfqn/pfqn_xzgsblow.h"
#include "line/api/pfqn/pfqn_xzgsbup.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/solvers/ba/solver_ba_bgt.h"
#include "line/solvers/ba/solver_ba_bpt.h"
#include "line/solvers/ba/solver_ba_snc.h"
#include "line/solvers/mva/sn_chain.h"

namespace line {
namespace ba {

using lang::SchedStrategy;

/** The options SolverBA reads. Defaults are `SolverBA.defaultOptions`. */
struct BaOptions {
    /** Bound method; `default` resolves to `gb.upper` in the runner. */
    std::string method = "default";
    /**
     * `options.level`: the hierarchy level of pbh/cbh/sib and the iteration
     * count k of pbk/bjbk. MATLAB's default is 2 and is NOT the pfqn default,
     * which is 1 for pbh and 3 for sib -- the solver overrides both.
     */
    int level = 2;
    /**
     * `options.config.qrf_alpha`, the (nstations x N) load-dependent scaling of
     * the two load-dependent QRF arms. Empty means all ones.
     */
    Matrix<double> qrf_alpha;
    /**
     * `options.config.qrf_params`, the blocking tables the BAS and RS-RD arms
     * need. `sn_to_qrf_params` in the reference refuses rather than defaulting:
     * assuming no blocking puts the bound ~31x farther from exact.
     */
    struct QrfParams {
        bool supplied = false;
        int f = 1;   ///< finite-capacity queue, 1-based as in the reference
        int MR = 1;  ///< number of blocking configurations
        std::vector<std::vector<int> > BB;   ///< (MR x M) blocking state
        std::vector<std::vector<int> > MM;   ///< (MR x 2) blocking order
        std::vector<std::vector<int> > MM1;  ///< (MR x M) extended order
        std::vector<int> ZZ;                 ///< (MR) blocked count per config
        std::vector<int> F;                  ///< (M) capacity; empty takes sn.cap
    };
    QrfParams qrf_params;
};

/** Class-level results, the [Q,U,R,T,C,X] of `solver_ba_analyzer`. */
template <class T>
struct BaSolution {
    Matrix<T> Q, U, R, Tp;
    std::vector<T> C, X;
    /** `-N log X`, the reference's approximate normalizing constant. */
    double lG = 0.0;
    int iter = 1;
};

namespace detail {

/** Station-indexed visit vector of the single chain of a single-class model. */
template <class T>
std::vector<T> ba_station_visits(const qn::NetworkStruct<T>& L) {
    // THE CONVERSION IS LOAD-BEARING, not a formality: `L.visits` is
    // STATEFUL-indexed (nstateful x nclasses) while `L.rates`, `L.stations` and
    // the sched are STATION-indexed. Every station is stateful, but not every
    // stateful node is a station -- a Transition is stateful without being one,
    // a Place is both -- so on a Petri net the two lengths differ (4 against 7
    // on the fork-join SPN). Walking the VISIT rows and indexing the station
    // arrays with the row is the bug the other three codebases carried into
    // their degeneracy predicate; do not "simplify" this loop into it.
    std::vector<T> V(L.nstations, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < L.nstations; ++i)
        V[i] = L.visits[0](L.stateful_of_station(i + 1) - 1, 0);
    return V;
}

template <class T>
bool ba_is_delay(const qn::NetworkStruct<T>& L, std::size_t i) {
    return L.stations[i].sched == SchedStrategy::INF;
}

/**
 * Port of the local `ba_sc_demands`: the visit vector, the aggregate think
 * time, the per-queue demand vector and the closed population.
 */
template <class T>
struct ScDemands {
    std::vector<T> V;  ///< (M) station-indexed visits
    std::vector<T> D;  ///< (Mq) demands of the queueing stations, in station order
    T Z;               ///< aggregate think time of the delay stations
    long N;            ///< closed population
};

template <class T>
ScDemands<T> ba_sc_demands(const qn::NetworkStruct<T>& L) {
    // No single-class test here: `method_refusal` owns that rule for every
    // caller, and a second copy is what would let the report and the run
    // disagree. `solver_ba_analyzer` has already asked it by the time we arrive.
    ScDemands<T> d;
    d.V = ba_station_visits(L);
    d.Z = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < L.nstations; ++i) {
        const T di = T(d.V[i] / L.rates(i, 0));
        if (ba_is_delay(L, i))
            d.Z += di;
        else
            d.D.push_back(di);
    }
    d.N = static_cast<long>(L.nclosedjobs());
    return d;
}

/**
 * The per-station metrics implied by a scalar chain-throughput bound and a
 * system response time, shared by every family (`ba_fill` in the reference,
 * inlined in each noniterative branch there).
 *
 * The optimistic side charges every queueing station the full-contention
 * response time N/mu and the pessimistic side the no-contention 1/mu; a delay
 * station takes 1/mu on both, having no contention to bound.
 */
template <class T>
BaSolution<T> ba_fill_xc(const qn::NetworkStruct<T>& L, const std::vector<T>& V, long N, const T& X,
                         const T& C, bool is_upper) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations;
    BaSolution<T> s;
    s.Q = Matrix<T>(M, 1, zero);
    s.U = Matrix<T>(M, 1, zero);
    s.R = Matrix<T>(M, 1, zero);
    s.Tp = Matrix<T>(M, 1, zero);
    for (std::size_t i = 0; i < M; ++i) {
        s.Tp(i, 0) = T(V[i] * X);
        const T svc = T(one / L.rates(i, 0));
        s.R(i, 0) = (is_upper && !ba_is_delay(L, i)) ? T(svc * num_traits<T>::from_int(N)) : svc;
        s.Q(i, 0) = T(s.Tp(i, 0) * s.R(i, 0));
        if (ba_is_delay(L, i)) {
            s.U(i, 0) = s.Q(i, 0);
        } else {
            // Utilization law per SERVER: without the nservers divisor a
            // multiserver station reports U > 1 (ssd/ldbcmp/auto reach here)
            const T c = num_traits<T>::from_double(std::max(1.0, L.stations[i].nservers));
            s.U(i, 0) = T(s.Tp(i, 0) / T(c * L.rates(i, 0)));
        }
    }
    s.C.assign(1, C);
    s.X.assign(1, X);
    s.lG = -static_cast<double>(N) * num_traits<T>::log_as_double(X);
    return s;
}

/**
 * Port of the local `ba_chain_qfill`, the chain-indexed queue length of the two
 * multiclass families.
 *
 * NEITHER RESIDENCE THE BOUND CARRIES YIELDS A QUEUE LENGTH ON THE DECLARED
 * SIDE: the optimistic no-contention residence understates it and the Theorem-1
 * residence of mwrbb overstates it, so the reference builds Q from the
 * utilizations instead. The pessimistic side takes U itself, which is the jobs
 * in service and so a lower bound on the jobs present; the optimistic side
 * charges the whole chain population to the station in proportion to its
 * saturation, capped at one.
 */
template <class T>
Matrix<T> ba_chain_qfill(const Matrix<T>& Uchain, const std::vector<T>& Xchain,
                         const Matrix<T>& Lchain, const std::vector<T>& Nchain,
                         const std::vector<bool>& isdelay, bool is_upper) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = Uchain.rows(), C = Uchain.cols();
    std::vector<T> Utot(M, zero);
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t c = 0; c < C; ++c) Utot[i] += Uchain(i, c);
        if (Utot[i] > one) Utot[i] = one;
    }
    Matrix<T> Qchain(M, C, zero);
    for (std::size_t c = 0; c < C; ++c)
        for (std::size_t i = 0; i < M; ++i) {
            if (isdelay[i])
                Qchain(i, c) = T(Xchain[c] * Lchain(i, c));
            else if (is_upper)
                Qchain(i, c) = T(Nchain[c] * Utot[i]);
            else
                Qchain(i, c) = Uchain(i, c);
        }
    return Qchain;
}

/** Port of the local `ba_fill`, whose cycle time is the ABA one at this side. */
template <class T>
BaSolution<T> ba_fill(const qn::NetworkStruct<T>& L, const std::vector<T>& V, long N, const T& X,
                      const T& Z, const std::vector<T>& D, bool is_upper) {
    T Dsum = num_traits<T>::from_int(0);
    for (const T& d : D) Dsum += d;
    const T C = is_upper ? T(Z + num_traits<T>::from_int(N) * Dsum) : T(Z + Dsum);
    return ba_fill_xc(L, V, N, X, C, is_upper);
}

/** min, max, sum and mean of the demand vector, which every family needs. */
template <class T>
struct DemandStats {
    T sum, max, mean;
};

template <class T>
DemandStats<T> ba_demand_stats(const std::vector<T>& D) {
    if (D.empty()) throw UnsupportedError("solver_ba_analyzer: the model has no queueing station");
    DemandStats<T> st;
    st.sum = num_traits<T>::from_int(0);
    st.max = D[0];
    for (const T& d : D) {
        st.sum += d;
        if (d > st.max) st.max = d;
    }
    st.mean = T(st.sum / num_traits<T>::from_int(static_cast<long>(D.size())));
    return st;
}

template <class T>
T ba_min(const T& a, const T& b) {
    return a < b ? a : b;
}

/**
 * Port of the local `mwrbb_disc_code`.
 *
 * The reference additionally maps the priority-scheduling variants of PS/DPS/GPS
 * and the preemptive-resume priority disciplines; neither family exists in this
 * port's SchedStrategy, so no input can reach those codes and the default
 * (ABA full-contention, code 4) covers every remaining work-conserving case.
 */
inline pfqn::MwrbbSched mwrbb_disc_code(SchedStrategy s) {
    switch (s) {
        case SchedStrategy::FCFS:
            return pfqn::MwrbbSched::Fifo;
        case SchedStrategy::PS:
        case SchedStrategy::DPS:
        case SchedStrategy::GPS:
            return pfqn::MwrbbSched::Ps;
        case SchedStrategy::HOL:
            return pfqn::MwrbbSched::PrioNonPreemptive;
        default:
            return pfqn::MwrbbSched::Aba;
    }
}

/** The chain-level scaffolding the two multiclass families share. */
template <class T>
struct ChainView {
    mva::ChainDemands<T> d;
    std::vector<bool> isdelay;
    std::vector<std::size_t> qstat;  ///< 0-based station index of each queueing station
    std::vector<T> Nv, Zc;
};

template <class T>
ChainView<T> ba_chain_view(const qn::NetworkStruct<T>& L) {
    ChainView<T> v;
    v.d = mva::sn_get_demands_chain(L);
    const std::size_t M = L.nstations, C = L.nchains;
    v.isdelay.assign(M, false);
    for (std::size_t i = 0; i < M; ++i) {
        v.isdelay[i] = ba_is_delay(L, i);
        if (!v.isdelay[i]) v.qstat.push_back(i);
    }
    v.Nv.resize(C);
    v.Zc.assign(C, num_traits<T>::from_int(0));
    for (std::size_t c = 0; c < C; ++c) {
        v.Nv[c] = num_traits<T>::from_double(v.d.Nchain[c]);
        for (std::size_t i = 0; i < M; ++i)
            if (v.isdelay[i]) v.Zc[c] += v.d.Lchain(i, c);
    }
    return v;
}


/**
 * The families whose bound is a function of the single-chain demand vector
 * D = V/rates, the think time Z and the population N. A multiclass or open model
 * simply does not have those, which is why the rule is total.
 */
inline bool ba_is_single_class_family(const std::string& fam) {
    // 'mapamva' is single-class for a different reason from the rest -- its LP
    // variables QN(i,k)/UN(i,k) are indexed by station and MAP phase, with no
    // class index at all -- but the premise it fails on is the same one.
    return fam == "auto" || fam == "aba" || fam == "bjb" || fam == "pb" || fam == "sb" ||
           fam == "gb" || fam == "harel" || fam == "lr" || fam == "pbh" || fam == "cbh" ||
           fam == "pbk" || fam == "bjbk" || fam == "ssd" || fam == "sib" || fam == "scb" ||
           fam == "ldbcmp" || fam == "mapamva";
}

/**
 * The multiclass families: they take a per-chain demand MATRIX and a population
 * VECTOR, so several classes are fine and an infinite population is not.
 */
inline bool ba_is_fully_closed_family(const std::string& fam) {
    return fam == "mwba" || fam == "cub" || fam == "mbjb" || fam == "looping";
}

/**
 * The single-server families whose alternative on a multiserver model IS 'ssd',
 * which is why the reason names it. 'ssd' is the multiserver bound itself,
 * 'ldbcmp' is parameterized by the limiting demand of a load-dependent station
 * and 'auto' composes whichever candidates survive, so all three are absent.
 */
inline bool ba_is_single_server_family(const std::string& fam) {
    return fam == "aba" || fam == "bjb" || fam == "pb" || fam == "sb" || fam == "gb" ||
           fam == "harel" || fam == "lr" || fam == "pbh" || fam == "cbh" || fam == "pbk" ||
           fam == "bjbk" || fam == "sib" || fam == "scb" || fam == "mapamva";
}

}  // namespace detail

/**
 * The STRUCTURAL premises of the SolverBA bound families, in one place: the
 * reason METHOD cannot bound the model L, or "" when it can.
 *
 * ONE PREDICATE, THREE CALLERS. `solver_ba_analyzer` asks it before dispatching
 * and throws what it returns, `solver_ba_run_analyzer` asks it on the way in, and
 * `list_valid_methods(L)` asks it so the name never reaches a report at all,
 * which is the route `auto_family_methods` takes. A second copy of any rule
 * below is how the report and the run drift apart: the report offers a pair that throws the moment it is run, which is the
 * defect this function exists to remove.
 *
 * WHAT BELONGS HERE AND WHAT DOES NOT. Only the rules the feature registry
 * cannot name. `qn::Feature` has no entry for "one class", for a server count or
 * for a station count, so those are structural and live here. Rules of the form
 * "this family does not accept a delay station" ARE nameable and belong in
 * `qn::ba_feature_set`, which unsets SchedStrategy_INF for the offending method
 * instead: a feature set can refuse a model for HAVING a construct, never for
 * lacking one.
 *
 * METHOD is taken as the caller spells it and resolved through
 * `qn::ba_resolve_method_name`, the copy of `resolve_method` that lives beside the
 * feature set (this header is below the runner and cannot reach back into it),
 * so 'default' is judged as the gb.upper it runs as and the reason names that.
 * The marking-parameterized spnlp and the QRF reduction bounds carry no rule
 * here: the QRF premise is the reducibility test `list_valid_methods(L)` already
 * applies. Of the three OPEN families, 'bpt' and 'bgt' carry none either -- a
 * closed model is refused by their feature set and their analyzers walk the
 * routing matrix for the rest -- while 'snc' carries one, the SERVICE law.
 *
 * WHY THE SNC SERVICE LAW IS HERE AND THE bpt/bgt ONE IS NOT. All three
 * analyzers refuse a non-exponential law at a queueing station. For bpt and bgt
 * that rule extends to the SOURCE and is registry-expressible, so it rides in
 * `qn::ba_feature_set` as a dropped law: both are invariant to the arrival law
 * beyond its mean, so a non-exponential source is not something they refuse, it
 * is something they silently bound as if it were Poisson. snc is the opposite:
 * it CONSUMES the arrival law and its analyzer branches on a non-exponential
 * source deliberately. Its rule is about the SERVICE only, and no feature name
 * can say "Erlang at a Queue but not at a Source", so it is structural.
 *
 * Mirrors `matlab/src/solvers/BA/ba_method_refusal.m` and its JAR and native
 * python twins.
 */
template <class T>
std::string method_refusal(const qn::NetworkStruct<T>& L, const std::string& method) {
    const std::string resolved = qn::ba_resolve_method_name(method);
    const std::string fam = qn::ba_family_of(resolved);
    bool anyOpen = false;
    for (std::size_t r = 0; r < L.classes.size(); ++r)
        if (!std::isfinite(L.classes[r].population)) anyOpen = true;
    const bool anyClosed = L.nclosedjobs() > 0.0;
    if (detail::ba_is_single_class_family(fam)) {
        if (L.nclasses != 1 || !anyClosed)
            return "Method '" + resolved + "' supports single-class closed networks only.";
    } else if (detail::ba_is_fully_closed_family(fam)) {
        if (!anyClosed || anyOpen)
            return "Method '" + resolved + "' supports fully closed networks only.";
    }
    bool multiserver = false;
    for (std::size_t i = 0; i < L.stations.size() && !multiserver; ++i) {
        if (L.stations[i].sched == SchedStrategy::INF) continue;
        if (L.stations[i].nservers > 1.0) multiserver = true;
    }
    if (multiserver) {
        if (detail::ba_is_single_server_family(fam))
            return "Method '" + resolved +
                   "' does not support multi-server stations (use 'ssd').";
        if (detail::ba_is_fully_closed_family(fam))
            return "Method '" + resolved + "' does not support multi-server stations.";
    }

    // The SNC service law. Judged over the pairs a station COULD serve rather
    // than over the ones that carry traffic: the analyzer restricts to the
    // latter, which needs the traffic equations solved, and this is their
    // conservative outer approximation -- it never admits a pair the analyzer
    // refuses, and can only differ on a pair given a service time at a station
    // its class never visits. A Source is skipped, which is the whole reason
    // this is not a feature-set delta.
    if (fam == "snc") {
        for (std::size_t i = 0; i < L.nstations; ++i) {
            if (L.stations[i].nodetype == qn::NodeType::Source) continue;
            for (std::size_t r = 0; r < L.nclasses; ++r) {
                const double mu = num_traits<T>::to_double(L.rates(i, r));
                if (!std::isfinite(mu) || mu <= 0.0) continue;
                if (L.procid(i + 1, r + 1) != lang::ProcessType::EXP)
                    return "Method '" + resolved + "' requires exponential service: station " +
                           std::to_string(i + 1) + " class " + std::to_string(r + 1) +
                           " is not exponential.";
            }
        }
    }
    return "";
}

/**
 * Whether METHOD APPLIES to L but its bound carries no information there, and
 * why. Empty when the bound is informative, and empty for every method that has
 * no such regime.
 *
 * THIS IS A DIFFERENT QUESTION FROM `method_refusal`, which is why it is a
 * different function. That one answers "is this model outside the method's
 * domain", and its answer is what the analyzer throws. This one answers "inside
 * the domain, does the formula still say anything", and its answer is NOT
 * thrown: a degenerate bound is a VALID bound, just a vacuous one, so an
 * analyzer asked for it by name is entitled to publish it -- which is also what
 * `tests/test_ba.cpp` pins at the regime boundary. What must not happen is
 * OFFERING it: `list_valid_methods(L)` names the pairs a caller can act on, and
 * a table of zeros over a network with jobs circulating in it is not one.
 *
 * THE ONE METHOD WITH SUCH A REGIME IS 'ldbcmp.lower'. The Anselmi-Cremonesi
 * bound is built from the population SURPLUS a = N - Qhat, where Qhat is the
 * occupancy the non-bottleneck stations and the think time would hold in the
 * open network fed at the bottleneck's saturation rate. `pfqn_ldbcmp` reports
 * `applicable = false` below the regime, which the analyzer already throws on;
 * AT the boundary a = 0 it returns Xlo = 0, which is formally the trivial bound
 * X >= 0 and propagates into a table whose queue lengths, utilizations and
 * throughputs are all zero. Every entry of that table is a true lower bound and
 * none of them is usable, and a caller cannot tell it from a real answer of
 * zero.
 *
 * Mirrors `matlab/src/solvers/BA/ba_method_degenerate.m` and its JAR and native
 * python twins.
 */
template <class T>
std::string method_degenerate(const qn::NetworkStruct<T>& L, const std::string& method) {
    if (qn::ba_resolve_method_name(method) != "ldbcmp.lower") return "";
    // The applicability rules come first and are not restated: a model this
    // method is outside the domain of has no bound to be degenerate about.
    if (!method_refusal(L, method).empty()) return "";

    const detail::ScDemands<T> sc = detail::ba_sc_demands(L);
    // NO QUEUEING STATION, NO BOUND. Every station is a delay (or a Place, on a
    // Petri net, which is an INF station too), so there is no bottleneck to
    // build Qhat on and `pfqn_ldbcmp` refuses an empty demand vector. A
    // PREDICATE MUST NOT THROW -- this one is asked once per name by
    // `list_valid_methods`, before the Petri sieve has had a chance to drop
    // anything -- so the case is answered rather than propagated.
    if (sc.D.empty())
        return "Method 'ldbcmp.lower' has no queueing station to bound here: every station is "
               "an infinite server, so the bottleneck the open-network occupancy is built on "
               "does not exist.";
    const T zero = num_traits<T>::from_int(0);
    const std::vector<T> cz(sc.D.size(), zero);
    const pfqn::LdBcmpBound<T> b = pfqn::pfqn_ldbcmp(
        sc.D, num_traits<T>::from_int(static_cast<int>(sc.N)), sc.Z, cz,
        num_traits<T>::from_double(1e-10));
    if (!b.applicable)
        return "Method 'ldbcmp.lower' does not apply here: the population is below the "
               "open-network occupancy Qhat the bound is built from.";
    const double xlo = num_traits<T>::to_double(b.Xlo);
    if (!std::isfinite(xlo) || xlo <= 0.0) {
        std::ostringstream os;
        os << "Method 'ldbcmp.lower' needs a population strictly above the open-network "
              "occupancy the bound is built from (Qhat="
           << num_traits<T>::to_double(b.Qhat) << ", N=" << sc.N
           << "): with no surplus it degenerates to the trivial bound X >= 0 and reports a "
              "table of zeros.";
        return os.str();
    }
    return "";
}

/**
 * Port of `solver_ba_analyzer`.
 *
 * @param L      the refreshed struct of a closed model
 * @param opt    the method and, for a hierarchical family, the level
 * @return the [Q,U,R,T,C,X] of the requested bound
 */
template <class T>
BaSolution<T> solver_ba_analyzer(const qn::NetworkStruct<T>& L, const BaOptions& opt) {
    const std::string& method = opt.method;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    // The STRUCTURAL premises of every family -- single-class closed, fully
    // closed, single-server -- are asked here and nowhere else, so this run and
    // the report give one answer whichever a caller meets first. Asking once at
    // the top is what replaced the per-branch copies: an applicability test
    // inlined in each family is how the gate and the run drift apart.
    //
    // IT ALSO RETIRED `ba_out_of_domain`. The five noniterative single-class
    // families and `mwba` used to answer a model outside their CLASS domain
    // with zeroed metrics and an empty C and X, reproducing a reference whose
    // own branches were wrapped in an applicability test with no else. That was
    // a deliberate alignment (user ruling 2026-07-25) and it is OVERTURNED
    // (user ruling 2026-09-04): the reference refuses by name since
    // `ba_method_refusal.m` landed, and a table a caller cannot tell from a
    // real answer of zero is the defect class this gate exists to remove.
    {
        const std::string refusal = method_refusal(L, method);
        if (!refusal.empty()) throw UnsupportedError("solver_ba_analyzer: " + refusal);
    }

    // AUTO composite: evaluate every noniterative bound and keep the tightest
    // side. Feasibility is probed by execution -- a candidate that rejects the
    // model (multiserver, delay station, regime gate) throws and is skipped --
    // so the list stays correct as families are added.
    if (method == "auto.upper" || method == "auto.lower") {
        const bool up = method == "auto.upper";
        const detail::ScDemands<T> sc = detail::ba_sc_demands(L);
        static const char* const cand_up[] = {"aba.upper", "bjb.upper",  "pb.upper",  "gb.upper",
                                              "sb.upper",  "mwba.upper", "ssd.upper", "cub.upper"};
        static const char* const cand_lo[] = {"aba.lower",  "bjb.lower", "pb.lower",
                                              "gb.lower",   "sb.lower",  "mwba.lower",
                                              "ssd.lower",  "mbjb.lower", "ldbcmp.lower"};
        const char* const* cand = up ? cand_up : cand_lo;
        const std::size_t ncand = up ? sizeof(cand_up) / sizeof(*cand_up)
                                     : sizeof(cand_lo) / sizeof(*cand_lo);
        bool have = false;
        T Xbest = zero;
        for (std::size_t ci = 0; ci < ncand; ++ci) {
            BaOptions oc = opt;
            oc.method = cand[ci];
            BaSolution<T> r;
            try {
                r = solver_ba_analyzer(L, oc);
            } catch (const std::exception&) {
                continue;
            }
            if (r.X.empty()) continue;
            const T Xc = r.X[0];
            if (!std::isfinite(num_traits<T>::to_double(Xc)) || !(Xc > zero)) continue;
            if (!have || (up && Xc < Xbest) || (!up && Xc > Xbest)) {
                Xbest = Xc;
                have = true;
            }
        }
        if (!have)
            throw UnsupportedError("solver_ba_analyzer: method '" + method +
                                   "' found no feasible bound for this model");
        return detail::ba_fill(L, sc.V, sc.N, Xbest, sc.Z, sc.D, up);
    }

    // The noniterative single-class families (aba, bjb, pb, sb, gb) share the
    // demand extraction; the reference inlines it in every branch.
    const bool single_class_family =
        method == "aba.upper" || method == "aba.lower" || method == "bjb.upper" ||
        method == "bjb.lower" || method == "pb.upper" || method == "pb.lower" ||
        method == "sb.upper" || method == "sb.lower" || method == "gb.upper" ||
        method == "gb.lower";

    if (single_class_family) {
        const detail::ScDemands<T> sc = detail::ba_sc_demands(L);
        const std::vector<T>& V = sc.V;
        const std::vector<T>& D = sc.D;
        const T& Z = sc.Z;
        const long N = sc.N;
        const T Nt = num_traits<T>::from_int(N);
        const detail::DemandStats<T> st = detail::ba_demand_stats(D);

        if (method == "aba.upper") {
            const T X = detail::ba_min(T(one / st.max), T(Nt / T(Z + st.sum)));
            return detail::ba_fill_xc(L, V, N, X, T(Z + Nt * st.sum), true);
        }
        if (method == "aba.lower") {
            const T X = T(Nt / T(Z + Nt * st.sum));
            return detail::ba_fill_xc(L, V, N, X, T(Z + st.sum), false);
        }
        // The balanced and proportional families evaluate the ABA bracket at
        // N-1: both are one exact MVA step taken from a bounded arrival-theorem
        // queue length, so the population they see is the one a job arriving to
        // the system leaves behind.
        const T Nm1 = num_traits<T>::from_int(N - 1);
        const T xup1 = detail::ba_min(T(one / st.max), T(Nm1 / T(Z + st.sum)));
        const T xlo1 = T(Nm1 / T(Z + Nm1 * st.sum));
        if (method == "bjb.upper") {
            const T C = T(Z + st.sum + st.max * T(Nm1 - Z * xlo1));
            const T X = detail::ba_min(
                T(one / st.max), T(Nt / T(Z + st.sum + st.mean * T(Nm1 - Z * xup1))));
            return detail::ba_fill_xc(L, V, N, X, C, true);
        }
        if (method == "bjb.lower") {
            const T C = T(Z + st.sum + st.mean * T(Nm1 - Z * xup1));
            const T X = T(Nt / T(Z + st.sum + st.max * T(Nm1 - Z * xlo1)));
            return detail::ba_fill_xc(L, V, N, X, C, false);
        }
        if (method == "pb.upper" || method == "pb.lower") {
            T d2 = zero, dN = zero, dNm1 = zero;
            for (const T& d : D) {
                d2 += T(d * d);
                dN += num_pow_int(d, static_cast<unsigned>(N));
                dNm1 += num_pow_int(d, static_cast<unsigned>(N - 1));
            }
            const T Dpb2 = T(d2 / st.sum);
            const T DpbN = T(dN / dNm1);
            if (method == "pb.upper") {
                const T C = T(Z + st.sum + DpbN * T(Nm1 - Z * xlo1));
                const T X = detail::ba_min(
                    T(one / st.max), T(Nt / T(Z + st.sum + Dpb2 * T(Nm1 - Z * xup1))));
                return detail::ba_fill_xc(L, V, N, X, C, true);
            }
            const T C = T(Z + st.sum + Dpb2 * T(Nm1 - Z * xup1));
            const T X = T(Nt / T(Z + st.sum + DpbN * T(Nm1 - Z * xlo1)));
            return detail::ba_fill_xc(L, V, N, X, C, false);
        }
        if (method == "sb.upper" || method == "sb.lower") {
            // The power sums of Harel-Namn-Sturm are written for a network with
            // no terminal population, and the reference refuses a delay station
            // rather than drop its think time.
            for (std::size_t i = 0; i < L.nstations; ++i)
                if (detail::ba_is_delay(L, i))
                    throw UnsupportedError("solver_ba_analyzer: method '" + method +
                                           "' does not support infinite-server stations");
            T A1 = zero, A2 = zero, A3 = zero;
            for (const T& d : D) {
                A1 += d;
                A2 += T(d * d);
                A3 += T(d * d * d);
            }
            if (method == "sb.upper") {
                const T C = T(Z + A1 + Nm1 * T(T(A1 * A2 + A3) / T(A1 * A1 + A2)));
                const T X = detail::ba_min(T(one / st.max), T(Nt / C));
                // The power-sum family leaves the per-station response time
                // undefined, and the reference pairs BOTH sides with the
                // no-contention residence rather than only the lower one.
                return detail::ba_fill_xc(L, V, N, X, C, false);
            }
            // The (N-1)-st root of the N-th power sum has no field expression.
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "solver_ba_analyzer: method 'sb.lower' needs an (N-1)-st root and is "
                    "unavailable in exact arithmetic");
            } else {
                using std::pow;
                T AN = zero;
                for (const T& d : D) AN += num_pow_int(d, static_cast<unsigned>(N));
                // pow is taken on T, not through double: a high-precision
                // backend would otherwise lose every digit past the 17th here.
                const T root = pow(T(AN / A1), T(one / Nm1));
                const T C = T(Z + A1 + Nm1 * root);
                return detail::ba_fill_xc(L, V, N, T(Nt / C), C, false);
            }
        }
        // gb: the geometric bounds solve a quadratic in the throughput.
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError("solver_ba_analyzer: the geometric bounds '" + method +
                                   "' solve a quadratic and are unavailable in exact arithmetic");
        } else {
            const bool up = (method == "gb.upper");
            const T Xup = pfqn::pfqn_xzgsbup(D, Nt, Z);
            const T Xlo = pfqn::pfqn_xzgsblow(D, Nt, Z);
            const T X = up ? detail::ba_min(T(one / st.max), Xup) : Xlo;
            const T C = T(Nt / (up ? Xlo : Xup));
            // The queue-length bound is the geometric one, not the residence
            // convention of ba_fill; the response time follows from it by the
            // OPPOSITE side of the throughput bracket, which is what keeps
            // Q = X R consistent with a bound rather than with a point estimate.
            const T Xden = up ? Xlo : Xup;
            BaSolution<T> s;
            s.Q = Matrix<T>(L.nstations, 1, zero);
            s.U = Matrix<T>(L.nstations, 1, zero);
            s.R = Matrix<T>(L.nstations, 1, zero);
            s.Tp = Matrix<T>(L.nstations, 1, zero);
            std::size_t k = 0;
            for (std::size_t i = 0; i < L.nstations; ++i) {
                s.Tp(i, 0) = T(V[i] * X);
                if (detail::ba_is_delay(L, i)) {
                    s.R(i, 0) = T(one / L.rates(i, 0));
                    s.Q(i, 0) = T(X * s.R(i, 0));
                } else {
                    s.Q(i, 0) = up ? pfqn::pfqn_qzgbup(D, Nt, Z, k) : pfqn::pfqn_qzgblow(D, Nt, Z, k);
                    s.R(i, 0) = T(T(s.Q(i, 0) / Xden) / V[i]);
                    ++k;
                }
                s.U(i, 0) = detail::ba_is_delay(L, i) ? s.Q(i, 0) : T(s.Tp(i, 0) / L.rates(i, 0));
            }
            s.C.assign(1, C);
            s.X.assign(1, X);
            s.lG = -static_cast<double>(N) * num_traits<T>::log_as_double(X);
            return s;
        }
    }

    if (method == "lr.upper" || method == "lr.lower") {
        // LP linear-reduction bound: one LP per station, each minimizing or
        // maximizing that station's utilization over a polytope that CONTAINS
        // the exact stationary solution.
        const detail::ScDemands<T> sc = detail::ba_sc_demands(L);
        for (std::size_t i = 0; i < L.nstations; ++i)
            if (detail::ba_is_delay(L, i))
                throw UnsupportedError("solver_ba_analyzer: method '" + method +
                                       "' does not support delay (infinite-server) stations");
        const bool up = (method == "lr.upper");
        const std::size_t M = L.nstations;
        const std::vector<T>& V = sc.V;

        mapqn::LrPfParams<T> par;
        par.M = static_cast<int>(M);
        // MATLAB reads sn.njobs(1), the class population, which for the
        // single-class model this branch admits is nclosedjobs.
        par.N = static_cast<int>(sc.N);
        par.mu.resize(M);
        T totV = zero;
        for (std::size_t i = 0; i < M; ++i) {
            par.mu[i] = L.rates(i, 0);
            totV += V[i];
        }
        // Every row of r is the visit vector normalized by its total: the
        // reference builds it with repmat, so a job leaves any station for j
        // with the same probability V_j / sum V.
        par.r = Matrix<T>(M, M, zero);
        if (totV == zero) throw NumericError("solver_ba_analyzer: the visit ratios are all zero");
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < M; ++j) par.r(i, j) = T(V[j] / totV);

        const mapqn::MapqnSense sense = up ? mapqn::MapqnSense::Max : mapqn::MapqnSense::Min;
        std::vector<T> U(M, zero);
        for (std::size_t ti = 0; ti < M; ++ti) {
            const mapqn::LrPfResult<T> r =
                mapqn::mapqn_bnd_lr_pf(par, static_cast<int>(ti) + 1, sense);
            if (!r.ok)
                throw NumericError("solver_ba_analyzer: the '" + method +
                                   "' linear program did not solve (" + r.status + ")");
            U[ti] = r.objective;
        }
        // The chain throughput implied by the bounded utilizations, U_i = X V_i / mu_i.
        // MATLAB drops the non-finite candidates (a station with no visit) and
        // takes the minimum, or zero when none survives.
        bool any = false;
        T X = zero;
        for (std::size_t i = 0; i < M; ++i) {
            if (V[i] == zero) continue;  // MATLAB: the candidate is Inf and isfinite drops it
            const T cand = T(U[i] * L.rates(i, 0) / V[i]);
            if (!any || cand < X) {
                X = cand;
                any = true;
            }
        }
        if (!any) X = zero;
        // Q, R, T and C follow the ABA residence convention of ba_fill at this
        // side, with no think time since a delay station is refused above; only
        // U comes from the LP rather than from the utilization law.
        BaSolution<T> s = detail::ba_fill(L, V, sc.N, X, zero, sc.D, up);
        for (std::size_t i = 0; i < M; ++i) s.U(i, 0) = U[i];
        s.lG = std::numeric_limits<double>::quiet_NaN();  // the reference leaves it unset
        return s;
    }

    if (method == "mapamva.upper" || method == "mapamva.lower") {
        // MAP-AMVA (Casale-Smirni, DSN 2009): the LP over the EXACT mean-value
        // balance equations of a closed MAP queueing network. It is the only
        // family here that consumes the CORRELATION between successive services
        // rather than the service mean alone -- its variables are the per-phase
        // queue lengths QN(i,k) and utilizations UN(i,k), so a workload whose
        // burstiness moves the bottleneck between stations is bounded rather
        // than averaged into a renewal process. That is why 'MAP' and 'MMPP2'
        // reach this family's feature set and no other.
        //
        // The LP carries phases at ONE queue and requires it to be the LAST --
        // q(i,j,k,h) reads the scalar muM(i) for i < M and the (D0,D1) pair
        // muMAP/v for i == M -- so a model whose phase-carrying station sits
        // elsewhere is PERMUTED rather than refused, and the results are
        // permuted back before they are returned.
        const detail::ScDemands<T> sc = detail::ba_sc_demands(L);
        for (std::size_t i = 0; i < L.nstations; ++i)
            if (detail::ba_is_delay(L, i))
                throw UnsupportedError(
                    "solver_ba_analyzer: method '" + method +
                    "' does not support delay (infinite-server) stations: the MAP-AMVA program "
                    "of Casale-Smirni (DSN 2009) is written for a network of queues and the "
                    "paper names the delay extension as open work. Use a QRF method, which "
                    "carries the load-dependent rate law");
        const bool up = (method == "mapamva.upper");
        const std::size_t M = L.nstations;
        const std::vector<T>& V = sc.V;

        // Phase order per station. One phase is an exponential server, which
        // enters the LP as the scalar rate muM(i); more than one is the (D0,D1)
        // pair, which only queue M can hold.
        std::vector<std::pair<Matrix<T>, Matrix<T>>> MAPs(M);
        std::vector<int> kph(M, 1);
        for (std::size_t i = 0; i < M; ++i) {
            const mam::Map<T> m = lang::dist_to_map(L.service[i][0]);
            if (m.D0.rows() == 0) {
                MAPs[i] = std::make_pair(Matrix<T>(1, 1, num_traits<T>::from_int(-1)),
                                         Matrix<T>(1, 1, num_traits<T>::from_int(1)));
                kph[i] = 1;
            } else {
                MAPs[i] = std::make_pair(m.D0, m.D1);
                kph[i] = static_cast<int>(m.D0.rows());
            }
        }
        std::vector<std::size_t> phased;
        for (std::size_t i = 0; i < M; ++i)
            if (kph[i] > 1) phased.push_back(i);
        if (phased.size() > 1) {
            std::string names;
            for (std::size_t j = 0; j < phased.size(); ++j)
                names += (j ? ", " : "") + std::to_string(phased[j] + 1);
            throw UnsupportedError(
                "solver_ba_analyzer: method '" + method +
                "' carries phases at ONE station: the LP gives queue M the (D0,D1) pair and "
                "every other queue a scalar rate. Stations " + names + " are all "
                "non-exponential. Use a QRF method, whose q carries a phase at every station");
        }
        // Every station exponential: the program is still the right one, it just
        // degenerates to K = 1, where the per-phase variables collapse and the
        // balances become the product-form ones of mapqn_bnd_lr_pf.
        const std::size_t map_idx = phased.empty() ? M - 1 : phased.front();
        std::vector<std::size_t> perm;
        perm.reserve(M);
        for (std::size_t i = 0; i < M; ++i)
            if (i != map_idx) perm.push_back(i);
        perm.push_back(map_idx);
        const int K = kph[map_idx];
        const std::size_t Ks = static_cast<std::size_t>(K);

        mapqn::LrMvaParams<T> par;
        par.M = static_cast<int>(M);
        par.N = static_cast<int>(sc.N);
        par.K = K;
        par.muM.resize(M - 1);
        for (std::size_t a = 0; a + 1 < M; ++a) par.muM[a] = L.rates(perm[a], 0);
        // muMAP(k,h) is the completion rate out of phase k landing in phase h,
        // i.e. D1(k,h); v(k,h) is the background phase change that completes no
        // job, i.e. D0 off the diagonal. Same (from,to) convention as
        // qrf_extract_mu_v -- writing either as its transpose is invisible for a
        // reversible D0 and silently reverses the phase order of an Erlang.
        par.muMAP = Matrix<T>(Ks, Ks, zero);
        par.v = Matrix<T>(Ks, Ks, zero);
        for (std::size_t a = 0; a < Ks; ++a)
            for (std::size_t b = 0; b < Ks; ++b) {
                par.muMAP(a, b) = MAPs[map_idx].second(a, b);
                par.v(a, b) = (a == b) ? zero : MAPs[map_idx].first(a, b);
            }
        par.r = Matrix<T>(M, M, zero);
        for (std::size_t a = 0; a < M; ++a)
            for (std::size_t b = 0; b < M; ++b) par.r(a, b) = L.rt(perm[a], perm[b]);

        std::vector<T> Vp(M, zero), Sp(M, zero);
        for (std::size_t a = 0; a < M; ++a) {
            const std::size_t i = perm[a];
            Vp[a] = V[i];
            mam::Map<T> mi;
            mi.D0 = MAPs[i].first;
            mi.D1 = MAPs[i].second;
            Sp[a] = mam::map_mean(mi);
        }

        // THREE SWEEPS OF THE SAME LP, and the utilization one runs in BOTH
        // senses on purpose. R_i = Q_i/(V_i*X) rises with Q_i and FALLS with X,
        // so an upper bound on the response time pairs Q_i^max with X^min;
        // dividing by X^max on both sides is what would report an upper R below
        // the exact value and break the bracket.
        const mapqn::MapqnSense sense = up ? mapqn::MapqnSense::Max : mapqn::MapqnSense::Min;
        std::vector<T> umax(M, zero), umin(M, zero), qbnd(M, zero);
        for (std::size_t a = 0; a < M; ++a) {
            const int ti = static_cast<int>(a);
            const mapqn::MapqnBndLrMvaResult<T> ru =
                mapqn::mapqn_bnd_lr_mva(par, ti, -1, mapqn::MapqnSense::Max,
                                        mapqn::MapqnObjectiveVar::UN);
            const mapqn::MapqnBndLrMvaResult<T> rl =
                mapqn::mapqn_bnd_lr_mva(par, ti, -1, mapqn::MapqnSense::Min,
                                        mapqn::MapqnObjectiveVar::UN);
            const mapqn::MapqnBndLrMvaResult<T> rq =
                mapqn::mapqn_bnd_lr_mva(par, ti, -1, sense, mapqn::MapqnObjectiveVar::QN);
            if (!ru.ok || !rl.ok || !rq.ok)
                throw NumericError("solver_ba_analyzer: the '" + method +
                                   "' linear program did not solve (" +
                                   (ru.ok ? (rl.ok ? rq.status : rl.status) : ru.status) + ")");
            umax[a] = ru.objective;
            umin[a] = rl.objective;
            qbnd[a] = rq.objective;
        }

        // Utilization law U_i = X*V_i*S_i, exact at a single server under ANY
        // service law, so each station turns its own utilization bound into a
        // throughput bound and the tightest of the M survives. A station with no
        // visits or no service time carries no information and is skipped.
        bool any = false;
        T x_up = zero, x_lo = zero;
        for (std::size_t a = 0; a < M; ++a) {
            const T load = T(Vp[a] * Sp[a]);
            if (load <= zero) continue;
            const T cu = T(umax[a] / load), cl = T(umin[a] / load);
            if (!any) {
                x_up = cu;
                x_lo = cl;
                any = true;
            } else {
                if (cu < x_up) x_up = cu;
                if (cl > x_lo) x_lo = cl;
            }
        }
        if (!any)
            throw NumericError("solver_ba_analyzer: method '" + method +
                               "' found no station with both a positive visit ratio and a "
                               "positive mean service time");
        const T xb = up ? x_up : x_lo;
        const T xopp = up ? x_lo : x_up;

        // Unpermute: the LP orders the stations with the phase-carrying one last.
        BaSolution<T> s;
        s.Q = Matrix<T>(M, 1, zero);
        s.U = Matrix<T>(M, 1, zero);
        s.R = Matrix<T>(M, 1, zero);
        s.Tp = Matrix<T>(M, 1, zero);
        for (std::size_t a = 0; a < M; ++a) {
            const std::size_t i = perm[a];
            s.U(i, 0) = up ? umax[a] : umin[a];
            s.Q(i, 0) = qbnd[a];
            s.Tp(i, 0) = T(Vp[a] * xb);
            s.R(i, 0) = (xopp > zero && Vp[a] > zero) ? T(qbnd[a] / (Vp[a] * xopp)) : zero;
        }
        // Delay stations are refused above, so the closed-network response time
        // is N/X exactly and the throughput bracket transfers to it directly.
        // Summing the per-station R bounds instead would add M separately
        // attained maxima and report a looser number.
        s.C.assign(1, xopp > zero
                          ? T(num_traits<T>::from_int(static_cast<long>(sc.N)) / xopp)
                          : zero);
        s.X.assign(1, xb);
        s.lG = -static_cast<double>(sc.N) * num_traits<T>::log_as_double(xb);
        return s;
    }

    if (method == "mwba.upper" || method == "mwba.lower") {
        // Majumdar-Woodside robust box bounds: multiclass, and the only family
        // that reads the scheduling discipline, since its lower bound is a
        // per-discipline residence guarantee.
        const detail::ChainView<T> v = detail::ba_chain_view(L);
        const std::size_t M = L.nstations, C = L.nchains, Kq = v.qstat.size();
        Matrix<T> Vq(Kq, C, zero), Sq(Kq, C, zero);
        std::vector<pfqn::MwrbbSched> schedq(Kq);
        for (std::size_t k = 0; k < Kq; ++k) {
            const std::size_t i = v.qstat[k];
            schedq[k] = detail::mwrbb_disc_code(L.stations[i].sched);
            for (std::size_t c = 0; c < C; ++c) {
                Vq(k, c) = v.d.Vchain(i, c);
                Sq(k, c) = v.d.STchain(i, c);
            }
        }
        // Chain priority: the reference class's when the chain has one, else the
        // highest priority (lowest value) any class of the chain carries.
        std::vector<int> prioc(C, 0);
        for (std::size_t c = 0; c < C; ++c) {
            if (L.refclass[c] > 0) {
                prioc[c] = L.classes[L.refclass[c] - 1].prio;
            } else {
                int p = L.classes[L.inchain[c][0] - 1].prio;
                for (std::size_t r : L.inchain[c]) p = std::min(p, L.classes[r - 1].prio);
                prioc[c] = p;
            }
        }
        const pfqn::MwrbbBounds<T> b = pfqn::pfqn_mwrbb(Vq, Sq, v.Nv, v.Zc, schedq, prioc);
        const bool up = (method == "mwba.upper");
        const std::vector<T>& Xchain = up ? b.Xup : b.Xlo;
        Matrix<T> Tchain(M, C, zero), Uchain(M, C, zero), Qchain(M, C, zero);
        for (std::size_t c = 0; c < C; ++c)
            for (std::size_t i = 0; i < M; ++i) {
                Tchain(i, c) = T(Xchain[c] * v.d.Vchain(i, c));
                Uchain(i, c) = T(Xchain[c] * v.d.Lchain(i, c));  // utilization law
            }
        Qchain = detail::ba_chain_qfill(Uchain, Xchain, v.d.Lchain, v.Nv, v.isdelay, up);
        const mva::ClassResults<T> cr = mva::sn_deaggregate_chain_results(
            L, v.d, Qchain, Uchain, Matrix<T>(), Tchain, Xchain);
        BaSolution<T> s;
        s.Q = cr.Q;
        s.U = cr.U;
        s.R = cr.R;
        s.Tp = cr.Tp;
        s.C = cr.C;
        s.X = cr.X;
        s.lG = std::numeric_limits<double>::quiet_NaN();
        return s;
    }

    // Eager Looping: the multiclass bracket that initializes the multiple-class
    // PBH. Pessimistic side from the heap-inflated response time, optimistic
    // side from the response-time lower bound.
    if (method == "looping.upper" || method == "looping.lower") {
        const detail::ChainView<T> v = detail::ba_chain_view(L);
        const std::size_t M = L.nstations, C = L.nchains, Kq = v.qstat.size();
        Matrix<T> Lq(Kq, C, zero);
        for (std::size_t k = 0; k < Kq; ++k)
            for (std::size_t c = 0; c < C; ++c) Lq(k, c) = v.d.Lchain(v.qstat[k], c);
        const pfqn::LoopingBounds<T> b = pfqn::pfqn_looping(Lq, v.Nv, v.Zc);
        const bool up = (method == "looping.upper");
        const std::vector<T>& Xchain = up ? b.Xup : b.Xlo;
        Matrix<T> Tchain(M, C, zero), Uchain(M, C, zero), Qchain(M, C, zero);
        for (std::size_t c = 0; c < C; ++c)
            for (std::size_t i = 0; i < M; ++i) {
                Tchain(i, c) = T(Xchain[c] * v.d.Vchain(i, c));
                Uchain(i, c) = T(Xchain[c] * v.d.Lchain(i, c));  // utilization law
            }
        Qchain = detail::ba_chain_qfill(Uchain, Xchain, v.d.Lchain, v.Nv, v.isdelay, up);
        const mva::ClassResults<T> cr = mva::sn_deaggregate_chain_results(
            L, v.d, Qchain, Uchain, Matrix<T>(), Tchain, Xchain);
        BaSolution<T> s;
        s.Q = cr.Q;
        s.U = cr.U;
        s.R = cr.R;
        s.Tp = cr.Tp;
        s.C = cr.C;
        s.X = cr.X;
        s.lG = std::numeric_limits<double>::quiet_NaN();
        return s;
    }

    // Achievable-region LP relaxation (Bertsimas-Paschalidis-Tsitsiklis 1994).
    // The only OPEN-network family here: it lower bounds the mean response
    // times attainable by ANY non-idling policy, so it is refused on the closed
    // models every other family requires.
    if (method == "bpt.lower") {
        BaSolution<T> s;
        solver_ba_bpt(L, s);
        return s;
    }

    // Piecewise-linear Lyapunov bound (Bertsimas-Gamarnik-Tsitsiklis 2001), the
    // OPEN-network upper side. A feasible gamma > 0 both certifies that every
    // work-conserving policy is stable and yields the (loose) queue-length bound.
    if (method == "bgt.upper") {
        BaSolution<T> s;
        solver_ba_bgt(L, s);
        return s;
    }

    // Stochastic network calculus (Fidler-Rizk 2015), the third OPEN-network
    // family and the only one whose native object is a TAIL: MGF envelopes are
    // propagated hop by hop and the delay bound is integrated into a mean.
    if (method == "snc.upper") {
        BaSolution<T> s;
        solver_ba_snc(L, s);
        return s;
    }

    if (method == "cub.upper" || method == "mbjb.lower") {
        const detail::ChainView<T> v = detail::ba_chain_view(L);
        const std::size_t M = L.nstations, C = L.nchains, Kq = v.qstat.size();
        Matrix<T> Lq(Kq, C, zero);
        for (std::size_t k = 0; k < Kq; ++k)
            for (std::size_t c = 0; c < C; ++c) Lq(k, c) = v.d.Lchain(v.qstat[k], c);
        const pfqn::McubBounds<T> b = pfqn::pfqn_mcub(Lq, v.Nv, v.Zc);
        const bool up = (method == "cub.upper");
        const std::vector<T>& Xchain = up ? b.Xub : b.Xlb;
        Matrix<T> Tchain(M, C, zero), Uchain(M, C, zero), Qchain(M, C, zero);
        for (std::size_t c = 0; c < C; ++c)
            for (std::size_t i = 0; i < M; ++i) {
                Tchain(i, c) = T(Xchain[c] * v.d.Vchain(i, c));
                Uchain(i, c) = T(Xchain[c] * v.d.Lchain(i, c));  // utilization law
            }
        Qchain = detail::ba_chain_qfill(Uchain, Xchain, v.d.Lchain, v.Nv, v.isdelay, up);
        const mva::ClassResults<T> cr = mva::sn_deaggregate_chain_results(
            L, v.d, Qchain, Uchain, Matrix<T>(), Tchain, Xchain);
        BaSolution<T> s;
        s.Q = cr.Q;
        s.U = cr.U;
        s.R = cr.R;
        s.Tp = cr.Tp;
        s.C = cr.C;
        s.X = cr.X;
        s.lG = std::numeric_limits<double>::quiet_NaN();
        return s;
    }

    // Sharp bounds of Harel-Namn-Sturm, distinct from the 'sb' family of the
    // same paper: they extrapolate from the EXACT normalizing constant at
    // populations n <= 7 instead of using the first three power sums.
    if (method == "harel.upper" || method == "harel.lower") {
        const detail::ScDemands<T> sch = detail::ba_sc_demands(L);
        if (!(sch.Z == zero))
            throw UnsupportedError("solver_ba_analyzer: method '" + method +
                                   "' does not support think times (infinite-server stations)");
        const int Nh = static_cast<int>(sch.N);
        const int maxUB = Nh < 7 ? Nh : 7;
        const pfqn::HarelBoundsResult<T> b = pfqn::pfqn_harel_bounds(sch.D, Nh, zero, maxUB);
        const bool uph = method == "harel.upper";
        T Xb = b.LB;
        if (uph) {
            T dmax = sch.D[0];
            for (const T& d : sch.D)
                if (d > dmax) dmax = d;
            const T cap = T(num_traits<T>::from_int(1) / dmax);
            const T ext = maxUB >= 2 ? b.UB[static_cast<std::size_t>(maxUB)] : b.TH[1];
            Xb = ext < cap ? ext : cap;
        }
        return detail::ba_fill(L, sch.V, sch.N, Xb, sch.Z, sch.D, uph);
    }

    // The hierarchical / iterative single-class families.
    const bool hier = method == "pbh.upper" || method == "pbh.lower" || method == "cbh.upper" ||
                      method == "cbh.lower" || method == "pbk.upper" || method == "pbk.lower" ||
                      method == "bjbk.upper" || method == "bjbk.lower" || method == "ssd.upper" ||
                      method == "ssd.lower" || method == "sib.upper" || method == "sib.lower" ||
                      method == "scb.upper" || method == "scb.lower" ||
                      method == "ldbcmp.lower";
    if (!hier)
        throw UnsupportedError("solver_ba_analyzer: unknown bound method '" + method + "'");

    const detail::ScDemands<T> sc = detail::ba_sc_demands(L);
    const bool up = method.size() > 6 && method.compare(method.size() - 6, 6, ".upper") == 0;
    const T Nt = num_traits<T>::from_int(sc.N);
    const int lvl = opt.level;

    // ssd and ldbcmp are the two families the reference does NOT gate on server
    // count: ssd bounds a multiserver station by disaggregating it into single
    // servers, and ldbcmp reads only the limiting demand.
    if (method == "ssd.upper" || method == "ssd.lower") {
        std::vector<T> cvec;
        for (std::size_t i = 0; i < L.nstations; ++i)
            if (!detail::ba_is_delay(L, i))
                cvec.push_back(num_traits<T>::from_double(L.stations[i].nservers));
        const pfqn::SsdBounds<T> b = pfqn::pfqn_ssd(sc.D, Nt, sc.Z, cvec);
        return detail::ba_fill(L, sc.V, sc.N, up ? b.Xhi : b.Xlo, sc.Z, sc.D, up);
    }
    if (method == "ldbcmp.lower") {
        // Fixed-rate parameterization (Heffes c = 0); the bound holds only in
        // the asymptotic regime N >= Qhat, and the reference errors below it
        // rather than return the NaN pfqn_ldbcmp produces there.
        const std::vector<T> cz(sc.D.size(), zero);
        const pfqn::LdBcmpBound<T> b =
            pfqn::pfqn_ldbcmp(sc.D, Nt, sc.Z, cz, num_traits<T>::from_double(1e-10));
        if (!b.applicable)
            throw UnsupportedError(
                "solver_ba_analyzer: method 'ldbcmp.lower' requires the asymptotic regime N >= "
                "Qhat");
        return detail::ba_fill(L, sc.V, sc.N, b.Xlo, sc.Z, sc.D, false);
    }
    if (method == "scb.upper" || method == "scb.lower") {
        // Single-class bounds of Dowdy et al. (1992). THE BRACKETED OBJECT IS NOT
        // THIS MODEL: scb brackets the multiclass system that this single-class
        // model aggregates, so scb.lower is the EXACT single-class throughput and
        // scb.upper adds the demand-free Expression-(3) gap. That is why scb is
        // absent from the auto candidate list -- mixing it with families that
        // bracket this model's own solution would compare two different quantities.
        if (!(sc.Z == zero))
            throw UnsupportedError("solver_ba_analyzer: method '" + method +
                                   "' supports Z=0 (no delay station) only; Theorem 3 rests on "
                                   "the delay-free balanced-network throughput");
        const pfqn::ScbBounds<T> b = pfqn::pfqn_scb(sc.D, sc.N);
        return detail::ba_fill(L, sc.V, sc.N, up ? b.Xhi : b.Xlo, sc.Z, sc.D, up);
    }
    if (method == "pbh.upper" || method == "pbh.lower") {
        const pfqn::PbhBounds<T> b = pfqn::pfqn_pbh(sc.D, static_cast<int>(sc.N), sc.Z, lvl);
        return detail::ba_fill(L, sc.V, sc.N, up ? b.Xhi : b.Xlo, sc.Z, sc.D, up);
    }
    if (method == "pbk.upper" || method == "pbk.lower") {
        const pfqn::PbhBounds<T> b = pfqn::pfqn_pbk(sc.D, static_cast<int>(sc.N), sc.Z, lvl);
        return detail::ba_fill(L, sc.V, sc.N, up ? b.Xhi : b.Xlo, sc.Z, sc.D, up);
    }
    if (method == "bjbk.upper" || method == "bjbk.lower") {
        const pfqn::PbhBounds<T> b = pfqn::pfqn_bjbk(sc.D, static_cast<int>(sc.N), sc.Z, lvl);
        return detail::ba_fill(L, sc.V, sc.N, up ? b.Xhi : b.Xlo, sc.Z, sc.D, up);
    }
    if (method == "cbh.upper" || method == "cbh.lower") {
        const pfqn::CbhBounds<T> b = pfqn::pfqn_cbh(sc.D, static_cast<int>(sc.N), sc.Z, lvl);
        return detail::ba_fill(L, sc.V, sc.N, up ? b.Xhi : b.Xlo, sc.Z, sc.D, up);
    }
    // sib
    if (!(sc.Z == zero))
        throw UnsupportedError("solver_ba_analyzer: method '" + method +
                               "' supports Z=0 (no delay station) only; delay needs the SIB "
                               "Section-3.2 extension");
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError("solver_ba_analyzer: method '" + method +
                               "' solves a quadratic and is unavailable in exact arithmetic");
    } else {
        const pfqn::SibBounds<T> b = pfqn::pfqn_sib(sc.D, static_cast<int>(sc.N), zero, lvl);
        return detail::ba_fill(L, sc.V, sc.N, up ? b.Xhi : b.Xlo, sc.Z, sc.D, up);
    }
}

}  // namespace ba
}  // namespace line

#endif  // LINE_SOLVERS_BA_SOLVER_BA_ANALYZER_H
