/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_MVA_DISPATCH_H
#define LINE_SOLVERS_MVA_MVA_DISPATCH_H

/**
 * Port of `@@SolverMVA/mvaDispatch.m`: one inner solve, choosing the analyzer
 * that fits the model.
 *
 * THE ORDER IS THE CONTRACT. The branches are NOT disjoint -- a single-class
 * open Source-Queue-Sink model with a size-based discipline satisfies two of
 * them, a cache model with load dependence satisfies two more -- so the first
 * match wins and reordering silently changes which algorithm a model gets.
 * The sequence below is the reference's, top to bottom:
 *
 *   0  closed models with a shortest-job-next station <- ported
 *   1  order-independent / pass-and-swap stations
 *   2  delayed-hit retrieval caches, open then closed
 *   3  size-based scheduling in an open Source-Queue-Sink model
 *   4  single-class open Source-Queue-Sink            <- ported
 *   5  multiclass open polling
 *   6  multiclass open HOL priority                   <- ported
 *   7  multiclass open DPS, at most three classes     <- ported
 *   8  non-reentrant cache (Source-Cache-Sink)
 *   9  integrated caching-queueing
 *  10  bound methods (moved to SolverBA in the reference)
 *  11  Marie's aggregation-decomposition
 *  12  load- or class-dependent scaling               <- ported
 *  13  everything else                                <- solver_mva_analyzer,
 *      of whose method switch 'mvac' is dispatched here (solver_mvac.h)
 *
 * WHAT IS NOT PORTED IS REFUSED BY NAME, never allowed to fall through to the
 * generic analyzer. A cache model solved as an ordinary queueing network
 * returns numbers -- they are simply not the model's -- and the same holds for
 * a polling system and for an OI station, whose rate function the AMVA cannot
 * represent at all.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "line/api/qsys/qsys_gg1.h"
#include "line/api/qsys/qsys_ggnm_diffusion.h"
#include "line/api/qsys/qsys_gig1_bnds_extremal.h"
#include "line/api/qsys/qsys_gig1_approx_allencunneen.h"
#include "line/api/qsys/qsys_gig1_approx_gelenbe.h"
#include "line/api/qsys/qsys_gig1_approx_heyman.h"
#include "line/api/qsys/qsys_gig1_approx_kimura.h"
#include "line/api/qsys/qsys_gig1_approx_klb.h"
#include "line/api/qsys/qsys_gig1_approx_kobayashi.h"
#include "line/api/qsys/qsys_gig1_approx_marchal.h"
#include "line/api/qsys/qsys_gig1_rq.h"
#include "line/api/qsys/qsys_gig1_ubnd_kingman.h"
#include "line/api/qsys/qsys_gigk_approx.h"
#include "line/api/qsys/qsys_gigk_approx_kingman.h"
#include "line/api/qsys/qsys_gigk_approx_whitt.h"
#include "line/api/qsys/qsys_gigk_rqt.h"
#include "line/api/qsys/qsys_gigk_rqt_gamma.h"
#include "line/api/qsys/qsys_gm1.h"
#include "line/api/qsys/qsys_mg1.h"
#include "line/api/qsys/qsys_mg1k_loss_mgs.h"
#include "line/api/qsys/qsys_mgisrgi_whitt.h"
#include "line/api/qsys/qsys_mmk_qed.h"
#include "line/api/qsys/qsys_mg1_fb.h"
#include "line/api/qsys/qsys_mg1_lrpt.h"
#include "line/api/qsys/qsys_mg1_prio.h"
#include "line/api/qsys/qsys_mg1_psjf.h"
#include "line/api/qsys/qsys_mg1_setf.h"
#include "line/api/qsys/qsys_mg1_srpt.h"
#include "line/api/qsys/qsys_mm1.h"
#include "line/api/qsys/qsys_mm1_dps.h"
#include "line/api/qsys/qsys_mmk.h"
#include "line/api/qsys/qsys_phm1.h"
#include "line/api/mam/map_count_idc.h"
#include "line/api/pfqn/pfqn_marie.h"
#include "line/lang/distribution.h"
#include "line/util/rootfind.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/cache_metrics.h"
#include "line/solvers/mva/solver_mva.h"
#include "line/solvers/mva/solver_mva_cache.h"
#include "line/solvers/mva/solver_mva_oi.h"
#include "line/solvers/mva/solver_mva_cacheqn.h"
#include "line/solvers/mva/solver_mva_cacheqn_retrieval.h"
#include "line/solvers/mva/solver_mva_polling.h"
#include "line/solvers/mva/solver_mva_retrieval.h"
#include "line/solvers/mva/solver_mva_sjn.h"
#include "line/solvers/mva/solver_mvac.h"
#include "line/solvers/mva/solver_qna.h"
#include "line/solvers/mva/solver_rqna.h"
#include "line/solvers/mva/solver_rqt.h"
#include "line/solvers/mva/solver_mapqn.h"
#include "line/api/sn/sn_has_bursty_arrival.h"
#include "line/api/sn/sn_patience_handles.h"

namespace line {
namespace mva {

/** What the dispatch returns: the metrics plus the algorithm that produced them. */
template <class T>
struct DispatchResult {
    MvaSolution<T> sol;
    /** The concrete algorithm, as the reference's `actualmethod`. */
    std::string actualmethod;
    /**
     * Set only by the integrated cacheqn branch: the converged struct whose
     * routing carries the actual hit/miss probabilities, from which the runner
     * derives ArvR and ResidT. Empty for every other model.
     */
    std::shared_ptr<qn::NetworkStruct<T>> refreshed_struct;
    /**
     * What a cache branch measured, EMPTY on every model without a Cache.
     *
     * A cache's hit / miss / delayed-hit split is a SOLVER RESULT, not model
     * state, and it is the only part of the answer the (station x class) tables
     * cannot carry. Dropping it here is what left `-a cache` refusing under
     * `-s mva` and, through CPPLINE, left MATLAB's Cache node holding whatever
     * the PREVIOUS solver wrote -- retrieval_simple printed the LDES table under
     * the MVA label. Assembled by `solvers::cache_metrics_of`, one rule for
     * every branch of every solver.
     */
    solvers::CacheMetrics<T> cache;
    /**
     * Non-empty when the analyzer that ran would have raised a reference
     * `line_warning` and still returned a usable answer, carrying its text.
     * Set only by the SJN branch today (solver_mva_sjn.h).
     *
     * KNOWN DIVERGENCE, DELIBERATE. This port has no general warning channel --
     * no line_warning equivalent exists under cpp/include/line/ -- so the text
     * stops here, at the dispatch boundary. A caller of `mva_dispatch` sees it;
     * a user going through SolverMVA does NOT, because solver_mva_runner.h does
     * not read the field. Closing that last hop needs a runner change and a
     * decision about where a solver-level warning should surface at all, which
     * is wider than this analyzer.
     */
    std::string warning;
};

namespace detail {

/** True when the model is exactly a Source, a Queue and a Sink. */
template <class T>
bool is_open_sqs(const qn::NetworkStruct<T>& L) {
    if (L.nof_nodes() != 3) return false;
    int src = 0, q = 0, snk = 0;
    for (const qn::NodeDef& nd : L.nodes) {
        if (nd.nodetype == qn::NodeType::Source) ++src;
        else if (nd.nodetype == qn::NodeType::Queue) ++q;
        else if (nd.nodetype == qn::NodeType::Sink) ++snk;
    }
    return src == 1 && q == 1 && snk == 1 && L.nclosedjobs() == 0.0;
}

/**
 * The method names `solver_mva_qsys_analyzer` has an arm for.
 *
 * ONE PREDICATE FOR THE INTERCEPTION AND THE RUN. The Source-Queue-Sink shape is
 * claimed by that analyzer, which answers this FIXED list of closed forms and
 * refuses every other name -- so each general network method the report offers
 * on an open model was advertised on the one open shape it could not run on and
 * threw "not available for a model with one station and one class" the moment it
 * was asked for: `mva`, `amva`, `sum`, `esum`, `lin`, `gflin`, `egflin`, `qli`,
 * `fli`, `qd` and `qdlin`, eleven of them. They are NETWORK methods, and the
 * general branch solves a one-queue network exactly as it solves a larger one,
 * so the interception stands aside for them rather than claiming a model it
 * cannot answer. `qna` was the first name found this way and used to be excluded
 * by hand at the call site; it needs no special case now, having no arm here
 * either.
 *
 * Judged on the ALIASED name, the same `amva_method_alias` the analyzer applies
 * on the way in. Mirrors `matlab/src/solvers/MVA/@SolverMVA/mvaDispatch.m` and
 * its JAR and native python twins.
 */
inline bool qsys_serves_method(const std::string& method) {
    static const char* const kServed[] = {
        "default",     "exact",         "erlanga",       "mgisrgi",     "gigk.diffusion",
        "mm1",         "mmk",           "mg1",           "mgi1",        "gigk",
        "gigk.kingman_approx",          "gigk.whitt",    "gig1",        "gig1.allen",
        "gig1.kingman", "gig1.heyman",  "gig1.kobayashi", "gig1.klb",   "gig1.marchal",
        "gig1.gelenbe", "gig1.kimura",  "gig1.extremal", "qed",         "rqna",
        "rqt",         "gm1",           "gim1"};
    for (std::size_t i = 0; i < sizeof(kServed) / sizeof(*kServed); ++i)
        if (method == kServed[i]) return true;
    return false;
}

/** The same, with a Cache in place of the Queue. */
template <class T>
bool is_open_scs(const qn::NetworkStruct<T>& L) {
    if (L.nof_nodes() != 3) return false;
    int src = 0, ca = 0, snk = 0;
    for (const qn::NodeDef& nd : L.nodes) {
        if (nd.nodetype == qn::NodeType::Source) ++src;
        else if (nd.nodetype == qn::NodeType::Cache) ++ca;
        else if (nd.nodetype == qn::NodeType::Sink) ++snk;
    }
    return src == 1 && ca == 1 && snk == 1 && L.nclosedjobs() == 0.0;
}

/** The size-based disciplines the reference routes to its own analyzer. */
inline bool is_size_based(qn::SchedStrategy s) {
    return s == qn::SchedStrategy::SRPT || s == qn::SchedStrategy::PSJF ||
           s == qn::SchedStrategy::FB || s == qn::SchedStrategy::LRPT ||
           s == qn::SchedStrategy::SETF;
}

/** 1-based station index of the single node of this type, 0 when absent. */
template <class T>
std::size_t station_of_type(const qn::NetworkStruct<T>& L, qn::NodeType ty) {
    for (std::size_t i = 0; i < L.nof_nodes(); ++i)
        if (L.nodes[i].nodetype == ty) return L.nodes[i].station;
    return 0;
}

/** True when every finite SCV of a station's row is 1, as the reference tests. */
template <class T>
bool row_is_exponential(const qn::NetworkStruct<T>& L, std::size_t ist) {
    for (std::size_t r = 0; r < L.nclasses; ++r) {
        if (L.disabled[ist - 1][r]) continue;
        const double v = num_traits<T>::to_double(L.scv(ist - 1, r));
        if (!std::isfinite(v)) continue;
        if (std::fabs(v - 1.0) >= 1e-6) return false;
    }
    return true;
}

}  // namespace detail

/**
 * Port of `solver_mva_qsys_analyzer.m`: the closed forms for a single-class
 * open Source-Queue-Sink model.
 *
 * `exact` resolves to M/M/1, M/M/k, M/G/1 or G/M/1 and REFUSES anything else,
 * since no closed form covers it; `default` additionally falls back to the
 * G/G/k approximation and to the KLB G/G/1 approximation.
 */
template <class T>
DispatchResult<T> solver_mva_qsys_analyzer(const qn::NetworkStruct<T>& L, const MvaOptions& opt) {
    // The open-queue closed forms evaluate transcendentals (square roots, LSTs,
    // Brent roots), so under exact/Rational arithmetic the whole body is
    // discarded and the model is refused by name -- the field-arithmetic ladder
    // branches stay exact. Guarding the body (not just a static_assert) is what
    // lets mva_dispatch<Rational> compile at all.
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_mva_qsys_analyzer: the open queueing-system closed forms need transcendental "
            "arithmetic; rerun this model with --arith double or --arith real");
    } else {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations;
    const std::size_t src = detail::station_of_type(L, qn::NodeType::Source);
    const std::size_t q = detail::station_of_type(L, qn::NodeType::Queue);
    if (src == 0 || q == 0) throw InputError("solver_mva_qsys_analyzer: no Source-Queue pair");

    // The visit ratio of the queue, which for a feedback model exceeds one and
    // separates the per-visit quantities from the per-job ones.
    const std::size_t qstateful = L.stateful_of_station(q);
    const T Vq = L.visits[0](qstateful - 1, 0);
    const T srcRate = L.rates(src - 1, 0);
    const T lambda = T(srcRate * Vq);
    const T mu = L.rates(q - 1, 0);
    const double k = L.stations[q - 1].nservers;
    const unsigned ku = std::isfinite(k) ? static_cast<unsigned>(std::llround(k)) : 1u;
    const T ca = num_traits<T>::from_double(std::sqrt(num_traits<T>::to_double(L.scv(src - 1, 0))));
    const T cs = num_traits<T>::from_double(std::sqrt(num_traits<T>::to_double(L.scv(q - 1, 0))));
    const bool ca1 = std::fabs(num_traits<T>::to_double(ca) - 1.0) < 1e-12;
    const bool cs1 = std::fabs(num_traits<T>::to_double(cs) - 1.0) < 1e-12;

    // Finite-capacity loss branch (M/M/1/K with tail drop), the ONE finite
    // buffer this solver honours: the moment-based MacGregor Smith loss
    // probability, exact only at scv=1, with the queue length taken from the
    // truncated M/M/1/K distribution. Being an approximation in general it is
    // not offered under method='exact'; the capacity gate in the runner exempts
    // exactly this shape for every other method, so the two must agree.
    if (sn_is_mm1k_loss(L)) {
        if (opt.method == "exact")
            throw UnsupportedError(
                "solver_mva_qsys_analyzer: M/M/1/K tail-drop is solved by the approximate "
                "'mg1k.mgs' method (MacGregor Smith); it is not available under method='exact'. "
                "Use the default method, or SolverCTMC/SolverNC for an exact result");
        const double Kcap = L.cap[q - 1];
        const T rho = T(lambda / mu);
        const T Ploss =
            qsys::qsys_mg1k_loss_mgs(lambda, mu, T(cs * cs), static_cast<unsigned>(std::llround(Kcap)))
                .lossProbability;
        const T Tq = T(lambda * (one - Ploss));  // carried throughput
        const double rhod = num_traits<T>::to_double(rho);
        T Lsys = zero;
        if (std::fabs(rhod - 1.0) < 1e-10) {
            Lsys = num_traits<T>::from_double(Kcap / 2.0);  // L'Hopital limit at rho = 1
        } else {
            const T Kp1 = num_traits<T>::from_double(Kcap + 1.0);
            const T rKp1 = qsys::detail::num_pow(rho, Kp1);
            Lsys = T(rho / (one - rho) - Kp1 * rKp1 / (one - rKp1));
        }
        DispatchResult<T> lout;
        MvaSolution<T>& ls = lout.sol;
        ls.Q = Matrix<T>(M, 1, zero);
        ls.U = Matrix<T>(M, 1, zero);
        ls.R = Matrix<T>(M, 1, zero);
        ls.Tp = Matrix<T>(M, 1, zero);
        ls.X.assign(1, zero);
        ls.C.assign(1, zero);
        ls.method = opt.method;
        ls.iter = 1;
        ls.R(q - 1, 0) = T(Lsys / Tq);  // per-visit response time, by Little
        ls.Q(q - 1, 0) = Lsys;
        ls.U(q - 1, 0) = T(Tq / mu);  // single-server utilization
        ls.Tp(q - 1, 0) = Tq;         // carried (effective) rate
        ls.Tp(src - 1, 0) = lambda;   // offered arrival rate
        ls.X[0] = Tq;                 // system throughput = carried rate
        ls.C[0] = T(ls.R(q - 1, 0) * Vq);
        lout.actualmethod = "mg1k.mgs";
        return lout;
    }

    // Empty unless the queue reneges.
    const api::PatienceHandles<T> hpat = api::sn_patience_handles(L, q - 1, 0);

    std::string method = amva_method_alias(opt.method);
    if (method == "exact") {
        if (ca1 && cs1 && ku == 1) method = "mm1";
        else if (ca1 && cs1 && ku > 1) method = "mmk";
        else if (ca1 && ku == 1) method = "mg1";
        else if (cs1 && ku == 1) method = "gm1";
        else
            throw UnsupportedError(
                "solver_mva_qsys_analyzer: no exact closed form for this queueing system "
                "(neither arrivals nor service are exponential, or it has several servers)");
    } else if (method == "default") {
        // A station customers walk away from is a different model, not a
        // correction to one: nothing in the G/G/k family below carries an
        // abandonment rate, so the choice is made here and not by ca/cs.
        if (hpat.present) method = hpat.isExponential ? "erlanga" : "mgisrgi";
        else if (ca1 && cs1 && ku == 1) method = "mm1";
        else if (ca1 && cs1 && ku > 1) method = "mmk";
        else if (ca1 && ku == 1) method = "mg1";
        else if (cs1 && ku == 1) method = "gm1";
        else if (ku > 1) method = "gigk";
        else method = "gig1.klb";
    }

    // Whitt family, full metric set. These methods answer a station whose
    // CARRIED throughput is below the offered rate -- customers abandon, or are
    // blocked -- so Little's law on lambda would silently overstate the queue
    // and the common tail below cannot be used.
    if (method == "erlanga" || method == "mgisrgi" || method == "gigk.diffusion") {
        const double cap = L.cap[q - 1];
        // waiting spaces, servers excluded
        const double room = std::isfinite(cap)
                                ? std::max(0.0, cap - static_cast<double>(ku))
                                : std::numeric_limits<double>::infinity();
        T Lsys = zero, Tq = zero, Uq = zero;
        if (method == "gigk.diffusion") {
            const qsys::QsysGgnmResult<T> d = qsys::qsys_ggnm_diffusion(lambda, mu, ku, room, ca, cs);
            Lsys = d.meanNumber;
            Tq = d.throughput;
            Uq = d.utilization;
        } else {
            if (!hpat.present)
                throw UnsupportedError("solver_mva_qsys_analyzer: method '" + method +
                                       "' needs a reneging patience law on the queue");
            const qsys::QsysAbandonResult<T> ab =
                (method == "erlanga" || hpat.isExponential)
                    // Exponential patience makes the state-dependent
                    // approximation exact, so take the exact chain either way.
                    ? qsys::qsys_erlanga(lambda, mu, hpat.rate, ku, room)
                    : qsys::qsys_mgisrgi_whitt(lambda, mu, ku, room, hpat.as_patience());
            Lsys = ab.meanNumber;
            Tq = ab.throughput;
            Uq = ab.utilization;
        }
        DispatchResult<T> aout;
        MvaSolution<T>& as = aout.sol;
        as.Q = Matrix<T>(M, 1, zero);
        as.U = Matrix<T>(M, 1, zero);
        as.R = Matrix<T>(M, 1, zero);
        as.Tp = Matrix<T>(M, 1, zero);
        as.X.assign(1, zero);
        as.C.assign(1, zero);
        as.method = opt.method;
        as.iter = 1;
        // Little's law on the CARRIED rate, as the loss branch above and
        // SolverCTMC report it.
        as.R(q - 1, 0) = Tq > zero ? T(Lsys / Tq) : zero;
        as.Q(q - 1, 0) = Lsys;
        as.U(q - 1, 0) = Uq;
        as.Tp(q - 1, 0) = Tq;         // carried rate
        as.Tp(src - 1, 0) = srcRate;  // offered arrival rate
        as.X[0] = Tq;
        as.C[0] = T(as.R(q - 1, 0) * Vq);
        aout.actualmethod = method;
        return aout;
    }

    T R = zero;
    if (method == "mm1") R = qsys::qsys_mm1(lambda, mu).W;
    else if (method == "mmk") R = qsys::qsys_mmk(lambda, mu, ku).W;
    else if (method == "mg1" || method == "mgi1") R = qsys::qsys_mg1(lambda, mu, cs).W;
    else if (method == "gigk") R = qsys::qsys_gigk_approx(lambda, mu, ca, cs, ku).W;
    else if (method == "gigk.kingman_approx")
        R = qsys::qsys_gigk_approx_kingman(lambda, mu, ca, cs, ku).W;
    else if (method == "gig1.kingman") R = qsys::qsys_gig1_ubnd_kingman(lambda, mu, ca, cs).W;
    else if (method == "gig1.heyman") R = qsys::qsys_gig1_approx_heyman(lambda, mu, ca, cs).W;
    else if (method == "gig1" || method == "gig1.allen")
        R = qsys::qsys_gig1_approx_allencunneen(lambda, mu, ca, cs).W;
    else if (method == "gig1.kobayashi") R = qsys::qsys_gig1_approx_kobayashi(lambda, mu, ca, cs).W;
    else if (method == "gig1.klb") R = qsys::qsys_gig1_approx_klb(lambda, mu, ca, cs).W;
    else if (method == "gig1.marchal") R = qsys::qsys_gig1_approx_marchal(lambda, mu, ca, cs).W;
    else if (method == "gig1.gelenbe") R = qsys::qsys_gig1_approx_gelenbe(lambda, mu, ca, cs).W;
    else if (method == "gig1.kimura") R = qsys::qsys_gig1_approx_kimura(lambda, mu, ca, cs).W;
    else if (method == "gigk.whitt") R = qsys::qsys_gigk_approx_whitt(lambda, mu, ca, cs, ku).W;
    // RQNA and RQT reach the SINGLE-QUEUE robust formulas, not the network
    // analyzers of solver_rqna.h / solver_rqt.h: this shape is one node, and
    // the reference answers it here (solver_mva_qsys_analyzer.m cases rqna
    // and rqt). Without these two arms the qsys interception below claimed the
    // model and then refused the method, so `rqna` and `rqt` were advertised
    // on every Source-Queue-Sink model and ran on none of them.
    else if (method == "rqna") {
        // The arrival flow enters through its index of dispersion for counts.
        const mam::Map<T> arvMAP = lang::dist_to_map(L.service[src - 1][0]);
        const T rho1 = lambda / mu;
        auto IaFun = [&](const T& x) { return mam::map_count_idc(arvMAP, x); };
        R = T(qsys::qsys_gig1_rq(rho1, mu, T(cs * cs), IaFun).W + one / mu);
    } else if (method == "rqt") {
        // Polyhedral uncertainty sets for the arrival and service flows.
        const T rho1 = lambda / (num_traits<T>::from_int(static_cast<int>(ku)) * mu);
        const T Gamma_a = ca / lambda;
        const T two = num_traits<T>::from_int(2);
        const T Gamma_s =
            qsys::qsys_gigk_rqt_gamma(rho1, mu, Gamma_a, T(cs / mu), ku, two);
        R = qsys::qsys_gigk_rqt(lambda, mu, Gamma_a, Gamma_s, ku, two, two).W;
    }
    else if (method == "qed") R = T(qsys::qsys_mmk_qed(lambda, mu, ku).meanWait + one / mu);
    // The upper end, gig1.kingman already reporting a bound.
    else if (method == "gig1.extremal")
        R = T(qsys::qsys_gig1_bnds_extremal(lambda, mu, ca, cs).upperBound + one / mu);
    else if (method == "gm1" || method == "gim1") {
        // The PH path is exact only when sn.proc holds the arrival law ITSELF.
        // For a non-Markovian arrival the refresh substitutes an Erlang-n SCV
        // fit, so the PH sigma-root would answer a different model; those cases
        // take the transform's sigma-root instead, as the reference does.
        const lang::Distrib<T>& arv = L.service[src - 1][0];
        const bool markovian = arv.has_map();
        bool done = false;
        if (markovian) {
            const mam::Map<T> m = lang::dist_to_map(arv);
            const std::vector<T> pie = mam::map_pie(m);
            R = qsys::qsys_phm1(pie, m.D0, mu, num_traits<T>::from_double(1e-16)).meanSojournTime;
            done = true;
        }
        if (!done) {
            // sigma solves A*(mu - mu sigma) = sigma on (0,1)
            //
            // THE UPPER BRACKET CANNOT SIT AT 1. sigma = 1 is ALWAYS a root of
            // this equation, and near it the transform is evaluated at s = mu(1-x)
            // -> 0, where A*(s) is a difference of two exponentials that both
            // tend to 1: at x = 1-1e-12 the cancellation leaves A* correct only to
            // about 1e-4, so f came out POSITIVE at both ends and Brent refused
            // the bracket outright ("the bracket endpoints do not straddle a
            // root"), which killed SolverMVA on every G/M/1 with a non-Markovian
            // arrival -- gallery_um1 (Uniform(1,2) -> Exp(2)) among them. f(0) is
            // A*(mu) > 0 always, so only the upper end has to be found: walk it
            // toward 1 until the sign turns, exactly as far as the sought root
            // requires and no further. The reference sidesteps the same trap by
            // seeding an UNBRACKETED fzero at 0.5 (solver_mva_qsys_analyzer.m).
            auto f = [&](const T& x) { return T(lang::dist_lst(arv, T(mu - mu * x)) - x); };
            double hi = 0.5;
            bool straddles = false;
            for (int i = 0; i < 8; ++i) {
                if (num_traits<T>::to_double(f(num_traits<T>::from_double(hi))) < 0) {
                    straddles = true;
                    break;
                }
                hi = 1.0 - (1.0 - hi) / 10.0;
            }
            const RootResult<T> rr =
                straddles ? root_brent(f, num_traits<T>::from_double(1e-12),
                                       num_traits<T>::from_double(hi),
                                       num_traits<T>::from_double(1e-14))
                          : RootResult<T>();
            if (straddles && rr.converged) {
                R = qsys::qsys_gm1(rr.root, mu);
            } else {
                R = qsys::qsys_gg1(lambda, mu, T(ca * ca), one).W;
            }
        }
    } else {
        throw UnsupportedError("solver_mva_qsys_analyzer: method '" + method +
                               "' is not available for a model with one station and one class");
    }

    DispatchResult<T> out;
    MvaSolution<T>& s = out.sol;
    s.Q = Matrix<T>(M, 1, zero);
    s.U = Matrix<T>(M, 1, zero);
    s.R = Matrix<T>(M, 1, zero);
    s.Tp = Matrix<T>(M, 1, zero);
    s.X.assign(1, zero);
    s.C.assign(1, zero);
    s.method = opt.method;
    s.iter = 1;
    // Per-visit response time and queue length at the queue; per-job cycle time
    // and system throughput. For a feedback model (Vq > 1) the two differ.
    s.R(q - 1, 0) = R;
    s.C[0] = T(R * Vq);
    s.X[0] = srcRate;
    s.U(q - 1, 0) = T(lambda / mu / num_traits<T>::from_double(std::isfinite(k) ? k : 1.0));
    s.Tp(src - 1, 0) = srcRate;
    s.Tp(q - 1, 0) = lambda;
    s.Q(q - 1, 0) = T(lambda * R);
    out.actualmethod = method;
    return out;
    }  // if constexpr has_transcendental
}

/**
 * Port of `solver_mva_qsys_prio_analyzer.m`: the exact Cobham formula for a
 * single open M/G/1 queue with non-preemptive priorities.
 *
 * The classes are handed to `qsys_mg1_prio` in PRIORITY order (lowest value
 * first), which is the order the formula's nested sums assume.
 */
template <class T>
DispatchResult<T> solver_mva_qsys_prio_analyzer(const qn::NetworkStruct<T>& L,
                                                const MvaOptions& opt) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses;
    const std::size_t src = detail::station_of_type(L, qn::NodeType::Source);
    const std::size_t q = detail::station_of_type(L, qn::NodeType::Queue);

    std::vector<std::size_t> order(K);
    for (std::size_t r = 0; r < K; ++r) order[r] = r;
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return L.classes[a].prio < L.classes[b].prio;
    });

    std::vector<T> lam(K, zero), mu(K, zero), cs(K, one);
    for (std::size_t j = 0; j < K; ++j) {
        const std::size_t r = order[j];
        lam[j] = L.rates(src - 1, r);
        mu[j] = L.rates(q - 1, r);
        const double v = num_traits<T>::to_double(L.scv(q - 1, r));
        cs[j] = (std::isfinite(v) && v > 0.0) ? num_traits<T>::from_double(std::sqrt(v)) : one;
    }
    const std::vector<T> W = qsys::qsys_mg1_prio(lam, mu, cs).W;

    DispatchResult<T> out;
    MvaSolution<T>& s = out.sol;
    s.Q = Matrix<T>(M, K, zero);
    s.U = Matrix<T>(M, K, zero);
    s.R = Matrix<T>(M, K, zero);
    s.Tp = Matrix<T>(M, K, zero);
    s.C.assign(K, zero);
    s.X.assign(K, zero);
    s.method = opt.method;
    s.iter = 0;
    for (std::size_t j = 0; j < K; ++j) {
        const std::size_t r = order[j];
        if (!(lam[j] > zero) || !std::isfinite(num_traits<T>::to_double(mu[j]))) continue;
        s.R(q - 1, r) = W[j];
        s.C[r] = W[j];
        s.X[r] = lam[j];
        s.U(q - 1, r) = T(lam[j] / mu[j]);
        s.Tp(q - 1, r) = lam[j];
        s.Tp(src - 1, r) = lam[j];
        s.Q(q - 1, r) = T(lam[j] * W[j]);
    }
    out.actualmethod = "mg1.prio";
    return out;
}

/**
 * The exact multiclass M/M/1-DPS solve the reference inlines in mvaDispatch:
 * the truncated multiclass chain of `qsys_mm1_dps`, which the AMVA cross-term
 * correction cannot reproduce (it violates the equal-rate conservation law).
 */
template <class T>
DispatchResult<T> solver_mva_dps_exact(const qn::NetworkStruct<T>& L, const MvaOptions& opt) {
    // The truncated-chain M/M/1-DPS closed form (qsys_mm1_dps) needs
    // transcendental arithmetic; refuse by name under exact/Rational.
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_mva_dps_exact: the exact DPS closed form needs transcendental arithmetic; "
            "rerun this model with --arith double or --arith real");
    } else {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses;
    const std::size_t src = detail::station_of_type(L, qn::NodeType::Source);
    const std::size_t q = detail::station_of_type(L, qn::NodeType::Queue);

    std::vector<T> lam(K, zero), mu(K, zero), w(K, one);
    for (std::size_t r = 0; r < K; ++r) {
        lam[r] = L.rates(src - 1, r);
        mu[r] = L.rates(q - 1, r);
        const std::vector<T>& sp = L.stations[q - 1].schedparam;
        if (sp.size() > r && sp[r] > zero) w[r] = sp[r];
    }
    const std::vector<T> Wd = qsys::qsys_mm1_dps(lam, mu, w).T_;

    DispatchResult<T> out;
    MvaSolution<T>& s = out.sol;
    s.Q = Matrix<T>(M, K, zero);
    s.U = Matrix<T>(M, K, zero);
    s.R = Matrix<T>(M, K, zero);
    s.Tp = Matrix<T>(M, K, zero);
    s.C.assign(K, zero);
    s.X.assign(K, zero);
    s.method = opt.method;
    s.iter = 0;
    for (std::size_t r = 0; r < K; ++r) {
        s.R(q - 1, r) = Wd[r];
        s.C[r] = Wd[r];
        s.X[r] = lam[r];
        if (mu[r] > zero) s.U(q - 1, r) = T(lam[r] / mu[r]);
        s.Tp(q - 1, r) = lam[r];
        s.Tp(src - 1, r) = lam[r];
        s.Q(q - 1, r) = T(lam[r] * Wd[r]);
    }
    out.actualmethod = "mm1.dps";
    return out;
    }  // if constexpr has_transcendental
}

/**
 * Port of `solver_mvald_analyzer.m`: the load-dependent branch.
 *
 * `exact` and `mva` take the exact load-dependent recursion; `default` takes it
 * too on a small closed product-form model, mirroring the default-to-exact
 * upgrade the non-LD analyzer applies; everything else goes to the AMVA, whose
 * lld and cd terms are now in place.
 */
template <class T>
DispatchResult<T> solver_mvald_analyzer(const qn::NetworkStruct<T>& L, const MvaOptions& opt,
                                        const Matrix<T>& init_sol) {
    const std::string method = amva_method_alias(opt.method);
    const bool has_cd = [&] {
        for (const auto& st : L.stations)
            if (st.cdscaling) return true;
        return false;
    }();

    if (method == "exact" || method == "mva") {
        if (has_cd)
            throw UnsupportedError(
                "solver_mvald_analyzer: there is no exact class-dependent solver in MVA");
        DispatchResult<T> out;
        out.sol = solver_mvald(L, opt);
        out.actualmethod = "exact";
        return out;
    }
    if (!(method == "default" || method == "amva" || method == "qd" || method == "lin" ||
          method == "qdlin"))
        throw UnsupportedError("solver_mvald_analyzer: the '" + method +
                               "' method is not supported by the load-dependent MVA solver");

    if (method == "default" && !has_cd) {
        double Nsum = 0.0;
        bool integral = true, finite = true;
        for (const qn::JobClass& c : L.classes) {
            if (!std::isfinite(c.population)) finite = false;
            else {
                Nsum += c.population;
                if (c.population != std::floor(c.population)) integral = false;
            }
        }
        if (finite && integral && L.nchains <= 4 && Nsum <= 20.0 && L.has_product_form()) {
            DispatchResult<T> out;
            out.sol = solver_mvald(L, opt);
            out.actualmethod = "exact";
            return out;
        }
    }
    bool converged = true;
    const ChainDemands<T> d = sn_get_demands_chain(L);
    MvaOptions o = opt;
    DispatchResult<T> out;
    out.sol = solver_amva(L, d, o, init_sol, converged);
    // solver_amva computes the outer residual and reports it here; dropping the
    // flag on the floor left every caller with only `iter`, which on this route
    // aggregates the nested sweeps and cannot decide convergence.
    out.sol.converged = converged;
    out.actualmethod = out.sol.method;
    return out;
}

/**
 * Port of `solver_mva_qsys_sizebased_analyzer.m`: the M/G/1 formulas for the
 * size-based disciplines (Wierman and Harchol-Balter, SIGMETRICS 2003).
 *
 * The arrival rate of each class is the source rate times the queue's VISIT
 * ratio, and the response time comes back per visit, so both are multiplied
 * out the same way the reference does.
 *
 * PER-VISIT VERSUS PER-JOB. `W` comes back per visit, so the response time is
 * `W * visits` while the queue length stays `lambda * W` with lambda already
 * carrying the visit ratio. That makes `Q != T * R` at the station whenever the
 * visit ratio differs from one: the reference inflates R by the visits but not
 * Q, so Little's law fails on a feedback routing. It is a defect of the
 * reference, reproduced here rather than silently corrected -- this branch is
 * only ever dispatched on a plain Source-Queue-Sink, where the visit ratio is
 * one and the discrepancy cannot arise.
 */
template <class T>
DispatchResult<T> solver_mva_qsys_sizebased_analyzer(const qn::NetworkStruct<T>& L,
                                                     const MvaOptions& opt) {
    // The M/G/1 size-based closed forms (SRPT/PSJF/FB/LRPT/SETF) evaluate
    // transcendentals; refuse by name under exact/Rational.
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_mva_qsys_sizebased_analyzer: the M/G/1 size-based closed forms need "
            "transcendental arithmetic; rerun this model with --arith double or --arith real");
    } else {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses;
    const std::size_t src = detail::station_of_type(L, qn::NodeType::Source);
    const std::size_t q = detail::station_of_type(L, qn::NodeType::Queue);
    const std::size_t qstateful = L.stateful_of_station(q);

    // sn.visits is indexed by CHAIN. The reference wrote `sn.visits{source_ist}`
    // -- the chain whose number happens to equal the Source's STATION index --
    // so every class outside chain 1 read a visit of 0 and the analyzer rejected
    // its own multiclass model. Fixed in MATLAB in the same change; take the
    // chain each class belongs to.
    std::vector<std::size_t> chain_of(K, 0);
    for (std::size_t c = 0; c < L.nchains; ++c)
        for (std::size_t r : L.inchain[c]) chain_of[r - 1] = c;

    std::vector<T> lambda(K, zero), mu(K, zero), cs(K, one), vis(K, one);
    for (std::size_t r = 0; r < K; ++r) {
        vis[r] = L.visits[chain_of[r]](qstateful - 1, r);
        lambda[r] = T(L.rates(src - 1, r) * vis[r]);
        mu[r] = L.rates(q - 1, r);
        const double v = num_traits<T>::to_double(L.scv(q - 1, r));
        cs[r] = (std::isfinite(v) && v > 0.0) ? num_traits<T>::from_double(std::sqrt(v)) : one;
        if (!(lambda[r] > zero) || !(mu[r] > zero))
            throw InputError(
                "solver_mva_qsys_sizebased_analyzer: the arrival and service rates must be "
                "positive");
    }

    const qn::SchedStrategy sched = L.stations[q - 1].sched;
    std::vector<T> W;
    std::string actual;
    switch (sched) {
        case qn::SchedStrategy::SRPT:
            W = qsys::qsys_mg1_srpt(lambda, mu, cs).W;
            actual = "mg1.srpt";
            break;
        case qn::SchedStrategy::PSJF:
            W = qsys::qsys_mg1_psjf(lambda, mu, cs).W;
            actual = "mg1.psjf";
            break;
        case qn::SchedStrategy::FB:
            W = qsys::qsys_mg1_fb(lambda, mu, cs).W;
            actual = "mg1.fb";
            break;
        case qn::SchedStrategy::LRPT:
            W = qsys::qsys_mg1_lrpt(lambda, mu, cs).W;
            actual = "mg1.lrpt";
            break;
        case qn::SchedStrategy::SETF:
            W = qsys::qsys_mg1_setf(lambda, mu, cs).W;
            actual = "mg1.setf";
            break;
        default:
            throw UnsupportedError(std::string("solver_mva_qsys_sizebased_analyzer: ") +
                                   lang::sched_to_text(sched) + " is not a size-based discipline");
    }

    DispatchResult<T> out;
    MvaSolution<T>& s = out.sol;
    s.Q = Matrix<T>(M, K, zero);
    s.U = Matrix<T>(M, K, zero);
    s.R = Matrix<T>(M, K, zero);
    s.Tp = Matrix<T>(M, K, zero);
    s.C.assign(K, zero);
    s.X.assign(K, zero);
    s.method = opt.method;
    s.iter = 1;
    for (std::size_t r = 0; r < K; ++r) {
        s.R(q - 1, r) = T(W[r] * vis[r]);
        s.C[r] = s.R(q - 1, r);
        s.X[r] = lambda[r];
        s.U(q - 1, r) = T(lambda[r] / mu[r]);
        s.Tp(q - 1, r) = lambda[r];
        s.Tp(src - 1, r) = lambda[r];
        s.Q(q - 1, r) = T(lambda[r] * W[r]);
    }
    out.actualmethod = actual;
    return out;
    }  // if constexpr has_transcendental
}

/**
 * Port of `solver_mva_marie_analyzer.m`: Marie's iterative
 * aggregation-decomposition for a CLOSED network with non-exponential FCFS
 * service.
 *
 * Only FCFS is service-time sensitive; PS and LCFSPR are insensitive
 * (product form) and are handed an SCV of one, which is what makes the
 * decomposition exact for them. Delay stations fold into the per-chain think
 * time rather than entering the isolation, and multiserver isolation is
 * single-chain only -- all three restrictions are the reference's, and each is
 * refused by name rather than approximated.
 */
template <class T>
DispatchResult<T> solver_mva_marie_analyzer(const qn::NetworkStruct<T>& L, const MvaOptions& opt) {
    // Marie's decomposition fits a Coxian per isolated station (square roots) and
    // stops on a tolerance; refuse by name under exact/Rational.
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_mva_marie_analyzer: Marie's aggregation-decomposition needs transcendental "
            "arithmetic; rerun this model with --arith double or --arith real");
    } else {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const ChainDemands<T> d = sn_get_demands_chain(L);
    const std::size_t M = L.nstations, C = L.nchains;

    for (std::size_t c = 0; c < C; ++c)
        if (!std::isfinite(d.Nchain[c]))
            throw UnsupportedError(
                "solver_mva_marie_analyzer: the 'marie' method supports closed models only; this "
                "model has open classes");
    for (const auto& nd : L.nodes)
        if (nd.nodetype == qn::NodeType::Source)
            throw UnsupportedError(
                "solver_mva_marie_analyzer: the 'marie' method supports closed models only; this "
                "model has a Source");

    std::vector<std::size_t> qrows, drows;
    for (std::size_t i = 0; i < M; ++i) {
        const qn::SchedStrategy sc = L.stations[i].sched;
        const bool delay = std::isinf(L.stations[i].nservers) || sc == qn::SchedStrategy::INF;
        if (delay) {
            drows.push_back(i);
            continue;
        }
        if (!(sc == qn::SchedStrategy::FCFS || sc == qn::SchedStrategy::PS ||
              sc == qn::SchedStrategy::LCFSPR))
            throw UnsupportedError(std::string("solver_mva_marie_analyzer: the 'marie' method "
                                               "supports FCFS, PS, LCFSPR and Delay stations "
                                               "only; station '") +
                                   L.stations[i].name + "' is " + lang::sched_to_text(sc));
        qrows.push_back(i);
    }

    std::vector<T> Z(C, zero);
    for (std::size_t c = 0; c < C; ++c)
        for (std::size_t i : drows) Z[c] += d.Lchain(i, c);

    const std::size_t Mq = qrows.size();
    Matrix<T> Lq(Mq, C, zero), SCV(Mq, C, one);
    std::vector<int> ns(Mq, 1);
    bool multiserver = false;
    for (std::size_t a = 0; a < Mq; ++a) {
        const std::size_t i = qrows[a];
        const double srv = L.stations[i].nservers;
        ns[a] = std::isfinite(srv) ? static_cast<int>(std::llround(srv)) : 1;
        if (ns[a] > 1) multiserver = true;
        for (std::size_t c = 0; c < C; ++c) {
            Lq(a, c) = d.Lchain(i, c);
            if (L.stations[i].sched != qn::SchedStrategy::FCFS) continue;
            const double v = num_traits<T>::to_double(d.SCVchain(i, c));
            if (std::isfinite(v) && v > 0.0) SCV(a, c) = num_traits<T>::from_double(v);
        }
    }
    if (C > 1 && multiserver)
        throw UnsupportedError(
            "solver_mva_marie_analyzer: the 'marie' method supports multiserver stations for "
            "single-chain models only");

    std::vector<int> N(C, 0);
    for (std::size_t c = 0; c < C; ++c) N[c] = static_cast<int>(std::llround(d.Nchain[c]));
    pfqn::MarieResult<T> mr;
    if (Mq == 0) {
        // Nothing to isolate: with every station an infinite server the
        // aggregation-decomposition degenerates to the exact delay solution
        // X_c = N_c / Z_c, and pfqn_marie would be handed a zero-row demand matrix.
        mr.X.assign(C, zero);
        mr.Q = Matrix<T>(0, C, zero);
        mr.U = Matrix<T>(0, C, zero);
        mr.it = 1;
        for (std::size_t c = 0; c < C; ++c)
            if (Z[c] > zero) mr.X[c] = T(num_traits<T>::from_double(d.Nchain[c]) / Z[c]);
    } else {
        mr = C == 1 ? pfqn::pfqn_marie(Lq, N, Z, SCV, 1e-8, 1000, ns)
                    : pfqn::pfqn_marie(Lq, N, Z, SCV, 1e-8, 1000, std::vector<int>());
    }

    Matrix<T> Qchain(M, C, zero), Uchain(M, C, zero), Rchain(M, C, zero), Tchain(M, C, zero);
    for (std::size_t a = 0; a < Mq; ++a)
        for (std::size_t c = 0; c < C; ++c) {
            Qchain(qrows[a], c) = mr.Q(a, c);
            Uchain(qrows[a], c) = mr.U(a, c);
        }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < C; ++c) Tchain(i, c) = T(mr.X[c] * d.Vchain(i, c));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < C; ++c)
            if (Tchain(i, c) > zero) Rchain(i, c) = T(Qchain(i, c) / Tchain(i, c));
    // A delay station holds T S jobs, all of them in service.
    for (std::size_t i : drows)
        for (std::size_t c = 0; c < C; ++c) {
            Qchain(i, c) = T(Tchain(i, c) * d.STchain(i, c));
            Uchain(i, c) = Qchain(i, c);
            Rchain(i, c) = d.STchain(i, c);
        }

    const ClassResults<T> cr =
        sn_deaggregate_chain_results(L, d, Qchain, Uchain, Rchain, Tchain, mr.X);
    DispatchResult<T> out;
    out.sol.Q = cr.Q;
    out.sol.U = cr.U;
    out.sol.R = cr.R;
    out.sol.Tp = cr.Tp;
    out.sol.C = cr.C;
    out.sol.X = cr.X;
    out.sol.method = opt.method;
    out.sol.iter = mr.it;
    out.actualmethod = "marie";
    return out;
    }  // if constexpr has_transcendental
}

/**
 * The ladder itself. `init_sol` is the warm start the outer solver carries
 * across iterations; it reaches the AMVA paths only, as it does in the
 * reference.
 */
template <class T>
DispatchResult<T> mva_dispatch(const qn::NetworkStruct<T>& L, const MvaOptions& opt,
                               const Matrix<T>& init_sol) {
    using qn::NodeType;
    using qn::SchedStrategy;
    const std::string method = amva_method_alias(opt.method);

    // 0. closed models with a shortest-job-next station. This sits ABOVE every
    // other branch because the reference tests hasSJN first and the branches
    // are not disjoint: an SJF station in a closed model with a delay would
    // otherwise reach the generic AMVA path, which reads only the mean service
    // time and would report the station as if it scheduled size-blind. The
    // OPEN case has no population to recur over and is refused by name rather
    // than solved as though the discipline did nothing.
    if (sn_has_sjn(L)) {
        bool open = false;
        for (const auto& c : L.classes)
            if (std::isinf(c.population)) open = true;
        if (open)
            throw UnsupportedError(
                "SolverMVA supports shortest-job-next (SJF) scheduling only in closed models, the "
                "conditional waiting time equation being a population recursion. Use SolverLDES, "
                "or SolverMVA with SRPT or PSJF for the preemptive size-based open queue");
        const SjnAnalyzerResult<T> sr = solver_mva_sjn_analyzer(L, opt);
        DispatchResult<T> sout;
        sout.sol = sr.sol;
        sout.actualmethod = sr.actualmethod;
        sout.warning = sr.warning;
        return sout;
    }

    // 1. order-independent / pass-and-swap stations. Detection is
    // `nc_is_oi_model`'s: PAS or OI scheduling, a rate function, and an
    // ALL-ZERO swap graph. A station with those disciplines that does not
    // qualify is a genuine pass-and-swap station, which is not product-form
    // and has no exact mean-value analyzer, so it is refused by name rather
    // than handed to an AMVA that cannot represent its rate function at all.
    // The gate is the reference's, and all three conjuncts matter: an OI
    // station must be present, the WHOLE model must be order-independent
    // (nc_is_oi_model: closed, and every other station product-form), and the
    // method must be `default` or `exact`. MVA reaches OI/PAS ONLY through the
    // exact analyzer; any other method would fall through to an AMVA that only
    // ever sees the single-job rates and would silently report a zero queue
    // length at the OI station.
    {
        const bool has_oi_station = std::any_of(
            L.stations.begin(), L.stations.end(), [](const qn::Station<T>& st) {
                return st.sched == SchedStrategy::OI || st.sched == SchedStrategy::PAS;
            });
        const bool exact_method = (opt.method == "default" || opt.method == "exact");
        if (!find_oi_stations(L).empty() && nc_is_oi_model(L) && exact_method) {
            DispatchResult<T> oout;
            oout.sol = solver_mva_oi_analyzer(L, opt);
            oout.actualmethod = "oi";
            return oout;
        }
        if (has_oi_station)
            throw UnsupportedError(
                "SolverMVA supports order-independent (OI) and pass-and-swap (PAS) stations only "
                "through its exact order-independent analyzer, which requires method 'default' or "
                "'exact' (got '" + opt.method +
                "'), an empty/zero swap graph at every OI/PAS station, a closed model, and every "
                "other station to be product-form (INF, PS, LCFS-PR, SIRO, or "
                "class-independent-rate FCFS). Use SolverCTMC or SolverLDES for this model");
    }

    // 2. delayed-hit retrieval caches (Cache.setRetrievalSystem). The two
    // variants are DIFFERENT METHODS, not one method on two topologies, which
    // is why the Source test decides between them rather than parameterizing
    // one call. The OPEN model has no population to recur over, so the hit /
    // miss / delayed-hit split comes from the open fixed point
    // (`retrieval_fpi`) and all three are reported. The CLOSED model has no
    // exogenous arrival rate to run that fixed point on, so the split emerges
    // from a decomposition-aggregation sweep instead, and the delayed-hit
    // fraction folds into the miss (see solver_mva_cacheqn_retrieval.h).
    for (const auto& kv : L.nodeparam) {
        if (kv.second.retrieval_capacity <= 0) continue;
        DispatchResult<T> rout;
        if (detail::station_of_type(L, NodeType::Source) == 0) {
            const MvaCacheqnRetrievalSolution<T> cr = solver_mva_cacheqn_retrieval_analyzer(L, opt);
            rout.sol = cr.sol;
            rout.cache = solvers::cache_metrics_of(L, cr.hitprob, cr.missprob, cr.delayedprob,
                                                   cr.latency, cr.hitproblist, Matrix<T>(),
                                                   std::vector<T>());
        } else {
            MvaRetrievalCacheOutputs<T> co;
            rout.sol = solver_mva_retrieval_analyzer(L, opt, &co);
            rout.cache = solvers::cache_metrics_of(L, co.hitprob, co.missprob, co.delayedprob,
                                                   co.latency, co.hitproblist, co.itemprob,
                                                   std::vector<T>());
        }
        rout.actualmethod = "fpi";
        return rout;
    }

    const bool sqs = detail::is_open_sqs(L);
    const std::size_t qst = detail::station_of_type(L, NodeType::Queue);
    const SchedStrategy qsched = qst ? L.stations[qst - 1].sched : SchedStrategy::NONE;
    const std::size_t srcst = detail::station_of_type(L, NodeType::Source);

    // 8. non-reentrant cache: a Source-Cache-Sink model
    if (detail::is_open_scs(L)) {
        DispatchResult<T> cout;
        const CacheResult<T> cr = solver_mva_cache_analyzer(L, opt);
        cout.sol = cr.sol;
        cout.cache = solvers::cache_metrics_of(L, cr.hitprob, cr.missprob, std::vector<T>(),
                                               std::vector<T>(), cr.hitproblist, cr.itemprob,
                                               std::vector<T>());
        cout.actualmethod = cr.actualmethod;
        return cout;
    }

    // 9. any OTHER cache -- one embedded in a queueing network -- is the
    // integrated caching-queueing analyzer (decomposition-aggregation).
    if (!L.nodeparam.empty()) {
        DispatchResult<T> qout;
        std::shared_ptr<qn::NetworkStruct<T>> refreshed(new qn::NetworkStruct<T>());
        MvaCacheqnCacheOutputs<T> co;
        qout.sol = solver_mva_cacheqn_analyzer(L, opt, refreshed.get(), &co);
        qout.refreshed_struct = refreshed;
        // The integrated branch reports the split per (cache, class) and no
        // delayed-hit fraction; the per-item law rides in separately because
        // only this branch has one per cache rather than one per model.
        qout.cache = solvers::cache_metrics_of_matrix(L, co.hitprob, co.missprob);
        for (std::size_t ci = 0; ci < qout.cache.caches.size() && ci < co.itemprob.size(); ++ci)
            qout.cache.caches[ci].itemprob = co.itemprob[ci];
        qout.actualmethod = (opt.method == "exact") ? "exact" : "default";
        return qout;
    }

    // 3. size-based scheduling
    if (sqs && detail::is_size_based(qsched)) return solver_mva_qsys_sizebased_analyzer(L, opt);

    // 4. single-class open Source-Queue-Sink.
    //
    // CLAIMED ONLY FOR THE NAMES THE ANALYZER ANSWERS, `detail::qsys_serves_method`.
    // Anything else -- 'qna', 'mva', 'amva', the linearizers, the qd family,
    // sum/esum -- falls through to the network branches below and is solved
    // there, as it is on every larger open network; this branch used to claim
    // the model for them and then refuse the method by name. 'rqna' and 'rqt'
    // are served here: the analyzer answers them with single-queue robust
    // formulas of their own.
    //
    // THE M/M/1/K LOSS SHAPE IS THE EXCEPTION and is claimed whatever the name.
    // Its branch is chosen by the SHAPE and not by the method, and the capacity
    // gate in the runner exempts exactly this shape for every method but
    // 'exact', so the two have to agree on which models come here.
    if (sqs && L.nclasses == 1 && (detail::qsys_serves_method(method) || sn_is_mm1k_loss(L)))
        return solver_mva_qsys_analyzer(L, opt);

    // 5. multiclass open polling
    if (sqs && L.nclasses > 1 && qsched == SchedStrategy::POLLING) {
        DispatchResult<T> pout;
        pout.sol = solver_mva_polling_analyzer(L, opt);
        pout.actualmethod = "stationtime";
        return pout;
    }

    // 6. multiclass open HOL priority, single server, Poisson arrivals
    if (sqs && L.nclasses > 1 && qsched == SchedStrategy::HOL &&
        L.stations[qst - 1].nservers == 1.0 && detail::row_is_exponential(L, srcst))
        return solver_mva_qsys_prio_analyzer(L, opt);

    // 7. multiclass open DPS, at most three classes, exponential everywhere
    if (sqs && L.nclasses > 1 && L.nclasses <= 3 && qsched == SchedStrategy::DPS &&
        L.stations[qst - 1].nservers == 1.0 && detail::row_is_exponential(L, srcst) &&
        detail::row_is_exponential(L, qst))
        return solver_mva_dps_exact(L, opt);

    // 10. the bound family, which the reference moved to SolverBA
    static const char* kBounds[] = {"aba.upper", "aba.lower", "bjb.upper",  "bjb.lower",
                                    "pb.upper",  "pb.lower",  "gb.upper",   "gb.lower",
                                    "sb.upper",  "sb.lower",  "mwba.upper", "mwba.lower"};
    for (const char* b : kBounds)
        if (method == b)
            throw UnsupportedError("SolverMVA: bound methods have moved to SolverBA; use "
                                   "SolverBA with method '" +
                                   method + "'");

    // 11. Marie's aggregation-decomposition
    if (method == "marie") return solver_mva_marie_analyzer(L, opt);

    // 12. load- or class-dependent scaling
    bool has_scaling = false;
    for (const auto& st : L.stations)
        if (!st.lldscaling.empty() || st.cdscaling) has_scaling = true;
    if (has_scaling) return solver_mvald_analyzer(L, opt, init_sol);

    // 13. an ordinary queueing network. MVAC is a case of the analyzer's method
    // switch (solver_mva_analyzer.m:21-23) and so belongs BELOW branch 12: a
    // load-dependent model reaches solver_mvald_analyzer in the reference and
    // never sees the 'mvac' case, which is why this is not tested higher up.
    if (method == "mvac") {
        DispatchResult<T> mout;
        mout.sol = solver_mvac_analyzer(L, opt);
        mout.actualmethod = "mvac";
        return mout;
    }
    // QNA is dispatched here rather than in solver_mva_analyzer only because
    // solver_qna.h cannot include solver_mva.h without a cycle; the reference
    // reaches it from the analyzer's method table, and no branch above claims
    // 'qna', so the model it sees is the same.
    if (method == "qna") {
        DispatchResult<T> qout;
        qout.sol = solver_qna(L, opt);
        qout.actualmethod = "qna";
        return qout;
    }
    // RQNA (robust queueing-network analyzer) is likewise reached by method name,
    // and by the resolveMethod default->rqna upgrade for a bursty single-class
    // open network; dispatched here for the same include-cycle reason as qna.
    if (method == "rqna") {
        DispatchResult<T> rout;
        rout.sol = solver_rqna(L, opt);
        rout.actualmethod = "rqna";
        return rout;
    }
    // RQT (robust queueing theory) is reached by method name only; its
    // uncertainty-set decomposition has no default-dispatch upgrade.
    if (method == "rqt") {
        DispatchResult<T> tout;
        tout.sol = solver_rqt(L, opt);
        tout.actualmethod = "rqt";
        return tout;
    }
    // The horizontal-cut MVA for one exponential delay and one FCFS MAP queue
    // (mapqn_amva), reached by name; its shape gate is mva_mapqn_reason.
    if (method == "amva.mapqn" || method == "mapqn") {
        DispatchResult<T> mout;
        mout.sol = solver_mapqn(L, opt);
        mout.actualmethod = "amva.mapqn";
        return mout;
    }
    // bursty single-class open network: 'default' upgrades to RQNA, exactly as
    // the terminal branch of solver_mva_analyzer.m does. resolveMethod computes
    // this only for the feature gate and never writes it back to options, so the
    // substitution belongs HERE, after the structural special cases -- a 3-node
    // MAP/M/1 is a G/M/1 the qsys analyzer (branch 4) already solved and must
    // keep method='default'. A fork-transformed model is exempt: it is forced to
    // amva just below.
    if (method == "default" && !opt.base_has_fork && !L.has_fork() && L.nclasses == 1) {
        bool all_open = true;
        for (std::size_t r = 0; r < L.nclasses; ++r)
            if (std::isfinite(L.classes[r].population)) all_open = false;
        if (all_open && api::sn_has_bursty_arrival(L)) {
            DispatchResult<T> rout;
            rout.sol = solver_rqna(L, opt);
            rout.actualmethod = "rqna";
            return rout;
        }
    }

    MvaOptions o = opt;
    if ((opt.base_has_fork || L.has_fork()) && method == "default") {
        // The fork transform yields a mixed model with auxiliary near-zero-rate
        // open classes. The default ladder would send it to exact mixed MVA,
        // which degenerates to zero there; the approximation needs the AMVA.
        //
        // The test is on the BASE model, as the reference writes it: by the
        // time this runs the transform has already removed the fork, so
        // `L.has_fork()` alone is false on exactly the models the rule is for.
        o.method = "amva";
    }
    DispatchResult<T> out;
    out.sol = solver_mva_analyzer(L, o, init_sol);
    out.actualmethod = out.sol.method;
    return out;
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_MVA_DISPATCH_H
