/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_BA_SOLVER_BA_RUNNER_H
#define LINE_SOLVERS_BA_SOLVER_BA_RUNNER_H

/**
 * The SolverBA class surface: `@@SolverBA/runAnalyzer.m`, `listValidMethods`,
 * `getBounds` and `getBoundsTable`.
 *
 * What sits here rather than in the analyzer is everything around one bound
 * evaluation: the closed-model gate, the method aliases, the whitelist, the
 * arrival-rate conversion and the metric filter `@@NetworkSolver/getAvg` applies
 * to every solver's output. The result shape is SolverMVA's `AvgResult`, so a
 * bound and a point estimate are directly comparable -- which is what a caller
 * checking that a bracket contains the exact answer needs.
 *
 * FINITE-BUFFER BLOCKING IS REFUSED, NOT BOUNDED. Needing only demands and a
 * population is the BCMP parameterization, which presumes UNBOUNDED buffers; a
 * buffer that binds couples the station occupancies and the resulting numbers
 * do not bracket the blocked model. `solver_ba_run_analyzer` therefore gates on
 * `api::sn_has_blocking` and `list_valid_methods(L)` drops every blocking-blind
 * method. The exceptions are `qrf.bas*`/`qrf.rsrd`, which carry the blocking
 * tables explicitly. Use SolverMVA method `sqd` for a point estimate.
 *
 * THE RESIDENCE TIME IS ZERO BY CONSTRUCTION. `runAnalyzer.m` sets `WN` to a
 * zero matrix and does not derive it from `RN`: the bound families define no
 * per-station residence, and the reference declines to invent one.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "line/api/sn/sn_predicates.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/solvers/ba/solver_ba_analyzer.h"
#include "line/api/sn/sn_to_qrf_blocking.h"
#include "line/api/sn/sn_to_qrf_alpha.h"
#include "line/solvers/ba/solver_ba_qrf_analyzer.h"
#include "line/solvers/ba/solver_ba_spnlp.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace line {
namespace ba {

using lang::GlobalConstants;

/**
 * Port of `runAnalyzer`'s method aliases: `default` is the geometric upper
 * bound, bare `auto` is the AUTO composite's upper side, and `qr`/`lr` are the
 * friendly names of two QRF sub-methods.
 */
inline std::string resolve_method(const std::string& method) {
    if (method == "default") return "gb.upper";
    if (method == "auto") return "auto.upper";
    if (method == "lr") return "lr.upper";
    if (method == "qr") return "qrf.mmi";
    return method;
}

/**
 * Whether METHOD bounds a model as if its buffers were unbounded.
 *
 * Every family but the QRF BLOCKING bounds is parameterized by demands (visits
 * x service time) and a population alone, which is the BCMP parameterization:
 * unbounded buffers, and an equilibrium distribution that factorizes. A finite
 * buffer that BINDS breaks both premises, so the numbers do not bracket the
 * blocked model -- on `cqn_bas_blocking` (Queue2 capped at 1, N = 2) `gb.upper`
 * reports QLen 1.28 at a station that can never hold more than one job.
 * `qrf.bas*` and `qrf.rsrd` carry the blocking tables (MM, MM1, ZZ, ZM, BB, F)
 * explicitly and are the exceptions.
 *
 * METHOD must already be resolved through `resolve_method`.
 */
inline bool ignores_blocking(const std::string& method) {
    // 'spnlp.*' is an exception alongside the QRF blocking bounds: its polytope
    // is indexed by the marking itself, so a bounded place enters it as a
    // variable upper bound and as the P-invariant equality that produced the
    // bound. The buffer is modelled, not assumed away.
    return !(method.rfind("qrf.bas", 0) == 0 || method.rfind("qrf.rsrd", 0) == 0 ||
             method.rfind("spnlp", 0) == 0);
}

/**
 * Port of `SolverBA.listValidMethods`.
 *
 * A LISTED NAME MUST ACTUALLY RUN, OR BE REFUSED BY NAME. The `lr` family IS
 * listed: its bound is a pure linear program, which `lp::simplex_solve` serves
 * exactly. The four no-blocking `qrf.*` arms and the supported BAS-blocking
 * bounds are listed too, served by `solver_ba_qrf_analyzer`.
 *
 * `qrf.bas.mmi` was the one name here that used to be always refused, and it is listed
 * because the reference lists it: `SolverBA.listAllMethods` names it and
 * `solver_ba_qrf_analyzer.m:119` then refuses it by name, having lost
 * `qrf_bas_mmi_simple.m` on 2026-08-01. Omitting it here made the C++ answer
 * "unsupported method" where the reference answers with the analyzer's own
 * reason, so a caller could not tell a missing port from a retired bound.
 */
inline std::vector<std::string> list_valid_methods() {
    return {"default",
            "auto.upper", "auto.lower",
            "aba.upper",  "aba.lower",
            "bjb.upper",  "bjb.lower",
            "pb.upper",   "pb.lower",
            "gb.upper",   "gb.lower",
            "sb.upper",   "sb.lower",
            "harel.upper", "harel.lower",
            "mwba.upper", "mwba.lower",
            "pbh.upper",  "pbh.lower",
            "pbk.upper",  "pbk.lower",
            "bjbk.upper", "bjbk.lower",
            "cbh.upper",  "cbh.lower",
            "ssd.upper",  "ssd.lower",
            "cub.upper",  "mbjb.lower",
            "looping.upper", "looping.lower",
            "bpt.lower", "bgt.upper", "snc.upper",
            "sib.upper",  "sib.lower",
            "scb.upper",  "scb.lower",
            "ldbcmp.lower",
            "lr",         "lr.upper",   "lr.lower",
            "mapamva.upper", "mapamva.lower",
            "qr",         "qrf.mmi",    "qrf.mem",    "qrf.bethe",
            "qrf.mmi.ld", "qrf.mmi.linear",
            "qrf.bas",    "qrf.bas.mmi", "qrf.bas.mem", "qrf.bas.bethe",
            "qrf.rsrd",
            "spnlp.upper", "spnlp.lower", "spnlp.op.upper", "spnlp.op.lower"};
}

/**
 * The same list, narrowed to what THIS model can run.
 *
 * The QR/LR/QRF reduction bounds share one premise -- a single-class closed
 * network of single-server stations -- which `solver_ba_qrf_analyzer` enforces
 * and the `lr` family shares. A listed name that always throws is a method a
 * caller is invited to ask for and cannot have, so the model-aware overload
 * drops them; `SolverMVA` gates `sqni` the same way. Mirrors `SolverBA.m`,
 * `SolverBA.java` and the python `list_valid_methods`.
 */
/**
 * What `default`/`auto` must mean on a model with finite buffers: the QRF BAS
 * bound, or nothing.
 *
 * `default` resolves to the geometric upper bound, which is parameterized by
 * demands and a population alone and therefore bounds a blocked model as if its
 * buffers were unbounded. `ignores_blocking` refuses that, which is right; but
 * refusing is not the whole answer, because `qrf.bas` DOES model the finite
 * buffer and `sn_to_qrf_blocking` derives its tables from the model, so there
 * is nothing left for the caller to supply. A blocked model of the right shape
 * therefore gets `qrf.bas` as its default, exactly as SolverMVA routes a BAS
 * model to `sqd`.
 *
 * The shape is the one `list_valid_methods(L)` calls "reducible" and the QRF
 * analyzer gates on. On top of it the tables must actually derive, which is
 * asked of `sn_to_qrf_blocking` rather than re-tested here -- it owns the
 * single-finite-buffer rule and the size guard, and a second copy of either is
 * how the two drift apart.
 *
 * Only the UPPER side is routed: the analyzer solves qrf.bas in the 'max'
 * direction alone, so `auto.lower` has no blocking counterpart and keeps
 * refusing rather than being answered with the wrong side.
 *
 * @return {method, why}; method is empty when the routing does not apply, and
 *         why then carries the reason (empty for an unblocked model).
 */
template <class T>
std::pair<std::string, std::string> blocking_default(const qn::NetworkStruct<T>& L) {
    typedef std::pair<std::string, std::string> R;
    if (!api::sn_has_blocking(L)) return R("", "");
    if (L.nclasses != 1)
        return R("", "the QRF blocking bounds are derived for a single-class closed network,"
                     " which this model is not.");
    for (std::size_t r = 0; r < L.classes.size(); ++r)
        if (!std::isfinite(L.classes[r].population))
            return R("", "the QRF blocking bounds are derived for a single-class closed"
                         " network, which this model is not.");
    for (std::size_t i = 0; i < L.stations.size(); ++i) {
        if (L.stations[i].sched == SchedStrategy::INF)
            return R("", "the QRF blocking bounds model every station as a single server and"
                         " have no infinite-server notion, so a delay station rules them out.");
        if (L.stations[i].nservers > 1.0)
            return R("", "the QRF blocking bounds model every station as a single server, so a"
                         " multiserver station rules them out.");
    }
    // Same phase count `solver_ba_qrf_analyzer` derives, so the size guard
    // here judges the LP the analyzer would actually build.
    int Ktot = 0;
    for (std::size_t i = 0; i < L.stations.size(); ++i) {
        const mam::Map<T> mp = lang::dist_to_map(L.service[i][0]);
        Ktot += mp.D0.rows() == 0 ? 1 : static_cast<int>(mp.D0.rows());
    }
    const sn::QrfBlocking blk = sn::sn_to_qrf_blocking(L, Ktot);
    if (!blk.msg.empty()) return R("", blk.msg);
    return R("qrf.bas", "");
}

/** True for the aliases whose meaning a blocked model is allowed to change. */
inline bool is_default_request(const std::string& method) {
    return method == "default" || method == "auto" || method == "auto.upper";
}

template <class T>
std::vector<std::string> list_valid_methods(const qn::NetworkStruct<T>& L) {
    std::vector<std::string> all = list_valid_methods();
    // The STRUCTURAL premises -- single-class closed, fully closed,
    // single-server -- come from `method_refusal`, the same predicate
    // `solver_ba_run_analyzer` throws on. Asked first so every later narrowing
    // works on names this model could actually run, and asked HERE because this
    // is the list `autosolver::auto_family_methods` reads: before it, a
    // two-class closed network was offered all 36 demand-parameterized bounds
    // and 30 of them threw the moment they were run.
    // ... and `method_degenerate` withholds the second kind of name: one whose
    // premises this model MEETS but whose formula says nothing here.
    // 'ldbcmp.lower' at N == Qhat is the only such case: it reports the trivial
    // X >= 0, which propagates into an all-zero table a caller cannot tell from
    // an answer. Offering is what stops; asking for it by NAME still runs and
    // still publishes it, since a vacuous bound is a valid one -- which is what
    // `tests/test_ba.cpp` pins at that boundary.
    {
        std::vector<std::string> structural;
        for (std::size_t i = 0; i < all.size(); ++i)
            if (method_refusal(L, all[i]).empty() && method_degenerate(L, all[i]).empty())
                structural.push_back(all[i]);
        all = structural;
    }
    // 'spnlp.*' is the only family indexed by a MARKING rather than by demands
    // and a population, and the split is total in both directions: on a Petri
    // net nothing else has a representation of the model, and off one spnlp has
    // nothing to read. Two of the gates below already half-cover this by
    // accident -- a Place is an INF station, so `reducible` and `bptOk` are both
    // false on any Petri net -- but the demand-parameterized families survive
    // them and must be dropped by name. Applied first so the open/closed
    // narrowing cannot reinstate one.
    bool isPetri = false;
    for (std::size_t i = 0; i < L.nodes.size() && !isPetri; ++i)
        isPetri = L.nodes[i].nodetype == lang::NodeType::Transition;
    {
        std::vector<std::string> sieved;
        for (std::size_t i = 0; i < all.size(); ++i)
            if (is_spnlp_method(all[i]) == isPetri) sieved.push_back(all[i]);
        all = sieved;
        if (isPetri) return all;
    }
    bool reducible = (L.nclasses == 1);
    for (std::size_t r = 0; reducible && r < L.classes.size(); ++r)
        if (!std::isfinite(L.classes[r].population)) reducible = false;
    bool closedSingleClass = reducible;
    for (std::size_t i = 0; reducible && i < L.stations.size(); ++i) {
        if (L.stations[i].sched == SchedStrategy::INF)
            reducible = false;
        else if (L.stations[i].nservers > 1.0)
            reducible = false;
    }
    // The LOAD-DEPENDENT arms survive where the rest of the family cannot run:
    // alpha(i,n) is the rate law of a delay (alpha = n), of a c-server station
    // (alpha = min(n,c)) and of limited load dependence alike, so
    // `qrf.mmi.ld` and `qrf.mmi.linear` answer those models on the model's own
    // chain. `sn_to_qrf_alpha` owns the one restriction that survives,
    // exponential service wherever a station serves several jobs at once.
    // Dropping them with the rest would hide from a caller enumerating the list
    // the only two bound methods such a model has.
    const bool ldReducible =
        !reducible && closedSingleClass && sn::sn_to_qrf_alpha(L).msg.empty();
    // 'bpt', 'bgt' and 'snc' are the mirror image of the reduction bounds: all
    // three are derived for an OPEN network of single-server exponential
    // stations, so every closed model, every delay station and every
    // multiserver station rules them out. Every other family rules OUT the open
    // model, so on an open network the list narrows to those three.
    bool fullyOpen = !L.classes.empty();
    for (std::size_t r = 0; fullyOpen && r < L.classes.size(); ++r)
        if (std::isfinite(L.classes[r].population)) fullyOpen = false;
    bool bptOk = fullyOpen;
    for (std::size_t i = 0; bptOk && i < L.stations.size(); ++i) {
        if (L.stations[i].nodetype == qn::NodeType::Source) continue;
        if (L.stations[i].sched == SchedStrategy::INF)
            bptOk = false;
        else if (L.stations[i].nservers > 1.0)
            bptOk = false;
    }
    std::vector<std::string> keep;
    if (fullyOpen) {
        if (bptOk) {
            // Each is re-added only if it SURVIVED the structural narrowing
            // above. Pushing them unconditionally is what silently undid it:
            // `snc.upper` is dropped there on non-exponential service, and a
            // blind push put it straight back, so an Erlang-service model was
            // still offered a bound that throws.
            const char* const open_fams[] = {"bpt.lower",
                                             // 'bgt.upper' additionally needs
                                             // deterministic non-merging routes, and
                                             // 'snc.upper' a feed-forward station
                                             // graph; both are checked by walking the
                                             // routing matrix in their own analyzers,
                                             // too expensive to repeat here, so they
                                             // stay listed and refuse by name.
                                             "bgt.upper", "snc.upper"};
            for (const char* nm : open_fams)
                if (std::find(all.begin(), all.end(), std::string(nm)) != all.end())
                    keep.push_back(nm);
        }
    } else {
        for (std::size_t i = 0; i < all.size(); ++i) {
            const std::string& m = all[i];
            if (m == "bpt.lower" || m == "bgt.upper" || m == "snc.upper") continue;
            if (reducible) {
                keep.push_back(m);
                continue;
            }
            if (ldReducible && (m == "qrf.mmi.ld" || m == "qrf.mmi.linear")) {
                keep.push_back(m);
                continue;
            }
            // 'mapamva' shares that premise exactly -- single-class closed,
            // single-server, no delay -- so it narrows with them rather than
            // being offered on a model it would refuse on contact. It is NOT
            // one of the load-dependent arms: its q carries no population
            // index, so a delay or a c-server station has nowhere to go.
            if (m == "qr" || m == "lr" || m.compare(0, 3, "lr.") == 0 ||
                m.compare(0, 4, "qrf.") == 0 || m.compare(0, 7, "mapamva") == 0)
                continue;
            keep.push_back(m);
        }
    }
    // A binding finite buffer rules out everything but the QRF blocking bounds:
    // the other families presume unbounded buffers, and `solver_ba_run_analyzer` refuses
    // them by name on such a model. The list can legitimately come back EMPTY
    // -- a blocked model that is not single-class closed single-server has no
    // bound method at all, and offering one would be the mis-selection this
    // gate exists to prevent.
    if (api::sn_has_blocking(L)) {
        std::vector<std::string> unblocked;
        // 'default' is offered back when it now MEANS one of the survivors:
        // `solver_ba_run_analyzer` routes it to 'qrf.bas' on a blocked model of the right
        // shape, so a caller enumerating the list would otherwise be told the
        // model's own default is invalid.
        if (!blocking_default(L).first.empty()) unblocked.push_back("default");
        for (std::size_t i = 0; i < keep.size(); ++i)
            if (!ignores_blocking(resolve_method(keep[i]))) unblocked.push_back(keep[i]);
        return unblocked;
    }
    return keep;
}

/** Port of `runAnalyzer`'s method gate: an unlisted name is refused by name. */
inline void check_method(const std::string& method) {
    // `qrf.bas.mmi` used to be refused here, mirroring the reference. The
    // reference serves it again since 2026-09-03, so nothing is refused by name
    // beyond an unlisted one.
    const std::string m = resolve_method(method);
    const std::vector<std::string> valid = list_valid_methods();
    if (std::find(valid.begin(), valid.end(), m) != valid.end()) return;
    throw UnsupportedError("SolverBA: unknown bound method '" + method + "'");
}

/**
 * Port of `@@SolverBA/runAnalyzer.m` for the `lang='matlab'` path.
 *
 * @param L   the refreshed struct of a closed model
 * @param opt_in the method and, for a hierarchical family, the level
 */
template <class T>
mva::AvgResult<T> solver_ba_run_analyzer(const qn::NetworkStruct<T>& L, const BaOptions& opt_in) {
    // Every bound family here is parameterized by demands and a closed
    // population, except the three OPEN-network families 'bpt', 'bgt' and
    // 'snc', which are refused on a closed model instead.
    // ... and the 'spnlp' family, which is parameterized by a MARKING: whether
    // that marking is bounded is a question about the P-invariants of the net
    // and not about nclosedjobs, so `spn_lpbnd` decides it rather than this gate.
    if (!(L.nclosedjobs() > 0.0) && opt_in.method.compare(0, 3, "bpt") != 0 &&
        opt_in.method.compare(0, 3, "bgt") != 0 && opt_in.method.compare(0, 3, "snc") != 0 &&
        !is_spnlp_method(opt_in.method))
        throw UnsupportedError("SolverBA: supports closed queueing networks only");
    check_method(opt_in.method);
    BaOptions opt = opt_in;
    opt.method = resolve_method(opt_in.method);

    // Finite-buffer BLOCKING is outside the premises of every family here
    // except the QRF blocking bounds: the rest are parameterized by demands and
    // a population alone, which presumes unbounded buffers and a product form
    // that the truncation destroys. Refusing is not conservatism -- on
    // `cqn_bas_blocking`, `gb.upper` reports QLen 1.28 at a station capped at 1
    // job. Gated after the aliases so `default` is judged as the `gb.upper` it
    // resolves to.
    std::string blocking_why;
    if (ignores_blocking(opt.method) && api::sn_has_blocking(L) &&
        is_default_request(opt_in.method)) {
        // A blocked model whose shape admits the QRF BAS bound gets it as the
        // DEFAULT rather than a refusal: 'qrf.bas' models the finite buffer,
        // and since `sn_to_qrf_blocking` derives its tables from the model
        // there is nothing left for the caller to supply.
        const std::pair<std::string, std::string> routed = blocking_default(L);
        blocking_why = routed.second;
        if (!routed.first.empty()) opt.method = routed.first;
    }

    if (ignores_blocking(opt.method) && api::sn_has_blocking(L))
        throw UnsupportedError(
            "SolverBA: method '" + opt_in.method +
            "' does not support finite-buffer blocking: every SolverBA bound family but the"
            " QRF blocking ones is parameterized by demands and a population alone, so it"
            " bounds the model as if its buffers were unbounded. Use SolverMVA with method"
            " 'sqd', an exact solver (CTMC, SSA, JMT, LDES), or the QRF blocking bounds"
            " 'qrf.bas'/'qrf.rsrd', which model the finite buffer" +
            (blocking_why.empty()
                 ? std::string()
                 : ". The QRF blocking bounds do not apply here either: " + blocking_why));

    // The STRUCTURAL premises -- single-class closed, fully closed,
    // single-server -- are NOT re-asked here: `solver_ba_analyzer` asks
    // `method_refusal` on the way in and throws it, so a second ask would only
    // duplicate the throw one frame earlier and lose the analyzer's own prefix.
    // The report is gated by `list_valid_methods(L)`, which projects the same
    // predicate plus `method_degenerate`.

    const bool qrf = detail::is_qrf_noblo_method(opt.method) ||
                     detail::is_qrf_lp_method(opt.method) ||
                     detail::is_qrf_bas_nlp_method(opt.method);
    const BaSolution<T> s = is_spnlp_method(opt.method)
                                ? solver_ba_spnlp_analyzer(L, opt)
                                : (qrf ? solver_ba_qrf_analyzer(L, opt) : solver_ba_analyzer(L, opt));
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = L.nstations, K = L.nclasses;

    std::vector<std::vector<bool>> mask(M, std::vector<bool>(K, false));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k)
            mask[i][k] = num_traits<T>::to_double(s.R(i, k)) < 10.0 * GlobalConstants::FineTol;

    mva::AvgResult<T> out;
    out.QN = mva::filter_metric(L, s.Q, mva::MetricKind::QLen, &mask);
    out.UN = mva::filter_metric(L, s.U, mva::MetricKind::Util, &mask);
    out.RN = mva::filter_metric(L, s.R, mva::MetricKind::RespT, nullptr);
    out.TN = mva::filter_metric(L, s.Tp, mva::MetricKind::Tput, nullptr);
    // The arrival rate comes from the throughputs through the same helper
    // SolverMVA and SolverNC use. The bound goldens were recorded under those
    // solvers before the family moved here and carry a real per station-class
    // arrival rate; returning zeros silently emptied one column. The Source mask
    // is getAvg's and is kept although no model that reaches here has a Source:
    // the gates admit a mixed model, and every family then refuses it.
    std::vector<std::vector<bool>> srcmask(M, std::vector<bool>(K, false));
    for (std::size_t i = 0; i < M; ++i)
        if (L.stations[i].nodetype == qn::NodeType::Source)
            for (std::size_t k = 0; k < K; ++k) srcmask[i][k] = true;
    out.AN = mva::filter_metric(L, mva::sn_get_arvr_from_tput(L, out.TN), mva::MetricKind::ArvR,
                                &srcmask);
    out.WN = Matrix<T>(M, K, zero);
    out.CN = s.C;
    out.XN = s.X;
    out.method = opt_in.method;
    out.actualmethod = opt.method;
    out.iter = s.iter;
    return out;
}

/** Port of `SolverBA.getBounds`: the {lower,upper} bracket of a family. */
template <class T>
struct BaBounds {
    Matrix<T> Qlower, Qupper;  ///< (M x K), all-NaN on a side the family lacks
    Matrix<T> Tlower, Tupper;
    bool has_lower = false, has_upper = false;
    /**
     * `getBoundsTable`'s row filter, (M x K): whether the (station, class) pair
     * earns a row.
     *
     * This is the only part of `getBoundsTable` that is behaviour rather than
     * formatting, so it is carried here and the label columns are not ported --
     * no solver in `cpp/` has a table layer, `getAvgTable` included, and the
     * CLI formats its own. The rule mirrors `getAvgTable`'s drop of disabled
     * pairs but is NaN-SAFE: a row survives when any value that is present is
     * nonzero, so an all-NaN side (a one-sided family) never removes it.
     */
    std::vector<std::vector<bool>> keep;
};

/**
 * Port of `SolverBA.getBounds`.
 *
 * The family is the method's prefix before the first dot, and BOTH sides are
 * re-run with the caller's full option set: constructing the re-run with only
 * the method would reset `level` to its default and a hierarchical family would
 * never tighten as the level is raised. A one-sided family (cub upper-only,
 * mbjb/ldbcmp lower-only) leaves its missing side as NaN, never as zero.
 */
template <class T>
BaBounds<T> ba_bounds(const qn::NetworkStruct<T>& L, const BaOptions& opt) {
    const std::string m = resolve_method(opt.method);
    const std::string fam = m.substr(0, m.find('.'));
    // `solver_ba_run_analyzer` routes 'default'/'auto' to 'qrf.bas' on a blocked model,
    // and ba_bounds has to say so rather than let the re-run below fail with
    // the generic gate -- but it still cannot BRACKET, because the analyzer
    // solves qrf.bas in the 'max' direction alone.
    if (ignores_blocking(m) && api::sn_has_blocking(L) && is_default_request(opt.method)) {
        const std::pair<std::string, std::string> routed = blocking_default(L);
        if (!routed.first.empty())
            throw UnsupportedError(
                "SolverBA: '" + opt.method + "' resolves to '" + routed.first +
                "' on this model, which has a binding finite buffer, and that bound is"
                " UPPER-only: there is no bracket to return. Call getAvg for the upper"
                " bound, or SolverMVA with method 'sqd' for a point estimate");
    }
    const std::vector<std::string> valid = list_valid_methods();
    auto listed = [&](const std::string& x) {
        return std::find(valid.begin(), valid.end(), x) != valid.end();
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();
    BaBounds<T> b;
    const std::size_t M = L.nstations, K = L.nclasses;
    const T nanT = num_traits<T>::from_double(nan);
    b.Qlower = Matrix<T>(M, K, nanT);
    b.Qupper = Matrix<T>(M, K, nanT);
    b.Tlower = Matrix<T>(M, K, nanT);
    b.Tupper = Matrix<T>(M, K, nanT);
    if (listed(fam + ".lower")) {
        BaOptions o = opt;
        o.method = fam + ".lower";
        const mva::AvgResult<T> r = solver_ba_run_analyzer(L, o);
        b.Qlower = r.QN;
        b.Tlower = r.TN;
        b.has_lower = true;
    }
    if (listed(fam + ".upper")) {
        BaOptions o = opt;
        o.method = fam + ".upper";
        const mva::AvgResult<T> r = solver_ba_run_analyzer(L, o);
        b.Qupper = r.QN;
        b.Tupper = r.TN;
        b.has_upper = true;
    }
    b.keep.assign(M, std::vector<bool>(K, true));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k) {
            const T v[4] = {b.Qlower(i, k), b.Qupper(i, k), b.Tlower(i, k), b.Tupper(i, k)};
            bool any_present = false, any_nonzero = false;
            for (const T& x : v) {
                if (std::isnan(num_traits<T>::to_double(x))) continue;
                any_present = true;
                if (x != num_traits<T>::from_int(0)) any_nonzero = true;
            }
            b.keep[i][k] = !any_present || any_nonzero;
        }
    return b;
}

}  // namespace ba
}  // namespace line

#endif  // LINE_SOLVERS_BA_SOLVER_BA_RUNNER_H
