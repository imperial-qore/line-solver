/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_MVA_TYPES_H
#define LINE_SOLVERS_MVA_MVA_TYPES_H

/**
 * The option and result types every MVA analyzer shares.
 *
 * They live in their own header so that an analyzer in a separate translation
 * unit (solver_qna.h, and any analyzer that follows it) can speak the same
 * contract without including solver_mva.h and forming a cycle.
 */

#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/util/matrix.h"

namespace line {
namespace mva {

/** The options SolverMVA reads. Defaults are SolverOptions('MVA'). */
struct MvaOptions {
    std::string method = "default";
    double tol = 1e-4;
    double iter_tol = 1e-6;
    int iter_max = 1000;
    std::string multiserver = "default";
    std::string highvar = "default";
    /**
     * `options.config.np_priority`: which non-preemptive priority correction a
     * HOL station takes. `default` and `cl` are Chandy-Lakshmi, which discounts
     * the higher-priority utilization by the Linearizer's throughput
     * difference; `shadow` is Sevcik's shadow server, which does not.
     */
    std::string np_priority = "default";
    /**
     * `options.config.fork_join`: which fork-join arm the fixed point takes.
     * `default`, `mmt` and `fjt` are the MMT transform of fj_mmt.h;
     * `ht` and `heidelberger-trivedi` are the Heidelberger-Trivedi transform of
     * fj_ht.h, which is closed-model only. Read only by a model with a Fork.
     */
    std::string fork_join = "default";
    /** Warm start, (nstations x nchains) chain-aggregated queue lengths; empty for none. */
    bool has_init_sol = false;
    /**
     * Whether the model this solve came from has a Fork, which the model handed
     * to the analyzer no longer does.
     *
     * `mvaDispatch.m` tests `self.model.hasFork` -- the BASE model -- to force
     * `default` onto AMVA, and by the time the analyzer runs the fork-join
     * transform has already replaced the fork with a router, so the struct in
     * hand cannot answer the question. Testing the transformed struct instead
     * silently sends the model to exact mixed MVA, which the auxiliary
     * near-zero-rate open classes degenerate.
     */
    bool base_has_fork = false;
    /**
     * `options.config.rqt_regime`: which service adaptation regime of the
     * Robust Queueing Theory Table 1 the RQT analyzer takes. Empty and
     * `independent` leave the service distribution unknown; `normal` and
     * `pareto` are the service-dependent fits.
     */
    std::string rqt_regime;
    /** `options.config.rqt_exact`: use the exact worst case over the uncertainty
     *  sets instead of the closed-form bound of Theorem 3. */
    bool rqt_exact = false;
    /** `options.config.rqt_alpha_a`: arrival tail coefficient in (1,2]; 0 = 2. */
    double rqt_alpha_a = 0.0;
    /** `options.config.rqt_alpha_s`: service tail coefficient in (1,2]; 0 = 2. */
    double rqt_alpha_s = 0.0;
    /**
     * `options.config.interlock`: the interlocked-flow matrix of Franks (1999),
     * Eq. (4.7), CLASS-indexed and row-major (nclasses x nclasses). interlock[r][s]
     * is the share of the class-s queue that a class-r arrival must not see, because
     * that work was itself caused by the class-r request. Empty for every model but
     * the layers of SolverLN. It is kept CLASS-indexed, not chain-indexed, so that a
     * later chain refresh cannot leave it stale: the analyzer aggregates it to chains
     * against the struct it is about to solve. Doubles rather than T because the
     * probabilities come from a fixed-point iterate and never enter an exact
     * computation of their own.
     */
    std::vector<std::vector<double>> interlock;
};

/** Class-level results, the [Q,U,R,T,C,X] of the MATLAB analyzers. */
template <class T>
struct MvaSolution {
    Matrix<T> Q, U, R, Tp;
    std::vector<T> C, X;
    std::string method;
    int iter = 0;
    /**
     * Whether the fixed point met its tolerance, or empty when the handler
     * reports none. NOT `iter < iter_max`: `solver_amva` increments `totiter`
     * inside the nested per-chain and inner loops, so it aggregates inner sweeps
     * and can reach the cap on a solve whose outer residual is exactly zero.
     * Only the residual decides there, which is why the flag has to travel with
     * the solution instead of being recomputed by a caller from the count.
     */
    std::optional<bool> converged;
    /**
     * log of the normalizing constant, the reference's `lG`.
     *
     * Only the exact MVA recursion carries one: an AMVA result has no G, and
     * the reference returns 0 there rather than a fabricated value, which is
     * why `getProbNormConstAggr` re-enters the analyzer at method='exact'
     * instead of reusing whatever the last solve produced.
     */
    double lG = 0.0;
};

/**
 * Per-method structural gates, shared by `list_valid_methods` and by the
 * analyzers themselves.
 *
 * ONE predicate with TWO callers. A rule kept in two places is how the report
 * comes to offer a (solver, method) pair that the analyzer then refuses -- or,
 * worse, answers with a table of zeros. What the feature registry CAN name lives
 * in `qn::mva_feature_set` instead; these are the rules it cannot: a product
 * form, a class count and a server count have no registry name.
 *
 * Each returns the refusal sentence, or an empty string when the method may run.
 */

/**
 * The closed-population AMVA family: an open chain, or a network outside strict
 * product form, has nothing for the recursion to work on.
 *
 * Open chains are also expressed in the registry (`mva_feature_set` drops
 * OpenClass for these), which is what keeps them off the report; they are
 * repeated here because the analyzer must refuse by name with a sentence rather
 * than fall through to `solver_amvald` and answer under a method nobody asked
 * for. Strict product form has no registry name at all, so this is its only home.
 */
template <class T>
std::string mva_closed_population_reason(const qn::NetworkStruct<T>& L,
                                         const std::string& method) {
    if (!qn::mva_is_closed_population_method(method)) return "";
    const std::string base = qn::mva_base_method(method);
    if (L.has_open_classes())
        return "solver_amva: the '" + base +
               "' method approximates the arrival-instant queue length as a function of the "
               "closed population vector N, so it is defined for closed models only; use "
               "'default', 'lin', 'qd' or 'qna' for a model with open classes";
    // ab, schmidt and schmidt-ext ARE the class-dependent FCFS algorithms, so
    // heterogeneous FCFS service means are their subject matter rather than a
    // disqualification.
    const bool check_means = !(base == "ab" || base == "schmidt" || base == "schmidt-ext");
    if (L.has_product_form_not_het_fcfs(check_means)) return "";
    return "solver_amva: the '" + base +
           "' method is defined for strict product-form, load-independent models; use 'default', "
           "'lin' or 'qd' for this model";
}

/**
 * RQNA and RQT decompose an open network into GI/G/1 queues and build one
 * uncertainty set per flow out of the first two moments of a SINGLE stream, so a
 * multiclass model has no counterpart in their equations. No registry feature
 * names a class count, so that half of the rule is structural.
 *
 * A FORK-JOIN model is refused too. A Join is a synchronisation node, not a
 * queue: it carries no service process, so the index-of-dispersion curve these
 * analyzers read off every station does not exist for it, and neither has a
 * synchronisation term to put in its place. That half IS nameable, so
 * `qn::mva_feature_set` drops Fork/Join for these two as well and this is the
 * analyzer's half of it.
 */
template <class T>
std::string mva_single_class_open_reason(const qn::NetworkStruct<T>& L,
                                         const std::string& method) {
    const std::string base = qn::mva_base_method(method);
    if (base != "rqna" && base != "rqt") return "";
    const std::string label = (base == "rqna") ? "RQNA" : "RQT";
    // Scanned over the NODES and not the stations: a Fork is not a station in
    // this port, so a station scan would see only the Join.
    for (const qn::NodeDef& nd : L.nodes)
        if (nd.nodetype == qn::NodeType::Fork || nd.nodetype == qn::NodeType::Join)
            return "solver_" + base + ": " + label +
                   " decomposes an open network into GI/G/1 queues and has no synchronisation "
                   "term; a Join carries no service process for its index of dispersion to be "
                   "read from; use the 'default' method for a fork-join model";
    if (L.nclasses == 1) return "";
    return "solver_" + base + ": " + label +
           " supports single-class open networks only; use method 'qna' for multiclass models";
}

/**
 * The extended Schmidt method needs a customer of every class to tag.
 *
 * Schmidt's EXTENSION over plain Schmidt is an alpha correction applied at an
 * FCFS station, computed from the network with ONE class-r customer TAGGED, that
 * is at population N - 1_r. A class holding no customer has none to tag: the
 * sub-problem is formed at a negative population, whose state lattice prod(N+1)
 * collapses to zero and the recursion indexes an empty array. Plain `schmidt`
 * forms no such sub-problem, which is why the requirement is the -ext arm's
 * alone.
 *
 * THE TEST IS STATED AT THE FCFS STATION AND NOT AT A CLASS-DEPENDENT ONE,
 * because the four kernels differ on when they form the correction: MATLAB, this
 * port and native python form it only where the station's demands differ by
 * class, the JAR forms it at every FCFS station. Stating the union is what keeps
 * one rule safe for all four; the case it costs -- an FCFS station whose demands
 * are identical across classes, one of them empty -- is one where the extension
 * reduces to plain `schmidt`, which stays offered.
 *
 * `njobs` and `fcfs` are the numbers the CALLER'S OWN arm passes: CHAIN-indexed
 * here as in MATLAB, CLASS-indexed in the JAR and native python, whose arms
 * aggregate no chains. That difference belongs to those arms, not to this rule.
 */
inline std::string mva_schmidt_ext_reason(const std::vector<double>& njobs,
                                          const std::vector<bool>& fcfs,
                                          const std::string& method) {
    if (qn::mva_base_method(method) != "schmidt-ext") return "";
    bool any_fcfs = false;
    for (std::size_t i = 0; i < fcfs.size(); ++i)
        if (fcfs[i]) any_fcfs = true;
    if (!any_fcfs) return "";
    for (std::size_t r = 0; r < njobs.size(); ++r)
        if (std::isfinite(njobs[r]) && njobs[r] < 1.0)
            return "solver_amva: the 'schmidt-ext' method corrects an FCFS station from the "
                   "network with one customer of that class tagged, so it needs every class to "
                   "hold at least one customer; class " +
                   std::to_string(r + 1) +
                   " holds none; use 'schmidt' for the uncorrected recursion";
    return "";
}

/**
 * MVAC is the exact chain recursion over single-server fixed-rate (SSFR) queues
 * and infinite-server centres of a product-form network, and it recurs on the
 * queueing centres, so it needs at least one. Neither the server count nor
 * product form has a registry feature name; the scheduling restriction IS
 * nameable and lives in `qn::mva_feature_set`.
 */
template <class T>
std::string mva_mvac_reason(const qn::NetworkStruct<T>& L, const std::string& method) {
    if (qn::mva_base_method(method) != "mvac") return "";
    if (!L.has_product_form()) return "solver_mvac_analyzer: MVAC requires a product-form model";
    std::size_t nq = 0;
    for (const qn::Station<T>& st : L.stations) {
        if (st.sched == qn::SchedStrategy::INF || st.sched == qn::SchedStrategy::EXT) continue;
        ++nq;
        // An infinite count is refused here too, as the reference does: a station
        // scheduled FCFS with infinitely many servers is not an IS centre to MVAC,
        // and `infSET` below is built from the DISCIPLINE, not the count.
        if (st.nservers != 1.0)
            return "solver_mvac_analyzer: MVAC supports single-server (SSFR) queues only; use "
                   "method 'exact' for multiserver stations";
    }
    if (nq == 0)
        return "solver_mvac_analyzer: MVAC recurs on the queueing centers and needs at least one; "
               "this model has only delay stations";
    return "";
}

/**
 * QNA's station update has an arm for INF, PS and FCFS and none for any other
 * discipline, so a SIRO, LCFS, LCFS-PR, HOL or priority station used to leave
 * its whole row of Q, U, R and T at zero and the table was returned as a
 * solution. The registry expresses this too (`mva_feature_set` drops the
 * disciplines from QNA's envelope); this is the analyzer's half of it.
 */
template <class T>
std::string mva_qna_scheduling_reason(const qn::NetworkStruct<T>& L) {
    for (std::size_t i = 0; i < L.stations.size(); ++i) {
        // A Join station carries no service and is skipped by the update below.
        if (L.stations[i].nodetype == qn::NodeType::Join) continue;
        const qn::SchedStrategy sched = L.stations[i].sched;
        if (sched == qn::SchedStrategy::EXT || sched == qn::SchedStrategy::INF ||
            sched == qn::SchedStrategy::PS || sched == qn::SchedStrategy::FCFS)
            continue;
        return std::string("solver_qna: QNA decomposes every station as a GI/G/m centre and has "
                           "no arm for ") +
               lang::sched_to_text(sched) + " scheduling; use the 'default' or 'lin' methods";
    }
    return "";
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_MVA_TYPES_H
