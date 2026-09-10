/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NC_RUNNER_H
#define LINE_SOLVERS_NC_SOLVER_NC_RUNNER_H

/**
 * The SolverNC class surface: `@@SolverNC/runAnalyzer.m` and the gates around it.
 *
 * What sits here rather than in the dispatch is everything that happens BEFORE
 * and AFTER one inner solve: the method whitelist, the structural gates, the
 * MULTISERVER-TO-LOAD-DEPENDENCE conversion, the conversions from response time
 * to residence time and from throughput to arrival rate, and the metric filter.
 *
 * THE MULTISERVER CONVERSION IS THE INTERESTING PART, and it belongs here and
 * not in an analyzer. A c-server station is rewritten as the rate lattice
 * mu(n) = min(n, c), which routes the model to `solver_ncld` and is EXACT,
 * whereas leaving it alone routes it to `solver_nc` and Seidmann's
 * approximation. The reference does the rewrite on 'exact' and 'is' always, and
 * on 'default' only for the two-station Delay-plus-multiserver shape (which is
 * every SolverLN layer submodel), and only when the model is product-form --
 * a non-product-form model has no exact load-dependent solution, so it is sent
 * to 'comom' instead. The server count is deliberately KEPT: utilization is the
 * fraction of the c servers busy, and c is not recoverable from min(1:Nt, c)
 * once the population is below it.
 *
 * The metric filter is shared with SolverMVA verbatim (`filter_metric`,
 * `sn_get_residt_from_respt`, `sn_get_arvr_from_tput` in `solver_mva_runner.h`)
 * because `@@NetworkSolver/getAvg` is solver-independent: it is the same code
 * path for both solvers in the reference.
 */

#include "line/util/line_console.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/solvers/cache_metrics.h"
#include "line/solvers/mva/fj_driver.h"
#include "line/solvers/mva/fj_ht.h"
#include "line/solvers/mva/fj_mmt.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/nc/nc_dispatch.h"
#include "line/solvers/nc/solver_nc_cache.h"
#include "line/solvers/nc/solver_nc_cacheqn.h"
#include "line/solvers/nc/solver_nc_cacheqn_retrieval.h"
#include "line/solvers/nc/solver_nc_lossn.h"
#include "line/solvers/nc/solver_nc_spn.h"
#include "line/solvers/nc/solver_nc_mem.h"
#include "line/solvers/nc/solver_nc_oi.h"
#include "line/solvers/nc/solver_nc_dps.h"
#include "line/solvers/nc/solver_nc_retrieval.h"
#include "line/solvers/nc/nc_types.h"
#include "line/util/error.h"

namespace line {
namespace nc {

/**
 * Port of `SolverNC.listValidMethods`.
 *
 * A LISTED NAME MUST ACTUALLY RUN, or be refused with a message that names the
 * missing analyzer. `erlangfp` and `mci` (the loss-network Erlang fixed point
 * and its Monte Carlo counterpart) DO run: a Finite Capacity Region is now
 * representable -- `sn.regions` and `Network::add_region` -- and
 * `solver_nc_lossn_analyzer` claims the model above the product-form and
 * multiserver gates. A model that is not a loss network still reaches the
 * refusal further down, which is the honest outcome and different from silently
 * solving a model with its region ignored. `comomld` is not in the reference's
 * list; it is accepted because the reference selects it internally from
 * 'default'. `ms` names the same lossn_manjunath transform as `exact` on a loss
 * network (solver_nc_lossn.h's own method name alias) and was missing here, which
 * blocked it before it ever reached the analyzer that already accepts it.
 */
inline std::vector<std::string> list_valid_methods() {
    // "rayint" and "spm" both name the SPM saddle point on a cache, which serves
    // cache_spm_size once the items carry storage costs. On a retrieval model
    // "rayint" is instead the ray/WKB delayed-hit expansion, admissible only with
    // an infinite-server fetch system; solver_nc_retrieval branches on the method name
    // and warns and falls back to "exact" anywhere else.
    return {"default", "exact",  "rayint",  "spm",  "erlangfp", "mci",      "imci", "ls",   "le",
            "ble",    "aghq",   "mmint2",  "gleint", "pana",  "panald", "ca",  "clw",  "kt",   "bkt", "lekt",
            "bk",      "bkue",    "lc",       "lc.ue",
            "sampling", "is",      "propfair", "comom",  "cub",      "rgf",
            // "divdiff" is the divided-difference closed form of Casale (SIGMETRICS
            // 2017); it needs no think time, since a delay would ask for the integral
            // form of Cor. 3.4, and pfqn_nc refuses one by name. Load-dependent rates
            // ARE served: pfqn_ncld substitutes the limited load-dependent kernel of
            // Casale-Harrison-Ong (Perform. Eval. 2021), Thm. 1, and reports itself as
            // "divdiff.ld/...".
            "divdiff",
            // Chen-O'Cinneide regularization; a Markov chain Monte Carlo estimator of the
            // throughput RATIOS G(N-e_r)/G(N), which supplies no constant of its own
            "mcmc",
            "ger", "rd",   "nrp",
            "nrl",     "nre",     "gm",      "mem",    "comomld",  "ms",
            // Krzesinski state-dependent routing: `solver_nc` intercepts an
            // sdr model whatever the method says, but the gate runs first, so
            // omitting the tokens refused the very name the other three
            // codebases list -- and `sdr.mva` is the only way to reach the
            // Section 4 MVA arm rather than the eq. (16) enumeration.
            "sdr", "sdr.mva",
            // "morrison" is the heavy-usage asymptotic expansion of the generating
            // function for a closed think+DPS network (npfqn_dps_morrison,
            // solver_nc_dps_analyzer). It is the DEFAULT on that shape and
            // inadmissible anywhere else, where the runner refuses it: nothing else
            // in NC can see the DPS weights. Non-product-form, so it returns no lG.
            "morrison",
            // "rec" is the MDD-rec route: the reachable set lives in a decision
            // diagram and the product form supplies the rates. It is the ONLY
            // method admissible on a stochastic Petri net (solver_nc_spn.h) and
            // names the exact loss-network constant in solver_nc_lossn.h. Both
            // routes sit BELOW check_method, so leaving the token out of this
            // list refused the very name the SPN branch's own error tells the
            // caller to use. The other three codebases list it.
            "rec"};
}

/**
 * Port of `SolverNC.isStochasticMethod`.
 *
 * NC is deterministic except for the Monte Carlo integrators, the logistic
 * sampler, the importance-sampling estimators and the Chen-O'Cinneide Markov
 * chain Monte Carlo method, whose answer depends on the seed. The name is TOKENIZED on `.` and `/` so that a runtime-resolved name
 * such as `default/imci` and a prefixed one such as `nc.ls` classify alike.
 */
inline bool is_stochastic_method(const std::string& method) {
    std::string tok;
    std::vector<std::string> toks;
    for (char ch : method) {
        if (ch == '.' || ch == '/') {
            toks.push_back(tok);
            tok.clear();
        } else {
            tok += static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
        }
    }
    toks.push_back(tok);
    for (const std::string& t : toks)
        if (t == "mci" || t == "imci" || t == "ls" || t == "sampling" || t == "is" ||
            t == "mcmc")
            return true;
    return false;
}

/**
 * Port of `SolverNC.resolveMethod`: the feature-driven resolution of
 * `method='default'`.
 *
 * An open network with any non-unit SCV that MEM can carry resolves to `mem`,
 * because the normalizing-constant path would silently exponentialize it.
 *
 * IT IS INERT WITH RESPECT TO THE ANALYZER, and deliberately so: MATLAB's
 * `@@SolverNC/runAnalyzer.m` never consults `resolveMethod`, it reads
 * `options.method` directly. Measured on an M/E2/1 (Source Exp(0.5), FCFS
 * Erlang mean 0.5 scv 0.5), `resolveMethod` returns `mem` while
 * `SolverNC(model).getAvg` reports `default/exact`. The function feeds the
 * feature gate and the AUTO dispatch only, and this port matches that: nothing
 * here calls it on the solve path. See register row N5.
 */
template <class T>
std::string resolve_method(const qn::NetworkStruct<T>& L, const std::string& method) {
    if (method != "default") return method;
    if (!solver_nc_mem_supports(L).supported) return method;
    for (std::size_t i = 0; i < L.nstations; ++i)
        for (std::size_t r = 0; r < L.nclasses; ++r) {
            if (L.disabled[i][r]) continue;
            const double v = num_traits<T>::to_double(L.scv(i, r));
            if (std::isfinite(v) && std::fabs(v - 1.0) > GlobalConstants::FineTol) return "mem";
        }
    return method;
}

/** Port of `runAnalyzerChecks`' method gate: an unlisted method is refused. */
inline void check_method(const std::string& method) {
    const std::vector<std::string> valid = list_valid_methods();
    if (std::find(valid.begin(), valid.end(), method) != valid.end()) return;
    throw UnsupportedError("SolverNC: the '" + method + "' method is unsupported by this solver");
}

namespace detail {

/** Is there a node of this type? */
template <class T>
bool has_node_type(const qn::NetworkStruct<T>& sn, qn::NodeType ty) {
    for (const qn::NodeDef& nd : sn.nodes)
        if (nd.nodetype == ty) return true;
    return false;
}

/**
 * The cache-metric assemblers now live in `line/solvers/cache_metrics.h`, beside
 * the struct they build, because SolverMVA's cache branches need the identical
 * rule. Re-exported into this namespace so the call sites below read unchanged.
 */
using solvers::cache_metrics_of;
using solvers::cache_metrics_of_matrix;

/**
 * Rewrite every finite multiserver station as the rate lattice mu(n)=min(n,c).
 *
 * @return false when the total closed population is not finite, in which case
 *         the reference leaves the model alone
 */
template <class T>
bool multiserver_to_lld(qn::NetworkStruct<T>& sn) {
    double Nt = 0.0;
    for (const qn::JobClass& c : sn.classes) {
        if (std::isinf(c.population)) return false;
        Nt += c.population;
    }
    const std::size_t n = static_cast<std::size_t>(std::llround(Nt));
    if (n < 1) return false;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        std::vector<T> lld(n, num_traits<T>::from_int(1));
        const double c = sn.stations[i].nservers;
        if (std::isfinite(c) && c > 1.0)
            for (std::size_t k = 1; k <= n; ++k)
                lld[k - 1] =
                    num_traits<T>::from_double(std::min<double>(static_cast<double>(k), c));
        sn.stations[i].lldscaling = lld;
    }
    return true;
}

/**
 * Resolves `NcSolverOptions::multiserver` into the handling this solver
 * implements: `default`, `seidmann` or `lld`.
 *
 * Unknown values -- the SolverMVA approximations NC has no counterpart for
 * (`softmin`, `conway`, `krzesinski`, `suri`, `erlang`) -- resolve to `default`
 * rather than throwing, because one options object is commonly reused across
 * solvers. They are NOT honoured silently: `warning`, when given, receives the
 * text the other three codebases pass to `line_warning`. The port has no
 * `line_warning` channel, so a warning travels on `NcSolution::warning`, which
 * `line-cli` prints to stderr; see `NcSolution::warning`.
 */
inline std::string nc_multiserver_policy(const std::string& requested,
                                         std::string* warning = nullptr) {
    if (requested.empty() || requested == "default") return "default";
    if (requested == "seidmann") return "seidmann";
    if (requested == "lld" || requested == "exact" || requested == "loaddep" ||
        requested == "load-dependent")
        return "lld";
    if (warning != nullptr)
        *warning = "SolverNC does not implement config.multiserver='" + requested +
                   "' (it is a SolverMVA approximation); using 'default'. SolverNC accepts "
                   "'default', 'seidmann' and 'lld'.";
    return "default";
}

/**
 * The per-chain population lattice `prod(1+Nchain)`, which prices the exact
 * enumeration `solver_ncld` performs.
 */
template <class T>
double nc_population_lattice(const qn::NetworkStruct<T>& sn) {
    double lattice = 1.0;
    if (!sn.chains.empty()) {
        // sn.chains is (nchains x nclasses) as vector<vector<bool>>, not a Matrix
        for (std::size_t c = 0; c < sn.chains.size(); ++c) {
            double popc = 0.0;
            for (std::size_t r = 0; r < sn.nclasses && r < sn.chains[c].size(); ++r)
                if (sn.chains[c][r]) {
                    const double v = sn.classes[r].population;
                    if (std::isfinite(v)) popc += v;
                }
            lattice *= (1.0 + popc);
        }
    } else {
        for (const qn::JobClass& c : sn.classes)
            if (std::isfinite(c.population)) lattice *= (1.0 + c.population);
    }
    return lattice;
}

}  // namespace detail

/**
 * Is the closed model in NORMAL USAGE, the domain of the Mitra-McKenna PANACEA
 * asymptotic expansion (J. ACM 33(3), 1986)?
 *
 * Normal usage asks that every queueing centre absorb the load the think
 * stations offer it: with rho_j0 = Ztot(j) the aggregate think demand of chain
 * j, r_ij = L_ij / rho_j0 and mu_i(Ntot) the saturation rate,
 *
 *     alpha_i = 1 - (sum_j N_j r_ij) / mu_i(Ntot) > 0     at every centre i.
 *
 * Outside it the {phi(n)} series DIVERGES, which is why `pfqn_panaceald` returns
 * NaN there and `pfqn_ncld` turns that NaN into a refusal rather than a warning.
 * It is a property of the DEMANDS and not of a declared construct, so it has no
 * feature-registry name and cannot live in `nc_feature_set`.
 *
 * The rates are the ones `solver_ncld` would build: one for an ordinary single
 * server, min(n, c) for a finite multiserver (the conversion the runner performs
 * on the 'panald' arm), and the declared lldscaling row when the model sets
 * one. An infinite server is a think station and feeds Ztot.
 */
template <class T>
bool nc_is_normal_usage(const qn::NetworkStruct<T>& sn) {
    if (sn.has_open_classes()) return true;  // no rho_j0 to expand around
    const std::size_t M = sn.nstations, C = sn.nchains;
    const std::vector<double> Nchain = detail::chain_population(sn);
    double NtD = 0.0;
    for (std::size_t c = 0; c < C; ++c)
        if (std::isfinite(Nchain[c])) NtD += Nchain[c];
    const std::size_t Nt = static_cast<std::size_t>(std::llround(NtD));
    if (Nt < 1) return true;  // the empty network: G = 1, nothing to expand

    const mva::ChainDemands<T> d = mva::sn_get_demands_chain(sn);

    std::vector<double> Ztot(C, 0.0);
    for (std::size_t i = 0; i < M; ++i)
        if (std::isinf(sn.stations[i].nservers))
            for (std::size_t c = 0; c < C; ++c)
                Ztot[c] += num_traits<T>::to_double(d.Lchain(i, c));
    for (std::size_t c = 0; c < C; ++c)
        if (Nchain[c] > 0.0 && !(Ztot[c] > 0.0))
            return false;  // no think station on a populated chain's route

    for (std::size_t i = 0; i < M; ++i) {
        const double nserv = sn.stations[i].nservers;
        if (std::isinf(nserv)) continue;
        const std::vector<T>& lld = sn.stations[i].lldscaling;
        double muK;
        if (!lld.empty())
            muK = num_traits<T>::to_double(lld[std::min(Nt, lld.size()) - 1]);
        else if (nserv > 1.0)
            muK = std::min(static_cast<double>(Nt), nserv);
        else
            muK = 1.0;
        if (!(muK > 0.0) || !std::isfinite(muK)) return false;
        double lambda = 0.0;
        for (std::size_t c = 0; c < C; ++c)
            if (Ztot[c] > 0.0)
                lambda += num_traits<T>::to_double(d.Lchain(i, c)) / Ztot[c] * Nchain[c];
        if (!(1.0 - lambda / muK > 0.0)) return false;
    }
    // Every centre cleared the test; a model with no queueing centre at all
    // reaches here too, and there the delay-only constant is exact.
    return true;
}

/**
 * How many queueing (non-infinite-server) stations carry demand from a CLOSED
 * chain?
 *
 * That is the row count L reaches `pfqn_nc` and `pfqn_comomrm_ld` with, once the
 * delay rows have been folded into Z and the zero-demand rows dropped. Zero when
 * the model has no closed population at all.
 */
template <class T>
std::size_t nc_closed_queueing_stations(const qn::NetworkStruct<T>& sn) {
    const std::size_t M = sn.nstations, C = sn.nchains;
    const std::vector<double> Nchain = detail::chain_population(sn);
    bool any_closed = false;
    std::vector<bool> closed(C, false);
    for (std::size_t c = 0; c < C; ++c) {
        closed[c] = std::isfinite(Nchain[c]) && Nchain[c] > 0.0;
        if (closed[c]) any_closed = true;
    }
    if (!any_closed) return 0;

    const mva::ChainDemands<T> d = mva::sn_get_demands_chain(sn);
    std::size_t nq = 0;
    for (std::size_t i = 0; i < M; ++i) {
        if (std::isinf(sn.stations[i].nservers)) continue;
        for (std::size_t c = 0; c < C; ++c)
            if (closed[c] &&
                std::fabs(num_traits<T>::to_double(d.Lchain(i, c))) > GlobalConstants::FineTol) {
                ++nq;
                break;
            }
    }
    return nq;
}

/**
 * May `method` run on this model? "" when it may, otherwise the reason it may
 * not, in the words the runner refuses with.
 *
 * ONE PREDICATE, TWO CALLERS. `solver_nc_solve` asks it once, ahead of the
 * dispatch, and throws on a non-empty answer; `auto_family_refusal` asks it so
 * that `auto_find_solver` never offers a (family, method) pair that would throw,
 * and so that the ranking never delegates to one. Two copies of these rules is
 * precisely how the report and the run drift apart, which is the failure this
 * function exists to prevent, so a new rule goes here and not at a call site.
 *
 * ONLY WHAT THE FEATURE REGISTRY CANNOT NAME LIVES HERE. A feature set declares
 * what the method ACCEPTS, so it can refuse a model for HAVING a construct and
 * never for lacking one: "closed population only" and "no think time" are said
 * in `qn::nc_feature_set` by dropping OpenClass and SchedStrategy_INF, while
 * "requires a cache", "requires state-dependent routing", "requires a loss
 * network", "requires exactly two stations" and "requires normal usage" have no
 * such form and are decided here.
 *
 * `for_report` says WHICH QUESTION IS BEING ASKED, and for two method names the
 * two questions have different answers:
 *
 *   true  -- "should `auto_find_solver` offer this pair?" A pair that comes back
 *            as a table of zeros must not be offered, so the answer is no.
 *   false -- "what does the reference DO when asked for it by name?" For 'mmint2'
 *            and 'gleint' outside their shape the reference deliberately WARNS AND
 *            RETURNS A ZERO TABLE (pfqn_nc.m, case {'mmint2','gleint'}: lG = [] and
 *            return, unconditionally), and a caller who names the method keeps that
 *            answer -- which is what `test_nc.cpp` pins.
 *
 * THE ASYMMETRY IS A RULING, NOT AN OVERSIGHT (2026-07-25, reaffirmed when this
 * gate was added): the report answers "should this be offered" and the run answers
 * "what does the reference do". 'comomld' is NOT in that bucket --
 * `pfqn_comomrm_ld` refuses "The solver accepts at most a single queueing station."
 * natively -- so it is refused on both paths.
 *
 * @param sn         the refreshed struct
 * @param method     the concrete method name
 * @param slotted    true on the discrete-time route, which answers for itself
 * @param for_report true when the caller is the report, false when it is the run
 */
template <class T>
std::string nc_method_refusal(const qn::NetworkStruct<T>& sn, const std::string& method,
                              bool slotted = false, bool for_report = true) {
    const std::string m = method.empty() ? std::string("default") : method;

    // The discrete-time route answers for itself: `solver_nc_dt` decides
    // admissibility on the slot lattice, and every gate below is written about a
    // continuous-time queueing network.
    if (slotted) return "";

    // -- discriminatory processor sharing ---------------------------------
    // Morrison's heavy-usage expansion is the ONLY NC route that can see the DPS
    // weights; every other method builds a product-form normalizing constant
    // that silently drops them and answers with the egalitarian-PS network,
    // which is a wrong number rather than a coarse one.
    if (nc_is_dps_model(sn)) {
        if (m != "default" && m != "morrison")
            return "SolverNC: method '" + m +
                   "' cannot represent the DPS weights of a discriminatory processor-sharing "
                   "station; it would return the egalitarian-PS network. Use method 'default' or "
                   "'morrison' (npfqn_dps_morrison), SolverMVA, SolverFLD or SolverCTMC.";
        return "";
    }
    if (sn_has_dps(sn))
        // A DPS station outside Morrison's shape. SchedStrategy_DPS is declared
        // in the feature set because a boolean feature cannot express "this shape
        // only"; this is that imperative half.
        return "SolverNC analyzes a discriminatory processor-sharing station only in the shape "
               "Morrison's expansion is derived for: a CLOSED network of exactly two stations, one "
               "infinite-server (think) station and one single-server DPS station, exponential "
               "service, each class visiting the two equally often. Use SolverMVA, SolverFLD or "
               "SolverCTMC for any other DPS model.";
    if (m == "morrison")
        // The method named on a model that is not the shape at all -- not even a
        // DPS station in it. Left ungated it reaches no route of its own and
        // falls through to the ordinary normalizing-constant path, which would
        // answer the product-form model UNDER THE CALLER'S LABEL.
        return "SolverNC: method 'morrison' is the heavy-usage expansion of a CLOSED network of "
               "exactly two stations, one infinite-server (think) station and one single-server DPS "
               "station with exponential service, which this model is not. Remove the method option "
               "to let SolverNC choose, or use SolverMVA, SolverFLD or SolverCTMC.";

    // -- Krzesinski state-dependent routing -------------------------------
    // An SDR model is intercepted by `solver_nc_sdr` whatever the method says,
    // so reaching the second test means the model declares none.
    if (!sn.sdr.branch.empty()) return "";
    if (m == "sdr" || m == "sdr.mva")
        return "SolverNC: method '" + m +
               "' requires state-dependent routing, which this model does not declare.";

    // -- stochastic Petri net ---------------------------------------------
    // A net is served only by the MDD-rec route, and none of the gates below --
    // written about stations, capacities and the queueing-network product form --
    // says anything about a net. `spn_pf` decides its product-form class, by name.
    for (std::size_t ind = 0; ind < sn.nodes.size(); ++ind)
        if (sn.nodes[ind].nodetype == lang::NodeType::Place) {
            if (m != "default" && m != "rec")
                return "solver_nc: a stochastic Petri net is solved by the MDD-rec route; method '" +
                       m + "' is a normalizing-constant algorithm for queueing networks. Use 'rec' "
                       "or 'default'";
            return "";
        }

    // -- order-independent stations ---------------------------------------
    // Every method other than the four listed reads the single-job rate mu([r])
    // of an OI station: the rank rate mu(n) is silently dropped and the answer is
    // that of an ordinary queue.
    if (nc_is_oi_model(sn)) {
        if (m != "default" && m != "exact" && m != "is" && m != "sampling")
            return "SolverNC: method '" + m +
                   "' cannot represent the rank rate mu(n) of an order-independent station; use "
                   "method 'default' or 'exact' (pfqn_ncoi), 'is', SolverMVA, or SolverCTMC.";
        return "";
    }

    // -- caches -------------------------------------------------------------
    // 'rayint' and 'spm' both name the SPM saddle point of a cache (and, on a
    // retrieval model, the ray/WKB delayed-hit expansion), so they are admissible
    // here and nowhere else.
    for (std::size_t ind = 0; ind < sn.nodes.size(); ++ind)
        if (sn.nodes[ind].nodetype == qn::NodeType::Cache) {
            if (m == "exact" && nc_is_noreentrant_cache(sn)) {
                const auto itp = sn.nodeparam.find(ind + 1);
                // cache_prob_erec is exact for the exchangeable (RR/FIFO) family
                // only; a recency-based policy would silently receive the
                // exchangeable answer, so the exact route refuses it.
                if (itp != sn.nodeparam.end() &&
                    itp->second.replacestrat != lang::ReplacementStrategy::RR &&
                    itp->second.replacestrat != lang::ReplacementStrategy::FIFO)
                    return "solver_nc_cache_analyzer: NC does not support the exact solution of "
                           "this cache replacement policy -- only RR and FIFO are exchangeable, "
                           "and a recency-based policy (LRU, h-LRU, q-LRU, CLIMB) would silently "
                           "receive the exchangeable answer. Use the default (approximate) method "
                           "or SolverCTMC";
            }
            return "";
        }
    if (m == "rayint" || m == "spm")
        return "SolverNC: method " + m +
               " names the SPM saddle point of a cache and, on a retrieval model, the ray/WKB "
               "delayed-hit expansion; this model declares no Cache node.";

    // -- loss networks and finite capacity regions --------------------------
    if (nc_is_lossn_model(sn)) return "";  // 'ms', 'erlangfp', 'rec' all run here
    if (nc_has_lossn_shape(sn))
        return "SolverNC: the Finite Capacity Region holds a single infinite server but does not "
               "apply DROP to every class; holding an arrival back (WAITQ) or blocking the server "
               "(BAS/BBS/RSRD) keeps the job in the region while it waits, which the Erlang loss "
               "model has no state for -- use DROP, or SolverCTMC/SolverJMT";
    if (m == "erlangfp")
        return "SolverNC: the 'erlangfp' Erlang fixed point applies only to a loss network, which "
               "is an open model whose single Delay sits inside a Finite Capacity Region under a "
               "DROP rule; this model declares no such region (see nc_is_lossn_model)";
    if (m == "ms")
        return "SolverNC: method 'ms' is admissible only on a loss network (open model, one DROP "
               "region holding a single Delay).";
    if (m == "rec")
        return "SolverNC: method rec is the MDD-rec route, admissible on a stochastic Petri net or "
               "on a loss network (open model, one DROP region holding a single Delay); this model "
               "is neither.";
    if (!sn.regions.empty())
        return "SolverNC: this model applies a Finite Capacity Region to queueing stations, whose "
               "aggregate population limit no normalizing-constant algorithm here enforces; only "
               "the loss network -- one region over a single infinite server, DROP on every "
               "class -- is solvable. Use SolverCTMC or SolverJMT, or setCapacity for a "
               "single-station limit";

    // -- PANACEA's domain ----------------------------------------------------
    // Normal usage is a property of the demands rather than of a declared
    // construct, so it has no feature name; an open chain is refused earlier by
    // the closed-population feature set of the load-dependent evaluators.
    //
    // BOTH TOKENS ARE GATED, because `pfqn_ncld` evaluates 'pana' and
    // 'panald' with the SAME `pfqn_panaceald` -- its case label covers both --
    // so on a model carrying a rate lattice the load-INDEPENDENT name reaches the
    // load-dependent expansion and throws with it. Off that lattice 'pana'
    // takes its own `pfqn_nc` arm, which warns and returns an empty constant
    // rather than throwing, so it is left alone there. CLASS-dependent scaling
    // diverts the whole model to `solver_nc_conv`, which never reads the method
    // at all -- and only class-dependent, because that is the one `solver_ncld`
    // diverts on in this port; a joint-dependent model still reaches the kernel.
    if ((m == "pana" || m == "panald") && !sn.has_open_classes()) {
        bool diverted_to_conv = false, any_lld = false;
        for (std::size_t i = 0; i < sn.nstations; ++i) {
            if (static_cast<bool>(sn.stations[i].cdscaling)) diverted_to_conv = true;
            if (!sn.stations[i].lldscaling.empty()) any_lld = true;
        }
        const bool reaches_ld_kernel = (m == "panald") || any_lld;
        if (!diverted_to_conv && reaches_ld_kernel && !nc_is_normal_usage(sn)) {
            const std::string why =
                "the model is not in normal usage, so the 'panald' asymptotic expansion does "
                "not apply. Use 'exact', 'clw' or an approximate load-dependent method instead.";
            if (m == "pana")
                return "SolverNC: method 'pana' reaches the load-dependent kernel on this "
                       "model, where pfqn_ncld evaluates it as 'panald', and " + why;
            return "SolverNC: " + why;
        }
    }

    // -- the single-queueing-station recursions --------------------------------
    // Two families are stated for a model with a delay and ONE queueing station,
    // and neither can say so with a feature name: it is a COUNT, and a feature
    // set has no arithmetic. `pfqn_nc` states it for 'mmint2'/'gleint' in those
    // words and `pfqn_comomrm_ld` refuses with "The solver accepts at most a
    // single queueing station."
    //
    // The count is taken over the CLOSED chains only, and the rule is inactive
    // without a closed population, because `pfqn_nc` answers an open network with
    // the exact open formulas BEFORE its method switch -- the method name is never read
    // there, so a purely open model with three queues runs these names correctly
    // today and must go on doing so.
    // 'mmint2' and 'gleint' are gated for the REPORT ONLY: `pfqn_nc` answers them
    // with an empty constant and the caller renders a table of zeros, which is a
    // pair the report must not offer and a run the reference nonetheless performs.
    // See `for_report` above.
    if (m == "comomld" || (for_report && (m == "mmint2" || m == "gleint"))) {
        const std::size_t nq = nc_closed_queueing_stations(sn);
        if (nq > 1) {
            if (m == "comomld")
                return "SolverNC: method 'comomld' is the load-dependent CoMoM recursion, and "
                       "pfqn_comomrm_ld accepts at most a single queueing station; this model has " +
                       std::to_string(nq) + ".";
            return "SolverNC: the '" + m +
                   "' method requires a model with a delay and a single queueing station; this "
                   "model has " + std::to_string(nq) + ".";
        }
    }

    // -- 'exact' outside its domain ------------------------------------------
    if (m == "exact" && sn.has_open_classes()) {
        for (std::size_t i = 0; i < sn.nstations; ++i)
            if (std::isfinite(sn.stations[i].nservers) && sn.stations[i].nservers > 1.0)
                return "solver_nc_analyzer: the NC solver cannot provide exact solutions for open "
                       "or mixed multiserver queueing networks. Remove the 'exact' option";
    }
    return "";
}


/**
 * The gates, the multiserver conversion and the dispatch of
 * `@@SolverNC/runAnalyzer.m`, without the metric filter.
 *
 * Kept separate so that the reported log normalizing constant comes from the
 * SAME model the metrics do: reaching the dispatch without the conversion above
 * would evaluate the constant of the Seidmann-approximated network while the
 * table reported the load-dependent one.
 *
 * A model with a Fork is solved through the shared fork-join fixed point
 * (`fj_driver.h`), which drives `nc_dispatch` as its inner solve on the
 * transformed model; a model without one runs the dispatch exactly once.
 */
template <class T>
NcSolution<T> solver_nc_solve(const qn::NetworkStruct<T>& L_in, const NcSolverOptions& opt_in) {
    check_method(opt_in.method);
    // runAnalyzerChecks' universal feature gate, AFTER the method gate. Gated on
    // `L_in`, before multiserver_to_lld rewrites lldscaling. NC DECLARES Region,
    // so a loss network passes straight through to the imperative split below.
    qn::feature_gate("SolverNC", qn::nc_feature_set(opt_in.method), L_in);

    // THE STRUCTURAL METHOD GATE, asked once and in one place.
    //
    // `nc_method_refusal` holds every rule of the form "this method has no route
    // on this model": the DPS shape, state-dependent routing, the Petri net, the
    // order-independent rank rate, the cache and loss-network tokens, PANACEA's
    // normal usage. `auto_family_refusal` asks the SAME function, which is what
    // keeps `auto_find_solver` from offering a pair that would throw here.
    {
        // for_report=false: this is the RUN, and it asks what the reference DOES
        // rather than what the report should offer. The two answers differ for
        // 'mmint2'/'gleint', which `pfqn_nc` answers with an empty constant and a
        // zero table; see `nc_method_refusal`.
        const std::string refusal =
            nc_method_refusal(L_in, opt_in.method, opt_in.slotted, /*for_report=*/false);
        if (!refusal.empty()) throw UnsupportedError(refusal);
    }

    NcSolverOptions opt = opt_in;
    qn::NetworkStruct<T> L = L_in;
    NcSolution<T> res;

    // Discrete time before every other gate: the slot lattice is a property of
    // the MODEL, not of a method, and the imperative gates below (finite
    // capacity, the multiserver-to-lldscaling rewrite, the load-dependence
    // routing) all assume a continuous time scale. `solver_nc_dt` refuses a
    // model outside the discrete-time product form rather than falling through.
    if (opt.slotted) return solver_nc_dt(L, opt);

    // THE STRUCTURAL FINITE-CAPACITY GATE, port of `@@SolverNC/SolverNC.m:180-188`.
    //
    // It was missing here, and its absence was not a missing message: a closed
    // two-queue model with setCapacity(1) on the second station SOLVED, and
    // reported the UNCONSTRAINED product-form answer (QLen 0.852459/1.147541,
    // the values of the same model with no buffer at all) where MATLAB, the JAR
    // and python all refuse by name. A wrong number, silently, is the failure
    // mode `check_binding_capacity` exists to stop -- SolverMVA has called it
    // since the port (`mva_check_finite_capacity`) and NC never did.
    //
    // AFTER THE SLOTTED RETURN, as in the reference: a finite buffer on a
    // Bernoulli server is the loss system of Daduna's corollary 2.8, which
    // `solver_nc_dt` solves exactly, so this gate must not see it.
    //
    // THE ONE EXEMPTION IS mem.blocking, which represents the buffer as a
    // censored GE/GE/c/0;N queue and therefore DOES honour it. The reference
    // exempts one more, the single-station M/M/1/K with tail drop that its
    // `qsys_mm1k_loss` branch answers exactly; this port has no such branch in
    // `nc_dispatch`, so that shape is refused here rather than answered
    // unconstrained. Add the exemption in the same change that ports the branch.
    //
    // IT READS `opt.method`, NOT `resolve_method`, and the difference is not
    // cosmetic. `nc_dispatch` reaches the MEM algorithm on the literal name --
    // `if (opt.method == "mem")` -- and nothing on the solve path resolves
    // `default` into it, so exempting a `default` run because `resolve_method`
    // would have called it `mem` skips the gate and then dispatches somewhere
    // that does NOT honour the buffer. That is the very failure this gate
    // exists to stop. (The reference gates on its resolved method and its
    // analyzer reads `options.method` directly, so MATLAB has the same seam;
    // do not copy it.)
    {
        bool mem_blocking = false;
        if (opt.method == "mem") {
            const MemSupport ms = solver_nc_mem_supports(L);
            mem_blocking = ms.supported && ms.blocking;
        }
        if (!mem_blocking) qn::check_binding_capacity("SolverNC", L);
    }

    if (L.has_fork()) {
        // The transform turns the fork into a router, the join into a delay and
        // the branches into auxiliary open classes. `nc_dispatch` then sees a
        // plain mixed network, which is why none of the specialised NC routes
        // has to know about forks.
        mva::FjMmt<T> tr = mva::fj_fork_join_transform(L, opt.fork_join);
        std::vector<T> lam(tr.V.classes.size() + 1,
                           num_traits<T>::from_double(GlobalConstants::FineTol));
        mva::MvaOptions mopt;
        mopt.method = opt.method;
        mopt.tol = opt.tol;
        mopt.iter_tol = opt.iter_tol;
        mopt.iter_max = opt.iter_max;
        mopt.fork_join = opt.fork_join;
        mopt.base_has_fork = true;
        NcSolverOptions inner = opt;
        inner.base_has_fork = true;
        std::string am;
        res.sol = mva::fj_fixed_point(L, tr, lam, mopt, [&inner, &am](qn::NetworkStruct<T>& V) {
            const NcSolution<T> d = nc_dispatch(V, inner);
            am = d.actualmethod;
            return d.sol;
        });
        res.actualmethod = am;
        return res;
    } else {
        // Closed think+DPS network -> Morrison's heavy-usage generating-function
        // expansion, the DEFAULT for that shape. Intercepted FIRST, ahead of every
        // other branch: "SchedStrategy_DPS" is now inside the solver's reach, and
        // each of the branches below would silently drop the weights and answer
        // with the egalitarian-PS network. Not a product-form route: lG is NaN.
        if (nc_is_dps_model(L)) return solver_nc_dps_analyzer(L, opt);
        // The four refusal arms that used to follow -- another method on a DPS
        // model, a DPS station outside Morrison's shape, "morrison" on a model with
        // no DPS station at all, and "sdr"/"sdr.mva" on a model declaring no
        // state-dependent routing -- moved into `nc_method_refusal` above, with
        // their wording unchanged.

        // Order-independent networks are intercepted BEFORE the
        // multiserver-to-lldscaling rewrite below: an OI station is multiserver
        // but is not a plain min(n,c) load-dependent station, so the rewrite
        // would send it down a path that cannot represent its rate function.
        if (nc_is_oi_model(L) && (opt.method == "default" || opt.method == "exact"))
            return solver_nc_oi_analyzer(L, opt);

        // A loss network is intercepted above the gates below: it is judged by
        // has_product_form on the 'exact' path and its infinite server would be
        // run through the multiserver-to-lldscaling rewrite, neither of which
        // applies to a model whose only dynamics are the region's admission rule.
        if (nc_is_lossn_model(L)) return solver_nc_lossn_analyzer(L, opt).sol;

        // A stochastic Petri net takes the MDD-rec route: the reachable set
        // lives in a decision diagram and the product form supplies the rates,
        // so none of the queueing-network gates below say anything about it.
        // `spn_pf` is where a net's product-form class is decided, by name.
        for (std::size_t ind = 0; ind < L.nodes.size(); ++ind)
            if (L.nodes[ind].nodetype == lang::NodeType::Place)
                return solver_nc_spn_analyzer(L, opt).sol;

        bool multiserver = false;
        for (const qn::Station<T>& st : L.stations)
            if (std::isfinite(st.nservers) && st.nservers > 1.0) multiserver = true;
        bool anyLld = false;
        for (const qn::Station<T>& st : L.stations)
            if (!st.lldscaling.empty()) anyLld = true;

        // How this model's finite multiserver stations are represented. The
        // shipped "default" reproduces the historical dispatch exactly, so no
        // result moves unless config.multiserver is set.
        std::string ms_warning;
        const std::string ms_policy =
            detail::nc_multiserver_policy(opt.multiserver, &ms_warning);

        if (opt.method == "default") {
            // The two-station Delay-plus-multiserver shape: every SolverLN layer
            // submodel. Exact through load-dependent CoMoM when product-form,
            // Seidmann through 'comom' when not.
            if (L.nstations == 2 && !detail::has_node_type(L, qn::NodeType::Cache) &&
                detail::has_node_type(L, qn::NodeType::Delay) && multiserver) {
                if (L.has_product_form() && !anyLld) {
                    if (detail::multiserver_to_lld(L)) anyLld = true;
                } else {
                    opt.method = "comom";
                }
            } else if (ms_policy == "lld" && multiserver && !anyLld && L.has_product_form() &&
                       detail::nc_population_lattice(L) <= 6000.0) {
                // config.multiserver="lld" generalizes the exact load-dependent
                // lattice of the branch above to any closed product-form model,
                // under the same 6000-state enumeration budget. Off unless asked
                // for: with the shipped "default" policy this branch never runs
                // and the model keeps Seidmann's approximation, as it always has.
                if (detail::multiserver_to_lld(L)) anyLld = true;
            }
        } else if (opt.method == "exact" || opt.method == "is" || opt.method == "panald") {
            // 'is' and 'panald' need the same model as 'exact', i.e. the same
            // multiserver conversion.
            //
            // 'is' on an OI / P&S model is the specialised sampler, not the
            // generic one: an OI station is multiserver but not min(n,c), and a
            // P&S tandem has only per-communicating-class product form, so
            // has_product_form is false and the gate below would reject it.
            // nc_is_oi_model too, not only nc_is_pas_model: the latter demands
            // that BOTH stations be OI/PAS, so a Delay + OI cycle -- the canonical
            // topology pfqn_oi_is exists to sample -- fell to the guard below.
            if (opt.method == "is" && (nc_is_pas_model(L) || nc_is_oi_model(L))) {
                // handled by solver_nc_analyzer (pfqn_pas_is / pfqn_oi_is)
            } else if (nc_is_lossn_model(L)) {
                // A LOSS NETWORK (open, one DROP region holding a single Delay)
                // IS product form -- the truncated Poisson law the residue
                // transform of solver_nc_lossn evaluates exactly under 'exact'
                // -- but any region reads as blocking, so has_product_form says
                // no. Exempted here as `nc_method_refusal` already exempts it.
            } else if (!L.has_product_form())
                throw UnsupportedError(
                    "SolverNC: the '" + opt.method +
                    "' method requires the model to have a product-form solution, and this model "
                    "does not");
            // Only a GENUINE multiserver is converted: filling lldscaling with
            // ones on an all-single-server model is a semantic no-op that forces
            // the load-dependent path, whose recovery mishandles a
            // single-station-confined closed chain.
            // config.multiserver="seidmann" asks for Seidmann's approximation on
            // this arm too, so the conversion is skipped. Off by default.
            if (!anyLld && multiserver && ms_policy != "seidmann") {
                if (detail::multiserver_to_lld(L)) anyLld = true;
            }
        }

        // Caches, in the reference's order: a delayed-hit retrieval system
        // first (open or closed), then the non-reentrant Source-Cache-Sink
        // model, then any other cache as an integrated caching-queueing network.
        // The order is the contract -- a retrieval cache is ALSO a cache, and a
        // Source-Cache-Sink model is also "a model with a cache node".
        if (nc_has_retrieval(L)) {
            if (detail::has_node_type(L, qn::NodeType::Source)) {
                NcRetrievalSolution<T> rr = solver_nc_retrieval_analyzer(L, opt);
                rr.sol.cache = detail::cache_metrics_of(L, rr.hitprob, rr.missprob, rr.delayedprob,
                                                        rr.latency, rr.hitproblist, rr.itemprob,
                                                        std::vector<T>());
                return rr.sol;
            }
            NcCacheqnRetrievalSolution<T> cr = solver_nc_cacheqn_retrieval_analyzer(L, opt);
            cr.sol.cache = detail::cache_metrics_of(L, cr.hitprob, cr.missprob, cr.delayedprob,
                                                    cr.latency, cr.hitproblist, Matrix<T>(),
                                                    std::vector<T>());
            return cr.sol;
        }
        if (nc_is_noreentrant_cache(L)) {
            NcCacheSolution<T> cr = solver_nc_cache_analyzer(L, opt);
            // The non-reentrant model reports the per-list hit fractions and the
            // per-item law but no scalar hit probability of its own: the hit and
            // miss columns are the row sums of `hitproblist` and its complement,
            // and `cache_metrics_of` derives them rather than leaving the total
            // row of getAvgCacheTable empty on the one branch that knows most.
            cr.sol.cache = detail::cache_metrics_of(L, cr.hitprob, cr.missprob,
                                                    std::vector<T>(), std::vector<T>(),
                                                    cr.hitproblist, cr.itemprob, cr.listcost);
            return cr.sol;
        }
        if (nc_is_cacheqn(L)) {
            const NcCacheqnSolution<T> cr = solver_nc_cacheqn_analyzer(L, opt);
            res = cr.sol;
            res.refreshed_struct.reset(new qn::NetworkStruct<T>(cr.refreshed));
            // hitprob/missprob are (ncaches x nclasses) here, one row per cache,
            // against the (K) vectors the retrieval branches return.
            res.cache = detail::cache_metrics_of_matrix(L, cr.hitprob, cr.missprob);
            // the per-item law rides in separately: this branch has one per cache
            // rather than one per model
            for (std::size_t ci = 0; ci < res.cache.caches.size() && ci < cr.itemprob.size(); ++ci)
                res.cache.caches[ci].itemprob = cr.itemprob[ci];
            return res;
        }

        // A region this solver cannot enforce must STOP the solve, not be
        // dropped from it. Everything below computes an unconstrained answer,
        // so a model that reaches it with a live region would be reported as
        // though the region were absent -- numbers that are not wrong for any
        // model the user described. The reference raises at the same two points
        // (`runAnalyzer.m:326,334-336`) for the same reason.
        //
        // The loss-network SHAPE with a rule other than all-DROP is named
        // separately from a region on ordinary queueing stations, because the
        // two have different remedies: the first is asking for blocking, where
        // switching the rule to DROP makes the model solvable here, while the
        // second needs a solver that carries the region as state.
        // The token gates that used to stand here -- the loss-network shape under
        // a rule other than all-DROP, a Finite Capacity Region on queueing
        // stations, "erlangfp"/"ms"/"rec" off a loss network and "rayint"/"spm"
        // with no Cache node -- moved into `nc_method_refusal`, which decides them
        // from the same struct before the dispatch begins and which
        // `auto_family_refusal` asks too; their wording is unchanged. The six
        // load-dependent evaluators on an OPEN chain are now refused by the
        // feature set instead (`nc_feature_set` drops OpenClass from them):
        // "closed population only" is a rule the registry CAN name, and naming it
        // there is what makes `auto_find_solver` drop the row rather than report
        // it runnable.

        if (anyLld || nc::detail::has_scaling(L)) {
            res = solver_ncld_analyzer(L, opt);
        } else if (opt.method == "rd" || opt.method == "nrp" || opt.method == "nrl" ||
                   opt.method == "nre" || opt.method == "comomld" || opt.method == "panald") {
            res = solver_ncld_analyzer(L, opt);
        } else {
            res = solver_nc_analyzer(L, opt);
        }
        // A config.multiserver value this solver does not implement was silently
        // ignored before; it now says so. Joined rather than overwriting, as the
        // analyzer warnings are.
        if (!ms_warning.empty())
            res.warning = res.warning.empty() ? ms_warning : res.warning + " " + ms_warning;
        return res;
    }
}

/**
 * Port of `@@SolverNC/runAnalyzer.m` for the `lang='matlab'` path: solve, then
 * apply the metric filter `@@NetworkSolver/getAvg` puts between the analyzer and
 * the caller.
 */
template <class T>
mva::AvgResult<T> solver_nc_run_analyzer(const qn::NetworkStruct<T>& L_in, const NcSolverOptions& opt_in) {
    const std::string origmethod = opt_in.method;
    const NcSolution<T> d = solver_nc_solve(L_in, opt_in);
    const mva::MvaSolution<T>& s = d.sol;
    const std::string& actualmethod = d.actualmethod;
    const qn::NetworkStruct<T>& L = L_in;

    const std::size_t M = L.nstations, K = L.nclasses;
    std::vector<std::vector<bool>> mask(M, std::vector<bool>(K, false));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k)
            mask[i][k] = num_traits<T>::to_double(s.R(i, k)) < 10.0 * GlobalConstants::FineTol;
    std::vector<std::vector<bool>> srcmask(M, std::vector<bool>(K, false));
    for (std::size_t i = 0; i < M; ++i)
        if (L.stations[i].nodetype == qn::NodeType::Source)
            for (std::size_t k = 0; k < K; ++k) srcmask[i][k] = true;

    mva::AvgResult<T> out;
    out.QN = mva::filter_metric(L, s.Q, mva::MetricKind::QLen, &mask);
    out.UN = mva::filter_metric(L, s.U, mva::MetricKind::Util, &mask);
    out.RN = mva::filter_metric(L, s.R, mva::MetricKind::RespT, nullptr);
    out.TN = mva::filter_metric(L, s.Tp, mva::MetricKind::Tput, nullptr);
    // Both ArvR and ResidT come from the refreshed struct when the cacheqn
    // analyzer supplied one: its cache self-switch is normalized at the ACTUAL
    // hit/miss split, so the visits back the carried rate without the over-route
    // inflation. filter_metric still keys on the base L for its topology masking.
    const qn::NetworkStruct<T>& refL = d.refreshed_struct ? *d.refreshed_struct : L;
    out.WN = mva::filter_metric(L, mva::sn_get_residt_from_respt(refL, out.RN),
                                mva::MetricKind::ResidT, nullptr);
    out.AN = mva::filter_metric(L, mva::sn_get_arvr_from_tput(refL, out.TN), mva::MetricKind::ArvR,
                                &srcmask);
    out.CN = s.C;
    out.XN = s.X;
    out.method = origmethod;
    // `runAnalyzer` reports 'default/<algorithm>' when the caller asked for the
    // default and something more specific ran, so the algorithm that produced
    // the numbers is never lost.
    out.actualmethod = (origmethod == "default" && !actualmethod.empty() &&
                        actualmethod != "default")
                           ? "default/" + actualmethod
                           : actualmethod;
    // Approximated, not refused, for the reason SolverMVA gives: a normalizing
    // constant counts a fed-back visit as an ordinary re-entry, having no way
    // to express a job that keeps its server. `@@SolverNC/runAnalyzer.m` warns
    // and solves, and this text is that warning verbatim.
    if (L.has_immediate_feedback())
        out.warning =
            "SolverNC does not handle immediate feedback (immfeed); the solver will treat "
            "self-loops as class-switching with re-queueing.";
    // A warning raised by an analyzer (the cache branch's cost-cap conditions)
    // must survive the runner; both can fire, so they are joined rather than
    // one overwriting the other.
    if (!d.warning.empty())
        out.warning = out.warning.empty() ? d.warning : out.warning + " " + d.warning;
    out.listcost = d.listcost;
    // What the cache branch measured, carried to `getAvgCacheTable` and
    // `getAvgItemTable`; empty on every model without a Cache node.
    out.cache = d.cache;
    // `runAnalyzer.m`'s last line: the log normalizing constant is stored beside
    // the average table rather than recomputed, so `getProbNormConstAggr` costs
    // the caller nothing once `getAvg` has run.
    out.lognormconst = s.lG;
    out.iter = s.iter;
    {
        const double lg = num_traits<T>::to_double(s.lG);
        if (std::isfinite(lg))
            line::util::LineConsole::step("normalizing constant obtained: log G = %.6g", lg);
    }
    return out;
}

/** Port of `@@SolverNC/getProbNormConstAggr.m`: the log normalizing constant. */
template <class T>
double solver_nc_lognormconst(const qn::NetworkStruct<T>& L, const NcSolverOptions& opt) {
    return solver_nc_solve(L, opt).sol.lG;
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_RUNNER_H
