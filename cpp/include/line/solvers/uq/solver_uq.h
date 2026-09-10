/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_UQ_SOLVER_UQ_H
#define LINE_SOLVERS_UQ_SOLVER_UQ_H

/**
 * SolverUQ: uncertainty quantification by expansion over a Prior.
 *
 * Port of `matlab/src/solvers/UQ/@@UQ/UQ.m`. The model carries one or more
 * `Prior` distributions (lang/prior.h); UQ reduces them to a set of weighted
 * DESIGN POINTS, each of which is a concrete model with every Prior replaced by
 * one alternative, solves each with an ordinary solver, and reports the
 * prior-weighted expectation of every metric together with the per-point
 * results the expectation was formed from.
 *
 * WHAT THE WEIGHTS MEAN. E[Q] = sum_l w_l Q(theta_l) is the unconditional
 * expectation of Trivedi and Bobbio (2017), Eq. (3.68): an average over MODELS,
 * not over jobs. Its spread -- `uq_moments`, `uq_credible_interval` -- is the
 * epistemic uncertainty in the answer, and it is the reason the per-point table
 * is kept rather than reduced away: a mean of 4.2 over points at 1.1 and 12.4
 * is a different statement from a mean of 4.2 over points at 4.1 and 4.3, and
 * only the design carries the difference.
 *
 * THE DESIGN, and why it is a tensor product. Each Prior is discretized on its
 * own and the design is the product of the per-Prior alternative sets, so the
 * joint weight is the product of the marginal weights. That is the
 * product-density case of f(theta_1, ..., theta_l) in Eq. (3.67) and it ASSUMES
 * THE PRIORS ARE INDEPENDENT; a joint prior over several parameters is not
 * expressible here, in this port or in the reference. The size is capped at
 * `kMaxDesignPoints` because a design point is a full solver run, and beyond
 * the cap the Monte Carlo design -- whose cost does not grow with the number of
 * Priors -- is the right tool. The cap REFUSES rather than truncating: a design
 * silently cut to 4096 of 20000 points would report an expectation against a
 * prior nobody wrote.
 *
 * WHAT SOLVES A DESIGN POINT is supplied by the caller as a `UqStageSolver`,
 * the C++ spelling of the reference's `solverFactory` argument (`UQ(model,
 * @@SolverMVA)`). There is no default: the inner solver decides both the
 * accuracy and the admissible feature set of every number reported here, and
 * choosing one silently would answer a question the caller did not ask.
 * `uq_dispatch.h` builds one from a solver name.
 *
 * THE FEATURE GATE IS THE INNER SOLVER'S. UQ itself declares only `Prior`
 * (`uq_feature_set`) and applies no feature gate of its own; what `UQ.supports`
 * decides is only whether the model carries a Prior at all. The real check
 * happens per design point, inside the stage solver, on a model from which the
 * Prior has already been removed; that is what makes "SolverMVA cannot solve
 * this model" reach the caller as SolverMVA's own refusal instead of a UQ
 * paraphrase of it.
 */

#include <algorithm>
#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "line/api/pfqn/pfqn_mva_interval.h"
#include "line/lang/lang_types.h"
#include "line/lang/prior.h"
#include "line/lang/qn/feature_set.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/util/error.h"

namespace line {
namespace uq {

/**
 * The cap on the tensor-product design, MATLAB `UQ.MaxDesignPoints`.
 *
 * A design point is one full solver run, so this bounds the cost of a
 * quadrature design over several Priors.
 */
inline constexpr std::size_t kMaxDesignPoints = 4096;

/** `UQ.defaultOptions` plus the stream the Monte Carlo design draws from. */
struct UqOptions {
    /**
     * `default` | `discrete` | `quadrature` | `montecarlo`, MATLAB
     * `UQ.listValidMethods`. The first three all resolve to the quadrature
     * design: a discrete Prior is expanded as given, a continuous one is placed
     * at stratum medians. Only `montecarlo` draws.
     */
    std::string method = "default";
    /**
     * Nodes per continuous Prior, or design points under `montecarlo`; the
     * reference's `options.samples`, defaulted to 11 rather than to the
     * simulation-oriented default of `Solver.defaultOptions` because each node
     * is a full solver run.
     */
    std::size_t samples = lang::kPriorDefaultNodes;
    /** The Monte Carlo stream; unread by a quadrature design, which draws nothing. */
    unsigned long seed = 23000;
};

/** Where a Prior sits in the model, MATLAB's `priorInfo` entry. */
template <class T>
struct PriorSite {
    std::size_t node = 0;     ///< 1-based node index
    std::size_t station = 0;  ///< 1-based station index
    std::size_t cls = 0;      ///< 1-based class index
    /** True at a Source, where the Prior is an ARRIVAL process, not a service one. */
    bool arrival = false;
    /** The Prior itself, copied out of the service table. */
    lang::Distrib<T> prior;
};

/** One design point: a concrete distribution for every Prior, and its weight. */
template <class T>
struct UqDesignPoint {
    T weight = num_traits<T>::from_int(1);
    std::vector<lang::Distrib<T>> dists;  ///< one per site, in site order
};

/** What `solver_uq_run_analyzer` returns. */
template <class T>
struct UqSolution {
    /** The prior-weighted expectation of every metric, (nstations x nclasses). */
    mva::AvgResult<T> avg;
    /** The result at each design point, in design order. */
    std::vector<mva::AvgResult<T>> points;
    /** The design weights, summing to 1. */
    std::vector<T> weights;
    /** The alternatives each point substituted, one per site. */
    std::vector<UqDesignPoint<T>> design;
    /** Where the Priors were found. */
    std::vector<PriorSite<T>> sites;
    /** The RESOLVED discretization method: `quadrature` or `montecarlo`. */
    std::string method;
    /** The options the design was built with; `uq_interval` reads `samples`. */
    UqOptions options;
};

/**
 * `UQ.getUQMethod`: resolve the discretization method.
 *
 * `default` and `discrete` are aliases of `quadrature` and not separate rules:
 * a discrete Prior is already exact, so expanding it as given IS the quadrature
 * design for it, and the name survives only because the reference lists it.
 */
inline std::string uq_resolve_method(const std::string& m) {
    if (m.empty() || m == "default" || m == "discrete" || m == "quadrature") return "quadrature";
    if (m == "montecarlo") return "montecarlo";
    throw InputError("SolverUQ: unknown method '" + m +
                     "'; the valid names are default, discrete, quadrature and montecarlo");
}

/** `UQ.listValidMethods`. */
inline std::vector<std::string> uq_list_valid_methods() {
    return std::vector<std::string>{"default", "discrete", "quadrature", "montecarlo"};
}

/**
 * `UQ.getFeatureSet`: the one construct UQ adds, and nothing else.
 *
 * IT IS NOT USED AS A FEATURE GATE, here or in the reference: UQ solves nothing
 * itself, so the set of models it admits is the inner solver's, applied per
 * design point once the Prior is gone. It is declared because the registry is
 * the vocabulary in which a capability is stated, and "SolverUQ is the solver
 * that understands Prior" is a statement worth being able to make.
 *
 * THE ONE THING `UQ.supports` DOES DECIDE is whether the model carries a Prior
 * at all: a model with no uncertain parameter is not a UQ model, and its
 * posterior is a single design point equal to the point estimate the inner
 * solver already returns. `has_prior_distribution()` is that test here. MATLAB
 * used to return true unconditionally -- which made SolverAUTO offer every
 * 'uq.*' method name on an ordinary network -- and now gates on `UQ.modelHasPrior`,
 * matching the JAR (`detectPrior() != null`) and native python
 * (`hasPriorDistribution`).
 */
inline qn::FeatureSet uq_feature_set() {
    qn::FeatureSet f;
    f.set(qn::Feature::Prior);
    return f;
}

/**
 * `UQ.detectPriors`: find every Prior, in node order and then class order.
 *
 * SERVICE AT A QUEUE OR DELAY, ARRIVAL AT A SOURCE, which is the reference's own
 * pair of branches. A Prior anywhere else -- a Cache's read process, a
 * Transition's firing law -- is REFUSED by name rather than skipped: the
 * reference's loop would ignore it and then solve a model in which the
 * uncertainty silently became the mixture moments, which is a confident answer
 * to a question nobody asked.
 */
template <class T>
std::vector<PriorSite<T>> uq_detect_priors(const qn::NetworkStruct<T>& sn) {
    std::vector<PriorSite<T>> sites;
    for (std::size_t nd = 1; nd <= sn.nodes.size(); ++nd) {
        const std::size_t ist = sn.nodes[nd - 1].station;
        if (ist == 0 || ist > sn.service.size()) continue;
        const qn::NodeType ty = sn.nodes[nd - 1].nodetype;
        const bool servicer = (ty == qn::NodeType::Queue || ty == qn::NodeType::Delay);
        const bool source = (ty == qn::NodeType::Source);
        for (std::size_t r = 1; r <= sn.service[ist - 1].size(); ++r) {
            const lang::Distrib<T>& d = sn.service[ist - 1][r - 1];
            if (!d.is_prior()) continue;
            if (!servicer && !source)
                throw UnsupportedError(
                    "SolverUQ: node '" + sn.nodes[nd - 1].name +
                    "' carries a Prior, but a Prior is expanded only where the reference expands "
                    "one: the service process of a Queue or a Delay, or the arrival process of a "
                    "Source");
            PriorSite<T> s;
            s.node = nd;
            s.station = ist;
            s.cls = r;
            s.arrival = source;
            s.prior = d;
            sites.push_back(s);
        }
    }
    return sites;
}

namespace detail {

/**
 * `UQ.unrankIndex`: the linear index i in [0, prod(counts)) as a subscript
 * vector over a mixed-radix grid, FIRST COORDINATE VARYING FASTEST.
 *
 * The order is the reference's and is kept because it is what makes design
 * point k the same model in both codebases; any other unranking would permute
 * the per-point table while leaving the expectation unchanged, which is the
 * hardest kind of divergence to notice.
 */
inline std::vector<std::size_t> unrank_index(std::size_t i,
                                             const std::vector<std::size_t>& counts) {
    std::vector<std::size_t> idx(counts.size(), 0);
    std::size_t rem = i;
    for (std::size_t l = 0; l < counts.size(); ++l) {
        idx[l] = rem % counts[l];
        rem /= counts[l];
    }
    return idx;
}

}  // namespace detail

/**
 * `UQ.buildDesign`: reduce the detected Priors to weighted design points.
 *
 * With no Prior there is ONE point of weight 1 substituting nothing, so the
 * original model is solved once and the expectation is that solve -- the
 * degenerate case the reference also carries, and the reason a model without a
 * Prior is not an error here.
 */
template <class T>
std::vector<UqDesignPoint<T>> uq_build_design(const std::vector<PriorSite<T>>& sites,
                                              const UqOptions& opt) {
    const std::string method = uq_resolve_method(opt.method);
    const std::size_t n = opt.samples;
    if (n < 1) throw InputError("SolverUQ: options.samples must be at least 1");
    std::vector<UqDesignPoint<T>> design;
    if (sites.empty()) {
        design.push_back(UqDesignPoint<T>());
        return design;
    }
    const std::size_t L = sites.size();
    lang::PriorRng rng(opt.seed);

    if (method == "montecarlo") {
        // ALL PRIORS ARE DRAWN JOINTLY at each point, so the cost is n runs
        // whatever L is; that independence from L is the whole reason the
        // method exists beside the tensor product.
        const T w = T(num_traits<T>::from_int(1) / num_traits<T>::from_int(static_cast<long>(n)));
        for (std::size_t i = 0; i < n; ++i) {
            UqDesignPoint<T> p;
            p.weight = w;
            for (std::size_t l = 0; l < L; ++l) {
                const lang::PriorDesign<T> g =
                    lang::prior_discretize(sites[l].prior, 1, "montecarlo", rng);
                p.dists.push_back(g.dists[0]);
            }
            design.push_back(p);
        }
        return design;
    }

    std::vector<lang::PriorDesign<T>> marg(L);
    std::vector<std::size_t> counts(L, 0);
    std::size_t total = 1;
    for (std::size_t l = 0; l < L; ++l) {
        marg[l] = lang::prior_discretize(sites[l].prior, n, "quadrature", rng);
        counts[l] = marg[l].dists.size();
        if (counts[l] == 0) throw NumericError("SolverUQ: a Prior discretized to no alternative");
        // SATURATE RATHER THAN OVERFLOW: the product of a dozen 11-node Priors
        // does not fit a size_t, and a wrapped one would pass the cap check.
        if (total > kMaxDesignPoints / counts[l]) {
            total = kMaxDesignPoints + 1;
            break;
        }
        total *= counts[l];
    }
    if (total > kMaxDesignPoints)
        throw UnsupportedError(
            "SolverUQ: the tensor-product design has " + std::to_string(total) +
            " points, above the limit of " + std::to_string(kMaxDesignPoints) +
            "; use method 'montecarlo', whose cost does not grow with the number of Priors, or "
            "lower options.samples");

    for (std::size_t i = 0; i < total; ++i) {
        const std::vector<std::size_t> idx = detail::unrank_index(i, counts);
        UqDesignPoint<T> p;
        p.weight = num_traits<T>::from_int(1);
        for (std::size_t l = 0; l < L; ++l) {
            p.dists.push_back(marg[l].dists[idx[l]]);
            p.weight = T(p.weight * marg[l].weights[idx[l]]);
        }
        design.push_back(p);
    }
    return design;
}

/**
 * What solves one design point: the C++ spelling of `@(m) SolverXXX(m)`.
 *
 * It takes the REFRESHED struct of the expanded model and returns the same
 * `AvgResult` every Network solver in this port returns, so the aggregation
 * below is one loop rather than one loop per solver.
 */
template <class T>
using UqStageSolver = std::function<mva::AvgResult<T>(const qn::NetworkStruct<T>&)>;

namespace detail {

/** M += w * A, sizing M from A on first use. */
template <class T>
void accumulate(Matrix<T>& M, const Matrix<T>& A, const T& w) {
    if (A.empty()) return;
    if (M.empty()) M = Matrix<T>(A.rows(), A.cols(), num_traits<T>::from_int(0));
    if (M.rows() != A.rows() || M.cols() != A.cols())
        throw NumericError("SolverUQ: two design points reported metrics of different shape");
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) M(i, j) = T(M(i, j) + w * A(i, j));
}

template <class T>
void accumulate(std::vector<T>& v, const std::vector<T>& a, const T& w) {
    if (a.empty()) return;
    if (v.empty()) v.assign(a.size(), num_traits<T>::from_int(0));
    if (v.size() != a.size())
        throw NumericError("SolverUQ: two design points reported vectors of different length");
    for (std::size_t i = 0; i < a.size(); ++i) v[i] = T(v[i] + w * a[i]);
}

}  // namespace detail

/**
 * `UQ.aggregateResults`: the prior-weighted expectation of the solved points.
 *
 * Free rather than private to the solver below because it is a pure function of
 * `points` and `weights`, and `post()` is not the only caller that has those: a
 * host that solved the design itself -- one point per process, one per machine --
 * aggregates with this and needs nothing else from the class.
 */
template <class T>
void uq_aggregate(UqSolution<T>& sol) {
    sol.avg = mva::AvgResult<T>();
    // The expectation, metric by metric. CN and XN are aggregated beside the
    // station matrices because the reference's field list carries them.
    for (std::size_t e = 0; e < sol.points.size(); ++e) {
        const mva::AvgResult<T>& r = sol.points[e];
        const T w = sol.weights[e];
        detail::accumulate(sol.avg.QN, r.QN, w);
        detail::accumulate(sol.avg.UN, r.UN, w);
        detail::accumulate(sol.avg.RN, r.RN, w);
        detail::accumulate(sol.avg.TN, r.TN, w);
        detail::accumulate(sol.avg.AN, r.AN, w);
        detail::accumulate(sol.avg.WN, r.WN, w);
        detail::accumulate(sol.avg.CN, r.CN, w);
        detail::accumulate(sol.avg.XN, r.XN, w);
    }
    if (sol.points.empty()) return;
    // THE AGGREGATE NAMES THE ALGORITHM ONLY IF EVERY POINT RAN THE SAME ONE. A
    // method that resolves per model -- SolverMVA's `default` picks `exact` on
    // one alternative and an AMVA on another -- would otherwise have the
    // expectation reported under one point's name.
    sol.avg.method = sol.points[0].method;
    sol.avg.actualmethod = sol.points[0].actualmethod;
    for (const mva::AvgResult<T>& r : sol.points)
        if (r.actualmethod != sol.avg.actualmethod) {
            sol.avg.actualmethod = "mixed";
            break;
        }
    for (std::size_t e = 0; e < sol.points.size(); ++e)
        if (!sol.points[e].warning.empty()) {
            sol.avg.warning =
                "design point " + std::to_string(e + 1) + ": " + sol.points[e].warning;
            break;
        }
}

/**
 * The ensemble surface of UQ: `@UQ`'s `EnsembleSolver` implementation.
 *
 * WHY A CLASS AND NOT ONLY `solver_uq_run_analyzer`. The reference is an
 * `EnsembleSolver`, and the lifecycle it implements -- `init`, `pre`, `analyze`,
 * `post`, `converged`, `finish` -- is not decoration: it is the surface a caller
 * uses to drive the design POINT BY POINT, to inspect the ensemble before
 * solving it, to substitute its own stage solver per point, or to distribute the
 * points and aggregate afterwards. A single `run` function can do none of those,
 * so the port carried the numbers of UQ without the way UQ is meant to be
 * driven. `SolverEnv` in `env/solver_env.h` keeps the same lifecycle for the
 * same reason.
 *
 * ONE DELIBERATE DEVIATION: the reference's `init` materializes the whole
 * ensemble (`self.ensemble{i}` is a deep model copy per design point) and this
 * one does not -- `expand(e)` builds the copy for point e on demand, because a
 * 4096-point design would otherwise hold 4096 Networks alive to solve them one
 * at a time. The expansion is a pure function of the design, so a caller that
 * wants the model of point e asks for it and gets exactly what the reference's
 * `self.ensemble{e}` holds.
 *
 * `net` is taken by non-const reference because the expansion COPIES it once per
 * design point and each copy is then refreshed -- `get_struct()` re-derives
 * rates, chains and visits, which is exactly what changing a service process
 * requires and is why the substitution is done on the model rather than on a
 * struct.
 */
template <class T>
class SolverUq {
public:
    SolverUq(qn::Network<T>& n, const UqStageSolver<T>& s, const UqOptions& o = UqOptions())
        : net(&n), stage(s) {
        if (!stage)
            throw InputError(
                "SolverUQ: no stage solver was given. UQ solves nothing itself; it needs the "
                "solver that runs at each design point, the C++ spelling of UQ(model, "
                "@SolverMVA)");
        sol.options = o;
        init();
    }

    /** `UQ.init`: resolve the design; the ensemble itself is expanded per point. */
    void init() {
        sol.method = uq_resolve_method(sol.options.method);
        sol.sites = uq_detect_priors(net->get_struct());
        sol.design = uq_build_design(sol.sites, sol.options);
        sol.points.clear();
        sol.weights.clear();
        sol.avg = mva::AvgResult<T>();
    }

    /** `UQ.pre`: nothing to seed -- the design points are independent models. */
    void pre(int /*it*/) {}

    /**
     * `UQ.analyze`: solve design point `e` (1-BASED, as the reference indexes it).
     *
     * The result is APPENDED to `points`, so driving the ensemble in order gives
     * the same `points` vector `iterate()` builds. Calling it out of order is a
     * caller's choice and the order of `points` then follows the calls, which is
     * why `weights` is appended beside it rather than indexed into.
     */
    const mva::AvgResult<T>& analyze(int /*it*/, std::size_t e) {
        if (e < 1 || e > sol.design.size())
            throw InputError("SolverUQ::analyze: design point " + std::to_string(e) +
                             " is outside the design of " + std::to_string(sol.design.size()) +
                             " points");
        qn::Network<T> copy = expand(e);
        sol.points.push_back(stage(copy.get_struct()));
        sol.weights.push_back(sol.design[e - 1].weight);
        return sol.points.back();
    }

    /** `UQ.post`: the prior-weighted expectation over the points solved so far. */
    void post(int /*it*/) { uq_aggregate(sol); }

    /** `UQ.finish`: nothing to release. */
    void finish() {}

    /** `UQ.converged`: one iteration is all UQ needs, the design being fixed. */
    bool converged(int it) const { return it >= 1; }

    /** `UQ.runAnalyzer`: the whole lifecycle, in the reference's order. */
    const UqSolution<T>& iterate() {
        init();
        int it = 1;
        pre(it);
        for (std::size_t e = 1; e <= sol.design.size(); ++e) analyze(it, e);
        post(it);
        if (!converged(it))
            throw NumericError("SolverUQ: the design did not converge in one iteration");
        finish();
        return sol;
    }

    /**
     * `self.ensemble{e}`: the model of design point `e` (1-based), Prior gone.
     *
     * The reference names each copy `<model>_alt<e>`; this port keeps the copy's
     * name because a `qn::Network` name is read back by the JSON writer and the
     * LQN bridge, and renaming it would put a name no model file carries into
     * whatever a caller does next with the expanded model.
     */
    qn::Network<T> expand(std::size_t e) const {
        if (e < 1 || e > sol.design.size())
            throw InputError("SolverUQ::expand: design point " + std::to_string(e) +
                             " is outside the design of " + std::to_string(sol.design.size()) +
                             " points");
        qn::Network<T> copy = *net;
        for (std::size_t l = 0; l < sol.sites.size(); ++l) {
            const PriorSite<T>& s = sol.sites[l];
            if (s.arrival)
                copy.set_arrival(s.node, s.cls, sol.design[e - 1].dists[l]);
            else
                copy.set_service(s.node, s.cls, sol.design[e - 1].dists[l]);
        }
        return copy;
    }

    /** `UQ.hasPriorDistribution`. */
    bool has_prior_distribution() const { return !sol.sites.empty(); }

    /** `UQ.getNumAlternatives`: design points, 1 when the model carries no Prior. */
    std::size_t get_num_alternatives() const { return sol.design.size(); }

    /** `UQ.getNumberOfModels`: the same count, under the EnsembleSolver's name. */
    std::size_t get_number_of_models() const { return sol.design.size(); }

    /** `UQ.getProbabilities`: the design weights, from the DESIGN not the solve. */
    std::vector<T> get_probabilities() const {
        std::vector<T> w;
        w.reserve(sol.design.size());
        for (std::size_t e = 0; e < sol.design.size(); ++e) w.push_back(sol.design[e].weight);
        return w;
    }

    /** `UQ.getUQNodes`: nodes per continuous Prior. */
    std::size_t get_uq_nodes() const { return sol.options.samples; }

    /** `UQ.getUQMethod`: the RESOLVED design name. */
    const std::string& get_uq_method() const { return sol.method; }

    /** `UQ.getEnsembleAvg`: the per-point results, in the order they were solved. */
    const std::vector<mva::AvgResult<T>>& get_ensemble_avg() const { return sol.points; }

    /** `UQ.getAvg`: the aggregate. Valid once `post()` has run. */
    const mva::AvgResult<T>& get_avg() const { return sol.avg; }

    /** Where the Priors were found, the reference's `priorInfo`. */
    const std::vector<PriorSite<T>>& get_prior_info() const { return sol.sites; }

    /** The design itself, one entry per point. */
    const std::vector<UqDesignPoint<T>>& get_design() const { return sol.design; }

    /** Everything above in one value, which is what the CLI and the tests read. */
    const UqSolution<T>& get_solution() const { return sol; }

private:
    qn::Network<T>* net;
    UqStageSolver<T> stage;
    UqSolution<T> sol;
};

/**
 * `UQ.runAnalyzer` as a free call: expand, solve every design point, aggregate.
 *
 * The one-shot spelling of `SolverUq::iterate`, kept because most callers want
 * the solution and not the lifecycle, and because it is what `line_cli.cpp` and
 * `uq_interval_run` reach for.
 */
template <class T>
UqSolution<T> solver_uq_run_analyzer(qn::Network<T>& net, const UqStageSolver<T>& stage,
                            const UqOptions& opt = UqOptions()) {
    SolverUq<T> s(net, stage, opt);
    return s.iterate();
}

// ---------------------------------------------------------------------------
// Posterior summaries
// ---------------------------------------------------------------------------

/** The metric matrix a name selects, MATLAB's @c res.Avg.(metric) field. */
template <class T>
const Matrix<T>& uq_metric_matrix(const mva::AvgResult<T>& r, const std::string& metric) {
    if (metric == "Q") return r.QN;
    if (metric == "U") return r.UN;
    if (metric == "R") return r.RN;
    if (metric == "T") return r.TN;
    if (metric == "A") return r.AN;
    if (metric == "W") return r.WN;
    throw InputError("SolverUQ: unknown metric '" + metric + "'; use Q, U, R, T, A or W");
}

/**
 * `UQ.getSamples`: the value of one metric at every design point, with weights.
 *
 * `ist` and `r` are 1-based, as everywhere in the readable surface of this port.
 */
template <class T>
std::vector<T> uq_samples(const UqSolution<T>& sol, const std::string& metric, std::size_t ist,
                          std::size_t r) {
    std::vector<T> vals;
    for (std::size_t e = 0; e < sol.points.size(); ++e) {
        const Matrix<T>& M = uq_metric_matrix(sol.points[e], metric);
        if (M.empty() || ist == 0 || r == 0 || ist > M.rows() || r > M.cols())
            throw InputError("SolverUQ: metric " + metric + " is unavailable at station " +
                             std::to_string(ist) + ", class " + std::to_string(r) +
                             " for design point " + std::to_string(e + 1));
        vals.push_back(M(ist - 1, r - 1));
    }
    return vals;
}

/** The weighted mean and variance of a metric over the design. */
template <class T>
struct UqMoments {
    T mean = num_traits<T>::from_int(0);
    T var = num_traits<T>::from_int(0);
};

/**
 * `UQ.getMoments`: the unconditional mean of Trivedi and Bobbio Eq. (3.68) and
 * the second moment of the same weighting.
 *
 * Both are EXACT for a discrete Prior and quadrature- or sample-approximate for
 * a continuous one, which is the only sense in which a variance over 11 stratum
 * medians is a variance.
 */
template <class T>
UqMoments<T> uq_moments(const UqSolution<T>& sol, const std::string& metric, std::size_t ist,
                        std::size_t r) {
    const std::vector<T> vals = uq_samples(sol, metric, ist, r);
    UqMoments<T> out;
    for (std::size_t e = 0; e < vals.size(); ++e) out.mean += T(sol.weights[e] * vals[e]);
    for (std::size_t e = 0; e < vals.size(); ++e) {
        const T dv = T(vals[e] - out.mean);
        out.var += T(sol.weights[e] * dv * dv);
    }
    return out;
}

/** The weighted empirical law of a metric, sorted ascending; MATLAB's `EmpiricalCDF`. */
template <class T>
struct UqEmpiricalCdf {
    std::vector<T> values;         ///< the metric at each design point, ascending
    std::vector<T> probabilities;  ///< the weight of each value, in the same order
    std::vector<T> cdf;            ///< the running sum of `probabilities`
};

/** `UQ.getPosteriorDist`: the posterior law of a metric across the design. */
template <class T>
UqEmpiricalCdf<T> uq_posterior_cdf(const UqSolution<T>& sol, const std::string& metric,
                                   std::size_t ist, std::size_t r) {
    const std::vector<T> vals = uq_samples(sol, metric, ist, r);
    std::vector<std::size_t> ord(vals.size());
    for (std::size_t i = 0; i < ord.size(); ++i) ord[i] = i;
    std::stable_sort(ord.begin(), ord.end(),
                     [&vals](std::size_t a, std::size_t b) { return vals[a] < vals[b]; });
    UqEmpiricalCdf<T> out;
    T acc = num_traits<T>::from_int(0);
    for (std::size_t k = 0; k < ord.size(); ++k) {
        out.values.push_back(vals[ord[k]]);
        out.probabilities.push_back(sol.weights[ord[k]]);
        acc += sol.weights[ord[k]];
        out.cdf.push_back(acc);
    }
    return out;
}

/**
 * `UQ.getCredibleInterval`: the equal-tailed interval of the weighted empirical
 * law at coverage `level`.
 *
 * The endpoints are DESIGN-POINT VALUES, not interpolations between them: the
 * design is a finite set of models and the interval names two of them, which is
 * what the reference's `find(cw >= alpha, 1)` returns. On a coarse design the
 * interval is therefore conservative rather than smooth.
 */
template <class T>
std::pair<T, T> uq_credible_interval(const UqSolution<T>& sol, const std::string& metric,
                                     std::size_t ist, std::size_t r, double level = 0.95) {
    if (!(level > 0.0) || !(level < 1.0))
        throw InputError("SolverUQ: the coverage level must lie strictly between 0 and 1");
    const UqEmpiricalCdf<T> ec = uq_posterior_cdf(sol, metric, ist, r);
    if (ec.values.empty()) throw InputError("SolverUQ: the design is empty");
    T tot = num_traits<T>::from_int(0);
    for (const T& w : ec.probabilities) tot += w;
    const double alpha = (1.0 - level) / 2.0;
    T lo = ec.values.front(), hi = ec.values.back();
    bool lo_set = false, hi_set = false;
    for (std::size_t k = 0; k < ec.values.size(); ++k) {
        const double cw = num_traits<T>::to_double(T(ec.cdf[k] / tot));
        if (!lo_set && cw >= alpha) {
            lo = ec.values[k];
            lo_set = true;
        }
        if (!hi_set && cw >= 1.0 - alpha) {
            hi = ec.values[k];
            hi_set = true;
        }
    }
    return std::make_pair(lo, hi);
}

// ---------------------------------------------------------------------------
// Support-only (interval) uncertainty
// ---------------------------------------------------------------------------

/**
 * `UQ.getInterval`: the RANGE of every metric over the support of the Priors.
 *
 * A different epistemic question from the expectation above, and the one to ask
 * when a parameter can be BOUNDED but not distributed: the weights are dropped
 * and only the endpoints are kept. `exact` says which of the two regimes
 * produced it, and the distinction is not a quality label but a change of
 * meaning -- see `uq_interval`.
 *
 * THE INTERVAL IS CONDITIONAL on the true parameters lying inside the Prior
 * supports. It is not a bound on the exact solution of the network, and it must
 * not be composed with the brackets of SolverBA, which bracket the exact
 * solution of a model whose parameters are known.
 */
template <class T>
struct UqInterval {
    /** (nstations x nclasses) lower and upper endpoints of each metric. */
    Matrix<T> Qlo, Qup, Ulo, Uup, Rlo, Rup, Tlo, Tup, Wlo, Wup;
    /** System throughput and total response time; the EXACT path only. */
    T Xlo = num_traits<T>::from_int(0), Xup = num_traits<T>::from_int(0);
    T Rtot_lo = num_traits<T>::from_int(0), Rtot_up = num_traits<T>::from_int(0);
    bool has_totals = false;
    /** True when the interval is the attained hull rather than a sampled range. */
    bool exact = false;
    /** `mvainterval` or `sampled`. */
    std::string method;
    /** On the sampled path, the condition that disqualified the exact one. */
    std::string why;
};

/**
 * `UQ.priorMeanRange`: the range of a Prior's MEAN over its alternatives.
 *
 * Exact for a discrete Prior, whose alternatives ARE the support. A continuous
 * Prior is discretized first, so the range is that of the discretized support:
 * an unbounded parameter density is never reached at its tails, which is
 * exactly why the interval built from it is an inner approximation.
 */
template <class T>
std::pair<T, T> uq_prior_mean_range(const lang::Distrib<T>& prior, std::size_t n) {
    lang::PriorRng rng(0);
    const lang::PriorDesign<T> g = lang::prior_discretize(prior, n, "quadrature", rng);
    if (g.dists.empty()) throw NumericError("uq_prior_mean_range: the Prior has no alternative");
    T lo = g.dists[0].mean, up = g.dists[0].mean;
    for (const lang::Distrib<T>& d : g.dists) {
        if (d.mean < lo) lo = d.mean;
        if (up < d.mean) up = d.mean;
    }
    return std::make_pair(lo, up);
}

/**
 * `UQ.qualifiesForIntervalMVA`: whether the monotonicity theorems behind
 * `pfqn_mva_interval` hold for this model.
 *
 * The returned string NAMES the first violated condition rather than reporting
 * a bare false, because "this model does not qualify" leaves the modeller
 * guessing which of six conditions to change.
 */
template <class T>
std::pair<bool, std::string> uq_qualifies_for_interval_mva(const qn::NetworkStruct<T>& sn,
                                                           const std::vector<PriorSite<T>>& sites) {
    for (const PriorSite<T>& s : sites)
        if (s.arrival)
            return std::make_pair(false, "a Prior sits on an arrival process, so the model is open");
    if (sn.nclasses != 1)
        return std::make_pair(false, "the theorems are proved for a single class only");
    if (!(sn.nclosedjobs() > 0.0)) return std::make_pair(false, "the class is not closed");
    if (sn.nodes.size() != sn.nstations)
        return std::make_pair(false, "the model has nodes that are not stations");
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        const bool inf = sn.stations[i].sched == lang::SchedStrategy::INF;
        if (!inf && sn.stations[i].nservers > 1.0)
            return std::make_pair(false, "a queueing station has more than one server");
        if (!inf && sn.stations[i].sched != lang::SchedStrategy::PS &&
            sn.stations[i].sched != lang::SchedStrategy::FCFS)
            return std::make_pair(false, "a station is neither delay, PS nor FCFS");
    }
    for (std::size_t a = 0; a < sites.size(); ++a)
        for (std::size_t b = a + 1; b < sites.size(); ++b)
            if (sites[a].station == sites[b].station)
                return std::make_pair(false, "two Priors sit on the same station");
    return std::make_pair(true, std::string());
}

/**
 * `UQ.intervalByMVA`: the exact hull through `pfqn_mva_interval`.
 *
 * The demand box is the nominal demand vector with the prior-carrying stations
 * widened to the range of mean service times over the Prior support. NO
 * ENSEMBLE RUN HAPPENS: 2*(m+2) MVA calls replace the whole tensor design, and
 * the answer is the attained range rather than the range of what was sampled.
 */
template <class T>
UqInterval<T> uq_interval_by_mva(const qn::NetworkStruct<T>& sn,
                                 const std::vector<PriorSite<T>>& sites, std::size_t nodes) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = sn.nstations;
    if (sn.visits.empty()) throw NumericError("uq_interval_by_mva: the model carries no visits");
    std::vector<T> V(M, zero), STlo(M, zero), STup(M, zero);
    std::vector<bool> inf(M, false);
    for (std::size_t i = 0; i < M; ++i) {
        V[i] = sn.visits[0](i, 0);
        inf[i] = sn.stations[i].sched == lang::SchedStrategy::INF;
        // A disabled pair has a NaN rate; MATLAB's `ST(isnan(ST)) = 0` makes it
        // a station the class does not visit rather than an error.
        const double rate = num_traits<T>::to_double(sn.rates(i, 0));
        const T st = (std::isnan(rate) || rate == 0.0) ? zero : T(one / sn.rates(i, 0));
        STlo[i] = st;
        STup[i] = st;
    }
    for (const PriorSite<T>& s : sites) {
        const std::pair<T, T> mr = uq_prior_mean_range(s.prior, nodes);
        STlo[s.station - 1] = mr.first;
        STup[s.station - 1] = mr.second;
    }

    std::vector<std::size_t> qidx;
    T zlo = zero, zup = zero;
    for (std::size_t i = 0; i < M; ++i) {
        if (inf[i]) {
            zlo += T(V[i] * STlo[i]);
            zup += T(V[i] * STup[i]);
        } else {
            qidx.push_back(i);
        }
    }
    if (qidx.empty())
        throw UnsupportedError(
            "uq_interval_by_mva: the model has no queueing station, so there is no MVA recursion "
            "to bound; a pure delay model's metrics are the demand intervals themselves");

    Matrix<T> L(qidx.size(), 2, zero);
    for (std::size_t k = 0; k < qidx.size(); ++k) {
        L(k, 0) = T(V[qidx[k]] * STlo[qidx[k]]);
        L(k, 1) = T(V[qidx[k]] * STup[qidx[k]]);
    }
    const int n = static_cast<int>(std::lround(sn.nclosedjobs()));
    const pfqn::MvaIntervalResult<T> iv = pfqn::pfqn_mva_interval(L, n, n, zlo, zup);

    UqInterval<T> out;
    out.Qlo = Matrix<T>(M, 1, zero);
    out.Qup = Matrix<T>(M, 1, zero);
    out.Ulo = Matrix<T>(M, 1, zero);
    out.Uup = Matrix<T>(M, 1, zero);
    out.Rlo = Matrix<T>(M, 1, zero);
    out.Rup = Matrix<T>(M, 1, zero);
    out.Wlo = Matrix<T>(M, 1, zero);
    out.Wup = Matrix<T>(M, 1, zero);
    out.Tlo = Matrix<T>(M, 1, zero);
    out.Tup = Matrix<T>(M, 1, zero);
    for (std::size_t k = 0; k < qidx.size(); ++k) {
        const std::size_t i = qidx[k];
        out.Qlo(i, 0) = iv.Q(k, 0);
        out.Qup(i, 0) = iv.Q(k, 1);
        out.Ulo(i, 0) = iv.U(k, 0);
        out.Uup(i, 0) = iv.U(k, 1);
        out.Wlo(i, 0) = iv.R(k, 0);
        out.Wup(i, 0) = iv.R(k, 1);
        // The RESIDENCE time is per visit; the response time divides it out.
        out.Rlo(i, 0) = V[i] > zero ? T(iv.R(k, 0) / V[i]) : zero;
        out.Rup(i, 0) = V[i] > zero ? T(iv.R(k, 1) / V[i]) : zero;
    }
    for (std::size_t i = 0; i < M; ++i) {
        if (!inf[i]) continue;
        // A delay station never queues, so its residence time IS its own demand
        // interval and its population is the throughput times that demand,
        // enclosed as a product of two intervals.
        out.Wlo(i, 0) = T(V[i] * STlo[i]);
        out.Wup(i, 0) = T(V[i] * STup[i]);
        out.Rlo(i, 0) = STlo[i];
        out.Rup(i, 0) = STup[i];
        out.Qlo(i, 0) = T(iv.Xlo * V[i] * STlo[i]);
        out.Qup(i, 0) = T(iv.Xup * V[i] * STup[i]);
        out.Ulo(i, 0) = out.Qlo(i, 0);
        out.Uup(i, 0) = out.Qup(i, 0);
    }
    for (std::size_t i = 0; i < M; ++i) {
        out.Tlo(i, 0) = T(V[i] * iv.Xlo);
        out.Tup(i, 0) = T(V[i] * iv.Xup);
    }
    out.Xlo = iv.Xlo;
    out.Xup = iv.Xup;
    out.Rtot_lo = iv.Rtot_lo;
    out.Rtot_up = iv.Rtot_up;
    out.has_totals = true;
    out.exact = true;
    out.method = "mvainterval";
    return out;
}

/**
 * `UQ.intervalBySampling`: the range of each metric across the design points
 * that were actually solved.
 *
 * EXACT FOR A DISCRETE PRIOR, whose design visits the whole support, and an
 * INNER approximation for a continuous one, since a quadrature node is a
 * stratum median and never an endpoint. It is therefore not an enclosure, and
 * `exact` is false to say so.
 */
template <class T>
UqInterval<T> uq_interval_by_sampling(const UqSolution<T>& sol) {
    UqInterval<T> out;
    out.method = "sampled";
    out.exact = false;
    auto range = [&](const Matrix<T>& (*pick)(const mva::AvgResult<T>&), Matrix<T>& lo,
                     Matrix<T>& up) {
        for (const mva::AvgResult<T>& r : sol.points) {
            const Matrix<T>& v = pick(r);
            if (v.empty()) continue;
            if (lo.empty()) {
                lo = v;
                up = v;
                continue;
            }
            for (std::size_t i = 0; i < v.rows(); ++i)
                for (std::size_t j = 0; j < v.cols(); ++j) {
                    if (v(i, j) < lo(i, j)) lo(i, j) = v(i, j);
                    if (up(i, j) < v(i, j)) up(i, j) = v(i, j);
                }
        }
    };
    range([](const mva::AvgResult<T>& r) -> const Matrix<T>& { return r.QN; }, out.Qlo, out.Qup);
    range([](const mva::AvgResult<T>& r) -> const Matrix<T>& { return r.UN; }, out.Ulo, out.Uup);
    range([](const mva::AvgResult<T>& r) -> const Matrix<T>& { return r.RN; }, out.Rlo, out.Rup);
    range([](const mva::AvgResult<T>& r) -> const Matrix<T>& { return r.TN; }, out.Tlo, out.Tup);
    range([](const mva::AvgResult<T>& r) -> const Matrix<T>& { return r.WN; }, out.Wlo, out.Wup);
    return out;
}

/**
 * `UQ.getInterval`: the exact hull where the monotonicity theorems apply, the
 * sampled range otherwise.
 *
 * THE TWO PATHS DO NOT MEAN THE SAME THING and the caller must read `exact`
 * before quoting the numbers. The MVA path returns the ATTAINED range over the
 * whole (continuous) demand box; the sampled path returns the range over the
 * points that happened to be solved, which for a continuous Prior lies strictly
 * inside the true range. `why` carries the condition that forced the fallback,
 * which is the reference's `line_warning` text made into a returned value --
 * this port has no warning channel, and a range that is not an enclosure must
 * not be silently indistinguishable from one that is.
 */
template <class T>
UqInterval<T> uq_interval(const UqSolution<T>& sol, const qn::NetworkStruct<T>& sn) {
    const std::pair<bool, std::string> q = uq_qualifies_for_interval_mva(sn, sol.sites);
    if (q.first) return uq_interval_by_mva(sn, sol.sites, sol.options.samples);
    UqInterval<T> out = uq_interval_by_sampling(sol);
    out.why = q.second;
    return out;
}

/**
 * `getInterval` from the model, running the ensemble ONLY when it is needed.
 *
 * The exact path costs 2*(m+2) MVA calls and reads no design point, so a caller
 * who wants the range and not the expectation should not pay for the tensor
 * design: `UQ.getInterval` calls `intervalByMVA` without touching `self.results`
 * and only the sampling fallback calls `iterate`. The stage solver is therefore
 * never invoked on a qualifying model, which also means a model whose stage
 * solver would refuse it still has a computable interval.
 */
template <class T>
UqInterval<T> uq_interval_run(qn::Network<T>& net, const UqStageSolver<T>& stage,
                              const UqOptions& opt = UqOptions()) {
    const std::vector<PriorSite<T>> sites = uq_detect_priors(net.get_struct());
    const std::pair<bool, std::string> q = uq_qualifies_for_interval_mva(net.get_struct(), sites);
    if (q.first) return uq_interval_by_mva(net.get_struct(), sites, opt.samples);
    UqInterval<T> out = uq_interval_by_sampling(solver_uq_run_analyzer(net, stage, opt));
    out.why = q.second;
    return out;
}

}  // namespace uq
}  // namespace line

#endif  // LINE_SOLVERS_UQ_SOLVER_UQ_H
