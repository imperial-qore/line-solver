/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_WF_WF_AUTO_INTEGRATION_H
#define LINE_API_WF_WF_AUTO_INTEGRATION_H

/**
 * Workflow-aware solver recommendation: the port of
 * jar/src/main/java/jline/api/wf/Wf_auto_integration.java.
 *
 * The routine takes the pattern analysis of `wf_analyzer.h`, turns it into a
 * flat feature vector, and lets those features amend the recommendation the
 * AUTO heuristic would give from the model alone. It is a HEURISTIC, not an
 * analysis: the numbers below (0.7 seed confidence, the 5-way parallelism
 * threshold, the 0.8 and 0.5 loop-probability bands, the 1.5 and 0.5 entropy
 * bands, the 50-node size gate) are the reference's constants and are ported
 * verbatim, because a caller comparing the two codebases compares the ADVICE,
 * and advice that differs by a tuned constant is a different answer.
 *
 * TWO DELIBERATE DEPARTURES FROM THE JAVA, both structural rather than
 * behavioural:
 *
 * 1. The feature bag is a TYPED STRUCT, not a `Map<String, Object>`. The Java
 *    reads each feature back with an `instanceof` test and substitutes 0 when
 *    the cast fails, so a mistyped or missing key degrades silently into a
 *    neutral value. Every one of those defaults is reproduced here by the
 *    struct's initializer, and the `present` flags below stand in for the
 *    Java's key-absent case -- the arms guarded by `!patterns.X.isEmpty()`
 *    never populate their keys otherwise.
 *
 * 2. `createOptimalSolver` returns a CHOICE, not a constructed solver. The C++
 *    solvers are free functions over a NetworkStruct rather than a runtime
 *    polymorphic family, so there is no `NetworkSolver` to hand back; the
 *    caller dispatches on the enum. Returning the choice also keeps `api/` from
 *    depending on `solvers/`, which is the tree's layering.
 *
 * The JAR reaches `branchNodes` by REFLECTION in generatePatternInsights, which
 * is a workaround for its own generics and silently yields 0 on any exception;
 * here it is a plain member read, so the insight fires when it should.
 *
 * ARITHMETIC: field, plus what the branch entropy needs (gated in
 * wf_branch_detector.h).
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/wf/wf_analyzer.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace wf {

/** The solvers the reference chooses among. */
enum class WfSolver { MVA = 0, NC, SSA, FLUID, JMT, CTMC, AUTO };

/** The reference's own spelling of each choice. */
inline std::string wf_solver_name(WfSolver s) {
    switch (s) {
        case WfSolver::MVA: return "MVA";
        case WfSolver::NC: return "NC";
        case WfSolver::SSA: return "SSA";
        case WfSolver::FLUID: return "FLUID";
        case WfSolver::JMT: return "JMT";
        case WfSolver::CTMC: return "CTMC";
        default: return "AUTO";
    }
}

/**
 * The flat feature vector the recommendation reads.
 *
 * The four `has*` flags decide which of the remaining groups carry meaning;
 * the rest keep the Java's zero defaults so an unpopulated group is neutral.
 */
template <class T>
struct WorkflowFeatures {
    bool hasSequencePatterns = false;
    bool hasParallelPatterns = false;
    bool hasLoopPatterns = false;
    bool hasBranchPatterns = false;
    std::size_t numSequences = 0, numParallels = 0, numLoops = 0, numBranches = 0;
    std::size_t originalNodeCount = 0, originalLinkCount = 0;
    std::size_t optimizedNodeCount = 0, optimizedLinkCount = 0;
    T avgSequenceLength = num_traits<T>::from_int(0);
    std::size_t maxSequenceLength = 0;
    T avgParallelism = num_traits<T>::from_int(0);
    std::size_t maxParallelism = 0;
    T avgLoopProbability = num_traits<T>::from_int(0);
    T maxLoopProbability = num_traits<T>::from_int(0);
    T avgBranches = num_traits<T>::from_int(0);
    std::size_t maxBranches = 0;
    T avgBranchEntropy = num_traits<T>::from_int(0);
};

/** What the recommendation returns: the choice and why. */
template <class T>
struct ExtendedSolverRecommendation {
    WfSolver recommendedSolver = WfSolver::MVA;
    double confidence = 0.0;
    std::vector<std::string> reasoning;
    WorkflowFeatures<T> workflowFeatures;
    std::vector<WfSolver> alternativeSolvers;  ///< at most three, the reference's cap
};

/** The four model facts the base recommendation reads, i.e. ModelAnalyzer. */
template <class T>
struct WfModelFacts {
    bool hasProductForm = false;
    bool hasSingleChain = false;
    bool hasMultiChain = false;
    long totalJobs = 0;
    double avgJobsPerChain = 0.0;
};

/**
 * Read the base facts off a NetworkStruct.
 *
 * The job total SKIPS the infinite entries, which is how the Java avoids
 * `(int) Double.POSITIVE_INFINITY` becoming MAX_VALUE and swamping the count:
 * an open class contributes nothing to a population-based threshold.
 */
template <class T>
WfModelFacts<T> wf_model_facts(const qn::NetworkStruct<T>& sn) {
    WfModelFacts<T> f;
    f.hasProductForm = sn.has_product_form();
    f.hasSingleChain = (sn.nchains == 1);
    f.hasMultiChain = (sn.nchains > 1);
    const std::vector<double> N = sn.njobs();
    for (std::size_t i = 0; i < N.size(); ++i)
        if (!std::isinf(N[i]) && !std::isnan(N[i])) f.totalJobs += static_cast<long>(N[i]);
    f.avgJobsPerChain =
        (sn.nchains > 0) ? static_cast<double>(f.totalJobs) / static_cast<double>(sn.nchains) : 0.0;
    return f;
}

/** The AUTO heuristic before the workflow features amend it. */
template <class T>
WfSolver wf_base_recommendation(const WfModelFacts<T>& f) {
    if (f.hasSingleChain) return WfSolver::NC;
    if (f.hasMultiChain && f.hasProductForm && f.totalJobs < 10) return WfSolver::NC;
    if (f.hasMultiChain && f.hasProductForm && f.avgJobsPerChain < 30) return WfSolver::MVA;
    if (f.hasMultiChain && f.avgJobsPerChain > 30) return WfSolver::FLUID;
    return WfSolver::MVA;
}

/** Flatten the analysis into the feature vector the heuristic reads. */
template <class T>
WorkflowFeatures<T> extract_workflow_features(const WorkflowAnalysis<T>& a) {
    WorkflowFeatures<T> f;
    const DetectedPatterns<T>& p = a.detectedPatterns;
    f.hasSequencePatterns = !p.sequences.empty();
    f.hasParallelPatterns = !p.parallels.empty();
    f.hasLoopPatterns = !p.loops.empty();
    f.hasBranchPatterns = !p.branches.empty();
    f.numSequences = p.sequences.size();
    f.numParallels = p.parallels.size();
    f.numLoops = p.loops.size();
    f.numBranches = p.branches.size();

    f.originalNodeCount = a.statistics.originalComplexity.totalNodes;
    f.originalLinkCount = a.statistics.originalComplexity.totalLinks;
    f.optimizedNodeCount = a.statistics.optimizedComplexity.totalNodes;
    f.optimizedLinkCount = a.statistics.optimizedComplexity.totalLinks;

    if (f.hasSequencePatterns) {
        T total = num_traits<T>::from_int(0);
        std::size_t mx = 0;
        for (std::size_t i = 0; i < p.sequences.size(); ++i) {
            total += num_traits<T>::from_int(static_cast<int>(p.sequences[i].size()));
            mx = std::max(mx, p.sequences[i].size());
        }
        f.avgSequenceLength = total / num_traits<T>::from_int(static_cast<int>(p.sequences.size()));
        f.maxSequenceLength = mx;
    }
    if (f.hasParallelPatterns) {
        T total = num_traits<T>::from_int(0);
        std::size_t mx = 0;
        for (std::size_t i = 0; i < p.parallels.size(); ++i) {
            total += num_traits<T>::from_int(static_cast<int>(p.parallels[i].size()));
            mx = std::max(mx, p.parallels[i].size());
        }
        f.avgParallelism = total / num_traits<T>::from_int(static_cast<int>(p.parallels.size()));
        f.maxParallelism = mx;
    }
    if (f.hasLoopPatterns) {
        f.avgLoopProbability = a.statistics.loopStats.avgLoopProbability;
        f.maxLoopProbability = a.statistics.loopStats.maxLoopProbability;
    }
    if (f.hasBranchPatterns) {
        f.avgBranches = a.statistics.branchStats.avgBranches;
        f.maxBranches = a.statistics.branchStats.maxBranches;
        T total = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < p.branches.size(); ++i)
            total += calculate_branch_diversity(p.branches[i]).entropy;
        f.avgBranchEntropy = total / num_traits<T>::from_int(static_cast<int>(p.branches.size()));
    }
    return f;
}

namespace detail {

/** Append a solver unless it is the recommendation or already listed. */
inline void push_alternative(std::vector<WfSolver>* alt, WfSolver rec, WfSolver s) {
    if (s == rec) return;
    for (std::size_t i = 0; i < alt->size(); ++i)
        if ((*alt)[i] == s) return;
    alt->push_back(s);
}

}  // namespace detail

/**
 * Amend the base recommendation with what the workflow analysis found.
 *
 * The reasoning strings are user-facing and reproduced verbatim, so a caller
 * grepping them sees the same text from either codebase.
 */
template <class T>
ExtendedSolverRecommendation<T> enhance_recommendation_with_workflow(
    WfSolver base, const WorkflowFeatures<T>& f, const WorkflowAnalysis<T>& a) {
    ExtendedSolverRecommendation<T> r;
    r.workflowFeatures = f;
    WfSolver rec = base;
    double confidence = 0.7;

    if (f.hasSequencePatterns) {
        r.reasoning.push_back("Detected sequence patterns - suitable for analytical methods");
        confidence += 0.1;
        if (f.maxSequenceLength > 10) {
            r.reasoning.push_back("Long sequences detected - consider FLUID approximation");
            if (rec == WfSolver::MVA) r.alternativeSolvers.push_back(WfSolver::FLUID);
        }
    }
    if (f.hasParallelPatterns) {
        r.reasoning.push_back("Detected parallel patterns - fork-join structures present");
        if (f.maxParallelism > 5) {
            r.reasoning.push_back(
                "High parallelism detected - exact methods may be computationally expensive");
            if (rec == WfSolver::NC || rec == WfSolver::MVA) {
                rec = WfSolver::SSA;
                r.reasoning.push_back("Switching to SSA for high-parallelism workflow");
            }
            r.alternativeSolvers.push_back(WfSolver::JMT);
        } else {
            confidence += 0.05;
        }
    }
    if (f.hasLoopPatterns) {
        const double avgLoop = num_traits<T>::to_double(f.avgLoopProbability);
        const double maxLoop = num_traits<T>::to_double(f.maxLoopProbability);
        r.reasoning.push_back("Detected loop patterns with avg probability " +
                              std::to_string(avgLoop));
        if (maxLoop > 0.8) {
            r.reasoning.push_back(
                "High loop probability detected - may cause numerical instability");
            confidence -= 0.1;
            if (rec == WfSolver::NC || rec == WfSolver::MVA) {
                r.alternativeSolvers.push_back(WfSolver::SSA);
                r.alternativeSolvers.push_back(WfSolver::FLUID);
            }
        } else if (maxLoop > 0.5) {
            r.reasoning.push_back("Moderate loop probability - analytical methods suitable");
            confidence += 0.05;
        }
    }
    if (f.hasBranchPatterns) {
        const double avgEntropy = num_traits<T>::to_double(f.avgBranchEntropy);
        r.reasoning.push_back("Detected branch patterns with avg entropy " +
                              std::to_string(avgEntropy));
        if (avgEntropy > 1.5) {
            r.reasoning.push_back("High branching entropy - complex decision structure");
            if (f.maxBranches > 5) {
                r.reasoning.push_back("Many branches detected - consider simulation methods");
                r.alternativeSolvers.push_back(WfSolver::JMT);
                r.alternativeSolvers.push_back(WfSolver::SSA);
            }
        }
        if (avgEntropy < 0.5) {
            r.reasoning.push_back("Low branching entropy - deterministic-like behavior");
            confidence += 0.1;
        }
    }

    if (f.originalNodeCount > f.optimizedNodeCount) {
        const double reduction = static_cast<double>(f.originalNodeCount - f.optimizedNodeCount) /
                                 static_cast<double>(f.originalNodeCount);
        r.reasoning.push_back("Workflow complexity reduced by " +
                              std::to_string(static_cast<int>(reduction * 100)) +
                              "% through pattern optimization");
        confidence += 0.1;
    }
    if (f.optimizedNodeCount > 50) {
        r.reasoning.push_back("Large optimized workflow - consider approximation methods");
        if (rec == WfSolver::NC) {
            rec = WfSolver::MVA;
            r.reasoning.push_back("Switching from NC to MVA for large workflow");
        }
        r.alternativeSolvers.push_back(WfSolver::FLUID);
    }
    (void)a;  // the reference keeps the analysis in scope without reading it here

    confidence = std::min(1.0, std::max(0.1, confidence));

    // The reference pads with the standard list IN ITS ORDER and then keeps the
    // first three, so the padding order is part of the answer.
    const WfSolver standard[6] = {WfSolver::MVA,   WfSolver::NC,  WfSolver::SSA,
                                  WfSolver::FLUID, WfSolver::JMT, WfSolver::CTMC};
    for (std::size_t i = 0; i < 6; ++i) detail::push_alternative(&r.alternativeSolvers, rec, standard[i]);
    if (r.alternativeSolvers.size() > 3) r.alternativeSolvers.resize(3);

    r.recommendedSolver = rec;
    r.confidence = confidence;
    return r;
}

/** The entry point: analyse the workflow, then let it amend the base choice. */
template <class T>
ExtendedSolverRecommendation<T> recommend_solver_with_workflow_analysis(
    const qn::NetworkStruct<T>& sn) {
    const WorkflowRepresentation<T> w = wf_from_struct(sn);
    const WorkflowAnalysis<T> a = analyze_workflow(w);
    const WorkflowFeatures<T> f = extract_workflow_features(a);
    return enhance_recommendation_with_workflow(wf_base_recommendation(wf_model_facts(sn)), f, a);
}

/** The chosen solver alone, i.e. the reference's createOptimalSolver. */
template <class T>
WfSolver create_optimal_solver(const qn::NetworkStruct<T>& sn) {
    return recommend_solver_with_workflow_analysis(sn).recommendedSolver;
}

/** The advisory text the reference's getOptimizationInsights assembles. */
struct OptimizationInsights {
    std::vector<std::string> recommendations;
    std::vector<std::string> patternInsights;
    std::vector<std::string> performancePredictions;
};

/** Pattern-level advice; the strings are the reference's, verbatim. */
template <class T>
std::vector<std::string> generate_pattern_insights(const DetectedPatterns<T>& p) {
    std::vector<std::string> out;
    if (!p.sequences.empty())
        out.push_back("Consider merging sequential services to reduce overhead");
    if (!p.parallels.empty())
        out.push_back("Parallel patterns can benefit from resource pooling strategies");
    if (!p.loops.empty())
        out.push_back("High-probability loops may benefit from caching or memoization");
    if (!p.branches.empty()) {
        double total = 0.0;
        for (std::size_t i = 0; i < p.branches.size(); ++i)
            total += static_cast<double>(p.branches[i].branchNodes.size());
        if (total / static_cast<double>(p.branches.size()) > 3.0)
            out.push_back("Complex branching patterns - consider load balancing strategies");
    }
    return out;
}

/** Solve-time advice keyed off the collapse ratio and the pattern mix. */
template <class T>
std::vector<std::string> generate_performance_predictions(const WorkflowAnalysis<T>& a) {
    std::vector<std::string> out;
    const double ratio = num_traits<T>::to_double(a.statistics.updateStats.reductionRatio);
    if (ratio > 0.1)
        out.push_back("Expected " + std::to_string(static_cast<int>(ratio * 100)) +
                      "% reduction in solve time");
    if (!a.detectedPatterns.parallels.empty())
        out.push_back("High potential for parallel execution optimization");
    if (!a.detectedPatterns.loops.empty())
        out.push_back("Loop patterns may affect solver convergence rates");
    return out;
}

/** All three advisory blocks for one model. */
template <class T>
OptimizationInsights get_optimization_insights(const qn::NetworkStruct<T>& sn) {
    const WorkflowAnalysis<T> a = analyze_workflow(wf_from_struct(sn));
    OptimizationInsights out;
    out.recommendations = get_optimization_recommendations(a);
    out.patternInsights = generate_pattern_insights(a.detectedPatterns);
    out.performancePredictions = generate_performance_predictions(a);
    return out;
}

/**
 * The reference's self-check: a usable recommendation over a valid analysis.
 *
 * The Java swallows every exception and returns false; here the analysis is
 * total over a well-formed struct, so a throw is a defect and is left to
 * propagate rather than being reported as a failed validation.
 */
template <class T>
bool validate_workflow_enhancement(const qn::NetworkStruct<T>& sn) {
    const ExtendedSolverRecommendation<T> r = recommend_solver_with_workflow_analysis(sn);
    const WorkflowAnalysis<T> a = analyze_workflow(wf_from_struct(sn));
    return r.confidence > 0.1 && !r.reasoning.empty() && validate_analysis(a);
}

}  // namespace wf
}  // namespace line

#endif  // LINE_API_WF_WF_AUTO_INTEGRATION_H
