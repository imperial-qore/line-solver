/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_WF_WORKFLOW_MANAGER_H
#define LINE_API_WF_WORKFLOW_MANAGER_H

/**
 * The workflow facade: the port of
 * jar/src/main/java/jline/api/wf/WorkflowManager.java.
 *
 * It composes `wf_analyzer.h` and `wf_auto_integration.h` into the things a
 * user asks for -- one analysis object, a complexity report, efficiency
 * metrics, a benchmark table and three export formats. There is no new
 * queueing content here; what must be faithful is the arithmetic of the scores
 * and the exact bytes of the exports, since both are what a caller compares
 * across codebases.
 *
 * THREE DEPARTURES FROM THE JAVA, each named:
 *
 * 1. `factorial` is computed in DOUBLE, not int. The Java's `int factorial(int)`
 *    silently overflows at n = 13 and goes NEGATIVE at n = 17, so a workflow
 *    with a wide fork would report a negative parallel complexity and a
 *    complexity score that DROPS as the model grows. That is a defect, not a
 *    convention, so it is not reproduced; the port saturates at infinity in
 *    double instead, which keeps the score monotone. Below 13 branches the two
 *    agree exactly.
 *
 * 2. `benchmarkSolvers` takes a CALLER-SUPPLIED runner. The C++ solvers are
 *    free functions over a NetworkStruct rather than a runtime polymorphic
 *    family, and `api/` does not depend on `solvers/` in this tree. Everything
 *    the Java method actually does -- time the call, record success, aggregate
 *    the queue lengths, turn a thrown exception into a failed row -- is here;
 *    only the construction of the solver is injected.
 *
 * 3. `validateWorkflow` reads the node count and the routing off the
 *    NetworkStruct rather than the Network object, which is where the port
 *    keeps that information.
 *
 * The Java's `calculateLoopEfficiency` ignores its own `linkMatrix` argument
 * and returns the constant 0.5, which its comment admits ("Simplified"). It is
 * ported as written, because changing it would change every metric a caller
 * has recorded; the constant is documented at the function instead.
 *
 * ARITHMETIC: field, plus what the branch entropy needs.
 */

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <functional>
#include <exception>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "line/api/wf/wf_analyzer.h"
#include "line/api/wf/wf_auto_integration.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace wf {

/** The four efficiency scores plus the two headline numbers. */
struct WorkflowPerformanceMetrics {
    double complexityReduction = 0.0;
    double solverConfidence = 0.0;
    double sequenceEfficiency = 1.0;
    double parallelEfficiency = 1.0;
    double loopEfficiency = 1.0;
    double branchEfficiency = 1.0;
};

/** Everything analyze_workflow_full returns, i.e. WorkflowAnalysisResult. */
template <class T>
struct WorkflowAnalysisResult {
    WorkflowAnalysis<T> patternAnalysis;
    ExtendedSolverRecommendation<T> solverRecommendation;
    OptimizationInsights optimizationInsights;
    WorkflowPerformanceMetrics performanceMetrics;
};

/** Per-pattern-family complexity, one entry of the report's patternComplexity. */
struct PatternComplexityEntry {
    std::size_t count = 0;
    std::size_t totalNodes = 0;  ///< totalBranches for the branch family
    double complexity = 0.0;
};

/** The complexity report, one band per pattern family plus the overall score. */
template <class T>
struct ComplexityReport {
    WorkflowComplexity<T> originalMetrics;
    WorkflowComplexity<T> optimizedMetrics;
    PatternComplexityEntry sequences, parallels, loops, branches;
    double overallComplexityScore = 0.0;
    std::string complexityLevel;  ///< LOW, MEDIUM, HIGH or VERY_HIGH
};

/** One row of the benchmark table. */
struct BenchmarkRow {
    bool success = false;
    double solveTimeMs = 0.0;
    bool hasResults = false;
    double totalQueueLength = 0.0;
    double avgQueueLength = 0.0;
    std::string error;
};

/** What validate_workflow reports. */
struct WorkflowValidation {
    bool isValid = false;
    std::vector<std::string> issues;
    std::size_t nodeCount = 0;
    bool hasRouting = false;
};

/** The export formats the facade offers. */
enum class WfExportFormat { Json = 0, Csv, Summary };

namespace detail {

/**
 * n! in double.
 *
 * The Java computes this in int and overflows at 13, going negative at 17; see
 * the header note. Double saturates to infinity instead, so the score stays
 * monotone in the fork width.
 */
inline double wf_factorial(std::size_t n) {
    double f = 1.0;
    for (std::size_t k = 2; k <= n; ++k) f *= static_cast<double>(k);
    return f;
}

}  // namespace detail

/** avg sequence length capped at 5, the reference's saturation point. */
inline double calculate_sequence_efficiency(const std::vector<std::vector<int>>& sequences) {
    if (sequences.empty()) return 1.0;
    std::size_t total = 0;
    for (std::size_t i = 0; i < sequences.size(); ++i) total += sequences[i].size();
    const double avg = static_cast<double>(total) / static_cast<double>(sequences.size());
    return std::min(1.0, avg / 5.0);
}

/** Efficiency falls linearly past two-way parallelism, floored at 0.1. */
inline double calculate_parallel_efficiency(const std::vector<std::vector<int>>& parallels) {
    if (parallels.empty()) return 1.0;
    std::size_t total = 0;
    for (std::size_t i = 0; i < parallels.size(); ++i) total += parallels[i].size();
    const double avg = static_cast<double>(total) / static_cast<double>(parallels.size());
    return std::max(0.1, 1.0 - (avg - 2.0) / 10.0);
}

/**
 * The reference's placeholder: any loop at all costs half the efficiency.
 *
 * It assumes a loop probability of 0.5 and never reads the link matrix, which
 * its own comment calls "Simplified". Ported as written -- a caller comparing
 * the two codebases compares this number.
 */
inline double calculate_loop_efficiency(const std::vector<int>& loops) {
    if (loops.empty()) return 1.0;
    return 1.0 - 0.5;
}

/** The mean NORMALIZED entropy over the branches, i.e. how evenly they split. */
template <class T>
double calculate_branch_efficiency(const std::vector<BranchPattern<T>>& branches) {
    if (branches.empty()) return 1.0;
    double sum = 0.0;
    for (std::size_t i = 0; i < branches.size(); ++i)
        sum += num_traits<T>::to_double(calculate_branch_diversity(branches[i]).normalizedEntropy);
    return sum / static_cast<double>(branches.size());
}

/** The weighted size-plus-pattern score behind the complexity level. */
template <class T>
double calculate_complexity_score(std::size_t nodes, std::size_t links,
                                  const DetectedPatterns<T>& p) {
    double score = static_cast<double>(nodes) + static_cast<double>(links) * 0.5;
    for (std::size_t i = 0; i < p.sequences.size(); ++i) {
        const double n = static_cast<double>(p.sequences[i].size());
        score += n * n * 0.1;
    }
    for (std::size_t i = 0; i < p.parallels.size(); ++i)
        score += detail::wf_factorial(p.parallels[i].size()) * 0.2;
    score += static_cast<double>(p.loops.size()) * 10.0;
    for (std::size_t i = 0; i < p.branches.size(); ++i)
        score += static_cast<double>(p.branches[i].branchNodes.size()) * 2.0;
    return score;
}

/** The six metrics of the analysis result. */
template <class T>
WorkflowPerformanceMetrics calculate_performance_metrics(
    const WorkflowAnalysis<T>& a, const ExtendedSolverRecommendation<T>& r) {
    WorkflowPerformanceMetrics m;
    m.complexityReduction = num_traits<T>::to_double(a.statistics.updateStats.reductionRatio);
    m.solverConfidence = r.confidence;
    m.sequenceEfficiency = calculate_sequence_efficiency(a.detectedPatterns.sequences);
    m.parallelEfficiency = calculate_parallel_efficiency(a.detectedPatterns.parallels);
    m.loopEfficiency = calculate_loop_efficiency(a.detectedPatterns.loops);
    m.branchEfficiency = calculate_branch_efficiency(a.detectedPatterns.branches);
    return m;
}

/** The facade's headline call: analysis, recommendation, insights, metrics. */
template <class T>
WorkflowAnalysisResult<T> analyze_workflow_full(const qn::NetworkStruct<T>& sn) {
    WorkflowAnalysisResult<T> out;
    out.patternAnalysis = analyze_workflow(wf_from_struct(sn));
    out.solverRecommendation = recommend_solver_with_workflow_analysis(sn);
    out.optimizationInsights = get_optimization_insights(sn);
    out.performanceMetrics =
        calculate_performance_metrics(out.patternAnalysis, out.solverRecommendation);
    return out;
}

/** The patterns alone. */
template <class T>
DetectedPatterns<T> get_pattern_analysis(const qn::NetworkStruct<T>& sn) {
    return analyze_workflow(wf_from_struct(sn)).detectedPatterns;
}

/** The recommendation strings alone. */
template <class T>
std::vector<std::string> get_workflow_recommendations(const qn::NetworkStruct<T>& sn) {
    return get_optimization_recommendations(analyze_workflow(wf_from_struct(sn)));
}

/** The complexity report, with the reference's four bands on the score. */
template <class T>
ComplexityReport<T> generate_complexity_report(const qn::NetworkStruct<T>& sn) {
    const WorkflowAnalysis<T> a = analyze_workflow(wf_from_struct(sn));
    const DetectedPatterns<T>& p = a.detectedPatterns;
    ComplexityReport<T> r;
    r.originalMetrics = a.statistics.originalComplexity;
    r.optimizedMetrics = a.statistics.optimizedComplexity;

    r.sequences.count = p.sequences.size();
    for (std::size_t i = 0; i < p.sequences.size(); ++i) {
        const double n = static_cast<double>(p.sequences[i].size());
        r.sequences.totalNodes += p.sequences[i].size();
        r.sequences.complexity += n * n;
    }
    r.parallels.count = p.parallels.size();
    for (std::size_t i = 0; i < p.parallels.size(); ++i) {
        r.parallels.totalNodes += p.parallels[i].size();
        r.parallels.complexity += detail::wf_factorial(p.parallels[i].size());
    }
    r.loops.count = p.loops.size();
    r.loops.complexity = static_cast<double>(p.loops.size()) * 10.0;
    r.branches.count = p.branches.size();
    for (std::size_t i = 0; i < p.branches.size(); ++i) {
        r.branches.totalNodes += p.branches[i].branchNodes.size();
        r.branches.complexity += static_cast<double>(p.branches[i].branchNodes.size());
    }

    r.overallComplexityScore = calculate_complexity_score(
        a.statistics.originalComplexity.totalNodes, a.statistics.originalComplexity.totalLinks, p);
    if (r.overallComplexityScore < 10.0) r.complexityLevel = "LOW";
    else if (r.overallComplexityScore < 50.0) r.complexityLevel = "MEDIUM";
    else if (r.overallComplexityScore < 100.0) r.complexityLevel = "HIGH";
    else r.complexityLevel = "VERY_HIGH";
    return r;
}

/**
 * Time each solver and aggregate its queue lengths.
 *
 * @param solvers the solvers to try, in order
 * @param runner  runs one solver and returns its QN; it may throw, and a throw
 *                becomes a failed row exactly as the Java's catch does
 */
template <class T>
std::vector<std::pair<WfSolver, BenchmarkRow>> benchmark_solvers(
    const std::vector<WfSolver>& solvers,
    const std::function<Matrix<T>(WfSolver)>& runner) {
    std::vector<std::pair<WfSolver, BenchmarkRow>> out;
    for (std::size_t i = 0; i < solvers.size(); ++i) {
        BenchmarkRow row;
        const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
        try {
            const Matrix<T> QN = runner(solvers[i]);
            const std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
            row.success = true;
            row.solveTimeMs =
                std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(t1 - t0)
                    .count();
            const std::size_t n = QN.rows() * QN.cols();
            row.hasResults = (n > 0);
            double sum = 0.0;
            for (std::size_t rr = 0; rr < QN.rows(); ++rr)
                for (std::size_t cc = 0; cc < QN.cols(); ++cc)
                    sum += num_traits<T>::to_double(QN(rr, cc));
            row.totalQueueLength = sum;
            row.avgQueueLength = (n > 0) ? sum / static_cast<double>(n) : 0.0;
        } catch (const std::exception& e) {
            row.success = false;
            const char* what = e.what();
            row.error = (what != nullptr && what[0] != '\0') ? what : "Unknown error";
        }
        out.push_back(std::make_pair(solvers[i], row));
    }
    return out;
}

/** The reference's default benchmark set. */
inline std::vector<WfSolver> default_benchmark_solvers() {
    std::vector<WfSolver> v;
    v.push_back(WfSolver::MVA);
    v.push_back(WfSolver::NC);
    v.push_back(WfSolver::SSA);
    v.push_back(WfSolver::FLUID);
    return v;
}

/** The model is well formed and both analyses validate. */
template <class T>
WorkflowValidation validate_workflow(const qn::NetworkStruct<T>& sn) {
    WorkflowValidation v;
    v.nodeCount = sn.nodes.size();
    v.hasRouting = (sn.rtnodes.rows() > 0 && sn.rtnodes.cols() > 0);
    if (sn.nodes.empty()) v.issues.push_back("Network has no nodes");
    if (!validate_analysis(analyze_workflow(wf_from_struct(sn))))
        v.issues.push_back("Workflow analysis validation failed");
    if (!validate_workflow_enhancement(sn))
        v.issues.push_back("Workflow-enhanced solver selection validation failed");
    v.isValid = v.issues.empty();
    return v;
}

namespace detail {

/** The Java's `String.format("%.2f", x)`. */
inline std::string fixed2(double x) {
    std::ostringstream s;
    s.setf(std::ios::fixed);
    s.precision(2);
    s << x;
    return s.str();
}

/** The Java's default double rendering, which drops a trailing `.0` never. */
inline std::string plain(double x) {
    std::ostringstream s;
    s << x;
    return s.str();
}

}  // namespace detail

/**
 * Render the analysis.
 *
 * The three layouts are byte-for-byte the reference's, including the header
 * rule, the blank lines and the one-based numbering of the reasoning list: a
 * consumer that parses this text is parsing a contract.
 */
template <class T>
std::string export_analysis(const WorkflowAnalysisResult<T>& a,
                            WfExportFormat format = WfExportFormat::Summary) {
    const DetectedPatterns<T>& p = a.patternAnalysis.detectedPatterns;
    const std::string rec = wf_solver_name(a.solverRecommendation.recommendedSolver);
    std::ostringstream o;
    if (format == WfExportFormat::Json) {
        o << "{\n";
        o << "  \"solver_recommendation\": \"" << rec << "\",\n";
        o << "  \"confidence\": " << detail::plain(a.solverRecommendation.confidence) << ",\n";
        o << "  \"patterns\": {\n";
        o << "    \"sequences\": " << p.sequences.size() << ",\n";
        o << "    \"parallels\": " << p.parallels.size() << ",\n";
        o << "    \"loops\": " << p.loops.size() << ",\n";
        o << "    \"branches\": " << p.branches.size() << "\n";
        o << "  }\n";
        o << "}";
    } else if (format == WfExportFormat::Csv) {
        o << "Metric,Value\n";
        o << "Recommended Solver," << rec << "\n";
        o << "Confidence," << detail::plain(a.solverRecommendation.confidence) << "\n";
        o << "Sequences," << p.sequences.size() << "\n";
        o << "Parallels," << p.parallels.size() << "\n";
        o << "Loops," << p.loops.size() << "\n";
        o << "Branches," << p.branches.size() << "\n";
    } else {
        o << "=== Workflow Analysis Summary ===\n\n";
        o << "Recommended Solver: " << rec << "\n";
        o << "Confidence: " << detail::fixed2(a.solverRecommendation.confidence) << "\n\n";
        o << "Detected Patterns:\n";
        o << "- Sequences: " << p.sequences.size() << "\n";
        o << "- Parallels: " << p.parallels.size() << "\n";
        o << "- Loops: " << p.loops.size() << "\n";
        o << "- Branches: " << p.branches.size() << "\n\n";
        o << "Reasoning:\n";
        for (std::size_t i = 0; i < a.solverRecommendation.reasoning.size(); ++i)
            o << (i + 1) << ". " << a.solverRecommendation.reasoning[i] << "\n";
    }
    return o.str();
}

/** The reference's one-call summary. */
template <class T>
std::string quick_analysis(const qn::NetworkStruct<T>& sn) {
    return export_analysis(analyze_workflow_full(sn), WfExportFormat::Summary);
}

/** The chosen solver without the rest of the report. */
template <class T>
WfSolver get_optimal_solver(const qn::NetworkStruct<T>& sn) {
    return create_optimal_solver(sn);
}

}  // namespace wf
}  // namespace line

#endif  // LINE_API_WF_WORKFLOW_MANAGER_H
