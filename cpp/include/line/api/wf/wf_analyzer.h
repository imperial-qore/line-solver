/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_WF_WF_ANALYZER_H
#define LINE_API_WF_WF_ANALYZER_H

/**
 * The workflow analyzer: detect every pattern, collapse them, report the two
 * complexities and the recommendations that follow.
 *
 * Templated port of the native Python `line_solver/api/wf/analyzer.py`,
 * cross-checked against jar/src/main/java/jline/api/wf/Wf_analyzer.java. There
 * is no MATLAB counterpart; `api/wf` exists in the JAR and in Python only, and
 * Python is the reference for the same reason it is in `wf_pattern_updater.h`
 * (the JAR's convolutions are stubs, so its "optimized" workflow carries the
 * first branch's service law).
 *
 * THE NETWORK CONVERSION IS THE ONE PART THAT CANNOT BE COPIED. Both references
 * walk their own object model -- Python asks `type(node).__name__` and reads
 * `getLinkedRoutingMatrix`, the JAR walks `jline.lang.Network` -- and neither
 * shape exists here. `wf_from_struct` does the same job from a `NetworkStruct`:
 * the link matrix is the class-aggregated `rtnodes` above the zero tolerance,
 * and the node classification follows `NodeType`, Queue and Delay being service,
 * Fork and Join their own kinds, Router and ClassSwitch control. The service
 * laws come from `sn.service` rather than from Python's placeholder unit
 * exponential, which is strictly more information and changes no structure.
 *
 * `analyze_workflow` on a representation the caller built by hand is the entry
 * both references really exercise, and it is byte-for-byte their algorithm.
 *
 * ARITHMETIC: field, plus whatever the branch entropy needs -- the diversity
 * report is gated inside wf_branch_detector.h, not here.
 */

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "line/api/wf/wf_branch_detector.h"
#include "line/api/wf/wf_link_matrix.h"
#include "line/api/wf/wf_loop_detector.h"
#include "line/api/wf/wf_parallel_detector.h"
#include "line/api/wf/wf_pattern_updater.h"
#include "line/api/wf/wf_sequence_detector.h"
#include "line/lang/distribution.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace wf {

/** The workflow in matrix form: the reference's WorkflowRepresentation. */
template <class T>
struct WorkflowRepresentation {
    Matrix<T> linkMatrix;
    std::vector<int> serviceNodes, forkNodes, joinNodes, routerNodes;
    std::map<int, ServiceParameters<T>> serviceParameters;
};

/** Everything the four detectors found. */
template <class T>
struct DetectedPatterns {
    std::vector<std::vector<int>> sequences;
    std::vector<std::vector<int>> parallels;
    std::vector<int> loops;
    std::vector<BranchPattern<T>> branches;
};

/** The reference's complexity map, for either the original or the collapsed graph. */
template <class T>
struct WorkflowComplexity {
    std::size_t totalNodes = 0;
    std::size_t totalLinks = 0;
    std::size_t serviceNodes = 0;
    std::size_t controlNodes = 0;
    std::size_t connectedNodes = 0;
    T avgDegree = num_traits<T>::from_int(0);
};

/** The statistics block of WorkflowAnalysis. */
template <class T>
struct WorkflowStatistics {
    SequenceStats<T> sequenceStats;
    ParallelStats<T> parallelStats;
    LoopStats<T> loopStats;
    BranchStats<T> branchStats;
    UpdateStats<T> updateStats;
    WorkflowComplexity<T> originalComplexity;
    WorkflowComplexity<T> optimizedComplexity;
};

/** What analyze_workflow returns. */
template <class T>
struct WorkflowAnalysis {
    WorkflowRepresentation<T> originalWorkflow;
    DetectedPatterns<T> detectedPatterns;
    UpdatedWorkflow<T> optimizedWorkflow;
    WorkflowStatistics<T> statistics;
};

/** Run the four detectors on one representation. */
template <class T>
DetectedPatterns<T> detect_all_patterns(const WorkflowRepresentation<T>& w) {
    DetectedPatterns<T> p;
    p.sequences = detect_sequences(w.linkMatrix, w.serviceNodes);
    p.parallels = detect_parallel(w.linkMatrix, w.serviceNodes, w.forkNodes, w.joinNodes);
    p.loops = detect_loops(w.linkMatrix, w.serviceNodes, w.routerNodes, w.joinNodes);
    p.branches = detect_branches(w.linkMatrix, w.serviceNodes, w.joinNodes);
    return p;
}

namespace detail {

/** Distinct node ids the link matrix mentions, and the total degree of each. */
template <class T>
void wf_degrees(const Matrix<T>& m, std::set<int>* nodes, std::map<int, std::size_t>* deg) {
    for (std::size_t i = 0; i < m.rows(); ++i) {
        const int a = wf_id(m, i, 0), b = wf_id(m, i, 1);
        nodes->insert(a);
        nodes->insert(b);
        if (deg != 0) {
            (*deg)[a] += 1;
            (*deg)[b] += 1;
        }
    }
}

}  // namespace detail

/** Complexity of the workflow as declared. */
template <class T>
WorkflowComplexity<T> workflow_complexity(const WorkflowRepresentation<T>& w) {
    WorkflowComplexity<T> c;
    c.serviceNodes = w.serviceNodes.size();
    c.controlNodes = w.forkNodes.size() + w.joinNodes.size() + w.routerNodes.size();
    c.totalNodes = c.serviceNodes + c.controlNodes;
    c.totalLinks = w.linkMatrix.rows();

    std::set<int> nodes;
    std::map<int, std::size_t> deg;
    detail::wf_degrees(w.linkMatrix, &nodes, &deg);
    c.connectedNodes = nodes.size();
    if (!deg.empty()) {
        T s = num_traits<T>::from_int(0);
        for (std::map<int, std::size_t>::const_iterator it = deg.begin(); it != deg.end(); ++it)
            s += num_traits<T>::from_int(static_cast<long>(it->second));
        c.avgDegree = T(s / num_traits<T>::from_int(static_cast<long>(deg.size())));
    }
    return c;
}

/**
 * Complexity of the collapsed workflow.
 *
 * The reference reports only three of the six fields here, so the rest stay at
 * their defaults rather than being invented: after the collapse there is no
 * service/control split left to report, the node kinds having been merged.
 */
template <class T>
WorkflowComplexity<T> optimized_complexity(const UpdatedWorkflow<T>& w) {
    WorkflowComplexity<T> c;
    c.totalNodes = w.serviceParameters.size();
    c.totalLinks = w.linkMatrix.rows();
    std::set<int> nodes;
    detail::wf_degrees(w.linkMatrix, &nodes, static_cast<std::map<int, std::size_t>*>(0));
    c.connectedNodes = nodes.size();
    return c;
}

/** Detect, collapse, and report. */
template <class T>
WorkflowAnalysis<T> analyze_workflow(const WorkflowRepresentation<T>& w) {
    WorkflowAnalysis<T> a;
    a.originalWorkflow = w;
    a.detectedPatterns = detect_all_patterns(w);
    a.optimizedWorkflow = update_patterns(w.linkMatrix, w.serviceNodes, w.forkNodes, w.joinNodes,
                                          w.routerNodes, w.serviceParameters);

    a.statistics.sequenceStats = get_sequence_stats<T>(a.detectedPatterns.sequences);
    a.statistics.parallelStats = get_parallel_stats<T>(a.detectedPatterns.parallels);
    a.statistics.loopStats = get_loop_stats(a.detectedPatterns.loops, w.linkMatrix, w.routerNodes);
    a.statistics.branchStats = get_branch_stats(a.detectedPatterns.branches);
    a.statistics.updateStats = get_update_stats(w.linkMatrix, a.optimizedWorkflow);
    a.statistics.originalComplexity = workflow_complexity(w);
    a.statistics.optimizedComplexity = optimized_complexity(a.optimizedWorkflow);
    return a;
}

/**
 * The reference's recommendation strings, in its order.
 *
 * They are user-facing text, so they are reproduced verbatim rather than
 * paraphrased: a caller that greps them would otherwise see different output
 * from the same analysis in two codebases.
 */
template <class T>
std::vector<std::string> get_optimization_recommendations(const WorkflowAnalysis<T>& a) {
    std::vector<std::string> out;
    const DetectedPatterns<T>& p = a.detectedPatterns;
    if (!p.sequences.empty())
        out.push_back("Found " + std::to_string(p.sequences.size()) +
                      " sequence patterns that can be simplified");
    if (!p.parallels.empty())
        out.push_back("Found " + std::to_string(p.parallels.size()) +
                      " parallel patterns for potential optimization");
    if (!p.loops.empty())
        out.push_back("Found " + std::to_string(p.loops.size()) +
                      " loop patterns - consider loop unrolling for performance");
    if (!p.branches.empty()) {
        out.push_back("Found " + std::to_string(p.branches.size()) +
                      " branch patterns - analyze probability distributions");
        std::size_t high = 0;
        for (std::size_t i = 0; i < p.branches.size(); ++i)
            if (num_traits<T>::to_double(calculate_branch_diversity(p.branches[i]).entropy) > 1.0)
                ++high;
        if (high > 0)
            out.push_back(std::to_string(high) +
                          " branches have high entropy - consider load balancing");
    }
    const double ratio = num_traits<T>::to_double(a.statistics.updateStats.reductionRatio);
    if (ratio > 0.1)
        out.push_back("Workflow complexity reduced by " +
                      std::to_string(static_cast<int>(ratio * 100)) +
                      "% through pattern optimization");
    return out;
}

/** The collapsed workflow is consistent and every detected pattern validates. */
template <class T>
bool validate_analysis(const WorkflowAnalysis<T>& a) {
    if (!validate_updated_workflow(a.optimizedWorkflow)) return false;
    const Matrix<T>& L = a.originalWorkflow.linkMatrix;
    for (std::size_t i = 0; i < a.detectedPatterns.sequences.size(); ++i)
        if (!validate_sequence(a.detectedPatterns.sequences[i], L)) return false;
    for (std::size_t i = 0; i < a.detectedPatterns.parallels.size(); ++i)
        if (!validate_parallel_pattern(a.detectedPatterns.parallels[i], L,
                                       a.originalWorkflow.forkNodes,
                                       a.originalWorkflow.joinNodes))
            return false;
    for (std::size_t i = 0; i < a.detectedPatterns.loops.size(); ++i)
        if (!validate_loop_pattern(a.detectedPatterns.loops[i], L, a.originalWorkflow.routerNodes))
            return false;
    for (std::size_t i = 0; i < a.detectedPatterns.branches.size(); ++i)
        if (!validate_branch_pattern(a.detectedPatterns.branches[i], L)) return false;
    return true;
}

/**
 * Build a workflow representation from a NetworkStruct.
 *
 * The link matrix is the class-aggregated node routing `sn.rtnodes`, one row per
 * (source, target) pair carrying positive probability, with the probability
 * summed over the class pairs and normalized per source. Node ids are 0-based,
 * matching the reference's node indices.
 */
template <class T>
WorkflowRepresentation<T> wf_from_struct(const qn::NetworkStruct<T>& sn) {
    const T zero = num_traits<T>::from_int(0);
    const double tol = lang::GlobalConstants::Zero;
    const std::size_t I = sn.nodes.size(), K = sn.nclasses;

    WorkflowRepresentation<T> w;
    for (std::size_t ind = 0; ind < I; ++ind) {
        switch (sn.nodes[ind].nodetype) {
            case qn::NodeType::Queue:
            case qn::NodeType::Delay:
                w.serviceNodes.push_back(static_cast<int>(ind));
                break;
            case qn::NodeType::Fork:
                w.forkNodes.push_back(static_cast<int>(ind));
                break;
            case qn::NodeType::Join:
                w.joinNodes.push_back(static_cast<int>(ind));
                break;
            case qn::NodeType::Router:
            case qn::NodeType::ClassSwitch:
                w.routerNodes.push_back(static_cast<int>(ind));
                break;
            default:
                break;
        }
    }

    std::vector<std::vector<T>> agg(I, std::vector<T>(I, zero));
    if (sn.rtnodes.rows() == I * K && sn.rtnodes.cols() == I * K) {
        for (std::size_t i = 0; i < I; ++i)
            for (std::size_t j = 0; j < I; ++j)
                for (std::size_t r = 0; r < K; ++r)
                    for (std::size_t s = 0; s < K; ++s)
                        agg[i][j] += sn.rtnodes(i * K + r, j * K + s);
    }
    std::vector<std::pair<std::pair<int, int>, T>> links;
    for (std::size_t i = 0; i < I; ++i) {
        T rowsum = zero;
        for (std::size_t j = 0; j < I; ++j) rowsum += agg[i][j];
        if (!(num_traits<T>::to_double(rowsum) > tol)) continue;
        for (std::size_t j = 0; j < I; ++j) {
            if (!(num_traits<T>::to_double(agg[i][j]) > tol)) continue;
            links.push_back(std::make_pair(std::make_pair(static_cast<int>(i), static_cast<int>(j)),
                                           T(agg[i][j] / rowsum)));
        }
    }
    w.linkMatrix = Matrix<T>(links.size(), 3, zero);
    for (std::size_t k = 0; k < links.size(); ++k) {
        w.linkMatrix(k, 0) = num_traits<T>::from_int(links[k].first.first);
        w.linkMatrix(k, 1) = num_traits<T>::from_int(links[k].first.second);
        w.linkMatrix(k, 2) = links[k].second;
    }

    // Service laws from sn.service, one per service node, taking the first
    // enabled class. The references install a placeholder unit exponential
    // here; reading the declared law instead is strictly more information and
    // leaves the structure untouched.
    for (std::size_t k = 0; k < w.serviceNodes.size(); ++k) {
        const std::size_t ind = static_cast<std::size_t>(w.serviceNodes[k]);
        const std::size_t ist = sn.nodes[ind].station;
        ServiceParameters<T> p;
        p.alpha.assign(1, num_traits<T>::from_int(1));
        p.T_ = Matrix<T>(1, 1, num_traits<T>::from_int(-1));
        if (ist >= 1 && ist <= sn.nstations) {
            for (std::size_t r = 0; r < K; ++r) {
                if (sn.disabled[ist - 1][r] || sn.service[ist - 1][r].disabled) continue;
                const mam::Map<T> m = lang::dist_to_map(sn.service[ist - 1][r]);
                p.alpha = mam::map_pie(m);
                p.T_ = m.D0;
                break;
            }
        }
        w.serviceParameters[w.serviceNodes[k]] = p;
    }
    return w;
}

}  // namespace wf
}  // namespace line

#endif  // LINE_API_WF_WF_ANALYZER_H
