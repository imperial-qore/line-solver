/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_WF_WF_SEQUENCE_DETECTOR_H
#define LINE_API_WF_WF_SEQUENCE_DETECTOR_H

/**
 * Sequence pattern detection in a workflow network.
 *
 * Templated port of jar/src/main/java/jline/api/wf/Wf_sequence_detector.java
 * (no MATLAB counterpart exists, so the JAR is the reference). A sequence is a
 * maximal chain of service nodes connected service-to-service; the chain is
 * grown from an unused edge in both directions until no edge extends it, and
 * the number of chains looked for is half the number of service nodes that
 * appear exactly once among the service-to-service edges, i.e. half the number
 * of chain endpoints.
 *
 * detect_sequences, validate_sequence and get_sequence_stats are the three
 * public methods of the Java class, kept 1:1.
 *
 * REFERENCE DEFECT (JAR), NOT reproduced: validateSequence builds a
 * HashSet<Pair<Integer,Integer>> of the edges and asks whether each consecutive
 * pair of the chain is in it, but jline.util.Pair implements neither equals nor
 * hashCode, so contains() falls back to reference identity and NEVER matches.
 * The Java method therefore returns false for every sequence of length >= 2,
 * including the chains its own detectSequences just produced (verified against
 * common/jline.jar: validateSequence([2,3,4]) on 1->2->3->4->9 returns false).
 * Fixing Pair is outside this port, so validate_sequence here implements the
 * intended semantics - std::pair has value equality - and returns true for a
 * genuinely connected chain. This is the one place where the port deliberately
 * disagrees with the reference.
 *
 * Everything here is graph traversal and counting; the only arithmetic is the
 * mean chain length, one division. Finite field computation, instantiates at
 * exact arithmetic, no transcendental gate: for T = Rational the mean length
 * is an exact rational and the structural counts are exact integers.
 */

#include <algorithm>
#include <cstddef>
#include <map>
#include <set>
#include <vector>

#include "line/api/wf/wf_link_matrix.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace wf {

/** Mirrors the Java getSequenceStats map. */
template <class T>
struct SequenceStats {
    std::size_t numSequences = 0;
    std::size_t totalNodes = 0;
    T avgLength = num_traits<T>::from_int(0);
    std::size_t maxLength = 0;
    std::size_t minLength = 0;
};

namespace detail {

/**
 * Grow one chain out of the first remaining edge and delete the edges it used,
 * as Java's buildSequenceChain does (including the fact that the search for an
 * extension starts at index 1, edge 0 being the seed).
 */
inline std::vector<int> wf_build_sequence_chain(std::vector<std::pair<int, int>>& connections) {
    std::vector<int> sequence;
    if (connections.empty()) return sequence;

    std::vector<std::size_t> used;
    int first = connections[0].first;
    int last = connections[0].second;
    sequence.push_back(first);
    sequence.push_back(last);
    used.push_back(0);

    bool foundExtension = true;
    while (foundExtension) {
        foundExtension = false;
        const std::size_t currentSize = sequence.size();
        for (std::size_t i = 1; i < connections.size(); ++i) {
            if (std::find(used.begin(), used.end(), i) != used.end()) continue;
            const int start = connections[i].first;
            const int end = connections[i].second;
            if (start == last) {
                last = end;
                sequence.push_back(end);
                used.push_back(i);
                foundExtension = true;
            } else if (end == first) {
                first = start;
                sequence.insert(sequence.begin(), start);
                used.push_back(i);
                foundExtension = true;
            }
        }
        foundExtension = foundExtension && sequence.size() > currentSize;
    }

    std::sort(used.begin(), used.end(), std::greater<std::size_t>());
    for (std::size_t idx : used) connections.erase(connections.begin() + static_cast<long>(idx));
    return sequence;
}

}  // namespace detail

/**
 * @param linkMatrix   (nedges x 3) edge list
 * @param serviceNodes ids of the service nodes
 * @return the detected chains, each as an ordered list of node ids
 */
template <class T>
std::vector<std::vector<int>> detect_sequences(const Matrix<T>& linkMatrix,
                                               const std::vector<int>& serviceNodes) {
    detail::wf_check(linkMatrix);
    std::vector<std::vector<int>> chains;
    const std::set<int> serviceSet(serviceNodes.begin(), serviceNodes.end());

    std::vector<std::pair<int, int>> connections;
    for (std::size_t i = 0; i < linkMatrix.rows(); ++i) {
        const int s = detail::wf_id(linkMatrix, i, 0);
        const int e = detail::wf_id(linkMatrix, i, 1);
        if (serviceSet.count(s) && serviceSet.count(e))
            connections.push_back(std::make_pair(s, e));
    }
    if (connections.empty()) return chains;

    std::map<int, std::size_t> counts;
    for (std::size_t i = 0; i < connections.size(); ++i) {
        counts[connections[i].first] += 1;
        counts[connections[i].second] += 1;
    }
    std::size_t countOnce = 0;
    for (std::map<int, std::size_t>::const_iterator it = counts.begin(); it != counts.end(); ++it)
        if (it->second == 1) ++countOnce;
    const std::size_t numSequences = countOnce / 2;

    for (std::size_t seq = 0; seq < numSequences; ++seq) {
        if (connections.empty()) break;
        std::vector<int> chain = detail::wf_build_sequence_chain(connections);
        if (!chain.empty()) chains.push_back(chain);
    }
    return chains;
}

/** Every consecutive pair of the chain must be an edge of the workflow. */
template <class T>
bool validate_sequence(const std::vector<int>& sequence, const Matrix<T>& linkMatrix) {
    detail::wf_check(linkMatrix);
    if (sequence.size() < 2) return false;
    std::set<std::pair<int, int>> edges;
    for (std::size_t i = 0; i < linkMatrix.rows(); ++i)
        edges.insert(std::make_pair(detail::wf_id(linkMatrix, i, 0),
                                    detail::wf_id(linkMatrix, i, 1)));
    for (std::size_t i = 0; i + 1 < sequence.size(); ++i)
        if (!edges.count(std::make_pair(sequence[i], sequence[i + 1]))) return false;
    return true;
}

/** Count, total, mean, maximum and minimum chain length. */
template <class T>
SequenceStats<T> get_sequence_stats(const std::vector<std::vector<int>>& sequences) {
    SequenceStats<T> stats;
    stats.numSequences = sequences.size();
    std::size_t total = 0;
    for (std::size_t i = 0; i < sequences.size(); ++i) total += sequences[i].size();
    stats.totalNodes = total;
    if (sequences.empty()) return stats;

    stats.avgLength = num_traits<T>::from_int(static_cast<long>(total)) /
                      num_traits<T>::from_int(static_cast<long>(sequences.size()));
    stats.maxLength = sequences[0].size();
    stats.minLength = sequences[0].size();
    for (std::size_t i = 1; i < sequences.size(); ++i) {
        if (sequences[i].size() > stats.maxLength) stats.maxLength = sequences[i].size();
        if (sequences[i].size() < stats.minLength) stats.minLength = sequences[i].size();
    }
    return stats;
}

}  // namespace wf
}  // namespace line

#endif  // LINE_API_WF_WF_SEQUENCE_DETECTOR_H
