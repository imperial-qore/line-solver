/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_WF_WF_LINK_MATRIX_H
#define LINE_API_WF_WF_LINK_MATRIX_H

/**
 * Shared conventions of the workflow pattern detectors.
 *
 * The four detectors ported from jar/src/main/java/jline/api/wf/ all read the
 * same "link matrix": one row per directed edge, column 0 the source node id,
 * column 1 the target node id, column 2 the routing probability of the edge.
 * Node ids are integers carried in the numeric type, exactly as the Java
 * Matrix does, so this header only holds the accessor that turns column 0 or 1
 * back into an id (truncating, like the Java (int) cast) and the adjacency
 * builders that three of the four detectors would otherwise duplicate.
 *
 * DIVERGENCE, deliberate: the Java detectors key their adjacency on HashMap and
 * collect results in HashSet, so the ORDER of the detected patterns, and in
 * findCommonJoinPoint even WHICH join point is returned, depend on Java's hash
 * iteration order. That is not a property of the workflow. The port uses
 * ordered containers throughout, so every result is in ascending node order and
 * is reproducible; where the Java picks an arbitrary element of a set the port
 * picks the smallest. Sets of size one - the only case the callers validate -
 * are unaffected.
 */

#include <cstddef>
#include <map>
#include <set>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace wf {

namespace detail {

/** Node id in column j of edge i, truncated as the Java (int) cast does. */
template <class T>
int wf_id(const Matrix<T>& linkMatrix, std::size_t i, std::size_t j) {
    return static_cast<int>(num_traits<T>::to_double(linkMatrix(i, j)));
}

template <class T>
void wf_check(const Matrix<T>& linkMatrix) {
    if (linkMatrix.cols() < 3)
        throw InputError("wf: the link matrix needs three columns (from, to, probability)");
}

/** node -> successors, in ascending order. */
template <class T>
std::map<int, std::vector<int>> wf_adjacency(const Matrix<T>& linkMatrix) {
    wf_check(linkMatrix);
    std::map<int, std::vector<int>> adj;
    for (std::size_t i = 0; i < linkMatrix.rows(); ++i)
        adj[wf_id(linkMatrix, i, 0)].push_back(wf_id(linkMatrix, i, 1));
    return adj;
}

/** node -> (successor, probability), in edge order. */
template <class T>
std::map<int, std::vector<std::pair<int, T>>> wf_adjacency_prob(const Matrix<T>& linkMatrix) {
    wf_check(linkMatrix);
    std::map<int, std::vector<std::pair<int, T>>> adj;
    for (std::size_t i = 0; i < linkMatrix.rows(); ++i)
        adj[wf_id(linkMatrix, i, 0)].push_back(
            std::make_pair(wf_id(linkMatrix, i, 1), linkMatrix(i, 2)));
    return adj;
}

/** node -> predecessors, in edge order. */
template <class T>
std::map<int, std::vector<int>> wf_reverse_adjacency(const Matrix<T>& linkMatrix) {
    wf_check(linkMatrix);
    std::map<int, std::vector<int>> radj;
    for (std::size_t i = 0; i < linkMatrix.rows(); ++i)
        radj[wf_id(linkMatrix, i, 1)].push_back(wf_id(linkMatrix, i, 0));
    return radj;
}

}  // namespace detail

}  // namespace wf
}  // namespace line

#endif  // LINE_API_WF_WF_LINK_MATRIX_H
