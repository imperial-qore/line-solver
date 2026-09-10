/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_WF_WF_PATTERN_UPDATER_H
#define LINE_API_WF_WF_PATTERN_UPDATER_H

/**
 * Collapse the detected workflow patterns and convolve their service laws.
 *
 * Templated port of the native Python `line_solver/api/wf/pattern_updater.py`,
 * cross-checked against jar/src/main/java/jline/api/wf/Wf_pattern_updater.java.
 * There is no MATLAB counterpart: `api/wf` exists in the JAR and in Python only.
 *
 * PYTHON IS THE REFERENCE HERE, NOT THE JAR. The JAR's four convolutions are
 * STUBS -- `convolveSequence`, `convolveParallel` and `convolveBranches` return
 * `params.get(0)` unchanged and `convolveLoop` returns its argument -- so the
 * Java class rewrites the link matrix and then reports the FIRST branch's
 * service law as the law of the collapsed pattern. Its `removeMatrixRows` is
 * also defective: the loop over the sorted row list `return`s inside its first
 * iteration, so at most one row is ever removed. Python implements the actual
 * phase-type algebra and removes every row, and that is what is ported.
 *
 * THE FOUR CONVOLUTIONS, on representations (alpha, T) that need not be
 * honest phase types -- alpha may sum below one, and the deficit 1 - alpha e is
 * treated as instantaneous completion, which is what makes the formulas below
 * carry the sub-stochastic entry vectors the detectors produce:
 *
 *  - SEQUENCE, the convolution of the two durations:
 *      alpha = [a1, (1 - a1 e1) a2],  T = [[T1, (-T1 e1) a2], [0, T2]].
 *  - PARALLEL, the MAXIMUM of the two durations, so the state is the pair of
 *    phases until one branch finishes and the surviving branch alone after:
 *      alpha = [kron(a1,a2), (1 - a2 e2) a1, (1 - a1 e1) a2],
 *      T     = [[T1 (x) I + I (x) T2, I (x) (-T2 e2), (-T1 e1) (x) I],
 *               [0, T1, 0], [0, 0, T2]].
 *  - LOOP, a geometric number of repetitions with probability p: the exit flow
 *    is fed back into the entry law, T <- T + p (-T e) alpha, which leaves
 *    alpha unchanged. p outside (0,1) is the identity.
 *  - BRANCH, a probabilistic choice: alpha = [p1 a1, p2 a2, ...] over a
 *    block-diagonal T, with the probabilities renormalized (and made uniform
 *    when they sum to zero).
 *
 * `find_fork_join_for_parallel` returns "none" in BOTH references, so the
 * parallel arm of `update_patterns` never fires. That is reproduced rather than
 * invented: supplying a fork/join search here would collapse patterns neither
 * reference collapses, and the resulting workflow would not be the one any
 * other codebase produces. The convolution itself is implemented and reachable
 * through `convolve_parallel`, which is what a caller with its own fork/join
 * pairing needs.
 *
 * ARITHMETIC: field. Block assembly, Kronecker products, and one division in
 * the probability renormalization, so it instantiates under Rational.
 */

#include <algorithm>
#include <cstddef>
#include <map>
#include <set>
#include <vector>

#include "line/api/wf/wf_branch_detector.h"
#include "line/api/wf/wf_link_matrix.h"
#include "line/api/wf/wf_loop_detector.h"
#include "line/api/wf/wf_parallel_detector.h"
#include "line/api/wf/wf_sequence_detector.h"
#include "line/lang/lang_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace wf {

/** A phase-type-shaped service law: entry vector and transient generator. */
template <class T>
struct ServiceParameters {
    std::vector<T> alpha;
    Matrix<T> T_;
};

/** The collapsed link matrix and the service law of every surviving node. */
template <class T>
struct UpdatedWorkflow {
    Matrix<T> linkMatrix;
    std::map<int, ServiceParameters<T>> serviceParameters;
};

/** Statistics of one update pass; the Java getUpdateStats map. */
template <class T>
struct UpdateStats {
    std::size_t originalLinks = 0;
    std::size_t updatedLinks = 0;
    long linksReduced = 0;
    std::size_t serviceNodes = 0;
    T reductionRatio = num_traits<T>::from_int(0);
};

namespace detail {

/** The unit-mass fallback both references return for an empty parameter list. */
template <class T>
ServiceParameters<T> wf_unit_params() {
    ServiceParameters<T> p;
    p.alpha.assign(1, num_traits<T>::from_int(1));
    p.T_ = Matrix<T>(1, 1, num_traits<T>::from_int(0));
    return p;
}

/** -T e, the exit rate out of each phase. */
template <class T>
std::vector<T> wf_exit_rate(const Matrix<T>& Tm) {
    const T zero = num_traits<T>::from_int(0);
    std::vector<T> v(Tm.rows(), zero);
    for (std::size_t i = 0; i < Tm.rows(); ++i) {
        T s = zero;
        for (std::size_t j = 0; j < Tm.cols(); ++j) s += Tm(i, j);
        v[i] = -s;
    }
    return v;
}

/** 1 - alpha e, the mass the entry law leaves for instantaneous completion. */
template <class T>
T wf_exit_prob(const std::vector<T>& alpha) {
    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < alpha.size(); ++i) s += alpha[i];
    return T(num_traits<T>::from_int(1) - s);
}

/** Copy `src` into `dst` with its top-left corner at (r0, c0). */
template <class T>
void wf_place(Matrix<T>& dst, const Matrix<T>& src, std::size_t r0, std::size_t c0) {
    for (std::size_t i = 0; i < src.rows(); ++i)
        for (std::size_t j = 0; j < src.cols(); ++j) dst(r0 + i, c0 + j) = src(i, j);
}

/** Drop the listed rows; duplicates and out-of-range indices are ignored. */
template <class T>
Matrix<T> wf_remove_rows(const Matrix<T>& m, const std::vector<std::size_t>& rows) {
    if (rows.empty() || m.rows() == 0) return m;
    std::vector<bool> keep(m.rows(), true);
    for (std::size_t i = 0; i < rows.size(); ++i)
        if (rows[i] < m.rows()) keep[rows[i]] = false;
    std::size_t n = 0;
    for (std::size_t i = 0; i < m.rows(); ++i)
        if (keep[i]) ++n;
    Matrix<T> out(n, m.cols(), num_traits<T>::from_int(0));
    std::size_t r = 0;
    for (std::size_t i = 0; i < m.rows(); ++i) {
        if (!keep[i]) continue;
        for (std::size_t j = 0; j < m.cols(); ++j) out(r, j) = m(i, j);
        ++r;
    }
    return out;
}

/** Rewrite every reference to `oldNode` in columns 0 and 1 as `newNode`. */
template <class T>
void wf_replace_node(Matrix<T>& m, int oldNode, int newNode) {
    const T nn = num_traits<T>::from_int(newNode);
    for (std::size_t i = 0; i < m.rows(); ++i) {
        if (wf_id(m, i, 0) == oldNode) m(i, 0) = nn;
        if (wf_id(m, i, 1) == oldNode) m(i, 1) = nn;
    }
}

/** Rows whose source or target is one of `nodes`. */
template <class T>
std::vector<std::size_t> wf_rows_involving(const Matrix<T>& m, const std::vector<int>& nodes) {
    const std::set<int> s(nodes.begin(), nodes.end());
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < m.rows(); ++i)
        if (s.count(wf_id(m, i, 0)) || s.count(wf_id(m, i, 1))) out.push_back(i);
    return out;
}

}  // namespace detail

/** Convolution of the durations, i.e. the service laws run one after another. */
template <class T>
ServiceParameters<T> convolve_sequence(const std::vector<ServiceParameters<T>>& params) {
    if (params.empty()) return detail::wf_unit_params<T>();
    if (params.size() == 1) return params[0];
    const T zero = num_traits<T>::from_int(0);

    ServiceParameters<T> acc = params[0];
    for (std::size_t k = 1; k < params.size(); ++k) {
        const std::vector<T>& a2 = params[k].alpha;
        const Matrix<T>& T2 = params[k].T_;
        const std::size_t n1 = acc.alpha.size(), n2 = a2.size();
        const T ex = detail::wf_exit_prob(acc.alpha);
        const std::vector<T> er = detail::wf_exit_rate(acc.T_);

        std::vector<T> na(n1 + n2, zero);
        for (std::size_t i = 0; i < n1; ++i) na[i] = acc.alpha[i];
        for (std::size_t j = 0; j < n2; ++j) na[n1 + j] = T(ex * a2[j]);

        Matrix<T> nt(n1 + n2, n1 + n2, zero);
        detail::wf_place(nt, acc.T_, 0, 0);
        detail::wf_place(nt, T2, n1, n1);
        for (std::size_t i = 0; i < n1; ++i)
            for (std::size_t j = 0; j < n2; ++j) nt(i, n1 + j) = T(er[i] * a2[j]);

        acc.alpha = na;
        acc.T_ = nt;
    }
    return acc;
}

/** Maximum of the durations, i.e. a fork whose join waits for every branch. */
template <class T>
ServiceParameters<T> convolve_parallel(const std::vector<ServiceParameters<T>>& params) {
    if (params.empty()) return detail::wf_unit_params<T>();
    if (params.size() == 1) return params[0];
    const T zero = num_traits<T>::from_int(0);

    ServiceParameters<T> acc = params[0];
    for (std::size_t k = 1; k < params.size(); ++k) {
        const std::vector<T>& a2 = params[k].alpha;
        const Matrix<T>& T2 = params[k].T_;
        const std::size_t n1 = acc.alpha.size(), n2 = a2.size(), np = n1 * n2;
        const T e1 = detail::wf_exit_prob(acc.alpha), e2 = detail::wf_exit_prob(a2);
        const std::vector<T> r1 = detail::wf_exit_rate(acc.T_), r2 = detail::wf_exit_rate(T2);

        std::vector<T> na(np + n1 + n2, zero);
        for (std::size_t i = 0; i < n1; ++i)
            for (std::size_t j = 0; j < n2; ++j) na[i * n2 + j] = T(acc.alpha[i] * a2[j]);
        for (std::size_t i = 0; i < n1; ++i) na[np + i] = T(e2 * acc.alpha[i]);
        for (std::size_t j = 0; j < n2; ++j) na[np + n1 + j] = T(e1 * a2[j]);

        Matrix<T> nt(np + n1 + n2, np + n1 + n2, zero);
        // T1 (x) I + I (x) T2 on the both-alive block.
        for (std::size_t i = 0; i < n1; ++i)
            for (std::size_t j = 0; j < n2; ++j) {
                const std::size_t r = i * n2 + j;
                for (std::size_t ii = 0; ii < n1; ++ii) nt(r, ii * n2 + j) += acc.T_(i, ii);
                for (std::size_t jj = 0; jj < n2; ++jj) nt(r, i * n2 + jj) += T2(j, jj);
            }
        // Branch 2 finishes first -> only branch 1 is left, in phase i.
        for (std::size_t i = 0; i < n1; ++i)
            for (std::size_t j = 0; j < n2; ++j) nt(i * n2 + j, np + i) = r2[j];
        // Branch 1 finishes first -> only branch 2 is left, in phase j.
        for (std::size_t i = 0; i < n1; ++i)
            for (std::size_t j = 0; j < n2; ++j) nt(i * n2 + j, np + n1 + j) = r1[i];
        detail::wf_place(nt, acc.T_, np, np);
        detail::wf_place(nt, T2, np + n1, np + n1);

        acc.alpha = na;
        acc.T_ = nt;
    }
    return acc;
}

/**
 * Geometric repetition: the exit flow re-enters through alpha with probability
 * `loopProb`. Outside (0,1) the law is returned unchanged, as in the reference.
 */
template <class T>
ServiceParameters<T> convolve_loop(const ServiceParameters<T>& params, const T& loopProb) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (!(loopProb > zero) || !(loopProb < one)) return params;
    const std::vector<T> er = detail::wf_exit_rate(params.T_);
    ServiceParameters<T> out = params;
    for (std::size_t i = 0; i < out.T_.rows(); ++i)
        for (std::size_t j = 0; j < out.T_.cols(); ++j)
            out.T_(i, j) = T(out.T_(i, j) + loopProb * er[i] * params.alpha[j]);
    return out;
}

/** Probabilistic choice among the alternatives, on a block-diagonal generator. */
template <class T>
ServiceParameters<T> convolve_branches(const std::vector<ServiceParameters<T>>& params,
                                       const std::vector<T>& probsIn) {
    if (params.empty()) return detail::wf_unit_params<T>();
    if (params.size() == 1) return params[0];
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    std::vector<T> probs(params.size(), zero);
    T total = zero;
    for (std::size_t i = 0; i < params.size(); ++i) {
        probs[i] = i < probsIn.size() ? probsIn[i] : zero;
        total += probs[i];
    }
    if (total > zero) {
        for (std::size_t i = 0; i < probs.size(); ++i) probs[i] = T(probs[i] / total);
    } else {
        const T u = T(one / num_traits<T>::from_int(static_cast<long>(params.size())));
        for (std::size_t i = 0; i < probs.size(); ++i) probs[i] = u;
    }

    std::size_t n = 0;
    for (std::size_t i = 0; i < params.size(); ++i) n += params[i].alpha.size();
    ServiceParameters<T> out;
    out.alpha.assign(n, zero);
    out.T_ = Matrix<T>(n, n, zero);
    std::size_t off = 0;
    for (std::size_t i = 0; i < params.size(); ++i) {
        for (std::size_t j = 0; j < params[i].alpha.size(); ++j)
            out.alpha[off + j] = T(probs[i] * params[i].alpha[j]);
        detail::wf_place(out.T_, params[i].T_, off, off);
        off += params[i].alpha.size();
    }
    return out;
}

/**
 * The fork and join bracketing a parallel pattern.
 *
 * BOTH references return "none" unconditionally, so the parallel arm of
 * `update_patterns` never fires. Reproduced deliberately; see the header note.
 */
template <class T>
bool find_fork_join_for_parallel(const Matrix<T>&, const std::vector<int>&, int*, int*) {
    return false;
}

/**
 * Collapse the four pattern families in the reference's order: sequences,
 * parallels, loops, branches. Each stage re-detects on the matrix the previous
 * stage produced.
 */
template <class T>
UpdatedWorkflow<T> update_patterns(const Matrix<T>& linkMatrix,
                                   const std::vector<int>& serviceNodes,
                                   const std::vector<int>& forkNodes,
                                   const std::vector<int>& joinNodes,
                                   const std::vector<int>& routerNodes,
                                   const std::map<int, ServiceParameters<T>>& serviceParams) {
    detail::wf_check(linkMatrix);
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    Matrix<T> m = linkMatrix;
    std::map<int, ServiceParameters<T>> params = serviceParams;
    const std::set<int> serviceSet(serviceNodes.begin(), serviceNodes.end());

    // ---- sequences ------------------------------------------------------
    const std::vector<std::vector<int>> sequences = detect_sequences(m, serviceNodes);
    if (!sequences.empty()) {
        std::vector<std::size_t> drop;
        for (std::size_t i = 0; i < m.rows(); ++i) {
            const int a = detail::wf_id(m, i, 0), b = detail::wf_id(m, i, 1);
            if (!serviceSet.count(a) || !serviceSet.count(b)) continue;
            for (std::size_t s = 0; s < sequences.size(); ++s)
                for (std::size_t j = 0; j + 1 < sequences[s].size(); ++j)
                    if (sequences[s][j] == a && sequences[s][j + 1] == b) {
                        drop.push_back(i);
                        j = sequences[s].size();  // break, as the reference does
                    }
        }
        m = detail::wf_remove_rows(m, drop);
        for (std::size_t s = 0; s < sequences.size(); ++s) {
            const std::vector<int>& seq = sequences[s];
            if (seq.size() < 2) continue;
            detail::wf_replace_node(m, seq.back(), seq.front());
            std::vector<ServiceParameters<T>> sp;
            for (std::size_t j = 0; j < seq.size(); ++j) {
                typename std::map<int, ServiceParameters<T>>::const_iterator it = params.find(seq[j]);
                if (it != params.end()) sp.push_back(it->second);
            }
            params[seq.front()] = convolve_sequence(sp);
            for (std::size_t j = 1; j < seq.size(); ++j) params.erase(seq[j]);
        }
    }

    // ---- parallels ------------------------------------------------------
    const std::vector<std::vector<int>> parallels =
        detect_parallel(m, serviceNodes, forkNodes, joinNodes);
    for (std::size_t p = 0; p < parallels.size(); ++p) {
        const std::vector<int>& par = parallels[p];
        if (par.size() < 2) continue;
        int forkNode = -1, joinNode = -1;
        if (!find_fork_join_for_parallel(m, par, &forkNode, &joinNode)) continue;
        m = detail::wf_remove_rows(m, detail::wf_rows_involving(m, par));
        detail::wf_replace_node(m, forkNode, par.front());
        detail::wf_replace_node(m, joinNode, par.front());
        std::vector<ServiceParameters<T>> sp;
        for (std::size_t j = 0; j < par.size(); ++j) {
            typename std::map<int, ServiceParameters<T>>::const_iterator it = params.find(par[j]);
            if (it != params.end()) sp.push_back(it->second);
        }
        params[par.front()] = convolve_parallel(sp);
        for (std::size_t j = 1; j < par.size(); ++j) params.erase(par[j]);
    }

    // ---- loops ----------------------------------------------------------
    const std::vector<int> loops = detect_loops(m, serviceNodes, routerNodes, joinNodes);
    const std::set<int> routerSet(routerNodes.begin(), routerNodes.end());
    for (std::size_t l = 0; l < loops.size(); ++l) {
        const int loopNode = loops[l];
        const T loopProb = get_loop_probability(loopNode, m, routerNodes);
        if (!(num_traits<T>::to_double(loopProb) > lang::GlobalConstants::Zero)) continue;

        std::vector<int> routers;
        for (std::size_t i = 0; i < m.rows(); ++i) {
            const int a = detail::wf_id(m, i, 0), b = detail::wf_id(m, i, 1);
            if (!((a == loopNode && routerSet.count(b)) || (b == loopNode && routerSet.count(a))))
                continue;
            if (routerSet.count(a)) routers.push_back(a);
            if (routerSet.count(b)) routers.push_back(b);
        }
        std::sort(routers.begin(), routers.end());
        routers.erase(std::unique(routers.begin(), routers.end()), routers.end());

        const std::set<int> rs(routers.begin(), routers.end());
        std::vector<std::size_t> drop;
        for (std::size_t i = 0; i < m.rows(); ++i) {
            const int a = detail::wf_id(m, i, 0), b = detail::wf_id(m, i, 1);
            if ((a == loopNode && rs.count(b)) || (b == loopNode && rs.count(a)))
                drop.push_back(i);
        }
        m = detail::wf_remove_rows(m, drop);
        for (std::size_t r = 0; r < routers.size(); ++r)
            detail::wf_replace_node(m, routers[r], loopNode);
        for (std::size_t i = 0; i < m.rows(); ++i)
            if (detail::wf_id(m, i, 0) == loopNode) m(i, 2) = one;

        typename std::map<int, ServiceParameters<T>>::iterator it = params.find(loopNode);
        if (it != params.end()) it->second = convolve_loop(it->second, loopProb);
    }

    // ---- branches -------------------------------------------------------
    const std::vector<BranchPattern<T>> branches = detect_branches(m, serviceNodes, joinNodes);
    for (std::size_t b = 0; b < branches.size(); ++b) {
        const BranchPattern<T>& br = branches[b];
        if (br.branchNodes.size() < 2 || br.forkNode < 0) continue;
        m = detail::wf_remove_rows(m, detail::wf_rows_involving(m, br.branchNodes));
        detail::wf_replace_node(m, br.forkNode, br.branchNodes.front());
        if (br.hasJoinNode) detail::wf_replace_node(m, br.joinNode, br.branchNodes.front());
        std::vector<ServiceParameters<T>> sp;
        for (std::size_t j = 0; j < br.branchNodes.size(); ++j) {
            typename std::map<int, ServiceParameters<T>>::const_iterator it =
                params.find(br.branchNodes[j]);
            if (it != params.end()) sp.push_back(it->second);
        }
        params[br.branchNodes.front()] = convolve_branches(sp, br.probabilities);
        for (std::size_t j = 1; j < br.branchNodes.size(); ++j) params.erase(br.branchNodes[j]);
    }

    UpdatedWorkflow<T> out;
    out.linkMatrix = m;
    out.serviceParameters = params;
    (void)zero;
    return out;
}

/**
 * Every node the collapsed matrix still references carries a service law.
 *
 * The two references test OPPOSITE implications: Python asks that every
 * referenced node have a law (ported here), the JAR that every law belong to a
 * referenced node. Python's is the one that catches the failure mode the
 * collapse can actually produce -- a node left in the matrix whose law was
 * erased with its pattern.
 */
template <class T>
bool validate_updated_workflow(const UpdatedWorkflow<T>& w) {
    std::set<int> referenced;
    for (std::size_t i = 0; i < w.linkMatrix.rows(); ++i) {
        referenced.insert(detail::wf_id(w.linkMatrix, i, 0));
        referenced.insert(detail::wf_id(w.linkMatrix, i, 1));
    }
    for (std::set<int>::const_iterator it = referenced.begin(); it != referenced.end(); ++it)
        if (w.serviceParameters.find(*it) == w.serviceParameters.end()) return false;
    return true;
}

/** How much the collapse shrank the link matrix. */
template <class T>
UpdateStats<T> get_update_stats(const Matrix<T>& originalMatrix, const UpdatedWorkflow<T>& w) {
    UpdateStats<T> s;
    s.originalLinks = originalMatrix.rows();
    s.updatedLinks = w.linkMatrix.rows();
    s.linksReduced = static_cast<long>(s.originalLinks) - static_cast<long>(s.updatedLinks);
    s.serviceNodes = w.serviceParameters.size();
    if (s.originalLinks > 0)
        s.reductionRatio = T(num_traits<T>::from_int(s.linksReduced) /
                             num_traits<T>::from_int(static_cast<long>(s.originalLinks)));
    return s;
}

}  // namespace wf
}  // namespace line

#endif  // LINE_API_WF_WF_PATTERN_UPDATER_H
