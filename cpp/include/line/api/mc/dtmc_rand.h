/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_DTMC_RAND_H
#define LINE_API_MC_DTMC_RAND_H

/**
 * Random DTMC kernels, trajectory simulation and the weak-component split.
 *
 * Templated port of matlab/lib/kpctoolbox/mc: dtmc_rand.m, dtmc_simulate.m and
 * weaklyconncomp.m.
 *
 * dtmc_rand is defined as the uniformization of a random generator, so its
 * kernel always has a nonzero diagonal and never mixes at rate one. Building
 * the rows directly from normalised uniforms would look equivalent and is not:
 * the self-loop probability of the uniformized chain is 1 + q_ii/q, which
 * concentrates near one on the fast states.
 *
 * weaklyconncomp goes through dmperm in MATLAB and through a graph traversal
 * here. Component LABELS are therefore not comparable across the two, only the
 * partition is; the labels here are assigned in order of the smallest member,
 * which is the one canonical choice that is stable under recompilation.
 */

#include <algorithm>
#include <cstddef>
#include <random>
#include <vector>

#include "line/api/mc/ctmc_rand.h"
#include "line/api/mc/ctmc_randomization.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/** Random stochastic matrix, the uniformization of a random generator. */
template <class T, class Gen>
Matrix<T> dtmc_rand(std::size_t n, Gen& gen) {
    return ctmc_randomization(ctmc_rand<T>(n, gen)).P;
}

/**
 * Sample path of a DTMC, n states starting from pi0.
 *
 * An absorbing state ends the path early, exactly as the reference does: it
 * returns the prefix rather than padding, so the returned length is at most n.
 */
template <class T, class Gen>
std::vector<std::size_t> dtmc_simulate(const Matrix<T>& P, const std::vector<T>& pi0,
                                       std::size_t n, Gen& gen) {
    const std::size_t m = P.rows();
    if (P.cols() != m) throw InputError("dtmc_simulate: P is not square");
    if (pi0.size() != m) throw InputError("dtmc_simulate: pi0 does not match the state space");
    const T zero = num_traits<T>::from_int(0);
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    std::size_t st = m;
    {
        const T r = num_traits<T>::from_double(unif(gen));
        T acc = zero;
        for (std::size_t i = 0; i < m; ++i) {
            acc += pi0[i];
            if (pi0[i] > zero && r < acc) {
                st = i;
                break;
            }
        }
        if (st == m)
            for (std::size_t i = 0; i < m; ++i)
                if (pi0[i] > zero) {
                    st = i;
                    break;
                }
        if (st == m) throw InputError("dtmc_simulate: pi0 puts no mass on any state");
    }

    std::vector<std::size_t> sts;
    sts.reserve(n);
    for (std::size_t k = 0; k < n; ++k) {
        sts.push_back(st);
        T rowsum = zero;
        for (std::size_t j = 0; j < m; ++j) rowsum += P(st, j);
        if (rowsum == zero || P(st, st) == num_traits<T>::from_int(1)) return sts;
        const T r = num_traits<T>::from_double(unif(gen));
        T acc = zero;
        std::size_t nxt = m;
        for (std::size_t j = 0; j < m; ++j) {
            acc += P(st, j);
            if (P(st, j) > zero && r < acc) {
                nxt = j;
                break;
            }
        }
        if (nxt == m)
            for (std::size_t j = 0; j < m; ++j)
                if (P(st, j) > zero) {
                    nxt = j;
                    break;
                }
        st = nxt;
    }
    return sts;
}

/** Number of weakly connected components and the per-node label (weaklyconncomp.m). */
template <class T>
struct WeakCompResult {
    std::size_t count;              ///< S, the number of components
    std::vector<std::size_t> comp;  ///< C, zero-based component label per node
};

/** Weakly connected components of the graph whose adjacency is the support of G. */
template <class T>
WeakCompResult<T> weaklyconncomp(const Matrix<T>& G) {
    const std::vector<std::vector<std::size_t>> parts = detail::weak_components(G);
    WeakCompResult<T> out;
    out.count = parts.size();
    out.comp.assign(G.rows(), 0);
    for (std::size_t c = 0; c < parts.size(); ++c)
        for (std::size_t k = 0; k < parts[c].size(); ++k) out.comp[parts[c][k]] = c;
    return out;
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_DTMC_RAND_H
