/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SPN_SPN_METRICS_H
#define LINE_API_SPN_SPN_METRICS_H

/**
 * Stationary measures of a product-form stochastic Petri net from the MDD-rec
 * masses.
 *
 * S. Balsamo, A. Marin, I. Stojic, FGCS 111 (2020) 475-490, Sec. 3.1 for the
 * definitions and Sec. 5.3 for the recursions they are read off.
 *
 *   n(P_j) = sum_k k P(m_j = k)              mean tokens
 *   u(P_j) = 1 - P(m_j = 0)                  place utilization
 *   u(T_j) = P(e_j >= 1)                     transition utilization
 *   x(T_j) = sum_k min(k, c_j) W(T_j) P(e_j = k)   throughput
 *   x(P_j) = sum_T I_j(T) x(T)               tokens removed per unit time
 *
 * ONE DEVIATION FROM THE PAPER'S x(T_j), AND IT IS A GENERALISATION. The paper
 * writes x(T_j) = sum_k k W(T_j) P(e_j = k), which is INFINITE-SERVER firing
 * semantics -- every enabling set fires in parallel. LINE's own rate law is
 * min(enabling degree, nmodeservers) * W(T), so `c_j` above is the mode's server
 * count: c_j = 1 recovers single-server semantics, x = W(T) P(e >= 1), and
 * c_j = infinity recovers the paper's formula exactly. Using the paper's form
 * for a single-server mode would report a throughput that grows with the token
 * population of a net whose transition can only fire one set at a time.
 *
 * The measures come out of ONE reachable set and ONE set of g_l, so they are
 * mutually consistent by construction: no per-measure fixed point, no iteration.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/mdd/mdd.h"
#include "line/api/mdd/mdd_rec.h"
#include "line/api/spn/spn_mdd.h"
#include "line/api/spn/spn_rec_enabled.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace spn {

/** The stationary measures of Sec. 3.1, per place level and per mode. */
template <class T>
struct SpnMetrics {
    /** The normalising constant G the measures are taken against. */
    T G;
    /** Mean tokens per place level. */
    std::vector<T> tokens;
    /** Place utilization, P(m_j > 0). */
    std::vector<T> place_util;
    /** Place throughput, tokens removed per unit time. */
    std::vector<T> place_tput;
    /** Transition (mode) utilization, P(e_j >= 1). */
    std::vector<T> mode_util;
    /** Transition (mode) throughput. */
    std::vector<T> mode_tput;
    /** marginal[l][k] = P(m_l = k). */
    std::vector<std::vector<T>> marginal;
};

/**
 * Every measure of Sec. 3.1 from one diagram and one product form.
 *
 * @param mdds the reachable set built by `spn_mdd`
 * @param g per-level product-form factors g_l(v)
 * @param info the metadata `spn_mdd` returned alongside the diagram
 */
template <class T>
SpnMetrics<T> spn_metrics(const mdd::MddStruct& mdds, const std::vector<std::vector<T>>& g,
                          const SpnInfo<T>& info) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t L = info.nplacelevels;

    SpnMetrics<T> out;
    out.G = mdd::mdd_rec(mdds, g);
    if (!(num_traits<T>::to_double(out.G) > 0))
        throw NumericError("spn_metrics: the normalising constant is not positive; the g_l passed "
                           "do not describe a product form over this reachable set");

    out.marginal.assign(L, std::vector<T>());
    out.tokens.assign(L, zero);
    out.place_util.assign(L, zero);
    out.place_tput.assign(L, zero);
    for (std::size_t l = 0; l < L; ++l) {
        const std::vector<T> mass = mdd::mdd_rec_marginal(mdds, g, l);
        out.marginal[l].assign(mass.size(), zero);
        for (std::size_t k = 0; k < mass.size(); ++k) {
            out.marginal[l][k] = T(mass[k] / out.G);
            out.tokens[l] += T(num_traits<T>::from_int(static_cast<long>(k)) * out.marginal[l][k]);
        }
        out.place_util[l] = T(one - out.marginal[l][0]);
    }

    out.mode_util.assign(info.modes.size(), zero);
    out.mode_tput.assign(info.modes.size(), zero);
    for (std::size_t e = 0; e < info.modes.size(); ++e) {
        const SpnMode<T>& mde = info.modes[e];
        const SpnEnabling<T> en = spn_rec_enabled(mdds, g, mde, L);
        out.mode_util[e] = T(en.ge[1] / out.G);
        // W(T) is the scalar firing rate of the mode; a phase-type firing time
        // has no single rate, so its throughput is left to the phase-level
        // marginal rather than reported through this formula.
        if (mde.nph > 1)
            throw UnsupportedError(
                "spn_metrics: mode " + std::to_string(mde.mode + 1) + " of node " +
                std::to_string(mde.trans) +
                " has a phase-type firing time, whose throughput is not W(T) times an enabling "
                "probability; read it from the phase-level marginal instead");
        // The formula below is W(T)*E[min(enabling degree, servers)], which is
        // the rate law only when no marking-dependent multiplier is in play.
        // With one, the firing rate is not a function of the enabling degree at
        // all, so the enabling-degree law is the wrong summary to take it from.
        // The MARGINALS above are unaffected -- they come from the product form,
        // not the rates.
        if (mde.dep)
            throw UnsupportedError(
                "spn_metrics: mode " + std::to_string(mde.mode + 1) + " of node " +
                std::to_string(mde.trans) +
                " has a marking-dependent firing rate, so its throughput is not W(T) times a "
                "function of the enabling degree and cannot be read from the enabling-degree "
                "law. The token marginals are still exact");
        const T rate = mde.D1(0, 0);
        T x = zero;
        for (std::size_t k = 1; k < en.eq.size(); ++k) {
            const double kd = static_cast<double>(k);
            const double served = std::isinf(mde.srv) ? kd : (kd < mde.srv ? kd : mde.srv);
            x += T(num_traits<T>::from_double(served) * rate * T(en.eq[k] / out.G));
        }
        out.mode_tput[e] = x;
        for (std::size_t l = 0; l < L; ++l)
            if (mde.enab[l] > 0)
                out.place_tput[l] += T(num_traits<T>::from_double(mde.enab[l]) * x);
    }
    return out;
}

}  // namespace spn
}  // namespace line

#endif  // LINE_API_SPN_SPN_METRICS_H
