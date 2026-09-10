/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SPN_SPN_REC_ENABLED_H
#define LINE_API_SPN_SPN_REC_ENABLED_H

/**
 * Enabling-degree distribution of one mode of a product-form stochastic Petri
 * net, by the masked MDD-rec recursion.
 *
 * S. Balsamo, A. Marin, I. Stojic, FGCS 111 (2020) 475-490, Sec. 5.3.
 *
 * The enabling degree of a mode in marking m is
 *
 *     e(m) = min_{l : I_l > 0} floor(m_l / I_l),
 *
 * zero when any inhibitor threshold is met. P(e >= k) is therefore the mass of
 * the marking subset in which EVERY input level holds at least k*I_l tokens and
 * no inhibitor fires, which is a per-level restriction and so exactly what
 * `mdd::mdd_rec_masked` computes: the paper's second modified recurrence is the
 * same walk under a different mask, not a second algorithm.
 *
 * The masses returned are UNNORMALISED, as in the paper; divide by G from
 * `mdd::mdd_rec` for probabilities. `spn_metrics` does that and turns them into
 * the transition measures.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mdd/mdd.h"
#include "line/api/mdd/mdd_rec.h"
#include "line/api/spn/spn_mdd.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace spn {

/** Unnormalised enabling-degree masses of one mode. */
template <class T>
struct SpnEnabling {
    /** ge[k] is the mass of {e >= k}; ge[0] is the whole reachable set. */
    std::vector<T> ge;
    /** eq[k] is the mass of {e == k}, i.e. ge[k] - ge[k+1]. */
    std::vector<T> eq;
    /** The largest enabling degree the place bounds permit, E_j in the paper. */
    std::size_t max_degree = 0;
};

/**
 * Enabling-degree masses of mode `mde` over the reachable set in `mdds`.
 *
 * @param mdds the reachable set built by `spn_mdd`
 * @param g per-level product-form factors, one vector per level
 * @param mde the mode, as returned in `SpnInfo::modes`
 * @param nplacelevels how many leading levels are place levels
 */
template <class T>
SpnEnabling<T> spn_rec_enabled(const mdd::MddStruct& mdds, const std::vector<std::vector<T>>& g,
                               const SpnMode<T>& mde, std::size_t nplacelevels) {
    if (nplacelevels > mdds.K)
        throw InputError("spn_rec_enabled: more place levels than diagram levels");
    const T zero = num_traits<T>::from_int(0);

    // E_j: the enabling degree cannot exceed what the tightest input place bound
    // allows. A mode with no input place has no bound and is refused rather than
    // silently truncated, matching spn_mdd's own refusal.
    std::size_t emax = 0;
    bool has_input = false;
    for (std::size_t l = 0; l < nplacelevels; ++l) {
        if (!(mde.enab[l] > 0)) continue;
        const double top = static_cast<double>(mdds.domain[l] - 1);
        const std::size_t cap = static_cast<std::size_t>(std::floor(top / mde.enab[l]));
        emax = has_input ? (cap < emax ? cap : emax) : cap;
        has_input = true;
    }
    if (!has_input)
        throw InputError("spn_rec_enabled: the mode consumes from no place, so its enabling "
                         "degree is unbounded");

    SpnEnabling<T> out;
    out.max_degree = emax;
    out.ge.assign(emax + 2, zero);
    out.eq.assign(emax + 2, zero);
    for (std::size_t k = 0; k <= emax; ++k) {
        mdd::MddMask mask(mdds.K);
        for (std::size_t j = 0; j < mdds.K; ++j)
            mask[j].assign(static_cast<std::size_t>(mdds.domain[j]), true);
        for (std::size_t l = 0; l < nplacelevels; ++l) {
            const double need = mde.enab[l] * static_cast<double>(k);
            for (int v = 0; v < mdds.domain[l]; ++v) {
                const bool short_of_tokens = static_cast<double>(v) < need;
                const bool inhibited = static_cast<double>(v) >= mde.inhib[l];
                // k = 0 asks only that the marking exist, so the inhibitor test
                // belongs to k >= 1: e = 0 covers the inhibited markings too.
                if (short_of_tokens || (k > 0 && inhibited))
                    mask[l][static_cast<std::size_t>(v)] = false;
            }
        }
        out.ge[k] = mdd::mdd_rec_masked(mdds, g, mask);
    }
    for (std::size_t k = 0; k <= emax; ++k) out.eq[k] = T(out.ge[k] - out.ge[k + 1]);
    return out;
}

}  // namespace spn
}  // namespace line

#endif  // LINE_API_SPN_SPN_REC_ENABLED_H
