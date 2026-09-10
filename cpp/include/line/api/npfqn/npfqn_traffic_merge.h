/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_NPFQN_TRAFFIC_MERGE_H
#define LINE_API_NPFQN_TRAFFIC_MERGE_H

/**
 * Superposition of several marked arrival flows into one.
 *
 * Templated port of matlab/src/api/npfqn/npfqn_traffic_merge.m and
 * matlab/src/api/npfqn/npfqn_traffic_merge_cs.m.
 *
 * Given n MMAPs carrying the same R classes, the merge is the class-by-class
 * superposition mmap_super(., ., 'match'), folded left to right, followed
 * optionally by a compression back to a small representation. The class-
 * switching variant first re-marks flow i with its own (R x R) switching
 * matrix, prob((i-1)R + r, s) being the probability that a class-r arrival
 * from flow i leaves as class s, and then superposes.
 *
 * WHAT MUST HOLD. Superposition of independent flows adds the rates: the
 * merged per-class rate is the sum of the per-class rates of the operands, as
 * an identity and not as an approximation. The class-switching variant carries
 * the rates through the switching matrix, so the merged class-s rate is
 * sum_i sum_r lambda_{i,r} P_i(r, s), and the TOTAL rate is preserved whenever
 * every P_i is stochastic. Both are asserted exactly in the rational
 * instantiation by the tests. Compression does NOT preserve them exactly: it
 * preserves the class probabilities p_c and the aggregate mean, hence the
 * per-class rates, only up to the APH(2) moment fit.
 *
 * REFERENCE DEFECTS in matlab/src/api/npfqn/npfqn_traffic_merge.m:
 *
 *  1. A ONE-ARGUMENT CALL WITH MORE THAN ONE NON-EMPTY FLOW ALWAYS ERRORS.
 *     Line 9 builds the default configuration as struct('merge','default'),
 *     with no compress field, and line 46 then reaches `switch config.compress`
 *     unconditionally. MATLAB raises "Reference to non-existent field
 *     'compress'". The single-flow case returns early at line 13 and is
 *     unaffected, which is why the defect survives: the n == 1 shortcut is by
 *     far the most common call. Reproduction, from matlab/:
 *         M = map_exponential(1); M{3} = M{2};
 *         npfqn_traffic_merge({M, M})
 *     -> Reference to non-existent field 'compress'.
 *     while npfqn_traffic_merge({M}) returns M.
 *     The INTENDED behaviour is the documented default of the compress switch,
 *     i.e. compression by mmap_compress, and that is what this port does with
 *     its default-constructed configuration. The defect is NOT propagated:
 *     there is no way to spell "merge configured but compression unset" here,
 *     since MergeConfig::compress is a value with a default.
 *
 *  2. The 'default'/'super' branch guards each operand with ~isempty(MMAP{j})
 *     although the empty operands were already deleted at line 5, so the guard
 *     is dead. Reproduced as written (the port drops empty flows up front and
 *     the loop then has nothing to skip), and harmless.
 *
 * ARITHMETIC. The merge itself is a Kronecker sum and stays exact at
 * T = Rational, so npfqn_traffic_merge_cs and the no-compression merge are
 * offered at every arithmetic. Compression pulls in aph2_fit, so
 * npfqn_traffic_merge with Compress::Default is gated on
 * num_traits<T>::has_transcendental through mmap_compress.
 *
 * NOT PORTED. Merge::Mixture needs mmap_mixture_fit_mmap, which is not in this
 * tree; it raises UnsupportedError naming the missing MATLAB function rather
 * than falling back to another merge. Merge::Interpos is served, through
 * m3pp2m_fitc_theoretical and m3pp2m_interleave; it is gated on transcendental
 * arithmetic like the rest of the counting-process fitters, so it is not
 * offered at T = Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/m3pp2m_interleave.h"
#include "line/api/mam/mmap_compress.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace npfqn {

/**
 * Superposition matching classes one to one (mmap_super.m, option 'match').
 *
 * Every component of the MATLAB cell, D0 and D1 included, is combined with the
 * same Kronecker sum, so the result carries the same class list as its
 * operands. This is NOT line::mam::mmap_super, which is the 'default' option
 * and CONCATENATES the two class lists; the merge needs the matching form,
 * because merging n flows of R classes must yield R classes and not n R.
 * MATLAB errors when the class counts differ; so does this.
 *
 * It lives in the npfqn namespace, next to its only caller, rather than in the
 * mam domain: the mam MMAP algebra is owned elsewhere in this tree and this
 * option is not part of it.
 */
template <class T>
mam::Mmap<T> mmap_super_match(const mam::Mmap<T>& a, const mam::Mmap<T>& b) {
    if (a.classes() != b.classes())
        throw InputError("mmap_super_match: class matching failed, the MMAPs have different "
                         "numbers of classes");
    mam::Mmap<T> s;
    s.D0 = mam::krons(a.D0, b.D0);
    s.D1 = mam::krons(a.D1, b.D1);
    s.Dc.reserve(a.classes());
    for (std::size_t c = 0; c < a.classes(); ++c) s.Dc.push_back(mam::krons(a.Dc[c], b.Dc[c]));
    return mam::mmap_normalize(s);
}

/**
 * Re-mark an MMAP's K types into R classes (mmap_mark.m).
 *
 * @param m    an MMAP with K marked types
 * @param prob (K x R); prob(k, r) is the probability that a type-k arrival is
 *             marked as class r
 *
 * D0 and D1 are untouched and D1^(r) = sum_k D1^(k) prob(k, r). MATLAB does not
 * normalize here, so neither does this: the caller's mmap_super does it.
 *
 * Distinct from line::mam::mmap_mark, which marks a plain MAP with
 * phase-dependent weights; MATLAB's mmap_mark is this one.
 */
template <class T>
mam::Mmap<T> mmap_mark_types(const mam::Mmap<T>& m, const Matrix<T>& prob) {
    const std::size_t K = prob.rows(), R = prob.cols();
    if (K != m.classes())
        throw InputError("mmap_mark_types: prob needs one row per marked type of the MMAP");
    const T zero = num_traits<T>::from_int(0);
    mam::Mmap<T> out;
    out.D0 = m.D0;
    out.D1 = m.D1;
    const std::size_t n = m.order();
    for (std::size_t r = 0; r < R; ++r) {
        Matrix<T> Dr(n, n, zero);
        for (std::size_t k = 0; k < K; ++k) {
            const T& p = prob(k, r);
            if (p == zero) continue;
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j) Dr(i, j) += m.Dc[k](i, j) * p;
        }
        out.Dc.push_back(Dr);
    }
    return out;
}

/** Merge rule, MATLAB's config.merge. */
enum class Merge {
    Default,  ///< 'default', identical to 'super'
    Super,    ///< 'super'
    Mixture,  ///< 'mixture', not ported
    Interpos  ///< 'interpos', lumped interleaving of per-flow M3PP(2, m) fits
};

/** Post-merge compression, MATLAB's config.compress. */
enum class Compress {
    Default,  ///< 'default', mmap_compress with its own default method
    None      ///< 'none'
};

/** MATLAB's config struct. The defaults are the documented intended ones. */
struct MergeConfig {
    Merge merge = Merge::Default;
    Compress compress = Compress::Default;
};

namespace detail {

/**
 * Compression, if the arithmetic can carry it. The exact instantiation cannot:
 * mmap_compress fits an APH(2) and needs square roots. Rather than making
 * npfqn_traffic_merge uninstantiable at T = Rational -- which would remove the
 * exact merge too, and that one is a pure Kronecker sum -- the branch is
 * selected at compile time and the exact build refuses only at run time, and
 * only when compression is actually requested.
 */
template <class T>
mam::Mmap<T> compress_or_refuse(const mam::Mmap<T>& s) {
    if constexpr (num_traits<T>::has_transcendental) {
        return mam::mmap_compress(s, mam::MmapCompressMethod::MixtureOrder1);
    } else {
        throw UnsupportedError("npfqn_traffic_merge: compression needs transcendental arithmetic "
                               "(aph2_fit); merge at exact arithmetic with Compress::None");
    }
}

}  // namespace detail

/**
 * Merge a list of MMAPs carrying the same classes.
 *
 * @param flows the MMAPs to superpose; empty ones are dropped, as in MATLAB
 * @param config merge rule and compression
 * @return the merged MMAP, normalized
 */
template <class T>
mam::Mmap<T> npfqn_traffic_merge(const std::vector<mam::Mmap<T>>& flows,
                                 const MergeConfig& config = MergeConfig()) {
    std::vector<const mam::Mmap<T>*> nonEmpty;
    for (const mam::Mmap<T>& f : flows)
        if (f.order() != 0) nonEmpty.push_back(&f);
    if (nonEmpty.empty()) throw InputError("npfqn_traffic_merge: no non-empty flow to merge");
    if (nonEmpty.size() == 1) return *nonEmpty.front();

    switch (config.merge) {
        case Merge::Default:
        case Merge::Super:
            break;
        case Merge::Mixture:
            throw UnsupportedError("npfqn_traffic_merge: merge 'mixture' needs "
                                   "mmap_mixture_fit_mmap, which is not ported");
        case Merge::Interpos: {
            // REFUSED AT RUNTIME, NOT AT COMPILE TIME. `m3pp2m_fitc_theoretical`
            // static_asserts on transcendental arithmetic, and a template
            // instantiates every branch of this switch whatever the runtime
            // `config.merge` is -- so calling it unguarded made the whole of
            // `npfqn_traffic_merge` uninstantiable at exact arithmetic, including
            // the Super and single-flow paths that never reach here.
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError("npfqn_traffic_merge: merge 'interpos' fits an "
                                       "M3PP(2,m), which is transcendental; use double or real");
            } else {
                // Every flow is first reduced to an M3PP(2, m) on its exact counting
                // characteristics, then the L of them are lumped onto one
                // birth-death phase process of order L + 1.
                std::vector<mam::Mmap<T>> flowFits;
                for (std::size_t j = 0; j < nonEmpty.size(); ++j)
                    flowFits.push_back(mam::m3pp2m_fitc_theoretical(
                        *nonEmpty[j], std::string("exact_delta"), num_traits<T>::from_int(1),
                        num_traits<T>::from_double(1e6)));
                mam::Mmap<T> lumped = mam::m3pp2m_interleave(flowFits);
                if (config.compress == Compress::Default)
                    lumped = detail::compress_or_refuse(lumped);
                return mam::mmap_normalize(lumped);
            }
        }
    }

    mam::Mmap<T> s = *nonEmpty.front();
    for (std::size_t j = 1; j < nonEmpty.size(); ++j) s = mmap_super_match(s, *nonEmpty[j]);

    if (config.compress == Compress::Default) s = detail::compress_or_refuse(s);
    return mam::mmap_normalize(s);
}

}  // namespace npfqn
}  // namespace line

#endif  // LINE_API_NPFQN_TRAFFIC_MERGE_H
