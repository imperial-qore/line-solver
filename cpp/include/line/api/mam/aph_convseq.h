/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_APH_CONVSEQ_H
#define LINE_API_MAM_APH_CONVSEQ_H

/**
 * Convolution of a sequence of matrix-exponential laws.
 *
 * Port of `matlab/lib/kpctoolbox/aph/aph_convseq.m`: fold `aph_simplify` with
 * the sequence pattern over the list, left to right. The reference takes a flat
 * cell array of alternating alpha and T entries and special-cases a list of one
 * pair by returning it unchanged; here the list is a vector of pairs, so the
 * one-element case falls out of the fold and needs no arm of its own.
 *
 * The composite order is the SUM of the orders, so a long activity sequence of
 * high-order fits produces a large generator. That is the reference's cost too
 * -- neither it nor this reduces the representation -- and it is why the LN
 * caller fits each term to a low-order APH before convolving.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/aph_simplify.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace mam {

/** Convolve the sequence, i.e. the law of the sum of independent terms. */
template <class T>
AphPair<T> aph_convseq(const std::vector<AphPair<T>>& seq) {
    if (seq.empty()) throw InputError("aph_convseq: the sequence is empty");
    const T one = num_traits<T>::from_int(1);
    AphPair<T> acc = seq[0];
    for (std::size_t i = 1; i < seq.size(); ++i)
        acc = aph_simplify(acc, seq[i], one, one, AphPattern::Sequence);
    return acc;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_APH_CONVSEQ_H
