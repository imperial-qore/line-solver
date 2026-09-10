/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFF_DELAY_H
#define LINE_API_PFQN_PFF_DELAY_H

/**
 * Product-form factor of a delay station.
 *
 * Templated port of jar/src/main/java/jline/api/pfqn/nc/Pfqn_pff_delay.java.
 * MATLAB has no standalone counterpart: the same expression is inlined wherever
 * a normalizing-constant recursion folds in the delay term.
 *
 *   F_Z(n) = prod_r Z_r^{n_r} / n_r!
 *
 * the contribution of the delay stations to the normalizing constant of a
 * closed product-form network. The empty population gives 1, and a class with a
 * positive population but no think time gives 0, the delay being unreachable
 * for it.
 *
 * The product is accumulated in the LOG DOMAIN and exponentiated once, so a
 * large population does not overflow through the intermediate powers even where
 * the value itself is representable.
 *
 * Arithmetic: TRANSCENDENTAL, because of that log-domain accumulation. The same
 * quantity in exact arithmetic is available from the delay balance function
 * inside pfqn_ca, which forms it as a ratio of exact rationals.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/**
 * @param Z (R) think times of the delay station
 * @param n (R) population of each class
 */
template <class T>
T pfqn_pff_delay(const std::vector<T>& Z, const std::vector<int>& n) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_pff_delay accumulates in the log domain and needs transcendental "
                  "arithmetic");
    using std::exp;
    using std::log;
    const std::size_t R = n.size();
    if (Z.size() != R)
        throw InputError("pfqn_pff_delay: Z and n disagree on the class count");
    const T zero = num_traits<T>::from_int(0);
    long total = 0;
    for (std::size_t r = 0; r < R; ++r) total += n[r];
    if (total == 0) return num_traits<T>::from_int(1);

    T f = zero;
    for (std::size_t r = 0; r < R; ++r) {
        if (Z[r] > zero) {
            f += log(Z[r]) * num_traits<T>::from_int(n[r]);
            f -= detail::num_factln<T>(num_traits<T>::from_int(n[r]));
        } else if (n[r] > 0) {
            return zero;
        }
    }
    return exp(f);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFF_DELAY_H
