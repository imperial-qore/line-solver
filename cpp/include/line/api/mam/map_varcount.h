/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_VARCOUNT_H
#define LINE_API_MAM_MAP_VARCOUNT_H

/**
 * Variance of the counts of a MAP over windows of length t, in the spelling of
 * matlab/lib/kpctoolbox/map/map_varcount.m.
 *
 * With Q = D0 + D1, pi the stationary phase vector, lambda = pi D1 e and
 * tmp = (e pi - Q)^-1,
 *
 *   Var[N(t)] = (lambda - 2 lambda^2 + 2 pi D1 tmp D1 e) t
 *               - 2 pi D1 (I - exp(Q t)) tmp^2 D1 e.
 *
 * This is algebraically identical to map_count_var.m (see map_count_var.h):
 * e pi and Q commute because Q e = 0 and pi Q = 0, so tmp commutes with
 * exp(Q t) and tmp (I - exp) tmp = (I - exp) tmp^2. Both spellings are ported
 * because both are called from the MATLAB tree -- map_acfc.m normalizes by
 * map_varcount, not by map_count_var -- and their agreement is a useful check
 * on the exponential.
 *
 * ARITHMETIC: transcendental, as map_count_var.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_count_var.h"
#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * @param m    the MAP (D0, D1)
 * @param tset window lengths
 * @return Var[N(t)] for each window length, in the order of tset
 */
template <class T>
std::vector<T> map_varcount(const Map<T>& m, const std::vector<T>& tset) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_varcount requires transcendental arithmetic");
    const std::size_t n = m.order();
    const Matrix<T> Q = map_infgen(m);
    const std::vector<T> piq = map_prob(m);
    const Matrix<T> tmp = detail::map_count_deviation(Q, piq);
    const Matrix<T> tmp2 = matmul(tmp, tmp);
    const std::vector<T> e = ones<T>(n);
    const std::vector<T> D1e = mulvec(m.D1, e);

    const std::vector<T> piD1 = vecmul(piq, m.D1);  // row vector pi D1
    T lam = num_traits<T>::from_int(0);
    for (const T& v : piD1) lam += v;

    const std::vector<T> piD1tmp = vecmul(piD1, tmp);
    T pre_corr = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) pre_corr += piD1tmp[i] * D1e[i];
    const T two = num_traits<T>::from_int(2);
    const T pre = lam - two * lam * lam + two * pre_corr;

    const std::vector<T> tail = mulvec(tmp2, D1e);  // tmp^2 D1 e

    std::vector<T> out;
    out.reserve(tset.size());
    for (std::size_t k = 0; k < tset.size(); ++k) {
        if (tset[k] < num_traits<T>::from_int(0))
            throw InputError("map_varcount: negative window length");
        const Matrix<T> E = expm(Q, tset[k]);
        const std::vector<T> piD1E = vecmul(piD1, E);
        T post = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < n; ++i) post += (piD1[i] - piD1E[i]) * tail[i];
        out.push_back(pre * tset[k] - two * post);
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_VARCOUNT_H
