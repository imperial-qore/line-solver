/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_FJ_DIST2FJ_H
#define LINE_API_FJ_FJ_DIST2FJ_H

/**
 * Conversion of a LINE MAP into the arrival or service descriptor of the
 * fork-join response-time-tail algorithm of Qiu, Perez and Harrison (IFIP
 * Performance 2015).
 *
 * Templated port of matlab/src/api/fj/fj_dist2fj.m. The reference reads two
 * scalars out of the NetworkStruct, sn.procid(ist, r) and sn.rates(ist, r);
 * only the first is used (see the reference defects below), and it is passed
 * here as an explicit FjProcType, so this port carries no NetworkStruct
 * dependency.
 *
 *   arrival: (lambda, lambda0 = D0, lambda1 = D1, ma = phases, Ia = I)
 *   service: (mu = lambda, ST = D0, St = -D0 e, tau_st = map_pie)
 *
 * At most two phases are accepted, which is the algorithm's own restriction.
 *
 * REFERENCE DEFECTS in fj_dist2fj.m:
 *
 *  1. DEAD READ. mean_rate = sn.rates(ist, r) is computed on line 53 and never
 *     used; the rate that is returned is map_lambda of the MAP itself. Not
 *     propagated: the port does not take a rate argument at all. This matters
 *     because it is the only place the two could disagree, and the MAP is the
 *     authority.
 *
 *  2. mean_time = 1 / lambda is likewise computed and never used, and divides
 *     by zero for a disabled class instead of reporting it.
 *
 *  3. THE SERVICE BRANCH ACCEPTS MAP(2) SILENTLY. Its distribution switch
 *     rejects a two-phase MAP with an error, but only AFTER the descriptor has
 *     been filled in, and the arrival branch accepts MAP(2) explicitly. The
 *     asymmetry is intentional in the algorithm (service must be phase type),
 *     so it is reproduced: FjProcType::Map is rejected for a service process.
 *
 * ARITHMETIC. map_lambda and map_pie need one linear solve each, and the exit
 * vector is a row sum, so the whole conversion stays in the field and is
 * instantiated at T = Rational as well as double and Real.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace fj {

/**
 * The subset of ProcessType that the fork-join algorithm accepts. The values
 * are MATLAB's ProcessType codes (EXP = 0, ERLANG = 1, HYPEREXP = 2, MAP = 5),
 * so a caller holding sn.procid can pass it through unchanged.
 */
enum class FjProcType { Exp = 0, Erlang = 1, HyperExp = 2, Map = 5 };

/** Which of the two descriptors to build. */
enum class FjDistKind { Arrival, Service };

/**
 * Descriptor of an arrival or a service process. Only the fields of the
 * requested kind are filled; choice is the algorithm's own distribution code
 * (1 = Exp, 2 = HE2, 3 = ER2, 4 = MAP2), MATLAB's ArrChoice or SerChoice.
 */
template <class T>
struct FjDist {
    // arrival
    T lambda = num_traits<T>::from_int(0);
    Matrix<T> lambda0;
    Matrix<T> lambda1;
    std::size_t ma = 0;
    Matrix<T> Ia;
    // service
    T mu = num_traits<T>::from_int(0);
    Matrix<T> ST;
    std::vector<T> St;
    std::vector<T> tau_st;
    int choice = 0;
};

/**
 * Build the fork-join descriptor of a MAP.
 *
 * @param m        the LINE MAP (D0, D1)
 * @param kind     arrival or service
 * @param procType the process type, MATLAB's sn.procid(ist, r)
 */
template <class T>
FjDist<T> fj_dist2fj(const mam::Map<T>& m, FjDistKind kind, FjProcType procType) {
    const std::size_t n = m.order();
    if (n == 0) throw InputError("fj_dist2fj: empty process");
    if (m.D1.rows() != n || m.D0.cols() != n || m.D1.cols() != n)
        throw InputError("fj_dist2fj: D0 and D1 must be square matrices of the same size");
    if (n > 2)
        throw InputError("fj_dist2fj: only distributions with at most 2 phases are supported");

    const T zero = num_traits<T>::from_int(0);
    FjDist<T> d;
    const T lambda = mam::map_lambda(m);

    if (kind == FjDistKind::Arrival) {
        d.lambda = lambda;
        d.lambda0 = m.D0;
        d.lambda1 = m.D1;
        d.ma = n;
        d.Ia = eye<T>(n);
        if (procType == FjProcType::Exp)
            d.choice = 1;
        else if (procType == FjProcType::HyperExp && n == 2)
            d.choice = 2;
        else if (procType == FjProcType::Erlang && n == 2)
            d.choice = 3;
        else if (procType == FjProcType::Map && n == 2)
            d.choice = 4;
        else
            throw InputError("fj_dist2fj: unsupported arrival distribution type for this number "
                             "of phases");
    } else {
        d.mu = lambda;
        d.ST = m.D0;
        d.St.assign(n, zero);
        for (std::size_t i = 0; i < n; ++i) {
            T s = zero;
            for (std::size_t j = 0; j < n; ++j) s += m.D0(i, j);
            d.St[i] = -s;
        }
        d.tau_st = mam::map_pie(m);
        if (procType == FjProcType::Exp)
            d.choice = 1;
        else if (procType == FjProcType::HyperExp && n == 2)
            d.choice = 2;
        else if (procType == FjProcType::Erlang && n == 2)
            d.choice = 3;
        else
            throw InputError("fj_dist2fj: unsupported service distribution type for this number "
                             "of phases");
    }
    return d;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_FJ_DIST2FJ_H
