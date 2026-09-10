/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FES_MAP_INTERDEPARTURE_H
#define LINE_API_FES_MAP_INTERDEPARTURE_H

/**
 * Inter-departure MAP of a closed subnetwork made of one MAP station and one MAP
 * flow-equivalent server.
 *
 * Templated port of matlab/src/api/fes/fes_map_interdeparture.m, mirrored by the
 * JAR and native Python. Implements the block bidiagonal construction of Casale,
 * Mi, Cherkasova and Smirni, IEEE Trans. Soft. Eng. 37(5), 2011, Section 5.2.2:
 *
 *   T0 = [ D0 (x) I        0              0            0
 *          I (x) F1^1   D0 (+) F0^1       0            0
 *            ...           ...           ...          ...
 *            0        I (x) F1^{n-1} D0 (+) F0^{n-1}   0
 *            0            0        I (x) F1^n     I (x) F0^n ]
 *
 *   T1 = superdiagonal blocks D1 (x) I, last block row zero
 *
 * Level k is the population of the flow-equivalent server, so the station holds
 * n - k jobs and both processes may be load dependent. Marked transitions are the
 * completions of the station, which are the departures fed to the rest of the
 * model.
 *
 * ARITHMETIC: field operations only, exact at T = Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/fes/fes_map_levels.h"
#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace fes {

/**
 * @param maps per-level service processes of the station, index j-1 holding j jobs
 * @param fes  per-level processes of the flow-equivalent server, index k-1 holding k jobs
 * @param n    number of jobs circulating in the subnetwork
 * @return the pair (T0, T1) of the inter-departure MAP
 */
template <class T>
mam::Map<T> fes_map_interdeparture(const std::vector<mam::Map<T>>& maps,
                                   const std::vector<mam::Map<T>>& fes, std::size_t n) {
    if (n < 1) throw InputError("fes_map_interdeparture: the subnetwork population must be at least 1");
    const std::vector<mam::Map<T>> mapsLev = fes_map_levels(maps, n);
    const std::vector<mam::Map<T>> fesLev = fes_map_levels(fes, n);

    const std::size_t ms = mapsLev[0].order();
    const std::size_t mf = fesLev[0].order();
    const std::size_t blk = ms * mf;
    const std::size_t dim = (n + 1) * blk;
    const T zero = num_traits<T>::from_int(0);

    mam::Map<T> out;
    out.D0 = Matrix<T>(dim, dim, zero);
    out.D1 = Matrix<T>(dim, dim, zero);

    for (std::size_t k = 0; k <= n; ++k) {
        const std::size_t off = k * blk;
        const std::size_t j = n - k;
        // D0 (x) I on the station phase, I (x) F0 on the flow-equivalent one
        if (j > 0) {
            const Matrix<T>& D0 = mapsLev[j - 1].D0;
            for (std::size_t a = 0; a < ms; ++a)
                for (std::size_t b = 0; b < ms; ++b)
                    for (std::size_t f = 0; f < mf; ++f)
                        out.D0(off + a * mf + f, off + b * mf + f) += D0(a, b);
        }
        if (k > 0) {
            const Matrix<T>& F0 = fesLev[k - 1].D0;
            for (std::size_t a = 0; a < ms; ++a)
                for (std::size_t f = 0; f < mf; ++f)
                    for (std::size_t g = 0; g < mf; ++g)
                        out.D0(off + a * mf + f, off + a * mf + g) += F0(f, g);
            const Matrix<T>& F1 = fesLev[k - 1].D1;
            for (std::size_t a = 0; a < ms; ++a)
                for (std::size_t f = 0; f < mf; ++f)
                    for (std::size_t g = 0; g < mf; ++g)
                        out.D0(off + a * mf + f, off - blk + a * mf + g) += F1(f, g);
        }
        if (j > 0) {
            const Matrix<T>& D1 = mapsLev[j - 1].D1;
            for (std::size_t a = 0; a < ms; ++a)
                for (std::size_t b = 0; b < ms; ++b)
                    for (std::size_t f = 0; f < mf; ++f)
                        out.D1(off + a * mf + f, off + blk + b * mf + f) += D1(a, b);
        }
    }
    return out;
}

}  // namespace fes
}  // namespace line

#endif  // LINE_API_FES_MAP_INTERDEPARTURE_H
