/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#pragma once

#include <cmath>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"

namespace line {
namespace solvers {

/** One [F(t), t] curve of the base-class response-time CDF fallback. */
struct DefaultCdfCurve {
    std::vector<double> t;  ///< time points
    std::vector<double> F;  ///< CDF values, F[j] = P(R <= t[j])
};

/**
 * The NetworkSolver base-class response-time CDF: an exponential law with the
 * right mean per (station, class), tabulated on 100 quantile points.
 *
 * This is `@@NetworkSolver/getCdfRespT.m` verbatim -- "a trivial approximation
 * that assumes exponential distributions everywhere with mean as RN(i,r)" --
 * and it is what the reference serves for every solver without a
 * distributional result of its own (MVA, QNS, BA, AG). It says nothing about
 * the tail; the solvers with a real law (CTMC, NC, Fluid, MAM, JMT, LDES)
 * never reach it, and SSA refuses instead of inheriting it.
 *
 * A Source station is left with empty curves, as the reference leaves its
 * cells empty; a served cell whose mean is non-finite or non-positive gets the
 * reference's degenerate [1, 0] point mass at zero.
 */
template <class T>
std::vector<std::vector<DefaultCdfCurve>> solver_default_cdf_respt(
    const qn::NetworkStruct<T>& sn, const Matrix<T>& RN) {
    const std::size_t M = sn.nstations;
    const std::size_t K = sn.nclasses;
    std::vector<std::vector<DefaultCdfCurve>> RD(M, std::vector<DefaultCdfCurve>(K));
    const std::size_t npts = 100;
    for (std::size_t i = 0; i < M; ++i) {
        if (sn.stations[i].nodetype == lang::NodeType::Source) continue;
        for (std::size_t c = 0; c < K; ++c) {
            const double rn = (i < static_cast<std::size_t>(RN.rows()) &&
                               c < static_cast<std::size_t>(RN.cols()))
                                  ? num_traits<T>::to_double(RN(i, c))
                                  : 0.0;
            DefaultCdfCurve& cell = RD[i][c];
            if (std::isfinite(rn) && rn > 0.0) {
                cell.t.reserve(npts);
                cell.F.reserve(npts);
                for (std::size_t j = 0; j < npts; ++j) {
                    const double q =
                        0.001 + (0.999 - 0.001) * static_cast<double>(j) / (npts - 1);
                    cell.F.push_back(q);
                    cell.t.push_back(-std::log(1.0 - q) * rn);
                }
            } else {
                cell.F.push_back(1.0);
                cell.t.push_back(0.0);
            }
        }
    }
    return RD;
}

}  // namespace solvers
}  // namespace line
