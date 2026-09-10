/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MATLAB_ILT_H
#define LINE_API_MAM_MATLAB_ILT_H

/**
 * Numerical inverse Laplace transform in the Abate-Whitt framework, the port of
 * `matlab/lib/thirdparty/iltcme/matlab_ilt.m`.
 *
 * Every variant evaluates the same quadrature
 *
 *     f(t) ~= (1/t) sum_k Re( eta_k F(beta_k / t) )
 *
 * and differs only in the (eta, beta) pair:
 *
 *  - **cme** (the default, and the only one `solver_mam_transient_qbd` uses):
 *    concentrated matrix exponential weights read from the vendored ILT-CME
 *    table. The entry chosen is the steepest -- smallest cv2 -- whose n+1 does
 *    not exceed the evaluation budget, which is exactly the reference's scan.
 *  - **euler**: binomial (Euler) weights, no table.
 *  - **gaver**: Gaver-Stehfest weights, no table.
 *
 * WHY cme IS NOT INTERCHANGEABLE WITH THE OTHER TWO, since it is tempting to
 * reach for a table-free variant: they are different quadratures with different
 * error behaviour, so substituting one for another changes the answer. The
 * reference defaults to cme; a port that quietly used euler would produce
 * numbers that are not the reference's while reporting the same method. That is
 * why the table is vendored rather than avoided.
 *
 * ARITHMETIC. Double only, and deliberately: the transform is evaluated at
 * COMPLEX arguments, and the whole route is transcendental. The signature takes
 * a `std::function` over `std::complex<double>` so a caller supplies its own
 * transform.
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <vector>

#include "line/api/mam/iltcme_table.h"
#include "line/util/error.h"

namespace line {
namespace mam {

/** Which Abate-Whitt weights to use. */
enum class IltMethod { Cme, Euler, Gaver };

/**
 * Invert a Laplace transform at the requested time points.
 *
 * @param fun        the transform F(s), evaluated at complex s
 * @param times      the points to invert at; each must be positive
 * @param maxFnEvals evaluation budget per point, which selects the CME entry
 * @param method     the weight family; the reference's default is Cme
 */
inline std::vector<double> matlab_ilt(
    const std::function<std::complex<double>(const std::complex<double>&)>& fun,
    const std::vector<double>& times, std::size_t maxFnEvals,
    IltMethod method = IltMethod::Cme) {
    if (times.empty()) return {};
    for (double t : times)
        if (!(t > 0.0))
            throw InputError("matlab_ilt: every evaluation time must be strictly positive");

    std::vector<std::complex<double>> eta, beta;

    if (method == IltMethod::Cme) {
        // The reference's scan: start at entry 0 REGARDLESS of its cost, then
        // take any entry that is steeper and affordable. Seeding with entry 0
        // before testing the bound is why it is always reachable.
        if (iltcme::kTableSize == 0)
            throw InputError("matlab_ilt: the vendored ILT-CME table is empty");
        const iltcme::CmeEntry* best = &iltcme::kTable[0];
        for (std::size_t i = 1; i < iltcme::kTableSize; ++i) {
            const iltcme::CmeEntry& e = iltcme::kTable[i];
            if (e.cv2 < best->cv2 && static_cast<std::size_t>(e.n) + 1 <= maxFnEvals) best = &e;
        }
        const double mu1 = best->mu1;
        eta.reserve(static_cast<std::size_t>(best->n) + 1);
        beta.reserve(static_cast<std::size_t>(best->n) + 1);
        eta.emplace_back(best->c * mu1, 0.0);
        beta.emplace_back(mu1, 0.0);
        for (int k = 0; k < best->n; ++k) {
            eta.emplace_back(best->a[k] * mu1, best->b[k] * mu1);
            beta.emplace_back(mu1, mu1 * static_cast<double>(k + 1) * best->omega);
        }
    } else {
        // The reference also carries table-free 'euler' and 'gaver' weights.
        // Neither is ported: `solver_mam_transient_qbd` is the only caller and
        // it takes the default, so those branches are unreachable here, and an
        // unreachable branch is an untested one. They are NOT a substitute for
        // cme in any case -- different quadratures with different error -- which
        // is exactly why the cme table was vendored rather than avoided.
        throw UnsupportedError(
            "matlab_ilt: only the CME weights are ported; the reference's 'euler' and 'gaver' "
            "variants are unreachable from this tree and are different quadratures, not "
            "substitutes for cme");
    }

    std::vector<double> out(times.size(), 0.0);
    for (std::size_t i = 0; i < times.size(); ++i) {
        const double t = times[i];
        double acc = 0.0;
        for (std::size_t k = 0; k < eta.size(); ++k)
            acc += std::real(eta[k] * fun(beta[k] / t));
        out[i] = acc / t;
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MATLAB_ILT_H
