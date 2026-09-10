/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_NC_SANITIZE_H
#define LINE_API_PFQN_NC_SANITIZE_H

/**
 * Preprocessing shared by the normalizing-constant solvers: drop the classes
 * that cannot contribute, rescale the demands per class, and order the classes
 * so that the zero-think-time ones come first.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_nc_sanitize.m.
 *
 * The transformation is a change of variables on G, not an approximation:
 * dividing every demand and think time of class r by Lmax_r divides G(N) by
 * exactly Lmax_r^{N_r}, and removing a class whose demands are all zero
 * factors out the delay term Z_r^{N_r}/N_r!. Both are returned, so the caller
 * recovers G(N) of the original model as
 *
 *   G_original(N) = Gremaind * G_sanitized(N_sanitized).
 *
 * MATLAB returns only the LOGARITHM lGremaind of that factor, which is what
 * forces every downstream CoMoM routine into log space. Here the factor itself
 * is returned as a value of T, so the rational path never leaves the field,
 * and lGremaind is provided alongside for the double callers.
 *
 * Two divergences from the reference, both deliberate.
 *
 * 1. REFERENCE DEFECT, corrected here. MATLAB's zero-demand branch reads
 *
 *      lGremaind = lGremaind + N(zerodemands)*log(Z(zerodemands))' ...
 *                            - sum(log(N(zerodemands)));
 *
 *    The delay balance function of a class confined to the think-time node is
 *    Z_r^{N_r}/N_r!, whose log is N_r log Z_r - log(N_r!), i.e. -factln(N_r).
 *    The reference subtracted log(N_r) instead of log(N_r!), which agrees only
 *    at N_r = 1 and 2 (1! = 1, 2! = 2) and diverges like log((N_r-1)!) after
 *    that. Measured on L = [0 0.5], N = [3 2], Z = [1 0]: MATLAB reconstructed
 *    lG = -2.484906649788 against the exact -3.178053830348, short by exactly
 *    log 2. FIXED IN MATLAB (commit "m fix: pfqn_nc_sanitize delay term uses
 *    factln and column selection"), which now reconstructs the exact constant
 *    to 0. `Pfqn_nc_sanitize.java` still carries the old form.
 *
 * 2. MATLAB indexes the zero-demand and zero-think-time tests with find() on a
 *    MATRIX, e.g. `zerodemands = find(L < atol)`, and then uses the resulting
 *    LINEAR indices as COLUMN indices. That is only equivalent to the intended
 *    per-class test when L has a single row, which is exactly the case for the
 *    repairman callers this routine was written for. This port applies the test
 *    to the column sums; MATLAB now does the same, and the M = 2 case
 *    reconstructs the exact constant to 0 where it previously could not.
 *
 * Arithmetic: EXACT-CAPABLE. Every operation is a comparison, a division or a
 * multiplication in the field of the inputs. Pass atol = 0 at T = Rational to
 * get the exact "is identically zero" tests; a positive atol reproduces the
 * reference's tolerant filtering in any arithmetic.
 */

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

template <class T>
struct NcSanitizeResult {
    std::vector<T> lambda;  ///< retained arrival rates, in the new class order
    Matrix<T> L;            ///< retained demands, rescaled and reordered
    std::vector<int> N;     ///< retained populations, reordered
    Matrix<T> Z;            ///< retained think times, rescaled and reordered
    T Gremaind;             ///< multiplicative factor removed from G
    double lGremaind;       ///< its logarithm, for the log-space callers
    std::vector<std::size_t> classIndex;  ///< original index of each retained class
};

/**
 * @param lambda (R) arrival rates; may be empty for a purely closed model
 * @param L      (M x R) service demands
 * @param N      (R) populations
 * @param Z      (K x R) think times; may be empty
 * @param atol   threshold below which a demand counts as zero; use 0 for exact
 */
template <class T>
NcSanitizeResult<T> pfqn_nc_sanitize(const std::vector<T>& lambda, const Matrix<T>& L,
                                     const std::vector<int>& N, const Matrix<T>& Z, const T& atol) {
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_nc_sanitize: L and N disagree on the class count");
    if (!Z.empty() && Z.cols() != R)
        throw InputError("pfqn_nc_sanitize: Z and N disagree on the class count");
    if (!lambda.empty() && lambda.size() != R)
        throw InputError("pfqn_nc_sanitize: lambda and N disagree on the class count");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t M = L.empty() ? 0 : L.rows();
    const std::size_t D = Z.empty() ? 0 : Z.rows();

    // Column sums, the per-class aggregates every test below is phrased on.
    const auto colsum = [](const Matrix<T>& A, std::size_t r, const T& z) {
        T s = z;
        for (std::size_t i = 0; i < A.rows(); ++i) s += A(i, r);
        return s;
    };

    NcSanitizeResult<T> res;
    res.Gremaind = one;

    // ---- keep only the classes that have jobs and a well-defined demand ----
    std::vector<std::size_t> keep;
    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] == 0) continue;                        // empty class: contributes 1
        const T tot = colsum(L, r, zero) + colsum(Z, r, zero);
        if (tot < atol) continue;                       // ill-defined class
        keep.push_back(r);
    }

    // ---- classes with no queueing demand at all live entirely in the delay ----
    // Their contribution factors out exactly as Z_r^{N_r} / N_r!.
    std::vector<std::size_t> retained;
    for (std::size_t k = 0; k < keep.size(); ++k) {
        const std::size_t r = keep[k];
        const T ldem = colsum(L, r, zero);
        if (M > 0 && ldem < atol) {
            const T zr = colsum(Z, r, zero);
            res.Gremaind *= num_pow_int(zr, static_cast<unsigned>(N[r])) /
                            num_factorial<T>(static_cast<unsigned>(N[r]));
            continue;
        }
        retained.push_back(r);
    }

    const std::size_t Rk = retained.size();

    // per-class rescaling rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    std::vector<T> scale(Rk, one);
    for (std::size_t k = 0; k < Rk; ++k) {
        const std::size_t r = retained[k];
        if (M == 0) continue;
        T mx = L(0, r);
        for (std::size_t i = 1; i < M; ++i)
            if (L(i, r) > mx) mx = L(i, r);
        if (mx > zero) scale[k] = mx;
    }
    for (std::size_t k = 0; k < Rk; ++k)
        res.Gremaind *= num_pow_int(scale[k], static_cast<unsigned>(N[retained[k]]));

    // ---- order: ascending total think time, then zero-think-time first ----
    std::vector<std::size_t> ord(Rk);
    std::iota(ord.begin(), ord.end(), static_cast<std::size_t>(0));
    std::vector<T> zsum(Rk, zero);
    for (std::size_t k = 0; k < Rk; ++k) zsum[k] = colsum(Z, retained[k], zero) / scale[k];
    // Stable, so classes with equal think time keep their relative order, as
    // MATLAB's sort does.
    std::stable_sort(ord.begin(), ord.end(),
                     [&](std::size_t a, std::size_t b) { return zsum[a] < zsum[b]; });
    // zero-think-time reordering rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    std::stable_partition(ord.begin(), ord.end(),
                          [&](std::size_t a) { return !(zsum[a] >= atol); });

    res.L = Matrix<T>(M, Rk);
    res.Z = Matrix<T>(D, Rk);
    res.N.assign(Rk, 0);
    res.classIndex.assign(Rk, 0);
    if (!lambda.empty()) res.lambda.assign(Rk, zero);
    for (std::size_t k = 0; k < Rk; ++k) {
        const std::size_t src = retained[ord[k]];
        const T& sc = scale[ord[k]];
        for (std::size_t i = 0; i < M; ++i) res.L(i, k) = L(i, src) / sc;
        for (std::size_t i = 0; i < D; ++i) res.Z(i, k) = Z(i, src) / sc;
        res.N[k] = N[src];
        res.classIndex[k] = src;
        if (!lambda.empty()) res.lambda[k] = lambda[src];
    }

    res.lGremaind = num_traits<T>::log_as_double(res.Gremaind);
    return res;
}

/** Overload with the exact (zero-tolerance) tests. */
template <class T>
NcSanitizeResult<T> pfqn_nc_sanitize(const Matrix<T>& L, const std::vector<int>& N,
                                     const Matrix<T>& Z) {
    return pfqn_nc_sanitize(std::vector<T>(), L, N, Z, num_traits<T>::from_int(0));
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_NC_SANITIZE_H
