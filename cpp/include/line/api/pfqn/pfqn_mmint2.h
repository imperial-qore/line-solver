/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_MMINT2_H
#define LINE_API_PFQN_PFQN_MMINT2_H

/**
 * McKenna-Mitra integral form of the normalizing constant of a repairman
 * model (one queueing station, R classes, per-class think time), in its three
 * MATLAB quadratures.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_mmint2.m,
 * pfqn_mmint2_gausslegendre.m and pfqn_mmint2_gausslaguerre.m. All three
 * evaluate
 *
 *   G = 1/prod_r N_r! * int_0^inf u^{m-1} e^{-u} prod_r (Z_r + L_r u)^{N_r} du
 *
 * and differ only in the rule: adaptive Gauss-Kronrod on a truncated interval,
 * a fixed Gauss-Legendre rule on [0, 1e6], or Gauss-Laguerre on [0, inf).
 * The integrand is a polynomial times e^{-u}, so Gauss-Laguerre with enough
 * nodes is exact up to rounding, which makes it the natural cross-check on the
 * other two and on pfqn_ca.
 *
 * NODES. MATLAB loads a Julia-generated table (gausslegendre-data.mat,
 * gausslaguerre-data.mat). A table cannot be carried across arithmetics -- it
 * would pin every instantiation to the precision it was generated at -- so the
 * rules are regenerated in T by pfqn_asympt_common.h.
 *
 * The Legendre form needs care. MATLAB's node count is
 * n = max(300, min(tablesize, 2(sum N + m - 1) - 1)) and it takes the FIRST n
 * entries of a 20000-point rule on [0, 1e6], which is not the same thing as a
 * fresh n-point rule: the 20000-point prefix spans [0.0036, 557.8] and
 * resolves the e^{-u} factor, whereas a genuine 300-point rule on [0, 1e6] has
 * its first node at u = 13.7 and misses the mass entirely (it returns
 * log G = -3.44 where the answer is 1.63). The port therefore generates the
 * 20000-point rule and takes the same prefix, computing only the prefix since
 * each Newton iteration is independent of the other nodes. `nodecap` is that
 * table length and defaults to 20000, the length of MATLAB's
 * gausslegendre-nodes.txt.
 *
 * TRUNCATION of the adaptive form. MATLAB integrates over
 * [0, -log(1 - (1 - 1e-12))] = [0, 27.63...], i.e. the 1 - 10^-order quantile
 * of the unit exponential, and asks for AbsTol 1e-12. The port keeps both
 * constants. That truncation is the dominant error for large populations,
 * where the polynomial factor pushes mass well beyond the cutoff; the tests
 * record where it starts to bite.
 *
 * ARITHMETIC. Quadrature, so all three are gated on
 * num_traits<T>::has_transcendental.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/api/qsys/qsys_quadrature.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/** Return value of the McKenna-Mitra quadratures, mirroring [G, lG]. */
template <class T>
struct MmintResult {
    T G;
    T lG;
};

/**
 * Adaptive form (MATLAB pfqn_mmint2): Gauss-Kronrod on [0, 27.63] with
 * absolute tolerance 1e-12.
 *
 * @param L (R) demand at the station, @param N (R) population,
 * @param Z (R) think times
 */
template <class T>
MmintResult<T> pfqn_mmint2(const std::vector<T>& L, const std::vector<T>& N,
                           const std::vector<T>& Z) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_mmint2 requires transcendental arithmetic (quadrature of e^{-u} p(u))");
    using std::exp;
    using std::log;
    const std::size_t R = L.size();
    if (N.size() != R || Z.size() != R)
        throw InputError("pfqn_mmint2: L, N and Z must have the same length");
    const T zero = num_traits<T>::from_int(0);

    // The reference restricts the product to the classes with N_r > 0.
    std::vector<std::size_t> nz;
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] != zero) nz.push_back(r);

    const auto f = [&](const T& u) {
        using std::exp;
        T p = exp(T(-u));
        for (std::size_t k = 0; k < nz.size(); ++k) {
            const std::size_t r = nz[k];
            const double nd = num_traits<T>::to_double(N[r]);
            p *= num_pow_int(T(Z[r] + L[r] * u), static_cast<unsigned>(nd));
        }
        return p;
    };

    const int order = 12;
    const T hi = num_traits<T>::from_double(-std::log(1.0 - (1.0 - std::pow(10.0, -order))));
    const T atol = num_traits<T>::from_double(std::pow(10.0, -order));
    const T I = qsys::detail::num_integral<T>(f, zero, hi, num_traits<T>::from_double(1e-12), atol);
    if (I <= zero) throw NumericError("pfqn_mmint2: non-positive integral");

    MmintResult<T> res;
    T lG = log(I);
    for (std::size_t r = 0; r < R; ++r) lG -= detail::num_factln<T>(N[r]);
    res.lG = lG;
    res.G = exp(lG);
    return res;
}

/**
 * Gauss-Legendre form on [0, 1e6] (MATLAB pfqn_mmint2_gausslegendre).
 *
 * @param m station multiplicity, contributing the u^{m-1} factor
 * @param nodecap size of the underlying tabulated rule, i.e. MATLAB's table
 *                length; the routine uses its first n nodes
 * @param L (M) service demands
 * @param N (1) population, single class
 * @param Z (1) think time
 */
template <class T>
MmintResult<T> pfqn_mmint2_gausslegendre(const std::vector<T>& L, const std::vector<T>& N,
                                         const std::vector<T>& Z, int m, std::size_t nodecap) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_mmint2_gausslegendre requires transcendental arithmetic (quadrature)");
    using std::exp;
    using std::log;
    const std::size_t R = L.size();
    if (N.size() != R || Z.size() != R)
        throw InputError("pfqn_mmint2_gausslegendre: L, N and Z must have the same length");
    if (m < 1) throw InputError("pfqn_mmint2_gausslegendre: multiplicity must be at least one");
    const T zero = num_traits<T>::from_int(0);
    T Ntot = zero;
    for (const T& v : N) Ntot += v;

    const long want = 2 * (static_cast<long>(num_traits<T>::to_double(Ntot)) + m - 1) - 1;
    std::size_t n = 300;
    const std::size_t capped = std::min<std::size_t>(nodecap, want > 0 ? static_cast<std::size_t>(want) : 1);
    if (capped > n) n = capped;
    if (n > nodecap) n = nodecap;

    // 20000-point table prefix rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    std::vector<T> x, w;
    detail::gauss_legendre<T>(nodecap, zero, num_traits<T>::from_double(1e6), x, w, n);

    std::vector<T> g(n);
    for (std::size_t i = 0; i < n; ++i) {
        T y = zero;
        for (std::size_t r = 0; r < R; ++r) {
            if (N[r] == zero) continue;
            y += N[r] * log(T(Z[r] + L[r] * x[i]));
        }
        g[i] = T(log(w[i]) - x[i] + y);
        if (m > 1) g[i] += num_traits<T>::from_int(m - 1) * log(x[i]);
    }
    T coeff = zero;
    for (std::size_t r = 0; r < R; ++r) coeff -= detail::num_factln<T>(N[r]);
    coeff -= detail::num_factln<T>(num_traits<T>::from_int(m - 1));

    MmintResult<T> res;
    // stable logsumexp rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    res.lG = T(detail::logsumexp(g) + coeff);
    res.G = exp(res.lG);
    return res;
}

template <class T>
MmintResult<T> pfqn_mmint2_gausslegendre(const std::vector<T>& L, const std::vector<T>& N,
                                         const std::vector<T>& Z) {
    // 20000 is the length of MATLAB's gausslegendre-nodes.txt.
    return pfqn_mmint2_gausslegendre(L, N, Z, 1, 20000);
}

/**
 * Gauss-Laguerre form (MATLAB pfqn_mmint2_gausslaguerre).
 *
 * @param npts node count; MATLAB uses the length of its tabulated rule
 * @param L (M) service demands
 * @param N (1) population, single class
 * @param Z (1) think time
 * @param m multiplicity of the queueing station
 */
template <class T>
MmintResult<T> pfqn_mmint2_gausslaguerre(const std::vector<T>& L, const std::vector<T>& N,
                                         const std::vector<T>& Z, int m, std::size_t npts) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_mmint2_gausslaguerre requires transcendental arithmetic (quadrature)");
    using std::exp;
    using std::log;
    const std::size_t R = L.size();
    if (N.size() != R || Z.size() != R)
        throw InputError("pfqn_mmint2_gausslaguerre: L, N and Z must have the same length");
    if (m < 1) throw InputError("pfqn_mmint2_gausslaguerre: multiplicity must be at least one");
    if (npts < 2) throw InputError("pfqn_mmint2_gausslaguerre: at least two nodes are required");
    const T zero = num_traits<T>::from_int(0);

    std::vector<T> x, w;
    detail::gauss_laguerre<T>(npts, x, w);
    std::vector<T> g(npts);
    for (std::size_t i = 0; i < npts; ++i) {
        T F = zero;
        if (m > 1) F += num_traits<T>::from_int(m - 1) * log(x[i]);
        for (std::size_t r = 0; r < R; ++r) {
            if (N[r] == zero) continue;
            F += N[r] * log(T(Z[r] + L[r] * x[i]));
        }
        g[i] = T(log(w[i]) + F);
    }
    T coeff = zero;
    for (std::size_t r = 0; r < R; ++r) coeff -= detail::num_factln<T>(N[r]);
    coeff -= detail::num_factln<T>(num_traits<T>::from_int(m - 1));

    MmintResult<T> res;
    res.lG = T(detail::logsumexp(g) + coeff);
    res.G = exp(res.lG);
    return res;
}

template <class T>
MmintResult<T> pfqn_mmint2_gausslaguerre(const std::vector<T>& L, const std::vector<T>& N,
                                         const std::vector<T>& Z) {
    return pfqn_mmint2_gausslaguerre(L, N, Z, 1, 90);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_MMINT2_H
