/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_JOINT_H
#define LINE_API_MOMENT_MOMENT_JOINT_H

/**
 * Joint moment conversions on the house of moments.
 *
 * Templated port of the matlab/src/api/moment/moment_joint_*.m family. Every
 * edge except the cumulant and the central ones is SEPARABLE: the joint table
 * is the Kronecker product of the univariate tables, so the conversion is one
 * mode product per dimension. The cumulant edges are not separable and carry
 * their own multivariate recurrence, and the central edges are separable only
 * once the mean vector is fixed.
 *
 * The joint cumulant recurrence picks j as the FIRST dimension with a positive
 * order and differentiates in that variable. Any other choice gives the same
 * answer, but not the same summation, so a port that picks the last dimension
 * will disagree in floating point even when it is algebraically right.
 *
 * Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 */

#include <algorithm>
#include <vector>

#include "line/api/moment/moment_housematrix.h"
#include "line/api/moment/moment_tensor.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace moment {

/** Joint binomial moments from joint factorial moments. */
template <class T>
MomentTensor<T> moment_joint_binomial_from_factorial(const MomentTensor<T>& f) {
    return moment_jointtrans<T>(f, MomentEdge::BinomialFromFactorial);
}

/** Joint binomial moments from joint negative-binomial moments. */
template <class T>
MomentTensor<T> moment_joint_binomial_from_negbinomial(const MomentTensor<T>& bm) {
    return moment_jointtrans<T>(bm, MomentEdge::BinomialFromNegbinomial);
}

/** Joint binomial moments from joint tail moments. */
template <class T>
MomentTensor<T> moment_joint_binomial_from_tail(const MomentTensor<T>& t) {
    return moment_jointtrans<T>(t, MomentEdge::BinomialFromTail);
}

/** Joint factorial moments from joint binomial moments. */
template <class T>
MomentTensor<T> moment_joint_factorial_from_binomial(const MomentTensor<T>& b) {
    return moment_jointtrans<T>(b, MomentEdge::FactorialFromBinomial);
}

/** Joint factorial moments from joint raw moments. */
template <class T>
MomentTensor<T> moment_joint_factorial_from_raw(const MomentTensor<T>& m) {
    return moment_jointtrans<T>(m, MomentEdge::FactorialFromRaw);
}

/** Joint factorial moments from joint upper-factorial moments. */
template <class T>
MomentTensor<T> moment_joint_factorial_from_upfactorial(const MomentTensor<T>& fp) {
    return moment_jointtrans<T>(fp, MomentEdge::FactorialFromUpfactorial);
}

/** Joint negative-binomial moments from joint binomial moments. */
template <class T>
MomentTensor<T> moment_joint_negbinomial_from_binomial(const MomentTensor<T>& b) {
    return moment_jointtrans<T>(b, MomentEdge::NegbinomialFromBinomial);
}

/** Joint negative-binomial moments from joint upper-factorial moments. */
template <class T>
MomentTensor<T> moment_joint_negbinomial_from_upfactorial(const MomentTensor<T>& fp) {
    return moment_jointtrans<T>(fp, MomentEdge::NegbinomialFromUpfactorial);
}

/** Joint raw moments from joint factorial moments. */
template <class T>
MomentTensor<T> moment_joint_raw_from_factorial(const MomentTensor<T>& f) {
    return moment_jointtrans<T>(f, MomentEdge::RawFromFactorial);
}

/** Joint raw moments from joint upper-factorial moments. */
template <class T>
MomentTensor<T> moment_joint_raw_from_upfactorial(const MomentTensor<T>& fp) {
    return moment_jointtrans<T>(fp, MomentEdge::RawFromUpfactorial);
}

/** Joint tail moments from joint binomial moments. */
template <class T>
MomentTensor<T> moment_joint_tail_from_binomial(const MomentTensor<T>& b) {
    return moment_jointtrans<T>(b, MomentEdge::TailFromBinomial);
}

/** Joint upper-factorial moments from joint factorial moments. */
template <class T>
MomentTensor<T> moment_joint_upfactorial_from_factorial(const MomentTensor<T>& f) {
    return moment_jointtrans<T>(f, MomentEdge::UpfactorialFromFactorial);
}

/** Joint upper-factorial moments from joint negative-binomial moments. */
template <class T>
MomentTensor<T> moment_joint_upfactorial_from_negbinomial(const MomentTensor<T>& bm) {
    return moment_jointtrans<T>(bm, MomentEdge::UpfactorialFromNegbinomial);
}

/** Joint upper-factorial moments from joint raw moments. */
template <class T>
MomentTensor<T> moment_joint_upfactorial_from_raw(const MomentTensor<T>& m) {
    return moment_jointtrans<T>(m, MomentEdge::UpfactorialFromRaw);
}

namespace detail {

/** Advances a column-major multi-index within the box sz, MATLAB odometer order. */
inline void moment_odometer(std::vector<std::size_t>& ord, const std::vector<std::size_t>& sz) {
    for (std::size_t l = 0; l < ord.size(); ++l) {
        ++ord[l];
        if (ord[l] < sz[l]) return;
        ord[l] = 0;
    }
}

/**
 * Shared body of the joint cumulant recurrence. With include_full = false the
 * term b = a is skipped, which turns the raw-from-cumulant sum into the
 * cumulant-from-raw correction.
 */
template <class T>
std::vector<T> moment_joint_cumulant_body(const std::vector<std::size_t>& sz,
                                          const std::vector<T>& src, bool cumulant_from_raw) {
    const std::size_t d = sz.size();
    std::size_t nel = 1;
    for (std::size_t l = 0; l < d; ++l) nel *= sz[l];
    std::vector<std::size_t> stride(d, 1);
    for (std::size_t l = 1; l < d; ++l) stride[l] = stride[l - 1] * sz[l - 1];
    std::vector<T> dst(nel, num_traits<T>::from_int(0));
    if (!cumulant_from_raw) dst[0] = num_traits<T>::from_int(1);
    std::vector<std::size_t> ord(d, 0), bord(d, 0);
    for (std::size_t ia = 0; ia < nel; ++ia) {
        std::size_t j = d;
        for (std::size_t l = 0; l < d; ++l)
            if (ord[l] > 0) {
                j = l;
                break;
            }
        if (j < d) {
            T acc = num_traits<T>::from_int(0);
            std::size_t nb = 1;
            for (std::size_t l = 0; l < d; ++l) nb *= ord[l] + 1;
            std::fill(bord.begin(), bord.end(), static_cast<std::size_t>(0));
            for (std::size_t ib = 0; ib < nb; ++ib) {
                bool any_pos = false, equal_a = true;
                for (std::size_t l = 0; l < d; ++l) {
                    if (bord[l] > 0) any_pos = true;
                    if (bord[l] != ord[l]) equal_a = false;
                }
                const bool keep =
                    any_pos && bord[j] > 0 && (!cumulant_from_raw || !equal_a);
                if (keep) {
                    T c = num_traits<T>::from_int(1);
                    for (std::size_t l = 0; l < d; ++l) {
                        const int alpha = static_cast<int>(ord[l]) - (l == j ? 1 : 0);
                        const int beta = static_cast<int>(bord[l]) - (l == j ? 1 : 0);
                        c *= num_nck<T>(alpha, beta);
                    }
                    std::size_t ib_lin = 0, ic_lin = 0;
                    for (std::size_t l = 0; l < d; ++l) {
                        ib_lin += bord[l] * stride[l];
                        ic_lin += (ord[l] - bord[l]) * stride[l];
                    }
                    // the cumulant factor is always the one indexed by b
                    acc += cumulant_from_raw ? c * dst[ib_lin] * src[ic_lin]
                                             : c * src[ib_lin] * dst[ic_lin];
                }
                // the inner odometer runs over the box 0..ord, not over sz
                for (std::size_t l = 0; l < d; ++l) {
                    ++bord[l];
                    if (bord[l] <= ord[l]) break;
                    bord[l] = 0;
                }
            }
            dst[ia] = cumulant_from_raw ? src[ia] - acc : acc;
        }
        moment_odometer(ord, sz);
    }
    return dst;
}

}  // namespace detail

/** Joint cumulants from joint raw moments. */
template <class T>
MomentTensor<T> moment_joint_cumulant_from_raw(const MomentTensor<T>& m) {
    MomentTensor<T> kappa(m.sz);
    kappa.data = detail::moment_joint_cumulant_body<T>(m.sz, m.data, true);
    return kappa;
}

/** Joint raw moments from joint cumulants. */
template <class T>
MomentTensor<T> moment_joint_raw_from_cumulant(const MomentTensor<T>& kappa) {
    MomentTensor<T> m(kappa.sz);
    m.data = detail::moment_joint_cumulant_body<T>(kappa.sz, kappa.data, false);
    return m;
}

/** Joint factorial cumulants from joint factorial moments. */
template <class T>
MomentTensor<T> moment_joint_factcumulant_from_factorial(const MomentTensor<T>& f) {
    return moment_joint_cumulant_from_raw<T>(f);
}

/** Joint factorial moments from joint factorial cumulants. */
template <class T>
MomentTensor<T> moment_joint_factorial_from_factcumulant(const MomentTensor<T>& kappa) {
    return moment_joint_raw_from_cumulant<T>(kappa);
}

/** Joint central moments from joint raw moments about a given mean vector. */
template <class T>
MomentTensor<T> moment_joint_central_from_raw_mean(const MomentTensor<T>& m,
                                                   const std::vector<T>& mu) {
    const std::size_t d = m.order();
    if (mu.size() != d)
        throw InputError(
            "moment_joint_central_from_raw_mean: mu must have one entry per dimension of m");
    MomentTensor<T> mc = m;
    for (std::size_t mode = 0; mode < d; ++mode) {
        const int n = static_cast<int>(m.sz[mode]) - 1;
        Matrix<T> Tm(static_cast<std::size_t>(n) + 1, static_cast<std::size_t>(n) + 1,
                     num_traits<T>::from_int(0));
        const T negmu = -mu[mode];
        for (int i = 0; i <= n; ++i)
            for (int k = 0; k <= i; ++k)
                Tm(i, k) = num_nck<T>(i, k) * num_pow_int(negmu, static_cast<unsigned>(i - k));
        mc = moment_tensortrans<T>(mc, Tm, mode);
    }
    return mc;
}

/** Joint central moments from joint raw moments, reading the means off m. */
template <class T>
MomentTensor<T> moment_joint_central_from_raw(const MomentTensor<T>& m) {
    const std::size_t d = m.order();
    for (std::size_t l = 0; l < d; ++l)
        if (m.sz[l] < 2)
            throw InputError(
                "moment_joint_central_from_raw: the means m_(e_j) are required, hence every "
                "dimension of m must have at least 2 elements");
    std::vector<T> mu(d);
    for (std::size_t l = 0; l < d; ++l) mu[l] = m.data[m.stride(l)];
    return moment_joint_central_from_raw_mean<T>(m, mu);
}

/** Joint raw moments from joint central moments about a given mean vector. */
template <class T>
MomentTensor<T> moment_joint_raw_from_central(const MomentTensor<T>& mc,
                                              const std::vector<T>& mu) {
    const std::size_t d = mc.order();
    if (mu.size() != d)
        throw InputError(
            "moment_joint_raw_from_central: mu must have one entry per dimension of mc");
    MomentTensor<T> m = mc;
    for (std::size_t mode = 0; mode < d; ++mode) {
        const int n = static_cast<int>(mc.sz[mode]) - 1;
        Matrix<T> Tm(static_cast<std::size_t>(n) + 1, static_cast<std::size_t>(n) + 1,
                     num_traits<T>::from_int(0));
        for (int i = 0; i <= n; ++i)
            for (int k = 0; k <= i; ++k)
                Tm(i, k) = num_nck<T>(i, k) * num_pow_int(mu[mode], static_cast<unsigned>(i - k));
        m = moment_tensortrans<T>(m, Tm, mode);
    }
    return m;
}

/** Joint central moments from joint tail moments. */
template <class T>
MomentTensor<T> moment_joint_central_from_tail(const MomentTensor<T>& t) {
    const MomentTensor<T> b = moment_joint_binomial_from_tail<T>(t);
    const MomentTensor<T> f = moment_joint_factorial_from_binomial<T>(b);
    return moment_joint_central_from_raw<T>(moment_joint_raw_from_factorial<T>(f));
}

/** Factorial moments of a total count from the joint factorial moments of its parts. */
template <class T>
std::vector<T> moment_joint_aggregate(const MomentTensor<T>& F) {
    const std::size_t d = F.order();
    std::size_t nmax = F.sz[0];
    for (std::size_t l = 1; l < d; ++l) nmax = F.sz[l] < nmax ? F.sz[l] : nmax;
    if (nmax == 0) throw InputError("moment_joint_aggregate: F must be nonempty");
    --nmax;
    std::vector<T> f(nmax + 1, num_traits<T>::from_int(0));
    std::vector<std::size_t> ord(d, 0);
    const std::size_t nel = F.numel();
    for (std::size_t ia = 0; ia < nel; ++ia) {
        std::size_t n = 0;
        for (std::size_t l = 0; l < d; ++l) n += ord[l];
        if (n <= nmax) {
            T c = num_factorial<T>(static_cast<unsigned>(n));
            for (std::size_t l = 0; l < d; ++l)
                c = c / num_factorial<T>(static_cast<unsigned>(ord[l]));
            f[n] += c * F.data[ia];
        }
        detail::moment_odometer(ord, F.sz);
    }
    return f;
}

/** Joint factorial moments of a multinomially marked count from the aggregate ones. */
template <class T>
MomentTensor<T> moment_joint_marking(const std::vector<T>& f, const std::vector<T>& p,
                                     const std::vector<std::size_t>& dims) {
    const std::size_t d = p.size();
    if (dims.size() != d)
        throw InputError("moment_joint_marking: p and dims must have the same length");
    std::size_t total = 0;
    for (std::size_t l = 0; l < d; ++l) total += dims[l];
    if (f.empty() || total > f.size() - 1)
        throw InputError(
            "moment_joint_marking: the aggregate factorial moments must reach order sum(dims)");
    std::vector<std::size_t> sz(d);
    for (std::size_t l = 0; l < d; ++l) sz[l] = dims[l] + 1;
    MomentTensor<T> F(sz);
    std::vector<std::size_t> ord(d, 0);
    const std::size_t nel = F.numel();
    for (std::size_t ia = 0; ia < nel; ++ia) {
        std::size_t n = 0;
        T c = num_traits<T>::from_int(1);
        for (std::size_t l = 0; l < d; ++l) {
            c *= num_pow_int(p[l], static_cast<unsigned>(ord[l]));
            n += ord[l];
        }
        F.data[ia] = c * f[n];
        detail::moment_odometer(ord, sz);
    }
    return F;
}

}  // namespace moment
}  // namespace line

#endif
