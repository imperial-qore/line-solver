/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMAP_STATS_H
#define LINE_API_MAM_MMAP_STATS_H

/**
 * Marked MAP statistics: embedded chains, class-transition probabilities,
 * forward and cross moments, counting means and covariances.
 *
 * Templated port of the M3A MMAP statistics in matlab/lib/m3a/m3a/mmap:
 * mmap_pie.m, mmap_embedded.m, mmap_maps.m, mmap_timereverse.m, mmap_sigma.m,
 * mmap_sigma2.m, mmap_count_mean.m, mmap_count_idc.m, mmap_count_mcov.m,
 * mmap_idc.m, mmap_cross_moment.m, mmap_forward_moment.m and mmap_sum.m.
 *
 * The embedded per-class matrix is E_c = (-D0)^-1 D1^(c). It is SUBstochastic,
 * not stochastic: its row sums give the probability that the next arrival is of
 * class c, which is what makes sum_c E_c the embedded chain of the aggregate
 * MAP and E_c on its own the class-c defective kernel. Reading E_c as a
 * transition matrix and renormalizing it is the classic way to get the
 * per-class moments wrong.
 *
 * mmap_issym and mmap_shorten have no C++ counterpart: the first asks whether
 * the MATLAB cell holds symbolic entries, which here is the template parameter,
 * and the second reorders a MATLAB cell into BUTools order, which the Mmap
 * struct already encodes by construction.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_count_var.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace detail {

/** E_c = (-D0)^-1 D1^(c), the defective embedded kernel of class c. */
template <class T>
Matrix<T> mmap_embedded_class(const Mmap<T>& mm, std::size_t c) {
    Matrix<T> negD0 = mm.D0;
    for (std::size_t i = 0; i < negD0.rows(); ++i)
        for (std::size_t j = 0; j < negD0.cols(); ++j) negD0(i, j) = -negD0(i, j);
    return matmul(inverse(negD0), mm.Dc[c]);
}

}  // namespace detail

/** Embedded per-class kernels E_c = (-D0)^-1 D1^(c). */
template <class T>
std::vector<Matrix<T>> mmap_embedded(const Mmap<T>& mm) {
    std::vector<Matrix<T>> Pc;
    Pc.reserve(mm.classes());
    for (std::size_t c = 0; c < mm.classes(); ++c) Pc.push_back(detail::mmap_embedded_class(mm, c));
    return Pc;
}

/** The C MAPs seen by each class, MAP_c = (D0 + D1 - D1^(c), D1^(c)). */
template <class T>
std::vector<Map<T>> mmap_maps(const Mmap<T>& mm) {
    std::vector<Map<T>> maps;
    maps.reserve(mm.classes());
    for (std::size_t c = 0; c < mm.classes(); ++c) {
        Matrix<T> A = mm.D0;
        for (std::size_t i = 0; i < A.rows(); ++i)
            for (std::size_t j = 0; j < A.cols(); ++j) A(i, j) += mm.D1(i, j) - mm.Dc[c](i, j);
        maps.push_back(Map<T>{A, mm.Dc[c]});
    }
    return maps;
}

/**
 * Stationary phase distribution seen just after a class-c arrival, one row per
 * class. The row is the left invariant vector of the STOCHASTIC matrix
 * P_c = (-D0 - D1 + D1^(c))^-1 D1^(c), which is the chain watched only at
 * class-c epochs, and is not E_c.
 */
template <class T>
Matrix<T> mmap_pie(const Mmap<T>& mm) {
    const std::size_t n = mm.order();
    const std::size_t C = mm.classes();
    Matrix<T> pie(C, n, num_traits<T>::from_int(0));
    for (std::size_t c = 0; c < C; ++c) {
        Matrix<T> B(n, n, num_traits<T>::from_int(0));
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                B(i, j) = -mm.D0(i, j) - mm.D1(i, j) + mm.Dc[c](i, j);
        const Matrix<T> Pc = matmul(inverse(B), mm.Dc[c]);
        Matrix<T> A(n, n, num_traits<T>::from_int(0));
        for (std::size_t i = 0; i + 1 < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                A(i, j) = Pc(j, i) - (i == j ? num_traits<T>::from_int(1)
                                             : num_traits<T>::from_int(0));
        for (std::size_t j = 0; j < n; ++j) A(n - 1, j) = num_traits<T>::from_int(1);
        std::vector<T> b(n, num_traits<T>::from_int(0));
        b[n - 1] = num_traits<T>::from_int(1);
        const std::vector<T> x = solve(A, b);
        for (std::size_t j = 0; j < n; ++j) pie(c, j) = x[j];
    }
    return pie;
}

/** Time-reversed MMAP, D^-1 M' D with D = diag(map_prob) applied to every matrix. */
template <class T>
Mmap<T> mmap_timereverse(const Mmap<T>& mm) {
    const std::size_t n = mm.order();
    const std::vector<T> piq = map_prob(mm.map());
    Mmap<T> out;
    out.D0 = Matrix<T>(n, n, num_traits<T>::from_int(0));
    out.D1 = Matrix<T>(n, n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            out.D0(i, j) = mm.D0(j, i) * piq[j] / piq[i];
            out.D1(i, j) = mm.D1(j, i) * piq[j] / piq[i];
        }
    out.Dc.reserve(mm.classes());
    for (std::size_t c = 0; c < mm.classes(); ++c) {
        Matrix<T> R(n, n, num_traits<T>::from_int(0));
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) R(i, j) = mm.Dc[c](j, i) * piq[j] / piq[i];
        out.Dc.push_back(R);
    }
    return out;
}

/** sigma(i,j) = pie E_i E_j 1, the probability that two consecutive marks are (i,j). */
template <class T>
Matrix<T> mmap_sigma(const Mmap<T>& mm) {
    const std::size_t C = mm.classes();
    const std::vector<T> alpha = map_pie(mm.map());
    const std::vector<Matrix<T>> E = mmap_embedded(mm);
    Matrix<T> sigma(C, C, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < C; ++i) {
        const std::vector<T> start = vecmul(alpha, E[i]);
        for (std::size_t j = 0; j < C; ++j) {
            const std::vector<T> v = vecmul(start, E[j]);
            T acc = num_traits<T>::from_int(0);
            for (std::size_t k = 0; k < v.size(); ++k) acc += v[k];
            sigma(i, j) = acc;
        }
    }
    return sigma;
}

/** sigma2(i,j,h) = pie E_i E_j E_h 1, indexed as sigma2[i][j][h]. */
template <class T>
std::vector<std::vector<std::vector<T>>> mmap_sigma2(const Mmap<T>& mm) {
    const std::size_t C = mm.classes();
    const std::vector<T> alpha = map_pie(mm.map());
    const std::vector<Matrix<T>> E = mmap_embedded(mm);
    std::vector<std::vector<std::vector<T>>> sigma(
        C, std::vector<std::vector<T>>(C, std::vector<T>(C, num_traits<T>::from_int(0))));
    for (std::size_t i = 0; i < C; ++i) {
        const std::vector<T> starti = vecmul(alpha, E[i]);
        for (std::size_t j = 0; j < C; ++j) {
            const std::vector<T> startj = vecmul(starti, E[j]);
            for (std::size_t h = 0; h < C; ++h) {
                const std::vector<T> v = vecmul(startj, E[h]);
                T acc = num_traits<T>::from_int(0);
                for (std::size_t k = 0; k < v.size(); ++k) acc += v[k];
                sigma[i][j][h] = acc;
            }
        }
    }
    return sigma;
}

/** Per-class mean of the counting process over a window of length t. */
template <class T>
std::vector<T> mmap_count_mean(const Mmap<T>& mm, const T& t) {
    const std::size_t n = mm.order();
    const std::size_t C = mm.classes();
    const std::vector<T> theta = map_prob(mm.map());
    std::vector<T> mk(C, num_traits<T>::from_int(0));
    for (std::size_t c = 0; c < C; ++c) {
        T acc = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) acc += theta[i] * mm.Dc[c](i, j);
        mk[c] = acc * t;
    }
    return mk;
}

/** Per-class index of dispersion of counts over a window of length t. */
template <class T>
std::vector<T> mmap_count_idc(const Mmap<T>& mm, const T& t) {
    const std::vector<T> m = mmap_count_mean(mm, t);
    const std::vector<T> v = mmap_count_var(mm, t);
    std::vector<T> idc(m.size(), num_traits<T>::from_int(0));
    for (std::size_t c = 0; c < m.size(); ++c) idc[c] = v[c] / m[c];
    return idc;
}

/** Asymptotic per-class index of dispersion, evaluated at t = 1e6 / sum_c lambda_c. */
template <class T>
std::vector<T> mmap_idc(const Mmap<T>& mm) {
    const std::vector<T> lam = mmap_lambda(mm);
    T total = num_traits<T>::from_int(0);
    for (std::size_t c = 0; c < lam.size(); ++c) total += lam[c];
    const T tinf = num_traits<T>::from_int(1000000) / total;
    return mmap_count_idc(mm, tinf);
}

/**
 * Covariance matrix of the per-class counts over a window of length t. The
 * off-diagonal entries come from the polarization identity on the pooled
 * classes, since only the variance of a single mark is available in closed form.
 */
template <class T>
Matrix<T> mmap_count_mcov(const Mmap<T>& mm, const T& t) {
    const std::size_t C = mm.classes();
    const std::vector<T> mV = mmap_count_var(mm, t);
    Matrix<T> S(C, C, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < C; ++i) S(i, i) = mV[i];
    const T half = num_traits<T>::from_int(1) / num_traits<T>::from_int(2);
    for (std::size_t i = 0; i < C; ++i)
        for (std::size_t j = 0; j < C; ++j) {
            if (i == j) continue;
            Mmap<T> pooled;
            pooled.D0 = mm.D0;
            pooled.D1 = mm.D1;
            Matrix<T> A = mm.Dc[i];
            for (std::size_t r = 0; r < A.rows(); ++r)
                for (std::size_t s = 0; s < A.cols(); ++s) A(r, s) += mm.Dc[j](r, s);
            Matrix<T> B = mm.D1;
            for (std::size_t r = 0; r < B.rows(); ++r)
                for (std::size_t s = 0; s < B.cols(); ++s) B(r, s) -= A(r, s);
            pooled.Dc.push_back(A);
            pooled.Dc.push_back(B);
            const std::vector<T> pV = mmap_count_var(pooled, t);
            S(i, j) = half * (pV[0] - mV[i] - mV[j]);
        }
    return S;
}

/**
 * Cross moments of order k: MC(i,j) is E[T^k] of the interval that FOLLOWS a
 * class-i arrival, conditioned on that next arrival being of class j.
 */
template <class T>
Matrix<T> mmap_cross_moment(const Mmap<T>& mm, unsigned k) {
    const std::size_t C = mm.classes();
    const std::vector<T> pie = map_pie(mm.map());
    const std::vector<Matrix<T>> E = mmap_embedded(mm);
    Matrix<T> negD0 = mm.D0;
    for (std::size_t i = 0; i < negD0.rows(); ++i)
        for (std::size_t j = 0; j < negD0.cols(); ++j) negD0(i, j) = -negD0(i, j);
    const Matrix<T> M = inverse(negD0);
    const Matrix<T> Mk1 = matpow(M, k + 1);
    std::vector<T> TG(C, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < C; ++i) {
        const std::vector<T> v = vecmul(pie, E[i]);
        for (std::size_t r = 0; r < v.size(); ++r) TG[i] += v[r];
    }
    Matrix<T> MC(C, C, num_traits<T>::from_int(0));
    const T fk = num_factorial<T>(k);
    for (std::size_t i = 0; i < C; ++i) {
        std::vector<T> start = vecmul(pie, E[i]);
        for (std::size_t r = 0; r < start.size(); ++r) start[r] = start[r] / TG[i];
        for (std::size_t j = 0; j < C; ++j) {
            const std::vector<T> num = vecmul(vecmul(start, Mk1), mm.Dc[j]);
            const std::vector<T> den = vecmul(start, E[j]);
            T sn = num_traits<T>::from_int(0), sd = num_traits<T>::from_int(0);
            for (std::size_t r = 0; r < num.size(); ++r) sn += num[r];
            for (std::size_t r = 0; r < den.size(); ++r) sd += den[r];
            MC(i, j) = fk * sn / sd;
        }
    }
    return MC;
}

/**
 * Forward moments: MOMENTS(a,h) is the order-orders[h] moment of the interval
 * ENDING with a class-a arrival. With normalize false the per-class probability
 * is left in, which returns the unnormalized contribution instead.
 */
template <class T>
Matrix<T> mmap_forward_moment(const Mmap<T>& mm, const std::vector<unsigned>& orders,
                              bool normalize) {
    const std::size_t C = mm.classes();
    const std::vector<T> pie = map_pie(mm.map());
    const std::vector<Matrix<T>> E = mmap_embedded(mm);
    Matrix<T> negD0 = mm.D0;
    for (std::size_t i = 0; i < negD0.rows(); ++i)
        for (std::size_t j = 0; j < negD0.cols(); ++j) negD0(i, j) = -negD0(i, j);
    const Matrix<T> M = inverse(negD0);
    Matrix<T> out(C, orders.size(), num_traits<T>::from_int(0));
    for (std::size_t a = 0; a < C; ++a) {
        const std::vector<T> start = vecmul(pie, E[a]);
        T pa = num_traits<T>::from_int(1);
        if (normalize) {
            pa = num_traits<T>::from_int(0);
            for (std::size_t r = 0; r < start.size(); ++r) pa += start[r];
        }
        for (std::size_t h = 0; h < orders.size(); ++h) {
            const std::vector<T> v = vecmul(start, matpow(M, orders[h]));
            T acc = num_traits<T>::from_int(0);
            for (std::size_t r = 0; r < v.size(); ++r) acc += v[r];
            out(a, h) = num_factorial<T>(orders[h]) / pa * acc;
        }
    }
    return out;
}

/** Default normalization, matching the two-argument MATLAB call. */
template <class T>
Matrix<T> mmap_forward_moment(const Mmap<T>& mm, const std::vector<unsigned>& orders) {
    return mmap_forward_moment(mm, orders, true);
}

/** MMAP of the sum of n independent copies, a block bidiagonal concatenation. */
template <class T>
Mmap<T> mmap_sum(const Mmap<T>& mm, unsigned n) {
    if (n == 0) throw InputError("mmap_sum: n must be positive");
    const std::size_t ns = mm.order();
    const std::size_t C = mm.classes();
    const std::size_t N = ns * n;
    const T zero = num_traits<T>::from_int(0);
    Mmap<T> out;
    out.D0 = Matrix<T>(N, N, zero);
    out.D1 = Matrix<T>(N, N, zero);
    out.Dc.assign(C, Matrix<T>(N, N, zero));
    for (unsigned b = 0; b < n; ++b) {
        const std::size_t off = b * ns;
        for (std::size_t i = 0; i < ns; ++i)
            for (std::size_t j = 0; j < ns; ++j) out.D0(off + i, off + j) = mm.D0(i, j);
        if (b + 1 < n) {
            for (std::size_t i = 0; i < ns; ++i)
                for (std::size_t j = 0; j < ns; ++j) out.D0(off + i, off + ns + j) = mm.D1(i, j);
        } else {
            for (std::size_t i = 0; i < ns; ++i)
                for (std::size_t j = 0; j < ns; ++j) {
                    out.D1(off + i, j) = mm.D1(i, j);
                    for (std::size_t c = 0; c < C; ++c) out.Dc[c](off + i, j) = mm.Dc[c](i, j);
                }
        }
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif
