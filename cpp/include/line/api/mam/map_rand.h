/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_RAND_H
#define LINE_API_MAM_MAP_RAND_H

/**
 * Random MAP, MMPP, MMAP, acyclic-PH and hyperexponential generators, plus the
 * hyperexponential reader.
 *
 * Templated port of matlab/lib/kpctoolbox/map/map_rand.m, map_randn.m,
 * matlab/lib/kpctoolbox/mmpp/mmpp_rand.m, matlab/lib/kpctoolbox/aph/aph_rand.m,
 * hyper_rand.m, ph2hyper.m and matlab/lib/m3a/m3a/mmap/mmap_rand.m.
 *
 * The generators take an explicit engine, following ctmc_rand, so a caller can
 * reproduce a draw; they therefore do NOT reproduce the MATLAB stream and are
 * only distributionally equivalent to it. Everything is drawn on the raw
 * matrices and then passed through map_normalize, which is what makes the
 * result a generator rather than a nonnegative matrix pair.
 *
 * mmap_rand draws its class split once. MATLAB redraws it inside a loop over
 * the phases and keeps only the last draw, so the two agree in law and differ
 * only in how far the stream is advanced.
 */

#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** Random MAP of order K with uniform [0,1) entries, normalized. */
template <class T, class Gen>
Map<T> map_rand(std::size_t K, Gen& gen) {
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    Matrix<T> D0(K, K, num_traits<T>::from_int(0)), D1(K, K, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t j = 0; j < K; ++j) {
            D0(i, j) = num_traits<T>::from_double(unif(gen));
            D1(i, j) = num_traits<T>::from_double(unif(gen));
        }
    return map_normalize(Map<T>{D0, D1});
}

/** Random MAP of order K with folded normal entries, normalized. */
template <class T, class Gen>
Map<T> map_randn(std::size_t K, double mu, double sigma, Gen& gen) {
    std::normal_distribution<double> nrm(mu, sigma);
    Matrix<T> D0(K, K, num_traits<T>::from_int(0)), D1(K, K, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t j = 0; j < K; ++j) {
            D0(i, j) = num_traits<T>::from_double(std::fabs(nrm(gen)));
            D1(i, j) = num_traits<T>::from_double(std::fabs(nrm(gen)));
        }
    return map_normalize(Map<T>{D0, D1});
}

/** Random MMPP of order K: the arrival matrix is diagonal, so arrivals do not switch phase. */
template <class T, class Gen>
Map<T> mmpp_rand(std::size_t K, Gen& gen) {
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    Matrix<T> D0(K, K, num_traits<T>::from_int(0)), D1(K, K, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t j = 0; j < K; ++j) D0(i, j) = num_traits<T>::from_double(unif(gen));
    for (std::size_t i = 0; i < K; ++i) {
        for (std::size_t j = 0; j < K; ++j) {
            const double v = unif(gen);
            if (i == j) D1(i, j) = num_traits<T>::from_double(v);
        }
    }
    return map_normalize(Map<T>{D0, D1});
}

/**
 * Random M3PP of order K with `classes` marks (m3a/m3pp/m3pp_rand.m).
 *
 * The underlying process is `mmpp_rand`; the marks are a class-independent
 * Bernoulli split of D1 with probabilities drawn from a uniform simplex, so
 * sum_c Dc = D1 by construction. The reference wraps the draw in two further
 * loops over the phases whose bodies overwrite the same matrices, so only the
 * last draw survives; that is what a single draw here reproduces.
 */
template <class T, class Gen>
Mmap<T> m3pp_rand(std::size_t K, std::size_t classes, Gen& gen) {
    if (classes == 0) throw InputError("m3pp_rand: at least one class is required");
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    const Map<T> base = mmpp_rand<T>(K, gen);
    std::vector<T> p(classes);
    T tot = num_traits<T>::from_int(0);
    for (std::size_t c = 0; c < classes; ++c) {
        p[c] = num_traits<T>::from_double(unif(gen));
        tot += p[c];
    }
    Mmap<T> out;
    out.D0 = base.D0;
    out.D1 = base.D1;
    for (std::size_t c = 0; c < classes; ++c) {
        Matrix<T> Dc(K, K, num_traits<T>::from_int(0));
        for (std::size_t i = 0; i < K; ++i)
            for (std::size_t j = 0; j < K; ++j) Dc(i, j) = base.D1(i, j) * p[c] / tot;
        out.Dc.push_back(Dc);
    }
    return out;
}

/** Random acyclic PH renewal process of order K, upper triangular in D0. */
template <class T, class Gen>
Map<T> aph_rand(std::size_t K, Gen& gen) {
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    Matrix<T> D0(K, K, num_traits<T>::from_int(0)), D1(K, K, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t j = 0; j < K; ++j) D1(i, j) = num_traits<T>::from_double(unif(gen));
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t j = 0; j < K; ++j) {
            const double v = unif(gen);
            if (j >= i) D0(i, j) = num_traits<T>::from_double(v);
        }
    return map_normalize(map_renewal(Map<T>{D0, D1}));
}

/** Random hyperexponential of order k, given as a MAP. */
template <class T, class Gen>
Map<T> hyper_rand(std::size_t k, Gen& gen) {
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::vector<T> v(k), alpha(k);
    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < k; ++i) {
        v[i] = num_traits<T>::from_double(unif(gen));
        alpha[i] = num_traits<T>::from_double(unif(gen));
        s += alpha[i];
    }
    for (std::size_t i = 0; i < k; ++i) alpha[i] = alpha[i] / s;
    Matrix<T> H0(k, k, num_traits<T>::from_int(0)), H1(k, k, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < k; ++i) {
        H0(i, i) = -v[i];
        for (std::size_t j = 0; j < k; ++j) H1(i, j) = v[i] * alpha[j];
    }
    return Map<T>{H0, H1};
}

/** Random MMAP of the given order with a random split of D1 across the classes. */
template <class T, class Gen>
Mmap<T> mmap_rand(std::size_t order, std::size_t classes, Gen& gen) {
    if (classes == 0) throw InputError("mmap_rand: at least one class is required");
    const Map<T> base = map_rand<T>(order, gen);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::vector<T> p(classes);
    T s = num_traits<T>::from_int(0);
    for (std::size_t c = 0; c < classes; ++c) {
        p[c] = num_traits<T>::from_double(unif(gen));
        s += p[c];
    }
    Mmap<T> out;
    out.D0 = base.D0;
    out.D1 = base.D1;
    for (std::size_t c = 0; c < classes; ++c) {
        Matrix<T> Dc = base.D1;
        for (std::size_t i = 0; i < Dc.rows(); ++i)
            for (std::size_t j = 0; j < Dc.cols(); ++j) Dc(i, j) = Dc(i, j) * p[c] / s;
        out.Dc.push_back(Dc);
    }
    return out;
}

/** Rates and branch probabilities of a hyperexponential given as a MAP. */
template <class T>
struct HyperParams {
    std::vector<T> lambda;
    std::vector<T> prob;
};

/** Reads a hyperexponential MAP back into rates and branch probabilities. */
template <class T>
HyperParams<T> ph2hyper(const Map<T>& ph) {
    const std::size_t n = ph.D0.rows();
    const T tol = num_traits<T>::from_double(1e-10);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (i != j && num_abs(ph.D0(i, j)) > tol)
                throw InputError("ph2hyper: the PH distribution is not hyperexponential");
    Matrix<T> negD0(n, n, num_traits<T>::from_int(0));
    HyperParams<T> out;
    out.lambda.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        out.lambda[i] = -ph.D0(i, i);
        negD0(i, i) = out.lambda[i];
    }
    out.prob = mc::dtmc_solve(matmul(inverse(negD0), ph.D1));
    return out;
}

}  // namespace mam
}  // namespace line

#endif
