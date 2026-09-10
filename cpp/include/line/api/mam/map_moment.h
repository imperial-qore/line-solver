/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_MOMENT_H
#define LINE_API_MAM_MAP_MOMENT_H

/**
 * Markovian arrival process descriptors: stationary vectors, rate, moments,
 * autocorrelation and the index of dispersion.
 *
 * Templated port of the kpctoolbox MAP primitives used across LINE
 * (matlab/lib/kpctoolbox/map/map_prob.m, map_pie.m, map_lambda.m, map_mean.m,
 * map_moment.m, map_var.m, map_scv.m, map_embedded.m, map_acf.m, map_idc.m,
 * map_infgen.m).
 *
 * A MAP is the pair (D0, D1): D0 carries the hidden transitions and the
 * negative diagonal, D1 the transitions that emit an arrival. The whole family
 * is rational in the entries of D0 and D1 -- stationary vectors are linear
 * solves, moments are i! pie (-D0)^-i e -- so all of it is exact in rational
 * arithmetic. That matters for fitting work, where the moments of a candidate
 * MAP are compared against targets and a rounding artifact is easy to mistake
 * for a fitting error.
 */

#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * Tolerance exponent shared by the KPC feasibility checks (map_feastol.m).
 *
 * It lives HERE, in the header every MAP user already includes, because it is
 * odr-used from `map_transform.h` (map_isfeasible) and from `map_dist.h`. A
 * declaration in one header and a definition in another that does not include
 * it links only while some third translation unit happens to emit the weak
 * symbol; a program built from the declaring header alone would not link.
 */
inline int map_feastol() { return 8; }

/** A MAP as the pair of matrices (D0, D1). */
template <class T>
struct Map {
    Matrix<T> D0;
    Matrix<T> D1;

    std::size_t order() const { return D0.rows(); }
};

/** Generator of the underlying phase process, D0 + D1. */
template <class T>
Matrix<T> map_infgen(const Map<T>& m) {
    if (m.D0.rows() != m.D1.rows() || m.D0.cols() != m.D1.cols())
        throw InputError("map_infgen: D0 and D1 have different shapes");
    Matrix<T> Q(m.D0.rows(), m.D0.cols());
    for (std::size_t i = 0; i < Q.rows(); ++i)
        for (std::size_t j = 0; j < Q.cols(); ++j) Q(i, j) = m.D0(i, j) + m.D1(i, j);
    return Q;
}

/** Stationary distribution of the phase process, pi (D0 + D1) = 0. */
template <class T>
std::vector<T> map_prob(const Map<T>& m) {
    return mc::ctmc_solve(map_infgen(m));
}

/** Stationary arrival rate, lambda = pi D1 e. */
template <class T>
T map_lambda(const Map<T>& m) {
    const std::vector<T> p = map_prob(m);
    const std::vector<T> pD1 = vecmul(p, m.D1);
    T s = num_traits<T>::from_int(0);
    for (const T& v : pD1) s += v;
    return s;
}

/** Phase distribution seen by an arriving job, pie = pi D1 / (pi D1 e). */
template <class T>
std::vector<T> map_pie(const Map<T>& m) {
    const std::vector<T> p = map_prob(m);
    std::vector<T> a = vecmul(p, m.D1);
    T s = num_traits<T>::from_int(0);
    for (const T& v : a) s += v;
    if (s == num_traits<T>::from_int(0)) throw NumericError("map_pie: the MAP has zero arrival rate");
    for (T& v : a) v /= s;
    return a;
}

/** Mean inter-arrival time, 1/lambda. */
template <class T>
T map_mean(const Map<T>& m) {
    const T lam = map_lambda(m);
    if (lam == num_traits<T>::from_int(0)) throw NumericError("map_mean: zero arrival rate");
    return num_traits<T>::from_int(1) / lam;
}

/** Embedded DTMC at arrival epochs, P = (-D0)^-1 D1. */
template <class T>
Matrix<T> map_embedded(const Map<T>& m) {
    Matrix<T> negD0 = m.D0;
    for (std::size_t i = 0; i < negD0.rows(); ++i)
        for (std::size_t j = 0; j < negD0.cols(); ++j) negD0(i, j) = -negD0(i, j);
    return matmul(inverse(negD0), m.D1);
}

/** Raw moment of order k of the inter-arrival time: k! pie (-D0)^-k e. */
template <class T>
T map_moment(const Map<T>& m, unsigned k) {
    if (k == 0) return num_traits<T>::from_int(1);
    Matrix<T> negD0 = m.D0;
    for (std::size_t i = 0; i < negD0.rows(); ++i)
        for (std::size_t j = 0; j < negD0.cols(); ++j) negD0(i, j) = -negD0(i, j);
    const Matrix<T> A = matpow(inverse(negD0), k);
    const std::vector<T> x = map_pie(m);
    const std::vector<T> xA = vecmul(x, A);
    T s = num_traits<T>::from_int(0);
    for (const T& v : xA) s += v;
    return num_factorial<T>(k) * s;
}

/** Variance of the inter-arrival time. */
template <class T>
T map_var(const Map<T>& m) {
    const T mu = map_mean(m);
    return map_moment(m, 2) - mu * mu;
}

/** Squared coefficient of variation. */
template <class T>
T map_scv(const Map<T>& m) {
    const T mu = map_mean(m);
    return map_var(m) / (mu * mu);
}

/**
 * Autocorrelation coefficients of the inter-arrival times at the given lags,
 *
 *     rho_k = (x P^k y - 1) / scv,   x = lambda pi,  y = (-D0)^-1 e,
 *
 * which is matlab/lib/kpctoolbox/map/map_acf.m in full, closing line included.
 *
 * The closing normalization USED TO BE MISSING here, and the raw kpctoolbox
 * quantity x P^k y was returned instead. That is not an autocorrelation
 * coefficient: it is 1, not 0, for a renewal process, and it is unbounded
 * above rather than confined to [-1, 1]. On the two-class MMAP of
 * tests/test_qbd_family.cpp the unnormalized form gave 1.13166 where MATLAB
 * gives 0.0739911109309197, and every renewal MAP read 1 where MATLAB reads 0
 * (checked directly: map_acf({[-1]},{[1]}, [1 2]) is [0 0] in MATLAB, and the
 * raw quantity is 1). The ratio rho_{k+1}/rho_k that amap2_adjust_gamma uses
 * as its decay characteristic was correspondingly wrong, since scv cancels in
 * the ratio but the -1 does not: MATLAB's FGAMMA is (raw_4 - 1)/(raw_3 - 1)
 * and the port was computing raw_4/raw_3.
 *
 * Exact in rational arithmetic: the normalization is one subtraction and one
 * division, so a renewal MAP returns identically zero rather than 1e-16.
 */
template <class T>
std::vector<T> map_acf(const Map<T>& m, const std::vector<unsigned>& lags) {
    const Matrix<T> P = map_embedded(m);
    const T lam = map_lambda(m);
    std::vector<T> x = map_prob(m);
    for (T& v : x) v *= lam;
    Matrix<T> negD0 = m.D0;
    for (std::size_t i = 0; i < negD0.rows(); ++i)
        for (std::size_t j = 0; j < negD0.cols(); ++j) negD0(i, j) = -negD0(i, j);
    const std::vector<T> y = mulvec(inverse(negD0), ones<T>(m.order()));

    const T scv = map_scv(m);
    if (scv == num_traits<T>::from_int(0))
        throw NumericError("map_acf: the inter-arrival time is deterministic, the autocorrelation "
                           "coefficient is undefined (zero variance)");

    std::vector<T> out;
    out.reserve(lags.size());
    for (unsigned lag : lags) {
        const std::vector<T> xP = vecmul(x, matpow(P, lag));
        T s = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < xP.size(); ++i) s += xP[i] * y[i];
        out.push_back(T((s - num_traits<T>::from_int(1)) / scv));
    }
    return out;
}

/** Index of dispersion for counts, I = 1 + 2(lambda - pie (Q + e pi)^-1 D1 e). */
template <class T>
T map_idc(const Map<T>& m) {
    const std::size_t n = m.order();
    const Matrix<T> Q = map_infgen(m);
    const std::vector<T> p = map_prob(m);
    Matrix<T> A(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) A(i, j) = Q(i, j) + p[j];  // Q + e pi
    const std::vector<T> pie = map_pie(m);
    const std::vector<T> t1 = vecmul(pie, inverse(A));
    const std::vector<T> t2 = vecmul(t1, m.D1);
    T s = num_traits<T>::from_int(0);
    for (const T& v : t2) s += v;
    return num_traits<T>::from_int(1) + num_traits<T>::from_int(2) * (map_lambda(m) - s);
}

/** Two-phase MAP constructor for a Poisson process of rate lambda. */
template <class T>
Map<T> map_exponential(const T& lambda) {
    Map<T> m;
    m.D0 = Matrix<T>(1, 1);
    m.D1 = Matrix<T>(1, 1);
    m.D0(0, 0) = -lambda;
    m.D1(0, 0) = lambda;
    return m;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_MOMENT_H
