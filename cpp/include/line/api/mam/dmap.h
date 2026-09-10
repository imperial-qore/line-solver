/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_DMAP_H
#define LINE_API_MAM_DMAP_H

/**
 * Discrete-time Markovian arrival processes (D-MAPs).
 *
 * Templated port of matlab/lib/kpctoolbox/dmap: dmap_pie.m, dmap_moment.m,
 * dmap_isfeasible.m, dmap_exp_mul_int.m, dmap_dist.m, dmap_geo_mul_sum.m,
 * dmap_dist_acf.m, dmap_dist_lag1.m and dmap_sample.m.
 *
 * A D-MAP is the pair (D0, D1) of SUBSTOCHASTIC matrices with D0 + D1
 * stochastic: D0 carries a slot with no arrival and D1 a slot with one. The
 * continuous-time analogue has D0 with a negative diagonal and rows of D0 + D1
 * summing to ZERO, so a routine written for one representation silently
 * produces nonsense on the other. Interarrival times are counted in SLOTS and
 * are at least one, which is why the first moment is alpha (I - D0)^-1 1 and
 * not the continuous alpha (-D0)^-1 1.
 *
 * The distance functionals all reduce to Stein equations A X B - X + C = 0,
 * MATLAB's three-argument dlyap. They are solved here through the Kronecker
 * form, which is exact in the rational backend where a Schur-based solver could
 * not be.
 */

#include <cstddef>
#include <random>
#include <vector>

#include "line/api/mam/mmap_lambda.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** A discrete-time MAP: substochastic D0 (no arrival) and D1 (one arrival). */
template <class T>
struct Dmap {
    Matrix<T> D0;
    Matrix<T> D1;

    std::size_t order() const { return D0.rows(); }
};

namespace detail {

/** Column-major vectorization, matching the MATLAB reshape convention. */
template <class T>
std::vector<T> vec_colmajor(const Matrix<T>& A) {
    std::vector<T> v(A.rows() * A.cols());
    for (std::size_t j = 0; j < A.cols(); ++j)
        for (std::size_t i = 0; i < A.rows(); ++i) v[j * A.rows() + i] = A(i, j);
    return v;
}

/** MATLAB dlyap(A, B, C): solves A X B - X + C = 0 through the Kronecker form. */
template <class T>
Matrix<T> stein_solve(const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C) {
    const std::size_t m = A.rows(), n = B.cols();
    if (A.cols() != m || B.rows() != n)
        throw InputError("stein_solve: A and B must be square");
    if (C.rows() != m || C.cols() != n)
        throw InputError("stein_solve: the right-hand side is not conformable");
    const std::size_t N = m * n;
    Matrix<T> M(N, N, num_traits<T>::from_int(0));
    // vec(A X B) = kron(B', A) vec(X) in the column-major layout
    for (std::size_t jb = 0; jb < n; ++jb)
        for (std::size_t ib = 0; ib < n; ++ib)
            for (std::size_t ia = 0; ia < m; ++ia)
                for (std::size_t ja = 0; ja < m; ++ja)
                    M(jb * m + ia, ib * m + ja) = -B(ib, jb) * A(ia, ja);
    for (std::size_t k = 0; k < N; ++k) M(k, k) += num_traits<T>::from_int(1);
    const std::vector<T> x = solve(M, vec_colmajor(C));
    Matrix<T> X(m, n, num_traits<T>::from_int(0));
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < m; ++i) X(i, j) = x[j * m + i];
    return X;
}

/** (I - D0)^-1 D1, the embedded chain at arrival epochs. */
template <class T>
Matrix<T> dmap_embedded_chain(const Matrix<T>& D0, const Matrix<T>& D1) {
    const std::size_t n = D0.rows();
    Matrix<T> A = eye<T>(n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) A(i, j) -= D0(i, j);
    return matmul(inverse(A), D1);
}

}  // namespace detail

/** Stationary phase distribution at arrival epochs. */
template <class T>
std::vector<T> dmap_pie(const Dmap<T>& d) {
    return mc::dtmc_solve(detail::dmap_embedded_chain(d.D0, d.D1));
}

/** Raw moments of the interarrival time in slots, for orders 1, 2 and 3 only. */
template <class T>
std::vector<T> dmap_moment(const Dmap<T>& d, const std::vector<unsigned>& orders) {
    const std::size_t n = d.order();
    const std::vector<T> al = dmap_pie(d);
    Matrix<T> IA = eye<T>(n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) IA(i, j) -= d.D0(i, j);
    const Matrix<T> A = inverse(IA);
    const std::vector<T> e = ones<T>(n);
    const std::vector<T> Ae = mulvec(A, e);
    T m1 = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) m1 += al[i] * Ae[i];
    const std::vector<T> alA = vecmul(al, A);
    const std::vector<T> alAA = vecmul(alA, A);
    T alAAe = num_traits<T>::from_int(0), alAe = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        alAe += alA[i] * Ae[i];
        alAAe += alAA[i] * Ae[i];
    }
    std::vector<T> out(orders.size(), num_traits<T>::from_int(0));
    for (std::size_t t = 0; t < orders.size(); ++t) {
        switch (orders[t]) {
            case 1:
                out[t] = m1;
                break;
            case 2:
                out[t] = num_traits<T>::from_int(2) * alAe - m1;
                break;
            case 3:
                out[t] = num_traits<T>::from_int(6) * alAAe -
                         num_traits<T>::from_int(6) * alAe + m1;
                break;
            default:
                throw InputError("dmap_moment: raw moments of order > 3 are not implemented");
        }
    }
    return out;
}

/** True when D0 and D1 are nonnegative and D0 + D1 is stochastic. */
template <class T>
bool dmap_isfeasible(const Dmap<T>& d) {
    const std::size_t n = d.D0.rows();
    if (d.D0.cols() != n || d.D1.rows() != n || d.D1.cols() != n) return false;
    const T negtol = -num_traits<T>::from_double(1e-10);
    const T sumtol = num_traits<T>::from_double(1e-6);
    for (std::size_t i = 0; i < n; ++i) {
        T rs = num_traits<T>::from_int(0), d1rs = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < n; ++j) {
            if (d.D0(i, j) < negtol || d.D1(i, j) < negtol) return false;
            rs += d.D0(i, j) + d.D1(i, j);
            d1rs += d.D1(i, j);
        }
        if (num_abs(rs - num_traits<T>::from_int(1)) > sumtol) return false;
        if (d1rs < negtol) return false;
    }
    return true;
}

/** Inner product of the two interarrival densities truncated at lag L. */
template <class T>
T dmap_exp_mul_int(const Dmap<T>& a, const Dmap<T>& b, unsigned L, const std::vector<T>& alA,
                   const std::vector<T>& alB) {
    const std::size_t NA = a.order(), NB = b.order();
    Matrix<T> C(NB, NA, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < NB; ++i)
        for (std::size_t j = 0; j < NA; ++j) C(i, j) = alB[i] * alA[j];
    if (L == 0) throw InputError("dmap_exp_mul_int: the truncation lag L must be positive");
    Matrix<T> Z = detail::stein_solve(b.D0.transpose(), a.D0, C);
    for (unsigned i = 1; i <= L - 1; ++i) {
        Matrix<T> R = matmul(matmul(b.D1.transpose(), Z), a.D1);
        Z = detail::stein_solve(b.D0.transpose(), a.D0, R);
    }
    std::vector<T> dA(NA, num_traits<T>::from_int(0)), dB(NB, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < NA; ++i) {
        dA[i] = num_traits<T>::from_int(1);
        for (std::size_t j = 0; j < NA; ++j) dA[i] -= a.D0(i, j);
    }
    for (std::size_t i = 0; i < NB; ++i) {
        dB[i] = num_traits<T>::from_int(1);
        for (std::size_t j = 0; j < NB; ++j) dB[i] -= b.D0(i, j);
    }
    const std::vector<T> t = vecmul(dB, Z);
    T out = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < NA; ++i) out += t[i] * dA[i];
    return out;
}

/** Default stationary vectors, matching the three-argument MATLAB call. */
template <class T>
T dmap_exp_mul_int(const Dmap<T>& a, const Dmap<T>& b, unsigned L) {
    return dmap_exp_mul_int(a, b, L, dmap_pie(a), dmap_pie(b));
}

/** Squared L2 distance between the interarrival densities truncated at lag L. */
template <class T>
T dmap_dist(const Dmap<T>& a, const Dmap<T>& b, unsigned L, const std::vector<T>& alA,
            const std::vector<T>& alB) {
    return dmap_exp_mul_int(a, a, L + 1, alA, alA) -
           num_traits<T>::from_int(2) * dmap_exp_mul_int(a, b, L + 1, alA, alB) +
           dmap_exp_mul_int(b, b, L + 1, alB, alB);
}

/** Default stationary vectors, matching the three-argument MATLAB call. */
template <class T>
T dmap_dist(const Dmap<T>& a, const Dmap<T>& b, unsigned L) {
    return dmap_dist(a, b, L, dmap_pie(a), dmap_pie(b));
}

/**
 * Geometrically weighted sum of the lagged joint moments, the building block of
 * the autocorrelation distance. MATLAB returns 1/rcond as a large sentinel when
 * the Stein operator is near singular; here the linear solve is exact in the
 * rational backend and raises otherwise, so no sentinel is produced.
 */
template <class T>
T dmap_geo_mul_sum(const Dmap<T>& a, const Dmap<T>& b, const std::vector<T>& alA,
                   const std::vector<T>& alB) {
    const std::size_t NA = a.order(), NB = b.order();
    Matrix<T> IA = eye<T>(NA), IB = eye<T>(NB);
    for (std::size_t i = 0; i < NA; ++i)
        for (std::size_t j = 0; j < NA; ++j) IA(i, j) -= a.D0(i, j);
    for (std::size_t i = 0; i < NB; ++i)
        for (std::size_t j = 0; j < NB; ++j) IB(i, j) -= b.D0(i, j);
    const Matrix<T> D0Ai = inverse(IA), D0Bi = inverse(IB);
    Matrix<T> PAh = matmul(D0Ai, a.D1), PBh = matmul(D0Bi, b.D1);
    for (std::size_t i = 0; i < NA; ++i)
        for (std::size_t j = 0; j < NA; ++j) PAh(i, j) -= alA[j];
    for (std::size_t i = 0; i < NB; ++i)
        for (std::size_t j = 0; j < NB; ++j) PBh(i, j) -= alB[j];
    std::vector<T> rowA(NA, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < NA; ++i)
        for (std::size_t j = 0; j < NA; ++j) rowA[i] += D0Ai(i, j);
    const std::vector<T> alBD0Bi = vecmul(alB, D0Bi);
    Matrix<T> C(NA, NB, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < NA; ++i)
        for (std::size_t j = 0; j < NB; ++j) C(i, j) = rowA[i] * alBD0Bi[j];
    const Matrix<T> X = detail::stein_solve(PAh, PBh, C);
    const std::vector<T> v = vecmul(vecmul(vecmul(alA, D0Ai), X), D0Bi);
    T out = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < NB; ++i) out += v[i];
    return out;
}

/** Default stationary vectors, matching the two-argument MATLAB call. */
template <class T>
T dmap_geo_mul_sum(const Dmap<T>& a, const Dmap<T>& b) {
    return dmap_geo_mul_sum(a, b, dmap_pie(a), dmap_pie(b));
}

/** Squared distance between the autocorrelation structures of two D-MAPs. */
template <class T>
T dmap_dist_acf(const Dmap<T>& a, const Dmap<T>& b, const std::vector<T>& alA,
                const std::vector<T>& alB) {
    std::vector<unsigned> ord;
    ord.push_back(1);
    ord.push_back(2);
    const std::vector<T> momA = dmap_moment(a, ord);
    const std::vector<T> momB = dmap_moment(b, ord);
    const T varA = momA[1] - momA[0] * momA[0];
    const T varB = momB[1] - momB[0] * momB[0];
    const T four = num_traits<T>::from_int(4);
    return (dmap_geo_mul_sum(a, a, alA, alA) - momA[1] * momA[1] / four) / (varA * varA) -
           num_traits<T>::from_int(2) *
               (dmap_geo_mul_sum(a, b, alA, alB) - momA[1] * momB[1] / four) / (varA * varB) +
           (dmap_geo_mul_sum(b, b, alB, alB) - momB[1] * momB[1] / four) / (varB * varB);
}

/** Default stationary vectors, matching the two-argument MATLAB call. */
template <class T>
T dmap_dist_acf(const Dmap<T>& a, const Dmap<T>& b) {
    return dmap_dist_acf(a, b, dmap_pie(a), dmap_pie(b));
}

/** Squared distance between the lag-1 joint densities of two D-MAPs. */
template <class T>
T dmap_dist_lag1(const Dmap<T>& a, const Dmap<T>& b, const std::vector<T>& alA,
                 const std::vector<T>& alB) {
    const std::size_t NA = a.order(), NB = b.order();
    std::vector<T> ea(NA, num_traits<T>::from_int(0)), eb(NB, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < NA; ++i) {
        ea[i] = num_traits<T>::from_int(1);
        for (std::size_t j = 0; j < NA; ++j) ea[i] -= a.D0(i, j);
    }
    for (std::size_t i = 0; i < NB; ++i) {
        eb[i] = num_traits<T>::from_int(1);
        for (std::size_t j = 0; j < NB; ++j) eb[i] -= b.D0(i, j);
    }
    Matrix<T> Cab(NA, NB, num_traits<T>::from_int(0)), Caa(NA, NA, num_traits<T>::from_int(0)),
        Cbb(NB, NB, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < NA; ++i)
        for (std::size_t j = 0; j < NB; ++j) Cab(i, j) = alA[i] * alB[j];
    for (std::size_t i = 0; i < NA; ++i)
        for (std::size_t j = 0; j < NA; ++j) Caa(i, j) = alA[i] * alA[j];
    for (std::size_t i = 0; i < NB; ++i)
        for (std::size_t j = 0; j < NB; ++j) Cbb(i, j) = alB[i] * alB[j];
    const Matrix<T> Z_AB = detail::stein_solve(a.D0.transpose(), b.D0, Cab);
    const Matrix<T> Z_AA = detail::stein_solve(a.D0.transpose(), a.D0, Caa);
    const Matrix<T> Z_BB = detail::stein_solve(b.D0.transpose(), b.D0, Cbb);
    Matrix<T> Rab(NA, NB, num_traits<T>::from_int(0)), Raa(NA, NA, num_traits<T>::from_int(0)),
        Rbb(NB, NB, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < NA; ++i)
        for (std::size_t j = 0; j < NB; ++j) Rab(i, j) = ea[i] * eb[j];
    for (std::size_t i = 0; i < NA; ++i)
        for (std::size_t j = 0; j < NA; ++j) Raa(i, j) = ea[i] * ea[j];
    for (std::size_t i = 0; i < NB; ++i)
        for (std::size_t j = 0; j < NB; ++j) Rbb(i, j) = eb[i] * eb[j];
    const Matrix<T> X_AB = detail::stein_solve(a.D0, b.D0.transpose(), Rab);
    const Matrix<T> X_AA = detail::stein_solve(a.D0, a.D0.transpose(), Raa);
    const Matrix<T> X_BB = detail::stein_solve(b.D0, b.D0.transpose(), Rbb);
    const std::vector<T> vA = detail::vec_colmajor(a.D1);
    const std::vector<T> vB = detail::vec_colmajor(b.D1);
    // v' kron(X, Z) w with the column-major vec, i.e. sum over the four indices
    T sBB = num_traits<T>::from_int(0), sAA = num_traits<T>::from_int(0),
      sAB = num_traits<T>::from_int(0);
    for (std::size_t p = 0; p < NB; ++p)
        for (std::size_t q = 0; q < NB; ++q)
            for (std::size_t r = 0; r < NB; ++r)
                for (std::size_t s = 0; s < NB; ++s)
                    sBB += vB[p * NB + q] * X_BB(p, r) * Z_BB(q, s) * vB[r * NB + s];
    for (std::size_t p = 0; p < NA; ++p)
        for (std::size_t q = 0; q < NA; ++q)
            for (std::size_t r = 0; r < NA; ++r)
                for (std::size_t s = 0; s < NA; ++s)
                    sAA += vA[p * NA + q] * X_AA(p, r) * Z_AA(q, s) * vA[r * NA + s];
    for (std::size_t p = 0; p < NA; ++p)
        for (std::size_t q = 0; q < NA; ++q)
            for (std::size_t r = 0; r < NB; ++r)
                for (std::size_t s = 0; s < NB; ++s)
                    sAB += vA[p * NA + q] * X_AB(p, r) * Z_AB(q, s) * vB[r * NB + s];
    return sBB + sAA - num_traits<T>::from_int(2) * sAB;
}

/** Default stationary vectors, matching the two-argument MATLAB call. */
template <class T>
T dmap_dist_lag1(const Dmap<T>& a, const Dmap<T>& b) {
    return dmap_dist_lag1(a, b, dmap_pie(a), dmap_pie(b));
}

/** n interarrival times in slots, drawn by walking the phase process. */
template <class T, class Gen>
std::vector<unsigned> dmap_sample(const Dmap<T>& d, std::size_t n, Gen& gen) {
    const std::size_t N = d.order();
    const std::vector<T> al = dmap_pie(d);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::size_t phase = 0;
    double u = unif(gen), acc = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        acc += num_traits<T>::to_double(al[i]);
        if (u <= acc) {
            phase = i;
            break;
        }
    }
    std::vector<unsigned> X(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        unsigned t = 0;
        while (true) {
            ++t;
            const double v = unif(gen);
            double c = 0.0;
            std::size_t next = 2 * N - 1;
            for (std::size_t j = 0; j < 2 * N; ++j) {
                c += num_traits<T>::to_double(j < N ? d.D0(phase, j) : d.D1(phase, j - N));
                if (v <= c) {
                    next = j;
                    break;
                }
            }
            if (next >= N) {
                phase = next - N;
                X[i] = t;
                break;
            }
            phase = next;
        }
    }
    return X;
}

}  // namespace mam
}  // namespace line

#endif
