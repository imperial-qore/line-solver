/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SUM_SUM_CLOSED_H
#define LINE_API_SUM_SUM_CLOSED_H

/**
 * Summation method (SUM) and its extension (ESUM) for closed queueing
 * networks, including non-product-form stations with generally distributed
 * service times.
 *
 * Templated port of matlab/src/api/sum/sum_closed.m, cross-checked against
 * jar/src/main/java/jline/api/sum/Sum_closed.java (the two agree step for
 * step, including the Gauss-Seidel multiclass variant and the incremental
 * Erlang-C evaluation).
 *
 * The method writes the mean queue length of a station as a function of its
 * throughput, K_i = f_i(lambda_i), and closes the model with the population
 * constraint sum_i K_i + lambda Z = K. A single class is solved by bisection
 * on the system throughput (Bolch et al., Sec. 9.2.1); several classes by
 * Gauss-Seidel sweeps of per-class bisections on the per-class constraints,
 * which is more robust than the successive substitution of Sec. 9.2.2 because
 * it cannot overshoot the saturation polytope.
 *
 * Node functions:
 *   - product-form stations (scv = 1, or an insensitive discipline for which
 *     the caller passes scv = 1): Eq. (9.15)/(9.19)
 *   - FCFS with general service (scv != 1): the ESUM corrections, Eq. (10.88)
 *     for -/G/1 and Eq. (10.89) for -/G/m, with a_i = (1+scv_i)/2 and the
 *     Erlang-C waiting probability
 *   - infinite-server stations and think time: K_i = lambda_i L_i
 *
 * Reference: G. Bolch, S. Greiner, H. de Meer, K.S. Trivedi, Queueing Networks
 * and Markov Chains, 2nd ed., Wiley, 2006, Secs. 9.2 and 10.1.4.4.
 *
 * ARITHMETIC: a bisection stopped on a tolerance, so the answer is the root of
 * the population constraint only to within tol whatever the arithmetic. Gated
 * on has_transcendental for that reason; the Erlang-C evaluation is itself a
 * finite field computation (the a^k/k! terms are built incrementally) and
 * needs no transcendental function.
 *
 * Infinity is carried by the explicit Servers::infinite flag rather than by a
 * floating infinity, so the same code compiles for a number type that has no
 * infinity at all. An infinite population is not accepted here: MATLAB's
 * N(r) = Inf is only ever produced by sum_closing, which substitutes the
 * closing population before calling in.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace sum {

/** Number of servers of a station; MATLAB's mi(i) = Inf becomes infinite. */
struct Servers {
    long m = 1;
    bool infinite = false;

    static Servers of(long m) {
        Servers s;
        s.m = m;
        s.infinite = false;
        return s;
    }
    static Servers inf() {
        Servers s;
        s.m = 0;
        s.infinite = true;
        return s;
    }
};

/** Mirrors the [XN, QN, UN, RN, it] return list of the MATLAB function. */
template <class T>
struct SumClosedResult {
    std::vector<T> XN;  ///< (R) class throughputs
    Matrix<T> QN;       ///< (M x R) mean queue lengths
    Matrix<T> UN;       ///< (M x R) utilizations, per server at queueing stations
    Matrix<T> RN;       ///< (M x R) residence times, QN/XN
    std::size_t it = 0; ///< iterations of the outer loop
};

/** Convergence controls, mirroring the trailing (tol, maxiter) arguments. */
struct SumOptions {
    double tol = 1e-6;
    std::size_t maxiter = 10000;
};

namespace detail {

/**
 * Erlang-C waiting probability of an M/M/m queue, Eq. (6.28), built from the
 * incremental terms a^k/k! so that no factorial and no real power is formed.
 */
template <class T>
T sum_erlangc(long m, const T& rho) {
    const T one = num_traits<T>::from_int(1);
    if (!(rho < one)) return one;
    const T a = num_traits<T>::from_int(m) * rho;
    T s = num_traits<T>::from_int(0);
    T term = one;  // a^k / k!
    for (long k = 0; k < m; ++k) {
        if (k > 0) term *= a / num_traits<T>::from_int(k);
        s += term;
    }
    const T last = term * a / num_traits<T>::from_int(m) / (one - rho);
    return last / (s + last);
}

/** Per-station per-class mean queue lengths K_ir = f_ir(lambda_r). */
template <class T>
Matrix<T> sum_node_qlen(const Matrix<T>& L, const std::vector<T>& XN,
                        const std::vector<Servers>& mi, const Matrix<T>& scv, long K) {
    const std::size_t M = L.rows(), R = L.cols();
    const T one = num_traits<T>::from_int(1);
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> Qir(M, R, zero);
    for (std::size_t i = 0; i < M; ++i) {
        if (mi[i].infinite) {
            for (std::size_t r = 0; r < R; ++r) Qir(i, r) = XN[r] * L(i, r);  // Type 3, Eq. (9.15)
            continue;
        }
        const long m = mi[i].m;
        std::vector<T> Uir(R);
        T Ui = zero, ci2num = zero;
        for (std::size_t r = 0; r < R; ++r) {
            Uir[r] = XN[r] * L(i, r);
            Ui += Uir[r];
            ci2num += Uir[r] * scv(i, r);
        }
        if (Ui == zero) continue;
        // per-server utilization; the correction factors below stay finite at
        // rho = 1 because the node function is capped by the population
        T rho = Ui / num_traits<T>::from_int(m);
        if (rho > one) rho = one;
        const T ci2 = ci2num / Ui;  // demand-weighted node service SCV
        const T ai = (one + ci2) / num_traits<T>::from_int(2);
        if (K <= m) {
            // never more than m jobs at an m-server node: no queueing at all
            for (std::size_t r = 0; r < R; ++r) Qir(i, r) = Uir[r];
            continue;
        }
        const T Kt = num_traits<T>::from_int(K);
        const T mt = num_traits<T>::from_int(m);
        if (m == 1) {
            if (ci2 == one || K <= 1) {
                // Type 1, 2, 4 with m = 1, Eq. (9.15)/(9.19)
                const T den = one - (Kt - one) / Kt * rho;
                for (std::size_t r = 0; r < R; ++r) Qir(i, r) = Uir[r] / den;
            } else {
                // -/G/1 FCFS, Eq. (10.88)
                const T den = one - (Kt - one - ai) / (Kt - one) * rho;
                for (std::size_t r = 0; r < R; ++r) Qir(i, r) = Uir[r] * (one + rho * ai / den);
            }
        } else {
            const T Pm = sum_erlangc(m, rho);
            if (ci2 == one) {
                // Type 1 with m > 1, Eq. (9.15)/(9.19)
                const T den = one - (Kt - mt - one) / (Kt - mt) * rho;
                for (std::size_t r = 0; r < R; ++r) Qir(i, r) = Uir[r] + (Uir[r] / mt) * Pm / den;
            } else {
                // -/G/m FCFS, Eq. (10.89)
                const T den = one - (Kt - mt - ai) / (Kt - mt) * rho;
                for (std::size_t r = 0; r < R; ++r)
                    Qir(i, r) = Uir[r] + (Uir[r] / mt) * ai * Pm / den;
            }
        }
    }
    return Qir;
}

}  // namespace detail

/**
 * @param L   (M x R) service demands, L(i,r) = e(i,r)/mu(i,r)
 * @param N   (R) population per class, all finite
 * @param Z   (R) think times
 * @param mi  (M) servers per station
 * @param scv (M x R) squared coefficients of variation of the service times
 * @param options tolerances, iteration caps and the closing method
 */
template <class T>
SumClosedResult<T> sum_closed(const Matrix<T>& L, const std::vector<long>& N,
                              const std::vector<T>& Z, const std::vector<Servers>& mi,
                              const Matrix<T>& scv, const SumOptions& options = SumOptions()) {
    static_assert(num_traits<T>::has_transcendental,
                  "sum_closed requires transcendental arithmetic: it locates the root of the "
                  "population constraint by bisection to within tol, so its answer is inexact "
                  "whatever the arithmetic");
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("sum_closed: L and N disagree on the class count");
    if (Z.size() != R) throw InputError("sum_closed: L and Z disagree on the class count");
    if (mi.size() != M) throw InputError("sum_closed: L and mi disagree on the station count");
    if (scv.rows() != M || scv.cols() != R) throw InputError("sum_closed: scv has the wrong shape");
    for (long v : N)
        if (v < 0) throw InputError("sum_closed: negative population");

    const T zero = num_traits<T>::from_int(0);
    const T two = num_traits<T>::from_int(2);
    SumClosedResult<T> out;
    out.XN.assign(R, zero);
    out.QN = Matrix<T>(M, R, zero);
    out.UN = Matrix<T>(M, R, zero);
    out.RN = Matrix<T>(M, R, zero);

    long K = 0;
    for (long v : N) K += v;
    if (K == 0) return out;

    if (R == 1) {
        // single class: bisection on the system throughput (Sec. 9.2.1)
        T lambda_l = zero, lambda_u = zero;
        bool ub_set = false;
        for (std::size_t i = 0; i < M; ++i) {
            if (L(i, 0) > zero) {
                const T cand = mi[i].infinite ? num_traits<T>::from_int(K) / L(i, 0)
                                              : num_traits<T>::from_int(mi[i].m) / L(i, 0);
                if (!ub_set || cand < lambda_u) {
                    lambda_u = cand;
                    ub_set = true;
                }
            }
        }
        if (Z[0] > zero) {
            const T cand = num_traits<T>::from_int(K) / Z[0];
            if (!ub_set || cand < lambda_u) {
                lambda_u = cand;
                ub_set = true;
            }
        }
        if (!ub_set) throw InputError("sum_closed: all service demands are zero");

        T lambda = lambda_u;
        for (out.it = 1; out.it <= options.maxiter; ++out.it) {
            lambda = (lambda_l + lambda_u) / two;
            std::vector<T> X(1, lambda);
            const Matrix<T> Qir = detail::sum_node_qlen(L, X, mi, scv, K);
            T g = lambda * Z[0];
            for (std::size_t i = 0; i < M; ++i) g += Qir(i, 0);
            const double gap = num_traits<T>::to_double(num_abs(T(g - num_traits<T>::from_int(K))));
            const double width = num_traits<T>::to_double(T(lambda_u - lambda_l));
            if (gap <= options.tol ||
                width <= options.tol * num_traits<T>::to_double(lambda_u))
                break;
            if (g > num_traits<T>::from_int(K))
                lambda_u = lambda;
            else
                lambda_l = lambda;
        }
        if (out.it > options.maxiter) out.it = options.maxiter;
        out.XN[0] = lambda;
    } else {
        // multiclass: Gauss-Seidel sweeps of per-class bisections
        for (out.it = 1; out.it <= options.maxiter; ++out.it) {
            double delta = 0.0;
            for (std::size_t r = 0; r < R; ++r) {
                if (N[r] == 0) continue;
                T ub = zero;
                bool ub_set = false;
                for (std::size_t i = 0; i < M; ++i) {
                    if (!(L(i, r) > zero)) continue;
                    T cand;
                    if (mi[i].infinite) {
                        cand = num_traits<T>::from_int(K) / L(i, r);
                    } else {
                        T rowload = zero;
                        for (std::size_t q = 0; q < R; ++q) rowload += out.XN[q] * L(i, q);
                        T rem = num_traits<T>::from_int(mi[i].m) - (rowload - out.XN[r] * L(i, r));
                        if (rem < zero) rem = zero;
                        cand = rem / L(i, r);
                    }
                    if (!ub_set || cand < ub) {
                        ub = cand;
                        ub_set = true;
                    }
                }
                if (Z[r] > zero) {
                    const T cand = num_traits<T>::from_int(N[r]) / Z[r];
                    if (!ub_set || cand < ub) {
                        ub = cand;
                        ub_set = true;
                    }
                }
                if (!ub_set) throw InputError("sum_closed: all service demands are zero");

                const T lambda_old = out.XN[r];
                T lambda_l = zero, lambda_u = ub;
                const double ubd = num_traits<T>::to_double(ub);
                const double stop = options.tol * (ubd > 1.0 ? ubd : 1.0) / 1e3;
                while (num_traits<T>::to_double(T(lambda_u - lambda_l)) > stop) {
                    const T lambda = (lambda_l + lambda_u) / two;
                    out.XN[r] = lambda;
                    const Matrix<T> Qir = detail::sum_node_qlen(L, out.XN, mi, scv, K);
                    T g = lambda * Z[r];
                    for (std::size_t i = 0; i < M; ++i) g += Qir(i, r);
                    if (g > num_traits<T>::from_int(N[r]))
                        lambda_u = lambda;
                    else
                        lambda_l = lambda;
                }
                out.XN[r] = (lambda_l + lambda_u) / two;
                const double d = num_traits<T>::to_double(num_abs(T(out.XN[r] - lambda_old)));
                if (d > delta) delta = d;
            }
            if (delta <= options.tol) break;
        }
        if (out.it > options.maxiter) out.it = options.maxiter;
    }

    out.QN = detail::sum_node_qlen(L, out.XN, mi, scv, K);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) {
            out.UN(i, r) = mi[i].infinite
                               ? T(out.XN[r] * L(i, r))
                               : T(out.XN[r] * L(i, r) / num_traits<T>::from_int(mi[i].m));
            if (out.XN[r] > zero) out.RN(i, r) = out.QN(i, r) / out.XN[r];
        }
    return out;
}

}  // namespace sum
}  // namespace line

#endif  // LINE_API_SUM_SUM_CLOSED_H
