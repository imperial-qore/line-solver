/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_PROPFAIR_H
#define LINE_API_PFQN_PFQN_PROPFAIR_H

/**
 * Proportionally fair allocation estimate of the normalizing constant
 * (Schweitzer 1979; Walton, "Proportional fairness and its relationship with
 * multi-class queueing networks", 2009).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_propfair.m. The asymptotic
 * throughput vector solves the convex program
 *
 *   maximize   sum_r (N_r - X_r Z_r) log(X_r + 1e-6)
 *   subject to L X <= 1,  X >= 0
 *
 * after which log G = sum_r (N_r - X_r Z_r) log(1/X_r) - sum_r factln(X_r Z_r).
 * The estimate is asymptotically exact for networks of single-server PS
 * queues; delay stations are handled by the heuristic above.
 *
 * OPTIMIZER. MATLAB calls fmincon. There is no fmincon here, and a generic
 * nonlinear programming solver is not something to invent, so the port solves
 * the same program with a primal log-barrier Newton method: the objective is
 * separable and strictly concave (d^2/dX_r^2 = -2 Z_r/(X_r+eps) -
 * (N_r - X_r Z_r)/(X_r+eps)^2 < 0 wherever the objective is defined), the
 * feasible set is a polytope, so the barrier path is well defined and the
 * method converges to the same maximizer fmincon reports. The centering
 * parameter is raised geometrically until the duality gap (M+R)/t falls below
 * the tolerance. What differs from MATLAB is only the path taken, not the
 * point reached; the accompanying test asserts agreement with the MATLAB
 * value.
 *
 * The reference starts fmincon at the origin, which is on the boundary of the
 * feasible set and outside the domain of the barrier. The port starts at the
 * strictly interior point X_r = 1/(2 max_m sum_s L_ms) instead, which is the
 * only deviation the barrier formulation forces.
 *
 * MATLAB DEAD CODE, noted rather than reproduced: pfqn_propfair.m accumulates
 * a first value of lG in a loop over the classes with Z_r > 0 and then
 * OVERWRITES it on the next line. The loop has no effect on the returned
 * value, so it is not ported.
 *
 * ARITHMETIC. Logarithms throughout, in the objective and in the barrier, so
 * gated on num_traits<T>::has_transcendental.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_propfair, mirroring [G, lG, Xasy]. */
template <class T>
struct PropfairResult {
    T G;
    T lG;
    std::vector<T> Xasy;
};

/**
 * @param L (M x R) demands, @param N (R) population, @param Z (R) think times
 */
template <class T>
PropfairResult<T> pfqn_propfair(const Matrix<T>& L, const std::vector<T>& N,
                                const std::vector<T>& Z) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_propfair requires transcendental arithmetic (logarithmic objective)");
    using std::log;
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_propfair: L and N disagree on the class count");
    if (!Z.empty() && Z.size() != R) throw InputError("pfqn_propfair: Z has the wrong length");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::vector<T> Zv = Z.empty() ? std::vector<T>(R, zero) : Z;
    const T eps = num_traits<T>::from_double(1e-6);  // MATLAB's log(X + 1e-6)

    // Strictly interior start.
    T rowmax = zero;
    for (std::size_t m = 0; m < M; ++m) {
        T s = zero;
        for (std::size_t r = 0; r < R; ++r) s += L(m, r);
        if (s > rowmax) rowmax = s;
    }
    if (rowmax <= zero) throw InputError("pfqn_propfair: every demand is zero");
    std::vector<T> X(R, T(one / T(num_traits<T>::from_int(2) * rowmax)));

    const auto objective = [&](const std::vector<T>& x) {
        T f = zero;
        for (std::size_t r = 0; r < R; ++r) f += T(N[r] - x[r] * Zv[r]) * log(T(x[r] + eps));
        return f;
    };
    const auto feasible = [&](const std::vector<T>& x) {
        for (std::size_t r = 0; r < R; ++r)
            if (!(x[r] > zero)) return false;
        for (std::size_t m = 0; m < M; ++m) {
            T s = zero;
            for (std::size_t r = 0; r < R; ++r) s += L(m, r) * x[r];
            if (!(s < one)) return false;
        }
        return true;
    };

    const T gaptol = num_traits<T>::from_double(1e-10);
    T t = num_traits<T>::from_int(1);
    for (int outer = 0; outer < 60; ++outer) {
        // Newton on the barrier subproblem: maximize t f(X) + barrier(X).
        for (int inner = 0; inner < 100; ++inner) {
            std::vector<T> slack(M);
            for (std::size_t m = 0; m < M; ++m) {
                T s = zero;
                for (std::size_t r = 0; r < R; ++r) s += L(m, r) * X[r];
                slack[m] = T(one - s);
            }
            std::vector<T> grad(R, zero);
            Matrix<T> H(R, R, zero);
            for (std::size_t r = 0; r < R; ++r) {
                const T d = T(X[r] + eps);
                const T num = T(N[r] - X[r] * Zv[r]);
                grad[r] = T(t * T(T(-Zv[r] * log(d)) + T(num / d)));
                H(r, r) = T(t * T(T(-num_traits<T>::from_int(2) * Zv[r] / d) - T(num / T(d * d))));
                // Barrier for X_r > 0.
                grad[r] += T(one / X[r]);
                H(r, r) -= T(one / T(X[r] * X[r]));
            }
            for (std::size_t m = 0; m < M; ++m) {
                for (std::size_t r = 0; r < R; ++r) {
                    grad[r] -= T(L(m, r) / slack[m]);
                    for (std::size_t s = 0; s < R; ++s)
                        H(r, s) -= T(L(m, r) * L(m, s) / T(slack[m] * slack[m]));
                }
            }
            // Newton step solves H dx = -grad; H is negative definite here.
            std::vector<T> rhs(R);
            for (std::size_t r = 0; r < R; ++r) rhs[r] = T(-grad[r]);
            std::vector<T> dx;
            try {
                dx = solve(H, rhs);
            } catch (const NumericError&) {
                break;  // singular Hessian: accept the current iterate
            }
            T dec = zero;  // Newton decrement, -grad' dx
            for (std::size_t r = 0; r < R; ++r) dec -= grad[r] * dx[r];
            if (num_traits<T>::to_double(num_abs(dec)) < 1e-14) break;

            T step = one;
            std::vector<T> Xn(R);
            bool ok = false;
            const T f0 = T(t * objective(X));
            for (int b = 0; b < 80; ++b) {
                for (std::size_t r = 0; r < R; ++r) Xn[r] = T(X[r] + step * dx[r]);
                if (feasible(Xn) && T(t * objective(Xn)) >= f0) {
                    ok = true;
                    break;
                }
                step = T(step / num_traits<T>::from_int(2));
            }
            if (!ok) break;
            X = Xn;
        }
        const T gap = T(num_traits<T>::from_int(static_cast<long>(M + R)) / t);
        if (gap < gaptol) break;
        t = T(t * num_traits<T>::from_int(10));
    }

    PropfairResult<T> res;
    res.Xasy = X;
    T lG = zero;
    for (std::size_t r = 0; r < R; ++r) {
        if (X[r] <= zero) throw NumericError("pfqn_propfair: non-positive asymptotic throughput");
        lG += T(N[r] - X[r] * Zv[r]) * log(T(one / X[r]));
        lG -= detail::num_factln<T>(T(X[r] * Zv[r]));
    }
    using std::exp;
    res.lG = lG;
    res.G = exp(lG);
    return res;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_PROPFAIR_H
