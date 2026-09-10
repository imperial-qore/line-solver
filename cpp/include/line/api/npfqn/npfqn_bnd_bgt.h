/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_NPFQN_NPFQN_BND_BGT_H
#define LINE_API_NPFQN_NPFQN_BND_BGT_H

/**
 * Piecewise-linear Lyapunov UPPER bound on the steady-state queue lengths of a
 * multitype (deterministic-routing) multiclass Markovian queueing network,
 * valid for EVERY work-conserving Markovian policy.
 *
 * Templated port of matlab/src/api/npfqn/npfqn_bnd_bgt.m, cross-checked against
 * jar/src/main/java/jline/api/npfqn/Npfqn_bnd_bgt.java.
 *
 * MODEL. J single-server stations; I customer types; type i arrives as a
 * Poisson stream of rate lambda(i) and passes through stages k = 0..Ji-1, stage
 * k being served at station sigma[i][k] at exponential rate mu[i][k]. Class
 * (i,k) is the buffer of type i at stage k; N = sum_i Ji is the class count.
 *
 * METHOD. Solve the Down-Meyn global-stability linear program GLP[dm], eq.
 * (25)-(28) of the reference, in the piecewise-linear Lyapunov function
 * Phi(x) = max_j L^j'x:
 *
 *   L^j(i,1) lambda_i + mu(i,k) (L^j(i,k+1) - L^j(i,k)) + V_j <= -gamma
 *                                                    for (i,k) in station j
 *   mu(i,k) (L^j(i,k+1) - L^j(i,k)) <= V_j           for (i,k) not in j
 *   (1/(J-1)) sum_{j' != j} L^j'(i,k) >= L^j(i,k)    for (i,k) not in j
 *   L, V, gamma >= 0
 *
 * with L^j(i,Ji+1) = 0. A feasible solution with gamma > 0 certifies that EVERY
 * work-conserving policy is stable, and a smoothed Phi is then a Lyapunov
 * function with drift gamma/4 and an explicit exception parameter, giving the
 * reference's Theorem 4 bound
 *
 *   E[L^j'Q] <= 16 N J^2 (J-1) (Lmax+gamma)^3/gamma^2
 *               + 8 (Lmax + gamma/2)^2/gamma  =: U
 *
 * for every j, whence E[Q(i,k)] <= U / max_j L^j(i,k).
 *
 * THE RATES ARE RESCALED so that sum_i lambda_i + sum_{i,k} mu(i,k) = 1, the
 * uniformization the reference imposes before Theorem 4. Queue lengths are
 * counts and are unaffected by the time scale.
 *
 * NORMALIZATION, WHICH THE REFERENCE LEAVES OPEN. GLP[dm] is homogeneous and so
 * is the bound, so this routine fixes L^j(i,k) <= 1 and MAXIMIZES gamma, then
 * breaks ties among gamma-optimal solutions by maximizing sum L: a degenerate
 * optimum can otherwise zero some L^j(i,k) and report an infinite bound for a
 * class for no reason.
 *
 * THE BOUND IS LOOSE, and knowingly so: the exception parameter carries
 * (Lmax+gamma)^3/gamma^2 and dominates as soon as J > 1. What is sharp is the
 * STABILITY CERTIFICATE gamma > 0 and the geometric tail RATE.
 *
 * ARITHMETIC. Rational-clean: the LP data and the bound are polynomial in the
 * rates and the dense simplex is exact.
 *
 * Reference: D. Bertsimas, D. Gamarnik, J. N. Tsitsiklis (2001). Performance of
 * multiclass Markovian queueing networks via piecewise linear Lyapunov
 * functions. Annals of Applied Probability 11(4), 1384-1428, Section 5.1
 * (GLP[dm] of Down and Meyn 1997, and Theorem 4).
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lp_highs.h"
#include "line/util/simplex.h"

namespace line {
namespace npfqn {

template <class T>
struct BndBgt {
    std::vector<std::vector<T> > Qub;  ///< per type and stage, the bound on E[Q(i,k)]
    std::vector<bool> finite;          ///< false where max_j L is 0 and Qub is meaningless
    T gamma = T();                     ///< the drift certificate, strictly positive on success
    T Lmax = T();                      ///< max over j and (i,k) of L
    std::vector<std::vector<T> > L;    ///< the Lyapunov coefficients, [J][N]
    std::vector<T> V;                  ///< the per-station slack
    T B = T();                         ///< exception parameter of the smoothed function
    T U = T();                         ///< the Theorem 4 bound on E[L^j'Q]
    T tailRatio = T();                 ///< geometric decay ratio of the tail bound
    T tailStep = T();                  ///< step of the tail bound, 2(Lmax+gamma/2)
    std::vector<T> rho;                ///< per-class nominal load
    std::vector<T> rhoStation;         ///< per-station nominal load
    T scale = T();                     ///< the uniformization divisor
    std::vector<std::size_t> classType, classStage, classStation;
};

/**
 * @param lambda Poisson arrival rate of each type
 * @param mu     mu[i][k] = service rate of stage k of type i
 * @param sigma  sigma[i][k] = zero-based station of stage k of type i
 * @param J      number of stations
 */
template <class T>
BndBgt<T> npfqn_bnd_bgt(const std::vector<T>& lambda, const std::vector<std::vector<T> >& mu,
                        const std::vector<std::vector<std::size_t> >& sigma, std::size_t J) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2), four = num_traits<T>::from_int(4);
    const std::size_t I = lambda.size();
    if (mu.size() != I || sigma.size() != I)
        throw InputError("npfqn_bnd_bgt: mu and sigma must have one entry per type");

    // ---- flatten (i,k) into a class index ----
    BndBgt<T> out;
    std::vector<T> muc;
    std::vector<std::size_t> firstOf(I, 0), nextOf;
    for (std::size_t i = 0; i < I; ++i) {
        if (mu[i].size() != sigma[i].size())
            throw InputError("npfqn_bnd_bgt: mu and sigma disagree on the stage count of a type");
        if (mu[i].empty()) throw InputError("npfqn_bnd_bgt: a type has no stage");
        if (!(lambda[i] > zero))
            throw InputError("npfqn_bnd_bgt: every type needs a strictly positive arrival rate");
        firstOf[i] = muc.size();
        for (std::size_t k = 0; k < mu[i].size(); ++k) {
            out.classType.push_back(i);
            out.classStage.push_back(k);
            out.classStation.push_back(sigma[i][k]);
            muc.push_back(mu[i][k]);
        }
    }
    const std::size_t N = muc.size();
    for (std::size_t c = 0; c < N; ++c)
        if (!(muc[c] > zero))
            throw InputError("npfqn_bnd_bgt: every stage needs a strictly positive service rate");
    nextOf.assign(N, N);  // N means "leaves the network"
    for (std::size_t c = 0; c + 1 < N; ++c)
        if (out.classType[c + 1] == out.classType[c]) nextOf[c] = c + 1;

    // ---- loads ----
    out.rho.assign(N, zero);
    out.rhoStation.assign(J, zero);
    for (std::size_t c = 0; c < N; ++c) {
        out.rho[c] = T(lambda[out.classType[c]] / muc[c]);
        out.rhoStation[out.classStation[c]] = T(out.rhoStation[out.classStation[c]] + out.rho[c]);
    }
    for (std::size_t j = 0; j < J; ++j)
        if (!(out.rhoStation[j] < one))
            throw UnsupportedError("npfqn_bnd_bgt: station " + std::to_string(j + 1) +
                                   " is saturated: the load condition of the reference fails");

    // ---- uniformization ----
    T scale = zero;
    for (std::size_t i = 0; i < I; ++i) scale = T(scale + lambda[i]);
    for (std::size_t c = 0; c < N; ++c) scale = T(scale + muc[c]);
    out.scale = scale;
    std::vector<T> lam(I, zero), mus(N, zero);
    for (std::size_t i = 0; i < I; ++i) lam[i] = T(lambda[i] / scale);
    for (std::size_t c = 0; c < N; ++c) mus[c] = T(muc[c] / scale);

    // ---- LP layout: L(j,c) -> j*N + c ; V(j) -> J*N + j ; gamma -> J*N+J ----
    const std::size_t oV = J * N, ig = J * N + J, nv = J * N + J + 1;
    lp::LpModel<T> m(nv);
    m.set_maximize(true);
    for (std::size_t v = 0; v < J * N; ++v) m.set_bounds(v, zero, one);

    for (std::size_t j = 0; j < J; ++j) {
        for (std::size_t c = 0; c < N; ++c) {
            m.row_clear();
            if (out.classStation[c] == j) {
                m.row_add(j * N + firstOf[out.classType[c]], lam[out.classType[c]]);
                m.row_add(j * N + c, T(-mus[c]));
                if (nextOf[c] < N) m.row_add(j * N + nextOf[c], mus[c]);
                m.row_add(oV + j, one);
                m.row_add(ig, one);
                m.emit(lp::LpSense::LE, zero);
            } else {
                m.row_add(j * N + c, T(-mus[c]));
                if (nextOf[c] < N) m.row_add(j * N + nextOf[c], mus[c]);
                m.row_add(oV + j, T(-one));
                m.emit(lp::LpSense::LE, zero);
                if (J > 1) {
                    m.row_clear();
                    m.row_add(j * N + c, one);
                    const T w = T(-(one / num_traits<T>::from_int(static_cast<int>(J - 1))));
                    for (std::size_t jp = 0; jp < J; ++jp)
                        if (jp != j) m.row_add(jp * N + c, w);
                    m.emit(lp::LpSense::LE, zero);
                }
            }
        }
    }

    m.set_cost(ig, one);
    lp::LpSolution<T> s = lp::lp_solve(m);
    if (!s.ok())
        throw UnsupportedError(std::string("npfqn_bnd_bgt: GLP[dm] did not solve to optimality (") +
                               lp::lp_status_name(s.status) + ")");
    T gamma = s.objective;
    if (!(gamma > zero))
        throw UnsupportedError(
            "npfqn_bnd_bgt: GLP[dm] has no solution with gamma > 0: this network is not certified "
            "globally stable, so no finite piecewise-linear Lyapunov bound exists");

    // Tie-break among gamma-optimal solutions: maximize sum L, so a degenerate
    // vertex does not report an infinite bound for a class it zeroed arbitrarily.
    {
        lp::LpModel<T> m2 = m;
        m2.row_clear();
        m2.row_add(ig, one);
        m2.emit(lp::LpSense::GE, gamma);
        m2.set_cost(ig, zero);
        for (std::size_t v = 0; v < J * N; ++v) m2.set_cost(v, one);
        const lp::LpSolution<T> s2 = lp::lp_solve(m2);
        if (s2.ok()) {
            s = s2;
            gamma = s.x[ig];
        }
    }

    out.gamma = gamma;
    out.L.assign(J, std::vector<T>(N, zero));
    out.Lmax = zero;
    for (std::size_t j = 0; j < J; ++j)
        for (std::size_t c = 0; c < N; ++c) {
            out.L[j][c] = s.x[j * N + c];
            if (out.L[j][c] > out.Lmax) out.Lmax = out.L[j][c];
        }
    out.V.assign(J, zero);
    for (std::size_t j = 0; j < J; ++j) out.V[j] = s.x[oV + j];

    const T Lg = T(out.Lmax + gamma);
    const T Lh = T(out.Lmax + gamma / two);
    out.B = T(num_traits<T>::from_int(16) * num_traits<T>::from_int(static_cast<int>(N)) *
              num_traits<T>::from_int(static_cast<int>(J * J)) *
              num_traits<T>::from_int(static_cast<int>(J - 1)) * Lg * Lg * Lg / (gamma * gamma));
    out.U = T(out.B + num_traits<T>::from_int(8) * Lh * Lh / gamma);
    out.tailStep = T(two * Lh);
    out.tailRatio = T(Lh / (out.Lmax + num_traits<T>::from_int(3) * gamma / four));

    out.Qub.assign(I, std::vector<T>());
    out.finite.assign(N, true);
    for (std::size_t i = 0; i < I; ++i) out.Qub[i].assign(mu[i].size(), zero);
    for (std::size_t c = 0; c < N; ++c) {
        T best = zero;
        for (std::size_t j = 0; j < J; ++j)
            if (out.L[j][c] > best) best = out.L[j][c];
        if (best > zero) {
            out.Qub[out.classType[c]][out.classStage[c]] = T(out.U / best);
        } else {
            out.finite[c] = false;
            out.Qub[out.classType[c]][out.classStage[c]] = zero;
        }
    }
    return out;
}

}  // namespace npfqn
}  // namespace line

#endif  // LINE_API_NPFQN_NPFQN_BND_BGT_H
