/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_ISM_GREEN_H
#define LINE_API_FJ_ISM_GREEN_H

/**
 * Green's independent server model of simultaneous server requests.
 *
 * Templated port of matlab/src/api/fj/fj_ism_green.m.
 *
 * A customer needs j servers at once with probability c(j) and then releases
 * them asynchronously as each of its j tasks completes at rate mu. Servers can
 * idle while customers wait, which is what separates the model from M/G/s, and
 * customer service ends with the last of its tasks, so its mean is H_j/mu.
 *
 *   E[B] = sum_j c(j) sum_{i=0..j-1} 1/((s-i) mu)
 *
 * is the interservice time, the j-th order statistic of s exponentials because
 * all s servers are busy whenever a customer enters service in a queueing
 * period, and
 *
 *   E[D] = sum_i sum_{k=1..i} [ sum_{m=0..k-1} 1/((i-m) mu) ] q(i) c(s-i+k)/p_d
 *
 * the initial delay of the customer that starts one. The busy-server
 * distribution q and the nonqueue length come from the embedded chain absorbed
 * when a queue forms, V = (I-T)^-1. The waiting-time transform of Eq. (61)
 * factors into the equilibrium transform of D and the Pollaczek-Khinchine
 * transform of an M/G/1 queue with service B, so
 *
 *   E[W] = (1-pi0) [ E[D^2]/(2 E[D]) + lambda E[B^2]/(2 (1-lambda E[B])) ].
 *
 * Eq. (65) of the survey prints the inner sum as starting at 1/(s mu) even
 * though only i servers are busy; it is started at 1/(i mu) here, which is what
 * the accompanying text prescribes and what makes E[D] reduce to E[B] at i = s.
 */

#include <cstddef>
#include <vector>

#include "line/api/fj/fj_harmonic.h"
#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace fj {

/** Everything Green's cycle decomposition produces. */
template <class T>
struct FJIsmGreenResult {
    T W;
    T R;
    T EB;
    T EB2;
    T ED;
    T ED2;
    T EQ;
    T EQbar;
    T pq;
    T pd;
    T pi0;
    T rho;
    T ES;
    std::vector<T> q;
};

/**
 * @param lambda customer arrival rate
 * @param mu     per-task service rate
 * @param s      number of servers
 * @param c      c[j-1] = P(a customer needs j servers), summing to one
 * @return       the waiting and response times with the cycle quantities
 */
template <class T>
FJIsmGreenResult<T> fj_ism_green(const T& lambda, const T& mu, unsigned s,
                                 const std::vector<T>& c) {
    if (s < 1) throw InputError("fj_ism_green: s must be a positive integer");
    if (c.size() != s)
        throw InputError("fj_ism_green: c must have one entry per server requirement");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1),
            two = num_traits<T>::from_int(2);
    if (!(lambda > zero) || !(mu > zero))
        throw InputError("fj_ism_green: lambda and mu must be positive");
    T tot = zero;
    for (std::size_t i = 0; i < s; ++i) {
        if (c[i] < zero)
            throw InputError("fj_ism_green: the server-requirement probabilities must be non-negative");
        tot += c[i];
    }
    if (!(tot - one < num_traits<T>::from_double(1e-9)) ||
        !(one - tot < num_traits<T>::from_double(1e-9)))
        throw InputError("fj_ism_green: the server-requirement probabilities must sum to one");

    FJIsmGreenResult<T> out;

    // Interservice time B: the j-th order statistic of s exponentials of rate mu
    out.EB = zero;
    out.EB2 = zero;
    for (unsigned j = 1; j <= s; ++j) {
        T mj = zero, vj = zero;
        for (unsigned i = 0; i < j; ++i) {
            const T st = one / (num_traits<T>::from_int(static_cast<long>(s - i)) * mu);
            mj += st;
            vj += st * st;
        }
        out.EB += c[j - 1] * mj;
        out.EB2 += c[j - 1] * (vj + mj * mj);
    }

    out.rho = lambda * out.EB;
    if (out.rho >= one)
        throw NumericError("fj_ism_green: unstable system, rho = lambda*E[B] >= 1");

    // Embedded chain over the busy-server count, absorbed when an arrival needs
    // more servers than are free
    const std::size_t n = s + 1;
    Matrix<T> A(n, n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        const T den = lambda + num_traits<T>::from_int(static_cast<long>(i)) * mu;
        if (i > 0) A(i, i - 1) = num_traits<T>::from_int(static_cast<long>(i)) * mu / den;
        for (std::size_t j = 1; j + i <= s; ++j)
            A(i, i + j) = A(i, i + j) + lambda * c[j - 1] / den;
    }
    // V = (I - T)^-1; only the row that starts with all s servers busy is used
    Matrix<T> M(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) M(i, j) = (i == j ? one : zero) - A(i, j);
    // Gaussian elimination on M^T x = e_s, which is the row s of (I-T)^-1
    std::vector<std::vector<T> > G(n, std::vector<T>(n + 1, zero));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) G[i][j] = M(j, i);
        G[i][n] = (i == s) ? one : zero;
    }
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t piv = col;
        for (std::size_t r = col + 1; r < n; ++r) {
            const T a = G[r][col] > zero ? G[r][col] : -G[r][col];
            const T b = G[piv][col] > zero ? G[piv][col] : -G[piv][col];
            if (a > b) piv = r;
        }
        if (G[piv][col] == zero)
            throw NumericError("fj_ism_green: the embedded chain is singular");
        G[col].swap(G[piv]);
        for (std::size_t r = 0; r < n; ++r) {
            if (r == col) continue;
            const T f = G[r][col] / G[col][col];
            for (std::size_t j = col; j <= n; ++j) G[r][j] -= f * G[col][j];
        }
    }
    std::vector<T> v(n);
    for (std::size_t i = 0; i < n; ++i) v[i] = G[i][n] / G[i][i];

    out.EQbar = zero;
    std::vector<T> hold(n);
    for (std::size_t i = 0; i < n; ++i) {
        hold[i] = one / (lambda + num_traits<T>::from_int(static_cast<long>(i)) * mu);
        out.EQbar += v[i] * hold[i];
    }
    out.q.resize(n);
    for (std::size_t i = 0; i < n; ++i) out.q[i] = v[i] * hold[i] / out.EQbar;

    // A customer arriving during a nonqueue period is delayed when it needs more
    // servers than the s-i free ones
    out.pd = zero;
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t free = s - i;
        for (std::size_t j = free + 1; j <= s; ++j) out.pd += out.q[i] * c[j - 1];
    }
    if (!(out.pd > zero))
        throw NumericError("fj_ism_green: no arrival can ever be delayed; the model is M/M/s");

    // Initial delay D: i busy, the customer needs k = j-(s-i) more to free
    out.ED = zero;
    out.ED2 = zero;
    for (unsigned i = 1; i <= s; ++i) {
        for (unsigned k = 1; k <= i; ++k) {
            const unsigned j = s - i + k;
            if (j < 1 || j > s) continue;
            const T wgt = out.q[i] * c[j - 1] / out.pd;
            if (wgt == zero) continue;
            T mk = zero, vk = zero;
            for (unsigned m = 0; m < k; ++m) {
                const T st = one / (num_traits<T>::from_int(static_cast<long>(i - m)) * mu);
                mk += st;
                vk += st * st;
            }
            out.ED += wgt * mk;
            out.ED2 += wgt * (vk + mk * mk);
        }
    }

    out.EQ = out.ED / (one - out.rho);
    out.pq = out.EQ / (out.EQ + out.EQbar);
    out.pi0 = (one - out.rho) / (one - lambda * (out.EB - out.ED));

    const T Weq = out.ED2 / (two * out.ED);
    const T Wmg1 = lambda * out.EB2 / (two * (one - out.rho));
    out.W = (one - out.pi0) * (Weq + Wmg1);

    // Customer service time: the maximum of the j tasks it holds
    out.ES = zero;
    for (unsigned j = 1; j <= s; ++j) out.ES += c[j - 1] * fj_harmonic<T>(j) / mu;
    out.R = out.W + out.ES;
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_ISM_GREEN_H
