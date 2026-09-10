/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_SENS_RESPT_H
#define LINE_API_PFQN_SENS_RESPT_H

/**
 * Exact raw moments E[W^t], t = 1..3, of the sojourn time of a job at an FCFS
 * b-server center of a closed product-form network.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_sens_respt.m, Theorem 4.1 of
 * Strelen (Performance Evaluation 11:127-142, 1990). By the arrival theorem a
 * class-l job finds j jobs at center i with probability p_i(j, N - e_l);
 * conditioning on j and inverting the Laplace transform of the conditional
 * density gives equation (4.5),
 *
 *   E[W_(i,l)^t] = t!/mu^t + sum_tau a_(t,tau)(0) E[Qt_i^tau]
 *                  - sum_{j<b_i} p_i(j,N-e_l) sum_tau a_(t,tau)(0) j^tau
 *
 * with mu = 1/S(i), Qt_i the total queue at center i at population N - e_l and
 * the coefficients a_(t,tau)(0) of Remark 4.3 depending only on b and mu. The
 * moments E[Qt_i^tau] up to tau = 3 need the second derivative of the
 * b-server recursion (4.1)-(4.2), which is carried here by a second-order
 * forward-mode pass along the population lattice, exactly as pfqn_sens_mom
 * does for the single-server recursion.
 *
 * Only FCFS centers are covered: the reference is explicit that the
 * sojourn-time distribution at PS and LCFS centers is in general not known.
 * FCFS in a BCMP network requires a class-independent exponential service
 * time, which is why the input is a per-station service time S and a separate
 * visit matrix V rather than a demand matrix: the sojourn time is per visit, so
 * mu = 1/S(i) must be known and cannot be recovered from L(i,l) = S(i) V(i,l).
 *
 * Arithmetic. Every step is a field operation, so the moments instantiate at
 * line::Rational and are exact rationals. The two skewnesses are the only
 * derived quantities that leave the field (they divide by a variance to the
 * power 1.5), so they are returned as doubles and the header carries no
 * transcendental gate.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_sens_mva.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

template <class T>
struct SensResptResult {
    std::vector<T> XN;  ///< (R) throughput at population N
    Matrix<T> QN;       ///< (M x R) mean queue length
    Matrix<T> UN;       ///< (M x R) utilization

    Matrix<T> W;                ///< (M x R) E[W_(i,l)], the mean sojourn time per visit
    std::vector<Matrix<T>> WM;  ///< (tmax) matrices M x R, WM[t-1](i,l) = E[W^t]
    Matrix<T> WVar;             ///< (M x R) Var[W], zero unless tmax >= 2
    Matrix<double> WSkew;       ///< (M x R) skewness of W, zero unless tmax >= 3

    std::vector<T> m;    ///< (M) E[Q_i] at population N
    std::vector<T> Var;  ///< (M) Var[Q_i] at population N
    /// (M x max(b)) p(i,j) = P[Q_i = j] for j = 0..b_i-1. RAGGED: row i is
    /// meaningful only up to column b_i-1 and is zero-padded out to max(b).
    Matrix<T> p;
    Matrix<T> Wresid;  ///< (M x R) the residence time w_i(l) of the recursion
};

/**
 * @param S    (M) service time of each station, common to all classes
 * @param V    (M x R) visit ratios; the demand is L(i,r) = S(i) V(i,r)
 * @param N    (R) population per class, closed only
 * @param Z    (R) think times, empty for none
 * @param b    (M) servers per station, empty for all ones
 * @param tmax highest sojourn-time moment, 1..3
 */
template <class T>
SensResptResult<T> pfqn_sens_respt(const std::vector<T>& S, const Matrix<T>& V,
                                   const std::vector<int>& N, const std::vector<T>& Z,
                                   const std::vector<int>& b, int tmax) {
    const std::size_t M = V.rows();
    const std::size_t R = N.size();
    if (V.cols() != R)
        throw InputError("pfqn_sens_respt: visit matrix and population vector disagree on the class count");
    if (S.size() != M) throw InputError("pfqn_sens_respt: service-time vector has the wrong length");
    if (!Z.empty() && Z.size() != R)
        throw InputError("pfqn_sens_respt: think-time vector has the wrong length");
    if (!b.empty() && b.size() != M)
        throw InputError("pfqn_sens_respt: server-count vector has the wrong length");
    if (tmax < 1 || tmax > 3)
        throw InputError(
            "pfqn_sens_respt: tmax must be 1, 2 or 3, the orders at which the coefficients "
            "a_{t,tau}(0) are tabulated in the reference");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    std::vector<int> bv = b.empty() ? std::vector<int>(M, 1) : b;
    std::size_t bmax = 1;
    for (std::size_t i = 0; i < M; ++i) {
        if (bv[i] < 1)
            throw InputError("pfqn_sens_respt: every station must have at least one server");
        if (static_cast<std::size_t>(bv[i]) > bmax) bmax = static_cast<std::size_t>(bv[i]);
        if (S[i] <= zero)
            throw InputError("pfqn_sens_respt: every FCFS station needs a strictly positive service time");
    }

    SensResptResult<T> res;
    res.XN.assign(R, zero);
    res.QN = Matrix<T>(M, R, zero);
    res.UN = Matrix<T>(M, R, zero);
    res.W = Matrix<T>(M, R, zero);
    res.WM.assign(static_cast<std::size_t>(tmax), Matrix<T>(M, R, zero));
    res.WVar = Matrix<T>(M, R, zero);
    res.WSkew = Matrix<double>(M, R, 0.0);
    res.m.assign(M, zero);
    res.Var.assign(M, zero);
    res.p = Matrix<T>(M, bmax, zero);
    res.Wresid = Matrix<T>(M, R, zero);

    bool anyPositive = false;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_sens_respt: negative population");
        if (v > 0) anyPositive = true;
    }
    if (!anyPositive || M == 0 || R == 0) return res;

    Matrix<T> rho(M, R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t l = 0; l < R; ++l) rho(i, l) = S[i] * V(i, l);
    std::vector<T> mu(M, zero);
    for (std::size_t i = 0; i < M; ++i) mu[i] = one / S[i];

    const auto Zr = [&](std::size_t r) -> T { return Z.empty() ? zero : Z[r]; };

    const std::vector<std::size_t> radix = sens_lattice_radix(N);
    const std::size_t totpop = population_count(N);

    // lattice-state rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    Matrix<T> Mrow(totpop, M, zero);
    std::vector<Matrix<T>> D1m(totpop, Matrix<T>(M, M, zero));
    std::vector<Matrix<T>> D2m(totpop, Matrix<T>(M, M, zero));
    std::vector<Matrix<T>> Prow(totpop, Matrix<T>(M, bmax, zero));
    std::vector<std::vector<Matrix<T>>> D1p(totpop, std::vector<Matrix<T>>(M, Matrix<T>(bmax, M, zero)));
    std::vector<std::vector<Matrix<T>>> D2p(totpop, std::vector<Matrix<T>>(M, Matrix<T>(bmax, M, zero)));
    for (std::size_t i = 0; i < M; ++i) Prow[0](i, 0) = one;  // empty population

    std::vector<std::size_t> rows(R, 0);
    Matrix<T> wv(M, R, zero);
    std::vector<Matrix<T>> d1w(M, Matrix<T>(R, M, zero)), d2w(M, Matrix<T>(R, M, zero));
    std::vector<T> lam(R, zero);
    Matrix<T> d1lam(R, M, zero), d2lam(R, M, zero);

    for (std::size_t k = 1; k < totpop; ++k) {
        const std::vector<int> n = sens_lattice_decode(k, N, radix);
        wv.fill(zero);
        for (std::size_t i = 0; i < M; ++i) {
            d1w[i].fill(zero);
            d2w[i].fill(zero);
        }
        // ---- residence times (4.1) and their first two derivatives --------
        for (std::size_t s = 0; s < R; ++s) {
            const std::size_t row = n[s] > 0 ? k - radix[s] : 0;
            rows[s] = row;
            if (n[s] == 0) continue;  // w and every derivative stay zero, as does X(s)
            for (std::size_t i = 0; i < M; ++i) {
                const T bi = num_traits<T>::from_int(bv[i]);
                T brk = one + Mrow(row, i);
                for (int j = 0; j + 2 <= bv[i]; ++j)
                    brk += num_traits<T>::from_int(bv[i] - 1 - j) * Prow[row](i, static_cast<std::size_t>(j));
                wv(i, s) = rho(i, s) / bi * brk;
                for (std::size_t h = 0; h < M; ++h) {
                    T dbrk = D1m[row](i, h);
                    T d2brk = D2m[row](i, h);
                    for (int j = 0; j + 2 <= bv[i]; ++j) {
                        const T c = num_traits<T>::from_int(bv[i] - 1 - j);
                        dbrk += c * D1p[row][i](static_cast<std::size_t>(j), h);
                        d2brk += c * D2p[row][i](static_cast<std::size_t>(j), h);
                    }
                    if (i == h) {
                        d1w[i](s, h) = rho(i, s) / bi * (brk + dbrk);
                        d2w[i](s, h) = rho(i, s) / bi * (num_traits<T>::from_int(2) * dbrk + d2brk);
                    } else {
                        d1w[i](s, h) = rho(i, s) / bi * dbrk;
                        d2w[i](s, h) = rho(i, s) / bi * d2brk;
                    }
                }
            }
        }

        // ---- throughputs and their derivatives ------------------------------
        for (std::size_t s = 0; s < R; ++s) {
            lam[s] = zero;
            for (std::size_t h = 0; h < M; ++h) {
                d1lam(s, h) = zero;
                d2lam(s, h) = zero;
            }
            if (n[s] == 0) continue;
            T sw = zero;
            for (std::size_t i = 0; i < M; ++i) sw += wv(i, s);
            const T den = Zr(s) + sw;
            const T nsT = num_traits<T>::from_int(n[s]);
            lam[s] = nsT / den;
            const T den2 = den * den;
            const T den3 = den2 * den;
            for (std::size_t h = 0; h < M; ++h) {
                T dden = zero, d2den = zero;
                for (std::size_t i = 0; i < M; ++i) {
                    dden += d1w[i](s, h);
                    d2den += d2w[i](s, h);
                }
                d1lam(s, h) = -nsT * dden / den2;
                d2lam(s, h) = -nsT * d2den / den2 +
                              num_traits<T>::from_int(2) * nsT * dden * dden / den3;
            }
        }

        // ---- mean queue lengths ----------------------------------------------
        for (std::size_t i = 0; i < M; ++i) {
            T acc = zero;
            for (std::size_t s = 0; s < R; ++s)
                if (n[s] > 0) acc += lam[s] * wv(i, s);
            Mrow(k, i) = acc;
            for (std::size_t h = 0; h < M; ++h) {
                T d1acc = zero, d2acc = zero;
                for (std::size_t s = 0; s < R; ++s) {
                    if (n[s] == 0) continue;
                    d1acc += d1lam(s, h) * wv(i, s) + lam[s] * d1w[i](s, h);
                    d2acc += d2lam(s, h) * wv(i, s) +
                             num_traits<T>::from_int(2) * d1lam(s, h) * d1w[i](s, h) +
                             lam[s] * d2w[i](s, h);
                }
                D1m[k](i, h) = d1acc;
                D2m[k](i, h) = d2acc;
            }
        }

        // ---- marginal probabilities (4.2) and their derivatives ---------------
        int nc = 0;
        for (std::size_t s = 0; s < R; ++s) nc += n[s];
        for (std::size_t i = 0; i < M; ++i) {
            for (int j = 1; j <= bv[i] - 1; ++j) {
                const std::size_t ju = static_cast<std::size_t>(j);
                if (j > nc) continue;  // stays zero
                const T jT = num_traits<T>::from_int(j);
                T acc = zero;
                for (std::size_t l = 0; l < R; ++l) {
                    if (n[l] == 0) continue;
                    acc += lam[l] * rho(i, l) * Prow[rows[l]](i, ju - 1);
                }
                Prow[k](i, ju) = acc / jT;
                for (std::size_t h = 0; h < M; ++h) {
                    T d1acc = zero, d2acc = zero;
                    for (std::size_t l = 0; l < R; ++l) {
                        if (n[l] == 0) continue;
                        const T pprev = Prow[rows[l]](i, ju - 1);
                        const T d1prev = D1p[rows[l]][i](ju - 1, h);
                        const T d2prev = D2p[rows[l]][i](ju - 1, h);
                        T v1, v2;
                        if (i == h) {
                            v1 = pprev + d1prev;
                            v2 = num_traits<T>::from_int(2) * d1prev + d2prev;
                        } else {
                            v1 = d1prev;
                            v2 = d2prev;
                        }
                        d1acc += rho(i, l) * (d1lam(l, h) * pprev + lam[l] * v1);
                        d2acc += rho(i, l) * (d2lam(l, h) * pprev +
                                              num_traits<T>::from_int(2) * d1lam(l, h) * v1 +
                                              lam[l] * v2);
                    }
                    D1p[k][i](ju, h) = d1acc / jT;
                    D2p[k][i](ju, h) = d2acc / jT;
                }
            }
            // mean number of busy servers
            T ui = zero;
            std::vector<T> d1ui(M, zero), d2ui(M, zero);
            for (std::size_t l = 0; l < R; ++l) {
                if (n[l] == 0) continue;
                ui += lam[l] * rho(i, l);
                for (std::size_t h = 0; h < M; ++h) {
                    if (i == h) {
                        d1ui[h] += rho(i, l) * (d1lam(l, h) + lam[l]);
                        d2ui[h] += rho(i, l) * (d2lam(l, h) + num_traits<T>::from_int(2) * d1lam(l, h));
                    } else {
                        d1ui[h] += rho(i, l) * d1lam(l, h);
                        d2ui[h] += rho(i, l) * d2lam(l, h);
                    }
                }
            }
            const T biT = num_traits<T>::from_int(bv[i]);
            T acc0 = ui;
            for (int j = 1; j <= bv[i] - 1; ++j)
                acc0 += num_traits<T>::from_int(bv[i] - j) * Prow[k](i, static_cast<std::size_t>(j));
            Prow[k](i, 0) = one - acc0 / biT;
            for (std::size_t h = 0; h < M; ++h) {
                T d1acc0 = d1ui[h], d2acc0 = d2ui[h];
                for (int j = 1; j <= bv[i] - 1; ++j) {
                    const T c = num_traits<T>::from_int(bv[i] - j);
                    d1acc0 += c * D1p[k][i](static_cast<std::size_t>(j), h);
                    d2acc0 += c * D2p[k][i](static_cast<std::size_t>(j), h);
                }
                D1p[k][i](0, h) = -d1acc0 / biT;
                D2p[k][i](0, h) = -d2acc0 / biT;
            }
        }

        // keep the measures of the last (full) population
        for (std::size_t s = 0; s < R; ++s) res.XN[s] = lam[s];
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t s = 0; s < R; ++s) {
                res.Wresid(i, s) = wv(i, s);
                res.QN(i, s) = lam[s] * wv(i, s);
                res.UN(i, s) = lam[s] * rho(i, s);
            }
    }

    const std::size_t last = totpop - 1;
    for (std::size_t i = 0; i < M; ++i) {
        res.m[i] = Mrow(last, i);
        for (std::size_t j = 0; j < bmax; ++j) res.p(i, j) = Prow[last](i, j);
        res.Var[i] = D1m[last](i, i);
    }

    // ---- sojourn-time moments (4.5) ---------------------------------------
    std::vector<std::size_t> rowsN(R, 0);
    for (std::size_t l = 0; l < R; ++l)
        if (N[l] > 0) rowsN[l] = last - radix[l];

    for (std::size_t i = 0; i < M; ++i) {
        // a_{t,tau}(0) of Remark 4.3: they depend only on b and mu.
        const T bT = num_traits<T>::from_int(bv[i]);
        const T mT = mu[i];
        Matrix<T> a(3, 4, zero);
        a(0, 0) = (one - bT) / (bT * mT);
        a(0, 1) = one / (bT * mT);
        if (tmax >= 2) {
            const T d = bT * bT * mT * mT;
            a(1, 0) = (num_traits<T>::from_int(2) - bT - bT * bT) / d;
            a(1, 1) = num_traits<T>::from_int(3) / d;
            a(1, 2) = one / d;
        }
        if (tmax >= 3) {
            const T d = bT * bT * bT * mT * mT * mT;
            a(2, 0) = (num_traits<T>::from_int(6) - num_traits<T>::from_int(5) * bT +
                       num_traits<T>::from_int(3) * bT * bT - num_traits<T>::from_int(4) * bT * bT * bT) / d;
            a(2, 1) = (num_traits<T>::from_int(11) - num_traits<T>::from_int(3) * bT +
                       num_traits<T>::from_int(3) * bT * bT) / d;
            a(2, 2) = num_traits<T>::from_int(6) / d;
            a(2, 3) = one / d;
        }
        for (std::size_t l = 0; l < R; ++l) {
            if (N[l] == 0 || V(i, l) <= zero) continue;
            const std::size_t rl = rowsN[l];
            const T mt = Mrow(rl, i);
            const T d1t = D1m[rl](i, i);
            const T d2t = D2m[rl](i, i);
            T EQ[4];
            EQ[0] = one;
            EQ[1] = mt;
            EQ[2] = d1t + mt * mt;
            EQ[3] = d2t + (one + num_traits<T>::from_int(3) * mt) * d1t + mt * mt * mt;
            for (int t = 1; t <= tmax; ++t) {
                T val = num_factorial<T>(static_cast<unsigned>(t)) /
                        num_pow_int(mu[i], static_cast<unsigned>(t));
                for (int tau = 0; tau <= t; ++tau)
                    val += a(static_cast<std::size_t>(t - 1), static_cast<std::size_t>(tau)) *
                           EQ[tau];
                // correction over the states in which a server is idle
                for (int j = 0; j <= bv[i] - 1; ++j) {
                    T inner = zero;
                    for (int tau = 0; tau <= t; ++tau) {
                        // j^tau with the convention 0^0 = 1 of the reference
                        const T jp = tau == 0 ? one
                                              : num_pow_int(num_traits<T>::from_int(j),
                                                            static_cast<unsigned>(tau));
                        inner += a(static_cast<std::size_t>(t - 1), static_cast<std::size_t>(tau)) * jp;
                    }
                    val -= Prow[rl](i, static_cast<std::size_t>(j)) * inner;
                }
                res.WM[static_cast<std::size_t>(t - 1)](i, l) = val;
            }
            res.W(i, l) = res.WM[0](i, l);
        }
    }

    if (tmax >= 2)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t l = 0; l < R; ++l)
                res.WVar(i, l) = res.WM[1](i, l) - res.WM[0](i, l) * res.WM[0](i, l);
    if (tmax >= 3)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t l = 0; l < R; ++l) {
                const T mu3 = res.WM[2](i, l) -
                              num_traits<T>::from_int(3) * res.WM[0](i, l) * res.WM[1](i, l) +
                              num_traits<T>::from_int(2) * res.WM[0](i, l) * res.WM[0](i, l) *
                                  res.WM[0](i, l);
                const double var = num_traits<T>::to_double(res.WVar(i, l));
                res.WSkew(i, l) = var > 0.0 ? num_traits<T>::to_double(mu3) / std::pow(var, 1.5)
                                            : std::nan("");
            }

    return res;
}

/** pfqn_sens_respt with single servers and moments up to order three. */
template <class T>
SensResptResult<T> pfqn_sens_respt(const std::vector<T>& S, const Matrix<T>& V,
                                   const std::vector<int>& N, const std::vector<T>& Z) {
    return pfqn_sens_respt(S, V, N, Z, std::vector<int>(), 3);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_SENS_RESPT_H
