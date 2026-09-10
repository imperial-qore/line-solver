/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_ME_ME_CQN_H
#define LINE_API_ME_ME_CQN_H

/**
 * Maximum-entropy algorithm for closed multiclass queueing networks.
 *
 * Templated port of matlab/src/api/me/me_cqn.m, cross-checked against
 * jar/src/main/java/jline/api/nc/Me_cqn.java. Implements the two-stage
 * algorithm of Kouvatsos (1994) Section 3.3 for networks of G/G/1 and
 * G/G/inf queues:
 *
 *   Stage 1 solves a PSEUDO-OPEN network at trial class throughputs X, using
 *     the GE-type fixed point of Section 3.2 on the class-composed streams,
 *     and moves X until sum_i L(i,r) = N(r). The update is damped and
 *     step-clamped, and X is capped below the saturation point of every
 *     single-server station, so the iteration cannot walk into an unstable
 *     pseudo-open network.
 *   Stage 2 builds the ME product form (3.8) from the Stage 1 Lagrangian
 *     coefficients, normalizes it by a multiclass convolution over the
 *     population lattice, and iterates the work-rate (flow) equations until
 *     the throughputs implied by the closed solution agree with those used to
 *     parametrize the building blocks.
 *
 * The coefficient functions f_i are evaluated in the log domain and rescaled
 * by their maximum before the convolution, because they are products of up to
 * sum(N) factors; the per-station scaling cancels in the marginals.
 *
 * ARITHMETIC: log, exp and a damped tolerance-stopped fixed point.
 *   static_assert(num_traits<T>::has_transcendental)
 * This is the algorithm in the port that most repays extra precision: the
 * convolution of the rescaled coefficients cancels heavily at high
 * population, and the population constraint sum_i L(i,r) = N(r) is the
 * observable that degrades first when it does.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/api/me/me_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace me {

namespace detail {

/** Stage 1 output: the pseudo-open decomposition at the current flows. */
template <class T>
struct PseudoOpen {
    Matrix<T> L, Cd, rho;
};

/**
 * Scales the class throughputs uniformly so that every single-server station
 * of the pseudo-open network stays below utilization 0.999.
 */
template <class T>
void me_cqn_capacity_cap(std::vector<T>& X, const Matrix<T>& V, const Matrix<T>& mu,
                         const std::vector<long>& c, std::size_t M, std::size_t R) {
    const T zero = num_traits<T>::from_int(0);
    T maxrho = zero;
    for (std::size_t i = 0; i < M; ++i) {
        if (is_is(c, i)) continue;
        T rho_i = zero;
        for (std::size_t r = 0; r < R; ++r)
            if (V(i, r) > zero && mu(i, r) > zero) rho_i += X[r] * V(i, r) / mu(i, r);
        if (rho_i > maxrho) maxrho = rho_i;
    }
    const T cap = num_traits<T>::from_rational(999, 1000);
    if (maxrho >= cap) {
        const T f = cap / maxrho;
        for (std::size_t r = 0; r < R; ++r) X[r] *= f;
    }
}

/**
 * GE-type fixed point of the open algorithm on the pseudo-open network: no
 * external arrivals, flows given by lambda. The scvs are computed on the
 * class-composed streams and disaggregated by thinning, which is what keeps
 * the closed solution consistent with the single-class one when the classes
 * are statistically identical.
 */
template <class T>
PseudoOpen<T> me_cqn_pseudoopen(std::size_t M, std::size_t R, const Matrix<T>& lambda,
                                const Matrix<T>& mu, const Matrix<T>& mueff,
                                const Matrix<T>& Cseff, const std::vector<Matrix<T>>& Peff,
                                const Matrix<T>& selfp, const std::vector<long>& c,
                                const std::vector<char>& insens, Matrix<T>& Ca,
                                const MeOptions& opt) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T tol = num_traits<T>::from_double(opt.tol);

    Matrix<T> lameff(M, R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) lameff(i, r) = lambda(i, r) * (one - selfp(i, r));

    PseudoOpen<T> po;
    po.rho = Matrix<T>(M, R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) {
            if (!(mu(i, r) > zero)) continue;
            po.rho(i, r) = is_is(c, i) ? lameff(i, r) / mueff(i, r) : lambda(i, r) / mu(i, r);
        }

    // Class composition per station
    std::vector<T> lam_a(M, zero), mu_a(M, zero), Cs_a(M, one);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) lam_a[i] += lameff(i, r);
    for (std::size_t i = 0; i < M; ++i) {
        if (!(lam_a[i] > zero)) continue;
        T ES = zero, ES2 = zero;
        for (std::size_t u = 0; u < R; ++u) {
            if (!(lameff(i, u) > zero && mueff(i, u) > zero)) continue;
            const T wu = lameff(i, u) / lam_a[i];
            ES += wu / mueff(i, u);
            ES2 += wu * (Cseff(i, u) + one) / (mueff(i, u) * mueff(i, u));
        }
        if (ES > zero) {
            mu_a[i] = one / ES;
            Cs_a[i] = ES2 / (ES * ES) - one;
        }
    }
    Matrix<T> Pa(M, M, zero);
    for (std::size_t j = 0; j < M; ++j) {
        if (!(lam_a[j] > zero)) continue;
        for (std::size_t i = 0; i < M; ++i) {
            T num = zero;
            for (std::size_t r = 0; r < R; ++r)
                if (lameff(j, r) > zero) num += lameff(j, r) * Peff[r](j, i);
            Pa(j, i) = num / lam_a[j];
        }
    }

    // Fixed point on the aggregate arrival scvs, warm started from Ca
    std::vector<T> Ca_a(M, one), Cd_a(M, one), L_a(M, zero);
    for (std::size_t i = 0; i < M; ++i) {
        if (!(lam_a[i] > zero)) continue;
        for (std::size_t r = 0; r < R; ++r)
            if (lameff(i, r) > zero) {
                Ca_a[i] = one + (Ca(i, r) - one) * lam_a[i] / lameff(i, r);
                break;
            }
    }
    for (long it = 1; it <= opt.maxiter; ++it) {
        const std::vector<T> Ca_old = Ca_a;
        for (std::size_t i = 0; i < M; ++i) {
            if (!(lam_a[i] > zero)) continue;
            T rho_i = zero;
            for (std::size_t r = 0; r < R; ++r) rho_i += po.rho(i, r);
            if (is_is(c, i)) {
                L_a[i] = lam_a[i] / mu_a[i];
                Cd_a[i] = Ca_a[i];
            } else if (rho_i < one) {
                if (insens[i])
                    L_a[i] = rho_i / (one - rho_i);
                else
                    L_a[i] = rho_i * (Ca_a[i] + one) / two +
                             rho_i * rho_i * (Ca_a[i] + Cs_a[i]) / (two * (one - rho_i));
                Cd_a[i] = two * L_a[i] * (one - rho_i) + Ca_a[i] * (one - two * rho_i);
            }
        }
        for (std::size_t i = 0; i < M; ++i) {
            if (!(lam_a[i] > zero)) continue;
            T sum_inv = zero;
            for (std::size_t j = 0; j < M; ++j) {
                if (!(Pa(j, i) > zero) || !(lam_a[j] > zero)) continue;
                const T Cdji = one + Pa(j, i) * (Cd_a[j] - one);
                sum_inv += (lam_a[j] * Pa(j, i) / lam_a[i]) / (Cdji + one);
            }
            if (sum_inv > zero) Ca_a[i] = -one + one / sum_inv;
        }
        T delta = zero;
        for (std::size_t i = 0; i < M; ++i) {
            const T d = num_abs(T(Ca_a[i] - Ca_old[i]));
            if (d > delta) delta = d;
        }
        if (delta < tol) break;
    }

    // Disaggregation by thinning, then per-class mean queue lengths
    po.L = Matrix<T>(M, R, zero);
    po.Cd = Matrix<T>(M, R, one);
    for (std::size_t i = 0; i < M; ++i) {
        T rho_i = zero;
        for (std::size_t r = 0; r < R; ++r) rho_i += po.rho(i, r);
        for (std::size_t r = 0; r < R; ++r) {
            if (!(lameff(i, r) > zero)) continue;
            const T pr = lameff(i, r) / lam_a[i];
            Ca(i, r) = one + pr * (Ca_a[i] - one);
            po.Cd(i, r) = one + pr * (Cd_a[i] - one);
        }
        if (is_is(c, i)) {
            for (std::size_t r = 0; r < R; ++r)
                if (lameff(i, r) > zero && mueff(i, r) > zero)
                    po.L(i, r) = lameff(i, r) / mueff(i, r);
        } else if (rho_i < one) {
            if (insens[i]) {
                for (std::size_t r = 0; r < R; ++r)
                    if (lameff(i, r) > zero && mueff(i, r) > zero)
                        po.L(i, r) = po.rho(i, r) / (one - rho_i);
            } else {
                T resid = zero;
                for (std::size_t u = 0; u < R; ++u)
                    if (lameff(i, u) > zero && mueff(i, u) > zero)
                        resid += lameff(i, u) * (Cseff(i, u) + Ca(i, u)) /
                                 (mueff(i, u) * mueff(i, u));
                for (std::size_t r = 0; r < R; ++r)
                    if (lameff(i, r) > zero && mueff(i, r) > zero)
                        po.L(i, r) = po.rho(i, r) * (Ca(i, r) + one) / two +
                                     lameff(i, r) * resid / (two * (one - rho_i));
            }
        }
    }
    return po;
}

/** Mixed-radix enumeration of the population lattice {0..N(1)} x ... */
inline std::vector<std::vector<long>> me_cqn_lattice(const std::vector<long>& N) {
    std::size_t PIdx = 1;
    for (std::size_t r = 0; r < N.size(); ++r) PIdx *= static_cast<std::size_t>(N[r] + 1);
    std::vector<std::vector<long>> Dec(PIdx, std::vector<long>(N.size(), 0));
    for (std::size_t p = 0; p < PIdx; ++p) {
        std::size_t q = p;
        for (std::size_t r = 0; r < N.size(); ++r) {
            const std::size_t sz = static_cast<std::size_t>(N[r] + 1);
            Dec[p][r] = static_cast<long>(q % sz);
            q /= sz;
        }
    }
    return Dec;
}

/**
 * Auxiliary functions f_i of the ME solution (3.8): the right-hand sides of
 * (3.2) and (3.4) with the (1-rho) factor removed, evaluated from the Stage 1
 * Lagrangian coefficients and rescaled by their maximum.
 */
template <class T>
Matrix<T> me_cqn_coefficients(std::size_t M, std::size_t R, const std::vector<long>& N,
                              const std::vector<std::vector<long>>& Dec, const Matrix<T>& Lpo,
                              const Matrix<T>& rho_po, const Matrix<T>& lambda,
                              const Matrix<T>& mueff, const Matrix<T>& Cseff, const Matrix<T>& Ca,
                              const std::vector<long>& c, const Matrix<T>& selfp) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t PIdx = Dec.size();
    Matrix<T> F(PIdx, M, zero);
    Matrix<T> lameff(M, R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) lameff(i, r) = lambda(i, r) * (one - selfp(i, r));

    for (std::size_t i = 0; i < M; ++i) {
        if (is_is(c, i)) {
            // GE/GE/inf: f(n) = prod_r prod_{k=1}^{n_r} g_r(k)
            std::vector<std::vector<T>> logg(R);
            std::vector<std::vector<char>> ok(R);
            for (std::size_t r = 0; r < R; ++r) {
                logg[r].assign(static_cast<std::size_t>(N[r]), zero);
                ok[r].assign(static_cast<std::size_t>(N[r]), 0);
                for (long j = 1; j <= N[r]; ++j) {
                    if (!(lameff(i, r) > zero && mueff(i, r) > zero)) continue;
                    const T den = num_traits<T>::from_int(j) * mueff(i, r) *
                                  (Ca(i, r) + Cseff(i, r));
                    if (!(den > zero)) continue;
                    const T gj = (lameff(i, r) * (one + Cseff(i, r)) +
                                  num_traits<T>::from_int(j - 1) * mueff(i, r) *
                                      (Ca(i, r) - one)) /
                                 den;
                    if (!(gj > zero)) continue;
                    logg[r][static_cast<std::size_t>(j - 1)] = num_log(gj);
                    ok[r][static_cast<std::size_t>(j - 1)] = 1;
                }
            }
            for (std::size_t p = 0; p < PIdx; ++p) {
                T val = zero;
                bool good = true;
                for (std::size_t r = 0; r < R && good; ++r)
                    for (long j = 1; j <= Dec[p][r]; ++j) {
                        if (!ok[r][static_cast<std::size_t>(j - 1)]) {
                            good = false;
                            break;
                        }
                        val += logg[r][static_cast<std::size_t>(j - 1)];
                    }
                F(p, i) = good ? num_exp(val) : zero;
            }
            F(0, i) = one;
        } else {
            // GE/GE/1: f(n) = ((|n|-1)!/prod_r n_r!)
            //   * sum_r n_r (g_r x_r) x_r^{n_r-1} prod_{s!=r} x_s^{n_s}
            T rho_i = zero, Li = zero;
            for (std::size_t r = 0; r < R; ++r) {
                rho_i += rho_po(i, r);
                Li += Lpo(i, r);
            }
            std::vector<T> x(R, zero), gx(R, zero);
            if (Li > zero && rho_i < one) {
                for (std::size_t r = 0; r < R; ++r) {
                    if (!(lambda(i, r) > zero)) continue;
                    const T d = Lpo(i, r) - rho_po(i, r);
                    x[r] = d > zero ? d / Li : zero;
                    gx[r] = rho_po(i, r) * rho_i / ((one - rho_i) * Li);
                }
            }
            for (std::size_t p = 0; p < PIdx; ++p) {
                long ntot = 0;
                bool absent = false;
                for (std::size_t r = 0; r < R; ++r) {
                    ntot += Dec[p][r];
                    if (Dec[p][r] > 0 && !(lambda(i, r) > zero)) absent = true;
                }
                if (ntot == 0) {
                    F(p, i) = one;
                    continue;
                }
                if (absent) {
                    F(p, i) = zero;  // a class that does not visit this station
                    continue;
                }
                T logmult = log_factorial<T>(ntot - 1);
                for (std::size_t r = 0; r < R; ++r) logmult -= log_factorial<T>(Dec[p][r]);
                T tot = zero;
                for (std::size_t r = 0; r < R; ++r) {
                    if (!(Dec[p][r] > 0) || !(gx[r] > zero)) continue;
                    T lterm = num_log(num_traits<T>::from_int(Dec[p][r])) + num_log(gx[r]);
                    bool good = true;
                    for (std::size_t s = 0; s < R; ++s) {
                        const long es = (s == r) ? Dec[p][s] - 1 : Dec[p][s];
                        if (es <= 0) continue;
                        if (!(x[s] > zero)) {
                            good = false;
                            break;
                        }
                        lterm += num_traits<T>::from_int(es) * num_log(x[s]);
                    }
                    if (good) tot += num_exp(T(logmult + lterm));
                }
                F(p, i) = tot;
            }
        }
        T fmax = zero;
        for (std::size_t p = 0; p < PIdx; ++p)
            if (F(p, i) > fmax) fmax = F(p, i);
        if (fmax > zero)
            for (std::size_t p = 0; p < PIdx; ++p) F(p, i) /= fmax;
        if (!(F(0, i) > zero)) F(0, i) = num_traits<T>::from_double(1e-300);
    }
    return F;
}

/** Convolution of a partial normalizing constant with one station term. */
template <class T>
std::vector<T> me_cqn_convpair(const std::vector<T>& G, const std::vector<T>& f,
                               const std::vector<std::vector<long>>& Dec,
                               const std::vector<long>& N, const std::vector<std::size_t>& rad) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t PIdx = Dec.size();
    std::vector<T> G2(PIdx, zero);
    for (std::size_t p = 0; p < PIdx; ++p) {
        if (f[p] == zero) continue;
        for (std::size_t q = 0; q < PIdx; ++q) {
            if (G[q] == zero) continue;
            std::size_t idx = 0;
            bool fits = true;
            for (std::size_t r = 0; r < N.size(); ++r) {
                const long t = Dec[p][r] + Dec[q][r];
                if (t > N[r]) {
                    fits = false;
                    break;
                }
                idx += static_cast<std::size_t>(t) * rad[r];
            }
            if (fits) G2[idx] += f[p] * G[q];
        }
    }
    return G2;
}

/** Return value of the convolution step: marginals and busy probabilities. */
template <class T>
struct MeCqnConv {
    Matrix<T> L;
    std::vector<T> U;
};

/**
 * Normalizing constant by convolving the f_i over the population lattice, and
 * the per-station marginals by prefix/suffix convolutions.
 */
template <class T>
MeCqnConv<T> me_cqn_convolve(std::size_t M, std::size_t R, const std::vector<long>& N,
                             const std::vector<std::vector<long>>& Dec, const Matrix<T>& F) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t PIdx = Dec.size();
    std::vector<std::size_t> rad(R, 1);
    for (std::size_t r = 1; r < R; ++r)
        rad[r] = rad[r - 1] * static_cast<std::size_t>(N[r - 1] + 1);

    std::vector<std::vector<T>> Fcol(M, std::vector<T>(PIdx, zero));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t p = 0; p < PIdx; ++p) Fcol[i][p] = F(p, i);

    std::vector<T> G0(PIdx, zero);
    G0[0] = one;
    std::vector<std::vector<T>> Gpre(M + 1), Gsuf(M + 1);
    Gpre[0] = G0;
    for (std::size_t k = 0; k < M; ++k)
        Gpre[k + 1] = me_cqn_convpair(Gpre[k], Fcol[k], Dec, N, rad);
    Gsuf[M] = G0;
    for (std::size_t k = M; k-- > 0;)
        Gsuf[k] = me_cqn_convpair(Gsuf[k + 1], Fcol[k], Dec, N, rad);

    const T Z = Gpre[M][PIdx - 1];
    if (!(Z > zero)) throw NumericError("me_cqn: the normalizing constant vanished");

    MeCqnConv<T> out;
    out.L = Matrix<T>(M, R, zero);
    out.U.assign(M, zero);
    for (std::size_t i = 0; i < M; ++i) {
        const std::vector<T> Grest = me_cqn_convpair(Gpre[i], Gsuf[i + 1], Dec, N, rad);
        for (std::size_t p = 0; p < PIdx; ++p) {
            if (!(F(p, i) > zero)) continue;
            std::size_t q = 0;
            for (std::size_t r = 0; r < R; ++r)
                q += static_cast<std::size_t>(N[r] - Dec[p][r]) * rad[r];
            const T pin = F(p, i) * Grest[q] / Z;
            if (p > 0) out.U[i] += pin;
            for (std::size_t r = 0; r < R; ++r)
                if (Dec[p][r] > 0) out.L(i, r) += num_traits<T>::from_int(Dec[p][r]) * pin;
        }
    }
    return out;
}

}  // namespace detail

/**
 * @param M       number of stations
 * @param R       number of classes
 * @param N       class populations (R)
 * @param mu      service rates (M x R)
 * @param Cs      service scvs (M x R)
 * @param P       routing, R matrices (M x M)
 * @param c       servers per station, 0 for an infinite-server station;
 *                finite values must be 1
 * @param refstat_in reference station per class, or -1 for the first station the
 *                class is served at
 * @param insens  insensitive discipline flags per station
 * @param opt     tolerance and iteration budget
 */
template <class T>
MeResult<T> me_cqn(std::size_t M, std::size_t R, const std::vector<long>& N, const Matrix<T>& mu,
                   const Matrix<T>& Cs, const std::vector<Matrix<T>>& P,
                   const std::vector<long>& c, const std::vector<long>& refstat_in,
                   const std::vector<char>& insens, const MeOptions& opt = MeOptions()) {
    static_assert(num_traits<T>::has_transcendental, "me_cqn requires transcendental arithmetic");
    detail::check_dims(M, R, mu, Cs, P, c, insens, "me_cqn");
    if (N.size() != R) throw InputError("me_cqn: one population per class");
    if (refstat_in.size() != R) throw InputError("me_cqn: one reference station per class");
    for (std::size_t i = 0; i < M; ++i)
        if (c[i] > 1) throw InputError("me_cqn: only single-server and IS stations are supported");
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] < 0) throw InputError("me_cqn: populations must be nonnegative");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T half = num_traits<T>::from_rational(1, 2);
    const T tol = num_traits<T>::from_double(opt.tol);

    // Feedback correction
    std::vector<Matrix<T>> Peff = P;
    Matrix<T> mueff = mu, Cseff = Cs, selfp(M, R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) {
            const T pii = P[r](i, i);
            if (!(pii > zero)) continue;
            if (pii >= one) throw InputError("me_cqn: a self-loop probability of one");
            selfp(i, r) = pii;
            mueff(i, r) = mu(i, r) * (one - pii);
            Cseff(i, r) = pii + (one - pii) * Cs(i, r);
            for (std::size_t j = 0; j < M; ++j) Peff[r](i, j) = P[r](i, j) / (one - pii);
            Peff[r](i, i) = zero;
        }

    // Visit ratios, normalized at the reference station of each class
    Matrix<T> V(M, R, zero);
    std::vector<long> refstat = refstat_in;
    for (std::size_t r = 0; r < R; ++r) {
        std::size_t ref;
        if (refstat[r] < 0) {
            std::size_t k = 0;
            while (k < M && !(mu(k, r) > zero)) ++k;
            if (k == M) throw InputError("me_cqn: a class is not served anywhere");
            ref = k;
        } else {
            ref = static_cast<std::size_t>(refstat[r]);
            if (ref >= M) throw InputError("me_cqn: the reference station is out of range");
        }
        refstat[r] = static_cast<long>(ref);
        Matrix<T> A(M, M, zero);
        std::vector<T> b(M, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < M; ++j) A(i, j) = (i == j ? one : zero) - P[r](j, i);
        for (std::size_t j = 0; j < M; ++j) A(ref, j) = zero;
        A(ref, ref) = one;
        b[ref] = one;
        const std::vector<T> v = detail::linear_solve(A, b);
        const T eps = num_traits<T>::from_double(1e-14);
        for (std::size_t i = 0; i < M; ++i) V(i, r) = num_abs(v[i]) < eps ? zero : v[i];
    }

    // Stage 1: initial throughputs at half the single-station capacity bound
    std::vector<T> X(R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        bool have = false;
        T capr = zero;
        for (std::size_t i = 0; i < M; ++i) {
            if (detail::is_is(c, i) || !(V(i, r) > zero) || !(mu(i, r) > zero)) continue;
            const T cand = mu(i, r) / V(i, r);
            if (!have || cand < capr) {
                capr = cand;
                have = true;
            }
        }
        if (!have) capr = one;  // IS-only class
        X[r] = half * capr / num_traits<T>::from_int(static_cast<long>(R));
    }

    MeResult<T> out;
    out.Ca = Matrix<T>(M, R, one);
    Matrix<T> lambda(M, R, zero);
    detail::PseudoOpen<T> po;
    po.L = Matrix<T>(M, R, zero);
    po.Cd = Matrix<T>(M, R, one);
    po.rho = Matrix<T>(M, R, zero);

    const long maxit1 = std::min<long>(opt.maxiter, 100);
    for (long it1 = 1; it1 <= maxit1; ++it1) {
        ++out.iter;
        detail::me_cqn_capacity_cap(X, V, mu, c, M, R);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) lambda(i, r) = V(i, r) * X[r];
        po = detail::me_cqn_pseudoopen(M, R, lambda, mu, mueff, Cseff, Peff, selfp, c, insens,
                                       out.Ca, opt);
        T err1 = zero;
        std::vector<T> Ltot(R, zero);
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t i = 0; i < M; ++i) Ltot[r] += po.L(i, r);
        for (std::size_t r = 0; r < R; ++r) {
            if (!(N[r] > 0) || !(Ltot[r] > zero)) continue;
            const T e = num_abs(T(Ltot[r] - num_traits<T>::from_int(N[r]))) /
                        num_traits<T>::from_int(N[r]);
            if (e > err1) err1 = e;
        }
        if (err1 < tol) break;
        const std::vector<T> Xold = X;
        for (std::size_t r = 0; r < R; ++r) {
            if (!(Ltot[r] > zero)) continue;
            T fac = detail::num_sqrt(T(num_traits<T>::from_int(N[r]) / Ltot[r]));
            const T lo = num_traits<T>::from_rational(1, 4), hi = num_traits<T>::from_int(4);
            if (fac < lo) fac = lo;
            if (fac > hi) fac = hi;
            X[r] = half * X[r] + half * X[r] * fac;
        }
        // Stall guard: the stability cap can bind before the population
        // target is met, and then X stops moving
        std::vector<T> Xc = X;
        detail::me_cqn_capacity_cap(Xc, V, mu, c, M, R);
        T move = zero;
        for (std::size_t r = 0; r < R; ++r) {
            const T den = Xold[r] > num_traits<T>::from_double(1e-12)
                              ? Xold[r]
                              : num_traits<T>::from_double(1e-12);
            const T d = num_abs(T(Xc[r] - Xold[r])) / den;
            if (d > move) move = d;
        }
        if (move < tol) break;
    }

    // Stage 2: closed ME solution by convolution, iterated on the flows
    const std::vector<std::vector<long>> Dec = detail::me_cqn_lattice(N);
    out.L = po.L;
    out.rho = po.rho;
    T err2 = zero;
    for (long it2 = 1; it2 <= opt.maxiter; ++it2) {
        ++out.iter;
        const Matrix<T> F = detail::me_cqn_coefficients(M, R, N, Dec, po.L, po.rho, lambda, mueff,
                                                        Cseff, out.Ca, c, selfp);
        const detail::MeCqnConv<T> conv = detail::me_cqn_convolve(M, R, N, Dec, F);
        out.L = conv.L;
        out.rho = Matrix<T>(M, R, zero);
        std::vector<T> Xhat(R, zero);
        for (std::size_t r = 0; r < R; ++r) {
            T num = zero, den = zero;
            for (std::size_t i = 0; i < M; ++i) {
                if (!(lambda(i, r) > zero)) continue;
                if (detail::is_is(c, i)) {
                    out.rho(i, r) = out.L(i, r);
                    num += out.L(i, r) * mueff(i, r) / (one - selfp(i, r));
                } else {
                    T rho_i = zero;
                    for (std::size_t u = 0; u < R; ++u) rho_i += po.rho(i, u);
                    if (rho_i > zero) out.rho(i, r) = conv.U[i] * po.rho(i, r) / rho_i;
                    num += out.rho(i, r) * mu(i, r);
                }
                den += V(i, r);
            }
            if (den > zero) Xhat[r] = num / den;
        }
        err2 = zero;
        for (std::size_t r = 0; r < R; ++r) {
            if (!(X[r] > zero)) continue;
            const T e = num_abs(T(Xhat[r] - X[r])) / X[r];
            if (e > err2) err2 = e;
        }
        if (err2 < tol) {
            out.converged = true;
            break;
        }
        const std::vector<T> Xold = X;
        for (std::size_t r = 0; r < R; ++r) X[r] = half * X[r] + half * Xhat[r];
        detail::me_cqn_capacity_cap(X, V, mu, c, M, R);
        T move = zero;
        for (std::size_t r = 0; r < R; ++r) {
            const T den = Xold[r] > num_traits<T>::from_double(1e-12)
                              ? Xold[r]
                              : num_traits<T>::from_double(1e-12);
            const T d = num_abs(T(X[r] - Xold[r])) / den;
            if (d > move) move = d;
        }
        if (move < tol) break;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) lambda(i, r) = V(i, r) * X[r];
        po = detail::me_cqn_pseudoopen(M, R, lambda, mu, mueff, Cseff, Peff, selfp, c, insens,
                                       out.Ca, opt);
    }

    // Response times by Little's law on the visit-inclusive throughputs
    out.lambda = Matrix<T>(M, R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) out.lambda(i, r) = V(i, r) * X[r];
    out.W = Matrix<T>(M, R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r)
            if (out.lambda(i, r) > zero) out.W(i, r) = out.L(i, r) / out.lambda(i, r);
    out.Cd = po.Cd;
    out.X = X;
    return out;
}

}  // namespace me
}  // namespace line

#endif  // LINE_API_ME_ME_CQN_H
