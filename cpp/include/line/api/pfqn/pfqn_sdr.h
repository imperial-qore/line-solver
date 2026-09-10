/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_SDR_H
#define LINE_API_PFQN_SDR_H

/**
 * Product-form state-dependent routing.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_sdr*.m, from A. E. Krzesinski,
 * "Multiclass Queueing Networks with State-Dependent Routing", Performance
 * Evaluation 7(2):125-143, 1987, the multiclass generalization of D. Towsley,
 * "Queuing Network Models with State-Dependent Routing", J. ACM 27(2):323-337,
 * 1980.
 *
 * A network is split into a subnetwork Q(V,V) under SDR and its complement
 * M-V. Q(V,V) has one entry centre e and one departure centre d, both outside
 * it, and is partitioned into disjoint branches arranged in a hierarchy of
 * nested subnetworks V_1 > ... > V_T. Each branch has one entry centre, one
 * departure centre, and may hold several centres between them.
 *
 * Branch index 1 denotes the complement M-V and is unused, in the paper and in
 * every codebase. The SDR branches are numbered 2..B, that is indices 1..B-1 of
 * the zero-based arrays here. Keeping the paper's numbering is what lets d(t,b)
 * be transcribed straight from the text.
 *
 * UNLIKE the MATLAB, JAR and Python twins this routine carries no logarithms.
 * The unnormalized weight of eq. (16) is a product of field operations on the
 * inputs, so the whole evaluation stays inside T and an exact rational backend
 * returns an exact normalizing constant. The price is that a double backend can
 * overflow where the log-domain twins would not; that bites only at populations
 * far beyond what state enumeration can reach anyway.
 */

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mc/dtmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/**
 * Topology and coefficients of a state-dependent routing subnetwork.
 *
 * Centre indices are whatever space the caller uses consistently: node indices
 * for the per-state routing probabilities, station indices for the product
 * form.
 */
struct SdrStruct {
    /** Entry centre e of Q(V,V). */
    std::size_t entry = 0;
    /** Departure centre d of Q(V,V); may equal `entry`. */
    std::size_t departure = 0;
    /** branch[b] holds the centres of branch b, b >= 1; branch[0] is unused. */
    std::vector<std::vector<std::size_t>> branch;
    /** entryOf[b] is the entry centre e(b) of branch b. */
    std::vector<std::size_t> entryOf;
    /** departureOf[b] is the departure centre d(b) of branch b. */
    std::vector<std::size_t> departureOf;
    /** level[b] is the unique t with B_b in V_t - V_{t+1}; level[0] is unused. */
    std::vector<std::size_t> level;
    /** Coefficients C_t of eq. (11), length T. */
    std::vector<double> C;
    /** Coefficients d_tb of eq. (11), T rows by B columns. */
    Matrix<double> d;

    bool empty() const { return branch.size() < 2; }
};

/** Derived coefficients of an SDR structure, eqs. (11)-(14). */
struct SdrCoeff {
    std::size_t T = 0;
    std::size_t B = 0;
    /** inA[t] holds the branch indices b with level[b] >= t, the set A_t. */
    std::vector<std::vector<std::size_t>> inA;
    /** D_tt = sum over A_t of d_tb; index 0 unused. */
    std::vector<double> Dtt;
    /** D_{t-1,t} = sum over A_t of d_{t-1,b}; indices 0 and 1 unused. */
    std::vector<double> Dprev;
    /** Largest branch population with delta nonnegative; -1 marks unbounded. */
    std::vector<double> mmax;
    /** Largest subnetwork population with omega nonnegative; -1 marks unbounded. */
    std::vector<double> vmax;
    const SdrStruct* sdr = nullptr;
};

/**
 * Validates an SDR structure and returns its derived coefficients.
 *
 * The population bounds are consequences of the coefficients, not independent
 * inputs: with C_t negative the routing enforces m_b <= d_tb/(-C_t) and
 * v_t <= D_tt/(-C_t) by itself, because the cumulative Delta hits a zero factor
 * exactly at the bound.
 */
inline SdrCoeff pfqn_sdrcoeff(const SdrStruct& sdr) {
    const std::size_t B = sdr.branch.size();
    const std::size_t T = sdr.C.size();
    if (B < 2)
        throw InputError("pfqn_sdrcoeff: an SDR structure must declare at least one branch "
                         "(branch indices start at 2)");
    if (T < 1)
        throw InputError("pfqn_sdrcoeff: an SDR structure must declare at least one level of "
                         "subnetwork nesting");
    if (sdr.level.size() != B)
        throw InputError("pfqn_sdrcoeff: level must have one entry per branch index");
    if (sdr.entryOf.size() != B || sdr.departureOf.size() != B)
        throw InputError("pfqn_sdrcoeff: entryOf and departureOf must have one entry per branch index");
    if (sdr.d.rows() < T || sdr.d.cols() < B)
        throw InputError("pfqn_sdrcoeff: the coefficient matrix d is too small for the declared "
                         "levels and branches");
    for (std::size_t b = 1; b < B; ++b)
        if (sdr.level[b] < 1 || sdr.level[b] > T)
            throw InputError("pfqn_sdrcoeff: SDR branch levels must lie in 1..T");

    // The nesting V_1 > ... > V_T must be strict, else two hierarchically
    // adjacent subnetworks coincide and the ratio of eq. (10) is not the paper's.
    for (std::size_t t = 1; t <= T; ++t) {
        bool found = false;
        for (std::size_t b = 1; b < B && !found; ++b) found = (sdr.level[b] == t);
        if (!found)
            throw InputError("pfqn_sdrcoeff: SDR level " + std::to_string(t) +
                             " carries no branch: the subnetwork nesting must be strict");
    }

    std::vector<std::size_t> seen;
    for (std::size_t b = 1; b < B; ++b) {
        if (sdr.branch[b].empty())
            throw InputError("pfqn_sdrcoeff: SDR branch " + std::to_string(b + 1) + " is empty");
        bool hasEntry = false, hasDeparture = false;
        for (std::size_t k = 0; k < sdr.branch[b].size(); ++k) {
            const std::size_t c = sdr.branch[b][k];
            for (std::size_t q = 0; q < seen.size(); ++q)
                if (seen[q] == c)
                    throw InputError("pfqn_sdrcoeff: SDR branches must be mutually disjoint");
            if (c == sdr.entry || c == sdr.departure)
                throw InputError("pfqn_sdrcoeff: the entry and departure centres of Q(V,V) must "
                                 "not belong to any branch");
            if (c == sdr.entryOf[b]) hasEntry = true;
            if (c == sdr.departureOf[b]) hasDeparture = true;
            seen.push_back(c);
        }
        if (!hasEntry || !hasDeparture)
            throw InputError("pfqn_sdrcoeff: the entry and departure centres of SDR branch " +
                             std::to_string(b + 1) + " must belong to that branch");
    }

    SdrCoeff c;
    c.T = T;
    c.B = B;
    c.sdr = &sdr;
    c.inA.assign(T + 1, std::vector<std::size_t>());
    c.Dtt.assign(T + 1, 0.0);
    c.Dprev.assign(T + 1, 0.0);
    for (std::size_t t = 1; t <= T; ++t) {
        for (std::size_t b = 1; b < B; ++b)
            if (sdr.level[b] >= t) {
                c.inA[t].push_back(b);
                c.Dtt[t] += sdr.d(t - 1, b);
                if (t > 1) c.Dprev[t] += sdr.d(t - 2, b);
            }
    }
    c.mmax.assign(B, -1.0);
    for (std::size_t b = 1; b < B; ++b) {
        const std::size_t t = sdr.level[b];
        if (sdr.C[t - 1] < 0) c.mmax[b] = std::floor(sdr.d(t - 1, b) / (-sdr.C[t - 1]));
    }
    c.vmax.assign(T + 1, -1.0);
    for (std::size_t t = 1; t <= T; ++t) {
        if (sdr.C[t - 1] < 0) c.vmax[t] = std::floor(c.Dtt[t] / (-sdr.C[t - 1]));
        if (t > 1 && sdr.C[t - 2] < 0) {
            const double alt = std::floor(c.Dprev[t] / (-sdr.C[t - 2]));
            c.vmax[t] = (c.vmax[t] < 0) ? alt : std::min(c.vmax[t], alt);
        }
    }
    return c;
}

/**
 * SDR routing probabilities of eq. (10).
 *
 * Entry b of the result is the probability of proceeding from the entry centre
 * e of Q(V,V) to the entry centre of branch b; entry 0 is zero because branch
 * index 1 denotes the complement M-V. The residual mass 1 - sum is the
 * probability of proceeding directly to the departure centre d, that is of
 * being denied entry into Q(V,V) and returned to e, which the paper calls the
 * busy form of waiting (Sec. 2.5).
 *
 * The probabilities are chain independent: they read the total branch and
 * subnetwork populations, not the per-chain ones. The chain-dependent form of
 * eq. (1) has no published product form and is not implemented. A branch
 * population beyond the bound SDR enforces itself is unreachable, and the
 * probability returned there is zero.
 *
 * @param c derived coefficients from pfqn_sdrcoeff
 * @param n per-centre total populations, indexed as the structure is
 */
inline std::vector<double> pfqn_sdrprob(const SdrCoeff& c, const std::vector<double>& n) {
    const SdrStruct& sdr = *c.sdr;
    std::vector<double> m(c.B, 0.0);
    for (std::size_t b = 1; b < c.B; ++b)
        for (std::size_t k = 0; k < sdr.branch[b].size(); ++k) m[b] += n[sdr.branch[b][k]];
    std::vector<double> v(c.T + 1, 0.0);
    for (std::size_t t = 1; t <= c.T; ++t)
        for (std::size_t k = 0; k < c.inA[t].size(); ++k) v[t] += m[c.inA[t][k]];

    std::vector<double> om(c.T + 1, 0.0), omprev(c.T + 1, 1.0);
    for (std::size_t s = 1; s <= c.T; ++s) {
        om[s] = sdr.C[s - 1] * v[s] + c.Dtt[s];
        if (s > 1) omprev[s] = sdr.C[s - 2] * v[s] + c.Dprev[s];
    }

    std::vector<double> P(c.B, 0.0);
    for (std::size_t b = 1; b < c.B; ++b) {
        const std::size_t t = sdr.level[b];
        bool closed = false;
        for (std::size_t s = 1; s <= t && !closed; ++s)
            closed = (om[s] <= 0.0);  // eq. (10): the branch is closed to new arrivals
        if (closed) continue;
        const double delta = sdr.C[t - 1] * m[b] + sdr.d(t - 1, b);
        if (delta <= 0.0) continue;
        double ratio = 1.0;
        for (std::size_t s = 1; s <= t; ++s) ratio *= omprev[s] / om[s];
        P[b] = delta * ratio;
    }
    return P;
}

/** Probability of being denied entry and routed straight to the departure centre. */
inline double pfqn_sdrped(const std::vector<double>& P) {
    double s = 0.0;
    for (std::size_t b = 0; b < P.size(); ++b) s += P[b];
    return 1.0 - s;
}

/** Mean performance measures returned by pfqn_sdr. */
template <class T>
struct SdrResult {
    Matrix<T> QN;  ///< mean queue lengths, centres by chains
    Matrix<T> XN;  ///< per-centre chain throughputs
    Matrix<T> UN;  ///< mean number in service, XN elementwise times S
    Matrix<T> RN;  ///< mean response times at the centre, QN elementwise over XN
    T G;           ///< normalizing constant
};

namespace detail {

/** All nonnegative integer m-vectors summing to n, appended to out. */
inline void sdr_compositions(std::size_t n, std::size_t m,
                             std::vector<std::vector<std::size_t>>& out) {
    if (m == 1) {
        out.push_back(std::vector<std::size_t>(1, n));
        return;
    }
    for (std::size_t k = 0; k <= n; ++k) {
        std::vector<std::vector<std::size_t>> tail;
        sdr_compositions(n - k, m - 1, tail);
        for (std::size_t t = 0; t < tail.size(); ++t) {
            std::vector<std::size_t> row;
            row.reserve(m);
            row.push_back(k);
            row.insert(row.end(), tail[t].begin(), tail[t].end());
            out.push_back(row);
        }
    }
}

}  // namespace detail

/**
 * Exact product form of eq. (16), by summation over the reachable state space.
 *
 * P(n) is G^-1 times the product over centres of f_i(n_i), the product over
 * levels of Omega_{t-1,t}(v_t)/Omega_tt(v_t), and the product over branches of
 * Delta_tb(m_b), with f_i(n_i) = [n_i!/beta_i(n_i)] times the product over
 * chains of gamma_ij^n_ij/n_ij! and gamma_ij = xi_ij/mu_ij.
 *
 * General in the branch topology: a branch may hold several interconnected
 * centres. Only the paper's Section 4 MVA and convolution algorithm, not ported
 * here, is restricted to single-centre branches.
 *
 * S and xi are required separately rather than as their product because under
 * SDR the xi are not visit ratios, so the per-centre throughputs cannot be
 * recovered from the demands alone.
 *
 * @param S     (M x J) mean service times 1/mu_ij
 * @param xi    (M x J) coefficients of Section 3.2, see pfqn_sdrvisits
 * @param N     (J) chain populations
 * @param sdr   the routing structure, in centre indices
 * @param alpha (M x sum(N)) load-dependent rate scalings, alpha(i, k-1) = alpha_i(k);
 *              empty for a fixed-rate centre. Use k for an infinite server and
 *              min(k, c) for a c-server centre
 */
template <class T>
SdrResult<T> pfqn_sdr(const Matrix<T>& S, const Matrix<T>& xi, const std::vector<std::size_t>& N,
                      const SdrStruct& sdr, const Matrix<T>& alpha = Matrix<T>()) {
    const std::size_t M = S.rows(), J = S.cols();
    if (xi.rows() != M || xi.cols() != J)
        throw InputError("pfqn_sdr: S and xi must have the same shape");
    if (N.size() != J)
        throw InputError("pfqn_sdr: the population vector must have one entry per chain");

    const SdrCoeff c = pfqn_sdrcoeff(sdr);
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    std::size_t Ntot = 0;
    for (std::size_t j = 0; j < J; ++j) Ntot += N[j];
    for (std::size_t b = 1; b < c.B; ++b)
        for (std::size_t k = 0; k < sdr.branch[b].size(); ++k)
            if (sdr.branch[b][k] >= M)
                throw InputError("pfqn_sdr: the SDR structure references a centre index beyond "
                                 "the number of centres");

    // beta_i(n) = alpha_i(n) beta_i(n-1), beta_i(0) = 1
    Matrix<T> beta(M, Ntot + 1, one);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 1; k <= Ntot; ++k) {
            const T a = (!alpha.empty() && i < alpha.rows() && (k - 1) < alpha.cols())
                            ? alpha(i, k - 1)
                            : one;
            beta(i, k) = beta(i, k - 1) * a;
        }

    Matrix<T> gamma(M, J, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < J; ++j) gamma(i, j) = xi(i, j) * S(i, j);

    // Cumulative Delta and Omega. A zero factor is the SDR population bound
    // closing the branch or the subnetwork, and it makes the weight vanish on
    // its own; no separate capacity test is needed anywhere.
    Matrix<T> Delta(c.B, Ntot + 1, one);
    for (std::size_t b = 1; b < c.B; ++b) {
        const std::size_t t = sdr.level[b];
        for (std::size_t n = 1; n <= Ntot; ++n) {
            const double f = sdr.C[t - 1] * static_cast<double>(n - 1) + sdr.d(t - 1, b);
            Delta(b, n) = (f <= 0.0) ? zero : Delta(b, n - 1) * num_traits<T>::from_double(f);
        }
    }
    Matrix<T> OmTT(c.T + 1, Ntot + 1, one), OmPrev(c.T + 1, Ntot + 1, one);
    for (std::size_t t = 1; t <= c.T; ++t)
        for (std::size_t n = 1; n <= Ntot; ++n) {
            const double f = sdr.C[t - 1] * static_cast<double>(n - 1) + c.Dtt[t];
            OmTT(t, n) = (f <= 0.0) ? zero : OmTT(t, n - 1) * num_traits<T>::from_double(f);
            if (t > 1) {
                const double g = sdr.C[t - 2] * static_cast<double>(n - 1) + c.Dprev[t];
                OmPrev(t, n) = (g <= 0.0) ? zero : OmPrev(t, n - 1) * num_traits<T>::from_double(g);
            }
        }

    std::vector<std::vector<std::size_t>> states;
    states.push_back(std::vector<std::size_t>());
    for (std::size_t j = 0; j < J; ++j) {
        std::vector<std::vector<std::size_t>> comps;
        detail::sdr_compositions(N[j], M, comps);
        std::vector<std::vector<std::size_t>> next;
        next.reserve(states.size() * comps.size());
        for (std::size_t a = 0; a < states.size(); ++a)
            for (std::size_t b = 0; b < comps.size(); ++b) {
                std::vector<std::size_t> row = states[a];
                row.insert(row.end(), comps[b].begin(), comps[b].end());
                next.push_back(row);
            }
        states.swap(next);
    }

    // Factorials as exact integers of T, so a rational backend stays exact
    std::vector<T> fact(Ntot + 1, one);
    for (std::size_t k = 1; k <= Ntot; ++k)
        fact[k] = fact[k - 1] * num_traits<T>::from_int(static_cast<long>(k));

    SdrResult<T> res;
    res.QN = Matrix<T>(M, J, zero);
    res.XN = Matrix<T>(M, J, zero);
    res.UN = Matrix<T>(M, J, zero);
    res.RN = Matrix<T>(M, J, zero);

    std::vector<T> w(states.size(), zero);
    T Gs = zero;
    for (std::size_t s = 0; s < states.size(); ++s) {
        const std::vector<std::size_t>& flat = states[s];
        std::vector<std::size_t> ni(M, 0);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < J; ++j) ni[i] += flat[j * M + i];

        T weight = one;
        bool alive = true;
        for (std::size_t i = 0; i < M && alive; ++i) {
            weight = weight * fact[ni[i]] / beta(i, ni[i]);
            for (std::size_t j = 0; j < J; ++j) {
                const std::size_t nij = flat[j * M + i];
                if (nij == 0) continue;
                if (!(gamma(i, j) > zero)) {
                    alive = false;
                    break;
                }
                for (std::size_t k = 0; k < nij; ++k) weight = weight * gamma(i, j);
                weight = weight / fact[nij];
            }
        }
        if (!alive) continue;

        std::vector<std::size_t> m(c.B, 0);
        for (std::size_t b = 1; b < c.B && alive; ++b) {
            for (std::size_t k = 0; k < sdr.branch[b].size(); ++k) m[b] += ni[sdr.branch[b][k]];
            if (Delta(b, m[b]) == zero) alive = false;
            weight = weight * Delta(b, m[b]);
        }
        if (!alive) continue;
        for (std::size_t t = 1; t <= c.T && alive; ++t) {
            std::size_t v = 0;
            for (std::size_t k = 0; k < c.inA[t].size(); ++k) v += m[c.inA[t][k]];
            if (OmTT(t, v) == zero) {
                alive = false;
                break;
            }
            weight = weight / OmTT(t, v);
            if (t > 1) weight = weight * OmPrev(t, v);
        }
        if (!alive) continue;
        w[s] = weight;
        Gs = Gs + weight;
    }
    if (Gs == zero)
        throw InputError("pfqn_sdr: the SDR network has no reachable state at the given "
                         "populations: the routing coefficients forbid every state");
    res.G = Gs;

    for (std::size_t s = 0; s < states.size(); ++s) {
        if (w[s] == zero) continue;
        const T p = w[s] / Gs;
        const std::vector<std::size_t>& flat = states[s];
        for (std::size_t i = 0; i < M; ++i) {
            std::size_t nitot = 0;
            for (std::size_t j = 0; j < J; ++j) nitot += flat[j * M + i];
            for (std::size_t j = 0; j < J; ++j) {
                const std::size_t nij = flat[j * M + i];
                if (nij > 0)
                    res.QN(i, j) =
                        res.QN(i, j) + p * num_traits<T>::from_int(static_cast<long>(nij));
                if (nitot > 0 && S(i, j) > zero) {
                    const T a = (!alpha.empty() && i < alpha.rows() && (nitot - 1) < alpha.cols())
                                    ? alpha(i, nitot - 1)
                                    : one;
                    res.XN(i, j) = res.XN(i, j) +
                                   p * a * num_traits<T>::from_int(static_cast<long>(nij)) /
                                       num_traits<T>::from_int(static_cast<long>(nitot)) / S(i, j);
                }
            }
        }
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < J; ++j) {
            res.UN(i, j) = res.XN(i, j) * S(i, j);
            if (res.XN(i, j) > zero) res.RN(i, j) = res.QN(i, j) / res.XN(i, j);
        }
    return res;
}

namespace detail {

/** Every population vector V with 0 <= V <= N, mixed-radix ordered. */
inline std::vector<std::vector<std::size_t>> sdr_lattice(const std::vector<std::size_t>& N) {
    std::size_t tot = 1;
    for (std::size_t j = 0; j < N.size(); ++j) tot *= N[j] + 1;
    std::vector<std::vector<std::size_t>> out(tot, std::vector<std::size_t>(N.size(), 0));
    for (std::size_t r = 0; r < tot; ++r) {
        std::size_t rem = r;
        for (std::size_t j = 0; j < N.size(); ++j) {
            out[r][j] = rem % (N[j] + 1);
            rem /= (N[j] + 1);
        }
    }
    return out;
}

/** Mixed-radix row of V in the lattice of sdr_lattice. */
inline std::size_t sdr_key(const std::vector<std::size_t>& V, const std::vector<std::size_t>& N) {
    std::size_t k = 0, mul = 1;
    for (std::size_t j = 0; j < N.size(); ++j) {
        k += V[j] * mul;
        mul *= N[j] + 1;
    }
    return k;
}

/** Omega_{t-1,t}(v)/Omega_tt(v), the cumulative ratio of eq. (16). */
inline double sdr_omega_cum(const SdrCoeff& c, std::size_t t, std::size_t v) {
    double num = 1.0, den = 1.0;
    for (std::size_t k = 0; k < v; ++k) {
        const double f = -static_cast<double>(k) + c.Dtt[t];
        if (f <= 0.0) return 0.0;
        den *= f;
        if (t > 1) {
            const double g = -static_cast<double>(k) + c.Dprev[t];
            if (g <= 0.0) return 0.0;
            num *= g;
        }
    }
    return num / den;
}

/**
 * omega_{t-1,t}(v)/omega_tt(v), the single-step ratio of eq. (10). Distinct
 * from sdr_omega_cum: the routing probability carries the lowercase omega, the
 * normalizing constant the uppercase one.
 */
inline double sdr_omega_step(const SdrCoeff& c, std::size_t t, std::size_t v) {
    const double den = -static_cast<double>(v) + c.Dtt[t];
    if (den <= 0.0) return 0.0;
    double num = 1.0;  // omega_{0,1} is one
    if (t > 1) {
        num = -static_cast<double>(v) + c.Dprev[t];
        if (num <= 0.0) return 0.0;
    }
    return num / den;
}

}  // namespace detail

/**
 * Section 4 mean value analysis and convolution.
 *
 * Same inputs and outputs as pfqn_sdr, which evaluates eq. (16) exactly by
 * state enumeration, so the two are directly comparable. This routine costs
 * O(J T M (V_1...V_J)^2) rather than the size of the state space, at the price
 * of two restrictions the paper itself imposes: every SDR branch must hold a
 * single centre, and every C_t must be negative. A C_t other than -1 is
 * rescaled internally, which leaves eqs. (10) and (16) unchanged because the
 * level factors telescope.
 *
 * Two formulas of Section 4 are corrected here, both verified against
 * pfqn_sdr. The initialise step of 4.2.2 divides by T_j(V-1_j,V_T) where the
 * convolution identity G(V)/G(V-1_j) = 1/T_j(V) gives T_j(V,V_T); this
 * implementation forms G = g_mva Omega_{T-1,T}/Omega_TT directly instead. And
 * 4.2.3's T_ij = xi_ij [d_1i - Q_i] T_j drops the state-dependent omega ratios
 * of eq. (10); the exact identity is T_ij = xi_ij T_j(N,M) E_{N-1_j}[P_{e,e(i)}].
 *
 * Unlike pfqn_sdr this routine divides by intermediate normalizing constants,
 * so it needs a field with division but no transcendentals beyond the final
 * logarithm of G.
 */
template <class T>
SdrResult<T> pfqn_sdrmva(const Matrix<T>& S, const Matrix<T>& xi, const std::vector<std::size_t>& N,
                         const SdrStruct& sdr, const Matrix<T>& alpha = Matrix<T>()) {
    const std::size_t M = S.rows(), J = S.cols();
    if (xi.rows() != M || xi.cols() != J)
        throw InputError("pfqn_sdrmva: S and xi must have the same shape");
    if (N.size() != J)
        throw InputError("pfqn_sdrmva: the population vector must have one entry per chain");

    const SdrCoeff c0 = pfqn_sdrcoeff(sdr);
    for (std::size_t b = 1; b < c0.B; ++b)
        if (sdr.branch[b].size() != 1)
            throw UnsupportedError("pfqn_sdrmva: every SDR branch must hold a single centre; the MVA "
                                   "and convolution of Krzesinski (1987) Section 4 is stated that way "
                                   "and its general case is in an unpublished technical report. Use "
                                   "pfqn_sdr, which evaluates eq. (16) exactly for any branch topology");
    for (std::size_t t = 0; t < c0.T; ++t)
        if (sdr.C[t] >= 0.0)
            throw UnsupportedError("pfqn_sdrmva: every C_t must be negative; Section 2.5 assumes it and "
                                   "Section 4 is written for C_t = -1");

    // Rescale each level to C_t = -1, which leaves eqs. (10) and (16) unchanged
    SdrStruct sdr1 = sdr;
    for (std::size_t t = 0; t < c0.T; ++t) {
        const double k = -sdr.C[t];
        sdr1.C[t] = -1.0;
        for (std::size_t b = 0; b < sdr1.d.cols(); ++b) sdr1.d(t, b) = sdr.d(t, b) / k;
    }
    const SdrCoeff c = pfqn_sdrcoeff(sdr1);

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::size_t Ntot = 0;
    for (std::size_t j = 0; j < J; ++j) Ntot += N[j];

    Matrix<T> gamma(M, J, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < J; ++j) gamma(i, j) = xi(i, j) * S(i, j);
    Matrix<T> alp(M, Ntot > 0 ? Ntot : 1, one);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < alp.cols(); ++k)
            if (!alpha.empty() && i < alpha.rows() && k < alpha.cols()) alp(i, k) = alpha(i, k);

    const std::vector<std::vector<std::size_t>> latt = detail::sdr_lattice(N);
    const std::size_t nl = latt.size();

    std::vector<bool> inV(M, false);
    std::vector<double> dvec(M, 0.0);
    std::vector<bool> hasd(M, false);
    std::vector<std::size_t> lvl(M, 0);
    for (std::size_t b = 1; b < c.B; ++b) {
        const std::size_t i = sdr1.branch[b][0];
        inV[i] = true;
        lvl[i] = sdr1.level[b];
        dvec[i] = sdr1.d(sdr1.level[b] - 1, b);
        hasd[i] = true;
    }
    std::vector<std::size_t> mv;
    for (std::size_t i = 0; i < M; ++i)
        if (!inV[i]) mv.push_back(i);

    // Sec. 4.2.1, also used for the complement with delta = 1
    struct SetOut {
        std::vector<Matrix<T>> Q;  // per lattice point, M x J
        std::vector<std::vector<T>> Tp;
        std::vector<T> g;
    };
    auto set_mva = [&](const std::vector<std::size_t>& cidx, bool sdrset) {
        SetOut o;
        o.Q.assign(nl, Matrix<T>(M, J, zero));
        o.Tp.assign(nl, std::vector<T>(J, zero));
        o.g.assign(nl, zero);
        const std::size_t nc = cidx.size();
        const std::size_t z = detail::sdr_key(std::vector<std::size_t>(J, 0), N);
        if (nc == 0) {
            for (std::size_t v = 0; v < nl; ++v) {
                std::size_t tot = 0;
                for (std::size_t j = 0; j < J; ++j) tot += latt[v][j];
                o.g[v] = (tot == 0) ? one : zero;
            }
            return o;
        }
        // Ps[k][n][v] is P_k(n : V)
        std::vector<std::vector<std::vector<T>>> Ps(
            nc, std::vector<std::vector<T>>(Ntot + 1, std::vector<T>(nl, zero)));
        o.g[z] = one;
        for (std::size_t k = 0; k < nc; ++k) Ps[k][0][z] = one;

        std::vector<std::size_t> ord(nl);
        for (std::size_t v = 0; v < nl; ++v) ord[v] = v;
        std::stable_sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) {
            std::size_t sa = 0, sb = 0;
            for (std::size_t j = 0; j < J; ++j) { sa += latt[a][j]; sb += latt[b][j]; }
            return sa < sb;
        });
        for (std::size_t oi = 0; oi < nl; ++oi) {
            const std::size_t v = ord[oi];
            const std::vector<std::size_t>& V = latt[v];
            std::size_t vv = 0;
            for (std::size_t j = 0; j < J; ++j) vv += V[j];
            if (vv == 0) continue;
            std::vector<std::vector<T>> A(nc, std::vector<T>(J, zero));
            std::vector<std::size_t> vm(J, 0);
            std::vector<bool> hasm(J, false);
            for (std::size_t j = 0; j < J; ++j) {
                if (V[j] == 0) continue;
                std::vector<std::size_t> Vm = V;
                Vm[j] -= 1;
                vm[j] = detail::sdr_key(Vm, N);
                hasm[j] = true;
                for (std::size_t k = 0; k < nc; ++k) {
                    T acc = zero;
                    for (std::size_t n = 1; n <= vv; ++n) {
                        const double df = sdrset ? (dvec[cidx[k]] - static_cast<double>(n - 1)) : 1.0;
                        if (df <= 0.0) break;
                        acc = acc + num_traits<T>::from_int(static_cast<long>(n)) *
                                        num_traits<T>::from_double(df) / alp(cidx[k], n - 1) *
                                        Ps[k][n - 1][vm[j]];
                    }
                    A[k][j] = acc;
                }
            }
            for (std::size_t j = 0; j < J; ++j) {
                if (!hasm[j]) continue;
                T den = zero;
                for (std::size_t k = 0; k < nc; ++k) den = den + gamma(cidx[k], j) * A[k][j];
                if (den > zero)
                    o.Tp[v][j] = num_traits<T>::from_int(static_cast<long>(V[j])) / den;
            }
            for (std::size_t k = 0; k < nc; ++k)
                for (std::size_t j = 0; j < J; ++j) {
                    if (!hasm[j]) continue;
                    o.Q[v](cidx[k], j) = gamma(cidx[k], j) * o.Tp[v][j] * A[k][j];
                }
            for (std::size_t k = 0; k < nc; ++k) {
                T tot = zero;
                for (std::size_t n = 1; n <= vv; ++n) {
                    const double df = sdrset ? (dvec[cidx[k]] - static_cast<double>(n - 1)) : 1.0;
                    if (df <= 0.0) break;
                    T acc = zero;
                    for (std::size_t j = 0; j < J; ++j) {
                        if (!hasm[j]) continue;
                        acc = acc + gamma(cidx[k], j) * o.Tp[v][j] * Ps[k][n - 1][vm[j]];
                    }
                    Ps[k][n][v] = num_traits<T>::from_double(df) / alp(cidx[k], n - 1) * acc;
                    tot = tot + Ps[k][n][v];
                }
                Ps[k][0][v] = one - tot;
            }
            for (std::size_t j = 0; j < J; ++j)
                if (hasm[j] && o.Tp[v][j] > zero) {
                    o.g[v] = o.g[vm[j]] / o.Tp[v][j];
                    break;
                }
        }
        return o;
    };

    const SetOut comp = set_mva(mv, false);

    std::vector<T> Gin(nl, zero);
    Gin[detail::sdr_key(std::vector<std::size_t>(J, 0), N)] = one;
    std::vector<Matrix<T>> Qin(nl, Matrix<T>(M, J, zero));
    std::vector<std::vector<T>> Ain(nl, std::vector<T>(M, zero));
    for (std::size_t tt = c.T; tt >= 1; --tt) {
        std::vector<std::size_t> St, inner;
        for (std::size_t i = 0; i < M; ++i) {
            if (inV[i] && lvl[i] == tt) St.push_back(i);
            else if (inV[i] && lvl[i] > tt) inner.push_back(i);
        }
        const SetOut lev = set_mva(St, true);
        std::vector<T> Gnew(nl, zero);
        std::vector<Matrix<T>> Qnew(nl, Matrix<T>(M, J, zero));
        std::vector<std::vector<T>> Tnew(nl, std::vector<T>(J, zero));
        std::vector<std::vector<T>> Anew(nl, std::vector<T>(M, zero));
        for (std::size_t v = 0; v < nl; ++v) {
            const std::vector<std::size_t>& V = latt[v];
            std::size_t vv = 0;
            for (std::size_t j = 0; j < J; ++j) vv += V[j];
            const double omr = detail::sdr_omega_cum(c, tt, vv);
            if (omr == 0.0) continue;
            std::vector<T> anum(M, zero);
            const std::vector<std::vector<std::size_t>> sub = detail::sdr_lattice(V);
            for (std::size_t q = 0; q < sub.size(); ++q) {
                std::vector<std::size_t> VL(J);
                for (std::size_t j = 0; j < J; ++j) VL[j] = V[j] - sub[q][j];
                const std::size_t iL = detail::sdr_key(sub[q], N), iVL = detail::sdr_key(VL, N);
                const T pb = num_traits<T>::from_double(omr) * lev.g[iVL] * Gin[iL];
                if (pb == zero) continue;
                Gnew[v] = Gnew[v] + pb;
                for (std::size_t k = 0; k < St.size(); ++k)
                    for (std::size_t j = 0; j < J; ++j)
                        Qnew[v](St[k], j) = Qnew[v](St[k], j) + lev.Q[iVL](St[k], j) * pb;
                for (std::size_t j = 0; j < J; ++j) Tnew[v][j] = Tnew[v][j] + lev.Tp[iVL][j] * pb;
                for (std::size_t q2 = 0; q2 < inner.size(); ++q2) {
                    const std::size_t i = inner[q2];
                    for (std::size_t j = 0; j < J; ++j)
                        Qnew[v](i, j) = Qnew[v](i, j) + Qin[iL](i, j) * pb;
                    anum[i] = anum[i] + Ain[iL][i] * pb;
                }
            }
            if (Gnew[v] > zero) {
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t j = 0; j < J; ++j) Qnew[v](i, j) = Qnew[v](i, j) / Gnew[v];
                for (std::size_t j = 0; j < J; ++j) Tnew[v][j] = Tnew[v][j] / Gnew[v];
                // eq. (10) carries, besides delta_ti, the single-step omega ratios
                // down to the centre's own level. Conditioning on the population of
                // Q(V,V_t) fixes v_t, so that ratio leaves the expectation and the
                // rest recurses through the same convolution as the queue lengths.
                const T st = num_traits<T>::from_double(detail::sdr_omega_step(c, tt, vv));
                for (std::size_t k = 0; k < St.size(); ++k) {
                    T qi = zero;
                    for (std::size_t j = 0; j < J; ++j) qi = qi + Qnew[v](St[k], j);
                    Anew[v][St[k]] = st * (num_traits<T>::from_double(dvec[St[k]]) - qi);
                }
                for (std::size_t q2 = 0; q2 < inner.size(); ++q2)
                    Anew[v][inner[q2]] = st * anum[inner[q2]] / Gnew[v];
            }
        }
        Gin = Gnew;
        Qin = Qnew;
        Ain = Anew;
        if (tt == 1) break;
    }

    // Sec. 4.3, at every population because the throughputs read N - 1_j
    std::vector<Matrix<T>> Qall(nl, Matrix<T>(M, J, zero));
    std::vector<std::vector<T>> Tall(nl, std::vector<T>(J, zero));
    std::vector<std::vector<T>> Aall(nl, std::vector<T>(M, zero));
    std::vector<T> Gall(nl, zero);
    for (std::size_t v = 0; v < nl; ++v) {
        const std::vector<std::size_t>& Np = latt[v];
        const std::vector<std::vector<std::size_t>> sub = detail::sdr_lattice(Np);
        for (std::size_t q = 0; q < sub.size(); ++q) {
            std::vector<std::size_t> C2(J);
            for (std::size_t j = 0; j < J; ++j) C2[j] = Np[j] - sub[q][j];
            const std::size_t iV = detail::sdr_key(sub[q], N), iC = detail::sdr_key(C2, N);
            const T pb = comp.g[iC] * Gin[iV];
            if (pb == zero) continue;
            Gall[v] = Gall[v] + pb;
            for (std::size_t k = 0; k < mv.size(); ++k)
                for (std::size_t j = 0; j < J; ++j)
                    Qall[v](mv[k], j) = Qall[v](mv[k], j) + comp.Q[iC](mv[k], j) * pb;
            for (std::size_t i = 0; i < M; ++i) {
                if (!inV[i]) continue;
                for (std::size_t j = 0; j < J; ++j)
                    Qall[v](i, j) = Qall[v](i, j) + Qin[iV](i, j) * pb;
                Aall[v][i] = Aall[v][i] + Ain[iV][i] * pb;
            }
            for (std::size_t j = 0; j < J; ++j) Tall[v][j] = Tall[v][j] + comp.Tp[iC][j] * pb;
        }
        if (Gall[v] > zero) {
            for (std::size_t i = 0; i < M; ++i) {
                for (std::size_t j = 0; j < J; ++j) Qall[v](i, j) = Qall[v](i, j) / Gall[v];
                Aall[v][i] = Aall[v][i] / Gall[v];
            }
            for (std::size_t j = 0; j < J; ++j) Tall[v][j] = Tall[v][j] / Gall[v];
        }
    }

    const std::size_t vN = detail::sdr_key(N, N);
    SdrResult<T> res;
    res.QN = Qall[vN];
    res.XN = Matrix<T>(M, J, zero);
    res.UN = Matrix<T>(M, J, zero);
    res.RN = Matrix<T>(M, J, zero);
    res.G = Gall[vN];
    for (std::size_t j = 0; j < J; ++j) {
        if (N[j] == 0) continue;
        std::vector<std::size_t> Nm = N;
        Nm[j] -= 1;
        const std::size_t vm = detail::sdr_key(Nm, N);
        for (std::size_t i = 0; i < M; ++i)
            res.XN(i, j) = inV[i] ? xi(i, j) * Tall[vN][j] * Aall[vm][i]
                                  : xi(i, j) * Tall[vN][j];
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < J; ++j) {
            res.UN(i, j) = res.XN(i, j) * S(i, j);
            if (res.XN(i, j) > zero) res.RN(i, j) = res.QN(i, j) / res.XN(i, j);
        }
    return res;
}

/**
 * Coefficients xi of Section 3.2.
 *
 * P holds one centre-by-centre state-independent routing matrix per chain.
 * Three rules fix the coefficients: the complement M-V obeys the ordinary
 * traffic equations with the whole SDR subnetwork collapsed into a single
 * e -> d arc of probability one; every branch obeys its own traffic equations
 * driven by an injection of xi_e at its entry centre; and xi_e is one.
 *
 * The paper states xi_ij = xi_ej for the branch entry and departure centres and
 * works out only single-centre branches. The traffic equations above are the
 * reading that extends it: they return xi at the branch departure equal to
 * xi_e because a customer leaves a branch only through it, and xi at the branch
 * entry equal to xi_e whenever that centre takes no internal feedback. They have
 * been checked against a brute-force CTMC on a branch that does take such
 * feedback, where the literal rule fails.
 *
 * These xi are NOT relative visit counts: the rate at which customers enter a
 * branch is state dependent, so a ratio of two xi carries no flow meaning.
 */
template <class T>
Matrix<T> pfqn_sdrvisits(const SdrStruct& sdr, const std::vector<Matrix<T>>& P) {
    const SdrCoeff c = pfqn_sdrcoeff(sdr);
    if (P.empty()) throw InputError("pfqn_sdrvisits: no routing matrix was supplied");
    const std::size_t M = P[0].rows(), J = P.size();
    if (P[0].cols() != M)
        throw InputError("pfqn_sdrvisits: the SIR routing matrices must be square");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    Matrix<T> xi(M, J, zero);

    std::vector<bool> inV(M, false);
    for (std::size_t b = 1; b < c.B; ++b)
        for (std::size_t k = 0; k < sdr.branch[b].size(); ++k) inV[sdr.branch[b][k]] = true;
    std::vector<std::size_t> mv;
    for (std::size_t i = 0; i < M; ++i)
        if (!inV[i]) mv.push_back(i);
    std::size_t ie = mv.size(), id = mv.size();
    for (std::size_t k = 0; k < mv.size(); ++k) {
        if (mv[k] == sdr.entry) ie = k;
        if (mv[k] == sdr.departure) id = k;
    }
    if (ie == mv.size() || id == mv.size())
        throw InputError("pfqn_sdrvisits: the entry and departure centres of Q(V,V) must lie "
                         "outside every branch");

    for (std::size_t j = 0; j < J; ++j) {
        const std::size_t nm = mv.size();
        // Complement M-V with the SDR subnetwork collapsed into the arc e -> d
        Matrix<T> Pmv(nm, nm, zero);
        for (std::size_t a = 0; a < nm; ++a) {
            if (a == ie) {
                Pmv(a, id) = one;
                continue;
            }
            T rs = zero;
            for (std::size_t b = 0; b < nm; ++b) {
                Pmv(a, b) = P[j](mv[a], mv[b]);
                rs = rs + Pmv(a, b);
            }
            if (!(rs > zero))
                throw InputError("pfqn_sdrvisits: the SIR routing does not keep customers inside "
                                 "the complement M-V");
        }
        std::vector<T> xmv = mc::dtmc_solve(Pmv);
        if (!(xmv[ie] > zero))
            throw InputError("pfqn_sdrvisits: the entry centre of Q(V,V) is unreachable");
        for (std::size_t a = 0; a < nm; ++a) xi(mv[a], j) = xmv[a] / xmv[ie];

        // Each branch, driven by an injection of xi_e at its entry centre:
        // xi_Bb = xi_e e_{e(b)} (I - P_bb)^-1, written as the transposed solve
        for (std::size_t b = 1; b < c.B; ++b) {
            const std::vector<std::size_t>& sb = sdr.branch[b];
            const std::size_t nb = sb.size();
            Matrix<T> Ab(nb, nb, zero);
            std::vector<T> rb(nb, zero);
            for (std::size_t a = 0; a < nb; ++a) {
                for (std::size_t k = 0; k < nb; ++k)
                    Ab(k, a) = ((a == k) ? one : zero) - P[j](sb[a], sb[k]);
                if (sb[a] == sdr.entryOf[b]) rb[a] = xi(sdr.entry, j);
            }
            std::vector<T> sol = solve(Ab, rb);
            for (std::size_t a = 0; a < nb; ++a) xi(sb[a], j) = sol[a];
        }
    }
    return xi;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_SDR_H
