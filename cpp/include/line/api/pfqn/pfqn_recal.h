/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_RECAL_H
#define LINE_API_PFQN_RECAL_H

/**
 * RECAL (REcursive CALculation) for the exact normalizing constant of a closed
 * product-form network (Conway and Georganas 1986).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_recal.m, cross-checked against
 * mp_pfqn's recal/recal-multi-exact.c for the exact path.
 *
 * Where convolution recurses on stations over the population lattice, RECAL
 * recurses on jobs over the space of station multiplicity vectors. Jobs are
 * added one at a time, class 0 first, and the state carried between steps is a
 * function g_n indexed by a multiplicity vector m with sum(m) = Ntot - n:
 *
 *   g_0(m)   = 1                                          for all sum(m) = Ntot
 *   g_n(m)   = ( Z_r g_{n-1}(m + e_delay)
 *                + sum_j (m_j + 1 + m0_j - 1) L(j,r) g_{n-1}(m + e_j) ) / n_r
 *
 * where r is the class of the n-th job, n_r its index within that class, and
 * m0(j) the multiplicity of station j. The answer is g_Ntot(0), the single
 * state with an empty multiplicity vector. Every operation is a field
 * operation (the only division is by the small integer n_r), so the recursion
 * is exact in rational arithmetic with no reformulation, which is what makes
 * RECAL usable as an exactness oracle for the rest of the pfqn family.
 *
 * Delay column: think time enters as one extra column of m, appended after the
 * M queueing stations. It is allocated only when some class has a non-zero
 * think time, following mp_pfqn. MATLAB always carries the column; with Z = 0
 * the column is a spectator (g never reads across it and the level-0 value is
 * 1 everywhere), so the returned constant is identical and only the size of
 * the intermediate arrays differs.
 *
 * Station consolidation: stations with identical demand rows are merged and
 * their multiplicities added, which is exactly what the (m_j + m0_j - 1)
 * coefficient is for. This mirrors pfqn_unique in the MATLAB reference and
 * shrinks the state space from multichoose(M+1, Ntot) to multichoose(M'+1,
 * Ntot). Rows are compared for exact equality rather than with the MATLAB
 * 1e-14 tolerance: merging rows that only nearly agree would perturb the
 * constant, which an exact-capable algorithm must not do.
 *
 * Scaling: as for convolution, in IEEE double the recursion can leave the
 * exponent range, so the same power-of-two rescaling of pfqn_ca is applied
 * there and nowhere else. G is homogeneous of degree Ntot in (L, Z), so
 * dividing every demand and think time by 2^k divides G by exactly 2^(k*Ntot),
 * and the recovery is an exponent adjustment rather than an exp().
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include "line/api/pfqn/pfqn_ca.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

namespace detail {

/** Largest number of multiplicity vectors the port will allocate g over. */
constexpr unsigned long long RECAL_MAX_STATES = 100000000ULL;

/**
 * Binomial coefficient as an exact integer. The running product is a binomial
 * at every step, so the division is exact and no rounding is possible; an
 * overflow of the accumulator is reported rather than wrapped.
 */
inline unsigned long long binom_exact(unsigned long long n, unsigned long long k) {
    if (k > n) return 0ULL;
    if (k > n - k) k = n - k;
    unsigned long long r = 1ULL;
    for (unsigned long long i = 1ULL; i <= k; ++i) {
        const unsigned long long a = n - k + i;
        if (r > std::numeric_limits<unsigned long long>::max() / a)
            throw NumericError("pfqn_recal: multiplicity state space overflows a 64-bit count");
        r = r * a / i;
    }
    return r;
}

/** Number of vectors of ncols non-negative integers summing to k. */
inline unsigned long long multichoose_count(std::size_t ncols, int k) {
    if (k < 0) return 0ULL;
    if (ncols == 0) return k == 0 ? 1ULL : 0ULL;
    return binom_exact(static_cast<unsigned long long>(ncols) + static_cast<unsigned long long>(k) - 1ULL,
                       static_cast<unsigned long long>(k));
}

/**
 * Position of m in the enumeration of all ncols-vectors summing to ksum,
 * ordered by the first component ascending and then recursively on the tail.
 * O(ncols + ksum), replacing the linear row scan of the MATLAB matchrow; the
 * enumeration it ranks is the one recal_next_composition walks.
 */
inline std::size_t recal_rank(const std::vector<int>& m, std::size_t ncols, int ksum) {
    unsigned long long idx = 0ULL;
    std::size_t pos = 0;
    while (ncols > 1) {
        const int c = m[pos];
        for (int i = 0; i < c; ++i) idx += multichoose_count(ncols - 1, ksum - i);
        ksum -= c;
        --ncols;
        ++pos;
    }
    return static_cast<std::size_t>(idx);
}

/**
 * Advance m to the next vector of the same length and sum, in the order that
 * recal_rank ranks. Returns false once the enumeration is exhausted.
 */
inline bool recal_next_composition(std::vector<int>& m) {
    const std::size_t n = m.size();
    if (n < 2) return false;
    const std::size_t last = n - 1;
    long rest = 0;
    for (long p = static_cast<long>(n) - 2; p >= 0; --p) {
        rest += m[static_cast<std::size_t>(p) + 1];
        if (rest > 0) {
            m[static_cast<std::size_t>(p)] += 1;
            for (std::size_t j = static_cast<std::size_t>(p) + 1; j < last; ++j) m[j] = 0;
            m[last] = static_cast<int>(rest - 1);
            return true;
        }
    }
    return false;
}

/**
 * Merge stations whose demand rows are identical, summing their
 * multiplicities. Mirrors pfqn_unique followed by the m0 accumulation of the
 * MATLAB reference, with exact row equality as the merge test.
 */
template <class T>
void consolidate_stations(const Matrix<T>& L, const std::vector<int>& m0, Matrix<T>& Lu,
                          std::vector<int>& m0u) {
    const std::size_t M = L.rows(), R = L.cols();
    std::vector<std::size_t> keep;
    std::vector<std::size_t> mapping(M, 0);
    for (std::size_t i = 0; i < M; ++i) {
        std::size_t hit = keep.size();
        for (std::size_t u = 0; u < keep.size(); ++u) {
            bool same = true;
            for (std::size_t r = 0; r < R; ++r)
                if (!(L(i, r) == L(keep[u], r))) {
                    same = false;
                    break;
                }
            if (same) {
                hit = u;
                break;
            }
        }
        if (hit == keep.size()) keep.push_back(i);
        mapping[i] = hit;
    }
    Lu = Matrix<T>(keep.size(), R);
    for (std::size_t u = 0; u < keep.size(); ++u)
        for (std::size_t r = 0; r < R; ++r) Lu(u, r) = L(keep[u], r);
    m0u.assign(keep.size(), 0);
    for (std::size_t i = 0; i < M; ++i) m0u[mapping[i]] += m0[i];
}

}  // namespace detail

/**
 * @param L  (M x R) service demands, M queueing stations, R classes
 * @param N  (R) population per class, non-negative
 * @param Z  (K x R) think times, summed over rows; may be empty
 * @param m0 (M) station multiplicities, each at least one; empty for all ones
 */
template <class T>
NcResult<T> pfqn_recal(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                       const std::vector<int>& m0) {
    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_recal: demand matrix and population vector disagree on the class count");
    if (!m0.empty() && m0.size() != M)
        throw InputError("pfqn_recal: multiplicity vector has the wrong length");
    for (std::size_t i = 0; i < m0.size(); ++i)
        if (m0[i] < 1) throw InputError("pfqn_recal: station multiplicity below one");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    // Z summed over its rows, so a per-node think-time matrix is accepted.
    std::vector<T> Zsum(R, zero);
    if (!Z.empty()) {
        if (Z.cols() != R) throw InputError("pfqn_recal: Z and N disagree on the class count");
        for (std::size_t k = 0; k < Z.rows(); ++k)
            for (std::size_t r = 0; r < R; ++r) Zsum[r] += Z(k, r);
    }

    // negative-population rejection rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    long Nt = 0;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_recal: negative population");
        Nt += v;
    }

    if (M == 0) {
        // Delay-only network: G = prod_r Z_r^{N_r} / N_r!.
        const T G = detail::pff_delay(Zsum, N);
        return {G, num_traits<T>::log_as_double(G)};
    }
    if (Nt == 0) return {one, 0.0};

    std::vector<int> mult(M, 1);
    for (std::size_t i = 0; i < m0.size(); ++i) mult[i] = m0[i];
    Matrix<T> Lc;
    std::vector<int> m0c;
    detail::consolidate_stations(L, mult, Lc, m0c);
    const std::size_t Mq = Lc.rows();

    const int kscale = detail::scale_exponent(Lc, N, Zsum);
    if constexpr (std::is_same<T, double>::value) {
        // exponent-only rescaling rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
        if (kscale != 0) {
            for (std::size_t i = 0; i < Mq; ++i)
                for (std::size_t r = 0; r < R; ++r) Lc(i, r) = std::ldexp(Lc(i, r), -kscale);
            for (std::size_t r = 0; r < R; ++r) Zsum[r] = std::ldexp(Zsum[r], -kscale);
        }
    }

    bool hasZ = false;
    for (std::size_t r = 0; r < R; ++r)
        if (!(Zsum[r] == zero)) {
            hasZ = true;
            break;
        }
    // Columns of the multiplicity vector: the queueing stations, plus one
    // delay column when any class thinks.
    const std::size_t Mz = hasZ ? Mq + 1 : Mq;
    const std::size_t delay = Mq;  // column index of the delay slot, if present

    const unsigned long long states = detail::multichoose_count(Mz, static_cast<int>(Nt));
    if (states > detail::RECAL_MAX_STATES)
        throw NumericError("pfqn_recal: multiplicity state space too large for this model");

    // g at the previous and the current job count. Level 0 is the largest, so a
    // single pair of buffers of that size serves every level.
    std::vector<T> gprev(static_cast<std::size_t>(states), one);
    std::vector<T> gcur(static_cast<std::size_t>(states), zero);

    std::vector<int> m(Mz, 0);
    int n = 0;
    for (std::size_t r = 0; r < R; ++r) {
        for (int nr = 1; nr <= N[r]; ++nr) {
            ++n;
            const int k = static_cast<int>(Nt) - n;  // sum of the current level
            const int kprev = k + 1;                 // sum of the previous level
            const unsigned long long ncfg = detail::multichoose_count(Mz, k);
            const T nrv = num_traits<T>::from_int(nr);
            const bool thinks = hasZ && !(Zsum[r] == zero);

            m.assign(Mz, 0);
            m[Mz - 1] = k;
            for (unsigned long long i = 0; i < ncfg; ++i) {
                T acc = zero;
                if (thinks) {
                    m[delay] += 1;
                    acc += Zsum[r] * gprev[detail::recal_rank(m, Mz, kprev)];
                    m[delay] -= 1;
                }
                for (std::size_t j = 0; j < Mq; ++j) {
                    if (Lc(j, r) == zero) continue;  // the term is exactly zero
                    m[j] += 1;
                    acc += num_traits<T>::from_int(m[j] + m0c[j] - 1) * Lc(j, r) *
                           gprev[detail::recal_rank(m, Mz, kprev)];
                    m[j] -= 1;
                }
                gcur[static_cast<std::size_t>(i)] = acc / nrv;
                if (i + 1 < ncfg) detail::recal_next_composition(m);
            }
            gprev.swap(gcur);
        }
    }

    // The final level holds the single empty multiplicity vector.
    const T raw = gprev[0];
    const double lG =
        num_traits<T>::log_as_double(raw) + static_cast<double>(Nt) * kscale * std::log(2.0);
    T G = raw;
    if constexpr (std::is_same<T, double>::value) {
        if (kscale != 0) G = std::ldexp(raw, static_cast<int>(Nt * kscale));
    }
    return {G, lG};
}

/** Overload with unit station multiplicities. */
template <class T>
NcResult<T> pfqn_recal(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z) {
    return pfqn_recal(L, N, Z, std::vector<int>());
}

/** Overload without think times. */
template <class T>
NcResult<T> pfqn_recal(const Matrix<T>& L, const std::vector<int>& N) {
    return pfqn_recal(L, N, Matrix<T>(), std::vector<int>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_RECAL_H
