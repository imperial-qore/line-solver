/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_RETRIEVAL_RETRIEVAL_FPI_LATENCY_H
#define LINE_API_RETRIEVAL_RETRIEVAL_FPI_LATENCY_H

/**
 * FPI-based approximation of the delayed-hit count and the expected latency of
 * a list-based cache with a phase-type retrieval system.
 *
 * Templated port of matlab/src/api/retrieval/retrieval_fpi_latency.m,
 * cross-checked against
 * jar/src/main/java/jline/api/retrieval/Retrieval_fpi_latency.java.
 *
 * From the fetch-period duration F_i of item i,
 *   d_i = phi_i lambda_i E0[F_i^2] / (2 E0[F_i]),
 * with the Palm moments of a reduced absorbing CTMC of item i's visits to the
 * retrieval stations, E0[F_i^k] = k! pi_e (-D0)^{-k} e. The steps are
 *   1. retrieval_fpi on the whole system            -> phi_i, pi_{i,0}
 *   2. retrieval_fpi without item i                 -> phitilde_s, the PS occupancy
 *   3. one PH block per station, shared stations slowed by 1/(1+phitilde_s),
 *      routed by R, absorbing on return to the cache
 *   4. the two moments, hence d_i
 *   5. Z = sum_i (phi_i + d_i) / sum_i lambda_i (phi_i + pi_{i,0})
 *
 * Routing convention, as in MATLAB: index 0 is the outside (entry on a miss,
 * return to the cache on completion) and indices 1..S are the retrieval
 * stations, R[i](a,b) being the probability of a -> b for item i.
 *
 * Station types: IS fetches are independent; PS and LCFSPR are symmetric
 * insensitive disciplines and admit general phase-type, class-dependent
 * service; SIRO and FCFS reduce to the same single-exponential sojourn only
 * with exponential service at a class-independent rate, and are rejected
 * otherwise. Any other discipline is rejected outright, matching MATLAB.
 *
 * ARITHMETIC: it calls retrieval_fpi, a tolerance-stopped successive
 * substitution, so it is gated on has_transcendental for the same reason. The
 * linear algebra around it (the visit-ratio solve and the two moment solves)
 * would itself be exact, but the phitilde it is fed is not.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/retrieval/retrieval_fpi.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace retrieval {

/** Scheduling of a retrieval station, mirroring the MATLAB station_type strings. */
enum class RetrievalStationType { IS, PS, SIRO, FCFS, LCFSPR };

/** True for the disciplines that carry the mean-field sharing slowdown. */
inline bool is_shared_station(RetrievalStationType t) {
    return t != RetrievalStationType::IS;
}

/** Mirrors the [Z, d, phi, pi0] return list of the MATLAB function. */
template <class T>
struct RetrievalFpiLatencyResult {
    T Z;                 ///< expected latency of the delayed-hit system
    std::vector<T> d;    ///< (n) mean delayed hits awaiting fetch, per item
    std::vector<T> phi;  ///< (n) delayed-hit ratio phi_i
    std::vector<T> pi0;  ///< (n) miss ratio pi_{i,0}
};

/**
 * Phase-type service of every item at one retrieval station.
 *
 * MATLAB carries this as alpha{s}(1,:,i) and T{s}(:,:,i); here alpha is one
 * (n x f) matrix whose row i is the entry vector of item i, and sub is one
 * (f x f) subgenerator per item.
 */
template <class T>
struct RetrievalStationPH {
    Matrix<T> alpha;             ///< (n x f) entry vectors, row per item
    std::vector<Matrix<T>> sub;  ///< (n) subgenerators, (f x f) each
    RetrievalStationType type = RetrievalStationType::PS;

    std::size_t phases() const { return alpha.cols(); }
};

namespace detail {

/** Mean phase-type service time -alpha inv(Tm) 1. */
template <class T>
T ph_mean(const Matrix<T>& alpha_row, const Matrix<T>& Tm) {
    const std::size_t f = Tm.rows();
    std::vector<T> e(f, num_traits<T>::from_int(1));
    Matrix<T> LU = Tm;
    const std::vector<std::size_t> piv = lu_factor(LU);
    lu_solve(LU, piv, e);  // e <- inv(Tm) 1
    T tau = num_traits<T>::from_int(0);
    for (std::size_t k = 0; k < f; ++k) tau += alpha_row(0, k) * e[k];
    return -tau;
}

}  // namespace detail

/**
 * @param m       (h) cache list capacities
 * @param lambda  (n) per-item arrival rates
 * @param gamma   (n x h) access factors
 * @param station (S) phase-type service and discipline of each retrieval station
 * @param R       (n) routing matrices, each (S+1) x (S+1), index 0 = outside
 * @param options fixed-point options (tolerance, iteration cap, damping)
 */
template <class T>
RetrievalFpiLatencyResult<T> retrieval_fpi_latency(const std::vector<int>& m,
                                                   const std::vector<T>& lambda,
                                                   const Matrix<T>& gamma,
                                                   const std::vector<RetrievalStationPH<T>>& station,
                                                   const std::vector<Matrix<T>>& R,
                                                   const FpiOptions& options = FpiOptions()) {
    static_assert(num_traits<T>::has_transcendental,
                  "retrieval_fpi_latency requires transcendental arithmetic: it is driven by "
                  "retrieval_fpi, a successive substitution stopped on a relative tolerance");
    const std::size_t n = lambda.size();
    const std::size_t S = station.size();
    if (S == 0) throw InputError("retrieval_fpi_latency: no retrieval station");
    if (gamma.rows() != n)
        throw InputError("retrieval_fpi_latency: gamma and lambda disagree on the item count");
    if (R.size() != n)
        throw InputError("retrieval_fpi_latency: one routing matrix per item is required");
    for (std::size_t i = 0; i < n; ++i)
        if (R[i].rows() != S + 1 || R[i].cols() != S + 1)
            throw InputError("retrieval_fpi_latency: routing matrices must be (S+1) x (S+1)");

    std::vector<std::size_t> fsz(S);
    for (std::size_t s = 0; s < S; ++s) {
        fsz[s] = station[s].phases();
        if (station[s].sub.size() != n)
            throw InputError("retrieval_fpi_latency: one subgenerator per item is required");
        if (station[s].alpha.rows() != n)
            throw InputError("retrieval_fpi_latency: alpha must have one row per item");
        for (std::size_t i = 0; i < n; ++i)
            if (station[s].sub[i].rows() != fsz[s] || station[s].sub[i].cols() != fsz[s])
                throw InputError("retrieval_fpi_latency: subgenerator size disagrees with alpha");
        // SIRO/FCFS PS-collapse rationale: see _kb/03-api-layer.md (cpp port notes: retrieval)
        if ((station[s].type == RetrievalStationType::SIRO ||
             station[s].type == RetrievalStationType::FCFS) &&
            fsz[s] > 1)
            throw UnsupportedError(
                "retrieval_fpi_latency: SIRO/FCFS retrieval stations require exponential "
                "(single-phase) service");
    }

    // per-item mean service times at each station
    Matrix<T> tau(n, S, num_traits<T>::from_int(0));
    for (std::size_t s = 0; s < S; ++s)
        for (std::size_t i = 0; i < n; ++i) {
            Matrix<T> arow(1, fsz[s]);
            for (std::size_t k = 0; k < fsz[s]; ++k) arow(0, k) = station[s].alpha(i, k);
            tau(i, s) = detail::ph_mean(arow, station[s].sub[i]);
        }

    // ... and the class-independence requirement at SIRO/FCFS stations
    for (std::size_t s = 0; s < S; ++s) {
        if (station[s].type != RetrievalStationType::SIRO &&
            station[s].type != RetrievalStationType::FCFS)
            continue;
        double lo = 0.0, hi = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double x = num_traits<T>::to_double(tau(i, s));
            if (i == 0 || x < lo) lo = x;
            if (i == 0 || x > hi) hi = x;
        }
        if (hi - lo > 1e-9 * hi)
            throw UnsupportedError(
                "retrieval_fpi_latency: SIRO/FCFS retrieval stations require class-independent "
                "mean service rates");
    }

    // IS stations aggregate into column 0 of eta; each shared station gets its own column.
    std::vector<std::size_t> is_idx, ps_idx;
    for (std::size_t s = 0; s < S; ++s) {
        if (station[s].type == RetrievalStationType::IS)
            is_idx.push_back(s);
        else
            ps_idx.push_back(s);
    }
    const std::size_t r = ps_idx.size();

    // eta_{s,i} = (visits per fetch) * (mean service time)
    Matrix<T> eta(n, r + 1, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i) {
        Matrix<T> ImP(S, S);
        for (std::size_t a = 0; a < S; ++a)
            for (std::size_t b = 0; b < S; ++b)
                ImP(a, b) = (a == b ? num_traits<T>::from_int(1) : num_traits<T>::from_int(0)) -
                            R[i](a + 1, b + 1);
        const Matrix<T> Vinv = inverse(ImP);
        std::vector<T> visits(S, num_traits<T>::from_int(0));
        for (std::size_t b = 0; b < S; ++b)
            for (std::size_t a = 0; a < S; ++a) visits[b] += R[i](0, a + 1) * Vinv(a, b);

        T is_sum = num_traits<T>::from_int(0);
        for (std::size_t s : is_idx) is_sum += visits[s] * tau(i, s);
        eta(i, 0) = is_sum;
        for (std::size_t p = 0; p < r; ++p) eta(i, 1 + p) = visits[ps_idx[p]] * tau(i, ps_idx[p]);
    }

    // step 1: FPI on the full system
    const RetrievalFpiResult<T> full = retrieval_fpi(m, lambda, eta, gamma, options);
    RetrievalFpiLatencyResult<T> out;
    out.pi0 = full.pmiss;
    out.phi.assign(n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t s = 0; s <= r; ++s) out.phi[i] += full.pdh(s, i);

    out.d.assign(n, num_traits<T>::from_int(0));
    const std::size_t Phi = [&]() {
        std::size_t t = 0;
        for (std::size_t s = 0; s < S; ++s) t += fsz[s];
        return t;
    }();
    std::vector<std::size_t> off(S, 0);
    for (std::size_t s = 1; s < S; ++s) off[s] = off[s - 1] + fsz[s - 1];

    for (std::size_t i = 0; i < n; ++i) {
        // step 2: FPI without item i, giving the occupancy left by the others
        std::vector<T> lambda_i;
        lambda_i.reserve(n - 1);
        Matrix<T> eta_i(n - 1, r + 1);
        Matrix<T> gamma_i(n - 1, gamma.cols());
        std::size_t q = 0;
        for (std::size_t k = 0; k < n; ++k) {
            if (k == i) continue;
            lambda_i.push_back(lambda[k]);
            for (std::size_t c = 0; c <= r; ++c) eta_i(q, c) = eta(k, c);
            for (std::size_t c = 0; c < gamma.cols(); ++c) gamma_i(q, c) = gamma(k, c);
            ++q;
        }
        std::vector<T> phitilde(S, num_traits<T>::from_int(0));
        if (n > 1) {
            const RetrievalFpiResult<T> without =
                retrieval_fpi(m, lambda_i, eta_i, gamma_i, options);
            for (std::size_t p = 0; p < r; ++p) {
                T acc = num_traits<T>::from_int(0);
                for (std::size_t k = 0; k + 1 < n; ++k) acc += without.pdh(1 + p, k);
                phitilde[ps_idx[p]] = acc;
            }
        }

        // step 3: the reduced absorbing CTMC of item i's fetch
        Matrix<T> D0(Phi, Phi, num_traits<T>::from_int(0));
        std::vector<T> pe(Phi, num_traits<T>::from_int(0));
        for (std::size_t s = 0; s < S; ++s) {
            const T scale = is_shared_station(station[s].type)
                                ? num_traits<T>::from_int(1) /
                                      (num_traits<T>::from_int(1) + phitilde[s])
                                : num_traits<T>::from_int(1);
            Matrix<T> blk(fsz[s], fsz[s]);
            for (std::size_t a = 0; a < fsz[s]; ++a)
                for (std::size_t b = 0; b < fsz[s]; ++b) blk(a, b) = scale * station[s].sub[i](a, b);
            for (std::size_t a = 0; a < fsz[s]; ++a)
                for (std::size_t b = 0; b < fsz[s]; ++b) D0(off[s] + a, off[s] + b) += blk(a, b);
            // completion rates out of each phase of station s
            std::vector<T> compl_(fsz[s], num_traits<T>::from_int(0));
            for (std::size_t a = 0; a < fsz[s]; ++a) {
                T acc = num_traits<T>::from_int(0);
                for (std::size_t b = 0; b < fsz[s]; ++b) acc += blk(a, b);
                compl_[a] = -acc;
            }
            for (std::size_t sp = 0; sp < S; ++sp)
                for (std::size_t a = 0; a < fsz[s]; ++a)
                    for (std::size_t b = 0; b < fsz[sp]; ++b)
                        D0(off[s] + a, off[sp] + b) +=
                            compl_[a] * R[i](s + 1, sp + 1) * station[sp].alpha(i, b);
            for (std::size_t a = 0; a < fsz[s]; ++a)
                pe[off[s] + a] = R[i](0, s + 1) * station[s].alpha(i, a);
        }

        // step 4: the two Palm moments, hence d_i
        Matrix<T> A(Phi, Phi);
        for (std::size_t a = 0; a < Phi; ++a)
            for (std::size_t b = 0; b < Phi; ++b) A(a, b) = -D0(a, b);
        Matrix<T> LU = A;
        const std::vector<std::size_t> piv = lu_factor(LU);
        std::vector<T> y(Phi, num_traits<T>::from_int(1));
        lu_solve(LU, piv, y);  // y = inv(A) 1
        std::vector<T> y2 = y;
        lu_solve(LU, piv, y2);  // y2 = inv(A) y
        T M1 = num_traits<T>::from_int(0), M2 = num_traits<T>::from_int(0);
        for (std::size_t a = 0; a < Phi; ++a) {
            M1 += pe[a] * y[a];
            M2 += pe[a] * y2[a];
        }
        M2 *= num_traits<T>::from_int(2);
        if (M1 == num_traits<T>::from_int(0))
            throw NumericError("retrieval_fpi_latency: zero mean fetch period");
        out.d[i] = out.phi[i] * lambda[i] * M2 / (num_traits<T>::from_int(2) * M1);
    }

    // step 5: the expected latency
    T num = num_traits<T>::from_int(0), den = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        num += out.phi[i] + out.d[i];
        den += lambda[i] * (out.phi[i] + out.pi0[i]);
    }
    if (den == num_traits<T>::from_int(0))
        throw NumericError("retrieval_fpi_latency: no fetching traffic, the latency is undefined");
    out.Z = num / den;
    return out;
}

}  // namespace retrieval
}  // namespace line

#endif  // LINE_API_RETRIEVAL_RETRIEVAL_FPI_LATENCY_H
