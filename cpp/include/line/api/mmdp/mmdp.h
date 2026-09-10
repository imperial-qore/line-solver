/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MMDP_MMDP_H
#define LINE_API_MMDP_MMDP_H

/**
 * Markov-modulated deterministic process (MMDP), for fluid queues.
 *
 * Port of python/line_solver/api/mmdp/__init__.py. PYTHON-ONLY: there is no
 * MATLAB or JAR twin, so native Python is the reference.
 *
 * MMDP is the DETERMINISTIC analogue of MMPP. An MMPP modulates a Poisson
 * ARRIVAL RATE by a background chain; an MMDP modulates a deterministic FLUID
 * FLOW RATE by one. The parameterization is BUTools': `Q` is the generator of
 * the modulating chain (rows summing to zero) and `R` is the DIAGONAL matrix of
 * per-state flow rates.
 *
 * WHY THE SCV HERE IS NOT AN INTERARRIVAL SCV. In an MMPP the variability of
 * interest is that of the interarrival time; in an MMDP the flow is
 * deterministic within a state, so all the variability lives in WHICH state the
 * chain occupies. `mmdp_scv` is therefore the SCV of the RATE under the
 * stationary law, `Var[r]/E[r]^2` with `r` the per-state rate -- not of any
 * holding time. A two-state process with equal rates has SCV zero however
 * fast it switches, which is the tell that this is the rate's dispersion and
 * not a time's.
 *
 * ARITHMETIC: field. The stationary solve is the only numeric step and it goes
 * through `ctmc_solve`, which is templated.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mmdp {

/**
 * True when (Q, R) is a valid MMDP.
 *
 * Q must be a generator -- square, non-positive diagonal, non-negative
 * off-diagonal, rows summing to zero -- and R must be DIAGONAL with
 * non-negative entries. The diagonality is not a storage convention: a
 * non-diagonal R would make the flow rate depend on a transition rather than
 * on the state, which is a different process.
 */
template <class T>
bool mmdp_isfeasible(const Matrix<T>& Q, const Matrix<T>& R, double tol = 1e-10) {
    const std::size_t n = Q.rows();
    if (n == 0 || Q.cols() != n) return false;
    if (R.rows() != n || R.cols() != n) return false;
    for (std::size_t i = 0; i < n; ++i) {
        if (num_traits<T>::to_double(Q(i, i)) > tol) return false;
        double rowsum = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            const double q = num_traits<T>::to_double(Q(i, j));
            if (i != j && q < -tol) return false;
            rowsum += q;
        }
        if (std::fabs(rowsum) > tol) return false;
        if (num_traits<T>::to_double(R(i, i)) < -tol) return false;
        for (std::size_t j = 0; j < n; ++j)
            if (i != j && std::fabs(num_traits<T>::to_double(R(i, j))) > tol) return false;
    }
    return true;
}

/** The per-state rates, i.e. the diagonal of R. */
template <class T>
std::vector<T> mmdp_rates(const Matrix<T>& R) {
    std::vector<T> r(R.rows(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < R.rows(); ++i) r[i] = R(i, i);
    return r;
}

/** Stationary mean flow rate, `pi diag(R)`. */
template <class T>
T mmdp_mean_rate(const Matrix<T>& Q, const Matrix<T>& R) {
    const std::size_t n = Q.rows();
    if (n == 0) throw InputError("mmdp_mean_rate: empty generator");
    if (R.rows() != n || R.cols() != n)
        throw InputError("mmdp_mean_rate: R must be the same size as Q");
    // One state cannot switch, so its rate IS the mean and no solve is needed.
    if (n == 1) return R(0, 0);
    const std::vector<T> pi = mc::ctmc_solve(Q);
    T m = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) m += pi[i] * R(i, i);
    return m;
}

/**
 * SCV of the RATE under the stationary law: `Var[r] / E[r]^2`.
 *
 * A one-state process has no variability and returns zero. A mean rate of zero
 * leaves the ratio undefined, and the reference returns infinity rather than
 * dividing -- that is a signal the caller can test, where a NaN is not.
 */
template <class T>
T mmdp_scv(const Matrix<T>& Q, const Matrix<T>& R) {
    const std::size_t n = Q.rows();
    if (n == 0) throw InputError("mmdp_scv: empty generator");
    if (R.rows() != n || R.cols() != n)
        throw InputError("mmdp_scv: R must be the same size as Q");
    const T zero = num_traits<T>::from_int(0);
    if (n == 1) return zero;

    const std::vector<T> pi = mc::ctmc_solve(Q);
    T m = zero, m2 = zero;
    for (std::size_t i = 0; i < n; ++i) {
        m += pi[i] * R(i, i);
        m2 += pi[i] * R(i, i) * R(i, i);
    }
    if (!(num_traits<T>::to_double(m) > 0.0))
        return num_traits<T>::from_double(std::numeric_limits<double>::infinity());
    return T((m2 - m * m) / (m * m));
}

/** The (Q, R) pair of an MMDP. */
template <class T>
struct MmdpPair {
    Matrix<T> Q, R;
};

/**
 * The MMDP of a MAP: `Q = D0 + D1`, `R = diag(row sums of D1)`.
 *
 * The row sum of D1 is the total ARRIVAL rate out of a phase, so the fluid
 * analogue flows at exactly the rate at which that phase would be generating
 * arrivals. The MAP's phase process is unchanged, which is why Q is its full
 * generator.
 */
template <class T>
MmdpPair<T> mmdp_from_map(const Matrix<T>& D0, const Matrix<T>& D1) {
    const std::size_t n = D0.rows();
    if (n == 0 || D0.cols() != n || D1.rows() != n || D1.cols() != n)
        throw InputError("mmdp_from_map: D0 and D1 must be square and the same size");
    MmdpPair<T> out;
    out.Q = Matrix<T>(n, n, num_traits<T>::from_int(0));
    out.R = Matrix<T>(n, n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i) {
        T row = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < n; ++j) {
            out.Q(i, j) = D0(i, j) + D1(i, j);
            row += D1(i, j);
        }
        out.R(i, i) = row;
    }
    return out;
}

/**
 * The two-state MMDP, in the reference's own parameterization.
 *
 * `sigma0` is the rate out of state 0 and `sigma1` the rate out of state 1, so
 * the generator is `[[-s0, s0], [s1, -s1]]` and the stationary law is
 * `(s1, s0)/(s0+s1)` -- the chain spends LONGER in the state it leaves more
 * slowly, which is why the weights look swapped.
 */
template <class T>
MmdpPair<T> mmdp2(const T& r0, const T& r1, const T& sigma0, const T& sigma1) {
    const T zero = num_traits<T>::from_int(0);
    if (!(num_traits<T>::to_double(sigma0) > 0.0) || !(num_traits<T>::to_double(sigma1) > 0.0))
        throw InputError("mmdp2: both switching rates must be positive");
    MmdpPair<T> out;
    out.Q = Matrix<T>(2, 2, zero);
    out.R = Matrix<T>(2, 2, zero);
    out.Q(0, 0) = -sigma0;
    out.Q(0, 1) = sigma0;
    out.Q(1, 0) = sigma1;
    out.Q(1, 1) = -sigma1;
    out.R(0, 0) = r0;
    out.R(1, 1) = r1;
    return out;
}

/** The closed-form mean rate of a two-state MMDP. */
template <class T>
T mmdp2_mean_rate(const T& r0, const T& r1, const T& sigma0, const T& sigma1) {
    return T((r0 * sigma1 + r1 * sigma0) / (sigma0 + sigma1));
}

/** The closed-form SCV of a two-state MMDP. */
template <class T>
T mmdp2_scv(const T& r0, const T& r1, const T& sigma0, const T& sigma1) {
    const T tot = T(sigma0 + sigma1);
    const T p0 = T(sigma1 / tot), p1 = T(sigma0 / tot);
    const T m = T(p0 * r0 + p1 * r1);
    if (!(num_traits<T>::to_double(m) > 0.0))
        return num_traits<T>::from_double(std::numeric_limits<double>::infinity());
    const T v = T(p0 * r0 * r0 + p1 * r1 * r1 - m * m);
    return T(v / (m * m));
}

}  // namespace mmdp
}  // namespace line

#endif  // LINE_API_MMDP_MMDP_H
