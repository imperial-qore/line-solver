/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_DTIME_H
#define LINE_API_MAM_DTIME_H

/**
 * Discrete-time (slotted) matrix-analytic primitives and queues.
 *
 * Port of matlab/src/api/mam/dph_from_dist.m, dph_to_dmap.m, dmap_to_dph.m,
 * dmap_is_renewal.m, dmap_lambda.m, dmap_super.m, dmap_thin.m and
 * mg1_dt_queue.m, plus the queue-length half of the Q-MAM discrete-time queues
 * Q_DT_MAP_MAP_1.m and Q_DT_PH_PH_1.m.
 *
 * Everything measures time in SLOTS and follows the late arrival system with
 * delayed access (LAS-DA): within a slot the service completion resolves first,
 * arrivals are appended at the end of the slot and cannot enter service before
 * the next one, and the level is read after both. The QBD blocks state it
 * directly, A1 = kron(C1,D0) being the claim that a job arriving at the end of a
 * slot is not served within it. The LDES slotted engine orders its intra-slot
 * events the same way, so the two are directly comparable.
 *
 * A discrete phase-type law is (alpha, A) with P[X=k] = alpha A^(k-1) a,
 * a = e - A e, k = 1,2,... A batch stream is a vector [A_0, A_1, ...] whose
 * A_k carries the slots delivering k events; a plain D-MAP is the two-entry case.
 *
 * Two solver differences from the MATLAB twin, both deliberate and measured
 * against it: the QBD fundamental matrix uses logarithmic reduction rather than
 * SMCSolver's cyclic reduction (same minimal solution, both iterate to 1e-14),
 * and the batch M/G/1-type chain is solved by level truncation rather than by
 * MG1_CR, because C++ carries neither SMCSolver nor BuTools. The truncation
 * level is chosen from the tail mass, so it is an accuracy knob, not a model
 * change.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/dmap.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/lang/lang_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** A discrete phase-type law: initial row vector alpha and transient A. */
template <class T>
struct Dph {
    std::vector<T> alpha;
    Matrix<T> A;

    std::size_t order() const { return A.rows(); }
};

/** A discrete batch arrival stream, entry k carrying the slots with k events. */
template <class T>
using DBatch = std::vector<Matrix<T>>;

namespace detail {

/** Kronecker product; C++ has no shared templated kron. */
template <class T>
Matrix<T> dt_kron(const Matrix<T>& A, const Matrix<T>& B) {
    Matrix<T> C(A.rows() * B.rows(), A.cols() * B.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) {
            const T& a = A(i, j);
            for (std::size_t k = 0; k < B.rows(); ++k)
                for (std::size_t l = 0; l < B.cols(); ++l)
                    C(i * B.rows() + k, j * B.cols() + l) = a * B(k, l);
        }
    return C;
}

template <class T>
Matrix<T> dt_eye(std::size_t n) {
    Matrix<T> I(n, n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i) I(i, i) = num_traits<T>::from_int(1);
    return I;
}

template <class T>
Matrix<T> dt_add(const Matrix<T>& A, const Matrix<T>& B) {
    Matrix<T> C = A;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) = C(i, j) + B(i, j);
    return C;
}

template <class T>
Matrix<T> dt_sub(const Matrix<T>& A, const Matrix<T>& B) {
    Matrix<T> C = A;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) = C(i, j) - B(i, j);
    return C;
}

template <class T>
Matrix<T> dt_mul(const Matrix<T>& A, const Matrix<T>& B) {
    Matrix<T> C(A.rows(), B.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t k = 0; k < A.cols(); ++k) {
            const T& a = A(i, k);
            for (std::size_t j = 0; j < B.cols(); ++j) C(i, j) = C(i, j) + a * B(k, j);
        }
    return C;
}

/** Solves X * M = N for X, i.e. M^T X^T = N^T column by column. */
template <class T>
Matrix<T> dt_right_divide(const Matrix<T>& N, const Matrix<T>& M) {
    std::size_t n = M.rows();
    Matrix<T> Mt(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) Mt(i, j) = M(j, i);
    Matrix<T> X(N.rows(), n);
    for (std::size_t r = 0; r < N.rows(); ++r) {
        std::vector<T> rhs(n);
        for (std::size_t j = 0; j < n; ++j) rhs[j] = N(r, j);
        std::vector<T> sol = solve(Mt, rhs);
        for (std::size_t j = 0; j < n; ++j) X(r, j) = sol[j];
    }
    return X;
}

/** Left-multiplies the inverse: returns M^-1 * N. */
template <class T>
Matrix<T> dt_left_solve(const Matrix<T>& M, const Matrix<T>& N) {
    std::size_t n = M.rows();
    Matrix<T> X(n, N.cols());
    for (std::size_t c = 0; c < N.cols(); ++c) {
        std::vector<T> rhs(n);
        for (std::size_t i = 0; i < n; ++i) rhs[i] = N(i, c);
        std::vector<T> sol = solve(M, rhs);
        for (std::size_t i = 0; i < n; ++i) X(i, c) = sol[i];
    }
    return X;
}

}  // namespace detail

/**
 * Exact discrete phase-type representation of a lattice-valued law.
 *
 * Geometric, Det and DiscreteUniform are represented EXACTLY, not
 * moment-matched: a fitted surrogate would leave the lattice the caller relies
 * on, so any other family is an error here.
 */
template <class T>
Dph<T> dph_from_dist(lang::ProcessType type, const T& mean_slots, const T& scv) {
    const T one = num_traits<T>::from_int(1);
    const T zero = num_traits<T>::from_int(0);
    const double tol = 1e-8;

    if (type == lang::ProcessType::GEOMETRIC) {
        T p = one / mean_slots;
        if (num_traits<T>::to_double(p) > 1 + tol || num_traits<T>::to_double(p) <= 0)
            throw InputError("dph_from_dist: Geometric mean is outside the support {1,2,...}");
        if (num_traits<T>::to_double(p) > 1) p = one;
        Dph<T> d;
        d.alpha.assign(1, one);
        d.A = Matrix<T>(1, 1, one - p);
        return d;
    }

    if (type == lang::ProcessType::DET) {
        double m = num_traits<T>::to_double(mean_slots);
        long k = static_cast<long>(m + 0.5);
        if (std::abs(m - static_cast<double>(k)) > tol * std::max(1.0, m) || k < 1)
            throw InputError("dph_from_dist: Det is not a positive integral number of slots");
        std::size_t n = static_cast<std::size_t>(k);
        Dph<T> d;
        d.alpha.assign(n, zero);
        d.alpha[0] = one;
        d.A = Matrix<T>(n, n, zero);
        for (std::size_t i = 0; i + 1 < n; ++i) d.A(i, i + 1) = one;
        return d;
    }

    if (type == lang::ProcessType::DUNIFORM) {
        double m = num_traits<T>::to_double(mean_slots);
        double v = num_traits<T>::to_double(scv) * m * m;
        double width = std::sqrt(std::max(0.0, 12 * v + 1)) - 1;
        long lo = static_cast<long>(m - width / 2 + 0.5);
        long hi = static_cast<long>(m + width / 2 + 0.5);
        if (lo < 1 || hi < lo)
            throw InputError("dph_from_dist: DiscreteUniform is outside the support {1,2,...}");
        std::size_t n = static_cast<std::size_t>(hi);
        Dph<T> d;
        d.alpha.assign(n, zero);
        d.alpha[0] = one;
        d.A = Matrix<T>(n, n, zero);
        for (long j = 1; j < hi; ++j) {
            // hazard of absorbing at step j, zero below the lower bound
            T h = j < lo ? zero : one / num_traits<T>::from_int(static_cast<int>(hi - j + 1));
            d.A(static_cast<std::size_t>(j - 1), static_cast<std::size_t>(j)) = one - h;
        }
        return d;
    }

    throw InputError(
        "dph_from_dist: the process type has no exact discrete phase-type representation. "
        "The discrete-time path accepts Geometric, Det on the slot lattice, DiscreteUniform "
        "and DMAP.");
}

/** Renewal D-MAP (A, a alpha) of a discrete phase-type law. */
template <class T>
Dmap<T> dph_to_dmap(const Dph<T>& d) {
    const T one = num_traits<T>::from_int(1);
    std::size_t m = d.A.rows();
    Dmap<T> out;
    out.D0 = d.A;
    out.D1 = Matrix<T>(m, m, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < m; ++i) {
        T rowsum = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < m; ++j) rowsum = rowsum + d.A(i, j);
        T a = one - rowsum;
        for (std::size_t j = 0; j < m; ++j) out.D1(i, j) = a * d.alpha[j];
    }
    return out;
}

/** True when D1 has rank one, i.e. the process renews at every event. */
template <class T>
bool dmap_is_renewal(const Dmap<T>& d) {
    std::size_t m = d.D0.rows();
    if (m == 1) return true;
    std::vector<T> row_mass(m, num_traits<T>::from_int(0));
    std::size_t pivot = 0;
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) row_mass[i] = row_mass[i] + d.D1(i, j);
        if (num_traits<T>::to_double(row_mass[i]) > num_traits<T>::to_double(row_mass[pivot]))
            pivot = i;
    }
    if (num_traits<T>::to_double(row_mass[pivot]) <= 1e-14) return false;
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) {
            T expected = row_mass[i] * d.D1(pivot, j) / row_mass[pivot];
            if (std::abs(num_traits<T>::to_double(d.D1(i, j) - expected)) > 1e-8) return false;
        }
    return true;
}

/** Discrete phase-type law underlying a renewal D-MAP. */
template <class T>
Dph<T> dmap_to_dph(const Dmap<T>& d) {
    if (!dmap_is_renewal(d))
        throw InputError("dmap_to_dph: the D-MAP does not renew at events, so it has no DPH form");
    std::size_t m = d.D0.rows();
    std::vector<T> row_mass(m, num_traits<T>::from_int(0));
    std::size_t pivot = 0;
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) row_mass[i] = row_mass[i] + d.D1(i, j);
        if (num_traits<T>::to_double(row_mass[i]) > num_traits<T>::to_double(row_mass[pivot]))
            pivot = i;
    }
    Dph<T> out;
    out.A = d.D0;
    out.alpha.resize(m);
    for (std::size_t j = 0; j < m; ++j) out.alpha[j] = d.D1(pivot, j) / row_mass[pivot];
    return out;
}

/**
 * Mean number of EVENTS per slot, pi sum_k k A_k e. A slot carrying a batch of
 * two counts twice, which is what Little's law consumes downstream.
 */
template <class T>
T dmap_lambda_batch(const DBatch<T>& A) {
    std::size_t m = A[0].rows();
    Matrix<T> P(m, m, num_traits<T>::from_int(0));
    for (std::size_t k = 0; k < A.size(); ++k) P = detail::dt_add(P, A[k]);
    std::vector<T> pi = mc::dtmc_solve(P);
    T lambda = num_traits<T>::from_int(0);
    for (std::size_t k = 1; k < A.size(); ++k)
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j)
                lambda = lambda + pi[i] * num_traits<T>::from_int(static_cast<int>(k)) * A[k](i, j);
    return lambda;
}

/**
 * Superposition, E_k = sum_{i+j=k} kron(A_i, B_j).
 *
 * NOT closed on D-MAPs: two slotted streams fire in the same slot with positive
 * probability, so the merged stream carries batches. Folding E_2 into E_1 would
 * conserve neither the arrival rate nor the slot in which the work appears, so
 * the batch dimension is kept and the station is solved as an M/G/1-type chain
 * instead of a QBD.
 */
template <class T>
DBatch<T> dmap_super(const DBatch<T>& A, const DBatch<T>& B) {
    std::size_t p = A.size() - 1, q = B.size() - 1;
    DBatch<T> E;
    E.reserve(p + q + 1);
    for (std::size_t k = 0; k <= p + q; ++k) {
        Matrix<T> Ek(A[0].rows() * B[0].rows(), A[0].cols() * B[0].cols(),
                     num_traits<T>::from_int(0));
        std::size_t lo = k > q ? k - q : 0;
        for (std::size_t i = lo; i <= (k < p ? k : p); ++i)
            Ek = detail::dt_add(Ek, detail::dt_kron(A[i], B[k - i]));
        E.push_back(Ek);
    }
    return E;
}

/**
 * Bernoulli thinning, B_k = sum_{n>=k} C(n,k) p^k (1-p)^(n-k) A_n. The phase
 * process is untouched, so this is exact for PROB/RAND routing.
 */
template <class T>
DBatch<T> dmap_thin(const DBatch<T>& A, const T& p) {
    double pd = num_traits<T>::to_double(p);
    if (pd < 0 || pd > 1) throw InputError("dmap_thin: the routing probability must lie in [0,1]");
    const T one = num_traits<T>::from_int(1);
    std::size_t n = A.size() - 1;
    DBatch<T> B;
    B.reserve(n + 1);
    for (std::size_t k = 0; k <= n; ++k) {
        Matrix<T> Bk(A[0].rows(), A[0].cols(), num_traits<T>::from_int(0));
        for (std::size_t j = k; j <= n; ++j) {
            T w = num_traits<T>::from_int(1);
            for (std::size_t i = 0; i < k; ++i)
                w = w * num_traits<T>::from_int(static_cast<int>(j - i)) /
                    num_traits<T>::from_int(static_cast<int>(i + 1));
            for (std::size_t i = 0; i < k; ++i) w = w * p;
            for (std::size_t i = 0; i < j - k; ++i) w = w * (one - p);
            if (num_traits<T>::to_double(w) > 0)
                for (std::size_t r = 0; r < Bk.rows(); ++r)
                    for (std::size_t c = 0; c < Bk.cols(); ++c)
                        Bk(r, c) = Bk(r, c) + w * A[j](r, c);
        }
        B.push_back(Bk);
    }
    // trailing zero batch levels carry no mass and only inflate the blocks
    while (B.size() > 2) {
        double worst = 0;
        for (std::size_t r = 0; r < B.back().rows(); ++r)
            for (std::size_t c = 0; c < B.back().cols(); ++c)
                worst = std::max(worst, std::abs(num_traits<T>::to_double(B.back()(r, c))));
        if (worst >= 1e-14) break;
        B.pop_back();
    }
    return B;
}

namespace detail {

/**
 * Raw moments over {1,2,...} of a two-phase acyclic discrete phase-type law
 * written as a mixture: with probability b the sum of two geometrics with
 * success probabilities p1 and p2, otherwise the second geometric alone.
 */
inline void dph2_moments(double b, double p1, double p2, double* m1, double* m2, double* m3) {
    const double A1 = 1 / p1, A2 = 1 / p2;
    const double S1 = (2 - p1) / (p1 * p1), S2 = (2 - p2) / (p2 * p2);
    const double C1 = (p1 * p1 - 6 * p1 + 6) / (p1 * p1 * p1);
    const double C2 = (p2 * p2 - 6 * p2 + 6) / (p2 * p2 * p2);
    *m1 = b * A1 + A2;
    *m2 = b * (S1 + 2 * A1 * A2) + S2;
    *m3 = b * (C1 + 3 * S1 * A2 + 3 * A1 * S2) + C2;
}

/**
 * Fits the mixture above to three raw moments by a damped Newton iteration on
 * (p1,p2), the mixing weight following from the first moment. Returns false
 * when no feasible triple is reached from any start.
 */
inline bool dph2_from_3moments(double m1, double m2, double m3, double* b_out, double* p1_out,
                               double* p2_out) {
    const double starts[4][2] = {{0.5, 0.5}, {0.9, 1 / std::max(1.0, m1)}, {0.2, 0.8}, {0.99, 0.3}};
    for (int s = 0; s < 4; ++s) {
        double p1 = starts[s][0], p2 = starts[s][1];
        for (int it = 0; it < 200; ++it) {
            const double b = (m1 - 1 / p2) * p1;
            double f1, f2, f3;
            dph2_moments(b, p1, p2, &f1, &f2, &f3);
            const double r2 = f2 - m2, r3 = f3 - m3;
            if (std::abs(r2) <= 1e-11 * std::max(1.0, std::abs(m2)) &&
                std::abs(r3) <= 1e-11 * std::max(1.0, std::abs(m3))) {
                if (b < -1e-9 || b > 1 + 1e-9 || p1 <= 0 || p1 > 1 || p2 <= 0 || p2 > 1) break;
                *b_out = std::min(1.0, std::max(0.0, b));
                *p1_out = p1;
                *p2_out = p2;
                return true;
            }
            const double h = 1e-7;
            double a2, a3, c2, c3, d2, d3, dummy;
            dph2_moments((m1 - 1 / p2) * (p1 + h), p1 + h, p2, &dummy, &a2, &a3);
            dph2_moments((m1 - 1 / (p2 + h)) * p1, p1, p2 + h, &dummy, &c2, &c3);
            const double j11 = (a2 - f2) / h, j12 = (c2 - f2) / h;
            const double j21 = (a3 - f3) / h, j22 = (c3 - f3) / h;
            const double det = j11 * j22 - j12 * j21;
            if (std::abs(det) < 1e-18) break;
            double dp1 = -(j22 * r2 - j12 * r3) / det;
            double dp2 = -(-j21 * r2 + j11 * r3) / det;
            // damping keeps the step inside (0,1], where the parameters live
            double step = 1.0;
            while (step > 1e-6 && (p1 + step * dp1 <= 1e-9 || p1 + step * dp1 > 1 ||
                                   p2 + step * dp2 <= 1e-9 || p2 + step * dp2 > 1))
                step /= 2;
            if (step <= 1e-6) break;
            p1 += step * dp1;
            p2 += step * dp2;
        }
    }
    return false;
}

}  // namespace detail

/**
 * Reduces the order of a D-MAP by matching interevent moments.
 *
 * Leaves the process untouched while its order stays within `max_order`, and
 * otherwise replaces it by the two-phase discrete phase-type law with the same
 * first three interevent moments, read back as a renewal D-MAP. Correlation is
 * NOT preserved, which is the reason the multi-station discrete-time path is an
 * approximation. When the moment triple is outside the DPH(2) region the
 * fallback keeps the exact mean with a Geometric, so the event rate of the
 * decomposition is conserved in every branch.
 *
 * The MATLAB twin reaches the same three moments through BuTools
 * (MGFromMoments then CanonicalFromDPH2), which admits the wider
 * matrix-geometric class; this port solves the acyclic DPH(2) directly, so a
 * triple that is MG(2)-feasible but not DPH(2)-feasible takes the Geometric
 * fallback here and the two-phase fit there.
 */
template <class T>
Dmap<T> dmap_compress(const Dmap<T>& d, std::size_t max_order) {
    if (d.D0.rows() <= max_order) return d;

    std::vector<T> moms = dmap_moment(d, std::vector<unsigned>{1u, 2u, 3u});
    const double m1 = num_traits<T>::to_double(moms[0]);
    const double m2 = num_traits<T>::to_double(moms[1]);
    const double m3 = num_traits<T>::to_double(moms[2]);

    double b = 0, p1 = 0, p2 = 0;
    if (detail::dph2_from_3moments(m1, m2, m3, &b, &p1, &p2)) {
        Dph<T> fit;
        fit.alpha.assign(2, num_traits<T>::from_int(0));
        fit.alpha[0] = num_traits<T>::from_double(b);
        fit.alpha[1] = num_traits<T>::from_double(1 - b);
        fit.A = Matrix<T>(2, 2, num_traits<T>::from_int(0));
        fit.A(0, 0) = num_traits<T>::from_double(1 - p1);
        fit.A(0, 1) = num_traits<T>::from_double(p1);
        fit.A(1, 1) = num_traits<T>::from_double(1 - p2);
        Dmap<T> cand = dph_to_dmap(fit);
        if (dmap_isfeasible(cand)) return cand;
    }

    // the moment triple is not DPH(2)-feasible: keep the rate, drop the shape
    double p = 1 / m1;
    p = std::min(1.0, std::max(1e-12, p));
    Dph<T> geo;
    geo.alpha.assign(1, num_traits<T>::from_int(1));
    geo.A = Matrix<T>(1, 1, num_traits<T>::from_double(1 - p));
    return dph_to_dmap(geo);
}

/**
 * Order reduction of a BATCH stream.
 *
 * The process of NONEMPTY SLOTS is compressed as a plain D-MAP and the
 * stationary batch-size distribution, conditional on the slot being nonempty,
 * is reattached to it. Compressing the phase alone would keep the slot process
 * and lose the batch sizes, which would silently rescale the event rate.
 *
 * MATLAB twin: dmap_compress_batch.m
 */
template <class T>
DBatch<T> dmap_compress_batch(const DBatch<T>& A, std::size_t max_order) {
    if (A[0].rows() <= max_order) return A;
    const std::size_t nb = A.size() - 1, m = A[0].rows();

    Matrix<T> Ptot = A[0];
    for (std::size_t k = 1; k < A.size(); ++k) Ptot = detail::dt_add(Ptot, A[k]);
    std::vector<T> pi_phase = mc::dtmc_solve(Ptot);

    std::vector<T> qraw(nb, num_traits<T>::from_int(0));
    T mass = num_traits<T>::from_int(0);
    for (std::size_t k = 0; k < nb; ++k) {
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) qraw[k] = qraw[k] + pi_phase[i] * A[k + 1](i, j);
        mass = mass + qraw[k];
    }
    if (num_traits<T>::to_double(mass) <= 0)
        throw InputError("dmap_compress_batch: the batch stream carries no events");

    Dmap<T> marked;
    marked.D0 = A[0];
    marked.D1 = detail::dt_sub(Ptot, A[0]);
    Dmap<T> markedC = dmap_compress(marked, max_order);

    DBatch<T> out;
    out.push_back(markedC.D0);
    for (std::size_t k = 0; k < nb; ++k) {
        Matrix<T> blk = markedC.D1;
        const T w = qraw[k] / mass;
        for (std::size_t i = 0; i < blk.rows(); ++i)
            for (std::size_t j = 0; j < blk.cols(); ++j) blk(i, j) = blk(i, j) * w;
        out.push_back(blk);
    }
    return out;
}

/**
 * G matrix of a discrete-time QBD by logarithmic reduction (Latouche and
 * Ramaswami). The blocks are stochastic, so no uniformization is needed; the
 * MATLAB twin reaches the same minimal solution through SMCSolver's cyclic
 * reduction.
 */
template <class T>
Matrix<T> qbd_dt_g(const Matrix<T>& A0, const Matrix<T>& A1, const Matrix<T>& A2,
                   int max_iter = 200) {
    std::size_t m = A1.rows();
    Matrix<T> I = detail::dt_eye<T>(m);
    Matrix<T> inv_local = detail::dt_left_solve(detail::dt_sub(I, A1), I);
    Matrix<T> B0 = detail::dt_mul(inv_local, A0);
    Matrix<T> B2 = detail::dt_mul(inv_local, A2);
    Matrix<T> G = B0;
    Matrix<T> PI = B2;

    for (int it = 0; it < max_iter; ++it) {
        Matrix<T> A1n = detail::dt_add(detail::dt_mul(B0, B2), detail::dt_mul(B2, B0));
        Matrix<T> A0n = detail::dt_mul(B0, B0);
        Matrix<T> A2n = detail::dt_mul(B2, B2);
        Matrix<T> inv_n = detail::dt_left_solve(detail::dt_sub(I, A1n), I);
        B0 = detail::dt_mul(inv_n, A0n);
        B2 = detail::dt_mul(inv_n, A2n);
        G = detail::dt_add(G, detail::dt_mul(PI, B0));
        PI = detail::dt_mul(PI, B2);

        double residual = 0;
        for (std::size_t i = 0; i < m; ++i) {
            T rowsum = num_traits<T>::from_int(0);
            for (std::size_t j = 0; j < m; ++j) rowsum = rowsum + G(i, j);
            residual = std::max(residual,
                                std::abs(1.0 - num_traits<T>::to_double(rowsum)));
        }
        if (residual < 1e-14) break;
    }
    return G;
}

/**
 * G matrix of an M/G/1-type chain by functional iteration on
 * G = sum_k A_k G^k, the blocks being stochastic.
 *
 * `A[0]` is the down-one block and `A[k]` the block raising the level by k-1.
 * SMCSolver reaches the same minimal solution by cyclic reduction; the JAR port
 * of that routine costs 19.8 s at order 17 where functional iteration costs
 * 25 ms and agrees to 1.6e-15, which is why neither this port nor the JAR one
 * takes it.
 */
template <class T>
Matrix<T> mg1_dt_g(const std::vector<Matrix<T>>& A, int max_iter = 5000,
                   double tol = 1e-14) {
    std::size_t m = A[0].rows();
    Matrix<T> G = A[0];
    for (int it = 0; it < max_iter; ++it) {
        Matrix<T> Gnew = A[0];
        Matrix<T> Gpow = G;
        for (std::size_t k = 1; k < A.size(); ++k) {
            Gnew = detail::dt_add(Gnew, detail::dt_mul(A[k], Gpow));
            if (k + 1 < A.size()) Gpow = detail::dt_mul(Gpow, G);
        }
        double delta = 0;
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j)
                delta = std::max(delta, std::abs(num_traits<T>::to_double(Gnew(i, j) - G(i, j))));
        G = Gnew;
        if (delta < tol) break;
    }
    return G;
}

/**
 * Stationary vector of an M/G/1-type chain by the stable Ramaswami formula.
 *
 * `A` repeats from level one and `B` is the boundary row, both as block lists
 * with `A[0]` the down-one block. The recursion is level-by-level, so it costs
 * O(levels * m^3) and not O((levels*m)^3): a dense solve of the truncated chain
 * is cubic in the WHOLE state space and stops being affordable at the second
 * station of a decomposition, where the arrival process already carries the
 * level space of the first.
 *
 * MATLAB twin: MG1_pi.m (SMCSolver, Van Houdt), default-boundary branch.
 */
template <class T>
std::vector<T> mg1_dt_pi(const std::vector<Matrix<T>>& B, const std::vector<Matrix<T>>& A,
                         std::size_t max_num_comp = 1000) {
    const T one = num_traits<T>::from_int(1);
    const T zero = num_traits<T>::from_int(0);
    std::size_t m = A[0].rows();
    std::size_t dega = A.size() - 1;
    std::size_t degb = B.size() - 1;
    Matrix<T> I = detail::dt_eye<T>(m);
    Matrix<T> G = mg1_dt_g(A);

    // hatA_i = sum_{v>=i} A_v G^(v-i), while sumA accumulates the originals
    std::vector<Matrix<T>> hatA = A;
    Matrix<T> sumA = A[dega];
    std::vector<T> beta(m, zero);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) beta[i] = beta[i] + sumA(i, j);
    for (std::size_t i = dega; i-- > 1;) {
        sumA = detail::dt_add(sumA, A[i]);
        hatA[i] = detail::dt_add(A[i], detail::dt_mul(hatA[i + 1], G));
        for (std::size_t r = 0; r < m; ++r)
            for (std::size_t c = 0; c < m; ++c) beta[r] = beta[r] + sumA(r, c);
    }
    sumA = detail::dt_add(sumA, A[0]);

    std::vector<T> theta = mc::dtmc_solve(sumA);
    T drift = zero;
    for (std::size_t i = 0; i < m; ++i) drift = drift + theta[i] * beta[i];
    if (num_traits<T>::to_double(drift) >= 1)
        throw InputError("mg1_dt_pi: the chain characterized by A is not positive recurrent");

    Matrix<T> invBarA1 = detail::dt_left_solve(detail::dt_sub(I, hatA[1]), I);

    // hatB_i = sum_{v>=i} B_v G^(v-i), sumBB0 = sum_{v>=1} B_v,
    // Bbeta = sum_{v>=1} (v-1) B_v e
    std::vector<Matrix<T>> hatB = B;
    Matrix<T> sumBB0 = B[degb];
    std::vector<T> Bbeta(m, zero);
    for (std::size_t i = degb; i-- > 1;) {
        for (std::size_t r = 0; r < m; ++r)
            for (std::size_t c = 0; c < m; ++c) Bbeta[r] = Bbeta[r] + sumBB0(r, c);
        sumBB0 = detail::dt_add(sumBB0, B[i]);
        hatB[i] = detail::dt_add(B[i], detail::dt_mul(hatB[i + 1], G));
    }

    Matrix<T> Kmat = detail::dt_add(B[0], detail::dt_mul(hatB[1], G));
    std::vector<T> kappa = mc::dtmc_solve(Kmat);
    std::vector<T> g = mc::dtmc_solve(G);

    // temp = rowsum(inv(I - sumA - (e - beta) g))
    Matrix<T> W = detail::dt_sub(I, sumA);
    for (std::size_t r = 0; r < m; ++r)
        for (std::size_t c = 0; c < m; ++c) W(r, c) = W(r, c) - (one - beta[r]) * g[c];
    Matrix<T> invW = detail::dt_left_solve(W, I);
    std::vector<T> temp(m, zero);
    for (std::size_t r = 0; r < m; ++r)
        for (std::size_t c = 0; c < m; ++c) temp[r] = temp[r] + invW(r, c);

    const T inv_slack = one / (one - drift);
    std::vector<T> psi1(m, zero), psi2(m, one);
    Matrix<T> M1 = detail::dt_sub(detail::dt_sub(I, A[0]), hatA[1]);
    Matrix<T> M2 = detail::dt_sub(sumBB0, hatB[1]);
    for (std::size_t r = 0; r < m; ++r) {
        T acc1 = zero, acc2 = zero, a0row = zero;
        for (std::size_t c = 0; c < m; ++c) {
            acc1 = acc1 + M1(r, c) * temp[c];
            acc2 = acc2 + M2(r, c) * temp[c];
            a0row = a0row + A[0](r, c);
        }
        psi1[r] = acc1 + inv_slack * a0row;
        psi2[r] = one + acc2 + inv_slack * Bbeta[r];
    }
    // tildekappa1 = psi2 + hatB_1 inv(I - hatA_1) psi1
    std::vector<T> tmp(m, zero), tilde(m, zero);
    for (std::size_t r = 0; r < m; ++r)
        for (std::size_t c = 0; c < m; ++c) tmp[r] = tmp[r] + invBarA1(r, c) * psi1[c];
    for (std::size_t r = 0; r < m; ++r) {
        T acc = zero;
        for (std::size_t c = 0; c < m; ++c) acc = acc + hatB[1](r, c) * tmp[c];
        tilde[r] = psi2[r] + acc;
    }
    T denom = zero;
    for (std::size_t r = 0; r < m; ++r) denom = denom + kappa[r] * tilde[r];

    std::vector<std::vector<T>> pi;
    std::vector<T> pi0(m, zero);
    for (std::size_t r = 0; r < m; ++r) pi0[r] = kappa[r] / denom;
    pi.push_back(pi0);

    double sumpi = 0;
    for (std::size_t r = 0; r < m; ++r) sumpi += num_traits<T>::to_double(pi0[r]);
    std::size_t numit = 1;
    while (sumpi < 1 - 1e-10 && numit < max_num_comp) {
        std::vector<T> pin(m, zero);
        if (numit <= degb)
            for (std::size_t r = 0; r < m; ++r)
                for (std::size_t c = 0; c < m; ++c)
                    pin[c] = pin[c] + pi0[r] * hatB[numit](r, c);
        for (std::size_t j = 1; j <= std::min(numit - 1, dega - 1); ++j)
            for (std::size_t r = 0; r < m; ++r)
                for (std::size_t c = 0; c < m; ++c)
                    pin[c] = pin[c] + pi[numit - j][r] * hatA[j + 1](r, c);
        std::vector<T> row(m, zero);
        for (std::size_t r = 0; r < m; ++r)
            for (std::size_t c = 0; c < m; ++c) row[c] = row[c] + pin[r] * invBarA1(r, c);
        pi.push_back(row);
        for (std::size_t r = 0; r < m; ++r) sumpi += num_traits<T>::to_double(row[r]);
        ++numit;
    }

    std::vector<T> out;
    out.reserve(pi.size() * m);
    for (std::size_t i = 0; i < pi.size(); ++i)
        for (std::size_t r = 0; r < m; ++r) out.push_back(pi[i][r]);
    return out;
}

/**
 * Queue length distribution of a discrete-time D-MAP/D-MAP/1/FCFS queue, the
 * queue-length half of Q_DT_MAP_MAP_1.
 *
 * Entry i is Prob[i customers in system] under LAS-DA. The waiting and sojourn
 * pmfs of the Q-MAM routine are deliberately not ported: no LINE caller consumes
 * them, and the discrete-time solver path reads the queue length alone.
 */
template <class T>
std::vector<T> q_dt_map_map_1(const Dmap<T>& arv, const Dmap<T>& svc,
                              std::size_t max_num_comp = 1000) {
    const T one = num_traits<T>::from_int(1);
    std::size_t ma = arv.D0.rows(), ms = svc.D0.rows(), mtot = ma * ms;

    std::vector<T> pi_a = mc::dtmc_solve(detail::dt_add(arv.D0, arv.D1));
    T avga = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < ma; ++i)
        for (std::size_t j = 0; j < ma; ++j) avga = avga + pi_a[i] * arv.D1(i, j);
    std::vector<T> pi_s = mc::dtmc_solve(detail::dt_add(svc.D0, svc.D1));
    T avgs = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < ms; ++i)
        for (std::size_t j = 0; j < ms; ++j) avgs = avgs + pi_s[i] * svc.D1(i, j);
    if (num_traits<T>::to_double(avga / avgs) >= 1)
        throw InputError("q_dt_map_map_1: the load of the system exceeds one");

    Matrix<T> Ims = detail::dt_eye<T>(ms);
    Matrix<T> Am1 = detail::dt_kron(arv.D0, svc.D1);
    Matrix<T> A0 = detail::dt_add(detail::dt_kron(arv.D0, svc.D0),
                                  detail::dt_kron(arv.D1, svc.D1));
    Matrix<T> A1 = detail::dt_kron(arv.D1, svc.D0);
    Matrix<T> B0 = detail::dt_kron(arv.D0, Ims);
    Matrix<T> B1 = detail::dt_kron(arv.D1, Ims);

    Matrix<T> G = qbd_dt_g(Am1, A0, A1);
    // R = A1 (I - A0 - A1 G)^-1
    Matrix<T> I = detail::dt_eye<T>(mtot);
    Matrix<T> denom = detail::dt_sub(detail::dt_sub(I, A0), detail::dt_mul(A1, G));
    Matrix<T> R = detail::dt_right_divide(A1, denom);

    // General boundary [B1; A0 + R Am1]: the empty system has its own local
    // block, so levels 0 and 1 are solved together (QBD_pi.m else-branch)
    Matrix<T> lower = detail::dt_add(A0, detail::dt_mul(R, Am1));
    Matrix<T> joint(2 * mtot, 2 * mtot, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < mtot; ++i)
        for (std::size_t j = 0; j < mtot; ++j) {
            joint(i, j) = B0(i, j);
            joint(mtot + i, j) = Am1(i, j);
            joint(i, mtot + j) = B1(i, j);
            joint(mtot + i, mtot + j) = lower(i, j);
        }
    std::vector<T> pi01 = mc::dtmc_solve(joint);

    Matrix<T> temp = detail::dt_left_solve(detail::dt_sub(I, R), I);
    std::vector<T> pi0(pi01.begin(), pi01.begin() + mtot);
    std::vector<T> pi1(pi01.begin() + mtot, pi01.end());
    T norm = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < mtot; ++i) norm = norm + pi0[i];
    for (std::size_t i = 0; i < mtot; ++i)
        for (std::size_t j = 0; j < mtot; ++j) norm = norm + pi1[i] * temp(i, j);
    for (std::size_t i = 0; i < mtot; ++i) {
        pi0[i] = pi0[i] / norm;
        pi1[i] = pi1[i] / norm;
    }

    std::vector<T> ql;
    T mass0 = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < mtot; ++i) mass0 = mass0 + pi0[i];
    ql.push_back(mass0);
    std::vector<T> cur = pi1;
    double acc = num_traits<T>::to_double(mass0);
    for (std::size_t it = 0; it < max_num_comp; ++it) {
        T mass = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < mtot; ++i) mass = mass + cur[i];
        ql.push_back(mass);
        acc += num_traits<T>::to_double(mass);
        if (acc > 1 - 1e-10) break;
        std::vector<T> nxt(mtot, num_traits<T>::from_int(0));
        for (std::size_t j = 0; j < mtot; ++j)
            for (std::size_t i = 0; i < mtot; ++i) nxt[j] = nxt[j] + cur[i] * R(i, j);
        cur = nxt;
    }

    T total = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < ql.size(); ++i) total = total + ql[i];
    for (std::size_t i = 0; i < ql.size(); ++i) ql[i] = ql[i] / total;
    return ql;
}

/** Queue length of a discrete-time DPH/DPH/1/FCFS queue, via the D-MAP route. */
template <class T>
std::vector<T> q_dt_ph_ph_1(const Dph<T>& arv, const Dph<T>& svc,
                            std::size_t max_num_comp = 1000) {
    return q_dt_map_map_1(dph_to_dmap(arv), dph_to_dmap(svc), max_num_comp);
}

/** Outcome of a slotted station solve. */
template <class T>
struct DtQueueResult {
    T QN;
    T UN;
    T TN;
    std::vector<T> ql;
    Dmap<T> dep;
};

/**
 * Discrete-time single-server queue with batch D-MAP arrivals, DBMAP/D-MAP/1.
 *
 * The chain is M/G/1-type because a slot may deliver a batch: with arrival
 * matrices A_k and service pair (S0,S1),
 * A^(-1) = kron(A_0,S1), A^(k) = kron(A_k,S0) + kron(A_{k+1},S1),
 * B^(k) = kron(A_k,I), the boundary row holding the empty system where no
 * service runs.
 *
 * Solved by the Ramaswami level recursion of `mg1_dt_pi`. The DEPARTURE process
 * is a different object: it is read off the chain truncated at the level where
 * the tail carries less than 1e-10, with arrivals that would cross the top held
 * there, which is a level cut and not a rate change.
 */
template <class T>
DtQueueResult<T> mg1_dt_queue(const DBatch<T>& arv, const Dmap<T>& svc,
                              std::size_t max_num_comp = 1000, bool want_departure = false) {
    const T one = num_traits<T>::from_int(1);
    std::size_t ms = svc.D0.rows(), ma = arv[0].rows(), K = arv.size() - 1;
    std::size_t m = ma * ms;
    Matrix<T> Ims = detail::dt_eye<T>(ms);

    T lambda = dmap_lambda_batch(arv);
    DBatch<T> svc_batch;
    svc_batch.push_back(svc.D0);
    svc_batch.push_back(svc.D1);
    T mu = dmap_lambda_batch(svc_batch);
    if (num_traits<T>::to_double(lambda / mu) >= 1)
        throw InputError("mg1_dt_queue: the discrete-time load of the station is not below one");

    // A^(-1) = kron(A_0,S1), A^(k) = kron(A_k,S0) + kron(A_{k+1},S1),
    // B^(k) = kron(A_k,I), the boundary row where no service runs
    std::vector<Matrix<T>> Ablocks, Bblocks;
    Ablocks.push_back(detail::dt_kron(arv[0], svc.D1));
    for (std::size_t k = 0; k <= K; ++k) {
        Matrix<T> blk = detail::dt_kron(arv[k], svc.D0);
        if (k + 1 <= K) blk = detail::dt_add(blk, detail::dt_kron(arv[k + 1], svc.D1));
        Ablocks.push_back(blk);
    }
    for (std::size_t k = 0; k <= K; ++k) Bblocks.push_back(detail::dt_kron(arv[k], Ims));

    std::vector<T> pi = mg1_dt_pi(Bblocks, Ablocks, max_num_comp);

    DtQueueResult<T> out;
    std::size_t nlev = pi.size() / m;
    out.ql.assign(nlev, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < nlev; ++i)
        for (std::size_t j = 0; j < m; ++j) out.ql[i] = out.ql[i] + pi[i * m + j];
    T total = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < nlev; ++i) total = total + out.ql[i];
    for (std::size_t i = 0; i < nlev; ++i) out.ql[i] = out.ql[i] / total;

    out.QN = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < nlev; ++i)
        out.QN = out.QN + num_traits<T>::from_int(static_cast<int>(i)) * out.ql[i];
    out.UN = one - out.ql[0];
    out.TN = lambda;

    if (want_departure) {
        double cum = 0;
        std::size_t L = nlev - 1;
        for (std::size_t i = 0; i < nlev; ++i) {
            cum += num_traits<T>::to_double(out.ql[i]);
            if (cum > 1 - 1e-10) {
                L = i;
                break;
            }
        }
        if (L < 1) L = 1;
        std::size_t nstates = (L + 1) * m;
        // departure process of the same truncated chain: D1 collects the
        // transitions carrying a completion, D0 the rest
        Matrix<T> D0(nstates, nstates, num_traits<T>::from_int(0));
        Matrix<T> D1(nstates, nstates, num_traits<T>::from_int(0));
        for (std::size_t i = 0; i <= L; ++i)
            for (std::size_t k = 0; k <= K; ++k) {
                if (i == 0) {
                    std::size_t tgt = std::min(L, k);
                    Matrix<T> blk = detail::dt_kron(arv[k], Ims);
                    for (std::size_t r = 0; r < m; ++r)
                        for (std::size_t c = 0; c < m; ++c)
                            D0(i * m + r, tgt * m + c) = D0(i * m + r, tgt * m + c) + blk(r, c);
                } else {
                    std::size_t tno = std::min(L, i + k);
                    Matrix<T> no = detail::dt_kron(arv[k], svc.D0);
                    for (std::size_t r = 0; r < m; ++r)
                        for (std::size_t c = 0; c < m; ++c)
                            D0(i * m + r, tno * m + c) = D0(i * m + r, tno * m + c) + no(r, c);
                    std::size_t tdep = std::min(L, i - 1 + k);
                    Matrix<T> dp = detail::dt_kron(arv[k], svc.D1);
                    for (std::size_t r = 0; r < m; ++r)
                        for (std::size_t c = 0; c < m; ++c)
                            D1(i * m + r, tdep * m + c) = D1(i * m + r, tdep * m + c) + dp(r, c);
                }
            }
        out.dep.D0 = D0;
        out.dep.D1 = D1;
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_DTIME_H
