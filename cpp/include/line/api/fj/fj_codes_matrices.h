/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_FJ_CODES_MATRICES_H
#define LINE_API_FJ_FJ_CODES_MATRICES_H

/**
 * The state-space construction of FJ_codes, the fork-join response-time-tail
 * approximation of Z. Qiu, J. F. Perez and P. Harrison, "Beyond the Mean in
 * Fork-Join Queues: Efficient Approximation for Response-Time Tails" (IFIP
 * Performance 2015). Third-party, BSD-3-Clause, Copyright 2015 Imperial College
 * London; see THIRD-PARTY-NOTICES.md and
 * `python/line_solver/lib/thirdparty/fj/LICENSE.txt`.
 *
 * Port of `build_index.m`, `vectmatch.m`, `build_Service_h.m`, `build_SA.m`,
 * `generateService.m`, `constructSRK.m` and `constructNotAllBusy.m` from
 * `matlab/lib/thirdparty/FJ_codes`. The solve layer is `fj_codes.h`.
 *
 * THE STATE SPACE. The algorithm analyses the TWO-node fork-join queue exactly
 * (Section 4 of the paper) and interpolates to K nodes (Section 6). Its phase
 * process tracks, for the job at the head of the two branches, the service
 * phase of the subtask in the LONGER queue and of the subtask in the SHORTER
 * one, together with c, the difference in queue length between the two
 * branches, truncated at C. A phase is therefore a pair of unit count vectors
 * `[e_long, e_short]` of length m each, and the level is c in 0..C; `S` is the
 * generator of phase changes at constant c, `A_jump` collects the transitions
 * on which the head job leaves and the next one enters.
 *
 * WHY COUNT VECTORS AND NOT PHASE INDICES. The reference indexes every phase by
 * a length-2m vector of counts summing to two (one per branch) and finds the
 * target of a transition with a linear search, `vectmatch`. For the two-node
 * queue this is a bijection with the pair (i, j), and a direct index would be
 * faster. It is reproduced because the transition builders multiply by
 * `countvect(i)`, the MULTIPLICITY of the source phase, and that factor is
 * where a general-K generalisation of the same code would differ from a pair
 * encoding; dropping the representation would silently drop the factor.
 *
 * ARITHMETIC. `double` only, in step with the solve layer: the T matrix needs
 * an ordered real Schur factorization and the two Sylvester equations need
 * Bartels-Stewart, all LAPACK. The construction itself is exact rational
 * arithmetic on the descriptors, but a `Matrix<T>` construction feeding a
 * double-only solve would only move the conversion one call inwards.
 */

#include <cstddef>
#include <vector>

#include "line/api/fj/fj_dist2fj.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace fj {

/** `build_Service_h.m`: the two-subtask phase process of one fork-join job. */
struct FjCodesServiceH {
    Matrix<double> service_phases;  ///< (m^2) x (2m) count vectors, [long, short]
    std::vector<double> beta;       ///< kron(tau_st, tau_st), length m^2
    Matrix<double> S;               ///< kronsum(ST, ST), (m^2) x (m^2)
};

/** `build_SA.m`: the level-constant generator and the head-of-line jump. */
struct FjCodesSA {
    Matrix<double> S;       ///< no job completes, ((C+1) m^2) square
    Matrix<double> A_jump;  ///< a job completes and the next enters service
};

/** `generateService.m`: the service process seen by a tagged job. */
struct FjCodesGenService {
    Matrix<double> T;               ///< (newdim + dim_notbusy) square
    std::size_t newdim = 0;         ///< the all-busy part
    std::size_t dim_notbusy = 0;    ///< the not-all-busy part
};

/** `constructSRK.m`: the extended generator and the busy/idle projectors. */
struct FjCodesSRK {
    Matrix<double> Se;      ///< busy and not-busy phases together
    Matrix<double> Sestar;  ///< the busy-to-not-busy block of Se, in place
    Matrix<double> R0;      ///< not-busy to busy, on an arrival
    Matrix<double> Ke;      ///< newdim x (newdim + dim_notbusy), the busy rows
    Matrix<double> Kc;      ///< (newdim + dim_notbusy) x newdim, the busy columns
};

/**
 * Port of `build_index.m`: the compositions of `cr` into `m` non-negative
 * parts, one per row, in the reference's own order.
 *
 * FJ_codes calls this with cr = 1 only, where it is the identity matrix and the
 * rows are the m single-subtask phases. The general case is ported because it
 * is what the reference computes, and because the row ORDER is what
 * `vectmatch` resolves against: any other enumeration of the same set would
 * permute every block of every matrix below.
 */
inline Matrix<double> fj_build_index(std::size_t m, std::size_t cr) {
    if (m == 0) throw InputError("fj_build_index: m must be positive");
    // nchoosek(cr + m - 1, cr), in exact integer arithmetic
    std::size_t total = 1;
    for (std::size_t k = 1; k <= cr; ++k) total = total * (cr + m - k) / k;
    Matrix<double> idx(total, m, 0.0);
    idx(0, 0) = static_cast<double>(cr);
    for (std::size_t row = 1; row < total; ++row) {
        std::size_t k = m;  // MATLAB find(..., 1) on an all-zero row yields empty
        for (std::size_t j = 0; j < m; ++j)
            if (idx(row - 1, j) > 0.0) {
                k = j;
                break;
            }
        if (k + 1 < m) {  // the reference's `k < m`, with k 1-based there
            for (std::size_t j = 0; j < m; ++j) idx(row, j) = idx(row - 1, j);
            idx(row, k + 1) += 1.0;
            idx(row, 0) = idx(row, k) - 1.0;
            for (std::size_t j = 1; j <= k; ++j) idx(row, j) = 0.0;
        }
    }
    return idx;
}

/**
 * Port of `vectmatch.m`: the row of `matrix` equal to `row`.
 *
 * Returns a 0-based index. The reference leaves its output unassigned when
 * there is no match, which MATLAB reports as an error one frame up; the same
 * condition is a construction bug here and is named as one.
 */
inline std::size_t fj_vectmatch(const std::vector<double>& row, const Matrix<double>& matrix) {
    if (row.size() != matrix.cols())
        throw InputError("fj_vectmatch: the probe and the table have different widths");
    for (std::size_t i = 0; i < matrix.rows(); ++i) {
        bool hit = true;
        for (std::size_t j = 0; j < matrix.cols(); ++j)
            if (matrix(i, j) != row[j]) {
                hit = false;
                break;
            }
        if (hit) return i;
    }
    throw NumericError("fj_vectmatch: the target phase is not in the phase table, so the "
                       "fork-join state space is not closed under its own transitions");
}

/** The reference's `A = -sum(ST, 2) * tau_st`: complete, then restart. */
inline Matrix<double> fj_restart_matrix(const FjDist<double>& service) {
    const std::size_t m = service.tau_st.size();
    Matrix<double> A(m, m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < m; ++j) s += service.ST(i, j);
        for (std::size_t j = 0; j < m; ++j) A(i, j) = -s * service.tau_st[j];
    }
    return A;
}

/** Port of `build_Service_h.m`. */
inline FjCodesServiceH fj_build_service_h(const FjDist<double>& service) {
    const std::size_t m = service.tau_st.size();
    if (m == 0) throw InputError("fj_build_service_h: the service process has no phases");
    const Matrix<double> single = fj_build_index(m, 1);

    FjCodesServiceH h;
    h.service_phases = Matrix<double>(m * m, 2 * m, 0.0);
    std::size_t k = 0;
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) {
            for (std::size_t c = 0; c < m; ++c) {
                h.service_phases(k, c) = single(i, c);          // longest queue
                h.service_phases(k, m + c) = single(j, c);      // shortest queue
            }
            ++k;
        }
    h.beta.assign(m * m, 0.0);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j)
            h.beta[i * m + j] = service.tau_st[i] * service.tau_st[j];
    h.S = mam::krons(service.ST, service.ST);
    return h;
}

/**
 * Port of `build_SA.m`.
 *
 * `C` is the truncation of the queue-length difference and must be at least
 * one: the reference indexes `A_jump((C-1)*dim+1 : C*dim)` unconditionally, so
 * C = 0 is an out-of-range write there rather than a degenerate model here.
 */
inline FjCodesSA fj_build_sa(const FjDist<double>& service, const FjCodesServiceH& h,
                             std::size_t C) {
    if (C < 1)
        throw InputError("fj_build_sa: the FJ_codes truncation C must be at least 1 (the "
                         "reference indexes the c = C - 1 block unconditionally)");
    const std::size_t dim = h.beta.size();
    const std::size_t m = service.tau_st.size();
    const std::size_t dim_C = C + 1;
    const std::size_t newdim = dim_C * dim;

    FjCodesSA out;
    out.S = Matrix<double>(newdim, newdim, 0.0);
    out.A_jump = Matrix<double>(newdim, newdim, 0.0);

    for (std::size_t b = 0; b < dim_C; ++b)
        for (std::size_t i = 0; i < dim; ++i)
            for (std::size_t j = 0; j < dim; ++j) out.S(b * dim + i, b * dim + j) = h.S(i, j);

    const Matrix<double> A = fj_restart_matrix(service);

    // The job in the LONGER queue completes: c drops by one.
    Matrix<double> S_Cminus1(dim, dim, 0.0);
    for (std::size_t row = 0; row < dim; ++row)
        for (std::size_t i = 0; i < m; ++i) {
            const double cnt = h.service_phases(row, i);
            if (cnt <= 0.0) continue;
            for (std::size_t j = 0; j < m; ++j) {
                std::vector<double> to(2 * m);
                for (std::size_t c = 0; c < 2 * m; ++c) to[c] = h.service_phases(row, c);
                to[i] -= 1.0;
                to[j] += 1.0;
                S_Cminus1(row, fj_vectmatch(to, h.service_phases)) += cnt * A(i, j);
            }
        }
    for (std::size_t b = 0; b + 1 < dim_C; ++b)
        for (std::size_t i = 0; i < dim; ++i)
            for (std::size_t j = 0; j < dim; ++j)
                out.S(b * dim + i, (b + 1) * dim + j) = S_Cminus1(i, j);

    // The job in the SHORTER queue completes: c grows by one. The reference
    // does NOT weight this one by countvect(i), unlike the two around it.
    Matrix<double> A_Cplus1(dim, dim, 0.0);
    for (std::size_t row = 0; row < dim; ++row)
        for (std::size_t i = m; i < 2 * m; ++i) {
            if (h.service_phases(row, i) <= 0.0) continue;
            for (std::size_t j = m; j < 2 * m; ++j) {
                std::vector<double> to(2 * m);
                for (std::size_t c = 0; c < 2 * m; ++c) to[c] = h.service_phases(row, c);
                to[i] -= 1.0;
                to[j] += 1.0;
                A_Cplus1(row, fj_vectmatch(to, h.service_phases)) += A(i - m, j - m);
            }
        }
    for (std::size_t b = 1; b + 1 < dim_C; ++b)
        for (std::size_t i = 0; i < dim; ++i)
            for (std::size_t j = 0; j < dim; ++j)
                out.A_jump(b * dim + i, (b - 1) * dim + j) = A_Cplus1(i, j);
    for (std::size_t i = 0; i < dim; ++i)
        for (std::size_t j = 0; j < dim; ++j) out.A_jump(i, j) = A_Cplus1(i, j);

    // Both queues have the same length: whichever subtask does not finish
    // becomes the longer one.
    Matrix<double> A_last(dim, dim, 0.0);
    for (std::size_t row = 0; row < dim; ++row)
        for (std::size_t k = 0; k < 2; ++k)
            for (std::size_t i = 0; i < m; ++i) {
                const double cnt = h.service_phases(row, k * m + i);
                if (cnt <= 0.0) continue;
                std::vector<double> other(m);
                for (std::size_t c = 0; c < m; ++c)
                    other[c] = h.service_phases(row, (1 - k) * m + c);
                for (std::size_t j = 0; j < m; ++j) {
                    std::vector<double> to(2 * m, 0.0);
                    for (std::size_t c = 0; c < m; ++c) to[c] = other[c];
                    to[m + j] = 1.0;
                    A_last(row, fj_vectmatch(to, h.service_phases)) += cnt * A(i, j);
                }
            }
    for (std::size_t i = 0; i < dim; ++i)
        for (std::size_t j = 0; j < dim; ++j)
            out.A_jump(C * dim + i, (C - 1) * dim + j) = A_last(i, j);

    return out;
}

/**
 * The busy-to-not-busy blocks `S_long` and `S_last`, shared verbatim by
 * `generateService.m` and `constructSRK.m`.
 *
 * S_long: the subtask in the SHORTER queue completes while c > 0, so the
 * remaining subtask is the one that was in the longer queue.
 * S_last: the two queues have equal length, so either subtask may be the one
 * that completes and the other one is what remains.
 */
inline void fj_busy_to_idle(const FjDist<double>& service, const FjCodesServiceH& h,
                            const Matrix<double>& idle_phases, Matrix<double>& S_long,
                            Matrix<double>& S_last) {
    const std::size_t dim = h.beta.size();
    const std::size_t m = service.tau_st.size();
    const std::size_t dim_NB = idle_phases.rows();
    S_long = Matrix<double>(dim, dim_NB, 0.0);
    S_last = Matrix<double>(dim, dim_NB, 0.0);

    for (std::size_t row = 0; row < dim; ++row)
        for (std::size_t i = m; i < 2 * m; ++i) {
            const double cnt = h.service_phases(row, i);
            if (cnt <= 0.0) continue;
            std::vector<double> to(m);
            for (std::size_t c = 0; c < m; ++c) to[c] = h.service_phases(row, c);
            S_long(row, fj_vectmatch(to, idle_phases)) += cnt * service.St[i - m];
        }

    for (std::size_t row = 0; row < dim; ++row)
        for (std::size_t k = 0; k < 2; ++k)
            for (std::size_t i = k * m; i < (k + 1) * m; ++i) {
                const double cnt = h.service_phases(row, i);
                if (cnt <= 0.0) continue;
                std::vector<double> to(m);
                for (std::size_t c = 0; c < m; ++c)
                    to[c] = h.service_phases(row, (1 - k) * m + c);
                S_last(row, fj_vectmatch(to, idle_phases)) += cnt * service.St[i - k * m];
            }
}

/** Port of `generateService.m`. */
inline FjCodesGenService fj_generate_service(const FjDist<double>& service,
                                             const FjCodesServiceH& h, std::size_t C,
                                             const Matrix<double>& S) {
    const std::size_t dim = h.beta.size();
    const std::size_t m = service.tau_st.size();
    const Matrix<double> idle = fj_build_index(m, 1);
    const std::size_t dim_NB = idle.rows();
    const std::size_t dim_C = C + 1;

    FjCodesGenService out;
    out.newdim = dim_C * dim;
    out.dim_notbusy = dim_C * dim_NB;
    const std::size_t n = out.newdim + out.dim_notbusy;
    out.T = Matrix<double>(n, n, 0.0);

    std::vector<double> t(n, 0.0);
    for (std::size_t i = 0; i < dim_NB; ++i) t[n - dim_NB + i] = service.St[i];

    for (std::size_t i = 0; i < out.newdim; ++i)
        for (std::size_t j = 0; j < out.newdim; ++j) out.T(i, j) = S(i, j);

    Matrix<double> S_long, S_last;
    fj_busy_to_idle(service, h, idle, S_long, S_last);

    for (std::size_t b = 0; b + 1 < dim_C; ++b)
        for (std::size_t i = 0; i < dim; ++i)
            for (std::size_t j = 0; j < dim_NB; ++j)
                out.T(b * dim + i, out.newdim + b * dim_NB + j) = S_long(i, j);
    for (std::size_t i = 0; i < dim; ++i)
        for (std::size_t j = 0; j < dim_NB; ++j)
            out.T((dim_C - 1) * dim + i, out.newdim + (dim_C - 1) * dim_NB + j) = S_last(i, j);

    for (std::size_t b = 0; b < dim_C; ++b)
        for (std::size_t i = 0; i < dim_NB; ++i)
            for (std::size_t j = 0; j < dim_NB; ++j)
                out.T(out.newdim + b * dim_NB + i, out.newdim + b * dim_NB + j) = service.ST(i, j);

    const Matrix<double> A = fj_restart_matrix(service);
    for (std::size_t b = 0; b + 1 < dim_C; ++b)
        for (std::size_t i = 0; i < dim_NB; ++i)
            for (std::size_t j = 0; j < dim_NB; ++j)
                out.T(out.newdim + b * dim_NB + i, out.newdim + (b + 1) * dim_NB + j) = A(i, j);

    for (std::size_t row = out.newdim; row < n; ++row) {
        out.T(row, row) = 0.0;
        double s = 0.0;
        for (std::size_t j = 0; j < n; ++j) s += out.T(row, j);
        out.T(row, row) = -s - t[row];
    }
    return out;
}

/** Port of `constructNotAllBusy.m`. */
inline Matrix<double> fj_construct_not_all_busy(std::size_t C, const FjDist<double>& service,
                                                const FjCodesServiceH& h) {
    (void)h;
    if (C < 1)
        throw InputError("fj_construct_not_all_busy: the FJ_codes truncation C must be at least 1");
    const std::size_t m = service.tau_st.size();
    const std::size_t dim_NB = fj_build_index(m, 1).rows();
    const std::size_t dim_C = C + 1;
    const std::size_t n = (dim_C - 1) * dim_NB + 1;
    Matrix<double> out(n, n, 0.0);

    for (std::size_t b = 0; b + 1 < dim_C; ++b)
        for (std::size_t i = 0; i < dim_NB; ++i)
            for (std::size_t j = 0; j < dim_NB; ++j)
                out(b * dim_NB + i, b * dim_NB + j) = service.ST(i, j);

    const Matrix<double> A = fj_restart_matrix(service);
    for (std::size_t b = 0; b + 2 < dim_C; ++b)
        for (std::size_t i = 0; i < dim_NB; ++i)
            for (std::size_t j = 0; j < dim_NB; ++j)
                out(b * dim_NB + i, (b + 1) * dim_NB + j) = A(i, j);

    // The last column is the empty state: the one remaining subtask completes.
    for (std::size_t i = 0; i < dim_NB; ++i)
        out((dim_C - 2) * dim_NB + i, n - 1) = service.St[i];

    for (std::size_t row = 0; row < n; ++row) {
        out(row, row) = 0.0;
        double s = 0.0;
        for (std::size_t j = 0; j < n; ++j) s += out(row, j);
        out(row, row) = -s;
    }
    return out;
}

/** Port of `constructSRK.m`. */
inline FjCodesSRK fj_construct_srk(std::size_t C, const FjDist<double>& service,
                                   const FjCodesServiceH& h, const Matrix<double>& S) {
    if (C < 1) throw InputError("fj_construct_srk: the FJ_codes truncation C must be at least 1");
    const std::size_t dim = h.beta.size();
    const std::size_t m = service.tau_st.size();
    const Matrix<double> idle = fj_build_index(m, 1);
    const std::size_t dim_NB = idle.rows();
    const std::size_t dim_C = C + 1;
    const std::size_t newdim = dim_C * dim;
    const std::size_t dim_notbusy = (dim_C - 1) * dim_NB + 1;
    const std::size_t n = newdim + dim_notbusy;

    FjCodesSRK out;
    out.Se = Matrix<double>(n, n, 0.0);
    out.Sestar = Matrix<double>(n, n, 0.0);
    for (std::size_t i = 0; i < newdim; ++i)
        for (std::size_t j = 0; j < newdim; ++j) out.Se(i, j) = S(i, j);

    Matrix<double> S_long, S_last;
    fj_busy_to_idle(service, h, idle, S_long, S_last);

    // Block row 0 lands on idle block 0; block row b >= 1 lands on idle block
    // b - 1. The offset by one is the reference's, and is what distinguishes
    // constructSRK's idle space (which carries the single empty state at the
    // end) from generateService's.
    for (std::size_t i = 0; i < dim; ++i)
        for (std::size_t j = 0; j < dim_NB; ++j) out.Se(i, newdim + j) = S_long(i, j);
    for (std::size_t b = 1; b + 1 < dim_C; ++b)
        for (std::size_t i = 0; i < dim; ++i)
            for (std::size_t j = 0; j < dim_NB; ++j)
                out.Se(b * dim + i, newdim + (b - 1) * dim_NB + j) = S_long(i, j);
    for (std::size_t i = 0; i < dim; ++i)
        for (std::size_t j = 0; j < dim_NB; ++j)
            out.Se((dim_C - 1) * dim + i, newdim + (dim_C - 2) * dim_NB + j) = S_last(i, j);

    for (std::size_t i = 0; i < newdim; ++i)
        for (std::size_t j = newdim; j < n; ++j) out.Sestar(i, j) = out.Se(i, j);

    for (std::size_t b = 0; b + 1 < dim_C; ++b)
        for (std::size_t i = 0; i < dim_NB; ++i)
            for (std::size_t j = 0; j < dim_NB; ++j)
                out.Se(newdim + b * dim_NB + i, newdim + b * dim_NB + j) = service.ST(i, j);

    const Matrix<double> A = fj_restart_matrix(service);
    for (std::size_t b = 0; b + 2 < dim_C; ++b)
        for (std::size_t i = 0; i < dim_NB; ++i)
            for (std::size_t j = 0; j < dim_NB; ++j)
                out.Se(newdim + b * dim_NB + i, newdim + (b + 1) * dim_NB + j) = A(i, j);

    for (std::size_t i = 0; i < dim_NB; ++i)
        out.Se(newdim + (dim_C - 2) * dim_NB + i, n - 1) = service.St[i];

    for (std::size_t row = newdim; row < n; ++row) {
        out.Se(row, row) = 0.0;
        double s = 0.0;
        for (std::size_t j = 0; j < n; ++j) s += out.Se(row, j);
        out.Se(row, row) = -s;
    }

    // R0: an arrival during a not-all-busy period restarts the second branch,
    // so the system re-enters the all-busy phase with the same long-queue
    // phase and a fresh short-queue phase drawn from tau_st.
    out.R0 = Matrix<double>(n, n, 0.0);
    Matrix<double> R_NB(dim_NB, dim, 0.0);
    for (std::size_t row = 0; row < dim_NB; ++row)
        for (std::size_t i = 0; i < m; ++i) {
            std::vector<double> to(2 * m, 0.0);
            for (std::size_t c = 0; c < m; ++c) to[c] = idle(row, c);
            to[m + i] = 1.0;
            R_NB(row, fj_vectmatch(to, h.service_phases)) += service.tau_st[i];
        }
    for (std::size_t b = 0; b + 1 < dim_C; ++b)
        for (std::size_t i = 0; i < dim_NB; ++i)
            for (std::size_t j = 0; j < dim; ++j)
                out.R0(newdim + b * dim_NB + i, b * dim + j) = R_NB(i, j);
    for (std::size_t j = 0; j < dim; ++j) out.R0(n - 1, newdim - dim + j) = h.beta[j];

    out.Ke = Matrix<double>(newdim, n, 0.0);
    for (std::size_t i = 0; i < newdim; ++i) out.Ke(i, i) = 1.0;
    out.Kc = Matrix<double>(n, newdim, 0.0);
    for (std::size_t i = 0; i < newdim; ++i) out.Kc(i, i) = 1.0;
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_FJ_CODES_MATRICES_H
