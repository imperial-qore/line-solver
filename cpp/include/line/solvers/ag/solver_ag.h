/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_AG_SOLVER_AG_H
#define LINE_SOLVERS_AG_SOLVER_AG_H

/**
 * Port of `solver_ag.m`: the RCAT (Reversed Compound Agent Theorem)
 * analyzers, reached by methods 'inap', 'inapplus', 'inapinf' and 'exact'.
 *
 * THE METHOD. Each (station, class) pair that carries jobs becomes an isolated
 * CTMC, and the pairs are coupled only through the REVERSED RATES x_l of the
 * synchronizing actions. A component is a QBD whose LEVEL is the queue length
 * and whose PHASE is the pair (arrival phase, service phase), laid out in the
 * Kronecker order of qbd_mapmap1.h: an arrival moves the level up carrying
 * kron(D1^a, I), a service completion moves it down carrying kron(I, D1^s),
 * the busy levels evolve under krons(D0^a, D0^s) and level zero under
 * kron(D0^a, I), because no server is running there. With exponential
 * processes every block is 1 x 1 and the QBD collapses to the scalar
 * birth-death chain this analyzer built before, entry for entry. A departure from one
 * component is an active transition there and a passive one at the destination;
 * RCAT says that if the reversed rate of every active label is
 * state-independent, the joint chain has a product form whose factors are the
 * isolated components solved with the passive rates set to those x_l. INAP is
 * the fixed point that looks for such an x: solve the components, re-estimate
 * each x_l from the resulting marginals, repeat.
 *
 * THESE ARE APPROXIMATIONS, and the reference is explicit about it. The
 * reversed rate is state-independent only for genuinely product-form models; on
 * everything else INAP converges to an x that is merely a good average, and the
 * marginals it returns are not the model's. The 'inapinf' variant reports the
 * RCAT residual of Remark 2, max_l ||pi (x_l I - T_l)||, which is zero exactly
 * when the product form is real, and it is exposed here as `rcat_residual` for
 * the same reason: it is the only honest indication of how far off the answer
 * is. Nothing in this header should be compared against an exact solver at a
 * tolerance that pretends otherwise.
 *
 * THE THREE VARIANTS differ only in how x_l is re-estimated and how the open
 * components are solved:
 *   inap      x_l = mean over the support of the state-wise reversed rate
 *             (pi A_l)_j / pi_j, which on a birth-death component is the
 *             entrywise mean of A_l(i,j) pi(i)/pi(j) term for term
 *   inapplus  x_l = sum over the support of A_l(i,j) pi(i), the rate-conserving
 *             estimator, which INAP also switches to on any component that is
 *             not a birth-death chain (a catastrophe or batch removal reaches
 *             beyond the neighbouring LEVEL, and the mean-of-ratios estimator
 *             is meaningless there) and on any component with more than one
 *             phase per level, where the same failure appears
 *   inapinf   x_l from the closed-form geometric tail, with each OPEN component
 *             solved on its infinite state space by the scalar QBD root instead
 *             of being truncated at maxStates (Marin, Rota Bulo and Balsamo,
 *             MASCOTS 2012)
 *
 * TWO REFERENCE QUIRKS THAT LOOK LIKE BUGS AND ARE HARMLESS, both worth knowing
 * before editing. First, `compute_equilibrium` subtracts unscaled row-sum
 * diagonals while adding x(c)-scaled off-diagonal blocks, which is inconsistent;
 * it does not matter because `ctmc_makeinfgen` discards the diagonal outright
 * and rebuilds it from the row sums. Second, and for the same reason, the
 * self-service term `L(n,n) = mu * P(self)` that `build_local_rates` writes on
 * the DIAGONAL is discarded too, so self-routing at a station has no effect on
 * the answer. Both are reproduced rather than corrected: correcting either
 * would change every number this analyzer returns.
 *
 * NOT PORTED, AND WHY. The reference's self-looping-class override reads
 * `sn.isslc`, which NetworkStruct does not carry. MATLAB guards that block with
 * `isfield(sn,'isslc')`, so a struct without the field skips it and the port
 * takes the same path; a model that genuinely has self-looping classes would
 * differ, and there is no way to detect one here without inventing the concept.
 *
 * ARITHMETIC. The fixed point stops on a tolerance and the QBD root takes a
 * square root, so the whole body is gated on transcendental arithmetic, as
 * solver_mam_basic.h and solver_mam_ldqbd.h are.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/qbd_r.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/ag/ag_exec.h"
#include "line/solvers/ag/ag_types.h"
#include "line/solvers/mam/mam_types.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace ag {

// Kronecker products and Neuts' logarithmic reduction are MAM primitives, not
// RCAT ones: they live in api/mam and are shared with every matrix-analytic
// analyzer. Named explicitly rather than pulled in wholesale, so the split
// between what AG owns and what it merely uses stays visible.
using mam::kron;
using mam::krons;
using mam::qbd_R_logred;

namespace ag_detail {

/** One synchronizing action: a departure from (from) that arrives at (to). */
struct RcatAction {
    std::size_t from_station = 0, from_class = 0;  ///< 0-based
    std::size_t to_station = 0, to_class = 0;      ///< 0-based
    double prob = 0.0;
    bool is_negative = false;
    bool is_catastrophe = false;
    std::size_t removal_class = 0;  ///< 0-based class whose signalremdist applies
    bool has_removal_dist = false;
};

/** The RCAT form of the network: components, actions, and their rate matrices. */
template <class T>
struct RcatModel {
    std::size_t num_processes = 0;
    /** 1-based process id per (station, class), 0 where the pair carries no jobs. */
    std::vector<std::vector<std::size_t>> process_map;
    std::vector<RcatAction> actions;
    std::vector<std::size_t> N;      ///< state count per process, nlev * mph
    std::vector<Matrix<T>> Aa, Pb;   ///< active and passive matrix per action
    std::vector<Matrix<T>> L;        ///< local (hidden) rate matrix per process
    std::vector<std::size_t> act, psv;  ///< 0-based active/passive process per action
    std::vector<bool> is_open_proc;
    std::vector<std::size_t> nlev;   ///< QBD levels per process
    std::vector<std::size_t> mph;    ///< phases per level per process
    std::vector<std::vector<std::size_t>> level;   ///< level index of every state
    /** Service completion rate out of every state; zero on level 0. */
    std::vector<std::vector<T>> svcrate;
    /** The same rate per phase at a busy level. */
    std::vector<std::vector<T>> svcdown;
};

/**
 * One component's Markovian processes: the service MAP of its station and the
 * arrival MAP of the external streams reaching it, plus the removal signals,
 * which stay scalar because a signal is a trigger with no service.
 */
template <class T>
struct RcatComponent {
    std::size_t ist = 0, r = 0;
    Matrix<T> Da0, Da1;   ///< arrival MAP of the external streams
    Matrix<T> Ds0, Ds1;   ///< service MAP
    Matrix<T> Dsvc;       ///< kron(I_na, D1^s): a completion, level down
    std::size_t na = 1, ns = 1, mph = 1, nlev = 0, N = 0;
    T lam_neg = num_traits<T>::from_int(0);
    T lam_cat = num_traits<T>::from_int(0);
    std::vector<std::pair<T, std::size_t>> batch;  ///< (rate, signal class)
};

/** Row/column offset of level N (0-based) in a component with MPH phases. */
inline std::size_t blk(std::size_t n, std::size_t mph) { return n * mph; }

/** Elementwise sum of two same-shaped matrices. */
template <class T>
Matrix<T> madd_local(const Matrix<T>& A, const Matrix<T>& B) {
    Matrix<T> C(A.rows(), A.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) = T(A(i, j) + B(i, j));
    return C;
}

/** Add SCALE times BLOCK into the (LI, LJ) level block of TARGET. */
template <class T>
void add_block(Matrix<T>& target, std::size_t li, std::size_t lj, std::size_t mph,
               const Matrix<T>& block, const T& scale) {
    const std::size_t r0 = blk(li, mph), c0 = blk(lj, mph);
    for (std::size_t i = 0; i < mph; ++i)
        for (std::size_t j = 0; j < mph; ++j) target(r0 + i, c0 + j) += T(scale * block(i, j));
}

/** Add SCALE times the identity into the (LI, LJ) level block of TARGET. */
template <class T>
void add_identity(Matrix<T>& target, std::size_t li, std::size_t lj, std::size_t mph,
                  const T& scale) {
    const std::size_t r0 = blk(li, mph), c0 = blk(lj, mph);
    for (std::size_t i = 0; i < mph; ++i) target(r0 + i, c0 + i) += scale;
}

/** Set the (LI, LJ) level block of TARGET to the identity. */
template <class T>
void set_identity(Matrix<T>& target, std::size_t li, std::size_t lj, std::size_t mph) {
    const std::size_t r0 = blk(li, mph), c0 = blk(lj, mph);
    const T one = num_traits<T>::from_int(1);
    for (std::size_t i = 0; i < mph; ++i) target(r0 + i, c0 + i) = one;
}

/** The (LI, LJ) level block of Q. */
template <class T>
Matrix<T> level_block(const Matrix<T>& Q, std::size_t li, std::size_t lj, std::size_t mph) {
    Matrix<T> out(mph, mph, num_traits<T>::from_int(0));
    const std::size_t r0 = blk(li, mph), c0 = blk(lj, mph);
    for (std::size_t i = 0; i < mph; ++i)
        for (std::size_t j = 0; j < mph; ++j) out(i, j) = Q(r0 + i, c0 + j);
    return out;
}

/**
 * True when (D0, D1) is a genuine MAP rather than a RAP or an ME process:
 * non-negative off-diagonal rates in D0, non-negative rates in D1, and
 * (D0 + D1) an infinitesimal generator. A CTMC assembled from anything else is
 * a rational generator whose stationary solution is a signed vector.
 */
template <class T>
bool is_markovian_map(const Matrix<T>& D0, const Matrix<T>& D1) {
    const std::size_t n = D0.rows();
    if (n == 0 || D0.cols() != n || D1.rows() != n || D1.cols() != n) return false;
    double scale = 1.0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            const double a = num_traits<T>::to_double(D0(i, j));
            const double b = num_traits<T>::to_double(D1(i, j));
            if (!std::isfinite(a) || !std::isfinite(b)) return false;
            scale = std::max(scale, std::max(std::fabs(a), std::fabs(b)));
        }
    const double tol = 1e-9 * scale;
    for (std::size_t i = 0; i < n; ++i) {
        double row = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            const double a = num_traits<T>::to_double(D0(i, j));
            const double b = num_traits<T>::to_double(D1(i, j));
            if (i != j && a < -tol) return false;
            if (b < -tol) return false;
            row += a + b;
        }
        if (std::fabs(row) > tol) return false;
    }
    return true;
}

/**
 * (D0,D1) of the process at (IST,R), or the exponential pair built from
 * `rates` when the station carries no usable matrix representation.
 *
 * A non-Markovian pair is refused here and answered as its mean rate; the
 * runner-level gate rejects those models before they reach this point.
 */
template <class T>
std::pair<Matrix<T>, Matrix<T>> proc_map(const qn::NetworkStruct<T>& L, std::size_t ist,
                                         std::size_t r) {
    const T zero = num_traits<T>::from_int(0);
    if (L.has_service_law(ist, r)) {
        try {
            const mam::Map<T> m = lang::dist_to_map(L.service[ist][r]);
            if (is_markovian_map(m.D0, m.D1)) return std::make_pair(m.D0, m.D1);
        } catch (const std::exception&) {
            // no usable representation: fall through to the exponential pair
        }
    }
    T rate = L.disabled[ist][r] ? zero : L.rates(ist, r);
    if (!(rate > zero) || !std::isfinite(num_traits<T>::to_double(rate))) rate = zero;
    Matrix<T> D0(1, 1, zero), D1(1, 1, zero);
    D0(0, 0) = T(zero - rate);
    D1(0, 0) = rate;
    return std::make_pair(D0, D1);
}

template <class T>
bool is_tridiagonal(const Matrix<T>& Q) {
    const std::size_t n = Q.rows();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            const std::size_t d = (i > j) ? (i - j) : (j - i);
            if (d > 1 && std::fabs(num_traits<T>::to_double(Q(i, j))) > 1e-14) return false;
        }
    return true;
}

/**
 * True when every transition of Q stays within the neighbouring level, LVL
 * being the level index of each state.
 */
template <class T>
bool is_block_tridiagonal(const Matrix<T>& Q, const std::vector<std::size_t>& lvl) {
    const std::size_t n = Q.rows();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            const std::size_t d = (lvl[i] > lvl[j]) ? (lvl[i] - lvl[j]) : (lvl[j] - lvl[i]);
            if (d > 1 && std::fabs(num_traits<T>::to_double(Q(i, j))) > 1e-14) return false;
        }
    return true;
}

/**
 * Equilibrium of a birth-death chain by the ratio recursion.
 *
 * Preferred over the null-space solve wherever the generator is tridiagonal
 * because the recursion is stable at the state-space sizes this analyzer builds
 * (maxStates is 100 by default), where the linear system is already
 * ill-conditioned.
 */
template <class T>
std::vector<T> birth_death_solve(const Matrix<T>& Q) {
    const std::size_t n = Q.rows();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (n <= 1) return std::vector<T>(1, one);

    std::vector<T> pi(n, zero);
    pi[0] = one;
    for (std::size_t i = 1; i < n; ++i) {
        const T birth = Q(i - 1, i);
        const T death = Q(i, i - 1);
        pi[i] = (death > zero) ? T(pi[i - 1] * birth / death) : zero;
    }
    T total = zero;
    for (const T& v : pi) total += v;
    if (total > zero) {
        for (T& v : pi) v /= total;
    } else {
        const T u = one / num_traits<T>::from_int(static_cast<long>(n));
        for (T& v : pi) v = u;
    }
    return pi;
}

/**
 * Stationary vector of a generator C, allowing a reducible one.
 *
 * Replaces the first balance equation by the normalization, which is the
 * equation it is redundant with (the columns of a generator sum to zero), and
 * solves the resulting square system. Unlike a null-space solve this stays well
 * posed when the chain is reducible with ONE closed class, which the level-0
 * chain of a phase-expanded component routinely is: a phase-type restarts in
 * the support of alpha, so every service phase outside that support is
 * unreachable once the queue has emptied at least once.
 */
template <class T>
std::vector<T> stat_vector(const Matrix<T>& C) {
    const std::size_t m = C.rows();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    // Solve v A = e_0, i.e. A^T v^T = e_0^T.
    Matrix<T> At(m, m, zero);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) At(j, i) = (j == 0) ? one : C(i, j);
    std::vector<T> rhs(m, zero);
    rhs[0] = one;
    const std::vector<std::size_t> piv = lu_factor(At);
    lu_solve(At, piv, rhs);
    return rhs;
}

/**
 * Stationary vector of a finite block-tridiagonal generator.
 *
 * Linear level reduction: censor the chain level by level from the top,
 *   C(nlev-1) = B(nlev-1),  C(n) = B(n) + F(n) (-C(n+1))^-1 D(n+1),
 * with B, F and D the diagonal, up and down blocks. C(0) is the generator of
 * the chain censored on level 0, so pi_0 is its stationary vector and the rest
 * follows from pi_(n+1) = pi_n F(n) (-C(n+1))^-1. This is the block form of
 * birth_death_solve and reduces to it entry for entry when m == 1.
 */
template <class T>
std::vector<T> qbd_finite_solve(const Matrix<T>& Q, std::size_t m, std::size_t nlev) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (nlev <= 1) return mc::ctmc_solve(Q);

    std::vector<Matrix<T>> C(nlev);
    C[nlev - 1] = level_block(Q, nlev - 1, nlev - 1, m);
    for (std::size_t t = nlev - 1; t-- > 0;) {
        const Matrix<T> F = level_block(Q, t, t + 1, m);
        const Matrix<T> D = level_block(Q, t + 1, t, m);
        Matrix<T> negC = C[t + 1];
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) negC(i, j) = T(zero - negC(i, j));
        const Matrix<T> step = matmul(F, inverse(negC));
        Matrix<T> Cn = level_block(Q, t, t, m);
        const Matrix<T> corr = matmul(step, D);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) Cn(i, j) += corr(i, j);
        C[t] = Cn;
    }

    std::vector<T> pi(nlev * m, zero);
    const std::vector<T> p0 = stat_vector(C[0]);
    for (std::size_t i = 0; i < m; ++i) pi[i] = p0[i];
    for (std::size_t n = 0; n + 1 < nlev; ++n) {
        Matrix<T> negC = C[n + 1];
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) negC(i, j) = T(zero - negC(i, j));
        const Matrix<T> step = matmul(level_block(Q, n, n + 1, m), inverse(negC));
        for (std::size_t j = 0; j < m; ++j) {
            T acc = zero;
            for (std::size_t i = 0; i < m; ++i) acc += T(pi[blk(n, m) + i] * step(i, j));
            pi[blk(n + 1, m) + j] = acc;
        }
    }

    T total = zero;
    for (const T& v : pi) total += v;
    if (total > zero) {
        for (T& v : pi) v /= total;
    } else {
        const T u = one / num_traits<T>::from_int(static_cast<long>(nlev * m));
        for (T& v : pi) v = u;
    }
    return pi;
}

/**
 * Stationary vector of one isolated component.
 *
 * A component with a single phase per level is the birth-death chain the
 * analyzer has always built, and the ratio recursion is both exact and stable
 * there; a phase-expanded component is block tridiagonal instead, and the matrix
 * analogue of that recursion keeps the same stability at the 100-level
 * truncation, where a null-space solve is already ill-conditioned. Anything that
 * reaches beyond the neighbouring level -- a catastrophe, a batch removal -- is
 * neither, and falls back to ctmc_solve.
 */
template <class T>
std::vector<T> solve_component(const Matrix<T>& Q, std::size_t mph, std::size_t nlev,
                               const std::vector<std::size_t>& lvl) {
    if (mph == 1) {
        if (is_tridiagonal(Q)) return birth_death_solve(Q);
    } else if (is_block_tridiagonal(Q, lvl)) {
        return qbd_finite_solve(Q, mph, nlev);
    }
    return mc::ctmc_solve(Q);
}

/** The matrix-geometric tail of one open component with more than one phase. */
template <class T>
struct QbdTail {
    Matrix<T> R;
    std::vector<T> pi0, pi1;
    std::vector<T> busy;   ///< sum_{n>=1} pi_n = pi_1 (I - R)^-1
    T qlen = num_traits<T>::from_int(0);
    bool ok = false;
};

/**
 * Neuts' matrix-geometric solution of one open component with MPH phases per
 * level: R from logarithmic reduction, then the boundary equations of levels 0
 * and 1,
 *   pi_0 B00 + pi_1 A2 = 0,   pi_0 A0 + pi_1 (A1 + R A2) = 0,
 * normalized by pi_0 e + pi_1 (I - R)^-1 e = 1. Returns ok = false when R has no
 * sub-unit spectral radius, i.e. when the isolated component is unstable and has
 * no stationary tail to report.
 *
 * Logarithmic reduction rather than successive substitutions: this runs once per
 * component per fixed-point sweep, and the quadratic convergence is what keeps
 * that affordable. The cap is small for the same reason -- an unstable component
 * must bail rather than grind to the library default of 1e5.
 */
template <class T>
QbdTail<T> qbd_matrix_tail(const Matrix<T>& Qk, const Matrix<T>& A0, const Matrix<T>& A1,
                           const Matrix<T>& A2, std::size_t mph) {
    QbdTail<T> out;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    Matrix<T> R;
    try {
        R = qbd_R_logred(A2, A1, A0, 500u, T(num_traits<T>::from_double(1e-14)));
    } catch (const std::exception&) {
        return out;
    }
    for (std::size_t i = 0; i < mph; ++i)
        for (std::size_t j = 0; j < mph; ++j) {
            const double v = num_traits<T>::to_double(R(i, j));
            // The minimal solution of a QBD is NON-NEGATIVE; anything else is
            // the iteration having failed rather than a rate matrix.
            if (!std::isfinite(v) || v < -1e-12) return out;
        }

    const Matrix<T> B00 = level_block(Qk, 0, 0, mph);
    Matrix<T> IR = eye<T>(mph);
    for (std::size_t i = 0; i < mph; ++i)
        for (std::size_t j = 0; j < mph; ++j) IR(i, j) = T(IR(i, j) - R(i, j));
    std::vector<T> tail_mass;
    try {
        tail_mass = mulvec(inverse(IR), ones<T>(mph));
    } catch (const std::exception&) {
        return out;
    }
    // STABILITY WITHOUT AN EIGENSOLVER. (I-R)^-1 = I + R + R^2 + ... converges
    // exactly when the spectral radius is below one, and every row of that
    // series is e_i plus non-negative terms, so (I-R)^-1 e >= 1 entrywise. When
    // the isolated component is unstable the series diverges and the inverse
    // picks up negative entries, so the test below is the spectral condition
    // without the LAPACK dependency an eigensolve would carry.
    for (std::size_t i = 0; i < mph; ++i) {
        const double w = num_traits<T>::to_double(tail_mass[i]);
        if (!std::isfinite(w) || w < 1.0 - 1e-9) return out;
    }

    // [pi_0 pi_1] Sys = 0, with the first column replaced by the normalization.
    const Matrix<T> lower_right = madd_local(A1, matmul(R, A2));
    const std::size_t n2 = 2 * mph;
    Matrix<T> SysT(n2, n2, zero);   // transposed: solve Sys^T v = e_0
    for (std::size_t i = 0; i < mph; ++i)
        for (std::size_t j = 0; j < mph; ++j) {
            SysT(j, i) = B00(i, j);
            SysT(mph + j, i) = A0(i, j);
            SysT(j, mph + i) = A2(i, j);
            SysT(mph + j, mph + i) = lower_right(i, j);
        }
    for (std::size_t i = 0; i < mph; ++i) {
        SysT(0, i) = one;
        SysT(0, mph + i) = tail_mass[i];
    }
    std::vector<T> v(n2, zero);
    v[0] = one;
    try {
        Matrix<T> lu = SysT;
        const std::vector<std::size_t> piv = lu_factor(lu);
        lu_solve(lu, piv, v);
    } catch (const std::exception&) {
        return out;
    }
    for (std::size_t i = 0; i < n2; ++i)
        if (!std::isfinite(num_traits<T>::to_double(v[i]))) return out;

    out.pi0.assign(v.begin(), v.begin() + static_cast<long>(mph));
    out.pi1.assign(v.begin() + static_cast<long>(mph), v.end());
    const Matrix<T> IRinv = inverse(IR);
    out.busy = vecmul(out.pi1, IRinv);
    const std::vector<T> qv = vecmul(out.busy, IRinv);
    T q = zero;
    for (const T& t : qv) q += t;
    out.qlen = q;
    out.R = R;
    out.ok = true;
    return out;
}

/**
 * Materialize the matrix-geometric tail over NLEV levels, so the block norm of
 * the fixed point and the RCAT residual read one vector shape for every
 * component. The metrics use the closed forms in the tail instead.
 */
template <class T>
std::vector<T> qbd_tail_expand(const QbdTail<T>& g, std::size_t nlev, std::size_t mph) {
    std::vector<T> pi(nlev * mph, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < mph; ++i) pi[i] = g.pi0[i];
    std::vector<T> v = g.pi1;
    for (std::size_t n = 1; n < nlev; ++n) {
        for (std::size_t i = 0; i < mph; ++i) pi[blk(n, mph) + i] = v[i];
        v = vecmul(v, g.R);
    }
    return pi;
}

/** P(batch = k) off the pmf indexed by batch size, zero past its end. */
template <class T>
T pmf_at(const std::vector<T>& pmf, std::size_t k) {
    return k < pmf.size() ? pmf[k] : num_traits<T>::from_int(0);
}

/**
 * Accumulate the n -> m block of a batch removal into `B`, scaled by `rate`.
 *
 * Landing on the empty state absorbs the whole upper tail of the pmf, which is
 * what keeps the block stochastic once the batch exceeds the queue length.
 */
template <class T>
void add_batch_removal(Matrix<T>& B, const std::vector<T>& pmf, const T& rate,
                       std::size_t nlev, std::size_t mph) {
    const T one = num_traits<T>::from_int(1), zero = num_traits<T>::from_int(0);
    for (std::size_t n = 1; n < nlev; ++n) {
        for (std::size_t m = 1; m <= n; ++m) {
            const T p = pmf_at(pmf, n - m);
            if (p > zero) add_identity(B, n, m, mph, T(rate * p));
        }
        T cdf = zero;
        for (std::size_t j = 0; j + 1 <= n; ++j) cdf += pmf_at(pmf, j);
        const T tail = T(one - cdf);
        if (tail > zero) add_identity(B, n, 0, mph, T(rate * tail));
    }
}

/** Sub-unit root of b rho^2 - (f+b+g) rho + f = 0, the block-size-1 Neuts R. */
inline double qbd_scalar_rho(double f, double b, double g) {
    const double inf = std::numeric_limits<double>::infinity();
    if (b <= 1e-14) {
        // Degenerate without a down transition: the catastrophe drain is the
        // only thing that keeps the chain stable.
        if (f + g <= 0.0) return inf;
        return f / (f + g);
    }
    const double c1 = -(f + b + g);
    const double disc = c1 * c1 - 4.0 * b * f;
    if (disc < 0.0) return inf;
    const double sq = std::sqrt(disc);
    double r1 = (-c1 - sq) / (2.0 * b), r2 = (-c1 + sq) / (2.0 * b);
    if (r1 > r2) std::swap(r1, r2);
    return (r1 > 0.0) ? r1 : r2;
}

/**
 * `sn.issignal(r)` under the struct's own convention: the five G-network arrays
 * are allocated together, and only when `set_signal` declares one, so an EMPTY
 * `issignal` means "no class is a signal" and is the ORDINARY case rather than a
 * malformed struct.
 *
 * IT MUST BE READ THROUGH A SIZE TEST, and not merely because the subscript is
 * out of range. `issignal` is a `std::vector<bool>`, whose empty state holds a
 * NULL word pointer, so `issignal[r]` dereferences null and takes the process
 * down -- it does not read a stray byte and carry on. That is what `build_rcat`
 * did on the first signal-free model ever handed to it, Source -> Q1 -> Sink,
 * and the SIGSEGV surfaced in a test whose subject is the M/M/1 marginal, with
 * nothing about it pointing at G-networks. Every other reader of these arrays in
 * the tree already guards: `solver_ctmc.h`, `state_events.h`, `tag_chain.h`.
 */
template <class T>
bool is_signal_class(const qn::NetworkStruct<T>& L, std::size_t r) {
    return r < L.issignal.size() && L.issignal[r];
}

/**
 * The declared signal type, defaulting exactly as `set_signal` initializes the
 * array. Only ever consulted for a class `is_signal_class` admits, so the
 * default is unreachable through this header; it is stated so that the pair of
 * arrays cannot disagree about how far they extend.
 */
template <class T>
lang::SignalType signal_type_of(const qn::NetworkStruct<T>& L, std::size_t r) {
    return r < L.signaltype.size() ? L.signaltype[r] : lang::SignalType::NEGATIVE;
}

/**
 * External (Source) streams reaching the component, as one arrival MAP for the
 * positive customers plus the scalar rates of the removal signals.
 *
 * Each stream is thinned by its routing probability -- a MAP thinned with
 * probability p is (D0 + (1-p) D1, p D1) -- and the streams are superposed by
 * the Kronecker sum, so several Poisson sources still collapse to the single
 * rate sum this analyzer used before. Removal signals stay scalar: a signal is a
 * trigger with no service, and its arrival process is required exponential.
 */
template <class T>
void arrival_map(const qn::NetworkStruct<T>& L, RcatComponent<T>& c,
                 const std::vector<std::size_t>& source_stations) {
    using lang::SignalType;
    const std::size_t K = L.nclasses;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    bool have_arrival = false;
    Matrix<T> Da0(1, 1, zero), Da1(1, 1, zero);

    for (std::size_t isrc : source_stations)
        for (std::size_t s = 0; s < K; ++s) {
            const bool is_signal = is_signal_class(L, s);
            T prob_src = zero;
            if (is_signal) {
                // A signal routes to itself, so its effect on this component is
                // its total probability of reaching this STATION in any class.
                for (std::size_t sd = 0; sd < K; ++sd)
                    prob_src += L.rt(isrc * K + s, c.ist * K + sd);
            } else {
                prob_src = L.rt(isrc * K + s, c.ist * K + c.r);
            }
            if (!(prob_src > zero) || L.disabled[isrc][s]) continue;
            const T src_rate = L.rates(isrc, s);
            const bool removal_signal =
                is_signal && (signal_type_of(L, s) == SignalType::NEGATIVE ||
                              signal_type_of(L, s) == SignalType::CATASTROPHE);
            if (removal_signal) {
                if (signal_type_of(L, s) == SignalType::CATASTROPHE) {
                    c.lam_cat += T(src_rate * prob_src);
                } else if (s < L.signalremdist.size() && !L.signalremdist[s].empty()) {
                    c.batch.emplace_back(T(src_rate * prob_src), s);
                } else {
                    c.lam_neg += T(src_rate * prob_src);
                }
                continue;
            }
            if (!(src_rate > zero)) continue;

            std::pair<Matrix<T>, Matrix<T>> sm = proc_map(L, isrc, s);
            Matrix<T> S0 = sm.first, S1 = sm.second;
            if (prob_src < one) {
                const T keep = prob_src, drop = T(one - prob_src);
                for (std::size_t i = 0; i < S0.rows(); ++i)
                    for (std::size_t j = 0; j < S0.cols(); ++j) S0(i, j) += T(drop * S1(i, j));
                for (std::size_t i = 0; i < S1.rows(); ++i)
                    for (std::size_t j = 0; j < S1.cols(); ++j) S1(i, j) = T(keep * S1(i, j));
            }
            if (have_arrival) {
                Da0 = krons(Da0, S0);
                Da1 = krons(Da1, S1);
            } else {
                Da0 = S0;
                Da1 = S1;
                have_arrival = true;
            }
        }
    c.Da0 = Da0;
    c.Da1 = Da1;
}

/** Build the local (hidden) rate matrix of the component at (station, class). */
template <class T>
Matrix<T> build_local_rates(const qn::NetworkStruct<T>& L, const RcatComponent<T>& c,
                            const std::vector<std::size_t>& sink_nodes) {
    const std::size_t K = L.nclasses;
    const std::size_t mph = c.mph, nlev = c.nlev;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    Matrix<T> Lm(c.N, c.N, zero);

    // Level-local blocks: the arrival phase always runs, the service phase only
    // while the server is busy (qbd_mapmap1's Lbar = kron(D0^a, I) at level 0
    // and L = krons(D0^a, D0^s) above it). With one phase each these are pure
    // diagonals, which ctmc_makeinfgen discards and rebuilds from the row sums.
    add_block(Lm, 0, 0, mph, kron(c.Da0, eye<T>(c.ns)), one);
    const Matrix<T> Lbusy = krons(c.Da0, c.Ds0);
    for (std::size_t n = 1; n < nlev; ++n) add_block(Lm, n, n, mph, Lbusy, one);

    // Positive arrivals: level n -> n+1, carrying kron(D1^a, I).
    const Matrix<T> Aup = kron(c.Da1, eye<T>(c.ns));
    for (std::size_t n = 0; n + 1 < nlev; ++n) add_block(Lm, n, n + 1, mph, Aup, one);
    // At the truncation the job is lost but the arrival process still moves on,
    // so the block stays on the top level. With a single arrival phase this is a
    // pure diagonal and is discarded, exactly as before.
    add_block(Lm, nlev - 1, nlev - 1, mph, Aup, one);

    // Catastrophe arrivals: every busy level drops to level 0.
    if (c.lam_cat > zero)
        for (std::size_t n = 1; n < nlev; ++n) add_identity(Lm, n, 0, mph, c.lam_cat);
    for (const std::pair<T, std::size_t>& ba : c.batch)
        add_batch_removal(Lm, L.signalremdist[ba.second], ba.first, nlev, mph);
    // Single-removal negative arrivals: level n -> n-1 (busy levels only).
    if (c.lam_neg > zero)
        for (std::size_t n = 1; n < nlev; ++n) add_identity(Lm, n, n - 1, mph, c.lam_neg);

    // Service completions that are not synchronizing actions: departures to a
    // Sink (level down) and self-routing (level unchanged, service restarted).
    if (L.disabled[c.ist][c.r]) return Lm;
    const T mu = L.rates(c.ist, c.r);
    if (!(mu > zero)) return Lm;

    // The node index IS needed here -- rtnodes is indexed by node, not by
    // station -- so a station the map does not cover has no row to read and
    // contributes no sink probability. Tested up front rather than relying on
    // the range check below to catch the wrapped `node - 1`.
    const std::size_t node = L.node_of_station(c.ist + 1);
    T prob_sink = zero;
    if (node != 0)
        for (std::size_t jsnk : sink_nodes)
            for (std::size_t s = 0; s < K; ++s) {
                const std::size_t from = (node - 1) * K + c.r, to = (jsnk - 1) * K + s;
                if (from < L.rtnodes.rows() && to < L.rtnodes.cols())
                    prob_sink += L.rtnodes(from, to);
            }
    const T prob_self = L.rt(c.ist * K + c.r, c.ist * K + c.r);

    if (prob_sink > zero)
        for (std::size_t n = 1; n < nlev; ++n) add_block(Lm, n, n - 1, mph, c.Dsvc, prob_sink);
    // With one service phase this lands on the diagonal, hence is discarded by
    // ctmc_makeinfgen and self-routing has no effect, exactly as the reference
    // records; with a phase block it restarts the service, which is the physics.
    if (prob_self > zero)
        for (std::size_t n = 1; n < nlev; ++n) add_block(Lm, n, n, mph, c.Dsvc, prob_self);
    return Lm;
}

/** Port of `build_rcat`: the network to its RCAT components and actions. */
template <class T>
RcatModel<T> build_rcat(const qn::NetworkStruct<T>& L, std::size_t max_states) {
    using lang::SignalType;
    const std::size_t M = L.nstations, K = L.nclasses;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    RcatModel<T> m;
    m.process_map.assign(M, std::vector<std::size_t>(K, 0));

    std::vector<std::size_t> source_stations, queue_stations, sink_nodes;
    for (std::size_t i = 0; i < L.nof_nodes(); ++i)
        if (L.nodes[i].nodetype == qn::NodeType::Sink) sink_nodes.push_back(i + 1);
    for (std::size_t ist = 0; ist < M; ++ist) {
        // FROM THE STATION, not through the node map. `add_station` copies the
        // type into both, so the two always agree -- but `node_of_station`
        // returns 0 for a station the map does not cover (a Layer built by the
        // LQN path carries fewer entries than stations), and `nodes[node - 1]`
        // then wraps to SIZE_MAX and reads off the end of the vector.
        const qn::NodeType ty = L.stations[ist].nodetype;
        if (ty == qn::NodeType::Source) source_stations.push_back(ist);
        else if (ty != qn::NodeType::Sink) queue_stations.push_back(ist);
    }

    // A signal class never gets a component of its own: it has no queue, it
    // only edits the state of the positive-customer components it reaches.
    std::size_t pidx = 0;
    for (std::size_t ist : queue_stations)
        for (std::size_t r = 0; r < K; ++r) {
            if (is_signal_class(L, r) || L.disabled[ist][r]) continue;
            if (!(L.rates(ist, r) > zero)) continue;
            m.process_map[ist][r] = ++pidx;
        }
    m.num_processes = pidx;
    if (m.num_processes == 0) return m;

    // The QBD shape of every component: both the service MAP of its station and
    // the arrival MAP of the external streams reaching it.
    const std::vector<double> njobs = L.njobs();
    std::vector<RcatComponent<T>> comp(m.num_processes);
    m.N.assign(m.num_processes, 0);
    m.nlev.assign(m.num_processes, 0);
    m.mph.assign(m.num_processes, 1);
    m.level.assign(m.num_processes, std::vector<std::size_t>());
    m.svcrate.assign(m.num_processes, std::vector<T>());
    m.svcdown.assign(m.num_processes, std::vector<T>());
    m.is_open_proc.assign(m.num_processes, false);
    for (std::size_t ist : queue_stations)
        for (std::size_t r = 0; r < K; ++r) {
            const std::size_t p = m.process_map[ist][r];
            if (p == 0 || m.N[p - 1] != 0) continue;
            RcatComponent<T>& c = comp[p - 1];
            c.ist = ist;
            c.r = r;
            const std::pair<Matrix<T>, Matrix<T>> svc = proc_map(L, ist, r);
            c.Ds0 = svc.first;
            c.Ds1 = svc.second;
            arrival_map(L, c, source_stations);
            c.ns = c.Ds0.rows();
            c.na = c.Da0.rows();
            c.mph = c.na * c.ns;
            if (std::isfinite(njobs[r])) {
                c.nlev = static_cast<std::size_t>(njobs[r]) + 1;
            } else {
                c.nlev = max_states;  // open class: truncate
                m.is_open_proc[p - 1] = true;
            }
            // Service completion: level down, arrival phase untouched.
            c.Dsvc = kron(eye<T>(c.na), c.Ds1);
            c.N = c.nlev * c.mph;

            m.N[p - 1] = c.N;
            m.nlev[p - 1] = c.nlev;
            m.mph[p - 1] = c.mph;
            m.level[p - 1].assign(c.N, 0);
            for (std::size_t n = 0; n < c.nlev; ++n)
                for (std::size_t j = 0; j < c.mph; ++j) m.level[p - 1][n * c.mph + j] = n;
            m.svcdown[p - 1].assign(c.mph, zero);
            for (std::size_t i = 0; i < c.mph; ++i) {
                T acc = zero;
                for (std::size_t j = 0; j < c.mph; ++j) acc += c.Dsvc(i, j);
                m.svcdown[p - 1][i] = acc;
            }
            m.svcrate[p - 1].assign(c.N, zero);
            for (std::size_t n = 1; n < c.nlev; ++n)
                for (std::size_t j = 0; j < c.mph; ++j)
                    m.svcrate[p - 1][n * c.mph + j] = m.svcdown[p - 1][j];
        }

    for (std::size_t ist : queue_stations)
        for (std::size_t r = 0; r < K; ++r) {
            if (m.process_map[ist][r] == 0) continue;
            const bool is_removal =
                is_signal_class(L, r) && (signal_type_of(L, r) == SignalType::NEGATIVE ||
                                          signal_type_of(L, r) == SignalType::CATASTROPHE);
            const bool is_cat = is_removal && signal_type_of(L, r) == SignalType::CATASTROPHE;
            const bool has_dist =
                is_removal && r < L.signalremdist.size() && !L.signalremdist[r].empty();
            for (std::size_t jst : queue_stations)
                for (std::size_t s = 0; s < K; ++s) {
                    if (m.process_map[jst][s] == 0) continue;
                    if (ist == jst && r == s) continue;
                    const T pr = L.rt(ist * K + r, jst * K + s);
                    if (!(pr > zero)) continue;
                    RcatAction a;
                    a.from_station = ist;
                    a.from_class = r;
                    a.to_station = jst;
                    a.to_class = s;
                    a.prob = num_traits<T>::to_double(pr);
                    a.is_negative = is_removal;
                    a.is_catastrophe = is_cat;
                    a.has_removal_dist = has_dist;
                    a.removal_class = r;
                    m.actions.push_back(a);
                }
        }

    m.L.reserve(m.num_processes);
    for (std::size_t p = 0; p < m.num_processes; ++p)
        m.L.push_back(build_local_rates(L, comp[p], sink_nodes));

    const std::size_t A = m.actions.size();
    m.Aa.reserve(A);
    m.Pb.reserve(A);
    m.act.reserve(A);
    m.psv.reserve(A);
    for (const RcatAction& a : m.actions) {
        const std::size_t pa = m.process_map[a.from_station][a.from_class] - 1;
        const std::size_t pp = m.process_map[a.to_station][a.to_class] - 1;
        m.act.push_back(pa);
        m.psv.push_back(pp);

        const RcatComponent<T>& ca = comp[pa];
        const T prob = num_traits<T>::from_double(a.prob);
        Matrix<T> Am(ca.N, ca.N, zero);
        for (std::size_t n = 1; n < ca.nlev; ++n)
            add_block(Am, n, n - 1, ca.mph, ca.Dsvc, prob);
        // The boundary self-loop is physical only for a closed class, where the
        // top state is the real population bound. On an open class the top state
        // is the artefact of the maxStates truncation, and a self-loop there
        // would feed the truncation bias straight into the reversed rate. It is
        // written on the DIAGONAL so it stays inert in the generator while still
        // contributing the pi(i)/pi(i) = 1 ratio the INAP estimators read off.
        if (std::isfinite(njobs[a.from_class])) {
            const std::size_t off = blk(ca.nlev - 1, ca.mph);
            for (std::size_t i = 0; i < ca.mph; ++i)
                Am(off + i, off + i) += T(m.svcdown[pa][i] * prob);
        }
        m.Aa.push_back(Am);

        const RcatComponent<T>& cp = comp[pp];
        Matrix<T> Bm(cp.N, cp.N, zero);
        if (a.is_negative) {
            if (a.is_catastrophe) {
                for (std::size_t n = 0; n < cp.nlev; ++n) set_identity(Bm, n, 0, cp.mph);
            } else if (a.has_removal_dist) {
                add_batch_removal(Bm, L.signalremdist[a.removal_class], one, cp.nlev, cp.mph);
                set_identity(Bm, 0, 0, cp.mph);  // an empty queue absorbs the signal
            } else {
                set_identity(Bm, 0, 0, cp.mph);
                for (std::size_t n = 1; n < cp.nlev; ++n) set_identity(Bm, n, n - 1, cp.mph);
            }
        } else {
            // The phase is untouched: a job joining does not restart the server,
            // and the service phase frozen at level 0 is the one the last
            // completion left behind, which for a phase-type is already its
            // entry distribution.
            for (std::size_t n = 0; n + 1 < cp.nlev; ++n) set_identity(Bm, n, n + 1, cp.mph);
            set_identity(Bm, cp.nlev - 1, cp.nlev - 1, cp.mph);
        }
        m.Pb.push_back(Bm);
    }
    return m;
}

/** Assemble a component's generator from the local rates and the current x. */
template <class T>
Matrix<T> assemble_generator(const RcatModel<T>& m, const std::vector<T>& x, std::size_t k) {
    Matrix<T> Qk = m.L[k];
    for (std::size_t c = 0; c < m.actions.size(); ++c) {
        if (m.psv[c] == k) {
            for (std::size_t i = 0; i < Qk.rows(); ++i)
                for (std::size_t j = 0; j < Qk.cols(); ++j) Qk(i, j) += T(x[c] * m.Pb[c](i, j));
        } else if (m.act[c] == k) {
            for (std::size_t i = 0; i < Qk.rows(); ++i)
                for (std::size_t j = 0; j < Qk.cols(); ++j) Qk(i, j) += m.Aa[c](i, j);
        }
    }
    // ctmc_makeinfgen discards whatever sits on the diagonal and rebuilds it
    // from the row sums, which is why the reference's unscaled diagonal
    // corrections here are inert. See the header note.
    return mc::ctmc_makeinfgen(Qk);
}

/** A matrix in double, whatever arithmetic the model carries. */
template <class T>
Matrix<double> ag_as_double(const Matrix<T>& m) {
    Matrix<double> out(m.rows(), m.cols(), 0.0);
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j)
            out(i, j) = num_traits<T>::to_double(m(i, j));
    return out;
}

/**
 * Agent k's static half, as the ag-worker protocol carries it. Double only,
 * because the wire is JSON; ag_sweep_cluster refuses any other arithmetic before
 * this is reached.
 */
template <class T>
AgWirePayload wire_payload(const RcatModel<T>& m, std::size_t k) {
    AgWirePayload a;
    a.k = static_cast<int>(k);
    a.n = static_cast<int>(m.N[k]);
    a.mph = static_cast<int>(m.mph[k]);
    a.nlev = static_cast<int>(m.nlev[k]);
    a.level.assign(m.level[k].begin(), m.level[k].end());
    a.L = ag_triplets(ag_as_double(m.L[k]));
    for (std::size_t c = 0; c < m.actions.size(); ++c) {
        if (m.psv[c] == k) {
            a.passive_c.push_back(static_cast<int>(c));
            a.passive_m.push_back(ag_triplets(ag_as_double(m.Pb[c])));
        } else if (m.act[c] == k) {
            a.active_c.push_back(static_cast<int>(c));
            a.active_m.push_back(ag_triplets(ag_as_double(m.Aa[c])));
        }
    }
    return a;
}

/**
 * Port of `compute_equilibrium`: every component solved in isolation.
 *
 * Agent k reads the rest of the model only through the scalar reversed rates x
 * and writes only its own slot, so this is a fan-out and not a recurrence. That
 * is what lets @p exec evaluate the agents on a thread pool or on remote workers
 * and still walk the same iterates as the serial loop; any cross-agent read
 * added here would silently make those backends race.
 */
template <class T>
void compute_equilibrium(const RcatModel<T>& m, const std::vector<T>& x,
                         std::vector<std::vector<T>>& pi, std::vector<Matrix<T>>& Q,
                         const AgOptions* opt = nullptr, AgWorkerPool* pool = nullptr) {
    pi.assign(m.num_processes, std::vector<T>());
    Q.assign(m.num_processes, Matrix<T>());

    auto gen = [&](std::size_t k) { return assemble_generator(m, x, k); };
    auto sol = [&](const Matrix<T>& Qk, std::size_t k) {
        return solve_component(Qk, m.mph[k], m.nlev[k], m.level[k]);
    };

    const std::string mode = (opt == nullptr) ? std::string(exec_serial()) : opt->exec;
    if (exec_is_parallel(mode)) {
        ag_sweep_parallel<T>(m.num_processes, opt->nworkers, gen, sol, Q, pi);
    } else if (mode == exec_cluster()) {
        auto payload = [&](std::size_t k) { return wire_payload(m, k); };
        ag_sweep_cluster<T>(m.num_processes, *pool, x, gen, sol, payload, Q, pi);
    } else {
        ag_sweep_serial<T>(m.num_processes, gen, sol, Q, pi);
    }
}


/**
 * Port of `compute_equilibrium_qbd`: open components on their infinite state
 * space, closed ones as before.
 */
template <class T>
void compute_equilibrium_qbd(const RcatModel<T>& m, const std::vector<T>& x,
                             std::vector<std::vector<T>>& pi, std::vector<Matrix<T>>& Q,
                             std::vector<double>& rho_proc, std::vector<bool>& is_geom,
                             std::vector<QbdTail<T>>& geom_data) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    pi.assign(m.num_processes, std::vector<T>());
    Q.assign(m.num_processes, Matrix<T>());
    rho_proc.assign(m.num_processes, 0.0);
    is_geom.assign(m.num_processes, false);
    geom_data.assign(m.num_processes, QbdTail<T>());

    for (std::size_t k = 0; k < m.num_processes; ++k) {
        const std::size_t Nk = m.N[k];
        Matrix<T> Off = m.L[k];
        for (std::size_t i = 0; i < Nk; ++i) Off(i, i) = zero;
        for (std::size_t c = 0; c < m.actions.size(); ++c) {
            if (m.psv[c] == k) {
                for (std::size_t i = 0; i < Nk; ++i)
                    for (std::size_t j = 0; j < Nk; ++j) Off(i, j) += T(x[c] * m.Pb[c](i, j));
            } else if (m.act[c] == k) {
                for (std::size_t i = 0; i < Nk; ++i)
                    for (std::size_t j = 0; j < Nk; ++j) Off(i, j) += m.Aa[c](i, j);
            }
        }
        for (std::size_t i = 0; i < Nk; ++i) Off(i, i) = zero;
        Q[k] = mc::ctmc_makeinfgen(Off);

        const std::size_t mph = m.mph[k], nlev = m.nlev[k];
        bool solved_geom = false;
        if (m.is_open_proc[k] && nlev >= 5) {
            if (mph == 1) {
                // Read the homogeneous interior one level below the truncation,
                // so the reflecting boundary of Off does not contaminate it.
                const std::size_t s0 = Nk - 2;
                const double f = num_traits<T>::to_double(Off(s0, s0 + 1));
                const double b = num_traits<T>::to_double(Off(s0, s0 - 1));
                const double g0 = num_traits<T>::to_double(Off(s0, 0));
                // A jump to a strictly interior lower level is batch removal onto
                // a non-empty state, which is not a scalar QBD; defer to the
                // finite solve rather than fit a geometric that does not hold.
                double inter_down = 0.0;
                for (std::size_t j = 1; j + 3 < Nk; ++j)
                    inter_down += num_traits<T>::to_double(Off(s0, j));
                if (inter_down <= 1e-11 && f > 0.0) {
                    const double rho = qbd_scalar_rho(f, b, g0);
                    if (std::isfinite(rho) && rho > 0.0 && rho < 1.0 - 1e-12) {
                        rho_proc[k] = rho;
                        is_geom[k] = true;
                        pi[k].assign(Nk, zero);
                        const T rt = num_traits<T>::from_double(rho);
                        T acc = T(one - rt);
                        for (std::size_t n = 0; n < Nk; ++n) {
                            pi[k][n] = acc;
                            acc = T(acc * rt);
                        }
                        solved_geom = true;
                    }
                }
            } else if (is_block_tridiagonal(Q[k], m.level[k])) {
                // Read the homogeneous interior BLOCKS one level below the
                // truncation, for the same reason the scalar branch reads the
                // interior row there.
                const std::size_t s0 = nlev - 2;
                const Matrix<T> A0 = level_block(Q[k], s0, s0 + 1, mph);   // up
                const Matrix<T> A1 = level_block(Q[k], s0, s0, mph);       // local
                const Matrix<T> A2 = level_block(Q[k], s0, s0 - 1, mph);   // down
                bool any_up = false;
                for (std::size_t i = 0; i < mph && !any_up; ++i) {
                    T rs = zero;
                    for (std::size_t j = 0; j < mph; ++j) rs += A0(i, j);
                    if (rs > zero) any_up = true;
                }
                if (any_up) {
                    const QbdTail<T> g = qbd_matrix_tail(Q[k], A0, A1, A2, mph);
                    if (g.ok) {
                        is_geom[k] = true;
                        geom_data[k] = g;
                        pi[k] = qbd_tail_expand(g, nlev, mph);
                        solved_geom = true;
                    }
                }
            }
        }
        if (!solved_geom) pi[k] = solve_component(Q[k], mph, nlev, m.level[k]);
    }
}

/** The reference's block norm: the largest 1-norm change over the components. */
template <class T>
double block_norm(const std::vector<std::vector<T>>& a, const std::vector<std::vector<T>>& b) {
    double e = 0.0;
    for (std::size_t k = 0; k < a.size() && k < b.size(); ++k) {
        const std::size_t n = std::min(a[k].size(), b[k].size());
        double s = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            s += std::fabs(num_traits<T>::to_double(T(a[k][i] - b[k][i])));
        if (s > e) e = s;
    }
    return e;
}

/** Port of `rcat_metrics`. */
template <class T>
mva::MvaSolution<T> rcat_metrics(const qn::NetworkStruct<T>& L, const RcatModel<T>& m,
                                 const std::vector<std::vector<T>>& pi,
                                 const std::vector<double>& rho_proc,
                                 const std::vector<bool>& is_geom,
                                 const std::vector<QbdTail<T>>& geom_data) {
    const std::size_t M = L.nstations, K = L.nclasses;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    mva::MvaSolution<T> s;
    s.Q = Matrix<T>(M, K, zero);
    s.U = Matrix<T>(M, K, zero);
    s.R = Matrix<T>(M, K, zero);
    s.Tp = Matrix<T>(M, K, zero);
    s.C.assign(K, zero);
    s.X.assign(K, zero);

    for (std::size_t ist = 0; ist < M; ++ist)
        for (std::size_t r = 0; r < K; ++r) {
            const std::size_t p = m.process_map.empty() ? 0 : m.process_map[ist][r];
            if (p == 0 || p > pi.size() || pi[p - 1].empty()) continue;
            const bool disabled = L.disabled[ist][r];
            const std::size_t mph = m.mph[p - 1];
            if (!is_geom.empty() && p <= is_geom.size() && is_geom[p - 1]) {
                if (mph == 1) {
                    // The infinite geometric marginal pi_n = (1-rho) rho^n, whose
                    // moments are closed forms; using the truncated vector here
                    // would put the truncation back into the answer.
                    const T mu = disabled ? zero : L.rates(ist, r);
                    const T rho = num_traits<T>::from_double(rho_proc[p - 1]);
                    s.Q(ist, r) = T(rho / (one - rho));
                    s.U(ist, r) = rho;
                    if (mu > zero) s.Tp(ist, r) = T(mu * rho);
                } else {
                    // The matrix-geometric tail pi_(n+1) = pi_n R.
                    const QbdTail<T>& g = geom_data[p - 1];
                    T busy = zero, tput = zero;
                    for (std::size_t i = 0; i < mph; ++i) {
                        busy += g.busy[i];
                        tput += T(g.busy[i] * m.svcdown[p - 1][i]);
                    }
                    s.Q(ist, r) = g.qlen;
                    s.U(ist, r) = busy;
                    s.Tp(ist, r) = tput;
                }
            } else {
                const std::vector<T>& v = pi[p - 1];
                T q = zero, tput = zero;
                for (std::size_t n = 0; n < v.size(); ++n) {
                    q += T(num_traits<T>::from_int(static_cast<long>(m.level[p - 1][n])) * v[n]);
                    // The rate of service completions, i.e. the phase-dependent
                    // departure rate averaged over the marginal. With one phase
                    // this is the mean rate times P(N>0).
                    tput += T(m.svcrate[p - 1][n] * v[n]);
                }
                s.Q(ist, r) = q;
                T level0 = zero;
                for (std::size_t i = 0; i < mph; ++i) level0 += v[i];
                s.U(ist, r) = T(one - level0);
                s.Tp(ist, r) = tput;
            }
        }

    for (std::size_t ist = 0; ist < M; ++ist)
        for (std::size_t r = 0; r < K; ++r)
            if (s.Tp(ist, r) > zero) s.R(ist, r) = T(s.Q(ist, r) / s.Tp(ist, r));

    // A SOURCE'S THROUGHPUT IS ITS ARRIVAL RATE, and the loop above cannot fill
    // it: a Source has no queue and no service process, so the RCAT process map
    // covers no station of that type and every one of its cells stayed zero. The
    // reference reports lambda there, and -- far more than a cosmetic row -- the
    // whole open network's ARRIVAL RATES are derived from the throughput vector
    // by `sn_get_arvr_from_tput`, so a zero at the Source propagated to zero
    // arrivals everywhere downstream: `ag_tandem_open` reported ArvR 0 at Queue1
    // against a golden of 0.5, and `ag_gnetwork` lost its Negative class's rows
    // entirely. It is the same fact `s.X[r]` below already reads off the Source
    // for an open class, written into the table it belongs in.
    for (std::size_t ist = 0; ist < M; ++ist) {
        if (L.stations[ist].nodetype != qn::NodeType::Source) continue;
        for (std::size_t r = 0; r < K; ++r)
            s.Tp(ist, r) = L.disabled[ist][r] ? zero : L.rates(ist, r);
    }

    const std::vector<double> njobs = L.njobs();
    for (std::size_t r = 0; r < K; ++r) {
        if (!std::isfinite(njobs[r])) {
            for (std::size_t ist = 0; ist < M; ++ist) {
                // From the STATION, for the reason `build_rcat` reads it there.
                if (L.stations[ist].nodetype == qn::NodeType::Source) {
                    s.X[r] = L.disabled[ist][r] ? zero : L.rates(ist, r);
                    break;
                }
            }
            T c = zero;
            for (std::size_t ist = 0; ist < M; ++ist) c += s.R(ist, r);
            s.C[r] = c;
        } else {
            const std::size_t refst = L.classes[r].refstat;
            if (refst > 0 && refst <= M) {
                s.X[r] = s.Tp(refst - 1, r);
                if (s.X[r] > zero)
                    s.C[r] = T(num_traits<T>::from_double(njobs[r]) / s.X[r]);
            }
        }
    }
    return s;
}

}  // namespace ag_detail

/**
 * The process types the RCAT construction can give a phase dimension to.
 *
 * After `sn_nonmarkov_toph` (which the RCAT methods run with the phase-type
 * fit and no Det preservation) each of these holds a genuine (D0,D1) pair with
 * non-negative off-diagonal rates and a single arrival per epoch. The list is
 * an ALLOW-list on purpose: a process type nobody has checked against this
 * construction must be refused, not answered. Refused are the laws whose
 * matrices are not a generator (ME, RAP), those that are not time-homogeneous
 * (NHPP, MAPt, PHt), those that are not continuous-time (DMAP), and those that
 * arrive in batches (BMAP, MMAP), since a batch moves the level by more than
 * one.
 */
inline bool rcat_supports_process(lang::ProcessType t) {
    switch (t) {
        case lang::ProcessType::EXP:
        case lang::ProcessType::ERLANG:
        case lang::ProcessType::HYPEREXP:
        case lang::ProcessType::PH:
        case lang::ProcessType::APH:
        case lang::ProcessType::COXIAN:
        case lang::ProcessType::COX2:
        case lang::ProcessType::MAP:
        case lang::ProcessType::MMPP2:
        case lang::ProcessType::DET:
        case lang::ProcessType::UNIFORM:
        case lang::ProcessType::GAMMA:
        case lang::ProcessType::PARETO:
        case lang::ProcessType::WEIBULL:
        case lang::ProcessType::LOGNORMAL:
        case lang::ProcessType::REPLAYER:
        case lang::ProcessType::IMMEDIATE:
        case lang::ProcessType::DISABLED:
            return true;
        default:
            return false;
    }
}

/** What the RCAT analyzer returns beyond the metrics. */
template <class T>
struct AgResult {
    mva::MvaSolution<T> sol;
    std::string actualmethod;
    /**
     * The RCAT product-form residual of Remark 2, max_l ||pi (x_l I - T_l)||_2,
     * computed by 'inapinf' only. Zero exactly when the reversed rates came out
     * state-independent, i.e. when the product form the method assumes is real;
     * anything else is the size of the modelling error, not of a numerical one.
     */
    double rcat_residual = 0.0;
};

/**
 * Port of `solver_ag.m`.
 *
 * `max_states` is the reference's `options.config.maxStates`, the truncation of
 * every OPEN component. It defaults to AgOptions::max_states, which carries the
 * reference's 100; the explicit parameter overrides it for a caller that has no
 * options object to hand.
 */
template <class T>
AgResult<T> solver_ag(const qn::NetworkStruct<T>& L, const AgOptions& opt,
                          std::size_t max_states = 0) {
    if (max_states == 0) max_states = opt.max_states;
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_ag: the RCAT analyzers stop on a tolerance and the matrix-geometric "
            "variant takes a square root, so they need transcendental arithmetic; rerun this "
            "model with --arith double or --arith real");
    } else {
    using namespace ag_detail;
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = L.nstations, K = L.nclasses;

    std::string method = opt.method;
    if (method == "default") method = "inap";
    if (method == "exact") {
        // The reference warns and falls back rather than erroring, because
        // autocat moved out of the tree; solver_mam_autocat.h records that.
        method = "inap";
    }
    if (method != "inap" && method != "inapplus" && method != "inapinf")
        throw UnsupportedError("solver_ag: unknown method '" + opt.method + "'");

    AgResult<T> out;
    out.actualmethod = method;

    if (opt.exec == exec_cluster() && method == "inapinf") {
        // The remote worker implements the FINITE agent solve. 'inapinf' replaces
        // it with the matrix-geometric treatment of an open agent -- Neuts' R
        // matrix and the scalar-tail detection that precedes it -- which the
        // worker does not carry, and answering with the finite solve instead
        // would silently change the method.
        throw UnsupportedError(
            "solver_ag: the 'cluster' execution backend does not carry the 'inapinf' agent "
            "solve (the matrix-geometric tail of an open agent runs on the coordinator "
            "only); use exec 'serial' or 'parallel' with 'inapinf', or method 'inap'/"
            "'inapplus' with 'cluster'");
    }
    if (opt.exec != exec_serial() && !exec_is_parallel(opt.exec)
            && opt.exec != exec_cluster()) {
        // 'threads' was this backend's name until 2026-08-19. Naming the rename
        // costs one branch and saves a caller with an old script from reading
        // "unknown backend" about a backend that still exists.
        const std::string hint = (opt.exec == "threads")
            ? "; 'threads' was renamed to 'parallel' (alias 'para')" : "";
        throw InputError("solver_ag: unknown execution backend '" + opt.exec +
                         "'; use 'serial', 'parallel' (alias 'para') or 'cluster'" + hint);
    }

    std::unique_ptr<AgWorkerPool> pool;
    if (opt.exec == exec_cluster()) {
        pool.reset(new AgWorkerPool(opt.endpoints, opt.worker_timeout));
    }

    const RcatModel<T> m = build_rcat(L, max_states);

    if (m.num_processes == 0) {
        // The reference warns and returns zeros rather than erroring: a model
        // with no serving station is degenerate, not invalid.
        out.sol.Q = Matrix<T>(M, K, zero);
        out.sol.U = Matrix<T>(M, K, zero);
        out.sol.R = Matrix<T>(M, K, zero);
        out.sol.Tp = Matrix<T>(M, K, zero);
        out.sol.C.assign(K, zero);
        out.sol.X.assign(K, zero);
        out.sol.method = opt.method;
        out.sol.iter = 0;
        return out;
    }

    const std::size_t A = m.actions.size();
    std::vector<std::vector<T>> pi;
    std::vector<Matrix<T>> Q;

    if (A == 0) {
        // No synchronizing action, so there is no fixed point to run: every
        // component is already closed under its local rates alone. This is the
        // single-queue G-network shape (Source -> Queue -> Sink).
        // Solved through the same dispatcher the fixed point uses: this branch
        // carries a whole M/PH/1 on its own, whose marginal spans tens of orders
        // of magnitude over the truncation, and the level recursions are stable
        // there where a null-space solve is not.
        pi.assign(m.num_processes, std::vector<T>());
        for (std::size_t k = 0; k < m.num_processes; ++k) {
            const Matrix<T> Qk = mc::ctmc_makeinfgen(m.L[k]);
            pi[k] = solve_component(Qk, m.mph[k], m.nlev[k], m.level[k]);
        }
        out.sol = rcat_metrics(L, m, pi, std::vector<double>(), std::vector<bool>(),
                               std::vector<QbdTail<T>>());
        out.sol.method = opt.method;
        out.sol.iter = 0;
        return out;
    }

    // A component reaching beyond the NEIGHBOURING LEVEL is not a birth-death
    // chain, and the mean-of-ratios estimator is meaningless there, so INAP
    // switches it to the rate-conserving one that INAP+ uses everywhere. The
    // within-level phase transitions of a PH sit far off the diagonal and are
    // NOT such a departure, which is why the test is on the level index.
    //
    // A PHASE-EXPANDED component takes the same estimator, for the same reason.
    // On a birth-death chain every state-wise reversed rate equals lambda, so
    // their mean is exact; with a phase block per level they do not, the deep
    // truncation levels dominate the unweighted mean, and the mean-of-ratios
    // overestimates the departure rate exactly as it does on a catastrophe
    // (measured on a tandem with Erlang(2) service at Q1: the reversed rate came
    // out 1.27 against the exact 0.5, so flow was not conserved).
    std::vector<bool> not_birth_death(m.num_processes, false);
    for (std::size_t k = 0; k < m.num_processes; ++k) {
        if (m.mph[k] > 1) not_birth_death[k] = true;
        for (std::size_t n = 0; n < m.N[k]; ++n)
            for (std::size_t j = 0; j < m.N[k]; ++j) {
                const std::size_t ln = m.level[k][n], lj = m.level[k][j];
                const std::size_t d = (ln > lj) ? (ln - lj) : (lj - ln);
                if (d > 1 && m.L[k](n, j) > zero) not_birth_death[k] = true;
            }
    }

    // Columns of each active matrix that carry any rate. Aa does not depend on
    // x, so this is fixed for the whole fixed point.
    std::vector<std::vector<std::size_t>> active_cols(A);
    for (std::size_t a = 0; a < A; ++a)
        for (std::size_t j = 0; j < m.Aa[a].cols(); ++j) {
            T colsum = zero;
            for (std::size_t i = 0; i < m.Aa[a].rows(); ++i) colsum += m.Aa[a](i, j);
            if (colsum > zero) active_cols[a].push_back(j);
        }

    // Deterministic initial guess, so the answer does not depend on a seed.
    std::vector<T> x(A, zero);
    for (std::size_t a = 0; a < A; ++a)
        x[a] = num_traits<T>::from_double(static_cast<double>(a + 1) /
                                          static_cast<double>(A + 1));

    std::vector<double> rho_proc;
    std::vector<bool> is_geom;
    std::vector<QbdTail<T>> geom_data;
    const bool qbd = (method == "inapinf");
    if (qbd) compute_equilibrium_qbd(m, x, pi, Q, rho_proc, is_geom, geom_data);
    else compute_equilibrium(m, x, pi, Q, &opt, pool.get());

    std::vector<std::vector<T>> a_rowsum(A);
    for (std::size_t a = 0; a < A; ++a) {
        a_rowsum[a].assign(m.Aa[a].rows(), zero);
        for (std::size_t i = 0; i < m.Aa[a].rows(); ++i)
            for (std::size_t j = 0; j < m.Aa[a].cols(); ++j) a_rowsum[a][i] += m.Aa[a](i, j);
    }

    // The fixed point is written out rather than driven by da_fpi: the
    // reference installs a custom `da_norm` (the per-component 1-norm above),
    // and the C++ da_fpi hard-codes the max-abs norm of a flat vector, which is
    // a strictly weaker stopping test and would stop earlier.
    const double tol = opt.tol;
    std::size_t iter = 0;
    bool converged = false;
    for (iter = 1; iter <= static_cast<std::size_t>(opt.iter_max); ++iter) {
        const std::vector<std::vector<T>> pi_ref = pi;
        for (std::size_t a = 0; a < A; ++a) {
            const std::size_t k = m.act[a];
            if (qbd) {
                if (is_geom[k]) {
                    if (m.mph[k] == 1) {
                        // On a geometric tail the active label fires only from an
                        // occupied state, so x = (per-state rate) * P(occupied).
                        const std::size_t idx = std::min<std::size_t>(1, m.N[k] - 1);
                        x[a] = T(a_rowsum[a][idx] * num_traits<T>::from_double(rho_proc[k]));
                    } else {
                        // Matrix-geometric tail: sum_{n>=1} pi_n = pi_1 (I-R)^-1,
                        // and the active label has the same row sums at every
                        // busy level.
                        T acc = zero;
                        for (std::size_t i = 0; i < m.mph[k]; ++i)
                            acc += T(geom_data[k].busy[i] * a_rowsum[a][m.mph[k] + i]);
                        x[a] = acc;
                    }
                } else {
                    T acc = zero;
                    for (std::size_t i = 0; i < pi[k].size() && i < a_rowsum[a].size(); ++i)
                        acc += T(pi[k][i] * a_rowsum[a][i]);
                    x[a] = acc;
                }
            } else {
                // pi Aa, the row vector both estimators are read off.
                const std::vector<T> v = vecmul(pi[k], m.Aa[a]);
                if (method == "inapplus" || not_birth_death[k]) {
                    // The departure rate of the active component.
                    T sum = zero;
                    for (std::size_t j = 0; j < v.size(); ++j) sum += v[j];
                    if (sum > zero) x[a] = sum;
                } else {
                    // The mean over the support of the STATE-WISE reversed rate
                    // (pi Aa)_j / pi_j, which RCAT requires to be independent of
                    // j. On a birth-death component every column of Aa holds one
                    // entry, so this is the reference's entrywise mean of
                    // Aa(i,j) pi(i) / pi(j) term for term.
                    T sum = zero;
                    std::size_t cnt = 0;
                    for (std::size_t j : active_cols[a])
                        if (pi[k][j] > zero && v[j] > zero) {
                            sum += T(v[j] / pi[k][j]);
                            ++cnt;
                        }
                    if (cnt > 0) x[a] = T(sum / num_traits<T>::from_int(static_cast<long>(cnt)));
                }
            }
        }
        if (qbd) compute_equilibrium_qbd(m, x, pi, Q, rho_proc, is_geom, geom_data);
        else compute_equilibrium(m, x, pi, Q, &opt, pool.get());

        if (block_norm<T>(pi, pi_ref) < tol) {
            converged = true;
            break;
        }
    }
    // The reference's legacy while-loop leaves the counter one past the cap
    // when it never converged, and callers read that as the "did not converge"
    // marker; reproduced so the iteration counts agree.
    if (!converged) iter = static_cast<std::size_t>(opt.iter_max) + 1;

    if (qbd) {
        double res = 0.0;
        for (std::size_t a = 0; a < A; ++a) {
            const std::size_t k = m.act[a];
            double s = 0.0;
            for (std::size_t j = 0; j < m.N[k]; ++j) {
                T acc = T(x[a] * pi[k][j]);
                for (std::size_t i = 0; i < m.N[k]; ++i) acc -= T(pi[k][i] * m.Aa[a](i, j));
                const double d = num_traits<T>::to_double(acc);
                s += d * d;
            }
            res = std::max(res, std::sqrt(s));
        }
        out.rcat_residual = res;
    }

    out.sol = rcat_metrics(L, m, pi, rho_proc, is_geom, geom_data);
    out.sol.method = opt.method;
    out.sol.iter = static_cast<int>(iter);
    return out;
    }  // if constexpr has_transcendental
}

}  // namespace ag
}  // namespace line

#endif  // LINE_SOLVERS_AG_SOLVER_AG_H
