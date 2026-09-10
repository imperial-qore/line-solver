/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_BMAP_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_BMAP_H

/**
 * The batch-arrival and batch-service queues of the MAM solver, and the two
 * finite-capacity helpers `solver_mam_basic.m` shares with them:
 * `solver_mam_bmap_map_1.m`, `solver_mam_map_bmap_1.m`, `mam_detect_mmck.m`
 * and `mam_truncate_renorm.m`.
 *
 * NEITHER SOLVER IS A CONVERSION. Both take the batch process AS a batch
 * process and build the exact level-transition blocks of the resulting chain;
 * neither one replaces a BMAP by a MAP of the same rate, and neither is
 * expressible as the other. What differs between them is the SHAPE of the
 * chain, and that shape is forced by which side carries the batches:
 *
 *   BMAP/MAP/1   the level rises by the batch size and falls by exactly one,
 *                so the chain is M/G/1-TYPE (skip-free to the left)
 *   MAP/BMAP/1   the level rises by exactly one and falls by the batch size,
 *                so the chain is GI/M/1-TYPE (skip-free to the right)
 *
 * Both block assemblies are EXACT: no batch is ever resolved into independent
 * single arrivals, which is the approximation that would destroy the whole
 * point. A batch of k is one epoch at which k jobs appear together, and the
 * queue-length law it induces is not the law induced by k Poisson epochs of
 * the same total rate -- the batch correlates the arrivals perfectly, so the
 * second moment of the queue length is strictly larger. That is why `A_{k+1}`
 * moves the level by k in one step rather than being folded into `A_2`.
 *
 * WHERE AN APPROXIMATION DOES ENTER, it is at the BOUNDARY of the GI/M/1-type
 * chain, and it is a clipping rather than a loss. A batch service of size k at
 * a level j < k cannot take the level to j - k, so `solver_mam_map_bmap_1.m`
 * routes every batch of size k >= j from level j to level ZERO:
 *
 *     B_{j+1} = sum_{k >= j} I (x) D_k
 *
 * and at level 0 it folds the whole service mass back onto the local block.
 * THE TAIL IS LUMPED, NOT DISCARDED, which is the same convention
 * `State.signalBatchPMF` uses for a negative signal (see `signal_batch_pmf` in
 * `lang/qn/state_events.h`: the whole tail P(B >= n) is assigned to "remove
 * all n"). The invariant it buys is that the boundary rows still sum to zero,
 * so no probability is created or destroyed by the clipping; the price is that
 * the boundary over-reports emptying events relative to a chain that could
 * represent the missing customers.
 *
 * `mam_truncate_renorm` uses the OPPOSITE convention on purpose, and the
 * contrast is the reason both are in this header. It truncates the marginal
 * queue length at the buffer capacity and RENORMALIZES, so the tail above capK
 * is deleted and its mass is redistributed over levels 0..capK in proportion,
 * rather than piled onto level capK. For an M/M/1 input that is not an
 * approximation at all -- the M/M/1/K law IS the truncated renormalized
 * geometric -- which is what makes it the right convention for a buffer, where
 * a blocked arrival leaves the system rather than joining at the top. Lumping
 * would instead report a boundary mass that no finite-buffer queue has.
 *
 * THE MEAN MEASURES ARE THE ETAQA SOLVE, and both reference files hand it to
 * third-party MAMSolver: `MG1_G_ETAQA`, `MG1_pi_ETAQA` and `MG1_qlen_ETAQA` on
 * the M/G/1-type side, `GIM1_R_ETAQA`, `GIM1_pi_ETAQA` and `GIM1_qlen_ETAQA` on
 * the GI/M/1-type one, routing in turn into `MG1_CR` (Bini-Meini cyclic
 * reduction with the shift technique and FFT polynomial products) and `GIM1_R`
 * (the Bright/Ramaswami dual plus functional iterations). All of it is now
 * ported, under `lib/smc/mg1.h` and `lib/smc/etaqa.h`, so these two entry
 * points answer instead of refusing. Read those headers before touching the
 * numbers: the ported code reproduces several reference defects verbatim and
 * says which.
 *
 * ARITHMETIC. The block assemblies are Kronecker products and one stationary
 * solve, so they are exact at Rational and are NOT gated. The ETAQA solve on
 * top of them is DOUBLE ONLY -- LAPACK eigenvalues in the caudal and decay
 * bisections, a complex FFT in cyclic reduction, an SVD in the rank test that
 * picks the redundant column -- so `solver_mam_bmap_map_1` and
 * `solver_mam_map_bmap_1` refuse at any other arithmetic rather than
 * down-convert behind the caller's back. `mam_truncate_renorm` is gated too,
 * because MMAP[K]/PH[K]/1 runs the ADDA doubling iteration.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <type_traits>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmapph1fcfs.h"
#include "line/api/mam/qbd_r.h"
#include "line/lang/qn/network_struct.h"
#include "line/lib/smc/etaqa.h"
#include "line/solvers/mam/solver_mam_basic.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace bmap_detail {

using qbd_detail::madd;

/** D1 + ... + DK, the epoch-marginal counterpart of a MAP's D1. */
template <class T>
Matrix<T> batch_d1(const std::vector<Matrix<T>>& D) {
    Matrix<T> tot(D[0].rows(), D[0].cols(), num_traits<T>::from_int(0));
    for (std::size_t k = 1; k < D.size(); ++k) tot = madd(tot, D[k]);
    return tot;
}

/** D0 + D1 + ... + DK, the generator of the batch process's phase chain. */
template <class T>
Matrix<T> batch_infgen(const std::vector<Matrix<T>>& D) {
    return madd(D[0], batch_d1(D));
}

/**
 * Customers per unit time, sum_k k theta D_k e.
 *
 * This is NOT the epoch rate theta (sum_k D_k) e: a batch of k counts k times.
 * Confusing the two is the error that makes a BMAP look like a MAP of the same
 * epoch rate, and it is the reason the two rates are computed by different
 * functions here rather than by one with a flag.
 */
template <class T>
T batch_customer_rate(const std::vector<Matrix<T>>& D) {
    const T zero = num_traits<T>::from_int(0);
    Map<T> phase;
    phase.D0 = D[0];
    phase.D1 = batch_d1(D);
    const std::vector<T> theta = map_prob(phase);
    T rate = zero;
    for (std::size_t k = 1; k < D.size(); ++k) {
        const std::vector<T> tD = vecmul(theta, D[k]);
        T s = zero;
        for (const T& v : tD) s += v;
        rate += T(num_traits<T>::from_int(static_cast<long>(k)) * s);
    }
    return rate;
}

/** Every matrix square, of one common order, and at least {D0, D1}. */
template <class T>
std::size_t check_batch_shape(const std::vector<Matrix<T>>& D, const std::string& who,
                              const std::string& what) {
    if (D.size() < 2)
        throw InputError(who + ": the " + what +
                         " must be given as {D0, D1, ..., DK} with at least D0 and D1");
    const std::size_t n = D[0].rows();
    for (std::size_t k = 0; k < D.size(); ++k)
        if (D[k].rows() != n || D[k].cols() != n)
            throw InputError(who + ": all " + what + " matrices must be " + std::to_string(n) +
                             "x" + std::to_string(n) + ", but D{" + std::to_string(k) + "} is " +
                             std::to_string(D[k].rows()) + "x" + std::to_string(D[k].cols()));
    return n;
}

/** The reference's `max(abs(sum(Q,2))) > 1e-10` generator test. */
template <class T>
void check_zero_rowsums(const Matrix<T>& Q, const std::string& msg) {
    for (std::size_t i = 0; i < Q.rows(); ++i) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < Q.cols(); ++j) s += Q(i, j);
        if (std::fabs(num_traits<T>::to_double(s)) > 1e-10) throw InputError(msg);
    }
}

}  // namespace bmap_detail

// ---------------------------------------------------------------------------
// mam_detect_mmck
// ---------------------------------------------------------------------------

/** What `mam_detect_mmck` returns; muRate is meaningful only when isMmck. */
template <class T>
struct MmckDetection {
    bool isMmck = false;
    T muRate = num_traits<T>::from_int(0);
};

/**
 * Port of `mam_detect_mmck.m`: is the exact M/M/c/K closed form legitimate at
 * this station?
 *
 * All three conditions are about losing nothing, not about convenience. A
 * multi-phase arrival MMAP is not Poisson, so the M/M/c/K birth-death chain
 * would answer a different arrival process; a non-exponential service breaks
 * the same chain's death rates; and per-class rates that differ leave the
 * aggregate service non-exponential even when each class is. A class the
 * station never serves is skipped rather than failing the test, because the
 * reference reads its NaN rate as "no inflow" -- the C++ `disabled` flag is
 * that same sentinel.
 *
 * @param ist 1-based station index
 * @param arv the assembled arrival stream at that station
 * @param L the refreshed struct whose station ist is being tested
 */
template <class T>
MmckDetection<T> mam_detect_mmck(const qn::NetworkStruct<T>& L, std::size_t ist,
                                 const Mmap<T>& arv) {
    MmckDetection<T> out;
    if (ist == 0 || ist > L.nstations)
        throw InputError("mam_detect_mmck: station index " + std::to_string(ist) +
                         " is out of range");
    const std::size_t i0 = ist - 1;
    if (arv.order() != 1) return out;

    bool any = false;
    double lo = 0.0, hi = 0.0;
    for (std::size_t r = 0; r < L.nclasses; ++r) {
        if (L.disabled[i0][r]) continue;
        if (L.service[i0][r].type != lang::ProcessType::EXP) return out;
        const double v = num_traits<T>::to_double(L.rates(i0, r));
        if (!(v > 0.0)) continue;
        if (!any) {
            lo = hi = v;
            out.muRate = L.rates(i0, r);
            any = true;
        } else {
            lo = std::min(lo, v);
            hi = std::max(hi, v);
        }
    }
    if (!any) return out;
    if (hi - lo > 1e-9 * std::max(1.0, hi)) return out;
    out.isMmck = true;
    return out;
}

// ---------------------------------------------------------------------------
// mam_truncate_renorm
// ---------------------------------------------------------------------------

/**
 * Port of `mam_truncate_renorm.m`: the finite-buffer marginal of an
 * MMAP[K]/PH[K]/1 FCFS queue, by truncation and renormalization.
 *
 * The body is `solver_mam_basic.h`'s `basic_detail::truncate_renorm`, which is
 * the same function under the name the analyzer that first needed it gave it.
 * It is re-exposed here under the REFERENCE's name so that a caller reading
 * `mam_truncate_renorm.m` finds it, and so that the tail convention documented
 * at the top of this file has one place to be documented.
 *
 * Multi-class input is aggregated first: `ncDistr` returns the per-class
 * marginal P(N_k = n) and the truncation needs the joint P(N_total = n).
 */
template <class T>
basic_detail::TruncRenorm<T> mam_truncate_renorm(const Mmap<T>& arv,
                                                 const std::vector<PhService<T>>& svc,
                                                 std::size_t capK) {
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)arv;
        (void)svc;
        (void)capK;
        throw UnsupportedError(
            "mam_truncate_renorm: the infinite-buffer marginal comes from MMAP[K]/PH[K]/1 FCFS, "
            "whose ADDA doubling iteration terminates on a tolerance; rerun this model with "
            "--arith double or --arith real");
    } else {
        if (capK < 1) throw InputError("mam_truncate_renorm: the buffer capacity must be positive");
        return basic_detail::truncate_renorm(arv, svc, capK);
    }
}

// ---------------------------------------------------------------------------
// BMAP/MAP/1, an M/G/1-type chain
// ---------------------------------------------------------------------------

/**
 * The level blocks of `solver_mam_bmap_map_1.m`.
 *
 * The phase is the pair (arrival phase, service phase) with the ARRIVAL phase
 * major, so every block is a Kronecker product in that order.
 */
template <class T>
struct BmapMap1Blocks {
    std::size_t ma = 0;  ///< BMAP phases
    std::size_t ms = 0;  ///< service MAP phases
    std::size_t m = 0;   ///< ma * ms
    std::size_t K = 0;   ///< largest batch size

    Matrix<T> A0;                ///< level -1: a service completion
    Matrix<T> A1;                ///< level 0: phase changes only
    std::vector<Matrix<T>> Aup;  ///< Aup[k-1]: level +k, a batch of k

    Matrix<T> B0;                ///< level 0 local block, at the empty queue
    std::vector<Matrix<T>> Bup;  ///< Bup[k-1]: level 0 -> +k

    T lambda = num_traits<T>::from_int(0);  ///< customers per unit time
    T mu = num_traits<T>::from_int(0);      ///< service completions per unit time
    T rho = num_traits<T>::from_int(0);
    /** The reference WARNS rather than errors when rho >= 1; recorded, not thrown. */
    bool stable = true;
};

/**
 * Port of the block assembly and the stability test of
 * `solver_mam_bmap_map_1.m`.
 *
 * B0 IS NOT `qbd_mapmap1_blocks`'s Lbar. The reference adds the service
 * completion back onto the level-zero local block, `kron(I, S0 + S1)`, so the
 * service phase process keeps running while the queue is empty and the
 * completion it would have made is absorbed as a self-loop. That is a
 * modelling choice about an idle server, not an oversight, and it makes
 * B0 = A1 + A0 exactly. `qbd_mapmap1.h` instead stops the service process at
 * level zero with `kron(D0, I)`. The two chains differ, and both are here
 * under their own names.
 */
template <class T>
BmapMap1Blocks<T> solver_mam_bmap_map_1_blocks(const std::vector<Matrix<T>>& D,
                                               const Map<T>& service) {
    using namespace bmap_detail;
    BmapMap1Blocks<T> b;
    b.ma = check_batch_shape(D, "solver_mam_bmap_map_1", "BMAP");
    b.K = D.size() - 1;
    b.ms = service.D0.rows();
    if (service.D0.cols() != b.ms || service.D1.rows() != b.ms || service.D1.cols() != b.ms)
        throw InputError("solver_mam_bmap_map_1: the service MAP matrices must be " +
                         std::to_string(b.ms) + "x" + std::to_string(b.ms));
    b.m = b.ma * b.ms;

    const Matrix<T> Ia = eye<T>(b.ma), Is = eye<T>(b.ms);
    b.A0 = kron(Ia, service.D1);
    b.A1 = madd(kron(D[0], Is), kron(Ia, service.D0));
    for (std::size_t k = 1; k <= b.K; ++k) b.Aup.push_back(kron(D[k], Is));

    b.B0 = madd(kron(D[0], Is), kron(Ia, madd(service.D0, service.D1)));
    b.Bup = b.Aup;

    b.lambda = batch_customer_rate(D);
    b.mu = map_lambda(service);
    const T zero = num_traits<T>::from_int(0);
    b.rho = (b.mu > zero) ? T(b.lambda / b.mu) : zero;
    b.stable = num_traits<T>::to_double(b.rho) < 1.0;
    return b;
}

// ---------------------------------------------------------------------------
// MAP/BMAP/1, a GI/M/1-type chain
// ---------------------------------------------------------------------------

/** The level blocks of `solver_mam_map_bmap_1.m`. */
template <class T>
struct MapBmap1Blocks {
    std::size_t ma = 0;  ///< arrival MAP phases
    std::size_t ms = 0;  ///< service BMAP phases
    std::size_t m = 0;
    std::size_t K = 0;  ///< largest service batch

    Matrix<T> A0;                  ///< level +1: an arrival
    Matrix<T> A1;                  ///< level 0: phase changes only
    std::vector<Matrix<T>> Adown;  ///< Adown[k-1]: level -k, a batch service of k

    Matrix<T> B1;                 ///< level 0 local block, service folded back
    std::vector<Matrix<T>> Bto0;  ///< Bto0[j-1]: level j -> level 0, j = 1..K

    T lambda = num_traits<T>::from_int(0);
    T mu = num_traits<T>::from_int(0);  ///< customers served per unit time
    T rho = num_traits<T>::from_int(0);
    bool stable = true;
};

/**
 * Port of the block assembly, the generator validation and the stability test
 * of `solver_mam_map_bmap_1.m`.
 *
 * THE BOUNDARY IS WHERE THE CLIPPING LIVES. `Bto0[j-1]` collects every batch
 * of size k >= j, so an oversized batch empties the queue instead of driving
 * the level negative, and `B1` collects the whole service mass at level zero
 * as a self-loop. Nothing is dropped: for every boundary level j the total
 * outflow A0 + A1 + sum_{k<j} Adown[k-1] + Bto0[j-1] is again
 * kron(C0+C1, I) + kron(I, sum_k D_k), whose rows sum to zero. That identity
 * is the whole justification for the convention and is asserted in the tests.
 *
 * ONE REFERENCE BRANCH IS UNREACHABLE HERE and is therefore not transcribed:
 * `K = bmapSvc.getNumberOfPhases() - 1` reads the PHASE count where the batch
 * count is meant, so the BMAP-object path mis-sizes D for any process whose
 * order differs from its largest batch. This port takes the matrices directly,
 * which is the reference's own cell-array path and the correct one.
 */
template <class T>
MapBmap1Blocks<T> solver_mam_map_bmap_1_blocks(const Map<T>& arrival,
                                               const std::vector<Matrix<T>>& D) {
    using namespace bmap_detail;
    MapBmap1Blocks<T> b;
    b.ma = arrival.D0.rows();
    if (arrival.D0.cols() != b.ma || arrival.D1.rows() != b.ma || arrival.D1.cols() != b.ma)
        throw InputError("solver_mam_map_bmap_1: the arrival MAP matrices must be " +
                         std::to_string(b.ma) + "x" + std::to_string(b.ma));
    b.ms = check_batch_shape(D, "solver_mam_map_bmap_1", "service BMAP");
    b.K = D.size() - 1;
    b.m = b.ma * b.ms;

    check_zero_rowsums(madd(arrival.D0, arrival.D1),
                       "solver_mam_map_bmap_1: MAP matrices C0 + C1 must have zero row sums");
    const Matrix<T> Dtot = batch_infgen(D);
    check_zero_rowsums(Dtot,
                       "solver_mam_map_bmap_1: BMAP matrices D0 + D1 + ... + DK must have zero "
                       "row sums");

    const Matrix<T> Ia = eye<T>(b.ma), Is = eye<T>(b.ms);
    b.A0 = kron(arrival.D1, Is);
    b.A1 = madd(kron(arrival.D0, Is), kron(Ia, D[0]));
    for (std::size_t k = 1; k <= b.K; ++k) b.Adown.push_back(kron(Ia, D[k]));

    b.B1 = b.A1;
    for (std::size_t k = 1; k <= b.K; ++k) b.B1 = madd(b.B1, kron(Ia, D[k]));
    for (std::size_t j = 1; j <= b.K; ++j) {
        Matrix<T> Bj(b.m, b.m, num_traits<T>::from_int(0));
        for (std::size_t k = j; k <= b.K; ++k) Bj = madd(Bj, kron(Ia, D[k]));
        b.Bto0.push_back(Bj);
    }

    b.lambda = map_lambda(arrival);
    b.mu = batch_customer_rate(D);
    const T zero = num_traits<T>::from_int(0);
    b.rho = (b.mu > zero) ? T(b.lambda / b.mu) : zero;
    b.stable = num_traits<T>::to_double(b.rho) < 1.0;
    return b;
}

// ---------------------------------------------------------------------------
// The mean measures, which are the third-party solve
// ---------------------------------------------------------------------------

/** What both reference files return once the ETAQA solve has run. */
template <class T>
struct BmapQueueResult {
    T QN = num_traits<T>::from_int(0), UN = num_traits<T>::from_int(0);
    T RN = num_traits<T>::from_int(0), TN = num_traits<T>::from_int(0);
    Matrix<T> piAgg;  ///< the ETAQA-aggregated stationary vector
    /** G for the M/G/1-type solve, R for the GI/M/1-type one. */
    Matrix<T> fund;
    /** Moments 1..nMoments of the queue length; `QN` is the first of them. */
    std::vector<T> qlenMoments;
};

namespace bmap_detail {

/** The ETAQA solve is double only; say so instead of down-converting. */
template <class T>
void require_double_arith(const std::string& who) {
    if (!std::is_same<T, double>::value)
        throw UnsupportedError(
            who +
            ": the ETAQA mean measures run in double precision only. The solve bisects on a "
            "Perron-Frobenius eigenvalue (LAPACK), evaluates the cyclic reduction at complex "
            "roots of unity (FFT) and picks the redundant balance equation by a numerical rank "
            "test (SVD), none of which this tree provides at exact or multiprecision arithmetic. "
            "The level blocks, the rates and the stability test ARE available at every "
            "arithmetic from " +
            who + "_blocks; rerun the measures with --arith double");
}

/** Copy of a templated matrix as doubles, for the third-party solve. */
template <class T>
Matrix<double> as_double(const Matrix<T>& A) {
    Matrix<double> out(A.rows(), A.cols(), 0.0);
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) out(i, j) = num_traits<T>::to_double(A(i, j));
    return out;
}

/** `[M0 M1 ... Mk]`, the wide block sequence MAMSolver's M/G/1 side takes. */
inline Matrix<double> hstack(const std::vector<Matrix<double>>& blk) {
    return smc::hcat(blk);
}

/** `[M0; M1; ...; Mk]`, the stacked sequence its GI/M/1 side takes. */
inline Matrix<double> vstack(const std::vector<Matrix<double>>& blk) {
    return smc::vcat(blk);
}

}  // namespace bmap_detail

/**
 * Port of `solver_mam_bmap_map_1.m`, mean measures included.
 *
 * The blocks and the stability test are assembled first, so a malformed input
 * is reported as such before any numerics run. The chain is then handed to
 * ETAQA exactly as the reference hands it: `A = [A0 A1 A2 ... A_{K+1}]` for the
 * repetitive levels and `B = [B0 B1 ... BK]` for the boundary, G from
 * `MG1_G_ETAQA`, the three aggregates from `MG1_pi_ETAQA`, and the moments from
 * `MG1_qlen_ETAQA`. `nMoments` matches the reference's default of 3; the mean
 * queue length is the first of them and the response time follows by Little.
 */
template <class T>
BmapQueueResult<T> solver_mam_bmap_map_1(const std::vector<Matrix<T>>& D, const Map<T>& service,
                                         std::size_t nMoments = 3) {
    const BmapMap1Blocks<T> b = solver_mam_bmap_map_1_blocks(D, service);
    bmap_detail::require_double_arith<T>("solver_mam_bmap_map_1");
    if (nMoments < 1) throw InputError("solver_mam_bmap_map_1: nMoments must be positive");

    std::vector<Matrix<double>> Ablk, Bblk;
    Ablk.push_back(bmap_detail::as_double(b.A0));
    Ablk.push_back(bmap_detail::as_double(b.A1));
    for (std::size_t k = 0; k < b.Aup.size(); ++k) Ablk.push_back(bmap_detail::as_double(b.Aup[k]));
    Bblk.push_back(bmap_detail::as_double(b.B0));
    for (std::size_t k = 0; k < b.Bup.size(); ++k) Bblk.push_back(bmap_detail::as_double(b.Bup[k]));

    const Matrix<double> A = bmap_detail::hstack(Ablk);
    const Matrix<double> B = bmap_detail::hstack(Bblk);

    const Matrix<double> G = smc::mg1_g_etaqa(A);
    const std::vector<double> pi = smc::mg1_pi_etaqa(B, A, G);

    BmapQueueResult<T> out;
    out.qlenMoments.reserve(nMoments);
    for (std::size_t n = 1; n <= nMoments; ++n)
        out.qlenMoments.push_back(
            num_traits<T>::from_double(smc::mg1_qlen_etaqa(B, A, pi, n)));

    out.QN = out.qlenMoments[0];
    out.UN = b.rho;
    out.TN = b.lambda;
    out.RN = T(out.QN / out.TN);
    out.piAgg = Matrix<T>(1, pi.size(), num_traits<T>::from_int(0));
    for (std::size_t j = 0; j < pi.size(); ++j) out.piAgg(0, j) = num_traits<T>::from_double(pi[j]);
    out.fund = Matrix<T>(G.rows(), G.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < G.rows(); ++i)
        for (std::size_t j = 0; j < G.cols(); ++j)
            out.fund(i, j) = num_traits<T>::from_double(G(i, j));
    return out;
}

/**
 * Port of `solver_mam_map_bmap_1.m`, mean measures included.
 *
 * The GI/M/1-type chain is stacked as `A = [A0; A1; ...; A_{K+1}]` and
 * `B = [B1; B2; ...; B_{K+1}; 0]`, which is the reference's own layout down to
 * the trailing zero block its `zeros(m*(K+1)+m, m)` allocation leaves unfilled.
 * R comes from `GIM1_R_ETAQA`, the aggregates and the mean queue length from
 * `GIM1_pi_ETAQA` and `GIM1_qlen_ETAQA`, both with A0 as the Boundary block.
 * The reference asks for the FIRST moment only on this side, and so does this.
 */
template <class T>
BmapQueueResult<T> solver_mam_map_bmap_1(const Map<T>& arrival, const std::vector<Matrix<T>>& D) {
    const MapBmap1Blocks<T> b = solver_mam_map_bmap_1_blocks(arrival, D);
    bmap_detail::require_double_arith<T>("solver_mam_map_bmap_1");

    std::vector<Matrix<double>> Ablk, Bblk;
    Ablk.push_back(bmap_detail::as_double(b.A0));
    Ablk.push_back(bmap_detail::as_double(b.A1));
    for (std::size_t k = 0; k < b.Adown.size(); ++k)
        Ablk.push_back(bmap_detail::as_double(b.Adown[k]));
    Bblk.push_back(bmap_detail::as_double(b.B1));
    for (std::size_t k = 0; k < b.Bto0.size(); ++k)
        Bblk.push_back(bmap_detail::as_double(b.Bto0[k]));
    Bblk.push_back(Matrix<double>(b.m, b.m, 0.0));  // the reference's unfilled tail block

    const Matrix<double> A = bmap_detail::vstack(Ablk);
    const Matrix<double> B = bmap_detail::vstack(Bblk);
    const Matrix<double> A0 = bmap_detail::as_double(b.A0);

    const Matrix<double> R = smc::gim1_r_etaqa(A);
    const std::vector<double> pi = smc::gim1_pi_etaqa(B, A, R, A0);
    const double QN = smc::gim1_qlen_etaqa(B, A, R, pi, 1, A0);

    BmapQueueResult<T> out;
    out.qlenMoments.push_back(num_traits<T>::from_double(QN));
    out.QN = out.qlenMoments[0];
    out.UN = b.rho;
    out.TN = b.lambda;
    out.RN = T(out.QN / out.TN);
    out.piAgg = Matrix<T>(1, pi.size(), num_traits<T>::from_int(0));
    for (std::size_t j = 0; j < pi.size(); ++j) out.piAgg(0, j) = num_traits<T>::from_double(pi[j]);
    out.fund = Matrix<T>(R.rows(), R.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < R.rows(); ++i)
        for (std::size_t j = 0; j < R.cols(); ++j)
            out.fund(i, j) = num_traits<T>::from_double(R(i, j));
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_BMAP_H
