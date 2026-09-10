/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_FJ_CODES_H
#define LINE_API_FJ_FJ_CODES_H

/**
 * FJ_codes, the fork-join response-time-tail approximation of Z. Qiu, J. F.
 * Perez and P. Harrison, "Beyond the Mean in Fork-Join Queues: Efficient
 * Approximation for Response-Time Tails" (IFIP Performance 2015). Third-party,
 * BSD-3-Clause, Copyright 2015 Imperial College London; see
 * THIRD-PARTY-NOTICES.md and `python/line_solver/lib/thirdparty/fj/LICENSE.txt`.
 *
 * Port of the solve layer of `matlab/lib/thirdparty/FJ_codes`: `computeT.m`,
 * `computeT_NARE.m`, `computePi.m`, `returnWait.m`, `returnPer.m`,
 * `returnRT1.m`, `returnRT2.m` and `mainFJ.m`. The state-space construction is
 * `fj_codes_matrices.h`.
 *
 * THE METHOD IN ONE PARAGRAPH. The response-time percentiles of the ONE-node
 * queue are exact (it is a MAP/PH/1 queue and its sojourn time is a phase-type
 * law). The TWO-node fork-join queue is solved by the approximation of Section
 * 4: the queue-length DIFFERENCE between the two branches is truncated at C, so
 * the all-busy period becomes a finite-phase Markov-modulated fluid whose
 * generator T solves a Riccati equation, and the waiting and response times
 * come out as phase-type laws over that phase space. A K-node queue is then
 * INTERPOLATED (Section 6) as RT_1 + (RT_2 - RT_1) * log(K) / log(2), one
 * percentile at a time. Nothing about K enters the matrices: K appears only in
 * that last line.
 *
 * WHAT `returnRT1` USES HERE. The reference calls `Q_CT_MAP_MAP_1` of the QMAM
 * toolbox for the one-node sojourn time. This port calls
 * `mmapph1fcfs_stdistr_ph`, which returns the same object -- the sojourn time
 * of the MAP/PH/1 queue as a phase-type pair (alpha, A) -- through BUTools'
 * MMAPPH1FCFS, the route the Java port also takes. Both are exact for this
 * queue, so this is a change of implementation and not of method. The Java port
 * additionally CATCHES a failure of that solve and substitutes the raw SERVICE
 * time PH; that is not reproduced, because a service time reported as a
 * response time is a wrong number rather than a degraded one.
 *
 * ARITHMETIC. `double` only. `computeT_NARE` needs an ORDERED real Schur
 * factorization (the stable invariant subspace of a 2m x 2m Hamiltonian-like
 * pencil), and both Sylvester equations in `computePi` are of order (C + 1) m^2
 * ma, which is in the hundreds to low thousands at the default C = 100 -- far
 * past what the field-generic Kronecker Sylvester solver can carry. Both are
 * LAPACK, exactly as `util/eig.h` is. Callers at another arithmetic must refuse
 * by name; `solver_mam_fj.h` does.
 *
 * REFERENCE DEFECTS, reproduced unless stated:
 *
 *  1. `mainFJ` LOOPS OVER `Cs` AND KEEPS ONLY THE LAST. `percentileRT_1` and
 *     `percentileRT_2` are overwritten each iteration and only the final C
 *     survives, so a vector of accuracies costs the full solve per entry and
 *     answers for one of them. Reproduced; LINE passes a scalar.
 *  2. `returnRT1` IS INSIDE THAT LOOP although it does not depend on C.
 *     Reproduced, and it is why a two-entry `Cs` doubles the MAP/PH/1 solve.
 *  3. `returnRT2` COMPUTES `percentileWait` AND DISCARDS IT. The waiting-time
 *     percentiles are a complete result of the same phase-type pair the
 *     response time is then built from. NOT reproduced: computing an unused
 *     percentile inversion is pure cost, and the inversion is the expensive
 *     step. `fj_return_rt2` returns the waiting-time PH pair instead, so a
 *     caller that wants those percentiles can invert it without a second solve.
 *  4. `returnPer` SCANS THE INVERSE CDF ON A FIXED 0.001 GRID downwards from
 *     3 * mean, which quantizes every percentile to a millisecond of model time
 *     regardless of the time scale of the model. Reproduced: the grid is the
 *     method's resolution and changing it would move every reported number.
 *  5. `returnPer` CAN FALL OFF ITS OWN SCAN. If the scan reaches t = 0 without
 *     the CDF dropping below the target, `temp_percentileRT` keeps the value
 *     from the PREVIOUS percentile, or is undefined on the first. NOT
 *     reproduced: this port reports it by name.
 *  6. `computeT`'s Sylvester iteration HAS NO ITERATION CAP. NOT reproduced:
 *     the loop is capped and a failure to converge is reported by name rather
 *     than hanging.
 *  7. THE RESIDUAL NORMS ARE PRINTED, from `computeT` and `computeT_NARE`, on
 *     every call. NOT reproduced: they are returned in `FjCodesT::residual`
 *     instead, because a solver that writes to stdout corrupts the CLI's JSON.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/fj/fj_codes_matrices.h"
#include "line/api/fj/fj_dist2fj.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmapph1fcfs.h"
#include "line/util/eig.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lstsq.h"
#include "line/util/matrix.h"
#include "line/util/sylvester.h"

namespace line {
namespace fj {

/** Which route `computeT.m` takes to the T matrix. */
enum class FjTMode { Nare, Sylvester };

/** Parse the reference's `T_Mode` string, whose default is 'NARE'. */
inline FjTMode fj_parse_tmode(const std::string& s) {
    // The reference tests `strfind(T_Mode, 'Sylvest') > 0` and defaults to NARE
    // for everything else, including an empty string and a misspelling.
    if (s.find("Sylvest") != std::string::npos) return FjTMode::Sylvester;
    return FjTMode::Nare;
}

/** What `computeT.m` returns. */
struct FjCodesT {
    Matrix<double> T;                ///< the all-busy generator, (newdim * ma) square
    Matrix<double> S;                ///< build_SA's S, newdim square
    Matrix<double> A_jump;           ///< build_SA's A_jump, newdim square
    Matrix<double> S_Arr;            ///< kron(S, I_ma)
    std::vector<double> sum_Ajump;   ///< row sums of kron(A_jump, I_ma)
    double residual = 0.0;           ///< inf-norm the reference prints
    std::size_t iterations = 0;      ///< Sylvester mode only
};

/** What `computePi.m` returns. */
struct FjCodesPi {
    std::vector<double> pi0;  ///< unnormalized, length newdim * ma
    double En1 = 0.0;         ///< mean number of arrivals in a not-all-busy period
};

/** What `returnWait.m` returns: the waiting time as a phase-type law. */
struct FjCodesWait {
    std::vector<double> wait_alpha;  ///< defective, mass prob_wait
    Matrix<double> wait_Smat;
    double prob_wait = 0.0;
    std::vector<double> alfa;  ///< -pi0 T^-1, the all-busy occupancy
};

/** One line of `mainFJ`'s output cell: the percentiles of a K-node queue. */
struct FjCodesPercentiles {
    std::size_t K = 0;
    std::vector<double> percentiles;  ///< in PERCENT, as the reference stores them
    std::vector<double> RTp;
};

namespace fjdetail {

/** MATLAB `A / B` for a row vector: x with x B = a. */
inline std::vector<double> rdivide_row(const std::vector<double>& a, const Matrix<double>& B) {
    const Matrix<double> iB = inverse(B);
    return vecmul(a, iB);
}

/** Elementwise infinity norm of A - B. */
inline double max_abs_diff(const Matrix<double>& A, const Matrix<double>& B) {
    double d = 0.0;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) {
            const double x = std::fabs(A(i, j) - B(i, j));
            if (x > d) d = x;
        }
    return d;
}

/** The inf-norm (maximum absolute row sum) the reference reports. */
inline double inf_norm(const Matrix<double>& A) {
    double best = 0.0;
    for (std::size_t i = 0; i < A.rows(); ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < A.cols(); ++j) s += std::fabs(A(i, j));
        if (s > best) best = s;
    }
    return best;
}

inline Matrix<double> madd(const Matrix<double>& A, const Matrix<double>& B) {
    Matrix<double> C(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) = A(i, j) + B(i, j);
    return C;
}

inline Matrix<double> msub(const Matrix<double>& A, const Matrix<double>& B) {
    Matrix<double> C(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) = A(i, j) - B(i, j);
    return C;
}

/** A sub-block of A, rows [r0, r1) and columns [c0, c1). */
inline Matrix<double> block(const Matrix<double>& A, std::size_t r0, std::size_t r1,
                            std::size_t c0, std::size_t c1) {
    Matrix<double> B(r1 - r0, c1 - c0, 0.0);
    for (std::size_t i = r0; i < r1; ++i)
        for (std::size_t j = c0; j < c1; ++j) B(i - r0, j - c0) = A(i, j);
    return B;
}

}  // namespace fjdetail

/**
 * Port of `computeT_NARE.m`: the T matrix as the stable invariant subspace of
 *
 *   H = [ I (x) D0        I (x) D1 ]
 *       [ -A_jump (x) I   -S       ]
 *
 * The m eigenvalues of SMALLEST real part are ordered to the front of the real
 * Schur form, X is read off the resulting basis as Q1(m+1:2m, 1:m) /
 * Q1(1:m, 1:m), and T = S + X (I (x) D1).
 *
 * A complex-conjugate pair astride the m/2m boundary would ask for an invariant
 * subspace that does not exist over the reals; LAPACK's dtrexc refuses to split
 * such a block, and the split is detected here first so that it is named as the
 * spectrum condition it is rather than as a reordering failure.
 */
inline Matrix<double> fj_compute_t_nare(const Matrix<double>& D0, const Matrix<double>& D1,
                                        const Matrix<double>& S, const Matrix<double>& A_jump,
                                        double* residual) {
    const std::size_t m = S.rows();
    const std::size_t ma = D0.rows();
    if (ma == 0 || m % ma != 0)
        throw InputError("fj_compute_t_nare: the phase space is not a multiple of the arrival "
                         "order");
    const std::size_t ms = m / ma;
    const Matrix<double> Ims = eye<double>(ms);
    const Matrix<double> Ima = eye<double>(ma);
    const Matrix<double> ImsD0 = mam::kron(Ims, D0);
    const Matrix<double> ImsD1 = mam::kron(Ims, D1);
    const Matrix<double> AjIma = mam::kron(A_jump, Ima);

    Matrix<double> H(2 * m, 2 * m, 0.0);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) {
            H(i, j) = ImsD0(i, j);
            H(i, m + j) = ImsD1(i, j);
            H(m + i, j) = -AjIma(i, j);
            H(m + i, m + j) = -S(i, j);
        }

    const RealSchur sc = schur_decomposition(H);
    // The real part of each diagonal entry: a 2 x 2 block contributes its own
    // trace/2 to both of its entries, which is what ordeig reports for the pair.
    const std::size_t n = 2 * m;
    std::vector<double> re(n, 0.0);
    std::vector<int> blk(n, 1);  // 1 for a 1 x 1 block, 2 for the first row of a 2 x 2
    for (std::size_t i = 0; i < n;) {
        const bool pair = (i + 1 < n) && sc.T(i + 1, i) != 0.0;
        if (pair) {
            const double r = 0.5 * (sc.T(i, i) + sc.T(i + 1, i + 1));
            re[i] = r;
            re[i + 1] = r;
            blk[i] = 2;
            blk[i + 1] = 0;
            i += 2;
        } else {
            re[i] = sc.T(i, i);
            blk[i] = 1;
            i += 1;
        }
    }
    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(),
                     [&re](std::size_t a, std::size_t b) { return re[a] < re[b]; });
    std::vector<double> key(n, 0.0);
    for (std::size_t i = 0; i < m; ++i) key[order[i]] = 1.0;
    for (std::size_t i = 0; i + 1 < n; ++i)
        if (blk[i] == 2 && key[i] != key[i + 1])
            throw NumericError(
                "fj_compute_t_nare: the m eigenvalues of smallest real part split a complex "
                "conjugate pair, so the stable invariant subspace of the Riccati pencil is not "
                "real and the T matrix of this model is not defined");

    const RealSchur ord = schur_reorder(sc, key);
    const Matrix<double> Q11 = fjdetail::block(ord.Z, 0, m, 0, m);
    const Matrix<double> Q21 = fjdetail::block(ord.Z, m, 2 * m, 0, m);
    const Matrix<double> X = matmul(Q21, inverse(Q11));
    const Matrix<double> T = fjdetail::madd(S, matmul(X, ImsD1));
    if (residual != nullptr)
        *residual = fjdetail::inf_norm(
            fjdetail::madd(fjdetail::madd(matmul(T, X), matmul(X, ImsD0)), AjIma));
    return T;
}

/**
 * Port of `computeT.m`.
 *
 * The Sylvester route is the paper's Section 5.1 and the NARE route its Section
 * 5.2; they solve the same equation and the reference defaults to NARE. Both
 * are offered because `options.config.fj_tmode` selects between them, and
 * because they fail on different models: the iteration converges linearly and
 * can stall, while the Schur route is direct but needs the stable subspace to
 * be real.
 */
inline FjCodesT fj_compute_t(const FjDist<double>& arrival, const FjDist<double>& service,
                             const FjCodesServiceH& h, std::size_t C, FjTMode mode) {
    const FjCodesSA sa = fj_build_sa(service, h, C);
    const std::size_t d0 = arrival.lambda0.rows();
    const Matrix<double> Id0 = eye<double>(d0);

    FjCodesT out;
    out.S = sa.S;
    out.A_jump = sa.A_jump;
    out.S_Arr = mam::kron(sa.S, Id0);
    const Matrix<double> A_jump_Arr = mam::kron(sa.A_jump, Id0);

    if (mode == FjTMode::Sylvester) {
        const std::size_t ms = sa.S.rows();
        const std::size_t m = ms * d0;
        const Matrix<double> ID0 = mam::kron(eye<double>(ms), arrival.lambda0);
        const Matrix<double> DS =
            matmul(mam::kron(eye<double>(ms), arrival.lambda1), A_jump_Arr);
        const Matrix<double> Im = eye<double>(m);
        Matrix<double> Tnew = out.S_Arr;
        Matrix<double> Told(m, m, 0.0);
        Matrix<double> L;
        // The reference iterates without a cap; 500 is far past the linear
        // convergence of every model the gate admits, and a stall is reported.
        const std::size_t max_iter = 500;
        std::size_t it = 0;
        while (fjdetail::max_abs_diff(Told, Tnew) > 1e-10) {
            if (++it > max_iter)
                throw NumericError(
                    "fj_compute_t: the Sylvester iteration of computeT.m did not reach 1e-10 in " +
                    std::to_string(max_iter) +
                    " steps; solve this model with the NARE route (config.fj_tmode = 'NARE')");
            Told = Tnew;
            L = lyap_schur(Tnew, ID0, Im);
            Tnew = fjdetail::madd(out.S_Arr, matmul(L, DS));
        }
        out.T = Tnew;
        out.iterations = it;
        if (it > 0)
            out.residual = fjdetail::inf_norm(
                fjdetail::madd(fjdetail::madd(matmul(out.T, L), matmul(L, ID0)), Im));
    } else {
        out.T = fj_compute_t_nare(arrival.lambda0, arrival.lambda1, out.S_Arr, sa.A_jump,
                                  &out.residual);
    }

    out.sum_Ajump.assign(A_jump_Arr.rows(), 0.0);
    for (std::size_t i = 0; i < A_jump_Arr.rows(); ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < A_jump_Arr.cols(); ++j) s += A_jump_Arr(i, j);
        out.sum_Ajump[i] = s;
    }
    return out;
}

/**
 * The boundary solve both branches of `computePi.m` end with:
 *
 *   pi0 [ pi0mat - I , T^-1 e ] = [ 0 ... 0 , -1 ]
 *
 * an OVERDETERMINED system by one column, which MATLAB's mrdivide answers in
 * the least-squares sense. It is consistent -- the extra column is the
 * normalization that fixes the scale of the null vector -- so least squares
 * returns the exact solution and not an approximation.
 */
inline std::vector<double> fj_boundary_solve(const Matrix<double>& pi0mat,
                                             const Matrix<double>& T) {
    const std::size_t nd = pi0mat.rows();
    if (pi0mat.cols() != nd || T.rows() != nd || T.cols() != nd)
        throw InputError("fj_boundary_solve: the boundary blocks are not conformable");
    const Matrix<double> iT = inverse(T);
    Matrix<double> M(nd, nd + 1, 0.0);
    for (std::size_t i = 0; i < nd; ++i) {
        for (std::size_t j = 0; j < nd; ++j) M(i, j) = pi0mat(i, j) - (i == j ? 1.0 : 0.0);
        double s = 0.0;
        for (std::size_t j = 0; j < nd; ++j) s += iT(i, j);
        M(i, nd) = s;
    }
    std::vector<double> b(nd + 1, 0.0);
    b[nd] = -1.0;
    const LstsqResult<double> r = lstsq(M.transpose(), b);
    return r.x;
}

/**
 * Port of `computePi.m`: the all-busy boundary vector and E[n1].
 *
 * The two branches are not two implementations of one formula. For EXPONENTIAL
 * service the not-all-busy space has exactly as many phases as the all-busy one
 * (m = 1 makes both (C + 1)-dimensional), so the return to the all-busy period
 * is a square map and `Q0` can be inverted directly. For phase-type service the
 * two spaces differ and the reference goes through `constructSRK`, whose Ke and
 * Kc project between them. Running the second branch on exponential service
 * would give the same answer; running the first on anything else is a shape
 * error, which is why the reference switches on `SerChoice`.
 */
inline FjCodesPi fj_compute_pi(const Matrix<double>& T, const FjDist<double>& arrival,
                               const FjDist<double>& service, const FjCodesServiceH& h,
                               std::size_t C, const Matrix<double>& S,
                               const Matrix<double>& A_jump) {
    const std::size_t da = arrival.lambda0.rows();
    FjCodesPi out;

    if (service.choice == 1) {
        const std::size_t ms = S.cols();
        const Matrix<double> S_notallbusy = fj_construct_not_all_busy(C, service, h);
        const Matrix<double> Q0 = mam::krons(S_notallbusy, arrival.lambda0);
        if (Q0.rows() != ms * da)
            throw NumericError(
                "fj_compute_pi: the not-all-busy space and the all-busy space have different "
                "dimensions, which the exponential branch of computePi.m assumes they do not");
        const Matrix<double> Igral =
            sylvester_schur(T, mam::kron(eye<double>(ms), arrival.lambda0), eye<double>(da * ms));
        const Matrix<double> pi0mat =
            matmul(matmul(matmul(Igral, mam::kron(A_jump, eye<double>(da))), inverse(Q0)),
                   mam::kron(eye<double>(ms), arrival.lambda1));
        out.pi0 = fj_boundary_solve(pi0mat, T);
        double sp = 0.0;
        for (double v : out.pi0) sp += v;
        const std::vector<double> pm = vecmul(out.pi0, pi0mat);
        double spm = 0.0;
        for (double v : pm) spm += v;
        out.En1 = spm / sp;
        return out;
    }

    const FjCodesSRK srk = fj_construct_srk(C, service, h, S);
    const std::size_t dtmat = T.rows() / arrival.lambda0.cols();
    const std::size_t dsexp = srk.Se.rows();
    const std::size_t dnb = dsexp - dtmat;

    const Matrix<double> Sedash = fjdetail::block(srk.Se, dtmat, dsexp, dtmat, dsexp);
    const Matrix<double> Rbusy = fjdetail::block(srk.R0, dtmat, dsexp, 0, dsexp);
    Matrix<double> Iidle_small(dsexp, dnb, 0.0);
    for (std::size_t i = 0; i < dnb; ++i) Iidle_small(dtmat + i, i) = 1.0;
    const Matrix<double> Ida = eye<double>(da);
    const Matrix<double> Iidle = mam::kron(Iidle_small, Ida);

    const Matrix<double> Qidle = mam::krons(Sedash, arrival.lambda0);
    const Matrix<double> Qbusy = mam::kron(Rbusy, arrival.lambda1);
    const Matrix<double> Bmap = [&] {
        const Matrix<double> M = matmul(matmul(Iidle, inverse(Qidle)), Qbusy);
        Matrix<double> N(M.rows(), M.cols());
        for (std::size_t i = 0; i < M.rows(); ++i)
            for (std::size_t j = 0; j < M.cols(); ++j) N(i, j) = -M(i, j);
        return N;
    }();
    const Matrix<double> Kemap = mam::kron(srk.Ke, Ida);
    const Matrix<double> Kcmap = mam::kron(srk.Kc, Ida);

    const Matrix<double> BB = mam::kron(eye<double>(srk.Sestar.cols()), arrival.lambda0);
    const Matrix<double> Igral =
        lyap_schur(T, BB, matmul(Kemap, mam::kron(srk.Sestar, Ida)));
    const Matrix<double> pi0mat = matmul(matmul(Igral, Bmap), Kcmap);
    out.pi0 = fj_boundary_solve(pi0mat, T);
    double sp = 0.0;
    for (double v : out.pi0) sp += v;

    // -pi0 Igral Iidle Qidle^-1 (I (x) D1), left to right as MATLAB reads it.
    std::vector<double> row = vecmul(out.pi0, Igral);
    for (double& v : row) v = -v;
    row = vecmul(row, Iidle);
    row = fjdetail::rdivide_row(row, Qidle);
    row = vecmul(row, mam::kron(eye<double>(Sedash.cols()), arrival.lambda1));
    double sr = 0.0;
    for (double v : row) sr += v;
    out.En1 = sr / sp;
    return out;
}

/** Port of `returnWait.m`: the stationary waiting time as a phase-type law. */
inline FjCodesWait fj_return_wait(double En1, const std::vector<double>& pi0,
                                  const Matrix<double>& T, const std::vector<double>& phi,
                                  const std::vector<double>& sum_Ajump) {
    const std::size_t ds = phi.size();
    FjCodesWait out;
    out.alfa = fjdetail::rdivide_row(pi0, T);
    for (double& v : out.alfa) v = -v;

    double ap = 0.0;
    for (std::size_t i = 0; i < ds; ++i) ap += out.alfa[i] * phi[i];
    std::vector<double> rhos(ds, 0.0);
    for (std::size_t i = 0; i < ds; ++i) rhos[i] = phi[i] * out.alfa[i] / ap;

    double sp = 0.0;
    for (double v : pi0) sp += v;
    double asum = 0.0;
    for (std::size_t i = 0; i < ds; ++i) asum += out.alfa[i] * sum_Ajump[i];
    const double En0 = asum / sp;

    out.prob_wait = (En0 - 1.0) / (En0 - 1.0 + En1);
    out.wait_alpha.assign(ds, 0.0);
    for (std::size_t i = 0; i < ds; ++i) out.wait_alpha[i] = out.prob_wait * rhos[i];

    out.wait_Smat = Matrix<double>(ds, ds, 0.0);
    for (std::size_t i = 0; i < ds; ++i)
        for (std::size_t j = 0; j < ds; ++j)
            out.wait_Smat(i, j) = out.alfa[j] * T(j, i) / out.alfa[i];
    return out;
}

/**
 * Port of `returnPer.m`: the percentiles of a (possibly defective) phase-type
 * law, by uniformization.
 *
 * The law is uniformized at c = max(-diag(A)) into P = A / c + I, the
 * absorption probability by time t is the Poisson mixture sum_k p_k(ct) a_k
 * with a_k = alpha P^k e, and the series is truncated where its partial sums
 * reach the total mass alpha (I - P)^-1 e. The percentile itself is then found
 * by scanning t downwards on a 0.001 grid from three times the mean, extending
 * the bracket by half a mean whenever the target is not yet reached. A target
 * below the defect 1 - sum(alpha) is answered with zero, which is the
 * reference's convention for a percentile the law never attains.
 */
inline std::vector<double> fj_return_per(const std::vector<double>& vec, const Matrix<double>& A,
                                         const std::vector<double>& pers) {
    using std::exp;
    const std::size_t m = A.cols();
    if (vec.size() != m || A.rows() != m)
        throw InputError("fj_return_per: the phase-type pair is not conformable");

    const std::vector<double> negvA = [&] {
        std::vector<double> v = vecmul(vec, inverse(A));
        for (double& x : v) x = -x;
        return v;
    }();
    double meanRT = 0.0;
    for (double v : negvA) meanRT += v;
    if (!(meanRT > 0.0))
        throw NumericError("fj_return_per: the phase-type law has a non-positive mean, so its "
                           "percentiles are not defined");

    double c = 0.0;
    for (std::size_t i = 0; i < m; ++i) c = std::max(c, -A(i, i));
    if (!(c > 0.0))
        throw NumericError("fj_return_per: the phase-type generator has no negative diagonal");
    Matrix<double> P(m, m, 0.0);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) P(i, j) = A(i, j) / c + (i == j ? 1.0 : 0.0);

    Matrix<double> ImP(m, m, 0.0);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) ImP(i, j) = (i == j ? 1.0 : 0.0) - P(i, j);
    const std::vector<double> vImP = vecmul(vec, inverse(ImP));
    double M = 0.0;
    for (double v : vImP) M += v;

    double a0 = 0.0;
    for (double v : vec) a0 += v;
    double sum_a = a0;
    std::vector<double> ak;
    std::vector<double> vP(m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < m; ++j) s += P(i, j);
        vP[i] = s;
    }
    while (std::fabs(sum_a - M) >= 1e-10) {
        double t = 0.0;
        for (std::size_t i = 0; i < m; ++i) t += vec[i] * vP[i];
        ak.push_back(t);
        sum_a += t;
        vP = mulvec(P, vP);
        if (ak.size() > 2000000)
            throw NumericError("fj_return_per: the uniformized Poisson series did not reach the "
                               "total absorption mass, so the percentile scan cannot terminate");
    }
    const std::size_t K1 = ak.size();

    // The absorption CDF at t, as the reference evaluates it.
    const auto cdf_at = [&](double t) {
        double pM = exp(-c * t);
        double F = pM * a0;
        for (std::size_t k = 1; k <= K1; ++k) {
            pM = c * t * pM / static_cast<double>(k);
            F += pM * ak[k - 1];
        }
        return 1.0 - F;
    };

    std::vector<double> out(pers.size(), 0.0);
    for (std::size_t p = 0; p < pers.size(); ++p) {
        if (pers[p] < 1.0 - a0) continue;  // the law never attains this percentile
        double MaxTime = 3.0 * meanRT;
        for (;;) {
            if (cdf_at(MaxTime) < pers[p]) {
                MaxTime += 0.5 * meanRT;
                continue;
            }
            bool found = false;
            // MATLAB's colon operator indexes the grid, it does not accumulate
            // a subtraction, so t is formed as MaxTime - k * 0.001.
            const std::size_t steps = static_cast<std::size_t>(MaxTime / 0.001);
            for (std::size_t k = 0; k <= steps; ++k) {
                const double t = MaxTime - static_cast<double>(k) * 0.001;
                if (cdf_at(t) < pers[p]) {
                    out[p] = t + 0.001;
                    found = true;
                    break;
                }
            }
            if (!found)
                throw NumericError(
                    "fj_return_per: the downward scan reached t = 0 without the response-time CDF "
                    "falling below " + std::to_string(pers[p]) +
                    ", so this percentile has no bracket; the reference silently reports the "
                    "previous percentile here");
            break;
        }
    }
    return out;
}

/**
 * Port of `returnRT1.m`: the response-time percentiles of the ONE-node queue,
 * which are exact.
 */
inline std::vector<double> fj_return_rt1(const FjDist<double>& arrival,
                                         const FjDist<double>& service,
                                         const std::vector<double>& pers) {
    mam::Mmap<double> mm;
    mm.D0 = arrival.lambda0;
    mm.D1 = arrival.lambda1;
    mm.Dc.push_back(arrival.lambda1);

    std::vector<mam::PhService<double>> svc(1);
    svc[0].sigma = service.tau_st;
    svc[0].S = service.ST;

    const std::vector<mam::StDistrPh<double>> ph = mam::mmapph1fcfs_stdistr_ph(mm, svc);
    if (ph.empty())
        throw NumericError("fj_return_rt1: MMAPPH1FCFS returned no sojourn-time law");
    return fj_return_per(ph[0].alpha, ph[0].A, pers);
}

/** What `returnRT2.m` produces, plus the waiting-time law it discards. */
struct FjCodesRT2 {
    std::vector<double> RTp;             ///< response-time percentiles
    std::vector<double> wait_alpha;      ///< the waiting-time law, defective
    Matrix<double> wait_Smat;
    double prob_wait = 0.0;
    double residual = 0.0;               ///< the T-matrix residual of computeT
};

/**
 * Port of `returnRT2.m`: the response-time percentiles of the TWO-node
 * fork-join queue, under the Section 4 approximation.
 *
 * The response time is assembled as one phase-type law over three blocks: the
 * service process of a job that arrives in a not-all-busy period, the TIME-
 * REVERSED service process of a job that arrives in an all-busy period, and the
 * waiting time. The reversal is what lets the job's own service be appended to
 * the waiting time it accrued: `stat_service_phase` is the stationary phase
 * occupancy, `tr_ST` is the reversed generator, and `tildeP` is the coupling
 * that hands the reversed process over to the waiting-time block. Phases of
 * zero stationary occupancy are dropped -- they are unreachable and the
 * reversal divides by their occupancy.
 */
inline FjCodesRT2 fj_return_rt2(const FjDist<double>& arrival, const FjDist<double>& service,
                                const std::vector<double>& pers, std::size_t C, FjTMode mode) {
    const FjCodesServiceH h = fj_build_service_h(service);
    const FjCodesT ct = fj_compute_t(arrival, service, h, C, mode);
    const std::size_t mWait = ct.A_jump.rows();
    const std::size_t n = ct.T.rows();

    std::vector<double> phi(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < n; ++j) s += ct.T(i, j) - ct.S_Arr(i, j);
        phi[i] = s;
    }
    const FjCodesPi pi = fj_compute_pi(ct.T, arrival, service, h, C, ct.S, ct.A_jump);
    const FjCodesWait w = fj_return_wait(pi.En1, pi.pi0, ct.T, phi, ct.sum_Ajump);

    const FjCodesGenService gs = fj_generate_service(service, h, C, ct.S);
    const Matrix<double> ST = mam::kron(gs.T, arrival.Ia);

    std::vector<double> pi0n = pi.pi0;
    double sp = 0.0;
    for (double v : pi0n) sp += v;
    for (double& v : pi0n) v /= sp;

    const std::size_t dim = arrival.ma * gs.newdim;
    const std::size_t dim_notbusy = arrival.ma * gs.dim_notbusy;
    const std::size_t dim_service = dim + dim_notbusy;
    const std::size_t Tr = ST.rows();
    const std::size_t Sc = w.wait_Smat.cols();
    if (Tr != dim_service)
        throw NumericError("fj_return_rt2: the tagged-job service space and the phase blocks "
                           "disagree in size");

    std::vector<double> notbusy_start(dim_service, 0.0);
    for (std::size_t i = 0; i < dim; ++i)
        notbusy_start[i] = (1.0 - w.prob_wait) * pi0n[i];

    // TS = T - S_Arr, the rates at which a new job enters service, then row
    // normalized into the jump kernel of the all-busy phase.
    Matrix<double> TS = fjdetail::msub(ct.T, ct.S_Arr);
    std::vector<double> busy_start(dim_service, 0.0);
    {
        const std::vector<double> aTS = vecmul(w.alfa, TS);
        if (aTS.size() != dim)
            throw NumericError("fj_return_rt2: the all-busy phase space and the tagged job's "
                               "all-busy service block have different sizes");
        double s = 0.0;
        for (double v : aTS) s += v;
        for (std::size_t i = 0; i < dim; ++i) busy_start[i] = w.prob_wait * aTS[i] / s;
    }
    for (std::size_t i = 0; i < TS.rows(); ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < TS.cols(); ++j) s += TS(i, j);
        for (std::size_t j = 0; j < TS.cols(); ++j) TS(i, j) /= s;
    }

    const std::vector<double> stat = [&] {
        Matrix<double> negST(Tr, Tr);
        for (std::size_t i = 0; i < Tr; ++i)
            for (std::size_t j = 0; j < Tr; ++j) negST(i, j) = -ST(i, j);
        return fjdetail::rdivide_row(busy_start, negST);
    }();
    std::vector<bool> nz(Tr, false);
    for (std::size_t i = 0; i < Tr; ++i) nz[i] = stat[i] > 0.0;

    std::vector<double> tr_start(Tr, 0.0);
    for (std::size_t j = 0; j < Tr; ++j) {
        double s = 0.0;
        for (std::size_t k = 0; k < Tr; ++k) s += ST(j, k);
        tr_start[j] = -s * stat[j];
    }

    Matrix<double> tr_ST(Tr, Tr, 0.0);
    for (std::size_t i = 0; i < Tr; ++i) {
        if (!nz[i]) continue;
        for (std::size_t j = 0; j < Tr; ++j)
            if (nz[j]) tr_ST(i, j) = ST(j, i) * stat[j] / stat[i];
    }
    std::vector<double> tr_exit(Tr, 0.0);
    for (std::size_t i = 0; i < Tr; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < Tr; ++j) s += tr_ST(i, j);
        tr_exit[i] = -s;
    }

    // TS2 = TS' diag(alfa), row normalized with the reference's 1e-11 floor.
    Matrix<double> TS2(Sc, Sc, 0.0);
    for (std::size_t i = 0; i < Sc; ++i)
        for (std::size_t j = 0; j < Sc; ++j) TS2(i, j) = TS(j, i) * w.alfa[j];
    for (std::size_t i = 0; i < Sc; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < Sc; ++j) s += TS2(i, j);
        if (std::fabs(s) < 10e-12) s = 1.0;
        for (std::size_t j = 0; j < Sc; ++j) TS2(i, j) /= s;
    }

    Matrix<double> tildeP(Tr, Sc, 0.0);
    for (std::size_t i = 0; i < Sc; ++i)
        for (std::size_t j = 0; j < Sc; ++j) tildeP(i, j) = tr_exit[i] * TS2(i, j);

    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < Tr; ++i)
        if (nz[i]) keep.push_back(i);
    const std::size_t m_tr = keep.size();

    const std::size_t total = Tr + m_tr + Sc;
    std::vector<double> gamma(total, 0.0);
    for (std::size_t i = 0; i < Tr; ++i) gamma[i] = notbusy_start[i];
    for (std::size_t i = 0; i < m_tr; ++i) gamma[Tr + i] = tr_start[keep[i]];

    Matrix<double> Cres(total, total, 0.0);
    for (std::size_t i = 0; i < Tr; ++i)
        for (std::size_t j = 0; j < Tr; ++j) Cres(i, j) = ST(i, j);
    for (std::size_t i = 0; i < m_tr; ++i) {
        for (std::size_t j = 0; j < m_tr; ++j) Cres(Tr + i, Tr + j) = tr_ST(keep[i], keep[j]);
        for (std::size_t j = 0; j < Sc; ++j) Cres(Tr + i, Tr + m_tr + j) = tildeP(keep[i], j);
    }
    for (std::size_t i = 0; i < Sc; ++i)
        for (std::size_t j = 0; j < Sc; ++j)
            Cres(Tr + m_tr + i, Tr + m_tr + j) = w.wait_Smat(i, j);

    FjCodesRT2 out;
    out.RTp = fj_return_per(gamma, Cres, pers);
    out.wait_alpha = w.wait_alpha;
    out.wait_Smat = w.wait_Smat;
    out.prob_wait = w.prob_wait;
    out.residual = ct.residual;
    (void)mWait;
    return out;
}

/**
 * Port of `mainFJ.m`: the response-time percentiles of a K-node fork-join
 * queue, interpolated between the exact one-node and the approximate two-node
 * results in log K.
 *
 * @param pers the target percentiles as PROBABILITIES in (0, 1)
 * @param K    one entry per fork-join width to report
 * @param Cs   the truncation levels; only the LAST one reaches the answer, as
 *             in the reference
 */
inline std::vector<FjCodesPercentiles> fj_main(const FjDist<double>& arrival,
                                               const FjDist<double>& service,
                                               const std::vector<double>& pers,
                                               const std::vector<std::size_t>& K,
                                               const std::vector<std::size_t>& Cs, FjTMode mode) {
    using std::log;
    if (!(arrival.lambda / service.mu < 1.0))
        throw InputError("mainFJ: the system is not stable, the mean arrival rate " +
                         std::to_string(arrival.lambda) +
                         " is not below the mean service rate " + std::to_string(service.mu));
    if (Cs.empty()) throw InputError("mainFJ: no truncation level C was given");
    if (pers.empty()) throw InputError("mainFJ: no percentile was requested");

    std::vector<double> rt1, rt2;
    for (std::size_t c = 0; c < Cs.size(); ++c) {
        rt1 = fj_return_rt1(arrival, service, pers);
        rt2 = fj_return_rt2(arrival, service, pers, Cs[c], mode).RTp;
    }

    std::vector<FjCodesPercentiles> out(K.size());
    for (std::size_t k = 0; k < K.size(); ++k) {
        out[k].K = K[k];
        out[k].percentiles.resize(pers.size());
        out[k].RTp.assign(pers.size(), 0.0);
        for (std::size_t p = 0; p < pers.size(); ++p) {
            out[k].percentiles[p] = 100.0 * pers[p];
            out[k].RTp[p] = rt1[p] + (rt2[p] - rt1[p]) * log(static_cast<double>(K[k])) / log(2.0);
        }
    }
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_FJ_CODES_H
