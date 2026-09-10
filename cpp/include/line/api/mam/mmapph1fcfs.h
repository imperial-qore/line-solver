/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMAPPH1FCFS_H
#define LINE_API_MAM_MMAPPH1FCFS_H

/**
 * The MMAP[K]/PH[K]/1 FCFS queue: per-class mean number in system and per-class
 * queue-length distribution.
 *
 * This is the workhorse `solver_mam_basic.m` calls at every FCFS station, where
 * MATLAB reaches BUTools' `MMAPPH1FCFS`. The algorithm is He's age process for
 * the SM[K]/PH[K]/1/FCFS queue (Qiming He, "Analysis of a continuous time
 * SM[K]/PH[K]/1/FCFS queue: age process, sojourn times, and queue lengths",
 * Journal of Systems Science and Complexity 25(1), 133-155, 2012), which is the
 * algorithm BUTools implements and which is what makes the two agree.
 *
 * WHY THE AGE PROCESS AND NOT A QBD. A level-independent QBD over (number in
 * system, arrival phase, in-service phase) is NOT Markovian here: under FCFS
 * with class-dependent service, the law of the next service depends on the
 * class of the job at the head of the queue, so the phase would have to carry
 * the whole waiting sequence of classes. `qsys_mapmap1.h` gets away with a QBD
 * precisely because it has one class. He's construction sidesteps this by
 * tracking the AGE of the job in service against the arrival process, which
 * closes as a fluid queue whose first-return matrix Psi is exactly what
 * `mfq_fundamental` computes.
 *
 * THE PIECES, and where each comes from in this port:
 *   Psi   -- the fluid first-return matrix of the age process, ADDA doubling
 *            (`mfq_solve.h`), the same routine BUTools reaches through
 *            FluidFundamentalMatrices(..., 'P').
 *   T     -- kron(I_N, Sa) + Psi iVec, the age generator.
 *   pi0   -- the age density at zero, a linear functional of the arrival
 *            stationary vector `theta` (`ctmc_solve`) and the per-class
 *            equilibrium service vectors `beta` (also `ctmc_solve`).
 *   the queue-length recursion -- a chain of Sylvester equations
 *            T X + X kron(D0+Da-Dk, I_Ns) + C = 0 that all share their left and
 *            right coefficient matrices, so `SylvesterFactor` factors the
 *            Kronecker operator once and reuses it (`util/sylvester.h`).
 *
 * ARITHMETIC. `mfq_fundamental` runs a tolerance-terminated doubling iteration
 * and is gated on `num_traits<T>::has_transcendental`, so this function is
 * {Double, Real} and refuses by name under exact/Rational. Everything else in
 * it -- the Kronecker assembly, the linear solves, the Sylvester chain -- is
 * field arithmetic and would be exact; the Riccati root is what is not.
 *
 * MEASURED AGAINST MATLAB'S BUTools MMAPPH1FCFS: see cpp/tests/test_mam.cpp,
 * which pins both entry points on a two-class MMAP/PH/1 fixture.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/mfq_solve.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"
#include "line/util/sylvester.h"

namespace line {
namespace mam {

/** One class's phase-type service law, He's (sigma_k, S_k). */
template <class T>
struct PhService {
    std::vector<T> sigma;  ///< initial probability row vector
    Matrix<T> S;           ///< transient generator
};

namespace mmapph1_detail {

/**
 * Everything the two entry points share: the age generator T, its density at
 * zero pi0, the per-class block layout of the stacked service generator, and
 * the aggregate load rho.
 */
template <class T>
struct AgeProcess {
    Matrix<T> T_;                    ///< age generator, (N*Ns) x (N*Ns)
    std::vector<T> pi0;              ///< age density at zero, length N*Ns
    Matrix<T> D0, Da;                ///< arrival D0 and the summed arrival matrix
    std::vector<Matrix<T>> Dk;       ///< per-class arrival matrices
    std::vector<std::size_t> Nsk;    ///< per-class service order
    std::size_t N = 0, Ns = 0;       ///< arrival order, total service order
    T rho = num_traits<T>::from_int(0);
};

/** Block-diagonal stacking of the per-class service generators. */
template <class T>
Matrix<T> blkdiag(const std::vector<PhService<T>>& svc, std::size_t Ns) {
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> Sa(Ns, Ns, zero);
    std::size_t off = 0;
    for (const PhService<T>& s : svc) {
        for (std::size_t i = 0; i < s.S.rows(); ++i)
            for (std::size_t j = 0; j < s.S.cols(); ++j) Sa(off + i, off + j) = s.S(i, j);
        off += s.S.rows();
    }
    return Sa;
}

/**
 * Equilibrium (stationary) vector of the PH renewal chain, MATLAB's
 * `CTMCSolve(S - sum(S,2) sigma)`.
 *
 * THE SIGN IS THE WHOLE FUNCTION. `sum(S,2)` is the ROW SUM of a subgenerator,
 * so it is the NEGATIVE exit rate, and the generator of the renewal chain is
 * `S - rowsum sigma` = `S + exitrate sigma`. Adding the row sum instead builds a
 * matrix that is not a generator at all: on an Erlang-2 it gives
 * [[-6 6],[-6 -6]], whose "stationary" vector makes the mean service rate
 * non-positive downstream, and every Sylvester system built on it is singular.
 */
template <class T>
std::vector<T> ph_equilibrium(const PhService<T>& s) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t n = s.S.rows();
    Matrix<T> G(n, n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        T rowsum = zero;
        for (std::size_t j = 0; j < n; ++j) rowsum += s.S(i, j);
        for (std::size_t j = 0; j < n; ++j) G(i, j) = s.S(i, j) - rowsum * s.sigma[j];
    }
    return mc::ctmc_solve(G);
}

/** Build the shared age process. */
template <class T>
AgeProcess<T> build(const Mmap<T>& arrival, const std::vector<PhService<T>>& svc) {
    static_assert(num_traits<T>::has_transcendental,
                  "mmapph1fcfs runs the ADDA doubling iteration of mfq_fundamental");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t K = svc.size();
    if (arrival.classes() != K)
        throw InputError("mmapph1fcfs: the arrival MMAP and the service list disagree on the "
                         "number of classes");
    AgeProcess<T> a;
    a.N = arrival.order();
    a.D0 = arrival.D0;
    a.Dk = arrival.Dc;
    a.Da = Matrix<T>(a.N, a.N, zero);
    for (std::size_t k = 0; k < K; ++k)
        for (std::size_t i = 0; i < a.N; ++i)
            for (std::size_t j = 0; j < a.N; ++j) a.Da(i, j) += a.Dk[k](i, j);

    Matrix<T> Q(a.N, a.N, zero);
    for (std::size_t i = 0; i < a.N; ++i)
        for (std::size_t j = 0; j < a.N; ++j) Q(i, j) = a.D0(i, j) + a.Da(i, j);
    const std::vector<T> theta = mc::ctmc_solve(Q);

    a.Nsk.resize(K);
    a.Ns = 0;
    for (std::size_t k = 0; k < K; ++k) {
        a.Nsk[k] = svc[k].S.rows();
        a.Ns += a.Nsk[k];
    }
    std::vector<T> lambda(K, zero), mu(K, zero);
    std::vector<std::vector<T>> beta(K);
    a.rho = zero;
    for (std::size_t k = 0; k < K; ++k) {
        const std::vector<T> td = vecmul(theta, a.Dk[k]);
        for (const T& v : td) lambda[k] += v;
        beta[k] = ph_equilibrium(svc[k]);
        for (std::size_t i = 0; i < a.Nsk[k]; ++i)
            for (std::size_t j = 0; j < a.Nsk[k]; ++j) mu[k] -= beta[k][i] * svc[k].S(i, j);
        if (!(mu[k] > zero)) throw NumericError("mmapph1fcfs: a service law has zero rate");
        a.rho += lambda[k] / mu[k];
    }

    const Matrix<T> Sa = blkdiag(svc, a.Ns);
    const Matrix<T> Ia = eye<T>(a.N);
    const Matrix<T> Is = eye<T>(a.Ns);

    // sa{q}, ba{q} and sv{q} place class q's vectors in its own block of the
    // stacked service space and leave the other blocks at zero.
    std::vector<std::vector<T>> sa(K, std::vector<T>(a.Ns, zero));
    std::vector<std::vector<T>> ba(K, std::vector<T>(a.Ns, zero));
    std::size_t off = 0;
    for (std::size_t k = 0; k < K; ++k) {
        for (std::size_t i = 0; i < a.Nsk[k]; ++i) {
            sa[k][off + i] = svc[k].sigma[i];
            ba[k][off + i] = beta[k][i];
        }
        off += a.Nsk[k];
    }

    // Fpp = kron(Ia,Sa), Fpm = kron(Ia,-Sa e), Fmp = iVec, Fmm = D0.
    const Matrix<T> Fpp = kron(Ia, Sa);
    Matrix<T> sexit(a.Ns, 1, zero);
    for (std::size_t i = 0; i < a.Ns; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < a.Ns; ++j) s += Sa(i, j);
        sexit(i, 0) = -s;
    }
    const Matrix<T> Fpm = kron(Ia, sexit);
    Matrix<T> iVec(a.N, a.N * a.Ns, zero);
    for (std::size_t k = 0; k < K; ++k) {
        Matrix<T> row(1, a.Ns, zero);
        for (std::size_t i = 0; i < a.Ns; ++i) row(0, i) = sa[k][i];
        const Matrix<T> blk = kron(a.Dk[k], row);
        for (std::size_t i = 0; i < iVec.rows(); ++i)
            for (std::size_t j = 0; j < iVec.cols(); ++j) iVec(i, j) += blk(i, j);
    }

    const FluidFundamental<T> ff =
        mfq_fundamental(Fpp, Fpm, iVec, a.D0, T(num_traits<T>::from_double(1e-14)), 150u,
                        RiccatiMethod::ADDA);
    const Matrix<T> Y0 = ff.Psi;

    a.T_ = matmul(Y0, iVec);
    for (std::size_t i = 0; i < a.T_.rows(); ++i)
        for (std::size_t j = 0; j < a.T_.cols(); ++j) a.T_(i, j) += Fpp(i, j);

    a.pi0.assign(a.N * a.Ns, zero);
    for (std::size_t k = 0; k < K; ++k) {
        const std::vector<T> tD = vecmul(theta, a.Dk[k]);
        for (std::size_t i = 0; i < a.N; ++i)
            for (std::size_t j = 0; j < a.Ns; ++j)
                a.pi0[i * a.Ns + j] += tD[i] * T(ba[k][j] / mu[k]);
    }
    const std::vector<T> tmp = vecmul(a.pi0, a.T_);
    for (std::size_t i = 0; i < a.pi0.size(); ++i) a.pi0[i] = -tmp[i];
    (void)one;
    return a;
}

/** The (Ns x 1) indicator of class k's block, or its complement. */
template <class T>
Matrix<T> class_indicator(const AgeProcess<T>& a, std::size_t k, bool complement) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    Matrix<T> jm(a.Ns, 1, complement ? one : zero);
    std::size_t off = 0;
    for (std::size_t q = 0; q < k; ++q) off += a.Nsk[q];
    for (std::size_t i = 0; i < a.Nsk[k]; ++i) jm(off + i, 0) = complement ? zero : one;
    return jm;
}

/** kron(ones(N,1), v) for an (Ns x 1) column, i.e. the vector pi0 contracts against. */
template <class T>
std::vector<T> lift(const AgeProcess<T>& a, const Matrix<T>& jm) {
    std::vector<T> out(a.N * a.Ns);
    for (std::size_t i = 0; i < a.N; ++i)
        for (std::size_t j = 0; j < a.Ns; ++j) out[i * a.Ns + j] = jm(j, 0);
    return out;
}

/** pi0 M v for a matrix M and a column v. */
template <class T>
T contract(const std::vector<T>& pi0, const Matrix<T>& M, const std::vector<T>& v) {
    const std::vector<T> row = vecmul(pi0, M);
    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < row.size(); ++i) s += row[i] * v[i];
    return s;
}

/** The (Ns x 1) exit-rate column -S_k e placed in class k's block, He's sv{k}. */
template <class T>
Matrix<T> class_exit(const AgeProcess<T>& a, const std::vector<PhService<T>>& svc, std::size_t k) {
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> sv(a.Ns, 1, zero);
    std::size_t off = 0;
    for (std::size_t q = 0; q < k; ++q) off += a.Nsk[q];
    for (std::size_t i = 0; i < a.Nsk[k]; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < a.Nsk[k]; ++j) s += svc[k].S(i, j);
        sv(off + i, 0) = -s;
    }
    return sv;
}

}  // namespace mmapph1_detail

/** A phase-type law (alpha, A) as BUTools' `'stDistrPH'` returns it. */
template <class T>
struct StDistrPh {
    std::vector<T> alpha;
    Matrix<T> A;
};

/**
 * Per-class SOJOURN TIME as a continuous phase-type law, BUTools' `'stDistrPH'`.
 *
 * The age process is already the sojourn-time engine: `T` generates it and
 * `pi0` starts it, so the sojourn time of a class-k job is the absorption time
 * of `T` with the class-k exit column as its closing vector. What this routine
 * does beyond that is the SIMILARITY TRANSFORM BUTools applies, which turns the
 * matrix-exponential pair into a genuine PH pair: it drops the states carrying
 * no probability (`vv > precision`), rescales by `delta = diag(vv(nz))`, and
 * TRANSPOSES the generator. Without the transpose the pair still has the right
 * transform but is not a subgenerator, and `map_cdf` on it returns values
 * outside [0,1].
 *
 * The result feeds `solver_mam_passage_time`, which reads it as the MAP
 * `{A, (-A e) alpha}` and evaluates the response-time CDF on a grid.
 */
template <class T>
std::vector<StDistrPh<T>> mmapph1fcfs_stdistr_ph(const Mmap<T>& arrival,
                                                 const std::vector<PhService<T>>& svc,
                                                 double precision = 1e-14) {
    using namespace mmapph1_detail;
    const T zero = num_traits<T>::from_int(0);
    const AgeProcess<T> a = build(arrival, svc);
    const std::size_t K = svc.size(), n = a.N * a.Ns;

    Matrix<T> negT(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) negT(i, j) = -a.T_(i, j);
    const Matrix<T> iT = inverse(negT);
    const std::vector<T> vv = vecmul(a.pi0, iT);

    // The support of the age density: states outside it contribute nothing and
    // would make delta singular.
    std::vector<std::size_t> nz;
    for (std::size_t i = 0; i < n; ++i)
        if (num_traits<T>::to_double(vv[i]) > precision) nz.push_back(i);
    if (nz.empty())
        throw NumericError("mmapph1fcfs_stdistr_ph: the age density has empty support");

    std::vector<StDistrPh<T>> out(K);
    for (std::size_t k = 0; k < K; ++k) {
        const std::vector<T> clo = mulvec(iT, lift(a, class_exit(a, svc, k)));
        T norm = zero;
        for (std::size_t i = 0; i < n; ++i) norm += a.pi0[i] * clo[i];
        if (!(num_traits<T>::to_double(norm) > 0.0))
            throw NumericError("mmapph1fcfs_stdistr_ph: class " + std::to_string(k + 1) +
                               " carries no sojourn-time mass");
        // cl = -T clo / (pi0 clo)
        const std::vector<T> Tclo = mulvec(a.T_, clo);
        std::vector<T> cl(n);
        for (std::size_t i = 0; i < n; ++i) cl[i] = T(-Tclo[i] / norm);

        out[k].alpha.assign(nz.size(), zero);
        for (std::size_t i = 0; i < nz.size(); ++i) out[k].alpha[i] = T(cl[nz[i]] * vv[nz[i]]);
        out[k].A = Matrix<T>(nz.size(), nz.size(), zero);
        for (std::size_t i = 0; i < nz.size(); ++i)
            for (std::size_t j = 0; j < nz.size(); ++j)
                // inv(delta) T(nz,nz)' delta, the transpose included
                out[k].A(i, j) = T(a.T_(nz[j], nz[i]) * vv[nz[j]] / vv[nz[i]]);
    }
    return out;
}

/**
 * Per-class mean number of customers in the system, BUTools' `'ncMoms', 1`.
 *
 * @param arrival the MMAP (D0, D1, D1^(1)..D1^(K)) of the arrival stream
 * @param svc     the per-class phase-type service laws, K of them
 */
template <class T>
std::vector<T> mmapph1fcfs_ncmean(const Mmap<T>& arrival, const std::vector<PhService<T>>& svc) {
    using namespace mmapph1_detail;
    const T zero = num_traits<T>::from_int(0);
    const AgeProcess<T> a = build(arrival, svc);
    const std::size_t K = svc.size();
    const Matrix<T> Is = eye<T>(a.Ns);
    Matrix<T> QA(a.N, a.N, zero);
    for (std::size_t i = 0; i < a.N; ++i)
        for (std::size_t j = 0; j < a.N; ++j) QA(i, j) = a.D0(i, j) + a.Da(i, j);
    const SylvesterFactor<T> F(a.T_, kron(QA, Is));
    const Matrix<T> EL1 = F.solve_lyap(eye<T>(a.N * a.Ns));

    std::vector<T> out(K, zero);
    for (std::size_t k = 0; k < K; ++k) {
        // n = 1 of the moment recursion: Btag is EL1 alone, and the moment is
        // pi0 EL2 e + pi0 Btag kron(e, jm).
        const Matrix<T> EL2 = F.solve_lyap(matmul(EL1, kron(a.Dk[k], Is)));
        const std::vector<T> row = vecmul(a.pi0, EL2);
        T s = zero;
        for (const T& v : row) s += v;
        const Matrix<T> jm = class_indicator(a, k, false);
        out[k] = T(s + contract(a.pi0, EL1, lift(a, jm)));
    }
    return out;
}

/**
 * Per-class queue-length distribution, BUTools' `'ncDistr', n`: P(N_k = 0..n-1).
 *
 * @param levels the number of probabilities per class, n
 * @param arrival the marked arrival process
 * @param svc per-class phase-type service processes
 */
template <class T>
std::vector<std::vector<T>> mmapph1fcfs_ncdistr(const Mmap<T>& arrival,
                                                const std::vector<PhService<T>>& svc,
                                                std::size_t levels) {
    using namespace mmapph1_detail;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (levels == 0) throw InputError("mmapph1fcfs_ncdistr: at least one level is required");
    const AgeProcess<T> a = build(arrival, svc);
    const std::size_t K = svc.size();
    const Matrix<T> Is = eye<T>(a.Ns);
    const Matrix<T> I = eye<T>(a.N * a.Ns);

    std::vector<std::vector<T>> out(K, std::vector<T>(levels, zero));
    for (std::size_t k = 0; k < K; ++k) {
        // The coefficient of the recursion excludes class k's own arrivals: a
        // class-k arrival advances the level, every other arrival does not.
        Matrix<T> B(a.N, a.N, zero);
        for (std::size_t i = 0; i < a.N; ++i)
            for (std::size_t j = 0; j < a.N; ++j)
                B(i, j) = a.D0(i, j) + a.Da(i, j) - a.Dk[k](i, j);
        const SylvesterFactor<T> F(a.T_, kron(B, Is));
        const Matrix<T> Dkl = kron(a.Dk[k], Is);
        const std::vector<T> ejm = lift(a, class_indicator(a, k, false));
        const std::vector<T> ejmc = lift(a, class_indicator(a, k, true));

        Matrix<T> LmCurr = F.solve_lyap(I);
        out[k][0] = T(one - a.rho + contract(a.pi0, LmCurr, ejmc));
        for (std::size_t i = 1; i < levels; ++i) {
            const Matrix<T> LmPrev = LmCurr;
            LmCurr = F.solve_lyap(matmul(LmPrev, Dkl));
            out[k][i] = T(contract(a.pi0, LmCurr, ejmc) + contract(a.pi0, LmPrev, ejm));
        }
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMAPPH1FCFS_H
