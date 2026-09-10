/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MFQ_PRIO_QUEUE_H
#define LINE_API_MAM_MFQ_PRIO_QUEUE_H

/**
 * Fluid priority queue: per-class fluid level and sojourn time of an
 * MMAP[K]/PH[K]/1-type continuous fluid queue served in priority order.
 *
 * Port of matlab/src/api/mam/mfq_prio_queue.m and the BUTools FluidPrioQueue it
 * wraps, which implements G. Horvath, "Efficient analysis of the
 * MMAP[K]/PH[K]/1 priority queue", EJOR 246(1):128-139, 2015. A background
 * chain with generator Q modulates the per-class fluid input rates R (one row
 * per class) and the server drains fluid at the constant rate d, higher
 * priority fluid first.
 *
 * PRIORITY ORDER: ROW K IS THE HIGHEST PRIORITY, ROW 1 THE LOWEST. This is the
 * opposite of the obvious reading and is worth stating first, because a port
 * that assumes row 1 is highest produces plausible and wrong numbers. Two
 * things in the reference fix it: the workload seen by class k is built from
 * sum(R(k:end,:)), i.e. class k together with every class ABOVE it, and the
 * k == K branch analyses class K as if no other class existed, which only the
 * highest priority class can be.
 *
 * METHOD. For each class k, the workload of classes k..K is a fluid queue with
 * net drift diag(sum(R(k:end,:)))/d - I, solved by mfq_general_solve. That
 * gives (mass0, ini, K, clo), which is put in the canonical similarity where
 * the closing vector sums to one. For the highest priority class the answer
 * follows directly from those matrices: its sojourn time is the matrix
 * exponential law they define, and its fluid level is the ordinary fluid queue
 * with constant service rate d (the needQL half of BUTools FluFluQueue, which
 * is a handful of lines given mfq_general_solve, so it is inlined here rather
 * than ported as a separate entry point).
 *
 * For a lower priority class the server is interrupted by everything above it,
 * so the measure is a BUSY PERIOD REWARD of the fluid queue formed by stacking
 * the class-k workload on top of the background chain:
 *   F = [K clo; 0 Q],  C = [I 0; 0 diag(sum(R(k+1:end,:)))/d - I],
 *   D = [0 0; 0 I]              for sojourn time, or
 *   D = [0 0; 0 diag(R(k,:))]   for fluid level.
 * Its moments come from a recursion in the derivatives of F(v) closed by a
 * Sylvester solve at each order; its distribution comes from ERLANGIZATION,
 * replacing the deterministic horizon t by an Erlang of order L and rate L/t,
 * which converges as O(1/L) and is why erlMaxOrder is an explicit option rather
 * than a hidden constant.
 *
 * AN EMPTY UP-DRIFT SET IS AN ANSWER HERE, NOT AN ERROR. The workload drift
 * diag(sum(R(k:end,:)))/d - I is entirely negative whenever classes k..K never
 * together exceed the service rate, and then that class simply never queues:
 * every moment is zero and every distribution is one at every point. MATLAB
 * returns exactly that, and it looks like a broken reference until the drift
 * signs are checked. mfq_general_solve returns the degenerate law (no density,
 * a point mass at zero equal to the stationary distribution) rather than
 * refusing, which is what makes this path work.
 *
 * DOUBLE ONLY. mfq_general_solve is gated on transcendental arithmetic and the
 * erlangization is a tolerance-controlled truncation, but neither of those is
 * the binding constraint: this function is written on Matrix<double> for the
 * same reason mfq_multiregime is, namely that it is only ever consumed
 * alongside them, and templating it would advertise a precision tier nobody
 * has asked for and no test covers. If a Real50 instantiation is ever wanted,
 * nothing in the algorithm forbids it -- unlike mfq_multiregime, there is no
 * eigendecomposition anywhere in this path.
 *
 * References: G. Horvath, EJOR 246(1):128-139, 2015.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/mfq_multiregime.h"
#include "line/api/mam/mfq_solve.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** Which measures to compute, and the numerical options. */
struct FluidPrioOptions {
    std::size_t flMoms = 0;            ///< number of fluid level moments, 0 = not wanted
    std::vector<double> flDistr;       ///< levels at which the fluid CDF is wanted
    std::size_t stMoms = 0;            ///< number of sojourn time moments
    std::vector<double> stDistr;       ///< times at which the sojourn CDF is wanted
    double prec = 1e-14;               ///< Riccati and matrix-quadratic tolerance
    std::size_t erlMaxOrder = 200;     ///< Erlang order of the erlangization
    std::vector<std::size_t> classes;  ///< 1-based classes to analyze, empty = all
};

/** One entry per analyzed class, in the order given by FluidPrioOptions::classes. */
struct FluidPrioResult {
    std::vector<std::size_t> classes;          ///< the classes analyzed, 1-based
    std::vector<std::vector<double>> flMoms;   ///< fluid level moments per class
    std::vector<std::vector<double>> flDistr;  ///< fluid level CDF per class
    std::vector<std::vector<double>> stMoms;   ///< sojourn time moments per class
    std::vector<std::vector<double>> stDistr;  ///< sojourn time CDF per class
};

namespace prio_detail {

using multiregime_detail::blk;
using multiregime_detail::expm0;
using multiregime_detail::sylvester;

/** MATLAB lyap(A, B, C) solves A X + X B + C = 0, i.e. sylvester(A, B, -C). */
inline Matrix<double> lyap(const Matrix<double>& A, const Matrix<double>& B,
                           const Matrix<double>& C) {
    Matrix<double> negC = C;
    for (std::size_t i = 0; i < negC.rows(); ++i)
        for (std::size_t j = 0; j < negC.cols(); ++j) negC(i, j) = -negC(i, j);
    return sylvester(A, B, negC);
}

/** Diagonal matrix from a vector. */
inline Matrix<double> dg(const std::vector<double>& v) {
    Matrix<double> m(v.size(), v.size(), 0.0);
    for (std::size_t i = 0; i < v.size(); ++i) m(i, i) = v[i];
    return m;
}

/** A with every entry negated. */
inline Matrix<double> neg(const Matrix<double>& A) {
    Matrix<double> B = A;
    for (std::size_t i = 0; i < B.rows(); ++i)
        for (std::size_t j = 0; j < B.cols(); ++j) B(i, j) = -B(i, j);
    return B;
}

/** Binomial coefficient, exact for the small orders reached here. */
inline double nchoosek(std::size_t n, std::size_t k) {
    double r = 1.0;
    for (std::size_t i = 0; i < k; ++i)
        r = r * static_cast<double>(n - i) / static_cast<double>(i + 1);
    return std::floor(r + 0.5);
}

/**
 * The n-th derivative of inv(v R - Q), valid even when R has zero diagonal
 * entries. BUTools DReward: the zero-rate states are censored out, the
 * derivative is taken on the remaining block in closed form, and the result is
 * lifted back through the censoring operator.
 */
inline Matrix<double> dreward(const Matrix<double>& Q, const Matrix<double>& R, std::size_t n,
                              double prec) {
    const std::size_t NQ = Q.rows();
    std::vector<std::size_t> ixz, ixp;
    for (std::size_t i = 0; i < NQ; ++i)
        if (std::fabs(R(i, i)) <= prec) ixz.push_back(i);
    // Both signs go into the non-zero group, in the reference's order: the
    // positive rates first, then the negative ones.
    for (std::size_t i = 0; i < NQ; ++i)
        if (R(i, i) > prec) ixp.push_back(i);
    for (std::size_t i = 0; i < NQ; ++i)
        if (R(i, i) < -prec) ixp.push_back(i);
    const std::size_t Nz = ixz.size(), Np = ixp.size();

    const Matrix<double> Rp = multiregime_detail::pick(R, ixp, ixp);
    const Matrix<double> Qpp = multiregime_detail::pick(Q, ixp, ixp);
    const Matrix<double> Qpz = multiregime_detail::pick(Q, ixp, ixz);
    const Matrix<double> Qzp = multiregime_detail::pick(Q, ixz, ixp);
    const Matrix<double> Qzz = multiregime_detail::pick(Q, ixz, ixz);
    const Matrix<double> iRp = inverse(Rp);
    Matrix<double> inner = neg(Qpp);
    if (Nz > 0)
        inner = mfq_detail::sub(inner, matmul(matmul(Qpz, inverse(neg(Qzz))), Qzp));
    // dXvn = (-1)^n n! (iRp inner)^-(n+1) iRp.
    const Matrix<double> base = inverse(matmul(iRp, inner));
    Matrix<double> pw = eye<double>(Np);
    for (std::size_t i = 0; i <= n; ++i) pw = matmul(pw, base);
    double fact = 1.0;
    for (std::size_t i = 2; i <= n; ++i) fact *= static_cast<double>(i);
    const double sgn = (n % 2 == 0) ? 1.0 : -1.0;
    Matrix<double> dXvn = matmul(pw, iRp);
    for (std::size_t i = 0; i < Np; ++i)
        for (std::size_t j = 0; j < Np; ++j) dXvn(i, j) *= sgn * fact;

    // Lift back: drpar is block [[Z Qzp dX Qpz Z, Z Qzp dX], [dX Qpz Z, dX]]
    // in the (z, p) ordering, then permuted to the original state order.
    Matrix<double> out(NQ, NQ, 0.0);
    if (Nz > 0) {
        const Matrix<double> Z = inverse(neg(Qzz));
        const Matrix<double> ZQzp = matmul(Z, Qzp);
        const Matrix<double> QpzZ = matmul(Qpz, Z);
        const Matrix<double> tl = matmul(matmul(ZQzp, dXvn), QpzZ);
        const Matrix<double> tr = matmul(ZQzp, dXvn);
        const Matrix<double> bl = matmul(dXvn, QpzZ);
        for (std::size_t i = 0; i < Nz; ++i) {
            for (std::size_t j = 0; j < Nz; ++j) out(ixz[i], ixz[j]) = tl(i, j);
            for (std::size_t j = 0; j < Np; ++j) out(ixz[i], ixp[j]) = tr(i, j);
        }
        for (std::size_t i = 0; i < Np; ++i)
            for (std::size_t j = 0; j < Nz; ++j) out(ixp[i], ixz[j]) = bl(i, j);
    }
    for (std::size_t i = 0; i < Np; ++i)
        for (std::size_t j = 0; j < Np; ++j) out(ixp[i], ixp[j]) = dXvn(i, j);
    return out;
}

/** The (z, p, n) partition of a fluid model by the sign of diag(C). */
struct SignSplit {
    std::vector<std::size_t> ixz, ixp, ixn;
};

inline SignSplit sign_split(const Matrix<double>& C, double prec) {
    SignSplit s;
    for (std::size_t i = 0; i < C.rows(); ++i) {
        if (std::fabs(C(i, i)) <= prec) s.ixz.push_back(i);
        else if (C(i, i) > prec) s.ixp.push_back(i);
        else s.ixn.push_back(i);
    }
    return s;
}

/**
 * Embed an Np x Nn block back into the full NF x NF state space at the (p, n)
 * position, undoing the sign permutation. This is the reference's
 * iPer [0; 0 X; 0] Per.
 */
inline Matrix<double> embed_pn(const Matrix<double>& X, const SignSplit& s, std::size_t NF) {
    Matrix<double> out(NF, NF, 0.0);
    for (std::size_t i = 0; i < s.ixp.size(); ++i)
        for (std::size_t j = 0; j < s.ixn.size(); ++j) out(s.ixp[i], s.ixn[j]) = X(i, j);
    return out;
}

/**
 * Moments of the busy period reward of the fluid model (F, C, D), BUTools
 * BusyPeriodRewardMoms. Returns numOfMoms + 1 matrices, the first being Psi.
 */
inline std::vector<Matrix<double>> busy_period_reward_moms(const Matrix<double>& F,
                                                           const Matrix<double>& C,
                                                           const Matrix<double>& D,
                                                           std::size_t numOfMoms, double prec) {
    const std::size_t NF = F.rows();
    const SignSplit s = sign_split(C, prec);
    const std::size_t Nz = s.ixz.size(), Np = s.ixp.size(), Nn = s.ixn.size();

    const Matrix<double> Fzz = multiregime_detail::pick(F, s.ixz, s.ixz);
    const Matrix<double> Fpz = multiregime_detail::pick(F, s.ixp, s.ixz);
    const Matrix<double> Fmz = multiregime_detail::pick(F, s.ixn, s.ixz);
    const Matrix<double> Fzp = multiregime_detail::pick(F, s.ixz, s.ixp);
    const Matrix<double> Fpp = multiregime_detail::pick(F, s.ixp, s.ixp);
    const Matrix<double> Fmp = multiregime_detail::pick(F, s.ixn, s.ixp);
    const Matrix<double> Fzm = multiregime_detail::pick(F, s.ixz, s.ixn);
    const Matrix<double> Fpm = multiregime_detail::pick(F, s.ixp, s.ixn);
    const Matrix<double> Fmm = multiregime_detail::pick(F, s.ixn, s.ixn);
    const Matrix<double> Cm = multiregime_detail::pick(C, s.ixn, s.ixn);
    const Matrix<double> Cp = multiregime_detail::pick(C, s.ixp, s.ixp);
    const Matrix<double> Dm = multiregime_detail::pick(D, s.ixn, s.ixn);
    const Matrix<double> Dp = multiregime_detail::pick(D, s.ixp, s.ixp);
    const Matrix<double> Dz = multiregime_detail::pick(D, s.ixz, s.ixz);
    const Matrix<double> iCp = inverse(Cp);
    const Matrix<double> iCm = inverse(neg(Cm));

    // Zeroth derivatives, censoring the zero-drift block.
    Matrix<double> Zf(Nz, Nz, 0.0);
    if (Nz > 0) Zf = inverse(neg(Fzz));
    auto cens = [&](const Matrix<double>& A, const Matrix<double>& Az,
                    const Matrix<double>& Zb) -> Matrix<double> {
        if (Nz == 0) return A;
        return mfq_detail::add(A, matmul(matmul(Az, Zf), Zb));
    };
    std::vector<Matrix<double>> Fppd(numOfMoms + 1), Fpmd(numOfMoms + 1), Fmpd(numOfMoms + 1),
        Fmmd(numOfMoms + 1);
    Fppd[0] = matmul(iCp, cens(Fpp, Fpz, Fzp));
    Fpmd[0] = matmul(iCp, cens(Fpm, Fpz, Fzm));
    Fmpd[0] = matmul(iCm, cens(Fmp, Fmz, Fzp));
    Fmmd[0] = matmul(iCm, cens(Fmm, Fmz, Fzm));
    for (std::size_t i = 1; i <= numOfMoms; ++i) {
        const Matrix<double> dr = (Nz > 0) ? dreward(Fzz, Dz, i, prec) : Matrix<double>(0, 0, 0.0);
        Fppd[i] = matmul(matmul(matmul(iCp, Fpz), dr), Fzp);
        Fpmd[i] = matmul(matmul(matmul(iCp, Fpz), dr), Fzm);
        Fmpd[i] = matmul(matmul(matmul(iCm, Fmz), dr), Fzp);
        Fmmd[i] = matmul(matmul(matmul(iCm, Fmz), dr), Fzm);
        if (i == 1) {
            Fppd[i] = mfq_detail::sub(Fppd[i], matmul(iCp, Dp));
            Fmmd[i] = mfq_detail::sub(Fmmd[i], matmul(iCm, Dm));
        }
    }

    const FluidFundamental<double> ff =
        mfq_fundamental(Fppd[0], Fpmd[0], Fmpd[0], Fmmd[0], prec, 150u, RiccatiMethod::ADDA);
    const Matrix<double>& Psi = ff.Psi;

    std::vector<Matrix<double>> BPM(numOfMoms + 1);
    BPM[0] = Psi;
    for (std::size_t i = 1; i <= numOfMoms; ++i) {
        Matrix<double> X = mfq_detail::add(neg(matmul(matmul(Psi, Fmpd[i]), Psi)), Fpmd[i]);
        for (std::size_t m = 0; m + 1 <= i; ++m) {
            const double c = nchoosek(i, m);
            Matrix<double> t = mfq_detail::add(
                matmul(mfq_detail::add(Fppd[i - m], matmul(Psi, Fmpd[i - m])), BPM[m]),
                matmul(BPM[m], mfq_detail::add(Fmmd[i - m], matmul(Fmpd[i - m], Psi))));
            X = mfq_detail::add(X, mfq_detail::scale(t, c));
        }
        for (std::size_t l = 1; l + 1 <= i; ++l)
            for (std::size_t m = 1; m + l <= i; ++m) {
                const double c = nchoosek(i, l) * nchoosek(i - l, m);
                const Matrix<double> t = matmul(matmul(BPM[l], Fmpd[i - l - m]), BPM[m]);
                X = mfq_detail::add(X, mfq_detail::scale(t, c));
            }
        BPM[i] = lyap(mfq_detail::add(Fppd[0], matmul(Psi, Fmpd[0])),
                      mfq_detail::add(Fmmd[0], matmul(Fmpd[0], Psi)), X);
    }
    for (std::size_t i = 0; i <= numOfMoms; ++i) BPM[i] = embed_pn(BPM[i], s, NF);
    (void)Np;
    (void)Nn;
    return BPM;
}

/** Result of the erlangized busy period reward distribution. */
struct BusyPeriodDistr {
    Matrix<double> pr;             ///< accumulated, embedded in the full state space
    std::vector<Matrix<double>> Pn;  ///< the per-order terms, likewise embedded
};

/**
 * Distribution of the busy period reward at horizon t, BUTools
 * BusyPeriodRewardDistr, by erlangization: the deterministic horizon t is
 * replaced by an Erlang of order L = erlMaxOrder and rate nu = L/t, which
 * converges as O(1/L). L is an option and not a constant precisely because it
 * is the accuracy knob.
 */
inline BusyPeriodDistr busy_period_reward_distr(const Matrix<double>& F, const Matrix<double>& C,
                                                const Matrix<double>& D, double t,
                                                std::size_t L, double prec) {
    const std::size_t NF = F.rows();
    const SignSplit s = sign_split(C, prec);
    const std::size_t Nz = s.ixz.size();

    const Matrix<double> Fzz = multiregime_detail::pick(F, s.ixz, s.ixz);
    const Matrix<double> Fpz = multiregime_detail::pick(F, s.ixp, s.ixz);
    const Matrix<double> Fmz = multiregime_detail::pick(F, s.ixn, s.ixz);
    const Matrix<double> Fzp = multiregime_detail::pick(F, s.ixz, s.ixp);
    const Matrix<double> Fpp = multiregime_detail::pick(F, s.ixp, s.ixp);
    const Matrix<double> Fmp = multiregime_detail::pick(F, s.ixn, s.ixp);
    const Matrix<double> Fzm = multiregime_detail::pick(F, s.ixz, s.ixn);
    const Matrix<double> Fpm = multiregime_detail::pick(F, s.ixp, s.ixn);
    const Matrix<double> Fmm = multiregime_detail::pick(F, s.ixn, s.ixn);
    const Matrix<double> Cm = multiregime_detail::pick(C, s.ixn, s.ixn);
    const Matrix<double> Cp = multiregime_detail::pick(C, s.ixp, s.ixp);
    const Matrix<double> Dm = multiregime_detail::pick(D, s.ixn, s.ixn);
    const Matrix<double> Dp = multiregime_detail::pick(D, s.ixp, s.ixp);
    const Matrix<double> Dz = multiregime_detail::pick(D, s.ixz, s.ixz);
    const Matrix<double> iCp = inverse(Cp);
    const Matrix<double> iCm = inverse(neg(Cm));

    const double nu = static_cast<double>(L) / t;
    Matrix<double> Z(Nz, Nz, 0.0);
    if (Nz > 0) Z = inverse(mfq_detail::sub(mfq_detail::scale(Dz, nu), Fzz));
    auto zc = [&](const Matrix<double>& A, const Matrix<double>& Az,
                  const Matrix<double>& Zb) -> Matrix<double> {
        if (Nz == 0) return A;
        return mfq_detail::add(A, matmul(matmul(Az, Z), Zb));
    };
    const Matrix<double> Fpp_e =
        matmul(iCp, zc(mfq_detail::sub(Fpp, mfq_detail::scale(Dp, nu)), Fpz, Fzp));
    const Matrix<double> Fpm_e = matmul(iCp, zc(Fpm, Fpz, Fzm));
    const Matrix<double> Fmp_e = matmul(iCm, zc(Fmp, Fmz, Fzp));
    const Matrix<double> Fmm_e =
        matmul(iCm, zc(mfq_detail::sub(Fmm, mfq_detail::scale(Dm, nu)), Fmz, Fzm));
    const Matrix<double> Psie =
        mfq_fundamental(Fpp_e, Fpm_e, Fmp_e, Fmm_e, prec, 150u, RiccatiMethod::ADDA).Psi;

    std::vector<Matrix<double>> Pn;
    Pn.push_back(Psie);
    Matrix<double> pr = Psie;
    const Matrix<double> AM = mfq_detail::add(Fpp_e, matmul(Psie, Fmp_e));
    const Matrix<double> BM = mfq_detail::add(Fmm_e, matmul(Fmp_e, Psie));
    const Matrix<double> iCpDp = mfq_detail::scale(matmul(iCp, Dp), nu);
    const Matrix<double> iCmDm = mfq_detail::scale(matmul(iCm, Dm), nu);
    const Matrix<double> iCmFmp = matmul(iCm, Fmp);
    const Matrix<double> iCpFpz = (Nz > 0) ? matmul(iCp, Fpz) : Matrix<double>(Fpp.rows(), 0, 0.0);
    const Matrix<double> iCmFmz = (Nz > 0) ? matmul(iCm, Fmz) : Matrix<double>(Fmm.rows(), 0, 0.0);
    const Matrix<double> nuDzZ = (Nz > 0) ? matmul(mfq_detail::scale(Dz, nu), Z)
                                          : Matrix<double>(0, 0, 0.0);

    for (std::size_t n = 1; n + 1 <= L; ++n) {
        Matrix<double> CM = mfq_detail::add(matmul(iCpDp, Pn[n - 1]), matmul(Pn[n - 1], iCmDm));
        for (std::size_t i = 1; i + 1 <= n; ++i)
            CM = mfq_detail::add(CM, matmul(matmul(Pn[i], iCmFmp), Pn[n - i]));
        if (Nz > 0) {
            // Z (nu Dz Z)^n, built once per n by one extra multiplication.
            Matrix<double> ZnuN = Z;
            for (std::size_t i = 0; i < n; ++i) ZnuN = matmul(ZnuN, nuDzZ);
            CM = mfq_detail::add(CM, matmul(matmul(iCpFpz, ZnuN), Fzm));
            CM = mfq_detail::sub(
                CM, matmul(matmul(matmul(matmul(Psie, iCmFmz), ZnuN), Fzp), Psie));
            for (std::size_t i = 0; i + 1 <= n; ++i) {
                Matrix<double> Zk = Z;
                for (std::size_t q = 0; q < n - i; ++q) Zk = matmul(Zk, nuDzZ);
                const Matrix<double> tail =
                    mfq_detail::add(Fzm, matmul(Fzp, Psie));
                CM = mfq_detail::add(CM, matmul(matmul(matmul(Pn[i], iCmFmz), Zk), tail));
                CM = mfq_detail::add(
                    CM, matmul(matmul(matmul(mfq_detail::add(iCpFpz, matmul(Psie, iCmFmz)), Zk),
                                      Fzp),
                               Pn[i]));
            }
            for (std::size_t i = 1; i + 1 <= n; ++i)
                for (std::size_t j = 1; j + i <= n; ++j) {
                    Matrix<double> Zk = Z;
                    for (std::size_t q = 0; q < n - i - j; ++q) Zk = matmul(Zk, nuDzZ);
                    CM = mfq_detail::add(
                        CM, matmul(matmul(matmul(matmul(Pn[i], iCmFmz), Zk), Fzp), Pn[j]));
                }
        }
        const Matrix<double> PM = lyap(AM, BM, CM);
        Pn.push_back(PM);
        pr = mfq_detail::add(pr, PM);
    }

    BusyPeriodDistr out;
    out.pr = embed_pn(pr, s, NF);
    out.Pn.reserve(Pn.size());
    for (const Matrix<double>& P : Pn) out.Pn.push_back(embed_pn(P, s, NF));
    return out;
}

}  // namespace prio_detail

/**
 * Fluid priority queue.
 *
 * @param Q   generator of the modulating chain, N x N
 * @param R   per-class fluid input rates, K x N; ROW K IS THE HIGHEST PRIORITY
 * @param d   constant fluid service rate, positive
 * @param opt which measures to compute and the numerical options
 */
inline FluidPrioResult mfq_prio_queue(const Matrix<double>& Q, const Matrix<double>& R, double d,
                                      const FluidPrioOptions& opt) {
    using namespace prio_detail;
    const std::size_t N = Q.rows();
    const std::size_t K = R.rows();
    if (Q.cols() != N) throw InputError("mfq_prio_queue: Q must be square");
    if (R.cols() != N)
        throw InputError("mfq_prio_queue: R must have one column per background state");
    if (K == 0) throw InputError("mfq_prio_queue: at least one fluid class is required");
    if (d <= opt.prec) throw InputError("mfq_prio_queue: the fluid service rate must be positive");
    for (std::size_t k = 0; k < K; ++k)
        for (std::size_t j = 0; j < N; ++j)
            if (R(k, j) < -opt.prec)
                throw InputError("mfq_prio_queue: the fluid arrival rate cannot be negative");
    if (opt.erlMaxOrder < 2)
        throw InputError("mfq_prio_queue: erlMaxOrder must be at least 2");

    std::vector<std::size_t> classes = opt.classes;
    if (classes.empty())
        for (std::size_t k = 1; k <= K; ++k) classes.push_back(k);
    for (std::size_t c : classes)
        if (c < 1 || c > K) throw InputError("mfq_prio_queue: class index out of range");

    const std::vector<double> pi = mc::ctmc_solve(Q);
    std::vector<double> lambda(K, 0.0);
    for (std::size_t k = 0; k < K; ++k)
        for (std::size_t j = 0; j < N; ++j) lambda[k] += pi[j] * R(k, j);

    FluidPrioResult out;
    out.classes = classes;
    for (std::size_t ci = 0; ci < classes.size(); ++ci) {
        const std::size_t k = classes[ci] - 1;  // 0-based
        if (lambda[k] <= 0.0)
            throw InputError("mfq_prio_queue: a requested class has zero mean input rate");

        // Workload of classes k..K-1 (0-based), i.e. class k and every class
        // above it in priority.
        std::vector<double> agg(N, 0.0);
        for (std::size_t kk = k; kk < K; ++kk)
            for (std::size_t j = 0; j < N; ++j) agg[j] += R(kk, j);
        std::vector<double> drift(N, 0.0);
        for (std::size_t j = 0; j < N; ++j) drift[j] = agg[j] / d - 1.0;
        const GeneralFluidSolution<double> gs =
            mfq_general_solve(Q, dg(drift), Matrix<double>(0, 0), opt.prec);
        const std::size_t KN = gs.K.rows();

        // clok, and the canonical similarity in which the closing vector sums
        // to one; kappa is the initial vector in that basis.
        Matrix<double> clok = matmul(gs.clo, dg(std::vector<double>(N, 1.0)));
        for (std::size_t i = 0; i < KN; ++i)
            for (std::size_t j = 0; j < N; ++j) clok(i, j) = gs.clo(i, j) * R(k, j) / lambda[k];
        std::vector<double> ini = gs.ini;
        Matrix<double> Km = gs.K;
        if (KN > 0) {
            const Matrix<double> iKm = inverse(neg(Km));
            std::vector<double> rowsum(KN, 0.0);
            for (std::size_t i = 0; i < KN; ++i)
                for (std::size_t j = 0; j < N; ++j) rowsum[i] += clok(i, j);
            const std::vector<double> delta = mulvec(iKm, rowsum);
            for (std::size_t i = 0; i < KN; ++i)
                if (delta[i] == 0.0)
                    throw NumericError(
                        "mfq_prio_queue: the canonical similarity is singular for this class");
            Matrix<double> K0(KN, KN, 0.0);
            for (std::size_t i = 0; i < KN; ++i)
                for (std::size_t j = 0; j < KN; ++j) K0(i, j) = Km(i, j) * delta[j] / delta[i];
            Matrix<double> K1(KN, N, 0.0);
            for (std::size_t i = 0; i < KN; ++i)
                for (std::size_t j = 0; j < N; ++j) K1(i, j) = clok(i, j) / delta[i];
            std::vector<double> kappa(KN, 0.0);
            for (std::size_t i = 0; i < KN; ++i) kappa[i] = ini[i] * delta[i];
            Km = K0;
            clok = K1;
            ini = kappa;
        }

        // mass0 weighted by this class's input, the probability that a drop
        // arrives to an empty higher-priority workload.
        double mass0R = 0.0;
        for (std::size_t j = 0; j < N; ++j) mass0R += gs.mass0[j] * R(k, j) / lambda[k];

        const bool highest = (k + 1 == K);
        if (highest) {
            // ---- highest priority class: no interruption ----
            if (opt.stMoms > 0) {
                std::vector<double> m(opt.stMoms, 0.0);
                if (KN > 0) {
                    const Matrix<double> iK = inverse(neg(Km));
                    Matrix<double> pw = iK;
                    double fact = 1.0;
                    for (std::size_t i = 1; i <= opt.stMoms; ++i) {
                        fact *= static_cast<double>(i);
                        pw = matmul(pw, iK);
                        const std::vector<double> v = vecmul(vecmul(ini, pw), clok);
                        double sum = 0.0;
                        for (double x : v) sum += x;
                        m[i - 1] = fact * sum;
                    }
                }
                out.stMoms.push_back(m);
            }
            if (!opt.stDistr.empty()) {
                std::vector<double> v(opt.stDistr.size(), 0.0);
                for (std::size_t x = 0; x < opt.stDistr.size(); ++x) {
                    double acc = mass0R;
                    if (KN > 0) {
                        const Matrix<double> iK = inverse(neg(Km));
                        Matrix<double> ImE = eye<double>(KN);
                        const Matrix<double> E = expm0(Km, opt.stDistr[x]);
                        for (std::size_t i = 0; i < KN; ++i)
                            for (std::size_t j = 0; j < KN; ++j) ImE(i, j) -= E(i, j);
                        const std::vector<double> w =
                            vecmul(vecmul(vecmul(ini, iK), ImE), clok);
                        for (double y : w) acc += y;
                    }
                    v[x] = acc;
                }
                out.stDistr.push_back(v);
            }
            if (opt.flMoms > 0 || !opt.flDistr.empty()) {
                // needQL/FluFluQueue rationale: see _kb/03-api-layer.md (cpp port notes: mam)
                std::vector<double> dr(N, 0.0);
                for (std::size_t j = 0; j < N; ++j) dr[j] = R(k, j) - d;
                const GeneralFluidSolution<double> fs =
                    mfq_general_solve(Q, dg(dr), Q, opt.prec);
                const std::size_t FN = fs.K.rows();
                if (opt.flMoms > 0) {
                    std::vector<double> m(opt.flMoms, 0.0);
                    if (FN > 0) {
                        const Matrix<double> iK = inverse(neg(fs.K));
                        Matrix<double> pw = iK;
                        double fact = 1.0;
                        for (std::size_t i = 1; i <= opt.flMoms; ++i) {
                            fact *= static_cast<double>(i);
                            pw = matmul(pw, iK);
                            const std::vector<double> v = vecmul(vecmul(fs.ini, pw), fs.clo);
                            double sum = 0.0;
                            for (double x : v) sum += x;
                            m[i - 1] = fact * sum;
                        }
                    }
                    out.flMoms.push_back(m);
                }
                if (!opt.flDistr.empty()) {
                    std::vector<double> v(opt.flDistr.size(), 0.0);
                    double m0 = 0.0;
                    for (double x : fs.mass0) m0 += x;
                    for (std::size_t x = 0; x < opt.flDistr.size(); ++x) {
                        double acc = m0;
                        if (FN > 0) {
                            const Matrix<double> iK = inverse(neg(fs.K));
                            Matrix<double> ImE = eye<double>(FN);
                            const Matrix<double> E = expm0(fs.K, opt.flDistr[x]);
                            for (std::size_t i = 0; i < FN; ++i)
                                for (std::size_t j = 0; j < FN; ++j) ImE(i, j) -= E(i, j);
                            const std::vector<double> w =
                                vecmul(vecmul(vecmul(fs.ini, ImE), iK), fs.clo);
                            for (double y : w) acc += y;
                        }
                        v[x] = acc;
                    }
                    out.flDistr.push_back(v);
                }
            }
            continue;
        }

        // busy period reward rationale: see _kb/03-api-layer.md (cpp port notes: mam)
        std::vector<double> lower(N, 0.0);
        for (std::size_t kk = k + 1; kk < K; ++kk)
            for (std::size_t j = 0; j < N; ++j) lower[j] += R(kk, j);
        const std::size_t NF = KN + N;
        Matrix<double> F(NF, NF, 0.0), Cmat(NF, NF, 0.0);
        for (std::size_t i = 0; i < KN; ++i) {
            for (std::size_t j = 0; j < KN; ++j) F(i, j) = Km(i, j);
            for (std::size_t j = 0; j < N; ++j) F(i, KN + j) = clok(i, j);
            Cmat(i, i) = 1.0;
        }
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) F(KN + i, KN + j) = Q(i, j);
            Cmat(KN + i, KN + i) = lower[i] / d - 1.0;
        }
        std::vector<double> inis(NF, 0.0);
        for (std::size_t i = 0; i < KN; ++i) inis[i] = ini[i];

        auto make_D = [&](bool fluid) {
            Matrix<double> D(NF, NF, 0.0);
            for (std::size_t i = 0; i < N; ++i)
                D(KN + i, KN + i) = fluid ? R(k, i) : 1.0;
            return D;
        };

        if (opt.stMoms > 0) {
            const std::vector<Matrix<double>> Tmp =
                busy_period_reward_moms(F, Cmat, make_D(false), opt.stMoms, opt.prec);
            std::vector<double> m(opt.stMoms, 0.0);
            for (std::size_t i = 1; i < Tmp.size(); ++i) {
                const std::vector<double> v = vecmul(inis, Tmp[i]);
                double sum = 0.0;
                for (double x : v) sum += x;
                m[i - 1] = ((i % 2 == 0) ? 1.0 : -1.0) * sum;
            }
            out.stMoms.push_back(m);
        }
        if (!opt.stDistr.empty()) {
            std::vector<double> v(opt.stDistr.size(), 0.0);
            for (std::size_t x = 0; x < opt.stDistr.size(); ++x) {
                const BusyPeriodDistr bp = busy_period_reward_distr(
                    F, Cmat, make_D(false), opt.stDistr[x], opt.erlMaxOrder, opt.prec);
                const std::vector<double> w = vecmul(inis, bp.pr);
                double acc = mass0R;
                for (double y : w) acc += y;
                v[x] = acc;
            }
            out.stDistr.push_back(v);
        }
        if (opt.flMoms > 0) {
            const std::vector<Matrix<double>> Tmp =
                busy_period_reward_moms(F, Cmat, make_D(true), opt.flMoms, opt.prec);
            // Keep only the background columns, then move from the law right
            // after a departure to the law at a random point in time.
            std::vector<std::vector<double>> FLDPn(Tmp.size(), std::vector<double>(N, 0.0));
            for (std::size_t i = 0; i < Tmp.size(); ++i) {
                const std::vector<double> row = vecmul(inis, Tmp[i]);
                for (std::size_t j = 0; j < N; ++j) FLDPn[i][j] = row[KN + j];
                if (i == 0)
                    for (std::size_t j = 0; j < N; ++j)
                        FLDPn[i][j] += gs.mass0[j] * R(k, j) / lambda[k];
            }
            // iTerm = inv(ones(N,1) pi - Q), the fundamental matrix of Q.
            Matrix<double> T2(N, N, 0.0);
            for (std::size_t i = 0; i < N; ++i)
                for (std::size_t j = 0; j < N; ++j) T2(i, j) = pi[j] - Q(i, j);
            const Matrix<double> iTerm = inverse(T2);
            std::vector<std::vector<double>> FLPn;
            FLPn.push_back(pi);
            std::vector<double> m(opt.flMoms, 0.0);
            for (std::size_t n = 1; n <= opt.flMoms; ++n) {
                const double nd = static_cast<double>(n);
                std::vector<double> a(N, 0.0);
                for (std::size_t j = 0; j < N; ++j)
                    a[j] = -FLDPn[n - 1][j] + FLPn[n - 1][j] * R(k, j) / lambda[k];
                const std::vector<double> aT = vecmul(a, iTerm);
                double sumP = 0.0;
                for (std::size_t j = 0; j < N; ++j) sumP += FLDPn[n][j];
                for (std::size_t j = 0; j < N; ++j) sumP += nd * aT[j] * R(k, j);
                std::vector<double> b(N, 0.0);
                for (std::size_t j = 0; j < N; ++j)
                    b[j] = -FLPn[n - 1][j] * R(k, j) + FLDPn[n - 1][j] * lambda[k];
                const std::vector<double> bT = vecmul(b, iTerm);
                std::vector<double> P(N, 0.0);
                for (std::size_t j = 0; j < N; ++j) P[j] = sumP * pi[j] + nd * bT[j];
                FLPn.push_back(P);
                double sum = 0.0;
                for (double x : P) sum += x;
                m[n - 1] = ((n % 2 == 0) ? 1.0 : -1.0) * sum;
            }
            out.flMoms.push_back(m);
        }
        if (!opt.flDistr.empty()) {
            std::vector<double> v(opt.flDistr.size(), 0.0);
            for (std::size_t x = 0; x < opt.flDistr.size(); ++x) {
                const double nu = static_cast<double>(opt.erlMaxOrder) / opt.flDistr[x];
                const BusyPeriodDistr bp = busy_period_reward_distr(
                    F, Cmat, make_D(true), opt.flDistr[x], opt.erlMaxOrder, opt.prec);
                Matrix<double> Wm(N, N, 0.0);
                for (std::size_t i = 0; i < N; ++i)
                    for (std::size_t j = 0; j < N; ++j)
                        Wm(i, j) = (i == j ? nu * R(k, i) : 0.0) - Q(i, j);
                const Matrix<double> iW = inverse(Wm);
                std::vector<double> Psiy(N, 0.0);
                {
                    const std::vector<double> row = vecmul(inis, bp.Pn[0]);
                    std::vector<double> a(N, 0.0);
                    for (std::size_t j = 0; j < N; ++j)
                        a[j] = gs.mass0[j] * R(k, j) / lambda[k] + row[KN + j];
                    const std::vector<double> t = vecmul(a, iW);
                    for (std::size_t j = 0; j < N; ++j) Psiy[j] = lambda[k] * nu * t[j];
                }
                for (std::size_t i = 1; i < bp.Pn.size(); ++i) {
                    const std::vector<double> row = vecmul(inis, bp.Pn[i]);
                    std::vector<double> a(N, 0.0);
                    for (std::size_t j = 0; j < N; ++j)
                        a[j] = lambda[k] * row[KN + j] + Psiy[j] * R(k, j);
                    const std::vector<double> t = vecmul(a, iW);
                    for (std::size_t j = 0; j < N; ++j) Psiy[j] = nu * t[j];
                }
                double acc = 0.0;
                for (double y : Psiy) acc += y;
                v[x] = acc;
            }
            out.flDistr.push_back(v);
        }
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MFQ_PRIO_QUEUE_H
