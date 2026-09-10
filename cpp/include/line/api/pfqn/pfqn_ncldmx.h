/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_NCLDMX_H
#define LINE_API_PFQN_NCLDMX_H

/**
 * Normalizing constant of a MIXED open/closed network with limited load
 * dependence.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_ncldmx.m.
 *
 * The closed-conditional constant of a mixed limited-load-dependent network is
 * the constant of a purely CLOSED load-dependent network in which every
 * station carries the Bruell-Balbo-Afshari effective-capacity rate
 *
 *   mu_i^eff(n) = 1 / EC_i(n),
 *
 * with EC from pfqn_ldmx_ec folding the open classes into the closed
 * subnetwork. The open classes contribute the separable prefactor
 *
 *   Gopen = prod_i E_i(0),
 *
 * which reduces to prod_i 1/(1 - rho_i) in the load-independent limit. The two
 * are returned separately, as in the reference, because callers use the closed
 * one for the closed-class ratios G(N - e_r)/G(N) and the open one only for
 * the joint state probabilities.
 *
 * OPEN CLASSES are marked by a NEGATIVE population, as everywhere in this port
 * (see pfqn_nc); MATLAB uses Inf, which has no counterpart in an exact field.
 *
 * Arithmetic: EXACT-CAPABLE. pfqn_ldmx_ec and pfqn_ncld are both exact, and
 * everything this routine adds is a reciprocal and a product. Gopen is
 * returned as a value of T, not only as its logarithm, which is what lets the
 * mixed constant stay in the rational field end to end.
 *
 * MEAN MEASURES come out of the same identification, without ever enumerating
 * the closed population lattice -- which is what makes this the
 * normalizing-constant counterpart of pfqn_mvaldmx rather than a rename of it:
 *
 *   - closed throughputs are the ratios X_r = G(N - e_r)/G(N);
 *   - closed queue lengths are the conditional normalizing-constant recursion
 *     of the load-dependent closed network (pfqn_mushift / pfqn_fnc), applied
 *     to the effective-capacity rates;
 *   - open queue lengths are the Bruell-Balbo-Afshari sum
 *     Q_ir = lambda_r D_ir sum_n (n+1) EC_i(n+1) P_i(n) with its SATURATED
 *     TAIL FOLDED ONTO THE CLOSED MEAN. EC_i(n) is constant for n >= b_i, the
 *     level where the rate row stops growing, so writing
 *     EC_i(n) = EC_i^inf + delta_i(n) with delta_i(n) = 0 for n >= b_i leaves
 *
 *       Q_ir = lambda_r D_ir [ EC_i^inf (Q_i^closed + 1)
 *                              + sum_{n=0}^{b_i-2} (n+1) delta_i(n+1) P_i(n) ]
 *
 *     using sum_n P_i(n) = 1 and sum_n n P_i(n) = Q_i^closed. Only the first
 *     b_i-1 marginals survive, and b_i is the SERVER COUNT, not the population:
 *     a single-server station needs none at all and the formula collapses to
 *     the classical lambda_r D_ir (1 + Q_i^closed)/(1 - rho_i).
 *
 * The ratios themselves are taken in the log domain and converted back with
 * from_double, exactly as solver_ncld's closed branch does: a mean measure is a
 * ratio of constants, not a constant, so it leaves the exact field there.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_fnc.h"
#include "line/api/pfqn/pfqn_ldmx_ec.h"
#include "line/api/pfqn/pfqn_mushift.h"
#include "line/api/pfqn/pfqn_ncld.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

template <class T>
struct NcldmxResult {
    T G;          ///< closed-conditional normalizing constant
    double lG;    ///< its logarithm
    T Gopen;      ///< the open-class prefactor prod_i E_i(0)
    double lGopen;///< its logarithm
    Matrix<T> EC; ///< the effective-capacity terms, for the caller's mean measures
    std::vector<T> XN; ///< (R) throughputs: G(N-e_r)/G(N) closed, lambda_r open
    Matrix<T> QN;      ///< (M x R) mean queue lengths
};

namespace detail {

/**
 * pfqn_ncld, with the empty-station residual network handled explicitly: a
 * network reduced to its think times alone has G(N) = prod_r Z_r^N_r / N_r!,
 * and no G at all when a class has jobs but neither demand nor think time.
 */
template <class T>
double ncldmx_lg(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                 const Matrix<T>& mu, NcldMethod method, const T& atol, const NcOptions& nopt) {
    if (L.rows() == 0) {
        bool empty = true;
        for (std::size_t r = 0; r < N.size(); ++r)
            if (N[r] > 0) { empty = false; break; }
        if (empty) return 0.0;
        double acc = 0.0;
        for (std::size_t r = 0; r < N.size(); ++r) {
            if (N[r] <= 0) continue;
            double zr = 0.0;
            for (std::size_t i = 0; i < Z.rows(); ++i) zr += num_traits<T>::to_double(Z(i, r));
            if (!(zr > 0.0)) return -std::numeric_limits<double>::infinity();
            acc += N[r] * std::log(zr) - std::lgamma(static_cast<double>(N[r]) + 1.0);
        }
        return acc;
    }
    return pfqn_ncld(L, N, Z, mu, method, atol, nopt).lG;
}

/** A copy of A with row `skip` removed. */
template <class T>
Matrix<T> ncldmx_drop_row(const Matrix<T>& A, std::size_t skip) {
    const std::size_t rows = A.rows() == 0 ? 0 : A.rows() - 1;
    Matrix<T> out(rows, A.cols(), num_traits<T>::from_int(0));
    std::size_t w = 0;
    for (std::size_t i = 0; i < A.rows(); ++i) {
        if (i == skip) continue;
        for (std::size_t j = 0; j < A.cols(); ++j) out(w, j) = A(i, j);
        ++w;
    }
    return out;
}

/** Non-negative integer vectors k with sum(k) == n and k <= cap. */
inline void ncldmx_compositions(int rem, const std::vector<int>& cap, std::size_t idx,
                                std::vector<int>& cur, std::vector<std::vector<int> >& out) {
    if (idx + 1 == cap.size()) {
        if (rem <= cap[idx]) {
            cur[idx] = rem;
            out.push_back(cur);
        }
        return;
    }
    const int hi = rem < cap[idx] ? rem : cap[idx];
    for (int v = 0; v <= hi; ++v) {
        cur[idx] = v;
        ncldmx_compositions(rem - v, cap, idx + 1, cur, out);
    }
    cur[idx] = 0;
}

/**
 * First column of the trailing constant run of a limited load-dependence row,
 * i.e. the level b with mu(n) = mu(b) for every n >= b. This is the level
 * pfqn_ldmx_ec infers, and hence the one past which EC is constant. One-based.
 */
template <class T>
std::size_t ncldmx_lld_level(const Matrix<T>& mu, std::size_t row) {
    std::size_t b = mu.cols();
    if (b == 0) return 1;
    while (b > 1 && mu(row, b - 2) == mu(row, b - 1)) --b;
    return b;
}

}  // namespace detail

/**
 * @param lambda (R) arrival rates; must be zero on the closed classes
 * @param D      (M x R) service demands
 * @param N      (R) populations; a NEGATIVE entry marks an open class
 * @param Z      (K x R) think times of the closed classes
 * @param mu     (M x >=Kc) load-dependent rate lattice
 */
template <class T>
NcldmxResult<T> pfqn_ncldmx(const std::vector<T>& lambda, const Matrix<T>& D,
                            const std::vector<int>& N, const Matrix<T>& Z, const Matrix<T>& mu,
                            NcldMethod method, const T& atol, const NcOptions& nopt) {
    const std::size_t M = D.rows();
    const std::size_t R = D.cols();
    if (N.size() != R) throw InputError("pfqn_ncldmx: D and N disagree on the class count");
    if (lambda.size() != R)
        throw InputError("pfqn_ncldmx: lambda and D disagree on the class count");
    if (mu.rows() != M) throw InputError("pfqn_ncldmx: mu and D disagree on the station count");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    std::vector<std::size_t> openCl, closedCl;
    for (std::size_t r = 0; r < R; ++r) (N[r] < 0 ? openCl : closedCl).push_back(r);
    for (std::size_t r : closedCl)
        if (lambda[r] != zero)
            throw InputError("pfqn_ncldmx: an arrival rate is specified on a closed class");

    std::vector<int> Nc(closedCl.size(), 0);
    long Kc = 0;
    for (std::size_t a = 0; a < closedCl.size(); ++a) {
        Nc[a] = N[closedCl[a]];
        Kc += Nc[a];
    }
    const std::size_t width = static_cast<std::size_t>(Kc > 0 ? Kc : 1);

    // PAD the rate lattice out to the closed population, then add one extra
    // column, matching how pfqn_mvaldmx is fed. Padding NEVER SHORTENS: the row
    // the caller handed in is kept in full and only extended, because
    // pfqn_ldmx_ec reads the saturation level b_i off this row -- the first
    // column equal to the last -- and a row cut at the closed population
    // declares a c-server station saturated at min(n,c) with n < c whenever c
    // exceeds it. Cutting to width+1 read the M/M/3 of test_nc as an M/M/2
    // (mean 3.4286 against the exact 1.7368) and, with no closed class at all,
    // as an M/M/1. The reference pads with `repmat(mup(:,end),...)` and appends,
    // never truncating; Pfqn_ncldmx.java and ncldmx.py already take the wider
    // of the two.
    const std::size_t padCols = (mu.cols() > width ? mu.cols() : width) + 1;
    Matrix<T> mup(M, padCols);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < padCols; ++k)
            mup(i, k) = mu.cols() == 0 ? one : mu(i, k < mu.cols() ? k : mu.cols() - 1);

    std::vector<T> lambdao(R, zero);
    for (std::size_t r : openCl) lambdao[r] = lambda[r];

    const LdmxEcResult<T> ec = pfqn_ldmx_ec(lambdao, D, mup);

    NcldmxResult<T> res;
    res.EC = ec.EC;
    res.Gopen = one;
    for (std::size_t i = 0; i < M; ++i) res.Gopen *= ec.E(i, 0);
    res.lGopen = num_traits<T>::log_as_double(res.Gopen);

    Matrix<T> Dc(M, closedCl.size());
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t a = 0; a < closedCl.size(); ++a) Dc(i, a) = D(i, closedCl[a]);
    Matrix<T> Zc(Z.empty() ? 0 : Z.rows(), closedCl.size());
    for (std::size_t i = 0; i < Zc.rows(); ++i)
        for (std::size_t a = 0; a < closedCl.size(); ++a) Zc(i, a) = Z(i, closedCl[a]);

    Matrix<T> muEff(M, width);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < width; ++k) {
            if (ec.EC(i, k) == zero)
                throw NumericError("pfqn_ncldmx: an effective capacity term is zero");
            muEff(i, k) = one / ec.EC(i, k);
        }

    if (Kc == 0) {
        res.G = one;
        res.lG = 0.0;
    } else {
        const NcldResult<T> nc = pfqn_ncld(Dc, Nc, Zc, muEff, method, atol, nopt);
        res.G = nc.G;
        res.lG = nc.lG;
    }

    // ---- mean measures ----
    const std::size_t C = closedCl.size();
    res.XN.assign(R, zero);
    res.QN = Matrix<T>(M, R, zero);
    for (std::size_t r : openCl) res.XN[r] = lambda[r];

    std::vector<double> lGr(C > 0 ? C : 1, 0.0);
    if (Kc > 0) {
        for (std::size_t a = 0; a < C; ++a) {
            if (Nc[a] <= 0) continue;
            std::vector<int> Ncr = Nc;
            --Ncr[a];
            lGr[a] = detail::ncldmx_lg(Dc, Ncr, Zc, muEff, method, atol, nopt);
            res.XN[closedCl[a]] = num_traits<T>::from_double(std::exp(lGr[a] - res.lG));
        }
        // closed queue lengths: the conditional normalizing-constant recursion of
        // the load-dependent closed network, on the effective-capacity rates
        for (std::size_t i = 0; i < M; ++i) {
            bool anyDemand = false;
            for (std::size_t a = 0; a < C; ++a)
                if (Dc(i, a) > zero) { anyDemand = true; break; }
            if (!anyDemand) continue;
            const Matrix<T> muhat = pfqn_mushift(muEff, i);
            Matrix<T> muhatRow(1, muhat.cols(), zero);
            for (std::size_t k = 0; k < muhat.cols(); ++k) muhatRow(0, k) = muhat(i, k);
            const FncResult<T> fnc = pfqn_fnc(muhatRow);
            const double cshift = num_traits<T>::to_double(fnc.c[0]);
            const Matrix<T> Dminus = detail::ncldmx_drop_row(Dc, i);
            const Matrix<T> muminus = detail::ncldmx_drop_row(muEff, i);
            Matrix<T> DcPlus(M + 1, C, zero);
            for (std::size_t i2 = 0; i2 < M; ++i2)
                for (std::size_t a = 0; a < C; ++a) DcPlus(i2, a) = Dc(i2, a);
            for (std::size_t a = 0; a < C; ++a) DcPlus(M, a) = Dc(i, a);
            Matrix<T> muhatPlus(M + 1, muhat.cols(), zero);
            for (std::size_t i2 = 0; i2 < M; ++i2)
                for (std::size_t k = 0; k < muhat.cols(); ++k) muhatPlus(i2, k) = muhat(i2, k);
            for (std::size_t k = 0; k < muhat.cols() && k < fnc.mu.cols(); ++k)
                muhatPlus(M, k) = fnc.mu(0, k);
            for (std::size_t a = 0; a < C; ++a) {
                if (Nc[a] <= 0 || !(Dc(i, a) > zero)) continue;
                std::vector<int> Ncr = Nc;
                --Ncr[a];
                const double lGhat = detail::ncldmx_lg(Dc, Ncr, Zc, muhat, method, atol, nopt);
                const double lGhatf =
                    detail::ncldmx_lg(DcPlus, Ncr, Zc, muhatPlus, method, atol, nopt);
                const double lGminus =
                    detail::ncldmx_lg(Dminus, Ncr, Zc, muminus, method, atol, nopt);
                const double CQ =
                    (std::exp(lGhatf - lGhat) - 1.0) + cshift * (std::exp(lGminus - lGhat) - 1.0);
                const double ldDemand = std::log(num_traits<T>::to_double(Dc(i, a))) + lGhat -
                                        std::log(num_traits<T>::to_double(muEff(i, 0))) - lGr[a];
                res.QN(i, closedCl[a]) = num_traits<T>::from_double(
                    std::exp(ldDemand) * num_traits<T>::to_double(res.XN[closedCl[a]]) * (1.0 + CQ));
            }
        }
    }

    // open queue lengths, with the saturated tail of EC folded onto the closed mean
    if (!openCl.empty()) {
        for (std::size_t i = 0; i < M; ++i) {
            double Qtot = 0.0;
            for (std::size_t a = 0; a < C; ++a)
                Qtot += num_traits<T>::to_double(res.QN(i, closedCl[a]));
            const std::size_t b = detail::ncldmx_lld_level(mup, i);
            const std::size_t bcap = b < ec.EC.cols() ? b : ec.EC.cols();
            const double ECinf = num_traits<T>::to_double(ec.EC(i, bcap - 1));
            double acc = ECinf * (Qtot + 1.0);
            if (b >= 2) {
                const Matrix<T> Dminus = detail::ncldmx_drop_row(Dc, i);
                const Matrix<T> muminus = detail::ncldmx_drop_row(muEff, i);
                Matrix<T> Drow(1, C, zero), murow(1, muEff.cols(), zero);
                for (std::size_t a = 0; a < C; ++a) Drow(0, a) = Dc(i, a);
                for (std::size_t k = 0; k < muEff.cols(); ++k) murow(0, k) = muEff(i, k);
                const Matrix<T> Zzero(1, C, zero);
                for (std::size_t n = 0; n + 2 <= b; ++n) {
                    const double delta = num_traits<T>::to_double(ec.EC(i, n)) - ECinf;
                    if (delta == 0.0) continue;
                    double Pn = 0.0;
                    if (Kc == 0) {
                        // WITH NO CLOSED POPULATION THE MARGINAL IS DEGENERATE, not
                        // absent: P_i(0) = 1 and P_i(n) = 0 above it, so only the n = 0
                        // term survives and acc collapses to EC_i(1), the exact open
                        // load-dependent mean. Skipping the loop instead left acc at
                        // EC_i^inf, i.e. read a c-server station as if every arrival found
                        // it saturated -- an M/M/3 at lambda = 1.5 came back with mean 1
                        // against the exact 1.7368.
                        Pn = n == 0 ? 1.0 : 0.0;
                    } else {
                        std::vector<std::vector<int> > ks;
                        std::vector<int> cur(C, 0);
                        if (C > 0)
                            detail::ncldmx_compositions(static_cast<int>(n), Nc, 0, cur, ks);
                        else if (n == 0)
                            ks.push_back(std::vector<int>());
                        for (std::size_t t = 0; t < ks.size(); ++t) {
                            std::vector<int> rest(C, 0);
                            for (std::size_t a = 0; a < C; ++a) rest[a] = Nc[a] - ks[t][a];
                            const double lF = n == 0 ? 0.0
                                                     : detail::ncldmx_lg(Drow, ks[t], Zzero, murow,
                                                                         method, atol, nopt);
                            const double lGbar =
                                detail::ncldmx_lg(Dminus, rest, Zc, muminus, method, atol, nopt);
                            Pn += std::exp(lF + lGbar - res.lG);
                        }
                    }
                    acc += static_cast<double>(n + 1) * delta * Pn;
                }
            }
            for (std::size_t r : openCl)
                res.QN(i, r) = num_traits<T>::from_double(
                    num_traits<T>::to_double(lambda[r]) * num_traits<T>::to_double(D(i, r)) * acc);
        }
    }

    return res;
}

/** Overload with the exact (zero-tolerance) filters and default sampling options. */
template <class T>
NcldmxResult<T> pfqn_ncldmx(const std::vector<T>& lambda, const Matrix<T>& D,
                            const std::vector<int>& N, const Matrix<T>& Z, const Matrix<T>& mu) {
    return pfqn_ncldmx(lambda, D, N, Z, mu, NcldMethod::Default, num_traits<T>::from_int(0),
                       NcOptions());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_NCLDMX_H
