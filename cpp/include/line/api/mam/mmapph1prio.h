/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMAPPH1PRIO_H
#define LINE_API_MAM_MMAPPH1PRIO_H

/**
 * @file
 * @ingroup api_mam
 * The MMAP[K]/PH[K]/1 priority queue, preemptive resume (`mmapph1prpr_*`) and
 * non-preemptive (`mmapph1nppr_*`).
 *
 * Port of BUTools' MMAPPH1PRPR.m and MMAPPH1NPPR.m (matlab/lib/thirdparty/
 * BUTools/queues), which `solver_mam_basic.m` calls at an FCFSPRPRIO or HOL
 * station with all-distinct class priorities, and `solver_mam_passage_time.m`
 * tabulates the sojourn CDF from. The algorithm is G. Horvath, "Efficient
 * analysis of the MMAP[K]/PH[K]/1 priority queue", European Journal of
 * Operational Research 246(1), 128-139, 2015: the workload of the classes at or
 * above k is a fluid queue whose first-return matrix Psi (ADDA doubling,
 * `mfq_fundamental`) yields the boundary vector, and the remaining sojourn time
 * of a class-k job is the first passage time of a second fluid model over the
 * higher classes. Moments come from a chain of Sylvester equations sharing
 * their coefficient matrices; the CDF from an Erlangized Laplace inversion of
 * order `erl_max_order`.
 *
 * CLASS ORDER IS BUTools' OWN: class 1 is the LOWEST priority and class K the
 * highest. The caller permutes LINE's classes (lower classprio value = higher
 * priority) into that order and maps the outputs back, as the reference does.
 *
 * Every entry point returns one row per class, in the order of the input.
 *
 * ARITHMETIC. The Riccati and QBD roots are tolerance-terminated iterations, so
 * these are {Double, Real} and refuse exact/Rational at compile time, like
 * `mmapph1fcfs`.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/mfq_solve.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmapph1fcfs.h"
#include "line/api/mam/qbd_r.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"
#include "line/util/sylvester.h"

namespace line {
namespace mam {

/** Options shared by both priority analyzers (the BUTools 'erlMaxOrder' and 'prec'). */
struct PrioQueueOptions {
    std::size_t erl_max_order = 200;  ///< Erlang order of the sojourn-time CDF inversion
    double precision = 1e-14;         ///< tolerance of the fluid and QBD fundamental matrices
};

namespace prio_detail {

enum class Measure { StMoms, StDistr, NcMoms, NcDistr };

template <class T>
Matrix<T> zeros(std::size_t r, std::size_t c) {
    return Matrix<T>(r, c, num_traits<T>::from_int(0));
}

template <class T>
Matrix<T> eye(std::size_t n) {
    Matrix<T> I = zeros<T>(n, n);
    for (std::size_t i = 0; i < n; ++i) I(i, i) = num_traits<T>::from_int(1);
    return I;
}

template <class T>
Matrix<T> ones(std::size_t r, std::size_t c) {
    return Matrix<T>(r, c, num_traits<T>::from_int(1));
}

template <class T>
Matrix<T> add(const Matrix<T>& A, const Matrix<T>& B) {
    Matrix<T> C = A;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) += B(i, j);
    return C;
}

template <class T>
Matrix<T> sub(const Matrix<T>& A, const Matrix<T>& B) {
    Matrix<T> C = A;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) -= B(i, j);
    return C;
}

template <class T>
Matrix<T> scal(const Matrix<T>& A, const T& s) {
    Matrix<T> C = A;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) = T(C(i, j) * s);
    return C;
}

template <class T>
Matrix<T> mul(const Matrix<T>& A, const Matrix<T>& B) {
    if (A.rows() == 0 || B.cols() == 0 || A.cols() == 0) return zeros<T>(A.rows(), B.cols());
    return matmul(A, B);
}

template <class T>
Matrix<T> mul(const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C) {
    return mul(mul(A, B), C);
}

/** inv(-A). */
template <class T>
Matrix<T> ninv(const Matrix<T>& A) {
    return inverse(scal(A, T(num_traits<T>::from_int(-1))));
}

template <class T>
Matrix<T> mpow(const Matrix<T>& A, std::size_t n) {
    Matrix<T> P = eye<T>(A.rows());
    for (std::size_t i = 0; i < n; ++i) P = mul(P, A);
    return P;
}

template <class T>
void setblk(Matrix<T>& A, std::size_t r0, std::size_t c0, const Matrix<T>& B) {
    for (std::size_t i = 0; i < B.rows(); ++i)
        for (std::size_t j = 0; j < B.cols(); ++j) A(r0 + i, c0 + j) = B(i, j);
}

template <class T>
Matrix<T> getblk(const Matrix<T>& A, std::size_t r0, std::size_t c0, std::size_t nr,
                 std::size_t nc) {
    Matrix<T> B = zeros<T>(nr, nc);
    for (std::size_t i = 0; i < nr; ++i)
        for (std::size_t j = 0; j < nc; ++j) B(i, j) = A(r0 + i, c0 + j);
    return B;
}

template <class T>
Matrix<T> hcat(const Matrix<T>& A, const Matrix<T>& B) {
    const std::size_t r = A.cols() > 0 ? A.rows() : B.rows();
    Matrix<T> C = zeros<T>(r, A.cols() + B.cols());
    setblk(C, 0, 0, A);
    setblk(C, 0, A.cols(), B);
    return C;
}

template <class T>
Matrix<T> vcat(const Matrix<T>& A, const Matrix<T>& B) {
    const std::size_t c = A.rows() > 0 ? A.cols() : B.cols();
    Matrix<T> C = zeros<T>(A.rows() + B.rows(), c);
    setblk(C, 0, 0, A);
    setblk(C, A.rows(), 0, B);
    return C;
}

template <class T>
Matrix<T> blkdiag(const Matrix<T>& A, const Matrix<T>& B) {
    Matrix<T> C = zeros<T>(A.rows() + B.rows(), A.cols() + B.cols());
    setblk(C, 0, 0, A);
    setblk(C, A.rows(), A.cols(), B);
    return C;
}

/** sum(A,2), as a column. */
template <class T>
Matrix<T> rowsum(const Matrix<T>& A) {
    Matrix<T> s = zeros<T>(A.rows(), 1);
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) s(i, 0) += A(i, j);
    return s;
}

/** sum(A(:)). */
template <class T>
T total(const Matrix<T>& A) {
    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) s += A(i, j);
    return s;
}

template <class T>
Matrix<T> rowvec(const std::vector<T>& v) {
    Matrix<T> r = zeros<T>(1, v.size());
    for (std::size_t j = 0; j < v.size(); ++j) r(0, j) = v[j];
    return r;
}

template <class T>
T factorial(std::size_t n) {
    T f = num_traits<T>::from_int(1);
    for (std::size_t i = 2; i <= n; ++i) f = T(f * num_traits<T>::from_int((int)i));
    return f;
}

/** BUTools MomsFromFactorialMoms: raw moments from factorial moments (signed Stirling numbers of the first kind). */
template <class T>
std::vector<T> moms_from_factorial_moms(const std::vector<T>& fm) {
    std::vector<T> m(fm.size());
    if (fm.empty()) return m;
    m[0] = fm[0];
    for (std::size_t i = 2; i <= fm.size(); ++i) {
        // coefficients of x(x-1)...(x-i+1), c[j] multiplying x^j
        std::vector<double> c(i + 1, 0.0);
        c[0] = 1.0;
        for (std::size_t r = 0; r < i; ++r) {
            std::vector<double> nc(i + 1, 0.0);
            for (std::size_t j = 0; j <= r; ++j) {
                nc[j + 1] += c[j];
                nc[j] -= static_cast<double>(r) * c[j];
            }
            c.swap(nc);
        }
        T v = fm[i - 1];
        for (std::size_t j = 1; j < i; ++j) v -= T(num_traits<T>::from_double(c[j]) * m[j - 1]);
        m[i - 1] = v;
    }
    return m;
}

/** The inputs of both analyzers in BUTools' layout: D{1} = D0, D{k+1} = class k. */
template <class T>
struct PrioInput {
    std::size_t N = 0, K = 0;
    Matrix<T> D0, sD, I;
    std::vector<Matrix<T>> D;      ///< D[k] is BUTools' D{k+1}
    std::vector<Matrix<T>> sigma;  ///< 1 x M[k]
    std::vector<Matrix<T>> S;      ///< M[k] x M[k]
    std::vector<Matrix<T>> s;      ///< M[k] x 1 exit vector
    std::vector<std::size_t> M;
    T prec;
    std::size_t erl = 200;

    std::size_t msum(std::size_t a, std::size_t b) const {  // sum M[a..b-1]
        std::size_t t = 0;
        for (std::size_t i = a; i < b; ++i) t += M[i];
        return t;
    }
};

template <class T>
PrioInput<T> prepare(const Mmap<T>& arrival, const std::vector<PhService<T>>& svc,
                     const PrioQueueOptions& opt, const char* who) {
    static_assert(num_traits<T>::has_transcendental,
                  "the MMAP/PH/1 priority analyzers run tolerance-terminated Riccati iterations");
    PrioInput<T> in;
    in.K = svc.size();
    if (in.K == 0) throw InputError(std::string(who) + ": at least one class is required");
    if (arrival.classes() != in.K)
        throw InputError(std::string(who) +
                         ": the arrival MMAP and the service list disagree on the number of classes");
    if (opt.erl_max_order < 1) throw InputError(std::string(who) + ": erl_max_order must be >= 1");
    in.N = arrival.order();
    in.D0 = arrival.D0;
    in.D = arrival.Dc;
    in.I = eye<T>(in.N);
    in.sD = in.D0;
    for (std::size_t k = 0; k < in.K; ++k) in.sD = add(in.sD, in.D[k]);
    for (std::size_t k = 0; k < in.K; ++k) {
        const std::size_t m = svc[k].S.rows();
        if (m == 0 || svc[k].S.cols() != m || svc[k].sigma.size() != m)
            throw InputError(std::string(who) + ": class " + std::to_string(k + 1) +
                             " has a malformed PH service representation");
        in.M.push_back(m);
        in.sigma.push_back(rowvec(svc[k].sigma));
        in.S.push_back(svc[k].S);
        in.s.push_back(scal(rowsum(svc[k].S), T(num_traits<T>::from_int(-1))));
    }
    in.prec = num_traits<T>::from_double(opt.precision);
    in.erl = opt.erl_max_order;
    return in;
}

template <class T>
FluidFundamental<T> fluid(const PrioInput<T>& in, const Matrix<T>& Fpp, const Matrix<T>& Fpm,
                          const Matrix<T>& Fmp, const Matrix<T>& Fmm) {
    return mfq_fundamental(Fpp, Fpm, Fmp, Fmm, in.prec, 150u, RiccatiMethod::ADDA);
}

/** The fluid model of the remaining sojourn (PRPR) or waiting (NPPR) time, started from `inis`. */
template <class T>
struct SojournModel {
    Matrix<T> Qspp, Qspm, Qsmp, Qsmm, inis, Psis;
};

/** P_n of the sojourn-moment recursion: C_n = -2n P_{n-1} + sum bino P_i Qsmp P_{n-i}. */
template <class T>
std::vector<Matrix<T>> st_moment_chain(const SojournModel<T>& sm, std::size_t n) {
    const Matrix<T> A = add(sm.Qspp, mul(sm.Psis, sm.Qsmp));
    const Matrix<T> B = add(sm.Qsmm, mul(sm.Qsmp, sm.Psis));
    const SylvesterFactor<T> F(A, B);
    std::vector<Matrix<T>> Pn{sm.Psis};
    std::vector<Matrix<T>> QP{mul(sm.Qsmp, sm.Psis)};
    for (std::size_t m = 1; m <= n; ++m) {
        Matrix<T> C = scal(Pn[m - 1], T(num_traits<T>::from_int(-2 * (int)m)));
        T bino = num_traits<T>::from_int(1);
        for (std::size_t i = 1; i + 1 <= m; ++i) {
            bino = T(bino * num_traits<T>::from_int((int)(m - i + 1)) /
                     num_traits<T>::from_int((int)i));
            C = add(C, scal(mul(Pn[i], QP[m - i]), bino));
        }
        Pn.push_back(F.solve_lyap(C));
        QP.push_back(mul(sm.Qsmp, Pn.back()));
    }
    return Pn;
}

/**
 * The Erlangized first-passage terms at one CDF point: Psie and P_1..P_{L-1}
 * of the lambda-shifted model, C_n = 2 lambda P_{n-1} + sum P_i Qsmp P_{n-i}.
 * Returns sum(inis*P_n) for n = 0..L-1.
 */
template <class T>
std::vector<T> erlang_terms(const PrioInput<T>& in, const SojournModel<T>& sm, const T& lambda) {
    const std::size_t L = in.erl;
    const Matrix<T> Ip = eye<T>(sm.Qspp.rows()), Im = eye<T>(sm.Qsmm.rows());
    const Matrix<T> Qpp = sub(sm.Qspp, scal(Ip, lambda));
    const Matrix<T> Qmm = sub(sm.Qsmm, scal(Im, lambda));
    const Matrix<T> Psie = fluid(in, Qpp, sm.Qspm, sm.Qsmp, Qmm).Psi;
    std::vector<T> out;
    out.push_back(total(mul(sm.inis, Psie)));
    const Matrix<T> A = add(Qpp, mul(Psie, sm.Qsmp));
    const Matrix<T> B = add(Qmm, mul(sm.Qsmp, Psie));
    const SylvesterFactor<T> F(A, B);
    std::vector<Matrix<T>> Pn{Psie};
    std::vector<Matrix<T>> QP{mul(sm.Qsmp, Psie)};
    const T two_l = T(num_traits<T>::from_int(2) * lambda);
    for (std::size_t n = 1; n + 1 <= L; ++n) {
        Matrix<T> C = scal(Pn[n - 1], two_l);
        for (std::size_t i = 1; i + 1 <= n; ++i) C = add(C, mul(Pn[i], QP[n - i]));
        Pn.push_back(F.solve_lyap(C));
        QP.push_back(mul(sm.Qsmp, Pn.back()));
        out.push_back(total(mul(sm.inis, Pn.back())));
    }
    return out;
}

/** Departure-instant factorial moments to random-time moments, shared by both analyzers. */
template <class T>
std::vector<T> random_time_moms(const PrioInput<T>& in, std::size_t k,
                                const std::vector<Matrix<T>>& QLDPn, std::size_t n) {
    const std::vector<T> piv = mc::ctmc_solve(in.sD);
    const Matrix<T> pi = rowvec(piv);
    const T lambdak = total(mul(pi, in.D[k]));
    const Matrix<T> iTerm = inverse(sub(mul(ones<T>(in.N, 1), pi), in.sD));
    const Matrix<T> dk1 = rowsum(in.D[k]);
    std::vector<Matrix<T>> QLPn{pi};
    std::vector<T> ql(n);
    for (std::size_t m = 1; m <= n; ++m) {
        const T mm = num_traits<T>::from_int((int)m);
        const Matrix<T> a = sub(QLDPn[m - 1], scal(mul(QLPn[m - 1], in.D[k]), T(num_traits<T>::from_int(1) / lambdak)));
        const T sumP = T(total(QLDPn[m]) + mm * total(mul(a, iTerm, dk1)));
        const Matrix<T> b = sub(mul(QLPn[m - 1], in.D[k]), scal(QLDPn[m - 1], lambdak));
        const Matrix<T> P = add(scal(pi, sumP), scal(mul(b, iTerm), mm));
        QLPn.push_back(P);
        ql[m - 1] = total(P);
    }
    return moms_from_factorial_moms(ql);
}

/** Departure-instant probabilities to random-time probabilities, shared by both analyzers. */
template <class T>
std::vector<T> random_time_probs(const PrioInput<T>& in, std::size_t k,
                                 const std::vector<Matrix<T>>& dql) {
    const Matrix<T> pi = rowvec(mc::ctmc_solve(in.sD));
    const T lambdak = total(mul(pi, in.D[k]));
    const Matrix<T> iTerm = ninv(sub(in.sD, in.D[k]));
    std::vector<T> out;
    Matrix<T> q = mul(scal(dql[0], lambdak), iTerm);
    out.push_back(total(q));
    for (std::size_t n = 1; n < dql.size(); ++n) {
        q = mul(add(mul(q, in.D[k]), scal(sub(dql[n], dql[n - 1]), lambdak)), iTerm);
        out.push_back(total(q));
    }
    return out;
}

/** The QBD tail both analyzers use for the number of jobs of the TOP class. */
template <class T>
std::vector<T> qbd_measure(const Matrix<T>& R, Matrix<T> p0, Measure what, std::size_t n) {
    const Matrix<T> IR = sub(eye<T>(R.rows()), R);
    const Matrix<T> iIR = inverse(IR);
    p0 = scal(p0, T(num_traits<T>::from_int(1) / total(mul(p0, iIR))));
    std::vector<T> out;
    if (what == Measure::NcMoms) {
        for (std::size_t i = 1; i <= n; ++i)
            out.push_back(total(scal(mul(p0, mpow(R, i), mpow(iIR, i + 1)), factorial<T>(i))));
        return moms_from_factorial_moms(out);
    }
    Matrix<T> v = p0;
    out.push_back(total(v));
    for (std::size_t i = 1; i < n; ++i) {
        v = mul(v, R);
        out.push_back(total(v));
    }
    return out;
}

template <class T>
Matrix<T> qbd_R_of(const PrioInput<T>& in, const Matrix<T>& B, const Matrix<T>& L,
                   const Matrix<T>& F) {
    return qbd_R_logred(B, L, F, 10000u, in.prec);
}

/** Sojourn moments / CDF of a class whose law is the PH (zeta, Z) with exit z. */
template <class T>
std::vector<T> ph_sojourn(const Matrix<T>& zeta, const Matrix<T>& Z, const Matrix<T>& z,
                          Measure what, std::size_t n, const std::vector<T>& pts) {
    const Matrix<T> iZ = ninv(Z);
    std::vector<T> out;
    if (what == Measure::StMoms) {
        for (std::size_t i = 1; i <= n; ++i)
            out.push_back(T(factorial<T>(i) * total(mul(zeta, mpow(iZ, i + 1), z))));
    } else {
        const Matrix<T> I = eye<T>(Z.rows());
        for (const T& t : pts)
            out.push_back(total(mul(mul(zeta, iZ), sub(I, expm(scal(Z, t))), z)));
    }
    return out;
}

// ---------------------------------------------------------------- PRPR

template <class T>
std::vector<T> prpr_class(const PrioInput<T>& in, std::size_t k, Measure what, std::size_t n,
                          const std::vector<T>& pts) {
    const std::size_t N = in.N, K = in.K;
    const Matrix<T>& I = in.I;
    // step 1. workload process of the classes k..K
    const std::size_t sM = in.msum(k, K);
    Matrix<T> Qwmm = in.D0;
    for (std::size_t i = 0; i < k; ++i) Qwmm = add(Qwmm, in.D[i]);
    Matrix<T> Qwpm = zeros<T>(N * sM, N), Qwmp = zeros<T>(N, N * sM), Qwpp = zeros<T>(N * sM, N * sM);
    std::size_t kix = 0;
    for (std::size_t i = k; i < K; ++i) {
        setblk(Qwmp, 0, kix, kron(in.D[i], in.sigma[i]));
        setblk(Qwpm, kix, 0, kron(I, in.s[i]));
        setblk(Qwpp, kix, kix, kron(I, in.S[i]));
        kix += N * in.M[i];
    }
    const FluidFundamental<T> fw = fluid(in, Qwpp, Qwpm, Qwmp, Qwmm);
    const Matrix<T>& Kw = fw.K;
    const Matrix<T> iKw = ninv(Kw);
    const Matrix<T> Ua =
        add(ones<T>(N, 1), scal(rowsum(mul(Qwmp, iKw)), T(num_traits<T>::from_int(2))));
    Matrix<T> UwUa = hcat(fw.U, Ua);  // N x (N+1); pm [Uw, Ua] = [0 .. 0, 1]
    Matrix<T> At = zeros<T>(N + 1, N);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j <= N; ++j) At(j, i) = UwUa(i, j);
    std::vector<T> rhs(N + 1, num_traits<T>::from_int(0));
    rhs[N] = num_traits<T>::from_int(1);
    const Matrix<T> pm = rowvec(mfq_detail::normal_equations_solve(At, rhs));

    Matrix<T> Bw = zeros<T>(N * sM, N);
    setblk(Bw, 0, 0, kron(I, in.s[k]));
    const Matrix<T> pmQ = mul(pm, Qwmp);
    const Matrix<T> kappa = scal(pmQ, T(num_traits<T>::from_int(1) / total(mul(pmQ, iKw, Bw))));

    if (k + 1 == K) {
        if (what == Measure::StMoms || what == Measure::StDistr)
            return ph_sojourn(kappa, Kw, rowsum(Bw), what, n, pts);
        const std::size_t Mk = in.M[k];
        const Matrix<T> IM = eye<T>(Mk);
        const Matrix<T> sDk = sub(in.sD, in.D[k]);
        const Matrix<T> L = add(kron(sDk, IM), kron(I, in.S[k]));
        const Matrix<T> B = kron(I, mul(in.s[k], in.sigma[k]));
        const Matrix<T> F = kron(in.D[k], IM);
        const Matrix<T> L0 = kron(sDk, IM);
        const Matrix<T> R = qbd_R_of(in, B, L, F);
        const Matrix<T> p0 = rowvec(mc::ctmc_solve(add(L0, mul(R, B))));
        return qbd_measure(R, p0, what, n);
    }

    // step 2. fluid model of the remaining sojourn time
    SojournModel<T> sm;
    sm.Qsmm = in.D0;
    for (std::size_t i = 0; i <= k; ++i) sm.Qsmm = add(sm.Qsmm, in.D[i]);
    const std::size_t Np = Kw.rows(), ext = N * in.msum(k + 1, K);
    sm.Qspm = zeros<T>(Np + ext, N);
    sm.Qsmp = zeros<T>(N, Np + ext);
    sm.Qspp = zeros<T>(Np + ext, Np + ext);
    setblk(sm.Qspp, 0, 0, Kw);
    setblk(sm.Qspm, 0, 0, Bw);
    kix = Np;
    for (std::size_t i = k + 1; i < K; ++i) {
        setblk(sm.Qsmp, 0, kix, kron(in.D[i], in.sigma[i]));
        setblk(sm.Qspm, kix, 0, kron(I, in.s[i]));
        setblk(sm.Qspp, kix, kix, kron(I, in.S[i]));
        kix += N * in.M[i];
    }
    sm.inis = hcat(kappa, zeros<T>(1, ext));
    sm.Psis = fluid(in, sm.Qspp, sm.Qspm, sm.Qsmp, sm.Qsmm).Psi;

    std::vector<T> out;
    if (what == Measure::StMoms) {
        const std::vector<Matrix<T>> Pn = st_moment_chain(sm, n);
        for (std::size_t m = 1; m <= n; ++m) {
            T v = T(total(mul(sm.inis, Pn[m])) / num_traits<T>::from_double(std::pow(2.0, (double)m)));
            if (m % 2 == 1) v = T(-v);
            out.push_back(v);
        }
    } else if (what == Measure::StDistr) {
        for (const T& t : pts) {
            const T lambda = T(num_traits<T>::from_int((int)in.erl) / t / num_traits<T>::from_int(2));
            T pr = num_traits<T>::from_int(0);
            for (const T& v : erlang_terms(in, sm, lambda)) pr += v;
            out.push_back(pr);
        }
    } else if (what == Measure::NcMoms) {
        const Matrix<T> A = add(sm.Qspp, mul(sm.Psis, sm.Qsmp));
        const Matrix<T> B = add(sm.Qsmm, mul(sm.Qsmp, sm.Psis));
        const SylvesterFactor<T> F(A, B);
        std::vector<Matrix<T>> P{sm.Psis};
        std::vector<Matrix<T>> QP{mul(sm.Qsmp, sm.Psis)};
        for (std::size_t m = 1; m <= n; ++m) {
            Matrix<T> C = scal(mul(P[m - 1], in.D[k]), T(num_traits<T>::from_int((int)m)));
            T bino = num_traits<T>::from_int(1);
            for (std::size_t i = 1; i + 1 <= m; ++i) {
                bino = T(bino * num_traits<T>::from_int((int)(m - i + 1)) /
                         num_traits<T>::from_int((int)i));
                C = add(C, scal(mul(P[i], QP[m - i]), bino));
            }
            P.push_back(F.solve_lyap(C));
            QP.push_back(mul(sm.Qsmp, P.back()));
        }
        std::vector<Matrix<T>> QLDPn;
        for (const Matrix<T>& p : P) QLDPn.push_back(mul(sm.inis, p));
        out = random_time_moms(in, k, QLDPn, n);
    } else {
        Matrix<T> sDk = in.D0;
        for (std::size_t i = 0; i < k; ++i) sDk = add(sDk, in.D[i]);
        const Matrix<T> Psid = fluid(in, sm.Qspp, sm.Qspm, sm.Qsmp, sDk).Psi;
        const Matrix<T> A = add(sm.Qspp, mul(Psid, sm.Qsmp));
        const Matrix<T> B = add(sDk, mul(sm.Qsmp, Psid));
        const SylvesterFactor<T> F(A, B);
        std::vector<Matrix<T>> P{Psid};
        std::vector<Matrix<T>> QP{mul(sm.Qsmp, Psid)};
        std::vector<Matrix<T>> dql{mul(sm.inis, Psid)};
        for (std::size_t m = 1; m < n; ++m) {
            Matrix<T> C = mul(P[m - 1], in.D[k]);
            for (std::size_t i = 1; i + 1 <= m; ++i) C = add(C, mul(P[i], QP[m - i]));
            P.push_back(F.solve_lyap(C));
            QP.push_back(mul(sm.Qsmp, P.back()));
            dql.push_back(mul(sm.inis, P.back()));
        }
        out = random_time_probs(in, k, dql);
    }
    return out;
}

// ---------------------------------------------------------------- NPPR

template <class T>
struct NpprModel {
    Matrix<T> pm;
    std::vector<Matrix<T>> Psiw, Qwmp, Qwzp, Qwpp, Qwmz, Qwpz, Qwzz, Qwmm, Qwpm, Qwzm;
    std::vector<Matrix<T>> q0, qL;  ///< q0[k], qL[k] for k = 1..K-1 (BUTools q0{k+1})
    std::vector<T> lambda;
};

template <class T>
NpprModel<T> nppr_build(const PrioInput<T>& in) {
    const std::size_t N = in.N, K = in.K;
    const Matrix<T>& I = in.I;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    NpprModel<T> md;

    // step 1. workload process of the joint queue
    const std::size_t sM = in.msum(0, K);
    Matrix<T> Qwpp = zeros<T>(N * sM, N * sM), Qwmp = zeros<T>(N, N * sM),
              Qwpm = zeros<T>(N * sM, N);
    std::size_t kix = 0;
    for (std::size_t i = 0; i < K; ++i) {
        setblk(Qwpp, kix, kix, kron(I, in.S[i]));
        setblk(Qwmp, 0, kix, kron(in.D[i], in.sigma[i]));
        setblk(Qwpm, kix, 0, kron(I, in.s[i]));
        kix += N * in.M[i];
    }
    const FluidFundamental<T> fw = fluid(in, Qwpp, Qwpm, Qwmp, in.D0);
    const Matrix<T> iKw = ninv(fw.K);
    const Matrix<T> Ua = add(ones<T>(N, 1), scal(rowsum(mul(Qwmp, iKw)), T(num_traits<T>::from_int(2))));
    const Matrix<T> UwUa = hcat(fw.U, Ua);
    Matrix<T> At = zeros<T>(N + 1, N);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j <= N; ++j) At(j, i) = UwUa(i, j);
    std::vector<T> rhs(N + 1, zero);
    rhs[N] = one;
    md.pm = rowvec(mfq_detail::normal_equations_solve(At, rhs));
    const T spm = total(md.pm);
    const T half_idle = T((one - spm) / num_traits<T>::from_int(2));
    const T ro = T(half_idle / (spm + half_idle));
    const Matrix<T> kappa = scal(md.pm, T(one / spm));

    const Matrix<T> pi = rowvec(mc::ctmc_solve(in.sD));
    md.lambda.resize(K);
    for (std::size_t i = 0; i < K; ++i) md.lambda[i] = total(mul(pi, in.D[i]));

    // step 2. workload process of the classes k..K, one per k
    for (std::size_t k = 0; k < K; ++k) {
        const std::size_t Mlo = in.msum(0, k), Mhi = in.msum(k, K);
        const std::size_t p = N * Mlo * Mhi + N * Mhi, z = N * Mlo;
        Matrix<T> Qkwpp = zeros<T>(p, p), Qkwpz = zeros<T>(p, z), Qkwpm = zeros<T>(p, N);
        Matrix<T> Qkwmz = zeros<T>(N, z), Qkwmp = zeros<T>(N, p);
        Matrix<T> Dlo = in.D0;
        for (std::size_t i = 0; i < k; ++i) Dlo = add(Dlo, in.D[i]);
        Matrix<T> Qkwzp = zeros<T>(z, p), Qkwzm = zeros<T>(z, N), Qkwzz = zeros<T>(z, z);
        kix = 0;
        for (std::size_t i = k; i < K; ++i) {
            std::size_t kix2 = 0;
            for (std::size_t j = 0; j < k; ++j) {
                const Matrix<T> IMj = eye<T>(in.M[j]);
                setblk(Qkwpp, kix, kix, kron(I, kron(IMj, in.S[i])));
                setblk(Qkwpz, kix, kix2, kron(I, kron(IMj, in.s[i])));
                setblk(Qkwzp, kix2, kix, kron(in.D[i], kron(IMj, in.sigma[i])));
                kix += N * in.M[j] * in.M[i];
                kix2 += N * in.M[j];
            }
        }
        for (std::size_t i = k; i < K; ++i) {
            setblk(Qkwpp, kix, kix, kron(I, in.S[i]));
            setblk(Qkwpm, kix, 0, kron(I, in.s[i]));
            setblk(Qkwmp, 0, kix, kron(in.D[i], in.sigma[i]));
            kix += N * in.M[i];
        }
        kix = 0;
        for (std::size_t j = 0; j < k; ++j) {
            setblk(Qkwzz, kix, kix, add(kron(Dlo, eye<T>(in.M[j])), kron(I, in.S[j])));
            setblk(Qkwzm, kix, 0, kron(I, in.s[j]));
            kix += N * in.M[j];
        }
        Matrix<T> Fpp = Qkwpp, Fpm = Qkwpm;
        if (z > 0) {
            const Matrix<T> iZ = ninv(Qkwzz);
            Fpp = add(Fpp, mul(Qkwpz, iZ, Qkwzp));
            Fpm = add(Fpm, mul(Qkwpz, iZ, Qkwzm));
        }
        md.Psiw.push_back(fluid(in, Fpp, Fpm, Qkwmp, Dlo).Psi);
        md.Qwzp.push_back(Qkwzp);
        md.Qwmp.push_back(Qkwmp);
        md.Qwpp.push_back(Qkwpp);
        md.Qwmz.push_back(Qkwmz);
        md.Qwpz.push_back(Qkwpz);
        md.Qwzz.push_back(Qkwzz);
        md.Qwmm.push_back(Dlo);
        md.Qwpm.push_back(Qkwpm);
        md.Qwzm.push_back(Qkwzm);
    }

    // step 3. the phi vectors; phi[0] is BUTools' phi{1}
    T lambdaS = zero;
    for (const T& l : md.lambda) lambdaS += l;
    const Matrix<T> iD0 = ninv(in.D0);
    std::vector<Matrix<T>> phi{scal(mul(kappa, scal(in.D0, T(-one))), T((one - ro) / lambdaS))};
    md.q0.assign(K, Matrix<T>());
    md.qL.assign(K, Matrix<T>());
    for (std::size_t k = 1; k < K; ++k) {  // BUTools k = 1..K-1, counting the first k classes
        Matrix<T> sDk = in.D0;
        for (std::size_t j = 0; j < k; ++j) sDk = add(sDk, in.D[j]);
        T lk = zero;
        for (std::size_t j = 0; j < k; ++j) lk += md.lambda[j];
        const T pk = T(lk / lambdaS - (one - ro) * total(mul(kappa, rowsum(sDk))) / lambdaS);
        const Matrix<T>& Qwzpk = md.Qwzp[k];
        std::vector<Matrix<T>> Ak(k), Gi(k);
        std::size_t vix = 0;
        for (std::size_t ii = 0; ii < k; ++ii) {
            const std::size_t bs = N * in.M[ii];
            const Matrix<T> V1 = getblk(Qwzpk, vix, 0, bs, Qwzpk.cols());
            Gi[ii] = ninv(add(kron(sDk, eye<T>(in.M[ii])), kron(I, in.S[ii])));
            Ak[ii] = mul(mul(kron(I, in.sigma[ii]), Gi[ii]),
                         add(kron(I, in.s[ii]), mul(V1, md.Psiw[k])));
            vix += bs;
        }
        const Matrix<T> Bk = mul(md.Qwmp[k], md.Psiw[k]);
        Matrix<T> ztag = mul(phi[0], add(sub(mul(iD0, in.D[k - 1], Ak[k - 1]), Ak[0]), mul(iD0, Bk)));
        for (std::size_t i = 0; i + 1 < k; ++i)
            ztag = add(ztag, add(mul(phi[i + 1], sub(Ak[i], Ak[i + 1])),
                                 mul(mul(phi[0], iD0), mul(in.D[i], Ak[i]))));
        Matrix<T> Mx = sub(eye<T>(N), Ak[k - 1]);
        for (std::size_t i = 0; i < N; ++i) Mx(i, 0) = one;
        Matrix<T> lhs = ztag;
        lhs(0, 0) = pk;
        phi.push_back(mul(lhs, inverse(Mx)));
        md.q0[k] = mul(phi[0], iD0);
        Matrix<T> qL;
        for (std::size_t ii = 0; ii < k; ++ii) {
            const Matrix<T> a = add(sub(phi[ii + 1], phi[ii]), mul(phi[0], iD0, in.D[ii]));
            const Matrix<T> q = mul(mul(a, kron(I, in.sigma[ii])), Gi[ii]);
            qL = (ii == 0) ? q : hcat(qL, q);
        }
        md.qL[k] = qL;
    }
    return md;
}

template <class T>
std::vector<T> nppr_class(const PrioInput<T>& in, const NpprModel<T>& md, std::size_t k,
                          Measure what, std::size_t n, const std::vector<T>& pts) {
    const std::size_t N = in.N, K = in.K;
    const Matrix<T>& I = in.I;
    const T one = num_traits<T>::from_int(1);
    Matrix<T> sD0k = in.D0;
    for (std::size_t i = 0; i < k; ++i) sD0k = add(sD0k, in.D[i]);
    const bool hasz = md.Qwzz[k].rows() > 0;
    const Matrix<T> iZz = hasz ? ninv(md.Qwzz[k]) : Matrix<T>();
    Matrix<T> Kw = add(md.Qwpp[k], mul(md.Psiw[k], md.Qwmp[k]));
    if (hasz) Kw = add(Kw, mul(md.Qwpz[k], iZz, md.Qwzp[k]));

    if (k + 1 == K) {
        if (what == Measure::StMoms || what == Measure::StDistr) {
            Matrix<T> AM, BM, CM;
            for (std::size_t i = 0; i < k; ++i) {
                const Matrix<T> am = kron(ones<T>(N, 1), kron(eye<T>(in.M[i]), in.s[k]));
                AM = (i == 0) ? am : blkdiag(AM, am);
                BM = (i == 0) ? in.S[i] : blkdiag(BM, in.S[i]);
                CM = (i == 0) ? in.s[i] : vcat(CM, in.s[i]);
            }
            const Matrix<T> AMz = vcat(AM, zeros<T>(N * in.M[k], AM.cols()));
            const Matrix<T> Z = vcat(hcat(Kw, AMz), hcat(zeros<T>(BM.rows(), Kw.cols()), BM));
            const Matrix<T> z =
                vcat(vcat(zeros<T>(AM.rows(), 1), kron(ones<T>(N, 1), in.s[k])), CM);
            const Matrix<T> iniw =
                hcat(add(mul(md.q0[k], md.Qwmp[k]), mul(md.qL[k], md.Qwzp[k])),
                     zeros<T>(1, BM.rows()));
            const Matrix<T> zeta = scal(iniw, T(one / total(mul(iniw, ninv(Z), z))));
            return ph_sojourn(zeta, Z, z, what, n, pts);
        }
        const std::size_t sM = in.msum(0, K), c0 = N * in.msum(0, k);
        Matrix<T> L = zeros<T>(N * sM, N * sM), B = zeros<T>(N * sM, N * sM),
                  F = zeros<T>(N * sM, N * sM);
        std::size_t kix = 0;
        for (std::size_t i = 0; i < K; ++i) {
            const Matrix<T> IMi = eye<T>(in.M[i]);
            setblk(F, kix, kix, kron(in.D[k], IMi));
            setblk(L, kix, kix, add(kron(sD0k, IMi), kron(I, in.S[i])));
            const Matrix<T> blk = kron(I, mul(in.s[i], in.sigma[k]));
            if (i + 1 < K)
                setblk(L, kix, c0, blk);
            else
                setblk(B, kix, c0, blk);
            kix += N * in.M[i];
        }
        const Matrix<T> R = qbd_R_of(in, B, L, F);
        const Matrix<T> p0 = hcat(md.qL[k], mul(md.q0[k], kron(I, in.sigma[k])));
        return qbd_measure(R, p0, what, n);
    }

    // step 4.1 workload right before the arrivals of class k
    Matrix<T> BM, CM, DM;
    for (std::size_t i = 0; i < k; ++i) {
        const Matrix<T> bm = kron(I, in.S[i]), cm = kron(I, in.s[i]),
                        dm = kron(in.D[k], eye<T>(in.M[i]));
        BM = (i == 0) ? bm : blkdiag(BM, bm);
        CM = (i == 0) ? cm : vcat(CM, cm);
        DM = (i == 0) ? dm : blkdiag(DM, dm);
    }
    Matrix<T> Kwu = Kw, Bwu = mul(md.Psiw[k], in.D[k]), iniw, pwu;
    if (k > 0) {
        const Matrix<T> top =
            mul(add(md.Qwpz[k], mul(md.Psiw[k], md.Qwmz[k])), iZz, DM);
        Kwu = vcat(hcat(Kw, top), hcat(zeros<T>(BM.rows(), Kw.cols()), BM));
        Bwu = vcat(Bwu, CM);
        iniw = hcat(add(mul(md.q0[k], md.Qwmp[k]), mul(md.qL[k], md.Qwzp[k])), mul(md.qL[k], DM));
        pwu = mul(md.q0[k], in.D[k]);
    } else {
        iniw = mul(md.pm, md.Qwmp[k]);
        pwu = mul(md.pm, in.D[k]);
    }
    const T nrm = T(total(pwu) + total(mul(iniw, ninv(Kwu), Bwu)));
    pwu = scal(pwu, T(one / nrm));
    iniw = scal(iniw, T(one / nrm));

    // step 4.2 fluid model whose first passage time is the WAITING time
    SojournModel<T> sm;
    const std::size_t KN = Kwu.rows(), ext = N * in.msum(k + 1, K);
    sm.Qspp = zeros<T>(KN + ext, KN + ext);
    sm.Qspm = zeros<T>(KN + ext, N);
    sm.Qsmp = zeros<T>(N, KN + ext);
    sm.Qsmm = add(sD0k, in.D[k]);
    std::size_t kix = KN;
    for (std::size_t i = k + 1; i < K; ++i) {
        setblk(sm.Qspp, kix, kix, kron(I, in.S[i]));
        setblk(sm.Qspm, kix, 0, kron(I, in.s[i]));
        setblk(sm.Qsmp, 0, kix, kron(in.D[i], in.sigma[i]));
        kix += N * in.M[i];
    }
    setblk(sm.Qspp, 0, 0, Kwu);
    setblk(sm.Qspm, 0, 0, Bwu);
    sm.inis = hcat(iniw, zeros<T>(1, ext));
    sm.Psis = fluid(in, sm.Qspp, sm.Qspm, sm.Qsmp, sm.Qsmm).Psi;

    const Matrix<T>& sig = in.sigma[k];
    const Matrix<T>& Sk = in.S[k];
    const std::size_t Mk = in.M[k];
    std::vector<T> out;
    if (what == Measure::StMoms) {
        const std::vector<Matrix<T>> Pn = st_moment_chain(sm, n);
        const Matrix<T> iSk = ninv(Sk);
        Matrix<T> Pnr = scal(sig, total(mul(sm.inis, Pn[0])));
        for (std::size_t m = 1; m <= n; ++m) {
            T w = T(total(mul(sm.inis, Pn[m])) / num_traits<T>::from_double(std::pow(2.0, (double)m)));
            if (m % 2 == 1) w = T(-w);
            const Matrix<T> P =
                add(scal(mul(Pnr, iSk), num_traits<T>::from_int((int)m)), scal(sig, w));
            Pnr = P;
            out.push_back(T(total(P) + total(pwu) * factorial<T>(m) * total(mul(sig, mpow(iSk, m)))));
        }
    } else if (what == Measure::StDistr) {
        const std::size_t L = in.erl;
        for (const T& t : pts) {
            const T lambdae = T(num_traits<T>::from_int((int)L) / t / num_traits<T>::from_int(2));
            const std::vector<T> terms = erlang_terms(in, sm, lambdae);
            // tail[m] = 1 - sum(sigma * inv(I - S/(2 lambdae))^m), m = 1..L
            const Matrix<T> G = inverse(sub(eye<T>(Mk), scal(Sk, T(one / (num_traits<T>::from_int(2) * lambdae)))));
            std::vector<T> tail(L + 1, num_traits<T>::from_int(0));
            Matrix<T> v = sig;
            for (std::size_t m = 1; m <= L; ++m) {
                v = mul(v, G);
                tail[m] = T(one - total(v));
            }
            T pr = T((total(pwu) + terms[0]) * tail[L]);
            for (std::size_t m = 1; m + 1 <= L; ++m) pr += T(terms[m] * tail[L - m]);
            out.push_back(pr);
        }
    } else {
        const Matrix<T> IMk = eye<T>(Mk);
        const Matrix<T> G = ninv(add(kron(sub(in.sD, in.D[k]), IMk), kron(I, Sk)));
        const Matrix<T> W = mul(G, kron(in.D[k], IMk));
        const Matrix<T> iW = inverse(sub(eye<T>(W.rows()), W));
        const Matrix<T> w = kron(I, sig);
        const Matrix<T> omega = mul(G, kron(I, in.s[k]));
        if (what == Measure::NcMoms) {
            const Matrix<T> A = add(sm.Qspp, mul(sm.Psis, sm.Qsmp));
            const Matrix<T> B = add(sm.Qsmm, mul(sm.Qsmp, sm.Psis));
            const SylvesterFactor<T> F(A, B);
            std::vector<Matrix<T>> Psii{sm.Psis};
            std::vector<Matrix<T>> QP{mul(sm.Qsmp, sm.Psis)};
            std::vector<Matrix<T>> QLDPn{mul(mul(sm.inis, sm.Psis), w, iW)};
            for (std::size_t m = 1; m <= n; ++m) {
                Matrix<T> C = scal(mul(Psii[m - 1], in.D[k]), T(num_traits<T>::from_int((int)m)));
                T bino = one;
                for (std::size_t i = 1; i + 1 <= m; ++i) {
                    bino = T(bino * num_traits<T>::from_int((int)(m - i + 1)) /
                             num_traits<T>::from_int((int)i));
                    C = add(C, scal(mul(Psii[i], QP[m - i]), bino));
                }
                Psii.push_back(F.solve_lyap(C));
                QP.push_back(mul(sm.Qsmp, Psii.back()));
                QLDPn.push_back(add(scal(mul(QLDPn[m - 1], iW, W), num_traits<T>::from_int((int)m)),
                                    mul(mul(sm.inis, Psii.back()), w, iW)));
            }
            for (std::size_t m = 0; m <= n; ++m)
                QLDPn[m] = mul(add(QLDPn[m], mul(mul(pwu, w), mpow(iW, m + 1), mpow(W, m))), omega);
            // random_time_moms uses the model's own class-k rate, lambda(k)
            out = random_time_moms(in, k, QLDPn, n);
        } else {
            const Matrix<T> Psid = fluid(in, sm.Qspp, sm.Qspm, sm.Qsmp, sD0k).Psi;
            const Matrix<T> A = add(sm.Qspp, mul(Psid, sm.Qsmp));
            const Matrix<T> B = add(sD0k, mul(sm.Qsmp, Psid));
            const SylvesterFactor<T> F(A, B);
            std::vector<Matrix<T>> P{Psid};
            std::vector<Matrix<T>> QP{mul(sm.Qsmp, Psid)};
            Matrix<T> XDn = mul(sm.inis, Psid, w);
            const Matrix<T> pw = mul(pwu, w);
            std::vector<Matrix<T>> dql{mul(add(XDn, pw), omega)};
            Matrix<T> Wn = eye<T>(W.rows());
            for (std::size_t m = 1; m < n; ++m) {
                Matrix<T> C = mul(P[m - 1], in.D[k]);
                for (std::size_t i = 1; i + 1 <= m; ++i) C = add(C, mul(P[i], QP[m - i]));
                P.push_back(F.solve_lyap(C));
                QP.push_back(mul(sm.Qsmp, P.back()));
                XDn = add(mul(XDn, W), mul(sm.inis, P.back(), w));
                Wn = mul(Wn, W);
                dql.push_back(mul(add(XDn, mul(pw, Wn)), omega));
            }
            out = random_time_probs(in, k, dql);
        }
    }
    return out;
}

template <class T>
std::vector<std::vector<T>> run(bool preemptive, const Mmap<T>& arrival,
                                const std::vector<PhService<T>>& svc, Measure what,
                                std::size_t n, const std::vector<T>& pts,
                                const PrioQueueOptions& opt) {
    const char* who = preemptive ? "mmapph1prpr" : "mmapph1nppr";
    const PrioInput<T> in = prepare(arrival, svc, opt, who);
    if (what == Measure::StDistr) {
        for (const T& t : pts)
            if (!(t > num_traits<T>::from_int(0)))
                throw InputError(std::string(who) +
                                 ": sojourn CDF points must be strictly positive (the Erlangization "
                                 "divides by t)");
    } else if (n == 0) {
        throw InputError(std::string(who) + ": at least one moment or level is required");
    }
    std::vector<std::vector<T>> out(in.K);
    if (preemptive) {
        for (std::size_t k = 0; k < in.K; ++k) out[k] = prpr_class(in, k, what, n, pts);
    } else {
        if (in.K < 2)
            throw InputError("mmapph1nppr: at least two classes are required (with one class the "
                             "queue is MMAPPH1FCFS)");
        const NpprModel<T> md = nppr_build(in);
        for (std::size_t k = 0; k < in.K; ++k) out[k] = nppr_class(in, md, k, what, n, pts);
    }
    return out;
}

}  // namespace prio_detail

/** Per-class moments 1..n of the number of jobs, MMAP[K]/PH[K]/1 preemptive resume priority. */
template <class T>
std::vector<std::vector<T>> mmapph1prpr_ncmoms(const Mmap<T>& arrival,
                                               const std::vector<PhService<T>>& svc,
                                               std::size_t n,
                                               const PrioQueueOptions& opt = PrioQueueOptions()) {
    return prio_detail::run(true, arrival, svc, prio_detail::Measure::NcMoms, n, std::vector<T>(), opt);
}

/** Per-class P(number of jobs = 0..nmax-1), preemptive resume priority. */
template <class T>
std::vector<std::vector<T>> mmapph1prpr_ncdistr(const Mmap<T>& arrival,
                                                const std::vector<PhService<T>>& svc,
                                                std::size_t nmax,
                                                const PrioQueueOptions& opt = PrioQueueOptions()) {
    return prio_detail::run(true, arrival, svc, prio_detail::Measure::NcDistr, nmax, std::vector<T>(), opt);
}

/** Per-class sojourn-time moments 1..n, preemptive resume priority. */
template <class T>
std::vector<std::vector<T>> mmapph1prpr_stmoms(const Mmap<T>& arrival,
                                               const std::vector<PhService<T>>& svc,
                                               std::size_t n,
                                               const PrioQueueOptions& opt = PrioQueueOptions()) {
    return prio_detail::run(true, arrival, svc, prio_detail::Measure::StMoms, n, std::vector<T>(), opt);
}

/** Per-class sojourn-time CDF at the (strictly positive) points, preemptive resume priority. */
template <class T>
std::vector<std::vector<T>> mmapph1prpr_stdistr(const Mmap<T>& arrival,
                                                const std::vector<PhService<T>>& svc,
                                                const std::vector<T>& points,
                                                const PrioQueueOptions& opt = PrioQueueOptions()) {
    return prio_detail::run(true, arrival, svc, prio_detail::Measure::StDistr, 0, points, opt);
}

/** Per-class moments 1..n of the number of jobs, MMAP[K]/PH[K]/1 non-preemptive priority. */
template <class T>
std::vector<std::vector<T>> mmapph1nppr_ncmoms(const Mmap<T>& arrival,
                                               const std::vector<PhService<T>>& svc,
                                               std::size_t n,
                                               const PrioQueueOptions& opt = PrioQueueOptions()) {
    return prio_detail::run(false, arrival, svc, prio_detail::Measure::NcMoms, n, std::vector<T>(), opt);
}

/** Per-class P(number of jobs = 0..nmax-1), non-preemptive priority. */
template <class T>
std::vector<std::vector<T>> mmapph1nppr_ncdistr(const Mmap<T>& arrival,
                                                const std::vector<PhService<T>>& svc,
                                                std::size_t nmax,
                                                const PrioQueueOptions& opt = PrioQueueOptions()) {
    return prio_detail::run(false, arrival, svc, prio_detail::Measure::NcDistr, nmax, std::vector<T>(), opt);
}

/** Per-class sojourn-time moments 1..n, non-preemptive priority. */
template <class T>
std::vector<std::vector<T>> mmapph1nppr_stmoms(const Mmap<T>& arrival,
                                               const std::vector<PhService<T>>& svc,
                                               std::size_t n,
                                               const PrioQueueOptions& opt = PrioQueueOptions()) {
    return prio_detail::run(false, arrival, svc, prio_detail::Measure::StMoms, n, std::vector<T>(), opt);
}

/** Per-class sojourn-time CDF at the (strictly positive) points, non-preemptive priority. */
template <class T>
std::vector<std::vector<T>> mmapph1nppr_stdistr(const Mmap<T>& arrival,
                                                const std::vector<PhService<T>>& svc,
                                                const std::vector<T>& points,
                                                const PrioQueueOptions& opt = PrioQueueOptions()) {
    return prio_detail::run(false, arrival, svc, prio_detail::Measure::StDistr, 0, points, opt);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMAPPH1PRIO_H
