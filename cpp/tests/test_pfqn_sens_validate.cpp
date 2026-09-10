/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Ports of the pfqn_sens_* SELF-CHECK HARNESSES, not of API functions.
 *
 * matlab/src/api/pfqn/pfqn_sens_validate.m and pfqn_sens_mva_validate.m are
 * `function foo()` with no inputs and no outputs: they generate models, compare
 * the routine under test against independent references, print a PASS/FAIL
 * table and raise once on a tolerance breach. Nothing a caller could consume,
 * so they belong here and are NOT registered.
 *
 * TWO DELIBERATE DEPARTURES FROM THE REFERENCE, both strengthening.
 *
 * 1. MODEL SET. MATLAB draws its models from rng(12345) / rng(0), a stream this
 *    port cannot reproduce and should not pretend to. The checks are identities
 *    and comparisons against internal references, not against MATLAB numbers, so
 *    ANY model set validates them; a deterministic LCG is used here and the
 *    models are fixed, so these remain hard regression guards. Every tolerance
 *    and every band is carried over verbatim.
 * 2. TOLERANCE, where the arithmetic allows. MATLAB compares the brute-force
 *    enumeration against pfqn_sens_mva at 1e-9 because both are computed in
 *    floating point. At Rational both are exact, so the port asserts EQUALITY of
 *    fractions with no tolerance at all. That is the sharpest form of the same
 *    statement, and it is the reason the enumerator was worth porting.
 */
#include <algorithm>
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_sens.h"
#include "line/api/pfqn/pfqn_linearizer.h"
#include "line/api/pfqn/pfqn_mvams.h"
#include "line/api/pfqn/pfqn_sens_mvaldmx.h"
#include "line/api/pfqn/pfqn_sens_respt.h"
#include "line/api/pfqn/pfqn_sens_linearizer.h"
#include "line/api/pfqn/pfqn_sens_mom.h"
#include "line/api/pfqn/pfqn_sens_mva.h"

using line::Matrix;
using line::Rational;
using line::num_traits;
using line::pfqn::pfqn_mva;
using line::pfqn::pfqn_sens;
using line::pfqn::pfqn_linearizer;
using line::pfqn::pfqn_sens_mvaldmx;
using line::pfqn::pfqn_sens_respt;
using line::pfqn::pfqn_sens_linearizer;
using line::pfqn::pfqn_sens_mom;
using line::pfqn::pfqn_sens_mva;

namespace {

/** Deterministic model generator, standing in for MATLAB's rng(seed). */
struct Lcg {
    unsigned s;
    explicit Lcg(unsigned seed) : s(seed) {}
    unsigned next() {
        s = s * 1103515245u + 12345u;
        return (s >> 16) & 0x7fffu;
    }
    /** A demand in [0.2, 1.2] on a 1/100 lattice, so it is exact at Rational. */
    template <class T>
    T demand() {
        return num_traits<T>::from_rational(20 + static_cast<long>(next() % 101u), 100);
    }
    int range(int lo, int hi) { return lo + static_cast<int>(next() % static_cast<unsigned>(hi - lo + 1)); }
};

/** MATLAB's relerr: max |a-b| / max(1,|a|,|b|), elementwise over a vector. */
double relerr_scaled(double a, double b) {
    double scale = 1.0;
    if (std::fabs(a) > scale) scale = std::fabs(a);
    if (std::fabs(b) > scale) scale = std::fabs(b);
    return std::fabs(a - b) / scale;
}

// -------------------------------------------------------------------------
// pfqn_sens_validate: Qplus and the central finite difference
// -------------------------------------------------------------------------

/**
 * Q_{i,s}^{+k}(N - 1_r): queue lengths at the ORIGINAL M stations when a
 * replica of station k is added, at population N with one class-r job removed.
 */
template <class T>
Matrix<T> Qplus(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z, std::size_t k,
                std::size_t r) {
    const std::size_t M = L.rows(), R = L.cols();
    std::vector<int> Np = N;
    Np[r] -= 1;
    Matrix<T> Lp(M + 1, R);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t s = 0; s < R; ++s) Lp(i, s) = L(i, s);
    for (std::size_t s = 0; s < R; ++s) Lp(M, s) = L(k, s);
    const line::pfqn::MvaResult<T> q = pfqn_mva(Lp, Np, Z);
    Matrix<T> out(M, R);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t s = 0; s < R; ++s) out(i, s) = q.QN(i, s);
    return out;
}

/** Central finite difference of Q_{k,r}(N) with respect to L(i,s). */
double fd_dQ(const Matrix<double>& L, const std::vector<int>& N, const Matrix<double>& Z,
             std::size_t k, std::size_t r, std::size_t i, std::size_t s) {
    const double h = 1e-6 * std::max(1.0, std::fabs(L(i, s)));
    Matrix<double> Lp = L, Lm = L;
    Lp(i, s) += h;
    Lm(i, s) -= h;
    return (pfqn_mva(Lp, N, Z).QN(k, r) - pfqn_mva(Lm, N, Z).QN(k, r)) / (2 * h);
}

/** Index of the 'L' parameter (i,s) in a SensResult's parameter list. */
template <class T>
long sens_param_L(const line::pfqn::SensResult<T>& s, std::size_t i, std::size_t cls) {
    for (std::size_t p = 0; p < s.params.size(); ++p)
        if (s.params[p].type == 'L' && s.params[p].station == static_cast<int>(i) &&
            s.params[p].cls == cls)
            return static_cast<long>(p);
    return -1;
}

// -------------------------------------------------------------------------
// pfqn_sens_mva_validate: the brute-force enumerator
// -------------------------------------------------------------------------

/** All nonnegative integer vectors of length M summing to at most n. */
std::vector<std::vector<int>> compositions_leq(int n, std::size_t M) {
    std::vector<std::vector<int>> out;
    if (M == 1) {
        for (int v = 0; v <= n; ++v) out.push_back(std::vector<int>(1, v));
        return out;
    }
    for (int first = 0; first <= n; ++first) {
        const std::vector<std::vector<int>> sub = compositions_leq(n - first, M - 1);
        for (std::size_t j = 0; j < sub.size(); ++j) {
            std::vector<int> row;
            row.reserve(M);
            row.push_back(first);
            row.insert(row.end(), sub[j].begin(), sub[j].end());
            out.push_back(row);
        }
    }
    return out;
}

/**
 * Exact moments of the closed product-form equilibrium distribution by
 * enumeration. Stations 1..M are single-server fixed rate; the think time is an
 * infinite server carrying no moment. Port of brute_moments in
 * pfqn_sens_mva_validate.m, with one change: the reference accumulates the
 * weight in LOG space (gammaln / exp) purely for floating-point safety, and the
 * weight is a plain ratio of products, so it is formed directly here and stays
 * in the field of the inputs.
 */
template <class T>
void brute_moments(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                   Matrix<T>& Q, std::vector<Matrix<T>>& QCov) {
    const std::size_t M = L.rows(), R = L.cols();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    std::vector<std::vector<std::vector<int>>> per(R);
    for (std::size_t r = 0; r < R; ++r) per[r] = compositions_leq(N[r], M);

    std::vector<Matrix<T>> states;  // each an M x R occupancy
    std::vector<T> w;
    std::vector<std::size_t> idx(R, 0);
    bool more = true;
    T wsum = zero;
    while (more) {
        Matrix<T> nir(M, R, zero);
        std::vector<std::vector<int>> ni(M, std::vector<int>(R, 0));
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t i = 0; i < M; ++i) {
                ni[i][r] = per[r][idx[r]][i];
                nir(i, r) = num_traits<T>::from_int(ni[i][r]);
            }
        T weight = one;
        bool ok = true;
        for (std::size_t i = 0; i < M && ok; ++i) {
            int tot = 0;
            for (std::size_t r = 0; r < R; ++r) tot += ni[i][r];
            weight *= line::num_factorial<T>(static_cast<unsigned>(tot));
            for (std::size_t r = 0; r < R; ++r) {
                if (ni[i][r] <= 0) continue;
                if (!(L(i, r) > zero)) {
                    ok = false;
                    break;
                }
                weight *= line::num_pow_int(L(i, r), static_cast<unsigned>(ni[i][r])) /
                          line::num_factorial<T>(static_cast<unsigned>(ni[i][r]));
            }
        }
        if (ok) {
            for (std::size_t r = 0; r < R; ++r) {
                int used = 0;
                for (std::size_t i = 0; i < M; ++i) used += ni[i][r];
                const int n0r = N[r] - used;
                if (n0r <= 0) continue;
                if (Z.empty() || !(Z[r] > zero)) {
                    ok = false;
                    break;
                }
                weight *= line::num_pow_int(Z[r], static_cast<unsigned>(n0r)) /
                          line::num_factorial<T>(static_cast<unsigned>(n0r));
            }
        }
        states.push_back(nir);
        w.push_back(ok ? weight : zero);
        wsum += w.back();

        std::size_t d = R;
        while (d-- > 0) {
            if (++idx[d] < per[d].size()) break;
            idx[d] = 0;
        }
        more = (d != static_cast<std::size_t>(-1));
    }

    Q = Matrix<T>(M, R, zero);
    QCov.assign(M, Matrix<T>(R, R, zero));
    for (std::size_t k = 0; k < states.size(); ++k) {
        const T pk = w[k] / wsum;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) Q(i, r) += pk * states[k](i, r);
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t s = 0; s < R; ++s) {
                T m2 = zero;
                for (std::size_t k = 0; k < states.size(); ++k)
                    m2 += w[k] / wsum * states[k](i, r) * states[k](i, s);
                QCov[i](r, s) = m2 - Q(i, r) * Q(i, s);
            }
}

}  // namespace

TEST_CASE("pfqn_sens_validate: the moment identities and the demand derivatives") {
    // Port of pfqn_sens_validate.m. The four models of the reference are kept
    // verbatim where they are literal; the two it draws from rand() are replaced
    // by fixed lattice demands.
    struct Model {
        Matrix<double> L;
        std::vector<int> N;
        Matrix<double> Z;
    };
    std::vector<Model> models;
    {
        Model m;
        m.L = Matrix<double>{{0.6, 0.3}, {0.4, 0.7}, {0.5, 0.2}};
        m.N = std::vector<int>{3, 2};
        m.Z = Matrix<double>{{0.0, 0.0}};
        models.push_back(m);
    }
    {
        Model m;
        m.L = Matrix<double>{{1.0, 0.5}, {0.3, 0.9}, {0.7, 0.6}};
        m.N = std::vector<int>{2, 3};
        m.Z = Matrix<double>{{0.4, 0.8}};
        models.push_back(m);
    }
    {
        Model m;
        m.L = Matrix<double>{{0.55, 0.80, 0.35}, {0.90, 0.25, 0.60}, {0.45, 0.70, 1.05}, {0.30, 0.95, 0.50}};
        m.N = std::vector<int>{2, 2, 2};
        m.Z = Matrix<double>{{0.0, 0.0, 0.0}};
        models.push_back(m);
    }
    {
        Model m;
        m.L = Matrix<double>{{0.35, 0.85}, {1.05, 0.20}, {0.50, 0.65}};
        m.N = std::vector<int>{4, 3};
        m.Z = Matrix<double>{{0.5, 0.0}};
        models.push_back(m);
    }

    // MATLAB's tolerances, carried over verbatim.
    const double tolIdent = 1e-8;  // pure MVA identities
    const double tolDeriv = 1e-6;  // against the analytic sensitivities
    const double tolFD = 1e-4;     // analytic against central finite differences
    double maxErr[7] = {0, 0, 0, 0, 0, 0, 0};

    for (std::size_t mm = 0; mm < models.size(); ++mm) {
        const Matrix<double>& L = models[mm].L;
        const std::vector<int>& N = models[mm].N;
        const Matrix<double>& Z = models[mm].Z;
        const std::size_t M = L.rows(), R = L.cols();
        std::vector<double> Zv(R);
        for (std::size_t r = 0; r < R; ++r) Zv[r] = Z(0, r);

        const Matrix<double> QN = pfqn_mva(L, N, Z).QN;
        const line::pfqn::SensResult<double> sens = pfqn_sens(L, N, Zv);

        // cache Q^{+k}(N - 1_r)
        std::vector<std::vector<Matrix<double>>> Qp(M, std::vector<Matrix<double>>(R));
        for (std::size_t k = 0; k < M; ++k)
            for (std::size_t r = 0; r < R; ++r)
                if (N[r] >= 1) Qp[k][r] = Qplus(L, N, Z, k, r);

        for (std::size_t r = 0; r < R; ++r) {
            if (N[r] < 1) continue;
            for (std::size_t s = 0; s < R; ++s) {
                for (std::size_t k = 0; k < M; ++k) {
                    for (std::size_t i = 0; i < M; ++i) {
                        if (N[s] >= 1) {
                            // identity 1: Q_{i,s}^{+k}(N-1_r) Q_{k,r}
                            //           = Q_{k,r}^{+i}(N-1_s) Q_{i,s}
                            maxErr[0] = std::max(
                                maxErr[0],
                                relerr_scaled(Qp[k][r](i, s) * QN(k, r), Qp[i][s](k, r) * QN(i, s)));
                            // identity 2 (the superscript is +k, forced by 1 and 3)
                            maxErr[1] = std::max(
                                maxErr[1],
                                relerr_scaled(Qp[k][r](i, s) * QN(k, r) * L(k, s) * L(i, r),
                                              Qp[k][s](i, r) * QN(k, s) * L(i, s) * L(k, r)));
                            // identity 3
                            maxErr[2] = std::max(
                                maxErr[2],
                                relerr_scaled(Qp[k][r](i, s) * QN(k, r) * L(k, s) * L(i, r),
                                              Qp[i][r](k, s) * QN(i, r) * L(i, s) * L(k, r)));
                        }
                        const long p = sens_param_L(sens, i, s);
                        REQUIRE(p >= 0);
                        const double dAna = L(i, s) * sens.dQ[static_cast<std::size_t>(p)](k, r);
                        const double dFD = L(i, s) * fd_dQ(L, N, Z, k, r, i, s);
                        maxErr[6] = std::max(maxErr[6], relerr_scaled(dAna, dFD));

                        if (i == k && r == s) {
                            maxErr[3] = std::max(
                                maxErr[3],
                                relerr_scaled(dAna, QN(k, r) * (1 + 2 * Qp[k][r](k, r) - QN(k, r))));
                        } else if (i == k) {
                            maxErr[4] = std::max(
                                maxErr[4],
                                relerr_scaled(dAna, QN(k, r) * (2 * Qp[k][r](k, s) - QN(k, s))));
                        } else {
                            maxErr[5] = std::max(
                                maxErr[5], relerr_scaled(dAna, QN(k, r) * (Qp[k][r](i, s) - QN(i, s))));
                        }
                    }
                }
            }
        }
    }

    INFO("pfqn_sens identities: ident " << maxErr[0] << "/" << maxErr[1] << "/" << maxErr[2]
                                        << ", deriv " << maxErr[3] << "/" << maxErr[4] << "/"
                                        << maxErr[5] << ", FD " << maxErr[6]);
    CHECK(maxErr[0] <= tolIdent);  // moment identity 1
    CHECK(maxErr[1] <= tolIdent);  // moment identity 2
    CHECK(maxErr[2] <= tolIdent);  // moment identity 3
    CHECK(maxErr[3] <= tolDeriv);  // derivative i=k, r=s
    CHECK(maxErr[4] <= tolDeriv);  // derivative i=k, r!=s
    CHECK(maxErr[5] <= tolDeriv);  // derivative i!=k
    CHECK(maxErr[6] <= tolFD);     // analytic against finite differences
}

TEST_CASE("pfqn_sens_validate: the three moment identities hold EXACTLY at Rational") {
    // The identities are polynomial in the demands, so at Rational they are not
    // "within 1e-8", they are equalities. MATLAB cannot make that statement.
    Matrix<Rational> L(3, 2);
    const long num[3][2] = {{3, 5}, {2, 7}, {1, 4}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t s = 0; s < 2; ++s) L(i, s) = num_traits<Rational>::from_rational(num[i][s], 10);
    std::vector<int> N(2);
    N[0] = 3;
    N[1] = 2;
    Matrix<Rational> Z(1, 2);
    Z(0, 0) = num_traits<Rational>::from_rational(1, 2);
    Z(0, 1) = num_traits<Rational>::from_rational(1, 4);

    const Matrix<Rational> QN = pfqn_mva(L, N, Z).QN;
    std::vector<std::vector<Matrix<Rational>>> Qp(3, std::vector<Matrix<Rational>>(2));
    for (std::size_t k = 0; k < 3; ++k)
        for (std::size_t r = 0; r < 2; ++r) Qp[k][r] = Qplus(L, N, Z, k, r);

    for (std::size_t r = 0; r < 2; ++r)
        for (std::size_t s = 0; s < 2; ++s)
            for (std::size_t k = 0; k < 3; ++k)
                for (std::size_t i = 0; i < 3; ++i) {
                    CHECK(Rational(Qp[k][r](i, s) * QN(k, r)) == Rational(Qp[i][s](k, r) * QN(i, s)));
                    CHECK(Rational(Qp[k][r](i, s) * QN(k, r) * L(k, s) * L(i, r)) ==
                          Rational(Qp[k][s](i, r) * QN(k, s) * L(i, s) * L(k, r)));
                    CHECK(Rational(Qp[k][r](i, s) * QN(k, r) * L(k, s) * L(i, r)) ==
                          Rational(Qp[i][r](k, s) * QN(i, r) * L(i, s) * L(k, r)));
                }
}

TEST_CASE("pfqn_sens_mva_validate: brute force, pfqn_sens Jacobian, pfqn_mva, symmetry") {
    // Port of pfqn_sens_mva_validate.m. Tolerances verbatim from the reference.
    const double tolBrute = 1e-9, tolSens = 1e-9, tolMva = 1e-10, tolSym = 1e-9;
    double errBrute = 0, errSens = 0, errMva = 0, errSym = 0;
    int nBrute = 0, nSens = 0;

    Lcg rng(20260722u);
    for (int trial = 0; trial < 40; ++trial) {
        const std::size_t M = static_cast<std::size_t>(rng.range(1, 3));
        const std::size_t R = static_cast<std::size_t>(rng.range(1, 3));
        Matrix<double> L(M, R);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) L(i, r) = rng.demand<double>();
        std::vector<int> N(R);
        bool any = false;
        for (std::size_t r = 0; r < R; ++r) {
            N[r] = rng.range(0, 3);
            if (N[r] > 0) any = true;
        }
        if (!any) N[0] = 2;
        std::vector<double> Z(R, 0.0);
        if (trial % 2 == 0)
            for (std::size_t r = 0; r < R; ++r) Z[r] = 0.3 + 0.01 * (rng.next() % 70u);
        // exercise a zero-demand column now and then
        if (trial % 5 == 0 && M > 1) L(0, 0) = 0.0;

        const line::pfqn::SensMvaResult<double> mom = pfqn_sens_mva(L, N, Z);

        // ---- C. base measures against pfqn_mva ---------------------------
        Matrix<double> Zm(1, R);
        for (std::size_t r = 0; r < R; ++r) Zm(0, r) = Z[r];
        const line::pfqn::MvaResult<double> mva = pfqn_mva(L, N, Zm);
        for (std::size_t r = 0; r < R; ++r)
            errMva = std::max(errMva, relerr_scaled(mom.XN[r], mva.XN[r]));
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                errMva = std::max(errMva, relerr_scaled(mom.QN(i, r), mva.QN(i, r)));
                errMva = std::max(errMva, relerr_scaled(mom.UN(i, r), mva.UN(i, r)));
                errMva = std::max(errMva, relerr_scaled(mom.CN(i, r), mva.CN(i, r)));
            }

        // ---- A. brute force ----------------------------------------------
        std::size_t lattice = 1;
        for (std::size_t r = 0; r < R; ++r) lattice *= static_cast<std::size_t>(N[r] + 1);
        if (lattice <= 64 && M <= 3) {
            Matrix<double> Qb;
            std::vector<Matrix<double>> QCovb;
            brute_moments(L, N, Z, Qb, QCovb);
            for (std::size_t i = 0; i < M; ++i) {
                for (std::size_t r = 0; r < R; ++r) {
                    errBrute = std::max(errBrute, relerr_scaled(mom.QN(i, r), Qb(i, r)));
                    for (std::size_t s = 0; s < R; ++s)
                        errBrute =
                            std::max(errBrute, relerr_scaled(mom.QCov[i](r, s), QCovb[i](r, s)));
                }
            }
            ++nBrute;
        }

        // ---- B. the raw pfqn_sens Jacobian, and the raw asymmetry ---------
        // Read the reference off the Jacobian, NOT off sens.QCov: pfqn_sens
        // sources its same-station blocks from pfqn_sens_mva, so comparing
        // against sens.QCov would compare the recursion with itself.
        const line::pfqn::SensResult<double> sens = pfqn_sens(L, N, Z);
        std::vector<Matrix<double>> CovRef(M, Matrix<double>(R, R, 0.0));
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r)
                for (std::size_t s = 0; s < R; ++s) {
                    const long p = sens_param_L(sens, i, s);
                    if (p >= 0) CovRef[i](r, s) = L(i, s) * sens.dQ[static_cast<std::size_t>(p)](i, r);
                }
        for (std::size_t i = 0; i < M; ++i) {
            double tot = 0;
            for (std::size_t r = 0; r < R; ++r)
                for (std::size_t s = 0; s < R; ++s) {
                    errSens = std::max(errSens, relerr_scaled(mom.QCov[i](r, s), CovRef[i](r, s)));
                    tot += CovRef[i](r, s);
                }
            // Theorem 3: the variance of the total queue length at a station
            errSens = std::max(errSens, relerr_scaled(mom.QTotVar[i], tot));
        }
        errSym = std::max(errSym, mom.QCovAsym);
        ++nSens;
    }

    INFO("pfqn_sens_mva: brute " << errBrute << " (" << nBrute << " models), sens " << errSens
                                 << " (" << nSens << "), mva " << errMva << ", asym " << errSym);
    CHECK(errBrute <= tolBrute);
    CHECK(errSens <= tolSens);
    CHECK(errMva <= tolMva);
    CHECK(errSym <= tolSym);
}

TEST_CASE("pfqn_sens_mva equals the brute-force enumeration EXACTLY at Rational") {
    // MATLAB holds this at 1e-9 because both sides are floating point. Both are
    // exact here, so the comparison needs no tolerance: the covariance recursion
    // of de Souza e Silva and Muntz and a direct sum over the state space return
    // the same fractions.
    Lcg rng(4242u);
    int checked = 0;
    for (int trial = 0; trial < 10; ++trial) {
        const std::size_t M = static_cast<std::size_t>(rng.range(1, 3));
        const std::size_t R = static_cast<std::size_t>(rng.range(1, 2));
        Matrix<Rational> L(M, R);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) L(i, r) = rng.demand<Rational>();
        std::vector<int> N(R);
        for (std::size_t r = 0; r < R; ++r) N[r] = rng.range(1, 3);
        std::vector<Rational> Z(R, num_traits<Rational>::from_int(0));
        if (trial % 2 == 0)
            for (std::size_t r = 0; r < R; ++r)
                Z[r] = num_traits<Rational>::from_rational(30 + static_cast<long>(rng.next() % 70u), 100);

        const line::pfqn::SensMvaResult<Rational> mom = pfqn_sens_mva(L, N, Z);
        Matrix<Rational> Qb;
        std::vector<Matrix<Rational>> QCovb;
        brute_moments(L, N, Z, Qb, QCovb);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                CHECK(mom.QN(i, r) == Qb(i, r));
                for (std::size_t s = 0; s < R; ++s) CHECK(mom.QCov[i](r, s) == QCovb[i](r, s));
            }
        ++checked;
    }
}

// ===========================================================================
// pfqn_sens_mom_validate: brute_totals / brute_perclass, and Strelen's table
// ===========================================================================

namespace {

/** States and normalized weights of the closed product form, shared by the
 *  three enumerators of pfqn_sens_mom_validate.m. Exact in the field of T. */
template <class T>
void brute_states(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                  std::vector<Matrix<T>>& states, std::vector<T>& w) {
    const std::size_t M = L.rows(), R = L.cols();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<std::vector<std::vector<int>>> per(R);
    for (std::size_t r = 0; r < R; ++r) per[r] = compositions_leq(N[r], M);
    states.clear();
    w.clear();
    std::vector<std::size_t> idx(R, 0);
    bool more = true;
    T wsum = zero;
    while (more) {
        Matrix<T> nir(M, R, zero);
        std::vector<std::vector<int>> ni(M, std::vector<int>(R, 0));
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t i = 0; i < M; ++i) {
                ni[i][r] = per[r][idx[r]][i];
                nir(i, r) = num_traits<T>::from_int(ni[i][r]);
            }
        T weight = one;
        bool ok = true;
        for (std::size_t i = 0; i < M && ok; ++i) {
            int tot = 0;
            for (std::size_t r = 0; r < R; ++r) tot += ni[i][r];
            weight *= line::num_factorial<T>(static_cast<unsigned>(tot));
            for (std::size_t r = 0; r < R; ++r) {
                if (ni[i][r] <= 0) continue;
                if (!(L(i, r) > zero)) {
                    ok = false;
                    break;
                }
                weight *= line::num_pow_int(L(i, r), static_cast<unsigned>(ni[i][r])) /
                          line::num_factorial<T>(static_cast<unsigned>(ni[i][r]));
            }
        }
        if (ok)
            for (std::size_t r = 0; r < R; ++r) {
                int used = 0;
                for (std::size_t i = 0; i < M; ++i) used += ni[i][r];
                const int n0r = N[r] - used;
                if (n0r <= 0) continue;
                if (Z.empty() || !(Z[r] > zero)) {
                    ok = false;
                    break;
                }
                weight *= line::num_pow_int(Z[r], static_cast<unsigned>(n0r)) /
                          line::num_factorial<T>(static_cast<unsigned>(n0r));
            }
        states.push_back(nir);
        w.push_back(ok ? weight : zero);
        wsum += w.back();
        std::size_t d = R;
        while (d-- > 0) {
            if (++idx[d] < per[d].size()) break;
            idx[d] = 0;
        }
        more = (d != static_cast<std::size_t>(-1));
    }
    for (std::size_t k = 0; k < w.size(); ++k) w[k] /= wsum;
}

/** Moments of the PER-STATION TOTAL queue lengths. Port of brute_totals. */
template <class T>
void brute_totals(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                  std::vector<T>& m, std::vector<T>& Var, Matrix<T>& Cov, std::vector<T>& M2,
                  std::vector<T>& M3) {
    const std::size_t M = L.rows(), R = L.cols();
    const T zero = num_traits<T>::from_int(0);
    std::vector<Matrix<T>> states;
    std::vector<T> w;
    brute_states(L, N, Z, states, w);
    const std::size_t K = states.size();
    std::vector<std::vector<T>> tot(K, std::vector<T>(M, zero));
    for (std::size_t k = 0; k < K; ++k)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) tot[k][i] += states[k](i, r);
    m.assign(M, zero);
    M2.assign(M, zero);
    M3.assign(M, zero);
    Var.assign(M, zero);
    Cov = Matrix<T>(M, M, zero);
    for (std::size_t k = 0; k < K; ++k)
        for (std::size_t i = 0; i < M; ++i) {
            m[i] += w[k] * tot[k][i];
            M2[i] += w[k] * tot[k][i] * tot[k][i];
            M3[i] += w[k] * tot[k][i] * tot[k][i] * tot[k][i];
        }
    for (std::size_t i = 0; i < M; ++i) Var[i] = M2[i] - m[i] * m[i];
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j) {
            T acc = zero;
            for (std::size_t k = 0; k < K; ++k) acc += w[k] * tot[k][i] * tot[k][j];
            Cov(i, j) = acc - m[i] * m[j];
        }
}

/** Per-class moments of n(i,r). Port of brute_perclass. */
template <class T>
void brute_perclass(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                    Matrix<T>& m, Matrix<T>& Var, Matrix<T>& M2, Matrix<T>& M3) {
    const std::size_t M = L.rows(), R = L.cols();
    const T zero = num_traits<T>::from_int(0);
    std::vector<Matrix<T>> states;
    std::vector<T> w;
    brute_states(L, N, Z, states, w);
    m = Matrix<T>(M, R, zero);
    M2 = Matrix<T>(M, R, zero);
    M3 = Matrix<T>(M, R, zero);
    Var = Matrix<T>(M, R, zero);
    for (std::size_t k = 0; k < states.size(); ++k)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                const T n = states[k](i, r);
                m(i, r) += w[k] * n;
                M2(i, r) += w[k] * n * n;
                M3(i, r) += w[k] * n * n * n;
            }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) Var(i, r) = M2(i, r) - m(i, r) * m(i, r);
}

}  // namespace

TEST_CASE("pfqn_sens_mom_validate: brute force, total variance, base measures, symmetry") {
    // Port of pfqn_sens_mom_validate.m checks A-D and F. Tolerances verbatim.
    const double tolBrute = 1e-9, tolMva = 1e-10, tolTot = 1e-9, tolSym = 1e-9;
    double errBrute = 0, errMva = 0, errTot = 0, errSym = 0, errGrp = 0;
    int nBrute = 0, nGrp = 0;

    Lcg rng(30303u);
    for (int trial = 0; trial < 40; ++trial) {
        const std::size_t M = static_cast<std::size_t>(rng.range(1, 3));
        const std::size_t R = static_cast<std::size_t>(rng.range(1, 2));
        Matrix<double> L(M, R);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) L(i, r) = rng.demand<double>();
        std::vector<int> N(R);
        bool any = false;
        for (std::size_t r = 0; r < R; ++r) {
            N[r] = rng.range(0, 3);
            if (N[r] > 0) any = true;
        }
        if (!any) N[0] = 2;
        std::vector<double> Z(R, 0.0);
        if (trial % 2 == 0)
            for (std::size_t r = 0; r < R; ++r) Z[r] = 0.3 + 0.01 * (rng.next() % 70u);
        std::vector<int> mi(M, 1);
        if (trial % 3 == 0)
            for (std::size_t i = 0; i < M; ++i) mi[i] = rng.range(1, 3);

        const line::pfqn::SensMomResult<double> mom =
            pfqn_sens_mom(L, N, Z, mi, std::vector<int>());

        // ---- C. base measures --------------------------------------------
        Matrix<double> Zm(1, R);
        for (std::size_t r = 0; r < R; ++r) Zm(0, r) = Z[r];
        const line::pfqn::MvaResult<double> mva = pfqn_mva(L, N, Zm, mi);
        for (std::size_t r = 0; r < R; ++r)
            errMva = std::max(errMva, relerr_scaled(mom.XN[r], mva.XN[r]));
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                errMva = std::max(errMva, relerr_scaled(mom.QN(i, r), mva.QN(i, r)));
                errMva = std::max(errMva, relerr_scaled(mom.UN(i, r), mva.UN(i, r)));
                errMva = std::max(errMva, relerr_scaled(mom.CN(i, r), mva.CN(i, r)));
            }

        // ---- D. symmetry --------------------------------------------------
        errSym = std::max(errSym, mom.CovAsym);

        // ---- B. total variance against the per-class covariances ----------
        const line::pfqn::SensMvaResult<double> ref = pfqn_sens_mva(L, N, Z, mi);
        for (std::size_t i = 0; i < M; ++i)
            errTot = std::max(errTot, relerr_scaled(mom.Var(i, 0), ref.QTotVar[i]));

        // ---- A / F. brute force, and the per-class grouping ---------------
        std::size_t lattice = 1;
        for (std::size_t r = 0; r < R; ++r) lattice *= static_cast<std::size_t>(N[r] + 1);
        bool unitmi = true;
        for (std::size_t i = 0; i < M; ++i)
            if (mi[i] != 1) unitmi = false;
        if (lattice <= 32 && M <= 3 && unitmi) {
            std::vector<double> mb, Varb, M2b, M3b;
            Matrix<double> Covb;
            brute_totals(L, N, Z, mb, Varb, Covb, M2b, M3b);
            for (std::size_t i = 0; i < M; ++i) {
                errBrute = std::max(errBrute, relerr_scaled(mom.m(i, 0), mb[i]));
                errBrute = std::max(errBrute, relerr_scaled(mom.Var(i, 0), Varb[i]));
                errBrute = std::max(errBrute, relerr_scaled(mom.M2(i, 0), M2b[i]));
                errBrute = std::max(errBrute, relerr_scaled(mom.M3(i, 0), M3b[i]));
                for (std::size_t j = 0; j < M; ++j)
                    errBrute = std::max(errBrute, relerr_scaled(mom.Cov(i, j), Covb(i, j)));
            }
            ++nBrute;

            // F. groups = 1:R scales one class at a time (Akyildiz-Strelen
            // Theorem 1 with T = {r}), so it must reproduce the per-class
            // moments of the brute-force distribution INCLUDING the third, and
            // its second moments must equal pfqn_sens_mva's.
            std::vector<int> perclass(R);
            for (std::size_t r = 0; r < R; ++r) perclass[r] = static_cast<int>(r);
            const line::pfqn::SensMomResult<double> momc = pfqn_sens_mom(L, N, Z, mi, perclass);
            Matrix<double> mc, Varc, M2c, M3c;
            brute_perclass(L, N, Z, mc, Varc, M2c, M3c);
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < R; ++r) {
                    errGrp = std::max(errGrp, relerr_scaled(momc.m(i, r), mc(i, r)));
                    errGrp = std::max(errGrp, relerr_scaled(momc.Var(i, r), Varc(i, r)));
                    errGrp = std::max(errGrp, relerr_scaled(momc.M2(i, r), M2c(i, r)));
                    errGrp = std::max(errGrp, relerr_scaled(momc.M3(i, r), M3c(i, r)));
                    errGrp = std::max(errGrp, relerr_scaled(momc.Var(i, r), ref.QVar(i, r)));
                }
            ++nGrp;

            // G. two classes in one group give the moments of their sum, i.e.
            // the station totals, so it must reproduce the default grouping.
            if (R == 2) {
                const line::pfqn::SensMomResult<double> momg =
                    pfqn_sens_mom(L, N, Z, mi, std::vector<int>(2, 0));
                for (std::size_t i = 0; i < M; ++i) {
                    errGrp = std::max(errGrp, relerr_scaled(momg.m(i, 0), mom.m(i, 0)));
                    errGrp = std::max(errGrp, relerr_scaled(momg.Var(i, 0), mom.Var(i, 0)));
                    errGrp = std::max(errGrp, relerr_scaled(momg.M3(i, 0), mom.M3(i, 0)));
                }
            }
        }
    }

    INFO("pfqn_sens_mom: brute " << errBrute << " (" << nBrute << "), tot " << errTot << ", mva "
                                 << errMva << ", asym " << errSym << ", grp " << errGrp << " ("
                                 << nGrp << ")");
    CHECK(errBrute <= tolBrute);
    CHECK(errTot <= tolTot);
    CHECK(errMva <= tolMva);
    CHECK(errSym <= tolSym);
    CHECK(errGrp <= tolBrute);
}

TEST_CASE("pfqn_sens_mom reproduces Strelen's published Example 3.4 table") {
    // Check E of pfqn_sens_mom_validate.m, and the most valuable one in the
    // batch: the expected values are numbers the AUTHOR printed, not numbers any
    // LINE code produced. Kobayashi central-server model, 12 type-1 queues, one
    // class; queues 1-9 x=0.0215 e=9.333, queues 10-11 x=0.104 e=10.5,
    // queue 12 x=0.019 e=105. Reference: J. C. Strelen, "Moment Analysis for
    // Closed Queuing Networks and its Linearizer", Perf. Eval. 11:127-142, 1990.
    const double tolPaper = 5e-5;  // the paper prints five significant digits
    Matrix<double> L(12, 1);
    for (std::size_t i = 0; i < 9; ++i) L(i, 0) = 0.0215 * 9.333;
    L(9, 0) = 0.104 * 10.5;
    L(10, 0) = 0.104 * 10.5;
    L(11, 0) = 0.019 * 105.0;

    // rows: queues 1-9, 10-11, 12; columns: n = 3, 2, 1
    const double paperM[3][3] = {{0.07606, 0.05835, 0.03353},
                                 {0.53316, 0.36327, 0.18246},
                                 {1.24917, 0.74835, 0.33334}};
    const double paperVar[3][3] = {{0.07893, 0.05873, 0.03240},
                                   {0.57689, 0.34341, 0.14917},
                                   {1.02546, 0.56250, 0.22222}};
    double errPaper = 0;
    for (int col = 0; col < 3; ++col) {
        const int nJobs = 3 - col;
        std::vector<int> N(1, nJobs);
        std::vector<double> Z(1, 0.0);
        const line::pfqn::SensMomResult<double> mk = pfqn_sens_mom(L, N, Z);
        const double got[3] = {mk.m(0, 0), mk.m(9, 0), mk.m(11, 0)};
        const double gotV[3] = {mk.Var(0, 0), mk.Var(9, 0), mk.Var(11, 0)};
        for (int row = 0; row < 3; ++row) {
            errPaper = std::max(errPaper, relerr_scaled(got[row], paperM[row][col]));
            errPaper = std::max(errPaper, relerr_scaled(gotV[row], paperVar[row][col]));
        }
        // the nine identical queues must be identical, and so must 10 and 11
        for (std::size_t i = 0; i < 9; ++i)
            errPaper = std::max(errPaper, relerr_scaled(mk.m(i, 0), mk.m(0, 0)));
        errPaper = std::max(errPaper, relerr_scaled(mk.m(9, 0), mk.m(10, 0)));
    }
    INFO("Strelen Example 3.4 published table: max relative error " << errPaper);
    CHECK(errPaper <= tolPaper);
}

TEST_CASE("pfqn_sens_mom equals the brute-force enumeration EXACTLY at Rational") {
    // As with pfqn_sens_mva: MATLAB compares at 1e-9 because both sides are
    // floating point, and both are exact here. Covers E[Q], Var, E[Q^2] and
    // E[Q^3] of the per-station totals.
    Lcg rng(5150u);
    int checked = 0;
    for (int trial = 0; trial < 8; ++trial) {
        const std::size_t M = static_cast<std::size_t>(rng.range(1, 3));
        const std::size_t R = static_cast<std::size_t>(rng.range(1, 2));
        Matrix<Rational> L(M, R);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) L(i, r) = rng.demand<Rational>();
        std::vector<int> N(R);
        for (std::size_t r = 0; r < R; ++r) N[r] = rng.range(1, 3);
        std::vector<Rational> Z(R, num_traits<Rational>::from_int(0));
        if (trial % 2 == 0)
            for (std::size_t r = 0; r < R; ++r)
                Z[r] = num_traits<Rational>::from_rational(30 + static_cast<long>(rng.next() % 70u), 100);

        const line::pfqn::SensMomResult<Rational> mom = pfqn_sens_mom(L, N, Z);
        std::vector<Rational> mb, Varb, M2b, M3b;
        Matrix<Rational> Covb;
        brute_totals(L, N, Z, mb, Varb, Covb, M2b, M3b);
        for (std::size_t i = 0; i < M; ++i) {
            CHECK(mom.m(i, 0) == mb[i]);
            CHECK(mom.Var(i, 0) == Varb[i]);
            CHECK(mom.M2(i, 0) == M2b[i]);
            CHECK(mom.M3(i, 0) == M3b[i]);
            for (std::size_t j = 0; j < M; ++j) CHECK(mom.Cov(i, j) == Covb(i, j));
        }
        ++checked;
    }
}

// ===========================================================================
// pfqn_sens_linearizer_validate: accuracy bands, degenerate exactness, invariants
// ===========================================================================

TEST_CASE("pfqn_sens_linearizer_validate: accuracy bands and structural invariants") {
    // Port of pfqn_sens_linearizer_validate.m. pfqn_sens_linearizer is an
    // APPROXIMATION, so it must NOT be held to machine precision. The bands are
    // the accuracy Strelen himself reports in Section 5 over his 51 networks,
    // carried over verbatim, so this asserts the port reproduces the paper's own
    // accuracy statement rather than some slacker figure.
    const double tolExact = 1e-9;  // B: the one-job case must be exact
    const double tolPop = 1e-8;    // C: population conservation
    const double tolLin = 5e-2;    // C: agreement with pfqn_linearizer on the means
    const double bandM = 0.021;    // r(E[Q])   < 2.1% in the reference
    const double bandM2 = 0.041;   // r(E[Q^2]) < 4.1% in the reference
    const double bandM3 = 0.062;   // r(E[Q^3]) < 6.2% in the reference
    // The paper reports no error on the variance. It is naturally larger than
    // the one on E[Q^2] because Var = E[Q^2] - E[Q]^2 is a difference of larger
    // numbers; banded here only to catch regressions.
    const double bandVar = 0.08;

    double errExact = 0, errPop = 0, errLin = 0;
    double worstM = 0, worstM2 = 0, worstM3 = 0, worstVar = 0;
    int nA = 0;

    Lcg rng(70707u);
    for (int trial = 0; trial < 60; ++trial) {
        const std::size_t M = static_cast<std::size_t>(rng.range(2, 4));
        const std::size_t R = static_cast<std::size_t>(rng.range(1, 2));
        Matrix<double> L(M, R);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) L(i, r) = rng.demand<double>();
        std::vector<int> N(R);
        for (std::size_t r = 0; r < R; ++r) N[r] = rng.range(1, 4);
        std::vector<double> Z(R, 0.0);
        if (trial % 2 == 0)
            for (std::size_t r = 0; r < R; ++r) Z[r] = 0.5 + 0.01 * (rng.next() % 100u);

        const line::pfqn::SensLinearizerResult<double> app = pfqn_sens_linearizer(L, N, Z);
        const line::pfqn::SensMomResult<double> ex = pfqn_sens_mom(L, N, Z);

        for (std::size_t i = 0; i < M; ++i) {
            worstM = std::max(worstM, relerr_scaled(app.m[i], ex.m(i, 0)));
            worstM2 = std::max(worstM2, relerr_scaled(app.M2[i], ex.M2(i, 0)));
            worstM3 = std::max(worstM3, relerr_scaled(app.M3[i], ex.M3(i, 0)));
            worstVar = std::max(worstVar, relerr_scaled(app.Var[i], ex.Var(i, 0)));
        }
        ++nA;

        // ---- C. population conservation -----------------------------------
        for (std::size_t r = 0; r < R; ++r) {
            double inNet = 0;
            for (std::size_t i = 0; i < M; ++i) inNet += app.QN(i, r);
            const double inDelay = app.XN[r] * Z[r];
            errPop = std::max(errPop, std::fabs(inNet + inDelay - N[r]) /
                                          std::max(1.0, static_cast<double>(N[r])));
        }

        // ---- C. means against LINE's own Linearizer -----------------------
        Matrix<double> Zm(1, R);
        for (std::size_t r = 0; r < R; ++r) Zm(0, r) = Z[r];
        const line::pfqn::LinearizerResult<double> QL =
            pfqn_linearizer(L, N, Zm, std::vector<line::pfqn::SchedStrategy>(), 1e-10, 500,
                            Matrix<double>());
        for (std::size_t i = 0; i < M; ++i) {
            double a = 0, b = 0;
            for (std::size_t r = 0; r < R; ++r) {
                a += app.QN(i, r);
                b += QL.Q(i, r);
            }
            errLin = std::max(errLin, relerr_scaled(a, b));
        }
    }

    // ---- B. one job: the Linearizer equations degenerate to the exact MVA ---
    // At a population of one the CORE estimate of the queue lengths one job down
    // is identically zero whatever the delta terms are, so every moment must
    // match pfqn_sens_mom to roundoff. This pins the derivative algebra
    // independently of the heuristic.
    for (int trial = 0; trial < 12; ++trial) {
        const std::size_t M = static_cast<std::size_t>(rng.range(2, 4));
        Matrix<double> L(M, 1);
        for (std::size_t i = 0; i < M; ++i) L(i, 0) = rng.demand<double>();
        std::vector<double> Z(1, (trial % 2 == 0) ? 0.4 + 0.01 * (rng.next() % 100u) : 0.0);
        std::vector<int> N(1, 1);
        const line::pfqn::SensLinearizerResult<double> app = pfqn_sens_linearizer(L, N, Z);
        const line::pfqn::SensMomResult<double> ex = pfqn_sens_mom(L, N, Z);
        for (std::size_t i = 0; i < M; ++i) {
            errExact = std::max(errExact, relerr_scaled(app.m[i], ex.m(i, 0)));
            errExact = std::max(errExact, relerr_scaled(app.Var[i], ex.Var(i, 0)));
            errExact = std::max(errExact, relerr_scaled(app.M2[i], ex.M2(i, 0)));
            errExact = std::max(errExact, relerr_scaled(app.M3[i], ex.M3(i, 0)));
            for (std::size_t j = 0; j < M; ++j)
                errExact = std::max(errExact, relerr_scaled(app.Cov(i, j), ex.Cov(i, j)));
        }
    }

    INFO("pfqn_sens_linearizer (" << nA << " models): E[Q] " << 100 * worstM << "% (band 2.1%), "
                                  << "Var " << 100 * worstVar << "% (band 8%), E[Q^2] "
                                  << 100 * worstM2 << "% (band 4.1%), E[Q^3] " << 100 * worstM3
                                  << "% (band 6.2%); one-job exact " << errExact << ", pop "
                                  << errPop << ", vs pfqn_linearizer " << errLin);
    CHECK(worstM <= bandM);
    CHECK(worstM2 <= bandM2);
    CHECK(worstM3 <= bandM3);
    CHECK(worstVar <= bandVar);
    CHECK(errExact <= tolExact);
    CHECK(errPop <= tolPop);
    CHECK(errLin <= tolLin);
}

// ===========================================================================
// pfqn_sens_mvaldmx_validate: brute_ldmx, finite differences, the LI limit
// ===========================================================================

namespace {

/** MATLAB's mu_at: limited load dependence saturates at the last tabulated rate. */
template <class T>
T mu_at(const Matrix<T>& mu, std::size_t i, std::size_t j) {
    // j is one-based, as in the reference
    const std::size_t col = (j <= mu.cols()) ? j - 1 : mu.cols() - 1;
    return mu(i, col);
}

/**
 * Exact moments of the MIXED load-dependent product form by enumeration. Port
 * of brute_ldmx in pfqn_sens_mvaldmx_validate.m:
 *
 *   p(n) ~ prod_i [ n_i! prod_r a(i,r)^n(i,r)/n(i,r)! prod_{j<=n_i} 1/mu(i,j) ]
 *          * prod_{closed c} Z(c)^n(0,c)/n(0,c)!
 *
 * with a(i,r) = D(i,r) on a closed class and lambda(r) D(i,r) on an open one.
 * Closed classes are enumerated exactly; OPEN classes are TRUNCATED at Kopen
 * jobs per station, which converges geometrically, so a mixed model is checked
 * at the reference's looser tolerance and is NOT an exact oracle. A closed
 * load-dependent model (no open class) IS exact and is asserted as such at
 * Rational.
 *
 * As in brute_moments, the reference accumulates in log space for
 * floating-point safety and the weight is a plain ratio of products, so it is
 * formed directly here and stays in the field of T.
 *
 * N encodes an open class as a NEGATIVE entry, matching the port's convention;
 * MATLAB writes Inf.
 */
template <class T>
void brute_ldmx(const std::vector<T>& lambda, const Matrix<T>& D, const std::vector<int>& N,
                const std::vector<T>& Z, const Matrix<T>& mu, int Kopen, Matrix<T>& Q,
                std::vector<Matrix<T>>& QCov) {
    const std::size_t M = D.rows(), R = D.cols();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    Matrix<T> a(M, R, zero);
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t i = 0; i < M; ++i)
            a(i, r) = (N[r] < 0) ? T(lambda[r] * D(i, r)) : D(i, r);

    std::vector<std::vector<std::vector<int>>> alloc(R);
    for (std::size_t r = 0; r < R; ++r)
        alloc[r] = compositions_leq(N[r] < 0 ? Kopen : N[r], M);

    std::vector<Matrix<T>> states;
    std::vector<T> w;
    std::vector<std::size_t> idx(R, 0);
    bool more = true;
    T wsum = zero;
    while (more) {
        Matrix<T> nir(M, R, zero);
        std::vector<std::vector<int>> ni(M, std::vector<int>(R, 0));
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t i = 0; i < M; ++i) {
                ni[i][r] = alloc[r][idx[r]][i];
                nir(i, r) = num_traits<T>::from_int(ni[i][r]);
            }
        T weight = one;
        bool ok = true;
        for (std::size_t i = 0; i < M && ok; ++i) {
            int tot = 0;
            for (std::size_t r = 0; r < R; ++r) tot += ni[i][r];
            weight *= line::num_factorial<T>(static_cast<unsigned>(tot));
            for (int j = 1; j <= tot; ++j) weight /= mu_at(mu, i, static_cast<std::size_t>(j));
            for (std::size_t r = 0; r < R; ++r) {
                if (ni[i][r] <= 0) continue;
                if (!(a(i, r) > zero)) {
                    ok = false;
                    break;
                }
                weight *= line::num_pow_int(a(i, r), static_cast<unsigned>(ni[i][r])) /
                          line::num_factorial<T>(static_cast<unsigned>(ni[i][r]));
            }
        }
        if (ok)
            for (std::size_t c = 0; c < R; ++c) {
                if (N[c] < 0) continue;  // open classes carry no delay term
                int used = 0;
                for (std::size_t i = 0; i < M; ++i) used += ni[i][c];
                const int n0c = N[c] - used;
                if (n0c <= 0) continue;
                if (!(Z[c] > zero)) {
                    ok = false;
                    break;
                }
                weight *= line::num_pow_int(Z[c], static_cast<unsigned>(n0c)) /
                          line::num_factorial<T>(static_cast<unsigned>(n0c));
            }
        states.push_back(nir);
        w.push_back(ok ? weight : zero);
        wsum += w.back();
        std::size_t d = R;
        while (d-- > 0) {
            if (++idx[d] < alloc[d].size()) break;
            idx[d] = 0;
        }
        more = (d != static_cast<std::size_t>(-1));
    }

    Q = Matrix<T>(M, R, zero);
    QCov.assign(M, Matrix<T>(R, R, zero));
    for (std::size_t k = 0; k < states.size(); ++k) {
        const T pk = w[k] / wsum;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) Q(i, r) += pk * states[k](i, r);
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t s = 0; s < R; ++s) {
                T m2 = zero;
                for (std::size_t k = 0; k < states.size(); ++k)
                    m2 += w[k] / wsum * states[k](i, r) * states[k](i, s);
                QCov[i](r, s) = m2 - Q(i, r) * Q(i, s);
            }
}

}  // namespace

TEST_CASE("pfqn_sens_mvaldmx_validate: base measures, finite differences, symmetry") {
    // Port of checks A, B and E of pfqn_sens_mvaldmx_validate.m. Tolerances
    // verbatim. Open classes are marked by a NEGATIVE population here where
    // MATLAB writes Inf.
    const double tolMva = 1e-12, tolFd = 1e-6, tolSym = 1e-8;
    double errMva = 0, errFd = 0, errSym = 0;
    int nFd = 0;

    Lcg rng(818181u);
    for (int trial = 0; trial < 12; ++trial) {
        const std::size_t M = static_cast<std::size_t>(rng.range(1, 2));
        const int Ropen = rng.range(0, 1);
        const std::size_t R = static_cast<std::size_t>(1 + Ropen);
        Matrix<double> D(M, R);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) D(i, r) = 0.2 + 0.006 * (rng.next() % 101u);
        std::vector<int> N(R, 0);
        N[0] = rng.range(1, 3);
        std::vector<double> lambda(R, 0.0);
        if (Ropen == 1) {
            N[1] = -1;  // open
            lambda[1] = 0.05 + 0.0015 * (rng.next() % 101u);
        }
        std::vector<double> Z(R, 0.0);
        Z[0] = 0.005 * (rng.next() % 101u);
        int NCtot = 0;
        for (std::size_t r = 0; r < R; ++r)
            if (N[r] > 0) NCtot += N[r];
        const int b = rng.range(1, 3);
        const std::size_t cols = static_cast<std::size_t>(std::max(NCtot, 1));
        Matrix<double> mu(M, cols);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t n = 1; n <= cols; ++n)
                mu(i, n - 1) = std::min<int>(static_cast<int>(n), b) *
                               (0.8 + 0.004 * (rng.next() % 101u));
        // keep the geometric tail of the limited load dependence stable
        bool unstable = false;
        for (std::size_t i = 0; i < M; ++i) {
            double Lo = 0;
            for (std::size_t r = 0; r < R; ++r) Lo += lambda[r] * D(i, r);
            if (Lo / mu(i, cols - 1) > 0.6) unstable = true;
        }
        if (unstable) continue;

        const line::pfqn::SensMvaldmxResult<double> mom =
            pfqn_sens_mvaldmx(lambda, D, N, Z, mu);

        // ---- A. base measures against pfqn_mvaldmx ------------------------
        Matrix<double> Zm(1, R);
        for (std::size_t r = 0; r < R; ++r) Zm(0, r) = Z[r];
        const line::pfqn::MvaResult<double> ref = line::pfqn::pfqn_mvaldmx(lambda, D, N, Zm, mu);
        for (std::size_t r = 0; r < R; ++r)
            errMva = std::max(errMva, relerr_scaled(mom.XN[r], ref.XN[r]));
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                errMva = std::max(errMva, relerr_scaled(mom.QN(i, r), ref.QN(i, r)));
                errMva = std::max(errMva, relerr_scaled(mom.UN(i, r), ref.UN(i, r)));
                errMva = std::max(errMva, relerr_scaled(mom.CN(i, r), ref.CN(i, r)));
            }

        // ---- E. symmetry ---------------------------------------------------
        errSym = std::max(errSym, mom.QCovAsym);

        // ---- B. finite differences: Cov = D d nbar / dD ---------------------
        const double h = 1e-6;
        for (std::size_t j = 0; j < M; ++j)
            for (std::size_t s = 0; s < R; ++s) {
                if (!(D(j, s) > 0)) continue;
                Matrix<double> Dp = D, Dm = D;
                Dp(j, s) = D(j, s) * (1 + h);
                Dm(j, s) = D(j, s) * (1 - h);
                const Matrix<double> QNp = line::pfqn::pfqn_mvaldmx(lambda, Dp, N, Zm, mu).QN;
                const Matrix<double> QNm = line::pfqn::pfqn_mvaldmx(lambda, Dm, N, Zm, mu).QN;
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t r = 0; r < R; ++r) {
                        const double fd = (QNp(i, r) - QNm(i, r)) / (2 * h);
                        const double an = mom.QCovFull(i * R + r, j * R + s);
                        errFd = std::max(errFd, relerr_scaled(an, fd));
                    }
                ++nFd;
            }
    }

    INFO("pfqn_sens_mvaldmx: mva " << errMva << ", fd " << errFd << " (" << nFd
                                   << " params), asym " << errSym);
    CHECK(errMva <= tolMva);
    CHECK(errFd <= tolFd);
    CHECK(errSym <= tolSym);
}

TEST_CASE("pfqn_sens_mvaldmx_validate: the closed load-independent limit is pfqn_sens_mva") {
    // Check C. With mu identically one and no open class the mixed
    // load-dependent recursion must collapse onto de Souza e Silva and Muntz.
    const double tolMom = 1e-9;
    double errMom = 0;
    int nMom = 0;
    Lcg rng(24680u);
    for (int trial = 0; trial < 12; ++trial) {
        const std::size_t M = static_cast<std::size_t>(rng.range(1, 3));
        const std::size_t R = static_cast<std::size_t>(rng.range(1, 2));
        Matrix<double> D(M, R);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) D(i, r) = rng.demand<double>();
        std::vector<int> N(R);
        int tot = 0;
        for (std::size_t r = 0; r < R; ++r) {
            N[r] = rng.range(1, 3);
            tot += N[r];
        }
        std::vector<double> Z(R), lambda(R, 0.0);
        for (std::size_t r = 0; r < R; ++r) Z[r] = 0.004 * (rng.next() % 101u);
        Matrix<double> mu(M, static_cast<std::size_t>(tot), 1.0);
        const line::pfqn::SensMvaldmxResult<double> mom =
            pfqn_sens_mvaldmx(lambda, D, N, Z, mu);
        const line::pfqn::SensMvaResult<double> ref = pfqn_sens_mva(D, N, Z);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                errMom = std::max(errMom, relerr_scaled(mom.QN(i, r), ref.QN(i, r)));
                errMom = std::max(errMom, relerr_scaled(mom.QVar(i, r), ref.QVar(i, r)));
                for (std::size_t s = 0; s < R; ++s)
                    errMom = std::max(errMom, relerr_scaled(mom.QCov[i](r, s), ref.QCov[i](r, s)));
            }
        for (std::size_t i = 0; i < M; ++i)
            errMom = std::max(errMom, relerr_scaled(mom.QTotVar[i], ref.QTotVar[i]));
        ++nMom;
    }
    INFO("pfqn_sens_mvaldmx closed LI limit (" << nMom << " models): " << errMom);
    CHECK(errMom <= tolMom);
}

TEST_CASE("pfqn_sens_mvaldmx_validate: brute force, closed load-dependent and mixed") {
    // Checks D1 and D2. D1 is an EXACT enumeration; D2 truncates the open
    // classes at 60 jobs per station and converges geometrically, so the
    // reference holds it at a looser tolerance and so does this port.
    const double tolBrt = 1e-8, tolBrtT = 5e-5;
    double errBrt = 0, errBrtT = 0;
    int nBrt = 0, nBrtT = 0;
    Lcg rng(13579u);

    // D1. closed load-dependent, exact enumeration
    for (int trial = 0; trial < 10; ++trial) {
        const std::size_t M = 2, R = 1;
        Matrix<double> D(M, R);
        for (std::size_t i = 0; i < M; ++i) D(i, 0) = 0.3 + 0.005 * (rng.next() % 101u);
        std::vector<int> N(1, rng.range(2, 4));
        std::vector<double> Z(1, 0.003 * (rng.next() % 101u));
        std::vector<double> lambda(1, 0.0);
        const int b = rng.range(2, 3);
        Matrix<double> mu(M, static_cast<std::size_t>(N[0]));
        for (std::size_t i = 0; i < M; ++i)
            for (int n = 1; n <= N[0]; ++n)
                mu(i, static_cast<std::size_t>(n - 1)) =
                    std::min(n, b) * (0.8 + 0.004 * (rng.next() % 101u));
        const line::pfqn::SensMvaldmxResult<double> mom =
            pfqn_sens_mvaldmx(lambda, D, N, Z, mu);
        Matrix<double> Qb;
        std::vector<Matrix<double>> QCovb;
        brute_ldmx(lambda, D, N, Z, mu, 0, Qb, QCovb);
        for (std::size_t i = 0; i < M; ++i) {
            errBrt = std::max(errBrt, relerr_scaled(mom.QN(i, 0), Qb(i, 0)));
            errBrt = std::max(errBrt, relerr_scaled(mom.QCov[i](0, 0), QCovb[i](0, 0)));
        }
        ++nBrt;
    }

    // D2. mixed load-dependent, truncated enumeration
    for (int trial = 0; trial < 6; ++trial) {
        const std::size_t M = 2, R = 2;
        Matrix<double> D(M, R);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) D(i, r) = 0.3 + 0.004 * (rng.next() % 101u);
        std::vector<int> N(2);
        N[0] = rng.range(1, 2);
        N[1] = -1;  // open
        std::vector<double> Z(2, 0.0);
        Z[0] = 0.003 * (rng.next() % 101u);
        std::vector<double> lambda(2, 0.0);
        lambda[1] = 0.05 + 0.001 * (rng.next() % 101u);
        const int b = rng.range(1, 2);
        Matrix<double> mu(M, static_cast<std::size_t>(N[0]));
        for (std::size_t i = 0; i < M; ++i)
            for (int n = 1; n <= N[0]; ++n)
                mu(i, static_cast<std::size_t>(n - 1)) =
                    std::min(n, b) * (1.0 + 0.003 * (rng.next() % 101u));
        bool unstable = false;
        for (std::size_t i = 0; i < M; ++i) {
            double Lo = 0;
            for (std::size_t r = 0; r < R; ++r) Lo += lambda[r] * D(i, r);
            if (Lo / mu(i, mu.cols() - 1) > 0.4) unstable = true;
        }
        if (unstable) continue;
        const line::pfqn::SensMvaldmxResult<double> mom =
            pfqn_sens_mvaldmx(lambda, D, N, Z, mu);
        Matrix<double> Qb;
        std::vector<Matrix<double>> QCovb;
        brute_ldmx(lambda, D, N, Z, mu, 60, Qb, QCovb);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                errBrtT = std::max(errBrtT, relerr_scaled(mom.QN(i, r), Qb(i, r)));
                for (std::size_t s = 0; s < R; ++s)
                    errBrtT = std::max(errBrtT, relerr_scaled(mom.QCov[i](r, s), QCovb[i](r, s)));
            }
        ++nBrtT;
    }

    INFO("pfqn_sens_mvaldmx brute: closed LD " << errBrt << " (" << nBrt << " models), mixed LD "
                                               << errBrtT << " (" << nBrtT << ")");
    CHECK(errBrt <= tolBrt);
    CHECK(errBrtT <= tolBrtT);
}

TEST_CASE("pfqn_sens_mvaldmx equals the closed load-dependent enumeration EXACTLY at Rational") {
    // Only the CLOSED case: the mixed enumeration truncates the open classes, so
    // it is an approximation and asserting equality there would be wrong.
    Lcg rng(97531u);
    int checked = 0;
    for (int trial = 0; trial < 6; ++trial) {
        const std::size_t M = 2;
        Matrix<Rational> D(M, 1);
        for (std::size_t i = 0; i < M; ++i)
            D(i, 0) = num_traits<Rational>::from_rational(30 + static_cast<long>(rng.next() % 51u), 100);
        std::vector<int> N(1, rng.range(2, 3));
        std::vector<Rational> Z(1, num_traits<Rational>::from_rational(
                                       static_cast<long>(rng.next() % 31u), 100));
        std::vector<Rational> lambda(1, num_traits<Rational>::from_int(0));
        Matrix<Rational> mu(M, static_cast<std::size_t>(N[0]));
        const int b = rng.range(2, 3);
        for (std::size_t i = 0; i < M; ++i)
            for (int n = 1; n <= N[0]; ++n)
                mu(i, static_cast<std::size_t>(n - 1)) =
                    num_traits<Rational>::from_int(std::min(n, b)) *
                    num_traits<Rational>::from_rational(80 + static_cast<long>(rng.next() % 41u), 100);
        const line::pfqn::SensMvaldmxResult<Rational> mom =
            pfqn_sens_mvaldmx(lambda, D, N, Z, mu);
        Matrix<Rational> Qb;
        std::vector<Matrix<Rational>> QCovb;
        brute_ldmx(lambda, D, N, Z, mu, 0, Qb, QCovb);
        for (std::size_t i = 0; i < M; ++i) {
            CHECK(mom.QN(i, 0) == Qb(i, 0));
            CHECK(mom.QCov[i](0, 0) == QCovb[i](0, 0));
        }
        ++checked;
    }
}

// ===========================================================================
// pfqn_sens_respt_validate: brute_respt, the W = w/V identity, Strelen's table
// ===========================================================================

namespace {

/**
 * E[(W|j)^t]: a job arriving to find j jobs at an FCFS b-server station waits an
 * Erlang(max(0, j-b+1), b mu) and is then served for an Exp(mu). Port of
 * cond_moment; exact in the field of T (binomials and factorials only).
 */
template <class T>
T cond_moment(int j, int b, const T& mu, int t) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const int k = std::max(0, j - b + 1);
    T v = zero;
    for (int s = 0; s <= t; ++s) {
        const T EX = line::num_factorial<T>(static_cast<unsigned>(s)) /
                     line::num_pow_int(mu, static_cast<unsigned>(s));
        const int p = t - s;
        T EY;
        if (k == 0) {
            EY = (p == 0) ? one : zero;
        } else {
            const T theta = T(num_traits<T>::from_int(b) * mu);
            EY = one;
            for (int A = 0; A < p; ++A) EY *= num_traits<T>::from_int(k + A);
            EY /= line::num_pow_int(theta, static_cast<unsigned>(p));
        }
        v += line::num_nck<T>(t, s) * EX * EY;
    }
    return v;
}

/**
 * P[Q_i = j] at every station, by enumerating the closed product form of a
 * network of FCFS b-server stations:
 *   f_i(q_i) = q_i! prod_l a(i,l)^{q_il}/q_il! prod_{j<=q_i} 1/min(j,b_i)
 * with a(i,l) = S(i) V(i,l), plus the delay term. Port of brute_marginals.
 * Returns an (M x 1+max(1,sum N)) matrix.
 */
template <class T>
Matrix<T> brute_marginals(const std::vector<T>& S, const Matrix<T>& V, const std::vector<int>& N,
                          const std::vector<T>& Z, const std::vector<int>& b) {
    const std::size_t M = V.rows(), R = V.cols();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    Matrix<T> a(M, R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) a(i, r) = S[i] * V(i, r);

    std::vector<std::vector<std::vector<int>>> per(R);
    for (std::size_t r = 0; r < R; ++r) per[r] = compositions_leq(N[r], M);

    int sumN = 0;
    for (std::size_t r = 0; r < R; ++r) sumN += N[r];
    const std::size_t maxj = static_cast<std::size_t>(std::max(sumN, 1));
    Matrix<T> pj(M, maxj + 1, zero);
    T wsum = zero;
    std::vector<std::vector<int>> tots;
    std::vector<T> w;

    std::vector<std::size_t> idx(R, 0);
    bool more = true;
    while (more) {
        std::vector<std::vector<int>> ni(M, std::vector<int>(R, 0));
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t i = 0; i < M; ++i) ni[i][r] = per[r][idx[r]][i];
        T weight = one;
        bool ok = true;
        std::vector<int> tot(M, 0);
        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t r = 0; r < R; ++r) tot[i] += ni[i][r];
        }
        for (std::size_t i = 0; i < M && ok; ++i) {
            weight *= line::num_factorial<T>(static_cast<unsigned>(tot[i]));
            for (int j = 1; j <= tot[i]; ++j)
                weight /= num_traits<T>::from_int(std::min(j, b[i]));
            for (std::size_t r = 0; r < R; ++r) {
                if (ni[i][r] <= 0) continue;
                if (!(a(i, r) > zero)) {
                    ok = false;
                    break;
                }
                weight *= line::num_pow_int(a(i, r), static_cast<unsigned>(ni[i][r])) /
                          line::num_factorial<T>(static_cast<unsigned>(ni[i][r]));
            }
        }
        if (ok)
            for (std::size_t r = 0; r < R; ++r) {
                int used = 0;
                for (std::size_t i = 0; i < M; ++i) used += ni[i][r];
                const int n0r = N[r] - used;
                if (n0r <= 0) continue;
                if (Z.empty() || !(Z[r] > zero)) {
                    ok = false;
                    break;
                }
                weight *= line::num_pow_int(Z[r], static_cast<unsigned>(n0r)) /
                          line::num_factorial<T>(static_cast<unsigned>(n0r));
            }
        w.push_back(ok ? weight : zero);
        wsum += w.back();
        tots.push_back(tot);
        std::size_t d = R;
        while (d-- > 0) {
            if (++idx[d] < per[d].size()) break;
            idx[d] = 0;
        }
        more = (d != static_cast<std::size_t>(-1));
    }
    for (std::size_t k = 0; k < w.size(); ++k)
        for (std::size_t i = 0; i < M; ++i)
            pj(i, static_cast<std::size_t>(tots[k][i])) += w[k] / wsum;
    return pj;
}

/**
 * Sojourn-time moments from first principles: the exact arrival-theorem
 * marginals p_i(j, N - e_l) mixed over the conditional sojourn-time moments.
 * Port of brute_respt. USES NONE OF THEOREM 4.1, which is what makes it an
 * independent oracle for pfqn_sens_respt rather than a restatement of it.
 */
template <class T>
void brute_respt(const std::vector<T>& S, const Matrix<T>& V, const std::vector<int>& N,
                 const std::vector<T>& Z, const std::vector<int>& b, int tmax,
                 std::vector<Matrix<T>>& WM, Matrix<T>& pN) {
    const std::size_t M = V.rows(), R = V.cols();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    int bmax = 1;
    for (std::size_t i = 0; i < M; ++i) bmax = std::max(bmax, b[i]);
    WM.assign(static_cast<std::size_t>(tmax), Matrix<T>(M, R, zero));
    for (std::size_t l = 0; l < R; ++l) {
        if (N[l] == 0) continue;
        std::vector<int> Nl = N;
        Nl[l] -= 1;
        const Matrix<T> pj = brute_marginals(S, V, Nl, Z, b);
        for (std::size_t i = 0; i < M; ++i) {
            if (!(V(i, l) > zero)) continue;
            const T mu = one / S[i];
            for (int t = 1; t <= tmax; ++t) {
                T acc = zero;
                for (std::size_t j = 0; j < pj.cols(); ++j)
                    acc += pj(i, j) * cond_moment(static_cast<int>(j), b[i], mu, t);
                WM[static_cast<std::size_t>(t - 1)](i, l) = acc;
            }
        }
    }
    const Matrix<T> pAll = brute_marginals(S, V, N, Z, b);
    pN = Matrix<T>(M, static_cast<std::size_t>(bmax), zero);
    for (std::size_t i = 0; i < M; ++i)
        for (int j = 0; j < bmax; ++j)
            if (static_cast<std::size_t>(j) < pAll.cols())
                pN(i, static_cast<std::size_t>(j)) = pAll(i, static_cast<std::size_t>(j));
}

}  // namespace

TEST_CASE("pfqn_sens_respt_validate: brute force, the W = w/V identity, base measures") {
    // Port of checks A, C and D of pfqn_sens_respt_validate.m. Tolerances verbatim.
    const double tolBrute = 1e-9, tolIdent = 1e-10, tolMva = 1e-10;
    double errBrute = 0, errIdent = 0, errMva = 0;
    int nBrute = 0;

    Lcg rng(556677u);
    for (int trial = 0; trial < 36; ++trial) {
        const std::size_t M = static_cast<std::size_t>(rng.range(1, 3));
        const std::size_t R = static_cast<std::size_t>(rng.range(1, 2));
        std::vector<double> S(M);
        for (std::size_t i = 0; i < M; ++i) S[i] = 0.2 + 0.008 * (rng.next() % 101u);
        Matrix<double> V(M, R);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) V(i, r) = 0.3 + 0.007 * (rng.next() % 101u);
        if (trial % 4 == 0 && M > 1) V(0, 0) = 0.0;  // a class that skips a station
        std::vector<int> N(R);
        for (std::size_t r = 0; r < R; ++r) N[r] = rng.range(1, 3);
        std::vector<double> Z(R, 0.0);
        if (trial % 2 == 0)
            for (std::size_t r = 0; r < R; ++r) Z[r] = 0.3 + 0.007 * (rng.next() % 101u);
        std::vector<int> b(M, 1);
        if (trial % 3 == 0)
            for (std::size_t i = 0; i < M; ++i) b[i] = rng.range(1, 3);

        const line::pfqn::SensResptResult<double> res = pfqn_sens_respt(S, V, N, Z, b, 3);

        // ---- C. the identity W(i,l) = w_i(l) / V(i,l) ---------------------
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t l = 0; l < R; ++l)
                if (V(i, l) > 0 && N[l] > 0)
                    errIdent = std::max(errIdent,
                                        relerr_scaled(res.W(i, l), res.Wresid(i, l) / V(i, l)));

        // ---- D. base measures, single-server only -------------------------
        bool single = true;
        for (std::size_t i = 0; i < M; ++i)
            if (b[i] != 1) single = false;
        if (single) {
            Matrix<double> L(M, R);
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < R; ++r) L(i, r) = S[i] * V(i, r);
            Matrix<double> Zm(1, R);
            for (std::size_t r = 0; r < R; ++r) Zm(0, r) = Z[r];
            const line::pfqn::MvaResult<double> mva = pfqn_mva(L, N, Zm);
            for (std::size_t r = 0; r < R; ++r)
                errMva = std::max(errMva, relerr_scaled(res.XN[r], mva.XN[r]));
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < R; ++r) {
                    errMva = std::max(errMva, relerr_scaled(res.QN(i, r), mva.QN(i, r)));
                    errMva = std::max(errMva, relerr_scaled(res.UN(i, r), mva.UN(i, r)));
                }
        }

        // ---- A. brute force -----------------------------------------------
        std::size_t lattice = 1;
        for (std::size_t r = 0; r < R; ++r) lattice *= static_cast<std::size_t>(N[r] + 1);
        if (lattice <= 24 && M <= 3) {
            std::vector<Matrix<double>> Wb;
            Matrix<double> pb;
            brute_respt(S, V, N, Z, b, 3, Wb, pb);
            for (int t = 0; t < 3; ++t)
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t l = 0; l < R; ++l)
                        errBrute = std::max(
                            errBrute, relerr_scaled(res.WM[static_cast<std::size_t>(t)](i, l),
                                                    Wb[static_cast<std::size_t>(t)](i, l)));
            // .p is RAGGED: station i defines only j = 0..b_i-1, the range the
            // b-server recursion needs, and the rest is zero padding. Comparing
            // the padding would be comparing against something the algorithm
            // never claims to compute.
            for (std::size_t i = 0; i < M; ++i)
                for (int j = 0; j < b[i]; ++j)
                    errBrute = std::max(errBrute,
                                        relerr_scaled(res.p(i, static_cast<std::size_t>(j)),
                                                      pb(i, static_cast<std::size_t>(j))));
            ++nBrute;
        }
    }

    INFO("pfqn_sens_respt: brute " << errBrute << " (" << nBrute << " models), identity "
                                   << errIdent << ", mva " << errMva);
    CHECK(errBrute <= tolBrute);
    CHECK(errIdent <= tolIdent);
    CHECK(errMva <= tolMva);
}

TEST_CASE("pfqn_sens_respt_validate: a single job sees an empty station, so W ~ Exp(mu)") {
    // Check E. With one job in the network the arriving job always finds the
    // station empty, so its sojourn time is exactly the service time and the
    // moments are 1/mu, 2/mu^2, 6/mu^3 with variance 1/mu^2. Held at 1e-12 in
    // MATLAB; at Rational it is an identity, so both are asserted.
    const double tolExp = 1e-12;
    SUBCASE("double") {
        std::vector<double> S(2);
        S[0] = 0.4;
        S[1] = 0.25;
        Matrix<double> V(2, 1);
        V(0, 0) = 1;
        V(1, 0) = 2;
        std::vector<int> N(1, 1), b(2, 1);
        std::vector<double> Z(1, 0.7);
        const line::pfqn::SensResptResult<double> r1 = pfqn_sens_respt(S, V, N, Z, b, 3);
        double errExp = 0;
        for (std::size_t i = 0; i < 2; ++i) {
            const double mu = 1.0 / S[i];
            errExp = std::max(errExp, relerr_scaled(r1.WM[0](i, 0), 1 / mu));
            errExp = std::max(errExp, relerr_scaled(r1.WM[1](i, 0), 2 / (mu * mu)));
            errExp = std::max(errExp, relerr_scaled(r1.WM[2](i, 0), 6 / (mu * mu * mu)));
            errExp = std::max(errExp, relerr_scaled(r1.WVar(i, 0), 1 / (mu * mu)));
        }
        CHECK(errExp <= tolExp);
    }
    SUBCASE("exact") {
        std::vector<Rational> S(2);
        S[0] = num_traits<Rational>::from_rational(2, 5);
        S[1] = num_traits<Rational>::from_rational(1, 4);
        Matrix<Rational> V(2, 1);
        V(0, 0) = num_traits<Rational>::from_int(1);
        V(1, 0) = num_traits<Rational>::from_int(2);
        std::vector<int> N(1, 1), b(2, 1);
        std::vector<Rational> Z(1, num_traits<Rational>::from_rational(7, 10));
        const line::pfqn::SensResptResult<Rational> r1 = pfqn_sens_respt(S, V, N, Z, b, 3);
        for (std::size_t i = 0; i < 2; ++i) {
            CHECK(r1.WM[0](i, 0) == S[i]);
            CHECK(r1.WM[1](i, 0) == Rational(num_traits<Rational>::from_int(2) * S[i] * S[i]));
            CHECK(r1.WM[2](i, 0) ==
                  Rational(num_traits<Rational>::from_int(6) * S[i] * S[i] * S[i]));
            CHECK(r1.WVar(i, 0) == Rational(S[i] * S[i]));
        }
    }
}

TEST_CASE("pfqn_sens_respt reproduces Strelen's published sojourn-time table") {
    // Check B: Example 3.4 continued, the Kobayashi central-server model at
    // n = 3, where the paper prints E(W_i) and sigma^2_(W_i). As with the
    // queue-length table these are the author's numbers, so this check cannot
    // be satisfied by a consistently wrong implementation.
    const double tolPaper = 5e-4;  // the paper prints five decimals on small times
    std::vector<double> xs(12), es(12);
    for (std::size_t i = 0; i < 9; ++i) {
        xs[i] = 0.0215;
        es[i] = 9.333;
    }
    xs[9] = xs[10] = 0.104;
    es[9] = es[10] = 10.5;
    xs[11] = 0.019;
    es[11] = 105.0;
    Matrix<double> V(12, 1);
    for (std::size_t i = 0; i < 12; ++i) V(i, 0) = es[i];
    std::vector<int> N(1, 3), b(12, 1);
    std::vector<double> Z(1, 0.0);
    const line::pfqn::SensResptResult<double> rk = pfqn_sens_respt(xs, V, N, Z, b, 3);

    const double paperW[3] = {0.02275, 0.14178, 0.03322};
    const double paperWV[3] = {0.00052, 0.01846, 0.00083};
    const double gotW[3] = {rk.W(0, 0), rk.W(9, 0), rk.W(11, 0)};
    const double gotWV[3] = {rk.WVar(0, 0), rk.WVar(9, 0), rk.WVar(11, 0)};
    double errPaper = 0;
    for (int i = 0; i < 3; ++i) {
        errPaper = std::max(errPaper, relerr_scaled(gotW[i], paperW[i]));
        errPaper = std::max(errPaper, relerr_scaled(gotWV[i], paperWV[i]));
    }
    INFO("Strelen Example 3.4 published sojourn times: max relative error "
         << errPaper << "; E(W_12) paper " << paperW[2] << " got " << gotW[2]
         << ", sigma2 paper " << paperWV[2] << " got " << gotWV[2]);
    CHECK(errPaper <= tolPaper);
}

TEST_CASE("pfqn_sens_respt equals the arrival-theorem enumeration EXACTLY at Rational") {
    // The enumerator uses none of Theorem 4.1, so this is an independent route,
    // and both sides are exact: no tolerance.
    Lcg rng(31415u);
    int checked = 0;
    for (int trial = 0; trial < 6; ++trial) {
        const std::size_t M = static_cast<std::size_t>(rng.range(1, 2));
        const std::size_t R = 1;
        std::vector<Rational> S(M);
        for (std::size_t i = 0; i < M; ++i)
            S[i] = num_traits<Rational>::from_rational(20 + static_cast<long>(rng.next() % 81u), 100);
        Matrix<Rational> V(M, R);
        for (std::size_t i = 0; i < M; ++i)
            V(i, 0) = num_traits<Rational>::from_rational(30 + static_cast<long>(rng.next() % 71u), 100);
        std::vector<int> N(1, rng.range(1, 3));
        std::vector<Rational> Z(1, num_traits<Rational>::from_rational(
                                       30 + static_cast<long>(rng.next() % 51u), 100));
        std::vector<int> b(M, 1);
        if (trial % 2 == 0)
            for (std::size_t i = 0; i < M; ++i) b[i] = rng.range(1, 2);

        const line::pfqn::SensResptResult<Rational> res = pfqn_sens_respt(S, V, N, Z, b, 3);
        std::vector<Matrix<Rational>> Wb;
        Matrix<Rational> pb;
        brute_respt(S, V, N, Z, b, 3, Wb, pb);
        for (int t = 0; t < 3; ++t)
            for (std::size_t i = 0; i < M; ++i)
                CHECK(res.WM[static_cast<std::size_t>(t)](i, 0) ==
                      Wb[static_cast<std::size_t>(t)](i, 0));
        for (std::size_t i = 0; i < M; ++i)
            for (int j = 0; j < b[i]; ++j)
                CHECK(res.p(i, static_cast<std::size_t>(j)) == pb(i, static_cast<std::size_t>(j)));
        ++checked;
    }
}
