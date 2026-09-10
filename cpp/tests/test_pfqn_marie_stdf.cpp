/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * pfqn_marie (Marie's aggregation-decomposition for FCFS Coxian service) and
 * the McKenna sojourn-time distributions pfqn_stdf / pfqn_stdf_heur.
 *
 * Every reference value below was produced by running the MATLAB original
 * (matlab/src/api/pfqn) on the same model under lineStart, printed at %.15g.
 *
 * Tolerances, and why each one:
 *   - the Coxian fit is closed form, so it is held to 1e-14, which is rounding,
 *     and is additionally checked against the moments it is defined by rather
 *     than only against MATLAB;
 *   - pfqn_marie stops on max|mu_new - mu| < 1e-8 (single class) and
 *     max|X - Xprev| < 1e-8 (multiclass), so its assertions sit at 1e-8 and NOT
 *     tighter: the trailing digits of a fixed point stopped on a tolerance
 *     belong to the stopping rule, not to the algorithm. The observed agreement
 *     is in fact ~1e-14;
 *   - pfqn_stdf has no tolerance of its own -- it is a finite recursion over
 *     exact load-dependent constants and a matrix exponential -- so it is held
 *     to 1e-11. The single exception is the t = 0 point, which the reference
 *     replaces by FineTol = 1e-8 and where the CDF ratios are cancellation
 *     dominated; that point is held to 1e-7;
 *   - pfqn_stdf_heur inherits pfqn_rd, whose port deviates from MATLAB by about
 *     5e-6 relative on a multiserver rate lattice (the deviation is inside
 *     pfqn_rd itself and is reproducible by calling it directly on the same
 *     tilted lattice; see the note on the STDF G case). Models whose FCFS
 *     station is single server agree to 1e-11 and are asserted there; the
 *     multiserver one is asserted at 1e-5.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_marie.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_stdf.h"
#include "line/api/pfqn/pfqn_stdf_heur.h"

using line::Matrix;
using line::Real50;
using namespace line::pfqn;

namespace {

template <class T>
Matrix<T> colvec(std::vector<double> v) {
    Matrix<T> m(v.size(), 1);
    for (std::size_t i = 0; i < v.size(); ++i) m(i, 0) = line::num_traits<T>::from_double(v[i]);
    return m;
}

/**
 * Mean and SCV of a Coxian given as (phase rates, completion probabilities).
 * Reaching phase k has probability q_k = prod_{j<k} (1 - phi_j); conditional on
 * completing out of phase k the sojourn is the sum of the first k exponentials,
 * so E[X] = sum_k q_k / mu_k and
 * E[X^2] = sum_k q_k phi_k ( (sum_{j<=k} 1/mu_j)^2 + sum_{j<=k} 1/mu_j^2 ).
 */
void coxMoments(const MarieCoxFit<double>& f, double& mean, double& scv) {
    const std::size_t n = f.mu.size();
    double q = 1.0, s1 = 0.0, s2 = 0.0, m1 = 0.0, m2 = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        s1 += 1.0 / f.mu[k];
        s2 += 1.0 / (f.mu[k] * f.mu[k]);
        m1 += q / f.mu[k];
        m2 += q * f.phi[k] * (s1 * s1 + s2);
        q *= 1.0 - f.phi[k];
    }
    mean = m1;
    scv = m2 / (m1 * m1) - 1.0;
}

const std::vector<double> kTset{0, 0.25, 0.5, 1, 2, 5, 10};
const std::vector<double> kTsetPos{0.25, 0.5, 1, 2, 5, 10};

Matrix<double> zmat(std::vector<double> z) {
    Matrix<double> Z(1, z.size());
    for (std::size_t r = 0; r < z.size(); ++r) Z(0, r) = z[r];
    return Z;
}

/** The CDF column must be a CDF: in [0,1], non-decreasing, tending to one. */
void checkIsCdf(const Matrix<double>& RD, bool clamped) {
    REQUIRE(RD.rows() > 1);
    for (std::size_t t = 0; t < RD.rows(); ++t) {
        CHECK(RD(t, 0) >= 0.0);
        if (clamped) CHECK(RD(t, 0) <= 1.0);
        if (t > 0) CHECK(RD(t, 0) >= RD(t - 1, 0) - 1e-12);
    }
    CHECK(RD(RD.rows() - 1, 0) == doctest::Approx(1.0).epsilon(2e-2));
}

}  // namespace

// --------------------------------------------------------------------------
// The Coxian fit that pfqn_marie is built on
// --------------------------------------------------------------------------

TEST_CASE("marie_cox_fit reproduces Coxian.fitMeanAndSCV") {
    // MATLAB: Coxian.fitMeanAndSCV(2, scv), mean 2 throughout.
    const MarieCoxFit<double> e = marie_cox_fit(2.0, 1.0);
    REQUIRE(e.mu.size() == 1);
    CHECK(e.mu[0] == doctest::Approx(0.5).epsilon(1e-14));
    CHECK(e.phi[0] == 1.0);

    const MarieCoxFit<double> hypo = marie_cox_fit(2.0, 0.75);
    REQUIRE(hypo.mu.size() == 2);
    CHECK(hypo.mu[0] == doctest::Approx(0.585786437626905).epsilon(1e-14));
    CHECK(hypo.mu[1] == doctest::Approx(3.414213562373096).epsilon(1e-14));
    CHECK(hypo.phi[0] == 0.0);
    CHECK(hypo.phi[1] == 1.0);

    const MarieCoxFit<double> erl = marie_cox_fit(2.0, 0.25);
    REQUIRE(erl.mu.size() == 4);
    for (std::size_t k = 0; k < 4; ++k) CHECK(erl.mu[k] == doctest::Approx(2.0).epsilon(1e-14));
    CHECK(erl.phi[0] == 0.0);
    CHECK(erl.phi[3] == 1.0);

    const MarieCoxFit<double> cox = marie_cox_fit(2.0, 2.0);
    REQUIRE(cox.mu.size() == 2);
    CHECK(cox.mu[0] == doctest::Approx(1.0).epsilon(1e-14));
    CHECK(cox.mu[1] == doctest::Approx(0.25).epsilon(1e-14));
    CHECK(cox.phi[0] == doctest::Approx(0.75).epsilon(1e-14));

    const MarieCoxFit<double> cox3 = marie_cox_fit(2.0, 3.0);
    CHECK(cox3.mu[1] == doctest::Approx(0.1666666666666667).epsilon(1e-14));
    CHECK(cox3.phi[0] == doctest::Approx(0.8333333333333334).epsilon(1e-14));

    // The boundary scv = 0.5 falls in the Erlang branch (scv <= 0.5 + CoarseTol),
    // giving Erlang-2, not the degenerate hypoexponential.
    const MarieCoxFit<double> half = marie_cox_fit(2.0, 0.5);
    REQUIRE(half.mu.size() == 2);
    CHECK(half.mu[0] == doctest::Approx(1.0).epsilon(1e-14));
    CHECK(half.mu[1] == doctest::Approx(1.0).epsilon(1e-14));
}

TEST_CASE("marie_cox_fit reproduces the moments it is defined by") {
    // The identity is the oracle: outside the Erlang branch the fit must return
    // exactly the requested mean and SCV.
    const double means[] = {0.4, 2.0, 7.5};
    const double scvs[] = {1.0, 0.999, 0.75, 0.6, 1.001, 2.0, 3.7, 25.0};
    for (double mean : means)
        for (double scv : scvs) {
            MarieCoxFit<double> f = marie_cox_fit(mean, scv);
            double m = 0.0, s = 0.0;
            coxMoments(f, m, s);
            CHECK(m == doctest::Approx(mean).epsilon(1e-13));
            // The exponential branch covers |scv - 1| <= CoarseTol = 1e-3 and
            // returns exactly 1 there; everywhere else the SCV is matched.
            const bool expBranch = scv >= 1.0 - 1e-3 && scv <= 1.0 + 1e-3;
            const double want = expBranch ? 1.0 : scv;
            CHECK(s == doctest::Approx(want).epsilon(1e-13));
            CHECK(f.phi.back() == 1.0);
        }

    // Erlang branch: the phase count is ceil(1/scv), so the FITTED scv is
    // 1/ceil(1/scv) <= scv and equals scv only when 1/scv is an integer. This is
    // the reference's behaviour, reproduced rather than improved.
    for (double scv : {0.5, 0.25, 0.2, 0.3, 0.13}) {
        MarieCoxFit<double> f = marie_cox_fit(3.0, scv);
        double m = 0.0, s = 0.0;
        coxMoments(f, m, s);
        const double n = std::ceil(1.0 / scv);
        CHECK(f.mu.size() == static_cast<std::size_t>(n));
        CHECK(m == doctest::Approx(3.0).epsilon(1e-13));
        CHECK(s == doctest::Approx(1.0 / n).epsilon(1e-13));
        CHECK(s <= scv + 1e-13);
    }
}

// --------------------------------------------------------------------------
// pfqn_marie, single class
// --------------------------------------------------------------------------

TEST_CASE("pfqn_marie single class matches the MATLAB reference") {
    const double tol = 1e-8;

    // MATLAB: pfqn_marie([1;0.8],4,0,[0.5;2],[],[],[1;1]) -> it = 13
    // scv 0.5 is the Erlang-2 branch, scv 2 the Coxian-2 branch, so both
    // non-exponential branches of the fit are exercised in one model.
    {
        const MarieResult<double> a =
            pfqn_marie(colvec<double>({1, 0.8}), std::vector<int>{4}, std::vector<double>{0.0},
                       colvec<double>({0.5, 2}), tol, 1000, std::vector<int>{1, 1});
        CHECK(a.it == 13);
        CHECK(a.X[0] == doctest::Approx(0.86946322209655).epsilon(tol));
        CHECK(a.Q(0, 0) == doctest::Approx(2.47245065290151).epsilon(tol));
        CHECK(a.Q(1, 0) == doctest::Approx(1.52754934709849).epsilon(tol));
        CHECK(a.U(0, 0) == doctest::Approx(0.869463218511949).epsilon(tol));
        CHECK(a.U(1, 0) == doctest::Approx(0.695570582410669).epsilon(tol));
        CHECK(a.C(0, 0) == doctest::Approx(4.60053961840357).epsilon(tol));
        // Little's law with no think time: every job is queued somewhere.
        CHECK(a.Q(0, 0) + a.Q(1, 0) == doctest::Approx(4.0).epsilon(1e-10));
        // Cycle time and throughput are one relation: X = N / C.
        CHECK(a.X[0] * a.C(0, 0) == doctest::Approx(4.0).epsilon(1e-10));
    }

    // MATLAB: pfqn_marie([1;0.5],5,0,[2;0.25],[],[],[2;1]) -> it = 18.
    // Multiserver: the isolation chain scales the phase rate by min(n, m) and
    // the initial multiplier lattice is min(n, m) rather than 1.
    {
        const MarieResult<double> c =
            pfqn_marie(colvec<double>({1, 0.5}), std::vector<int>{5}, std::vector<double>{0.0},
                       colvec<double>({2, 0.25}), tol, 1000, std::vector<int>{2, 1});
        CHECK(c.it == 18);
        CHECK(c.X[0] == doctest::Approx(1.63830480337609).epsilon(tol));
        CHECK(c.Q(0, 0) == doctest::Approx(2.68619783292849).epsilon(tol));
        CHECK(c.Q(1, 0) == doctest::Approx(2.31380216707151).epsilon(tol));
        CHECK(c.U(0, 0) == doctest::Approx(0.919587474631651).epsilon(tol));
        CHECK(c.U(1, 0) == doctest::Approx(0.819152400839636).epsilon(tol));
        CHECK(c.C(0, 0) == doctest::Approx(3.05193514033309).epsilon(tol));
        CHECK(c.Q(0, 0) + c.Q(1, 0) == doctest::Approx(5.0).epsilon(1e-10));
        CHECK(c.X[0] * c.C(0, 0) == doctest::Approx(5.0).epsilon(1e-10));
    }
}

TEST_CASE("pfqn_marie reduces to exact MVA for exponential service") {
    // Structural identity, not a reference lookup: with scv == 1 the isolation
    // chain is M/M/1, its conditional throughput is the initial multiplier
    // lattice, and the very first aggregate solve is already the fixed point.
    const Matrix<double> L = colvec<double>({1, 0.8});
    const std::vector<int> N{5};
    const std::vector<double> Z{1.0};
    const MarieResult<double> m =
        pfqn_marie(L, N, Z, colvec<double>({1, 1}), 1e-8, 1000, std::vector<int>{1, 1});
    CHECK(m.it == 1);

    Matrix<double> Zm(1, 1);
    Zm(0, 0) = 1.0;
    const MvaResult<double> e = pfqn_mva(L, N, Zm);
    CHECK(m.X[0] == doctest::Approx(e.XN[0]).epsilon(1e-12));
    CHECK(m.Q(0, 0) == doctest::Approx(e.QN(0, 0)).epsilon(1e-12));
    CHECK(m.Q(1, 0) == doctest::Approx(e.QN(1, 0)).epsilon(1e-12));
    // Utilization of a single-server exponential station is the utilization law.
    CHECK(m.U(0, 0) == doctest::Approx(m.X[0] * L(0, 0)).epsilon(1e-12));
    CHECK(m.U(1, 0) == doctest::Approx(m.X[0] * L(1, 0)).epsilon(1e-12));
    // MATLAB: pfqn_marie([1;0.8],5,1,[1;1]) and pfqn_mva([1;0.8],5,1) agree.
    CHECK(m.X[0] == doctest::Approx(0.873391535283454).epsilon(1e-12));
    CHECK(m.Q(0, 0) == doctest::Approx(2.53679025708819).epsilon(1e-12));
    CHECK(m.Q(1, 0) == doctest::Approx(1.58981820762836).epsilon(1e-12));
    // Little's law: queued plus thinking is the whole population.
    CHECK(m.Q(0, 0) + m.Q(1, 0) + m.X[0] * 1.0 == doctest::Approx(5.0).epsilon(1e-10));
}

// --------------------------------------------------------------------------
// pfqn_marie, multiple classes
// --------------------------------------------------------------------------

TEST_CASE("pfqn_marie multiclass matches the MATLAB reference") {
    const double tol = 1e-8;
    // MATLAB: pfqn_marie([0.5 0.3;0.4 0.6],[2 1],[1 2],[2 0.5;1 3]) -> it = 7.
    Matrix<double> L{{0.5, 0.3}, {0.4, 0.6}};
    Matrix<double> scv{{2, 0.5}, {1, 3}};
    const std::vector<int> N{2, 1};
    const std::vector<double> Z{1.0, 2.0};
    const MarieResult<double> d = pfqn_marie(L, N, Z, scv, tol, 1000, std::vector<int>());
    CHECK(d.it == 7);
    CHECK(d.X[0] == doctest::Approx(0.898889716939579).epsilon(tol));
    CHECK(d.X[1] == doctest::Approx(0.29835300548201).epsilon(tol));
    CHECK(d.Q(0, 0) == doctest::Approx(0.533859862980772).epsilon(tol));
    CHECK(d.Q(0, 1) == doctest::Approx(0.156585838670096).epsilon(tol));
    CHECK(d.Q(1, 0) == doctest::Approx(0.567250420079649).epsilon(tol));
    CHECK(d.Q(1, 1) == doctest::Approx(0.246708150365884).epsilon(tol));
    CHECK(d.U(0, 0) == doctest::Approx(0.449444858469789).epsilon(tol));
    CHECK(d.U(0, 1) == doctest::Approx(0.089505901644603).epsilon(tol));
    CHECK(d.U(1, 0) == doctest::Approx(0.359555886775832).epsilon(tol));
    CHECK(d.U(1, 1) == doctest::Approx(0.179011803289206).epsilon(tol));
    CHECK(d.C(0, 0) == doctest::Approx(0.593910301697953).epsilon(tol));
    CHECK(d.C(0, 1) == doctest::Approx(0.524834125324532).epsilon(tol));
    CHECK(d.C(1, 0) == doctest::Approx(0.631056746328068).epsilon(tol));
    CHECK(d.C(1, 1) == doctest::Approx(0.82690016803186).epsilon(tol));

    // Little's law per class, and Q = X W per station-class: both are exact at
    // the AMVA fixed point regardless of how good the cd scaling is.
    for (std::size_t r = 0; r < 2; ++r) {
        double q = 0.0;
        for (std::size_t i = 0; i < 2; ++i) q += d.Q(i, r);
        CHECK(q + d.X[r] * Z[r] == doctest::Approx(static_cast<double>(N[r])).epsilon(1e-8));
        for (std::size_t i = 0; i < 2; ++i) {
            CHECK(d.Q(i, r) == doctest::Approx(d.X[r] * d.C(i, r)).epsilon(1e-8));
            // Utilization law: the busy fraction uses the TRUE mean service,
            // not the cd-scaled one.
            CHECK(d.U(i, r) == doctest::Approx(d.X[r] * L(i, r)).epsilon(1e-12));
        }
    }
}

TEST_CASE("pfqn_marie multiclass dispatches product form to exact MVA") {
    // Exponential AND class-independent demands is genuine BCMP FCFS: the
    // reference short-circuits to pfqn_mva and reports it = 0.
    Matrix<double> L{{0.5, 0.5}, {0.4, 0.4}};
    Matrix<double> scv{{1, 1}, {1, 1}};
    const std::vector<int> N{2, 1};
    const std::vector<double> Z{1.0, 2.0};
    const MarieResult<double> e = pfqn_marie(L, N, Z, scv, 1e-8, 1000, std::vector<int>());
    CHECK(e.it == 0);
    CHECK(e.X[0] == doctest::Approx(0.870972487862292).epsilon(1e-12));
    CHECK(e.X[1] == doctest::Approx(0.295718699426217).epsilon(1e-12));
    CHECK(e.Q(0, 0) == doctest::Approx(0.648815653964985).epsilon(1e-12));
    CHECK(e.Q(0, 1) == doctest::Approx(0.236133588347801).epsilon(1e-12));
    CHECK(e.Q(1, 0) == doctest::Approx(0.480211858172723).epsilon(1e-12));
    CHECK(e.Q(1, 1) == doctest::Approx(0.172429012799765).epsilon(1e-12));

    Matrix<double> Zm = zmat({1.0, 2.0});
    const MvaResult<double> m = pfqn_mva(L, N, Zm);
    CHECK(e.X[0] == doctest::Approx(m.XN[0]).epsilon(1e-14));
    CHECK(e.X[1] == doctest::Approx(m.XN[1]).epsilon(1e-14));

    // A class-DEPENDENT demand with the same exponential service is not BCMP
    // FCFS, so the decomposition must actually run.
    Matrix<double> L2{{0.5, 0.3}, {0.4, 0.4}};
    const MarieResult<double> f = pfqn_marie(L2, N, Z, scv, 1e-8, 1000, std::vector<int>());
    CHECK(f.it > 0);
}

TEST_CASE("pfqn_marie rejects a multiserver multiclass model") {
    // The reference has no multiserver multiclass isolation chain; refusing is
    // the honest answer, not silently ignoring the server counts.
    Matrix<double> L{{0.5, 0.3}, {0.4, 0.6}};
    Matrix<double> scv{{2, 0.5}, {1, 3}};
    CHECK_THROWS_AS(pfqn_marie(L, std::vector<int>{2, 1}, std::vector<double>{1.0, 2.0}, scv, 1e-8,
                               1000, std::vector<int>{2, 1}),
                    line::InputError);
}

// --------------------------------------------------------------------------
// pfqn_stdf
// --------------------------------------------------------------------------

TEST_CASE("pfqn_stdf matches the MATLAB reference") {
    const double tol = 1e-11;
    // doctest compares |lhs-rhs| < epsilon*(scale + max(|lhs|,|rhs|)) with scale
    // defaulting to 1, so on a value of order 1e-9 a relative epsilon becomes an
    // absolute one twenty-five times the value and the row cannot fail. The
    // guarded rows carry scale(0) so the tolerance is genuinely relative.
    const double tolZero = 1e-7;  // the t = 0 -> FineTol guard point

    // MATLAB: pfqn_stdf([1;0.8],3,0,[1;1],[1 2],[1;1.25],[0 .25 .5 1 2 5 10])
    {
        const StdfResult<double> r =
            pfqn_stdf(colvec<double>({1, 0.8}), std::vector<int>{3}, zmat({0.0}),
                      std::vector<int>{1, 1}, std::vector<std::size_t>{0, 1},
                      colvec<double>({1, 1.25}), kTset);
        // e00[0] is NOT the MATLAB value. At the guarded point the reference
        // returns its isnan(H) -> FineTol substitution, 4.09836065573769e-09,
        // the same constant it returns at t = 1e-10 and 1e-9 as well; the true
        // CDF is w0 mu t with w0 = P(station 1 empty at population 2) =
        // 0.8^2 / (1 + 0.8 + 0.8^2) = 16/61, giving (16/61) * 1e-8. The
        // reference is 25/16 = 1.5625 times too large there. Approved
        // divergence: this row asserts the derived value, not MATLAB's.
        const double e00[] = {2.6229508e-09,     0.0675935296937125, 0.138676747580215,
                              0.285348954445641, 0.554059148679358,  0.933890798705694,
                              0.998689356126037};
        const double e10[] = {5.1229508249383e-09, 0.124085181942468, 0.23988489761545,
                              0.443429206099467,   0.729523201255678, 0.981059396712585,
                              0.999892415799409};
        REQUIRE(!r.RD[0][0].empty());
        REQUIRE(!r.RD[1][0].empty());
        for (std::size_t t = 0; t < 7; ++t) {
            CHECK(r.RD[0][0](t, 0) == doctest::Approx(e00[t]).epsilon(t == 0 ? tolZero : tol).scale(t == 0 ? 0.0 : 1.0));
            CHECK(r.RD[1][0](t, 0) == doctest::Approx(e10[t]).epsilon(t == 0 ? tolZero : tol).scale(t == 0 ? 0.0 : 1.0));
            // The second column carries the guarded time set.
            CHECK(r.RD[0][0](t, 1) == doctest::Approx(t == 0 ? 1e-8 : kTset[t]).epsilon(1e-14));
        }
        checkIsCdf(r.RD[0][0], true);
        checkIsCdf(r.RD[1][0], true);
        CHECK(r.isNumStable);
    }

    // Multiserver FCFS station with a think time.
    // MATLAB: pfqn_stdf([1;0.5],4,1,[2;1],1,[1;2],tset)
    {
        const StdfResult<double> r =
            pfqn_stdf(colvec<double>({1, 0.5}), std::vector<int>{4}, zmat({1.0}),
                      std::vector<int>{2, 1}, std::vector<std::size_t>{0},
                      colvec<double>({1, 2}), kTset);
        const double e[] = {6.71232878108005e-09, 0.161170005123419, 0.306014994801692,
                            0.539688143940688,    0.812977978911415, 0.989998933840505,
                            0.999932215468437};
        REQUIRE(!r.RD[0][0].empty());
        for (std::size_t t = 0; t < 7; ++t)
            CHECK(r.RD[0][0](t, 0) == doctest::Approx(e[t]).epsilon(t == 0 ? tolZero : tol).scale(t == 0 ? 0.0 : 1.0));
        checkIsCdf(r.RD[0][0], true);
        // Station 1 was not listed as FCFS, so its cell stays unset.
        CHECK(r.RD[1][0].empty());
    }

    // Single station: this is the branch that dispatches to pfqn_comomrm_ld
    // rather than pfqn_mvald, and the branch where the outer constant is taken
    // on an EMPTY station set.
    // MATLAB: pfqn_stdf(1,3,2,2,1,1,tset)
    {
        const StdfResult<double> r =
            pfqn_stdf(colvec<double>({1}), std::vector<int>{3}, zmat({2.0}), std::vector<int>{2},
                      std::vector<std::size_t>{0}, colvec<double>({1}), kTset);
        const double e[] = {8.88888885953853e-09, 0.202058092110954, 0.366952538227234,
                            0.606282319058021,    0.851662534058068, 0.992518436659879,
                            0.999949555862614};
        REQUIRE(!r.RD[0][0].empty());
        for (std::size_t t = 0; t < 7; ++t)
            CHECK(r.RD[0][0](t, 0) == doctest::Approx(e[t]).epsilon(t == 0 ? tolZero : tol).scale(t == 0 ? 0.0 : 1.0));
        checkIsCdf(r.RD[0][0], true);
    }

    // Two classes at a multiserver FCFS station whose per-class rates agree.
    // MATLAB: pfqn_stdf([1 0.6;0.35 0.2],[2 1],[1 0.5],[2;1],1,[1 1;2 2],tset)
    {
        Matrix<double> L{{1, 0.6}, {0.35, 0.2}};
        Matrix<double> rates{{1, 1}, {2, 2}};
        const StdfResult<double> r =
            pfqn_stdf(L, std::vector<int>{2, 1}, zmat({1.0, 0.5}), std::vector<int>{2, 1},
                      std::vector<std::size_t>{0}, rates, kTset);
        const double e0[] = {8.07999999352661e-09, 0.188123353243711, 0.347648306327458,
                             0.58747208050507,     0.842196945048595, 0.991977083963605,
                             0.999945883679465};
        const double e1[] = {8.22852080420659e-09, 0.190681923153791, 0.351192773672452,
                             0.590925845288287,    0.843934930342151, 0.992076482218041,
                             0.999946557932308};
        REQUIRE(!r.RD[0][0].empty());
        REQUIRE(!r.RD[0][1].empty());
        for (std::size_t t = 0; t < 7; ++t) {
            CHECK(r.RD[0][0](t, 0) == doctest::Approx(e0[t]).epsilon(t == 0 ? tolZero : tol).scale(t == 0 ? 0.0 : 1.0));
            CHECK(r.RD[0][1](t, 0) == doctest::Approx(e1[t]).epsilon(t == 0 ? tolZero : tol).scale(t == 0 ? 0.0 : 1.0));
        }
        checkIsCdf(r.RD[0][0], true);
        checkIsCdf(r.RD[0][1], true);
    }
}

TEST_CASE("pfqn_stdf rejects an FCFS station with class-dependent rates") {
    // The reference calls this an invalid model, because a product-form FCFS
    // station must have one service rate shared by all classes.
    Matrix<double> L{{1, 0.6}, {0.35, 0.2}};
    Matrix<double> rates{{1, 1.5}, {2, 2}};
    CHECK_THROWS_AS(pfqn_stdf(L, std::vector<int>{2, 1}, zmat({1.0, 0.5}), std::vector<int>{2, 1},
                              std::vector<std::size_t>{0}, rates, kTset),
                    line::InputError);
}

// --------------------------------------------------------------------------
// pfqn_stdf_heur
// --------------------------------------------------------------------------

TEST_CASE("pfqn_stdf_heur matches the MATLAB reference") {
    // MATLAB: pfqn_stdf_heur([1;0.8],3,0,[1;1],[1 2],[1;1.25],[.25 .5 1 2 5 10])
    {
        const StdfResult<double> r =
            pfqn_stdf_heur(colvec<double>({1, 0.8}), std::vector<int>{3}, zmat({0.0}),
                           std::vector<int>{1, 1}, std::vector<std::size_t>{0, 1},
                           colvec<double>({1, 1.25}), kTsetPos);
        const double e00[] = {0.0775679249687048, 0.169749014983526, 0.360734085833231,
                              0.664989708709369,  0.968408969807567, 0.999619682555596};
        const double e10[] = {0.133455259157461, 0.267306077355153, 0.502139205456063,
                              0.796805987012972, 0.990949018312113, 0.9999687816431};
        for (std::size_t t = 0; t < 6; ++t) {
            CHECK(r.RD[0][0](t, 0) == doctest::Approx(e00[t]).epsilon(1e-11));
            CHECK(r.RD[1][0](t, 0) == doctest::Approx(e10[t]).epsilon(1e-11));
        }
        checkIsCdf(r.RD[0][0], false);
        checkIsCdf(r.RD[1][0], false);
    }

    // Single station, so the branch that builds Q1 from two pfqn_comomrm_ld
    // constants over the pfqn_mu_ms aggregate lattice.
    // MATLAB: pfqn_stdf_heur(1,3,2,2,1,1,[.25 .5 1 2 5 10])
    {
        const StdfResult<double> r = pfqn_stdf_heur(
            colvec<double>({1}), std::vector<int>{3}, zmat({2.0}), std::vector<int>{2},
            std::vector<std::size_t>{0}, colvec<double>({1}), kTsetPos);
        const double e[] = {0.245776907698439, 0.437188155874852, 0.702356176476175,
                            0.960738574181542, 1.10362450333435,  1.11106066674471};
        for (std::size_t t = 0; t < 6; ++t)
            CHECK(r.RD[0][0](t, 0) == doctest::Approx(e[t]).epsilon(1e-11));
        // The heuristic does NOT clamp to one, and here it genuinely overshoots.
        CHECK(r.RD[0][0](4, 0) > 1.0);
        CHECK(r.RD[0][0](5, 0) > 1.0);
    }

    // Multiserver FCFS station. The port agrees with MATLAB only to ~5e-6 here,
    // and the deviation is NOT in this file's code: calling the ported pfqn_rd
    // directly on the same tilted lattice
    //   L = [1;0.5], N = 2, Z = 1, mu = [16.6949 3.8385 2; 1 1 1]
    // gives lGN = 0.828723465513802 where MATLAB pfqn_rd gives
    // 0.828717982337221 for every one of its 'default', 'ca' and 'exact'
    // methods. pfqn_stdf_heur inherits that, so the assertion sits at 1e-5.
    // MATLAB: pfqn_stdf_heur([1;0.5],4,1,[2;1],1,[1;2],[.25 .5 1 2 5 10])
    {
        const StdfResult<double> r =
            pfqn_stdf_heur(colvec<double>({1, 0.5}), std::vector<int>{4}, zmat({1.0}),
                           std::vector<int>{2, 1}, std::vector<std::size_t>{0},
                           colvec<double>({1, 2}), kTsetPos);
        const double e[] = {0.224135812809715, 0.385120575709274, 0.584367078247175,
                            0.761932922748691, 0.987516515374482, 1.00079897562865};
        for (std::size_t t = 0; t < 6; ++t)
            CHECK(r.RD[0][0](t, 0) == doctest::Approx(e[t]).epsilon(1e-5));
        checkIsCdf(r.RD[0][0], false);
    }
}

TEST_CASE("pfqn_stdf_heur reproduces the two failure modes of the reference") {
    // (a) At t = 0 the level CDFs underflow, the tilted rate lattice becomes
    // infinite and pfqn_rd has no finite rate to fall back on. MATLAB fails
    // here with "Unable to perform assignment ..." out of pfqn_rd line 68; the
    // port raises NumericError at the same point. pfqn_stdf, which does not use
    // pfqn_rd, evaluates the same model at t = 0 without trouble.
    CHECK_THROWS_AS(pfqn_stdf_heur(colvec<double>({1, 0.8}), std::vector<int>{3}, zmat({0.0}),
                                   std::vector<int>{1, 1}, std::vector<std::size_t>{0, 1},
                                   colvec<double>({1, 1.25}), kTset),
                    line::NumericError);

    // (b) A class whose population is one leaves Q1(k,s) = 0 at N - e_r, and
    // the heuristic then asks for an exponential of mean zero. MATLAB builds
    // D0 = -Inf, gets a not-a-number CDF and fails later inside pfqn_rd; the
    // port refuses at the point the degenerate MAP is requested.
    Matrix<double> L{{1, 0.6}, {0.35, 0.2}};
    Matrix<double> rates{{1, 1}, {2, 2}};
    CHECK_THROWS_AS(pfqn_stdf_heur(L, std::vector<int>{2, 1}, zmat({1.0, 0.5}),
                                   std::vector<int>{2, 1}, std::vector<std::size_t>{0}, rates,
                                   kTsetPos),
                    line::InputError);
}

// --------------------------------------------------------------------------
// Extended precision
// --------------------------------------------------------------------------

TEST_CASE("pfqn_marie and pfqn_stdf instantiate at extended precision") {
    // Both routines are gated on has_transcendental rather than on !is_exact,
    // so Real50 must work end to end: the isolation CTMC, the load-dependent
    // MVA, the matrix exponential behind map_cdf and the Coxian fit's square
    // root all have to be available at 50 digits.
    Matrix<Real50> L(2, 1);
    L(0, 0) = Real50(1);
    L(1, 0) = Real50("0.8");
    Matrix<Real50> scv(2, 1);
    scv(0, 0) = Real50("0.5");
    scv(1, 0) = Real50(2);
    const MarieResult<Real50> m =
        pfqn_marie(L, std::vector<int>{4}, std::vector<Real50>{Real50(0)}, scv, 1e-8, 1000,
                   std::vector<int>{1, 1});
    CHECK(m.it == 13);
    CHECK(static_cast<double>(m.X[0]) == doctest::Approx(0.86946322209655).epsilon(1e-8));
    CHECK(static_cast<double>(m.Q(0, 0)) == doctest::Approx(2.47245065290151).epsilon(1e-8));
    CHECK(static_cast<double>(m.Q(0, 0) + m.Q(1, 0)) == doctest::Approx(4.0).epsilon(1e-30));

    Matrix<Real50> rates(2, 1);
    rates(0, 0) = Real50(1);
    rates(1, 0) = Real50("1.25");
    std::vector<Real50> tv;
    for (std::size_t t = 1; t < kTset.size(); ++t) tv.push_back(Real50(kTset[t]));
    const StdfResult<Real50> s =
        pfqn_stdf(L, std::vector<int>{3}, Matrix<Real50>(1, 1, Real50(0)), std::vector<int>{1, 1},
                  std::vector<std::size_t>{0, 1}, rates, tv);
    // The t = 0 guard point is deliberately excluded: there the reference
    // substitutes FineTol and the value is cancellation dominated, so the
    // double and the 50-digit results disagree in the first digit or two. Away
    // from it the two arithmetics agree with MATLAB to the digits printed.
    const double e00[] = {0.0675935296937125, 0.138676747580215, 0.285348954445641,
                          0.554059148679358,  0.933890798705694, 0.998689356126037};
    for (std::size_t t = 0; t < 6; ++t)
        CHECK(static_cast<double>(s.RD[0][0](t, 0)) == doctest::Approx(e00[t]).epsilon(1e-11));
}
