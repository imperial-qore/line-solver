/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The BMAP/PH/N/N bufferless retrial queue. Oracles, in order of strength:
 *   1. A CONSERVATION LAW. With no batch rejection and no abandonment nothing
 *      is ever lost, so the throughput must equal the arrival rate and the
 *      mean number of busy servers must be lambda b1 EXACTLY, whatever the
 *      retrial rate or the admission threshold. That is an identity of the
 *      model, independent of the algorithm.
 *   2. AN INDEPENDENT CONSTRUCTION. For a single-phase service, a single BMAP
 *      state and unit batches the chain is a two-dimensional birth-death
 *      process that can be written down directly. The test builds it from the
 *      model description, solves it in EXACT rational arithmetic, and requires
 *      the port's stationary law to agree entry by entry. Two independent
 *      generators, one exact answer.
 *   3. MATLAB qsys_bmapphnn_retrial, on eight instances covering the retrial
 *      rate, abandonment, batch rejection, the admission threshold, PH service
 *      and a genuine two-state BMAP with batches.
 *   4. The truncation diagnostics, which are what the reference lacks.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/qsys/qsys_bmapphnn_retrial.h"
#include "line/util/lu.h"

using line::Matrix;
using line::Rational;
using line::qsys::BmapPhNnRetrialOptions;
using line::qsys::BmapPhNnRetrialResult;
using line::qsys::qsys_bmapphnn_retrial;

namespace {

constexpr double TOL = 1e-10;

/** Poisson(2) arrivals as a one-state BMAP with unit batches. */
template <class T>
std::vector<Matrix<T>> poisson2() {
    std::vector<Matrix<T>> D;
    D.push_back(Matrix<T>(1, 1, line::num_traits<T>::from_int(-2)));
    D.push_back(Matrix<T>(1, 1, line::num_traits<T>::from_int(2)));
    return D;
}

/** Exp(1) service as a one-phase PH. */
template <class T>
void exp1(std::vector<T>& beta, Matrix<T>& S) {
    beta.assign(1, line::num_traits<T>::from_int(1));
    S = Matrix<T>(1, 1, line::num_traits<T>::from_int(-1));
}

/**
 * The M/M/N/N retrial chain written down directly: state (i, n) with i in
 * orbit and n busy servers, level i occupying N+1 consecutive entries.
 *
 *   arrival   lambda        n < N -> (i, n+1);  n = N -> (i+1, N) at rate
 *                           lambda (1 - p), the rest being lost
 *   service   n mu          -> (i, n-1)
 *   retrial   i alpha       -> (i-1, n+1), but ONLY while n <= R
 *   abandon   i gamma       -> (i-1, n)
 *
 * Transitions out of the top level are dropped, which is exactly what the
 * reference's truncation does. Returns the stationary law of the truncated
 * chain, level-major.
 */
template <class T>
std::vector<T> direct_chain(const T& lambda, const T& mu, int N, const T& alpha, const T& gamma,
                            const T& p, long R, std::size_t levels) {
    using nt = line::num_traits<T>;
    const T zero = nt::from_int(0), one = nt::from_int(1);
    const std::size_t w = static_cast<std::size_t>(N) + 1;
    const std::size_t total = (levels + 1) * w;
    Matrix<T> Q(total, total, zero);
    for (std::size_t i = 0; i <= levels; ++i) {
        const T iT = nt::from_int(static_cast<long>(i));
        for (int n = 0; n <= N; ++n) {
            const std::size_t from = i * w + static_cast<std::size_t>(n);
            if (n < N)
                Q(from, i * w + static_cast<std::size_t>(n + 1)) += lambda;
            else if (i < levels)
                Q(from, (i + 1) * w + static_cast<std::size_t>(N)) += T(lambda * T(one - p));
            if (n > 0) Q(from, i * w + static_cast<std::size_t>(n - 1)) += nt::from_int(n) * mu;
            if (i >= 1) {
                if (static_cast<long>(n) <= R && n < N)
                    Q(from, (i - 1) * w + static_cast<std::size_t>(n + 1)) += iT * alpha;
                Q(from, (i - 1) * w + static_cast<std::size_t>(n)) += iT * gamma;
            }
        }
    }
    for (std::size_t i = 0; i < total; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < total; ++j)
            if (j != i) s += Q(i, j);
        Q(i, i) = -s;
    }
    for (std::size_t i = 0; i < total; ++i) Q(i, total - 1) = one;
    std::vector<T> rhs(total, zero);
    rhs[total - 1] = one;
    std::vector<T> pi = line::solve(Q.transpose(), rhs);
    T mass = zero;
    for (const T& v : pi) mass += v;
    for (T& v : pi) v = v / mass;
    return pi;
}

}  // namespace

TEST_CASE("qsys_bmapphnn_retrial matches MATLAB on the documented example") {
    std::vector<double> beta;
    Matrix<double> S;
    exp1(beta, S);
    const BmapPhNnRetrialResult<double> r =
        qsys_bmapphnn_retrial(poisson2<double>(), beta, S, 3, 0.5, 0.0, 0.0, 2L);
    CHECK(r.truncLevel == 150);  // max(100, ceil(50/(1 - 2/3)))
    CHECK(r.L_orbit == doctest::Approx(3.24306844547676).epsilon(TOL));
    CHECK(r.N_server == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(r.L_system == doctest::Approx(5.24306844547676).epsilon(TOL));
    CHECK(r.utilization == doctest::Approx(0.666666666666667).epsilon(1e-12));
    CHECK(r.throughput == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(r.P_idle == doctest::Approx(0.0809747243921922).epsilon(TOL));
    CHECK(r.P_empty_orbit == doctest::Approx(0.203430267326525).epsilon(TOL));
    CHECK(r.P_empty_system == doctest::Approx(0.0352113075849004).epsilon(TOL));
    CHECK(!r.clipped);
    CHECK(r.topLevelMass < 1e-12);  // this instance really has converged
}

TEST_CASE("qsys_bmapphnn_retrial matches MATLAB across the parameter set") {
    std::vector<double> beta;
    Matrix<double> S;
    exp1(beta, S);
    const std::vector<Matrix<double>> D = poisson2<double>();
    struct Case {
        const char* name;
        double alpha, gamma, p;
        long R;
        double L_orbit, N_server, P_idle;
    };
    const std::vector<Case> cases{
        {"alpha = 5", 5.0, 0.0, 0.0, 2, 1.1516474221525, 2.00000000000002, 0.101995350851283},
        {"alpha = 100", 100.0, 0.0, 0.0, 2, 0.903101396820811, 1.99999999999998,
         0.110570258728197},
        {"gamma = 0.3", 0.5, 0.3, 0.0, 2, 0.919091004937286, 1.72427269851879,
         0.128928333772151},
        {"p = 0.4", 0.5, 0.0, 0.4, 2, 1.17093318593371, 1.77418719108606, 0.118059943272096},
    };
    for (const Case& c : cases) {
        INFO(c.name);
        const BmapPhNnRetrialResult<double> r =
            qsys_bmapphnn_retrial(D, beta, S, 3, c.alpha, c.gamma, c.p, c.R);
        CHECK(r.L_orbit == doctest::Approx(c.L_orbit).epsilon(TOL));
        CHECK(r.N_server == doctest::Approx(c.N_server).epsilon(TOL));
        CHECK(r.P_idle == doctest::Approx(c.P_idle).epsilon(TOL));
    }

    // PH service: Erlang-2 of mean 1/3 at two servers. MATLAB values.
    std::vector<double> b2{1.0, 0.0};
    Matrix<double> S2{{-6.0, 6.0}, {0.0, -6.0}};
    const BmapPhNnRetrialResult<double> rp =
        qsys_bmapphnn_retrial(D, b2, S2, 2, 0.5, 0.0, 0.0, 1L);
    CHECK(rp.truncLevel == 100);
    CHECK(rp.L_orbit == doctest::Approx(0.753387862239275).epsilon(TOL));
    CHECK(rp.N_server == doctest::Approx(0.666666666666666).epsilon(1e-12));
    CHECK(rp.P_empty_system == doctest::Approx(0.283088688819676).epsilon(TOL));

    // A genuine BMAP: two states, batches of size one and two, with
    // abandonment and batch rejection at once. MATLAB values.
    std::vector<Matrix<double>> DD;
    DD.push_back(Matrix<double>{{-3.0, 0.5}, {0.2, -2.0}});
    DD.push_back(Matrix<double>{{1.0, 0.5}, {0.3, 0.5}});
    DD.push_back(Matrix<double>{{0.5, 0.5}, {0.5, 0.5}});
    const BmapPhNnRetrialResult<double> rb =
        qsys_bmapphnn_retrial(DD, beta, Matrix<double>(1, 1, -1.5), 3, 0.7, 0.1, 0.2, 2L);
    CHECK(rb.truncLevel == 159);
    CHECK(rb.L_orbit == doctest::Approx(1.92015976584849).epsilon(TOL));
    CHECK(rb.N_server == doctest::Approx(1.72539220773172).epsilon(TOL));
    CHECK(rb.L_system == doctest::Approx(3.64555197358021).epsilon(TOL));
    CHECK(rb.utilization == doctest::Approx(0.575130735910573).epsilon(TOL));
    CHECK(rb.throughput == doctest::Approx(2.58808831159758).epsilon(TOL));
    CHECK(rb.P_empty_orbit == doctest::Approx(0.320663784474886).epsilon(TOL));
}

TEST_CASE("qsys_bmapphnn_retrial conserves customers when nothing is lost") {
    // With p = 0 and gamma = 0 no customer ever leaves without being served,
    // so throughput = lambda and N_server = lambda b1, whatever alpha and R
    // are. This is an identity of the model, not a property of the solver.
    std::vector<double> beta;
    Matrix<double> S;
    exp1(beta, S);
    const std::vector<Matrix<double>> D = poisson2<double>();
    for (double alpha : {0.5, 5.0, 100.0}) {
        for (long R : {0L, 1L, 2L}) {
            INFO("alpha = ", alpha, ", R = ", R);
            BmapPhNnRetrialOptions opt;
            opt.maxLevel = 60;
            const BmapPhNnRetrialResult<double> r = qsys_bmapphnn_retrial(
                D, beta, S, 3, alpha, 0.0, 0.0, std::vector<long>(1, R), opt);
            if (r.topLevelMass < 1e-10) {
                // only meaningful where the truncation has converged; R = 0 is
                // the unstable case covered by its own test below
                CHECK(r.throughput == doctest::Approx(2.0).epsilon(1e-9));
                CHECK(r.N_server == doctest::Approx(2.0).epsilon(1e-9));
            }
        }
    }
    // With abandonment the throughput must fall strictly below lambda.
    const BmapPhNnRetrialResult<double> g =
        qsys_bmapphnn_retrial(D, beta, S, 3, 0.5, 0.3, 0.0, 2L);
    CHECK(g.throughput < 2.0);
    CHECK(g.throughput == doctest::Approx(1.72427269851879).epsilon(TOL));
}

TEST_CASE("qsys_bmapphnn_retrial agrees exactly with a chain built by hand") {
    // Single-phase service, single BMAP state, unit batches: the generator is
    // small enough to write down directly, and rational data makes the
    // comparison exact rather than approximate.
    const Rational lambda(1), mu(1), alpha(1, 2), gamma(1, 5), p(1, 4);
    const int N = 2;
    const long R = 1;
    const std::size_t levels = 4;

    std::vector<Matrix<Rational>> D;
    D.push_back(Matrix<Rational>(1, 1, -lambda));
    D.push_back(Matrix<Rational>(1, 1, lambda));
    std::vector<Rational> beta(1, Rational(1));
    Matrix<Rational> S(1, 1, -mu);
    BmapPhNnRetrialOptions opt;
    opt.maxLevel = levels;
    const BmapPhNnRetrialResult<Rational> r =
        qsys_bmapphnn_retrial(D, beta, S, N, alpha, gamma, p, std::vector<long>(1, R), opt);

    const std::vector<Rational> pi = direct_chain(lambda, mu, N, alpha, gamma, p, R, levels);
    REQUIRE(r.pi.rows() == levels + 1);
    REQUIRE(r.pi.cols() == static_cast<std::size_t>(N) + 1);
    for (std::size_t i = 0; i <= levels; ++i)
        for (std::size_t n = 0; n <= static_cast<std::size_t>(N); ++n) {
            INFO("state (", i, ", ", n, ")");
            CHECK(r.pi(i, n) == pi[i * (static_cast<std::size_t>(N) + 1) + n]);
        }
    // and the derived means follow from the same law, exactly
    Rational orbit(0), busy(0);
    for (std::size_t i = 0; i <= levels; ++i)
        for (std::size_t n = 0; n <= static_cast<std::size_t>(N); ++n) {
            orbit += Rational(static_cast<long>(i)) * r.pi(i, n);
            busy += Rational(static_cast<long>(n)) * r.pi(i, n);
        }
    CHECK(r.L_orbit == orbit);
    CHECK(r.N_server == busy);
    CHECK(r.L_system == Rational(orbit + busy));
}

TEST_CASE("qsys_bmapphnn_retrial exposes the truncation the reference hides") {
    // REFERENCE DEFECT 1, qsys_bmapphnn_retrial.m line 202. With R = 0 a
    // retrial succeeds only when the system is completely empty and the orbit
    // is unstable, but the truncated chain is always positive recurrent and is
    // renormalized, so MATLAB returns a finite number with no warning. What
    // gives it away is that the answer tracks the truncation level and that
    // the mass at the top level does not decay. MATLAB reproduction, with
    // MaxLevel 100/200/400/800/1600: L_orbit 94.5, 194.0, 393.7, 793.6,
    // 1593.5, top-level mass about 0.18 throughout.
    std::vector<double> beta;
    Matrix<double> S;
    exp1(beta, S);
    const std::vector<Matrix<double>> D = poisson2<double>();

    std::vector<double> orbit, top;
    for (std::size_t lv : {20u, 40u, 80u}) {
        BmapPhNnRetrialOptions opt;
        opt.maxLevel = lv;
        const BmapPhNnRetrialResult<double> r =
            qsys_bmapphnn_retrial(D, beta, S, 3, 0.5, 0.0, 0.0, std::vector<long>(1, 0L), opt);
        orbit.push_back(r.L_orbit);
        top.push_back(r.topLevelMass);
    }
    // the orbit mean sits just under the truncation level and doubles with it
    CHECK(orbit[0] > 12.0);
    CHECK(orbit[1] > 1.8 * orbit[0]);
    CHECK(orbit[2] > 1.8 * orbit[1]);
    // and the diagnostic that says so
    for (double t : top) CHECK(t > 0.1);

    // The same model with R = 2 is stable, and the diagnostic tracks the
    // convergence instead of staying flat: the top-level mass falls by orders
    // of magnitude with the level and the orbit mean settles. At level 20 it
    // has NOT converged yet (mass 5.8e-4, L_orbit 3.21869 against the
    // converged 3.24307), which is the case a caller needs to be able to see.
    std::vector<double> stable, stableTop;
    for (std::size_t lv : {20u, 40u, 80u, 160u}) {
        BmapPhNnRetrialOptions opt;
        opt.maxLevel = lv;
        const BmapPhNnRetrialResult<double> r =
            qsys_bmapphnn_retrial(D, beta, S, 3, 0.5, 0.0, 0.0, std::vector<long>(1, 2L), opt);
        stable.push_back(r.L_orbit);
        stableTop.push_back(r.topLevelMass);
    }
    for (std::size_t i = 1; i < stableTop.size(); ++i) {
        INFO("top mass ", stableTop[i - 1], " -> ", stableTop[i]);
        CHECK(stableTop[i] < 0.01 * stableTop[i - 1]);
    }
    CHECK(std::abs(stable[0] - 3.24306844547676) > 1e-4);   // level 20, not converged
    CHECK(stable[2] == doctest::Approx(3.24306844547676).epsilon(1e-9));
    CHECK(stable[3] == doctest::Approx(3.24306844547676).epsilon(1e-12));
}

TEST_CASE("qsys_bmapphnn_retrial accepts a per-state threshold and rejects bad input") {
    std::vector<double> beta;
    Matrix<double> S;
    exp1(beta, S);
    std::vector<Matrix<double>> DD;
    DD.push_back(Matrix<double>{{-3.0, 0.5}, {0.2, -2.0}});
    DD.push_back(Matrix<double>{{2.5, 0.0}, {0.0, 1.8}});
    BmapPhNnRetrialOptions opt;
    opt.maxLevel = 30;
    // a threshold per BMAP state, between the two uniform ones it brackets
    const BmapPhNnRetrialResult<double> lo =
        qsys_bmapphnn_retrial(DD, beta, S, 3, 0.5, 0.1, 0.0, std::vector<long>(1, 0L), opt);
    const BmapPhNnRetrialResult<double> hi =
        qsys_bmapphnn_retrial(DD, beta, S, 3, 0.5, 0.1, 0.0, std::vector<long>(1, 2L), opt);
    std::vector<long> mixed;
    mixed.push_back(0);
    mixed.push_back(2);
    const BmapPhNnRetrialResult<double> mid =
        qsys_bmapphnn_retrial(DD, beta, S, 3, 0.5, 0.1, 0.0, mixed, opt);
    // a looser admission control can only shorten the orbit
    CHECK(hi.L_orbit < lo.L_orbit);
    CHECK(mid.L_orbit < lo.L_orbit);
    CHECK(mid.L_orbit > hi.L_orbit);

    CHECK_THROWS_AS(qsys_bmapphnn_retrial(std::vector<Matrix<double>>(1, Matrix<double>(1, 1, -1.0)),
                                          beta, S, 3, 0.5, 0.0, 0.0, 2L),
                    line::InputError);
    CHECK_THROWS_AS(qsys_bmapphnn_retrial(poisson2<double>(), beta, S, 0, 0.5, 0.0, 0.0, 2L),
                    line::InputError);
    std::vector<long> wrong(3, 1L);
    CHECK_THROWS_AS(qsys_bmapphnn_retrial(DD, beta, S, 3, 0.5, 0.0, 0.0, wrong, opt),
                    line::InputError);
}
