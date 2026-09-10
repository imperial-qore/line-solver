/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The Lindley-recursion family of api/qsys plus the multiclass M/M/1-PS
 * sojourn moments and the two ForkTail entry points of api/fj.
 *
 * ORACLES.
 *  (a) The values printed in the reference .m docstrings, which MATLAB
 *      produces and which every assertion below naming a literal reproduces.
 *  (b) Closed forms that hold independently of the implementation: the M/M/1
 *      Lindley step at Wn = 0 is E[max(S-A,0)] = mu/(lambda(lambda+mu)) +
 *      (lambda-mu)/(lambda mu); a degenerate hyperexponential mixture is the
 *      M/M/1 step; single-class M/M/1-PS reduces to the Coffman-Muntz-Trotter
 *      second moment 4/(mu^2 (1-rho)^2 (2-rho)); and the ForkTail fit of a
 *      branch with SCV 1 is the exponential alpha = 1, beta = E[T].
 *  (c) Sample-path identities for the tandem recursion, which no distribution
 *      enters: an empty tandem never waits, and every departure epoch is the
 *      arrival epoch plus the accumulated sojourn.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/fj/fj_mg1_respt_moments.h"
#include "line/api/fj/fj_tail_forktail.h"
#include "line/api/qsys/qsys_hh1_lindley.h"
#include "line/api/qsys/qsys_lindley_moment.h"
#include "line/api/qsys/qsys_mm1_lindley.h"
#include "line/api/qsys/qsys_mm1_ps.h"
#include "line/api/qsys/qsys_mm1_tandem_lindley.h"
#include "line/api/qsys/qsys_tandem_lindley.h"

using line::Matrix;
namespace qsys = line::qsys;
namespace fj = line::fj;

TEST_CASE("qsys_mm1_lindley reproduces the reference values") {
    // Oracle (a), read from the FUNCTION and not from its docstring: MATLAB
    // returns mean 1.8902059153 and var 1.7016068097 at Wn = 2. The docstring
    // literal 1.7025 was stale and is corrected in qsys_mm1_lindley.m; a port
    // "fixed" to reach it would be wrong (see _kb/03-api-layer.md).
    const qsys::LindleyResult<double> r = qsys::qsys_mm1_lindley(0.8, 1.0, 2.0, 3u);
    CHECK(r.mean[0] == doctest::Approx(1.8902059153).epsilon(1e-9));
    CHECK(r.var[0] == doctest::Approx(1.7016068097).epsilon(1e-9));
    CHECK(r.mmax == 3);
    CHECK(r.moments.rows() == 1);
    CHECK(r.moments.cols() == 3);

    // the explicit mean and moments(:,1) are the same quantity by two routes
    CHECK(r.moments(0, 0) == doctest::Approx(r.mean[0]).epsilon(1e-12));
    // and the variance is the second moment less the squared FIRST moment
    CHECK(r.var[0] ==
          doctest::Approx(r.moments(0, 1) - r.moments(0, 0) * r.moments(0, 0)).epsilon(1e-12));

    // Oracle (a): the empty-queue step, docstring value 0.4444
    const qsys::LindleyResult<double> r0 = qsys::qsys_mm1_lindley(0.8, 1.0, 0.0);
    CHECK(r0.mean[0] == doctest::Approx(0.4444).epsilon(1e-4));
    // Oracle (b): at Wn = 0 the closed form is (lambda-mu)/(lambda mu) + mu/(lambda(lambda+mu))
    CHECK(r0.mean[0] == doctest::Approx((0.8 - 1.0) / 0.8 + 1.0 / (0.8 * 1.8)).epsilon(1e-12));

    // mmax below 2 is raised to 2, as in MATLAB, so var always exists
    CHECK(qsys::qsys_mm1_lindley(1.5, 0.7, 1.0, 1u).mmax == 2);

    // vector Wn is evaluated elementwise
    std::vector<double> w;
    w.push_back(0.0);
    w.push_back(0.5);
    w.push_back(3.25);
    const qsys::LindleyResult<double> rv = qsys::qsys_mm1_lindley(1.5, 0.7, w);
    REQUIRE(rv.mean.size() == 3);
    for (std::size_t i = 0; i < 3; ++i)
        CHECK(rv.mean[i] ==
              doctest::Approx(qsys::qsys_mm1_lindley(1.5, 0.7, w[i]).mean[0]).epsilon(1e-14));
    // the step is monotone in the current wait and the variance is positive
    CHECK(rv.mean[0] < rv.mean[1]);
    CHECK(rv.mean[1] < rv.mean[2]);
    for (std::size_t i = 0; i < 3; ++i) CHECK(rv.var[i] > 0.0);
}

TEST_CASE("qsys_lindley_moment is the shared kernel of the two Lindley entries") {
    std::vector<double> w;
    w.push_back(0.0);
    w.push_back(1.75);
    for (unsigned m = 1; m <= 4; ++m) {
        const std::vector<double> v = qsys::qsys_lindley_moment(0.9, 1.3, w, m);
        REQUIRE(v.size() == 2);
        CHECK(v[0] > 0.0);
        CHECK(v[1] > v[0]);  // a longer current wait gives a larger raw moment
    }
    // at Wn = 0 the m-th moment is lambda mu/(lambda+mu) * m!/mu^(m+1)
    // + (-1)^m m! (1 - 1)/lambda^(m+1), i.e. the positive-support half alone
    const double lambda = 0.9, mu = 1.3;
    for (unsigned m = 1; m <= 3; ++m) {
        double fact = 1.0;
        for (unsigned k = 2; k <= m; ++k) fact *= k;
        const double expected = lambda * mu / (lambda + mu) * fact / std::pow(mu, m + 1.0);
        CHECK(qsys::qsys_lindley_moment(lambda, mu, 0.0, m) ==
              doctest::Approx(expected).epsilon(1e-12));
    }
}

TEST_CASE("qsys_hh1_lindley mixes over the arrival and service phases") {
    std::vector<double> lam, pa, mu, ps;
    lam.push_back(0.5);
    lam.push_back(2.0);
    pa.push_back(0.4);
    pa.push_back(0.6);
    mu.push_back(1.0);
    mu.push_back(4.0);
    ps.push_back(0.7);
    ps.push_back(0.3);
    // Oracle (a): the docstring reports 1.0484
    const qsys::LindleyResult<double> h = qsys::qsys_hh1_lindley(lam, pa, mu, ps, 1.0);
    CHECK(h.mean[0] == doctest::Approx(1.0484).epsilon(1e-4));
    CHECK(h.var[0] > 0.0);

    // Oracle (b): a degenerate mixture is the M/M/1 step, docstring value 1.8902
    std::vector<double> l1(1, 0.8), p1(1, 1.0), m1(1, 1.0), q1(1, 1.0);
    const qsys::LindleyResult<double> hd = qsys::qsys_hh1_lindley(l1, p1, m1, q1, 2.0);
    CHECK(hd.mean[0] == doctest::Approx(1.8902).epsilon(1e-4));
    const qsys::LindleyResult<double> ref = qsys::qsys_mm1_lindley(0.8, 1.0, 2.0);
    CHECK(hd.mean[0] == doctest::Approx(ref.moments(0, 0)).epsilon(1e-12));
    CHECK(hd.var[0] == doctest::Approx(ref.var[0]).epsilon(1e-12));

    // the mixed variance EXCEEDS the mixture of the per-phase variances, because
    // the phase is itself random and adds the between-phase spread of the means
    double mixvar = 0.0;
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            mixvar += pa[i] * ps[j] * qsys::qsys_mm1_lindley(lam[i], mu[j], 1.0).var[0];
    CHECK(h.var[0] > mixvar);
}

TEST_CASE("qsys_tandem_lindley replays the sample path exactly") {
    std::vector<double> A;
    A.push_back(1.0);
    A.push_back(0.4);
    A.push_back(2.2);
    A.push_back(0.7);
    A.push_back(1.5);
    Matrix<double> S(5, 2, 0.0);
    const double sv[5][2] = {{0.6, 0.9}, {1.1, 0.3}, {0.2, 1.4}, {0.8, 0.5}, {0.35, 0.65}};
    for (std::size_t n = 0; n < 5; ++n)
        for (std::size_t k = 0; k < 2; ++k) S(n, k) = sv[n][k];

    const qsys::TandemLindleyResult<double> r = qsys::qsys_tandem_lindley(A, S);
    REQUIRE(r.W.rows() == 5);
    REQUIRE(r.W.cols() == 2);
    // the first customer starts from W0 = 0 at both stations
    CHECK(r.W(0, 0) == 0.0);
    CHECK(r.W(0, 1) == 0.0);
    // Oracle (c): station 1 is the plain Lindley recursion on A and S(:,1)
    double w = 0.0;
    for (std::size_t n = 0; n + 1 < 5; ++n) {
        w = std::max(w + sv[n][0] - A[n], 0.0);
        CHECK(r.W(n + 1, 0) == doctest::Approx(w).epsilon(1e-14));
    }
    // the gap at station 1 is the interarrival time itself
    for (std::size_t n = 0; n + 1 < 5; ++n) CHECK(r.G(n, 0) == doctest::Approx(A[n]));
    // the last customer has no successor, so its gap row stays unset
    CHECK(std::isnan(r.G(4, 0)));
    CHECK(std::isnan(r.G(4, 1)));
    // Oracle (c): departures accumulate the sojourn times along the tandem
    double epoch = 0.0;
    for (std::size_t n = 0; n < 5; ++n) {
        if (n > 0) epoch += A[n - 1];
        CHECK(r.departure(n, 0) == doctest::Approx(epoch + r.T_(n, 0)).epsilon(1e-14));
        CHECK(r.departure(n, 1) ==
              doctest::Approx(r.departure(n, 0) + r.T_(n, 1)).epsilon(1e-14));
    }
    // Oracle (c): a tandem whose gaps dominate every service never waits
    std::vector<double> Aslow(4, 100.0);
    Matrix<double> Sfast(4, 3, 0.5);
    const qsys::TandemLindleyResult<double> q = qsys::qsys_tandem_lindley(Aslow, Sfast);
    for (std::size_t n = 0; n < 4; ++n)
        for (std::size_t k = 0; k < 3; ++k) CHECK(q.W(n, k) == 0.0);
}

TEST_CASE("qsys_mm1_tandem_lindley couples the two stations through the departures") {
    std::vector<double> Wk, Wk1;
    Wk.push_back(0.0);
    Wk.push_back(2.5);
    Wk1.push_back(1.0);
    Wk1.push_back(0.25);
    const qsys::Mm1TandemLindleyResult<double> r =
        qsys::qsys_mm1_tandem_lindley(0.6, 1.1, 0.9, Wk, Wk1);
    REQUIRE(r.mean.size() == 2);
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(r.mean[i] > 0.0);
        CHECK(r.idleProb[i] > 0.0);
        CHECK(r.idleProb[i] < 1.0);
        // the interdeparture mean is 1/mu1 + q/lambda, between 1/mu1 and 1/mu1+1/lambda
        CHECK(r.interdepMean[i] > 1.0 / 1.1);
        CHECK(r.interdepMean[i] < 1.0 / 1.1 + 1.0 / 0.6);
        CHECK(r.interdepMean[i] ==
              doctest::Approx(1.0 / 1.1 + r.idleProb[i] / 0.6).epsilon(1e-14));
    }
    // a longer wait at station 1 makes an idle period less likely
    CHECK(r.idleProb[1] < r.idleProb[0]);
    CHECK(r.idleProb[0] == doctest::Approx(1.1 / (0.6 + 1.1)).epsilon(1e-14));

    // the coincident-rate branch is continuous in lambda: approaching mu1 from
    // either side must meet the Erlang(2,mu1) form the branch evaluates
    const double mu1 = 1.3, mu2 = 0.75, wk = 0.5, wk1 = 0.5;
    const double at = qsys::qsys_mm1_tandem_lindley(mu1, mu1, mu2, wk, wk1).mean[0];
    const double near = qsys::qsys_mm1_tandem_lindley(mu1 * (1.0 + 1e-6), mu1, mu2, wk, wk1).mean[0];
    CHECK(at == doctest::Approx(near).epsilon(1e-5));
}

TEST_CASE("qsys_mm1_ps returns the Mitra-Morrison moments") {
    // Oracle (b): one class reduces to Coffman, Muntz and Trotter
    for (int i = 1; i <= 4; ++i) {
        const double mu = 1.0 + 0.5 * i, rho = 0.15 * i;
        std::vector<double> lam(1, rho * mu), muv(1, mu);
        const qsys::Mm1PsResult<double> r = qsys::qsys_mm1_ps(lam, muv);
        CHECK(r.alpha == doctest::Approx(1.0 - rho).epsilon(1e-14));
        CHECK(r.W[0] == doctest::Approx(1.0 / (mu * (1.0 - rho))).epsilon(1e-14));
        CHECK(r.W2[0] == doctest::Approx(4.0 / (mu * mu * (1.0 - rho) * (1.0 - rho) *
                                                (2.0 - rho)))
                             .epsilon(1e-12));
    }

    std::vector<double> lam, mu;
    lam.push_back(0.2);
    lam.push_back(0.3);
    mu.push_back(1.0);
    mu.push_back(2.5);
    const qsys::Mm1PsResult<double> r = qsys::qsys_mm1_ps(lam, mu);
    CHECK(r.alpha == doctest::Approx(1.0 - 0.2 / 1.0 - 0.3 / 2.5).epsilon(1e-14));
    // the mean is 1/(alpha mu_r), so the slower class waits longer in proportion
    CHECK(r.W[0] / r.W[1] == doctest::Approx(mu[1] / mu[0]).epsilon(1e-14));
    // every second moment exceeds the square of its mean
    for (std::size_t k = 0; k < 2; ++k) CHECK(r.W2[k] > r.W[k] * r.W[k]);
    // an unstable mix is refused rather than returned negative
    std::vector<double> lbad(1, 2.0), mbad(1, 1.0);
    CHECK_THROWS(qsys::qsys_mm1_ps(lbad, mbad));
}

TEST_CASE("fj_mg1_respt_moments closes the ForkTail white-box route") {
    // M/M/1 branch: ES = 1/mu, ES2 = 2/mu^2, ES3 = 6/mu^3, so SCV = 1 and the
    // response time is Exp(mu-lambda) with mean and variance 1/(mu-lambda)^k
    const double lambda = 0.5, mu = 1.0;
    const fj::Mg1ResptMoments<double> r =
        fj::fj_mg1_respt_moments(lambda, 1.0 / mu, 2.0 / (mu * mu), 6.0 / (mu * mu * mu));
    CHECK(r.ET == doctest::Approx(1.0 / (mu - lambda)).epsilon(1e-14));
    CHECK(r.VT == doctest::Approx(1.0 / ((mu - lambda) * (mu - lambda))).epsilon(1e-14));

    // an unstable branch has no response time moments at all
    CHECK_THROWS(fj::fj_mg1_respt_moments(1.5, 1.0, 2.0, 6.0));
    // and neither has one whose third moment diverges
    CHECK_THROWS(fj::fj_mg1_respt_moments(0.5, 1.0, 2.0,
                                          std::numeric_limits<double>::infinity()));

    // a general branch: the mean is Pollaczek-Khinchine plus the service time
    const double ES = 1.5, ES2 = 3.6, ES3 = 14.0, lam = 0.3;
    const fj::Mg1ResptMoments<double> g = fj::fj_mg1_respt_moments(lam, ES, ES2, ES3);
    const double EW = lam * ES2 / (2.0 * (1.0 - lam * ES));
    CHECK(g.ET == doctest::Approx(EW + ES).epsilon(1e-12));
}

TEST_CASE("fj_tail_forktail fits and inverts the generalized exponential") {
    // Oracle (b): psi(1) = -gamma and psi'(1) = pi^2/6
    CHECK(fj::detail::digamma(1.0) == doctest::Approx(-0.5772156649015329).epsilon(1e-13));
    CHECK(fj::detail::trigamma(1.0) == doctest::Approx(1.6449340668482264).epsilon(1e-13));
    // psi(x+1) = psi(x) + 1/x and psi'(x+1) = psi'(x) - 1/x^2, at a point the
    // recurrence has to walk to reach the asymptotic regime
    CHECK(fj::detail::digamma(4.5) ==
          doctest::Approx(fj::detail::digamma(3.5) + 1.0 / 3.5).epsilon(1e-13));
    CHECK(fj::detail::trigamma(4.5) ==
          doctest::Approx(fj::detail::trigamma(3.5) - 1.0 / (3.5 * 3.5)).epsilon(1e-13));
    CHECK(fj::detail::digamma(0.5) == doctest::Approx(-1.9635100260214235).epsilon(1e-12));
    CHECK(fj::detail::trigamma(0.5) == doctest::Approx(4.934802200544679).epsilon(1e-12));

    // Oracle (b): SCV 1 is the exponential, alpha = 1 and beta = E[T] exactly,
    // and the K-fold maximum then inverts to -beta log(1 - p^(1/K))
    const fj::ForkTailResult<double> e1 = fj::fj_tail_forktail(2.0, 4.0, 1.0, 99.0);
    CHECK(e1.alpha[0] == 1.0);
    CHECK(e1.beta[0] == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(e1.xp == doctest::Approx(-2.0 * std::log(1.0 - 0.99)).epsilon(1e-10));

    const fj::ForkTailResult<double> e8 = fj::fj_tail_forktail(2.0, 4.0, 8.0, 99.0);
    CHECK(e8.xp == doctest::Approx(-2.0 * std::log1p(-std::pow(0.99, 1.0 / 8.0))).epsilon(1e-10));
    CHECK(e8.xp > e1.xp);  // more branches, a longer join

    // a percentage and the matching fraction are the same request
    CHECK(fj::fj_tail_forktail(3.0, 20.0, 4.0, 99.0).xp ==
          doctest::Approx(fj::fj_tail_forktail(3.0, 20.0, 4.0, 0.99).xp).epsilon(1e-12));

    // the fit reproduces the mean and variance it was matched on
    const fj::ForkTailResult<double> h = fj::fj_tail_forktail(3.0, 2.0, 5.0, 95.0);
    const double a = h.alpha[0], b = h.beta[0];
    CHECK(b * (fj::detail::digamma(a + 1.0) - fj::detail::digamma(1.0)) ==
          doctest::Approx(3.0).epsilon(1e-9));
    CHECK(b * b * (fj::detail::trigamma(1.0) - fj::detail::trigamma(a + 1.0)) ==
          doctest::Approx(2.0).epsilon(1e-9));
    CHECK(a > 1.0);  // SCV below 1 needs a shape above the exponential

    // heterogeneous branches: the request percentile solves the product of the
    // branch CDFs and must exceed every single-branch percentile
    std::vector<double> ET, VT;
    ET.push_back(2.0);
    ET.push_back(3.0);
    ET.push_back(1.5);
    VT.push_back(4.0);
    VT.push_back(2.0);
    VT.push_back(3.0);
    const fj::ForkTailResult<double> het = fj::fj_tail_forktail(ET, VT, std::vector<double>(),
                                                               99.0);
    REQUIRE(het.alpha.size() == 3);
    double resid = 0.0;
    for (std::size_t i = 0; i < 3; ++i)
        resid += het.alpha[i] * std::log1p(-std::exp(-het.xp / het.beta[i]));
    CHECK(resid == doctest::Approx(std::log(0.99)).epsilon(1e-9));
    for (std::size_t i = 0; i < 3; ++i) {
        const double single = -het.beta[i] * std::log1p(-std::pow(0.99, 1.0 / het.alpha[i]));
        CHECK(het.xp >= single - 1e-9);
    }

    // random fanout: the mixture percentile lies between the two extreme fanouts
    std::vector<double> K, P;
    K.push_back(2.0);
    K.push_back(5.0);
    K.push_back(9.0);
    P.push_back(0.2);
    P.push_back(0.5);
    P.push_back(0.3);
    const fj::ForkTailResult<double> mix =
        fj::fj_tail_forktail(std::vector<double>(1, 2.0), std::vector<double>(1, 4.0), K, 99.0, P);
    CHECK(mix.xp > fj::fj_tail_forktail(2.0, 4.0, 2.0, 99.0).xp);
    CHECK(mix.xp < fj::fj_tail_forktail(2.0, 4.0, 9.0, 99.0).xp);
    double cdf = 0.0;
    for (std::size_t i = 0; i < 3; ++i)
        cdf += P[i] * std::pow(1.0 - std::exp(-mix.xp / mix.beta[0]), K[i] * mix.alpha[0]);
    CHECK(cdf == doctest::Approx(0.99).epsilon(1e-9));

    // a percentile outside (0,1) after the percentage rescale is refused
    CHECK_THROWS(fj::fj_tail_forktail(2.0, 4.0, 1.0, 0.0));
}
