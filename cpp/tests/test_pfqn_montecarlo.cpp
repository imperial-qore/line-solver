/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Monte Carlo normalizing-constant estimators and the perfect sampler.
 *
 * Oracles, in order of strength:
 *  - the EXACT constant from pfqn_ca / pfqn_ncld, against which every estimator
 *    is checked both for accuracy at a stated sample count and for the
 *    1/sqrt(n) decay of its error. The decay is the property that separates a
 *    correct estimator from a biased one: a biased estimator is just as likely
 *    to look accurate at one sample count, and stops improving at the next.
 *    The error is measured as a root-mean-square over eight independent seeds,
 *    because a single seed's error is itself a random variable and its
 *    trajectory in n is not monotone.
 *  - closed forms that make an estimator's variance zero, which turn a Monte
 *    Carlo routine into a deterministic identity (constant OI rank rates; a
 *    swap graph admitting exactly one ordering).
 *  - exact marginal queue lengths from pfqn_mva / MATLAB pfqn_mvald for the
 *    perfect sampler.
 *
 * MATLAB comparability. These estimators are comparable to MATLAB only IN
 * DISTRIBUTION: the generators, the uniform-to-integer maps and the normal
 * transforms all differ, so no fixed seed reproduces a MATLAB sample path.
 * What is compared against MATLAB here is the ESTIMAND and the CONVERGENCE
 * RATE. For pfqn_ls, whose Z > 0 branch converges markedly more slowly than
 * 1/sqrt(n) would suggest at small n, the MATLAB rate was measured on the same
 * model (relative error 1.8e-1 at I = 200 falling to 2.0e-2 at I = 200000) and
 * the port reproduces it (2.5e-1 falling to 1.9e-2), so the slow rate is the
 * reference's and not a porting defect.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_cftp.h"
#include "line/api/pfqn/pfqn_is.h"
#include "line/api/pfqn/pfqn_ld_is.h"
#include "line/api/pfqn/pfqn_ls.h"
#include "line/api/pfqn/pfqn_mci.h"
#include "line/api/pfqn/pfqn_mmsample2.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_ncld.h"
#include "line/api/pfqn/pfqn_oi_is.h"
#include "line/api/pfqn/pfqn_pas_is.h"

using line::Matrix;
using namespace line::pfqn;

namespace {

/** Two single-server queues, two classes; the model every estimator is run on. */
Matrix<double> mcModel() {
    Matrix<double> L(2, 2);
    L(0, 0) = 0.5; L(0, 1) = 0.3;
    L(1, 0) = 0.2; L(1, 1) = 0.4;
    return L;
}

const std::vector<int> mcN{3, 2};
const std::vector<double> mcZ{1.0, 1.0};

/**
 * Root-mean-square relative error of an estimator over eight independent
 * seeds. One seed is not enough to see the 1/sqrt(n) law: its error is itself
 * random and moves non-monotonically in n.
 */
template <class F>
double rmsRelErr(F estimate, double exact, std::size_t n) {
    double acc = 0.0;
    for (unsigned k = 1; k <= 8; ++k) {
        McRng g(k * 7919u);
        const double e = (estimate(n, g) - exact) / exact;
        acc += e * e;
    }
    return std::sqrt(acc / 8.0);
}

/** Two-server capacities at station 1, single server at station 2. */
Matrix<double> mcCapacities() {
    Matrix<double> mu(2, 5);
    for (std::size_t k = 0; k < 5; ++k) {
        mu(0, k) = k + 1 < 2 ? 1.0 : 2.0;
        mu(1, k) = 1.0;
    }
    return mu;
}

}  // namespace

// ---------------------------------------------------------------------------
// pfqn_is / pfqn_ld_is
// ---------------------------------------------------------------------------

TEST_CASE("importance sampling converges to the exact constant at the 1/sqrt(n) rate") {
    const double exact = pfqn_ca(mcModel(), mcN, Matrix<double>::row({1.0, 1.0})).G;
    CHECK(exact == doctest::Approx(2.98312666666667).epsilon(1e-12));  // MATLAB pfqn_ca

    auto est = [](std::size_t n, McRng& g) { return pfqn_is(mcModel(), mcN, mcZ, n, g).G; };
    const double r100 = rmsRelErr(est, exact, 100);
    const double r10k = rmsRelErr(est, exact, 10000);
    // Observed: 6.93e-2 at n = 100, 2.69e-3 at n = 10^4, a factor of 25.8 for a
    // hundredfold sample count (the law predicts 10).
    CHECK(r100 < 1.0e-1);
    CHECK(r10k < 5.0e-3);
    CHECK(r10k < r100 / 4.0);
}

TEST_CASE("pfqn_is is exactly pfqn_ld_is with unit capacities") {
    McRng g1(4242), g2(4242);
    const double a = pfqn_is(mcModel(), mcN, mcZ, std::size_t(500), g1).G;
    const double b = pfqn_ld_is(mcModel(), mcN, mcZ, Matrix<double>(), std::size_t(500), g2).G;
    CHECK(a == b);  // same code path, same generator state: bit for bit
}

TEST_CASE("load-dependent importance sampling converges to the exact load-dependent constant") {
    const Matrix<double> mu = mcCapacities();
    // THE ORACLE MUST NAME 'exact'. On this model `default` clears the
    // Choudhury-Leung-Whitt cost gate (2 stations, 2 classes, N summing to 5)
    // and answers by contour inversion, which is accurate to ~1e-10 and not to
    // the last bit -- MATLAB's pfqn_ncld does the same and reports method
    // 'clw'. The exact ladder returns 824399/480000 = 1.7174979166666666, the
    // brute-force convolution of the load-dependent product form.
    const double exact = pfqn_ncld(mcModel(), mcN, Matrix<double>::row({1.0, 1.0}), mu,
                                   NcldMethod::Exact, 0.0, NcOptions())
                             .G;
    CHECK(exact == doctest::Approx(1.7174979166666666).epsilon(1e-12));
    // The default route is the same constant to the accuracy CLW gives it.
    CHECK(pfqn_ncld(mcModel(), mcN, Matrix<double>::row({1.0, 1.0}), mu).G ==
          doctest::Approx(exact).epsilon(1e-9));

    auto est = [&](std::size_t n, McRng& g) {
        return pfqn_ld_is(mcModel(), mcN, mcZ, mu, n, g).G;
    };
    const double r100 = rmsRelErr(est, exact, 100);
    const double r10k = rmsRelErr(est, exact, 10000);
    // Observed: 6.45e-2 -> 2.15e-3, a factor of 30.
    CHECK(r100 < 1.0e-1);
    CHECK(r10k < 5.0e-3);
    CHECK(r10k < r100 / 4.0);
}

TEST_CASE("an empty closed population makes the importance sampler return one exactly") {
    McRng g(1);
    const auto r = pfqn_is(mcModel(), std::vector<int>{0, 0}, mcZ, std::size_t(10), g);
    CHECK(r.G == 1.0);
    CHECK(r.lG == 0.0);
}

// ---------------------------------------------------------------------------
// pfqn_mci
// ---------------------------------------------------------------------------

TEST_CASE("Monte Carlo integration converges to the exact constant at the 1/sqrt(n) rate") {
    const double exact = pfqn_ca(mcModel(), mcN, Matrix<double>::row({1.0, 1.0})).G;

    auto imci = [](std::size_t n, McRng& g) {
        return pfqn_mci(mcModel(), mcN, mcZ, n, MciVariant::Imci, g).G;
    };
    const double r100 = rmsRelErr(imci, exact, 100);
    const double r10k = rmsRelErr(imci, exact, 10000);
    // Observed: 5.17e-2 -> 4.73e-3, a factor of 10.9, the textbook rate.
    CHECK(r100 < 1.0e-1);
    CHECK(r10k < 1.0e-2);
    CHECK(r10k < r100 / 4.0);
}

TEST_CASE("the two MCI proposal rules agree while every station is below the saturation clamp") {
    // On this model the Bard-Schweitzer utilizations are both under 0.9 and
    // 1 - U exceeds 0.01, so gamma = 1 - U under both rules and the estimates
    // are identical. The rules can only differ once a station saturates.
    McRng g1(31337), g2(31337);
    const double a = pfqn_mci(mcModel(), mcN, mcZ, std::size_t(2000), MciVariant::Imci, g1).G;
    const double b = pfqn_mci(mcModel(), mcN, mcZ, std::size_t(2000), MciVariant::Mci, g2).G;
    CHECK(a == b);
}

TEST_CASE("the repairman MCI variant is accepted for one station and refused for several") {
    Matrix<double> one(1, 2);
    one(0, 0) = 0.5;
    one(0, 1) = 0.3;
    const double exact = pfqn_ca(one, mcN, Matrix<double>::row({1.0, 1.0})).G;
    auto rm = [&](std::size_t n, McRng& g) {
        return pfqn_mci(one, mcN, mcZ, n, MciVariant::Rm, g).G;
    };
    const double r100 = rmsRelErr(rm, exact, 100);
    const double r10k = rmsRelErr(rm, exact, 10000);
    CHECK(r10k < r100 / 4.0);

    // For M > 1 the reference's rate expression is dimensionally inconsistent
    // (a (1 x R) row broadcast against an (M x R) matrix), so nothing faithful
    // can be computed and the port refuses rather than inventing a rule.
    McRng g(1);
    CHECK_THROWS_AS(pfqn_mci(mcModel(), mcN, mcZ, std::size_t(10), MciVariant::Rm, g),
                    line::UnsupportedError);
}

TEST_CASE("a demand-free model is answered by the delay term alone") {
    Matrix<double> empty(2, 2, 0.0);
    McRng g(1);
    const auto r = pfqn_mci(empty, mcN, mcZ, std::size_t(10), MciVariant::Imci, g);
    // G = prod_r Z_r^{N_r} / N_r! = 1/(3! 2!) = 1/12.
    CHECK(r.G == doctest::Approx(1.0 / 12.0).epsilon(1e-12));
}

// ---------------------------------------------------------------------------
// pfqn_ls
// ---------------------------------------------------------------------------

TEST_CASE("logistic sampling converges to the exact constant, Z = 0 branch") {
    const std::vector<double> Nd{3.0, 2.0};
    const double exact = pfqn_ca(mcModel(), mcN, Matrix<double>()).G;
    CHECK(exact == doctest::Approx(0.33586).epsilon(1e-10));  // MATLAB pfqn_ca

    auto est = [&](std::size_t n, McRng& g) {
        return pfqn_ls(mcModel(), Nd, std::vector<double>(), n, g).G;
    };
    const double r100 = rmsRelErr(est, exact, 100);
    const double r10k = rmsRelErr(est, exact, 10000);
    // Observed: 5.19e-1 -> 8.89e-3, a factor of 58. The proposal is a Gaussian
    // fitted at the mode, so the weights are heavy tailed at very small n.
    CHECK(r100 < 1.0);
    CHECK(r10k < 3.0e-2);
    CHECK(r10k < r100 / 4.0);
}

TEST_CASE("logistic sampling converges to the exact constant, Z > 0 branch") {
    const std::vector<double> Nd{3.0, 2.0};
    const double exact = pfqn_ca(mcModel(), mcN, Matrix<double>::row({1.0, 1.0})).G;

    auto est = [&](std::size_t n, McRng& g) { return pfqn_ls(mcModel(), Nd, mcZ, n, g).G; };
    const double r100 = rmsRelErr(est, exact, 100);
    const double r10k = rmsRelErr(est, exact, 10000);
    // Observed: 8.85e-1 -> 3.14e-2. MATLAB on the same model goes from 1.76e-1
    // at I = 200 to 2.03e-2 at I = 200000, i.e. the same slow rate, so this is
    // the reference's behaviour and not a porting defect.
    CHECK(r100 < 1.5);
    CHECK(r10k < 1.0e-1);
    CHECK(r10k < r100 / 4.0);
}

TEST_CASE("logistic sampling refuses a model with fewer than two loaded stations") {
    Matrix<double> one(1, 2);
    one(0, 0) = 0.5;
    one(0, 1) = 0.3;
    McRng g(1);
    CHECK_THROWS_AS(pfqn_ls(one, std::vector<double>{3.0, 2.0}, mcZ, std::size_t(10), g),
                    line::InputError);
}

// ---------------------------------------------------------------------------
// pfqn_oi_is / pfqn_pas_is
// ---------------------------------------------------------------------------

TEST_CASE("OI importance sampling is exact when the estimator has no variance") {
    // Two ingredients make the per-sample value deterministic here. Constant
    // rank rates make every ordering carry the same weight S(c), and a
    // population with one job per class makes the branching factors R, R-1,
    // ..., 1 the same along every ordering, so 1/p(c) is the same too. The
    // estimator then returns its mean at any sample count, and that mean is
    //   G = (ell! / prod_r N_r!) sum_{k=0}^{ell} mu_1^{-k} mu_2^{-(ell-k)}.
    // Note that this needs BOTH: with N = [2,1] the branching factors differ
    // between orderings and the estimator regains a variance, which the next
    // case exercises.
    const std::vector<int> N{1, 1};
    const double m1 = 2.0, m2 = 3.0;
    std::vector<OiRateFun<double>> mu{[&](const std::vector<int>&) { return m1; },
                                      [&](const std::vector<int>&) { return m2; }};
    double S = 0.0, Sk = 0.0;
    for (long k = 0; k <= 2; ++k) {
        const double w = std::pow(m1, -static_cast<double>(k)) *
                         std::pow(m2, -static_cast<double>(2 - k));
        S += w;
        Sk += k * w;
    }
    const double exact = 2.0 * S;  // ell! / (1! 1!) = 2 orderings
    for (unsigned seed : {11u, 12u, 13u}) {
        McRng g(seed);
        const auto r = pfqn_oi_is(N, mu, std::size_t(37), g);
        CHECK(r.G == doctest::Approx(exact).epsilon(1e-13));
        // The TOTAL queue length at station 1 is exact for the same reason: it
        // is the cut position, which does not depend on the ordering. Its SPLIT
        // between the classes does, since it is the prefix's class composition,
        // so the per-class values still carry sampling error and only converge.
        CHECK(r.Q(0, 0) + r.Q(0, 1) == doctest::Approx(Sk / S).epsilon(1e-13));
        CHECK(r.Q(0, 0) == doctest::Approx(Sk / S / 2.0).epsilon(5e-2));
        CHECK(r.Q(0, 1) == doctest::Approx(Sk / S / 2.0).epsilon(5e-2));
        for (std::size_t c = 0; c < 2; ++c)
            CHECK(r.Q(0, c) + r.Q(1, c) == doctest::Approx(N[c]).epsilon(1e-12));
    }
    // With enough samples the split settles on the symmetric answer.
    McRng gbig(11);
    const auto rbig = pfqn_oi_is(N, mu, std::size_t(200000), gbig);
    CHECK(rbig.Q(0, 0) == doctest::Approx(Sk / S / 2.0).epsilon(1e-2));
    CHECK(rbig.Q(0, 1) == doctest::Approx(Sk / S / 2.0).epsilon(1e-2));
}

TEST_CASE("OI importance sampling stays unbiased when the branching factors vary") {
    // N = [2,1] with constant rates: the weight S(c) is the same for every
    // ordering but 1/p(c) is 4 for the two orderings that start with class 0
    // and 2 for the one that starts with class 1, so the estimator has a
    // genuine variance and must be checked by convergence, not by identity.
    const std::vector<int> N{2, 1};
    const double m1 = 2.0, m2 = 3.0;
    std::vector<OiRateFun<double>> mu{[&](const std::vector<int>&) { return m1; },
                                      [&](const std::vector<int>&) { return m2; }};
    double S = 0.0;
    for (long k = 0; k <= 3; ++k)
        S += std::pow(m1, -static_cast<double>(k)) * std::pow(m2, -static_cast<double>(3 - k));
    const double exact = 3.0 * S;  // 3! / (2! 1!) = 3 orderings
    CHECK(exact == doctest::Approx(0.902777777777778).epsilon(1e-12));

    auto est = [&](std::size_t n, McRng& g) { return pfqn_oi_is(N, mu, n, g).G; };
    const double r100 = rmsRelErr(est, exact, 100);
    const double r10k = rmsRelErr(est, exact, 10000);
    CHECK(r10k < r100 / 4.0);
    CHECK(r10k < 5.0e-3);
}

TEST_CASE("OI importance sampling converges when the rank rates do depend on the support") {
    const std::vector<int> N{2, 1};
    auto r1 = [](const std::vector<int>& s) { return 1.0 + (s[0] ? 1.0 : 0.0) + (s[1] ? 0.5 : 0.0); };
    auto r2 = [](const std::vector<int>& s) { return 2.0 + (s[1] ? 1.0 : 0.0); };
    std::vector<OiRateFun<double>> mu{r1, r2};
    // Brute force over the three orderings of the multiset {0,0,1}.
    const std::vector<std::vector<int>> ords{{0, 0, 1}, {0, 1, 0}, {1, 0, 0}};
    double G = 0.0, Q0 = 0.0;
    for (const auto& c : ords) {
        std::vector<double> P1(4, 1.0), P2(4, 1.0);
        std::vector<int> sp(2, 0), occ(2, 0);
        std::vector<std::vector<int>> cnt(4, std::vector<int>(2, 0));
        double ph = 1.0;
        for (int k = 0; k < 3; ++k) {
            sp[c[k]] = 1;
            occ[c[k]] += 1;
            ph /= r1(sp);
            P1[k + 1] = ph;
            cnt[k + 1] = occ;
        }
        sp.assign(2, 0);
        ph = 1.0;
        for (int k = 3; k >= 1; --k) {
            sp[c[k - 1]] = 1;
            ph /= r2(sp);
            P2[k - 1] = ph;
        }
        for (int k = 0; k <= 3; ++k) {
            const double w = P1[k] * (k >= 3 ? 1.0 : P2[k]);
            G += w;
            if (k > 0) Q0 += w * cnt[k][0];
        }
    }
    CHECK(G == doctest::Approx(1.08481481481).epsilon(1e-10));

    auto estG = [&](std::size_t n, McRng& g) { return pfqn_oi_is(N, mu, n, g).G; };
    auto estQ = [&](std::size_t n, McRng& g) { return pfqn_oi_is(N, mu, n, g).Q(0, 0); };
    const double g100 = rmsRelErr(estG, G, 100);
    const double g10k = rmsRelErr(estG, G, 10000);
    const double q100 = rmsRelErr(estQ, Q0 / G, 100);
    const double q10k = rmsRelErr(estQ, Q0 / G, 10000);
    // Observed: G 8.07e-3 -> 7.09e-4; Q 4.86e-2 -> 3.32e-3.
    CHECK(g10k < g100 / 4.0);
    CHECK(q10k < q100 / 4.0);
    CHECK(g10k < 3.0e-3);
    CHECK(q10k < 1.0e-2);
}

TEST_CASE("the placement closure makes the swap graph transitive") {
    Matrix<int> H(3, 3, 0);
    H(0, 1) = 1;
    H(1, 2) = 1;
    const Matrix<int> P = pas_placement(H);
    CHECK(P(0, 1) == 1);
    CHECK(P(1, 2) == 1);
    CHECK(P(0, 2) == 1);  // the closure, not present in H
    CHECK(P(2, 0) == 0);
    CHECK(pas_placement(Matrix<int>()).empty());
}

TEST_CASE("pass-and-swap with a total-order swap graph is a deterministic identity") {
    // H forces class 0 before class 1, so with N = [1,1] exactly ONE ordering
    // is feasible and every sample carries the same value. With constant rates
    // mu_1 = 2, mu_2 = 4 the ordering (0,1) gives
    //   G = sum_{k=0}^{2} 2^{-k} 4^{-(2-k)} = 1/16 + 1/8 + 1/4 = 7/16,
    // and the prefix counts give Q_1 = (0.125 + 0.25)/G, 0.25/G.
    const std::vector<int> N{1, 1};
    Matrix<int> H(2, 2, 0);
    H(0, 1) = 1;
    std::vector<OiRateFun<double>> mu{[](const std::vector<int>&) { return 2.0; },
                                      [](const std::vector<int>&) { return 4.0; }};
    McRng g(5);
    const auto r = pfqn_pas_is(N, mu, H, std::size_t(50), g);
    CHECK(r.G == doctest::Approx(0.4375).epsilon(1e-14));
    CHECK(r.Q(0, 0) == doctest::Approx(0.375 / 0.4375).epsilon(1e-14));
    CHECK(r.Q(0, 1) == doctest::Approx(0.25 / 0.4375).epsilon(1e-14));
    CHECK(r.Q(1, 0) == doctest::Approx(1.0 - 0.375 / 0.4375).epsilon(1e-12));
}

TEST_CASE("pass-and-swap with an empty swap graph is the OI estimator") {
    const std::vector<int> N{2, 1};
    std::vector<OiRateFun<double>> mu{[](const std::vector<int>& s) { return 1.0 + s[0]; },
                                      [](const std::vector<int>& s) { return 2.0 + s[1]; }};
    McRng g1(808), g2(808);
    const auto a = pfqn_pas_is(N, mu, Matrix<int>(), std::size_t(300), g1);
    const auto b = pfqn_oi_is(N, mu, std::size_t(300), g2);
    CHECK(a.G == b.G);
    CHECK(a.Q(0, 0) == b.Q(0, 0));
}

TEST_CASE("a cyclic swap graph admits no ordering and is reported, not worked around") {
    const std::vector<int> N{1, 1};
    Matrix<int> H(2, 2, 0);
    H(0, 1) = 1;
    H(1, 0) = 1;
    std::vector<OiRateFun<double>> mu{[](const std::vector<int>&) { return 1.0; },
                                      [](const std::vector<int>&) { return 1.0; }};
    McRng g(1);
    CHECK_THROWS_AS(pfqn_pas_is(N, mu, H, std::size_t(1), g), line::NumericError);
}

// ---------------------------------------------------------------------------
// pfqn_cftp
// ---------------------------------------------------------------------------

TEST_CASE("coupling from the past samples the exact product-form stationary distribution") {
    const std::vector<double> L{0.5, 0.4, 0.3};
    const int K = 6;
    // Exact single-class marginal means from pfqn_mva.
    Matrix<double> Lc(3, 1);
    Lc(0, 0) = 0.5; Lc(1, 0) = 0.4; Lc(2, 0) = 0.3;
    const auto exact = pfqn_mva(Lc, std::vector<int>{K});
    CHECK(exact.QN(0, 0) == doctest::Approx(3.13500980782672).epsilon(1e-12));  // MATLAB
    CHECK(exact.QN(1, 0) == doctest::Approx(1.83259805947827).epsilon(1e-12));
    CHECK(exact.QN(2, 0) == doctest::Approx(1.03239213269501).epsilon(1e-12));

    McRng g(7);
    const auto r = pfqn_cftp(L, K, std::vector<int>(), std::size_t(50000), CftpMethod::Cftp, g);
    for (std::size_t i = 0; i < 3; ++i)
        CHECK(r.Q[i] == doctest::Approx(exact.QN(i, 0)).epsilon(0.02));

    // Every draw is a feasible state, and the coalescence horizon is a power of
    // two because the search doubles it.
    for (std::size_t s = 0; s < r.X.rows(); ++s) {
        int tot = 0;
        for (std::size_t i = 0; i < 3; ++i) {
            CHECK(r.X(s, i) >= 0);
            tot += r.X(s, i);
        }
        CHECK(tot == K);
        const long h = r.horizon[s];
        CHECK(h >= 1);
        CHECK((h & (h - 1)) == 0);
    }
}

TEST_CASE("the sampling error of coupling from the past shrinks at the 1/sqrt(n) rate") {
    const std::vector<double> L{0.5, 0.4, 0.3};
    Matrix<double> Lc(3, 1);
    Lc(0, 0) = 0.5; Lc(1, 0) = 0.4; Lc(2, 0) = 0.3;
    const double exact = pfqn_mva(Lc, std::vector<int>{6}).QN(0, 0);
    auto est = [&](std::size_t n, McRng& g) {
        return pfqn_cftp(L, 6, std::vector<int>(), n, CftpMethod::Cftp, g).Q[0];
    };
    const double r100 = rmsRelErr(est, exact, 100);
    const double r10k = rmsRelErr(est, exact, 10000);
    CHECK(r10k < r100 / 4.0);
}

TEST_CASE("coupling from the past reproduces the exact multiserver marginals") {
    const std::vector<double> L{0.5, 0.4, 0.3};
    // MATLAB pfqn_mvald with mu(i,k) = min(k, S_i), S = [2 1 1], K = 6.
    McRng g(2025);
    const auto r = pfqn_cftp(L, 6, std::vector<int>{2, 1, 1}, std::size_t(50000), CftpMethod::Cftp, g);
    CHECK(r.Q[0] == doctest::Approx(1.45256475682458).epsilon(0.03));
    CHECK(r.Q[1] == doctest::Approx(2.99423083638476).epsilon(0.03));
    CHECK(r.Q[2] == doctest::Approx(1.55320440679067).epsilon(0.03));

    // An infinite-server station, mu(1,k) = k.
    McRng g2(2026);
    const auto ri =
        pfqn_cftp(L, 6, std::vector<int>{cftp_inf_servers, 1, 1}, std::size_t(50000),
                  CftpMethod::Cftp, g2);
    CHECK(ri.Q[0] == doctest::Approx(1.14439318386371).epsilon(0.03));
    CHECK(ri.Q[1] == doctest::Approx(3.20975319042934).epsilon(0.03));
    CHECK(ri.Q[2] == doctest::Approx(1.64585362570694).epsilon(0.03));
}

TEST_CASE("the approximate sampler targets the same distribution as the exact one") {
    const std::vector<double> L{0.5, 0.4, 0.3};
    McRng g(2027);
    const auto r = pfqn_cftp(L, 6, std::vector<int>(), std::size_t(50000), CftpMethod::Approx, g);
    CHECK(r.Q[0] == doctest::Approx(3.13500980782672).epsilon(0.02));
    CHECK(r.Q[1] == doctest::Approx(1.83259805947827).epsilon(0.02));
    CHECK(r.Q[2] == doctest::Approx(1.03239213269501).epsilon(0.02));
    // Its step count is deterministic, unlike the coalescence horizon.
    for (long h : r.horizon) CHECK(h == r.horizon[0]);
}

TEST_CASE("coupling from the past rejects models it cannot sample") {
    McRng g(1);
    CHECK_THROWS_AS(pfqn_cftp(std::vector<double>{0.5}, 4, std::vector<int>(), std::size_t(1),
                              CftpMethod::Cftp, g),
                    line::InputError);
    CHECK_THROWS_AS(pfqn_cftp(std::vector<double>{0.5, 0.0}, 4, std::vector<int>(), std::size_t(1),
                              CftpMethod::Cftp, g),
                    line::InputError);
}

// ---------------------------------------------------------------------------
// pfqn_mmsample2
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_mmsample2 reproduces the reference's max-of-grid form, defect included") {
    // The exact repairman constant for L = [0.5 0.3], N = [3 2], Z = [1 1] is
    // (sum N)!/prod N_r! ... no: with a single queue the McKenna-Mitra integral
    // gives lG = log((sum N)!) + sum_r N_r log(Z_r + L_r) - sum_r log N_r!,
    // which pfqn_ca confirms as lG = 0.26300508473934.
    Matrix<double> L(1, 2);
    L(0, 0) = 0.5;
    L(0, 1) = 0.3;
    const double exactLG = pfqn_ca(L, mcN, Matrix<double>::row({1.0, 1.0})).lG;
    CHECK(exactLG == doctest::Approx(0.26300508473934).epsilon(1e-11));

    McRng g1(8), g2(8);
    const double lg1 = pfqn_mmsample2(L, mcN, mcZ, std::size_t(10000), g1).lG;
    const double lg2 = pfqn_mmsample2(L, mcN, mcZ, std::size_t(1000000), g2).lG;

    // The reference takes the MAXIMUM of du + f(v) over an UNSORTED grid rather
    // than summing the quadrature contributions, so du (a sign-indefinite
    // difference of consecutive grid points, the largest of which is of order
    // 1e5) is added in the LOG domain. The estimate therefore overshoots the
    // exact constant by a wide margin and, being an argmax over a grid whose
    // extent does not change with the sample count, does NOT converge: raising
    // the grid a hundredfold leaves the error essentially where it was. Both
    // properties are asserted here because they are what the reference does;
    // MATLAB on the same model returns lG = 2.3149216118351 against the same
    // exact 0.26300508473934.
    CHECK(lg1 > exactLG);
    CHECK(lg2 > exactLG);
    const double e1 = std::fabs(lg1 - exactLG);
    const double e2 = std::fabs(lg2 - exactLG);
    CHECK(e1 > 1.0);
    CHECK(e2 > e1 / 2.0);  // no 1/sqrt(n) decay: the error does not shrink
}

// ---------------------------------------------------------------------------
// generator contract
// ---------------------------------------------------------------------------

TEST_CASE("every estimator is reproducible from the generator state and from nothing else") {
    McRng a(2718), b(2718), c(3141);
    const double x1 = pfqn_is(mcModel(), mcN, mcZ, std::size_t(200), a).G;
    const double x2 = pfqn_is(mcModel(), mcN, mcZ, std::size_t(200), b).G;
    const double x3 = pfqn_is(mcModel(), mcN, mcZ, std::size_t(200), c).G;
    CHECK(x1 == x2);
    CHECK(x1 != x3);
    // The generator is advanced by the call, so a second call on the same
    // object does not repeat the first.
    const double x4 = pfqn_is(mcModel(), mcN, mcZ, std::size_t(200), a).G;
    CHECK(x1 != x4);
}
