/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * First passage times in Markov and semi-Markov chains.
 *
 * THE ORACLES SHARE NOTHING WITH THE FORMULAS. The CDF is checked against a
 * directly assembled 1 - alpha expm(St) 1 built in the test; the moments against
 * the dense closed form n! alpha (-S)^{-n} 1, which the implementation
 * deliberately does NOT use; the semi-Markov recursions against the Markov ones
 * on a chain that is both. The deterministic and Erlang holding-time cases have
 * exact analytic answers.
 *
 * Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time Distributions
 * in Large Markov Chains", 2002. Cross-checked against the MATLAB and native
 * Python twins, which agree with these numbers.
 */

#include <cmath>
#include <complex>
#include <vector>

#include "doctest.h"
#include "line/api/mc/ctmc_passage.h"
#include "line/util/expm.h"

namespace {

using line::Matrix;
namespace mc = line::mc;

/** M/M/1/K generator, K+1 states, arrival lambda and service mu. */
Matrix<double> mm1k(std::size_t K, double lambda, double mu) {
    const std::size_t n = K + 1;
    Matrix<double> Q(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i + 1 < n) Q(i, i + 1) = lambda;
        if (i > 0) Q(i, i - 1) = mu;
    }
    for (std::size_t i = 0; i < n; ++i) {
        double r = 0.0;
        for (std::size_t j = 0; j < n; ++j)
            if (j != i) r += Q(i, j);
        Q(i, i) = -r;
    }
    return Q;
}

}  // namespace

TEST_CASE("the passage law is the phase-type absorption law, built independently") {
    const std::size_t K = 6, n = K + 1;
    const Matrix<double> Q = mm1k(K, 1.0, 1.5);
    std::vector<double> pi0(n, 0.0);
    pi0[0] = 1.0;
    const std::vector<std::size_t> target{n - 1};

    // The oracle: assemble S by hand and exponentiate at each point.
    Matrix<double> S(n - 1, n - 1);
    for (std::size_t i = 0; i + 1 < n; ++i)
        for (std::size_t j = 0; j + 1 < n; ++j) S(i, j) = Q(i, j);
    std::vector<double> s0(n - 1, 0.0);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        double r = 0.0;
        for (std::size_t j = 0; j + 1 < n; ++j) r += S(i, j);
        s0[i] = -r;
    }

    std::vector<double> tset;
    for (double t = 0.0; t <= 20.0 + 1e-9; t += 0.25) tset.push_back(t);
    const mc::PassageCurve<double> c = mc::ctmc_passage_time(Q, pi0, target, tset);

    CHECK(c.atom == doctest::Approx(0.0));
    for (std::size_t it = 0; it < tset.size(); ++it) {
        Matrix<double> St = S;
        for (std::size_t i = 0; i + 1 < n; ++i)
            for (std::size_t j = 0; j + 1 < n; ++j) St(i, j) = S(i, j) * tset[it];
        const Matrix<double> E = line::expm(St);
        double sF = 0.0, sf = 0.0;
        for (std::size_t j = 0; j + 1 < n; ++j) {
            double acc = E(0, j);  // alpha is the first unit vector here
            sF += acc;
            sf += acc * s0[j];
        }
        CHECK(c.F[it] == doctest::Approx(1.0 - sF).epsilon(1e-12));
        CHECK(c.f[it] == doctest::Approx(sf).epsilon(1e-12));
    }
}

TEST_CASE("Eq. 3 reproduces n! alpha (-S)^-n 1, which it deliberately does not use") {
    const std::size_t K = 6, n = K + 1;
    const Matrix<double> Q = mm1k(K, 1.0, 1.5);
    std::vector<double> pi0(n, 0.0);
    pi0[0] = 1.0;
    const std::vector<std::size_t> target{n - 1};
    const mc::PassageMoments<double> pm = mc::ctmc_passage_moments(Q, pi0, target, 4);

    // Oracle: the dense inverse, formed once and powered.
    const std::size_t nA = n - 1;
    Matrix<double> A(nA, nA);
    for (std::size_t i = 0; i < nA; ++i)
        for (std::size_t j = 0; j < nA; ++j) A(i, j) = -Q(i, j);
    Matrix<double> Ainv(nA, nA);
    for (std::size_t c = 0; c < nA; ++c) {
        std::vector<double> e(nA, 0.0);
        e[c] = 1.0;
        const std::vector<double> col = line::solve(A, e);
        for (std::size_t r = 0; r < nA; ++r) Ainv(r, c) = col[r];
    }
    Matrix<double> P(nA, nA);
    for (std::size_t i = 0; i < nA; ++i) P(i, i) = 1.0;
    double fact = 1.0;
    for (std::size_t k = 1; k <= 4; ++k) {
        Matrix<double> Pn(nA, nA);
        for (std::size_t i = 0; i < nA; ++i)
            for (std::size_t j = 0; j < nA; ++j) {
                double acc = 0.0;
                for (std::size_t q = 0; q < nA; ++q) acc += P(i, q) * Ainv(q, j);
                Pn(i, j) = acc;
            }
        P = Pn;
        fact *= double(k);
        double m = 0.0;
        for (std::size_t j = 0; j < nA; ++j) m += P(0, j);
        CHECK(pm.m[k - 1] == doctest::Approx(fact * m).epsilon(1e-11));
    }
    // The value the MATLAB and Python twins report for the mean.
    CHECK(pm.m[0] == doctest::Approx(50.34375).epsilon(1e-12));
}

TEST_CASE("the atom at zero is reported, not dropped") {
    const std::size_t K = 6, n = K + 1;
    const Matrix<double> Q = mm1k(K, 1.0, 1.5);
    std::vector<double> pi0(n, 0.0);
    pi0[0] = 0.7;
    pi0[n - 1] = 0.3;
    const std::vector<std::size_t> target{n - 1};
    const std::vector<double> tset{0.0, 1.0, 5.0};
    const mc::PassageCurve<double> c = mc::ctmc_passage_time(Q, pi0, target, tset);
    CHECK(c.atom == doctest::Approx(0.3));
    // F(0) IS the atom: a passage started inside the target completed at once.
    CHECK(c.F[0] == doctest::Approx(0.3).epsilon(1e-12));
}

TEST_CASE("ctmc_hitting_time is the first moment, and says infinity when it must") {
    const std::size_t K = 6, n = K + 1;
    const Matrix<double> Q = mm1k(K, 1.0, 1.5);
    const std::vector<std::size_t> target{n - 1};
    std::vector<double> pi0(n, 1.0 / double(n));
    const mc::PassageMoments<double> pm = mc::ctmc_passage_moments(Q, pi0, target, 1);
    const std::vector<double> h = mc::ctmc_hitting_time(Q, target);
    for (std::size_t i = 0; i < n; ++i) CHECK(h[i] == doctest::Approx(pm.mall(i, 0)));
    CHECK(h[0] == doctest::Approx(50.34375).epsilon(1e-12));

    // An isolated absorbing state cannot reach state 0.
    Matrix<double> Q2(n + 1, n + 1);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) Q2(i, j) = Q(i, j);
    const std::vector<double> h2 = mc::ctmc_hitting_time(Q2, std::vector<std::size_t>{0});
    CHECK(std::isinf(h2[n]));
}

TEST_CASE("the transform route agrees with the exponential one") {
    const std::size_t K = 6, n = K + 1;
    const Matrix<double> Q = mm1k(K, 1.0, 1.5);
    std::vector<double> pi0(n, 0.0);
    pi0[0] = 1.0;
    const std::vector<std::size_t> target{n - 1};
    std::vector<double> tset;
    for (double t = 0.5; t <= 10.0 + 1e-9; t += 0.5) tset.push_back(t);
    const mc::PassageCurve<double> a = mc::ctmc_passage_time(Q, pi0, target, tset, "expm");
    const mc::PassageCurve<double> b = mc::ctmc_passage_time(Q, pi0, target, tset, "lt");
    for (std::size_t i = 0; i < tset.size(); ++i) {
        CHECK(b.F[i] == doctest::Approx(a.F[i]).epsilon(1e-8));
        CHECK(b.f[i] == doctest::Approx(a.f[i]).epsilon(1e-8));
    }
    // L(0) = 1 for a proper law.
    const std::vector<std::complex<double>> s{std::complex<double>(0.0, 0.0)};
    const std::vector<std::complex<double>> L = mc::ctmc_passage_lst(Q, pi0, target, s);
    CHECK(L[0].real() == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("the semi-Markov recursion reproduces the Markov one on a chain that is both") {
    const std::size_t K = 6, n = K + 1;
    const Matrix<double> Q = mm1k(K, 1.0, 1.5);
    std::vector<double> pi0(n, 0.0);
    pi0[0] = 1.0;
    const std::vector<std::size_t> target{n - 1};
    const std::size_t nmax = 4;

    Matrix<double> P(n, n);
    std::vector<double> rate(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        rate[i] = -Q(i, i);
        if (rate[i] > 0.0) {
            for (std::size_t j = 0; j < n; ++j) P(i, j) = (j == i) ? 0.0 : Q(i, j) / rate[i];
        } else {
            P(i, i) = 1.0;
        }
    }
    Matrix<double> hmom(n, nmax);
    for (std::size_t i = 0; i < n; ++i) {
        double f = 1.0, p = 1.0;
        for (std::size_t r = 1; r <= nmax; ++r) {
            f *= double(r);
            p *= rate[i];
            hmom(i, r - 1) = (rate[i] > 0.0) ? f / p : 0.0;
        }
    }

    const mc::PassageMoments<double> mC = mc::ctmc_passage_moments(Q, pi0, target, nmax);
    const mc::PassageMoments<double> mS = mc::smp_passage_moments(P, hmom, pi0, target, nmax);
    for (std::size_t q = 0; q < nmax; ++q)
        CHECK(mS.m[q] == doctest::Approx(mC.m[q]).epsilon(1e-11));

    // and the transform, through the same equivalence
    std::vector<std::function<std::complex<double>(std::complex<double>)>> hlst(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double r = rate[i];
        hlst[i] = [r](std::complex<double> s) {
            return std::complex<double>(r, 0.0) / (s + std::complex<double>(r, 0.0));
        };
    }
    const std::vector<std::complex<double>> sv{std::complex<double>(0.05, 0.0),
                                               std::complex<double>(0.2, 0.0)};
    const std::vector<std::complex<double>> Ls = mc::smp_passage_lst(P, hlst, pi0, target, sv);
    const std::vector<std::complex<double>> Lc = mc::ctmc_passage_lst(Q, pi0, target, sv);
    for (std::size_t i = 0; i < sv.size(); ++i)
        CHECK(Ls[i].real() == doctest::Approx(Lc[i].real()).epsilon(1e-11));
}

TEST_CASE("non-Markov holding times give their exact analytic moments") {
    // 1 -> 2 with a deterministic sojourn of 2: the passage IS 2.
    Matrix<double> P(2, 2);
    P(0, 1) = 1.0;
    P(1, 1) = 1.0;
    Matrix<double> hmom(2, 3);
    hmom(0, 0) = 2.0;
    hmom(0, 1) = 4.0;
    hmom(0, 2) = 8.0;
    const std::vector<double> pi0{1.0, 0.0};
    const mc::PassageMoments<double> m =
        mc::smp_passage_moments(P, hmom, pi0, std::vector<std::size_t>{1}, 3);
    CHECK(m.m[0] == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(m.m[1] == doctest::Approx(4.0).epsilon(1e-12));
    CHECK(m.m[2] == doctest::Approx(8.0).epsilon(1e-12));

    // 1 -> 2 -> 3 with Exp(2) sojourns: the passage 1 -> 3 is Erlang(2,2).
    Matrix<double> P3(3, 3);
    P3(0, 1) = 1.0;
    P3(1, 2) = 1.0;
    P3(2, 2) = 1.0;
    Matrix<double> h3(3, 3);
    for (std::size_t i = 0; i < 2; ++i) {
        h3(i, 0) = 0.5;
        h3(i, 1) = 0.5;
        h3(i, 2) = 0.75;
    }
    const std::vector<double> p3{1.0, 0.0, 0.0};
    const mc::PassageMoments<double> m3 =
        mc::smp_passage_moments(P3, h3, p3, std::vector<std::size_t>{2}, 3);
    CHECK(m3.m[0] == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(m3.m[1] == doctest::Approx(1.5).epsilon(1e-12));
    CHECK(m3.m[2] == doctest::Approx(3.0).epsilon(1e-12));
}

TEST_CASE("the refusals are by name") {
    const Matrix<double> Q = mm1k(3, 1.0, 1.5);
    std::vector<double> pi0(4, 0.0);
    pi0[0] = 1.0;
    CHECK_THROWS_AS(mc::ctmc_passage_ph(Q, pi0, std::vector<std::size_t>{}), line::InputError);
    CHECK_THROWS_AS(mc::ctmc_passage_ph(Q, pi0, std::vector<std::size_t>{99}), line::InputError);
    CHECK_THROWS_AS(mc::ctmc_passage_ph(Q, std::vector<double>(2, 0.5),
                                        std::vector<std::size_t>{3}),
                    line::InputError);
    // A matrix whose rows do not sum to zero is not a generator, and is refused
    // rather than silently repaired.
    Matrix<double> bad = Q;
    bad(0, 0) += 1.0;
    CHECK_THROWS_AS(mc::ctmc_passage_ph(bad, pi0, std::vector<std::size_t>{3}), line::InputError);
    CHECK_THROWS_AS(
        mc::ctmc_passage_time(Q, pi0, std::vector<std::size_t>{3}, std::vector<double>{1.0}, "nope"),
        line::InputError);
}
