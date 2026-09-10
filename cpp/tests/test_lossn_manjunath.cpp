/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The Manjunath-Sikdar transform, against oracles that do not come out of the
 * implementation: Erlang's loss formula by its rational recursion, exact box
 * enumeration of the admissible set, and the truncated Poisson tail.
 *
 * THE ENUMERATION IS THE PRIMARY ORACLE. `enumerate` sums the product form over
 * the whole box {0..N_r} and keeps the feasible states, which is the definition
 * of g(C); the transform is an algebraic identity for that sum, so agreement to
 * 1e-12 is the correct expectation and not a tolerance chosen to pass. It is
 * written independently here (a nested odometer, no series, no residues) rather
 * than shared with the implementation.
 *
 * Ported from matlab/src/api/lossn/test_lossn_manjunath.m; the enumeration oracle is
 * its `lossn_enum`, and the assertions are the same eight.
 */
#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/lossn/lossn_manjunath.h"
#include "line/api/lossn/lossn_mci.h"
#include "line/num/number.h"
#include "line/util/error.h"

using line::Matrix;
using line::Rational;
using line::lossn::LossnManjunathOptions;
using line::lossn::lossn_manjunath;

namespace {

/** Erlang B by the rational recursion B_k = nu B_{k-1} / (k + nu B_{k-1}). */
double erlangB(double nu, int C) {
    double b = 1.0;
    for (int k = 1; k <= C; ++k) b = nu * b / (k + nu * b);
    return b;
}

struct EnumResult {
    double g = 0.0;
    std::vector<double> loss;
};

/**
 * Exact g(C) and per-route blocking by box enumeration, the port of the
 * reference test's `lossn_enum`. The product form is accumulated directly (no
 * log-sum-exp) because every case here is small enough that nu^n/n! stays in
 * range, which keeps the oracle as simple as possible.
 */
EnumResult enumerate(const std::vector<double>& nu, const Matrix<double>& A,
                     const std::vector<double>& C) {
    const std::size_t R = nu.size(), J = C.size();
    std::vector<long> N(R, 0);
    for (std::size_t r = 0; r < R; ++r) {
        long v = -1;
        for (std::size_t j = 0; j < J; ++j) {
            if (A(j, r) <= 0.0) continue;
            const long cap = static_cast<long>(std::floor(C[j] / A(j, r)));
            if (v < 0 || cap < v) v = cap;
        }
        N[r] = v < 0 ? 0 : v;  // a route in no row is unbounded; none appear here
    }

    // q(n) = prod_r nu_r^{n_r} / n_r!, tabulated per route.
    std::vector<std::vector<double>> term(R);
    for (std::size_t r = 0; r < R; ++r) {
        term[r].assign(static_cast<std::size_t>(N[r]) + 1, 1.0);
        for (std::size_t n = 1; n < term[r].size(); ++n)
            term[r][n] = term[r][n - 1] * nu[r] / static_cast<double>(n);
    }

    EnumResult out;
    out.loss.assign(R, 0.0);
    std::vector<double> gShift(R, 0.0);
    std::vector<bool> impossible(R, false);
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t j = 0; j < J; ++j)
            if (C[j] - A(j, r) < 0.0) impossible[r] = true;

    std::vector<long> n(R, 0);
    while (true) {
        double q = 1.0;
        for (std::size_t r = 0; r < R; ++r) q *= term[r][static_cast<std::size_t>(n[r])];
        std::vector<double> load(J, 0.0);
        bool feas = true;
        for (std::size_t j = 0; j < J; ++j) {
            for (std::size_t r = 0; r < R; ++r) load[j] += A(j, r) * static_cast<double>(n[r]);
            if (load[j] > C[j] + 1e-9) feas = false;
        }
        if (feas) {
            out.g += q;
            // The class r acceptance region is A n <= C - A e_r.
            for (std::size_t r = 0; r < R; ++r) {
                if (impossible[r]) continue;
                bool ok = true;
                for (std::size_t j = 0; j < J && ok; ++j)
                    if (load[j] > C[j] - A(j, r) + 1e-9) ok = false;
                if (ok) gShift[r] += q;
            }
        }
        std::size_t k = 0;
        for (; k < R; ++k) {
            if (++n[k] <= N[k]) break;
            n[k] = 0;
        }
        if (k == R) break;
    }
    for (std::size_t r = 0; r < R; ++r)
        out.loss[r] = impossible[r] ? 1.0 : 1.0 - gShift[r] / out.g;
    return out;
}

Matrix<double> mat(std::size_t J, std::size_t R, const std::vector<double>& rowmajor) {
    Matrix<double> A(J, R, 0.0);
    for (std::size_t j = 0; j < J; ++j)
        for (std::size_t r = 0; r < R; ++r) A(j, r) = rowmajor[j * R + r];
    return A;
}

}  // namespace

TEST_CASE("lossn_manjunath: a single link with unit demand is Erlang's loss system") {
    // T1 of the reference test. nu = 8, C = 10.
    const double nu = 8.0;
    const int C = 10;
    Matrix<double> A(1, 1, 1.0);
    const auto r = lossn_manjunath(std::vector<double>{nu}, A, std::vector<double>{double(C)});

    CHECK(r.Loss[0] == doctest::Approx(erlangB(nu, C)).epsilon(1e-13));
    CHECK(r.QLen[0] == doctest::Approx(nu * (1.0 - r.Loss[0])).epsilon(1e-13));
    CHECK(r.iterations == 1);

    // g(C) = sum_{n=0}^{C} nu^n / n!, accumulated independently.
    double g = 0.0, t = 1.0;
    for (int n = 0; n <= C; ++n) {
        if (n > 0) t *= nu / n;
        g += t;
    }
    CHECK(r.lG == doctest::Approx(std::log(g)).epsilon(1e-13));
}

TEST_CASE("lossn_manjunath: two-link multirate network matches box enumeration") {
    // T2 of the reference test.
    const std::vector<double> nu{3.0, 1.5}, C{4.0, 5.0};
    const Matrix<double> A = mat(2, 2, {1, 1, 1, 2});
    const auto r = lossn_manjunath(nu, A, C);
    const EnumResult e = enumerate(nu, A, C);

    for (std::size_t k = 0; k < 2; ++k) {
        INFO("route ", k);
        CHECK(r.Loss[k] == doctest::Approx(e.loss[k]).epsilon(1e-12));
        CHECK(r.QLen[k] == doctest::Approx(nu[k] * (1.0 - r.Loss[k])).epsilon(1e-13));
    }
    CHECK(r.lG == doctest::Approx(std::log(e.g)).epsilon(1e-12));
}

TEST_CASE("lossn_manjunath: job cap, memory budget and two linear rows match enumeration") {
    // T3 of the reference test: the constraint shape a finite capacity region
    // produces from setGlobalMaxJobs / setGlobalMaxMemory with setClassSize and
    // setConstraint. Four rows over three classes, every row sharing classes
    // with every other, so nothing here decouples.
    const std::vector<double> nu{4.0, 2.5, 2.0}, C{9.0, 14.0, 10.0, 12.0};
    const Matrix<double> A = mat(4, 3, {1, 1, 1, 1, 2, 3, 2, 1, 1, 1, 3, 2});
    const auto r = lossn_manjunath(nu, A, C);
    const EnumResult e = enumerate(nu, A, C);

    for (std::size_t k = 0; k < 3; ++k) {
        INFO("route ", k);
        CHECK(r.Loss[k] == doctest::Approx(e.loss[k]).epsilon(1e-12));
    }
    CHECK(r.lG == doctest::Approx(std::log(e.g)).epsilon(1e-12));
    // Four rows live at once would be 10*15*11*13 coefficients; the elimination
    // order never holds more than the induced width, and every row shares a
    // class with route 0, so here that IS all four.
    CHECK(r.peak_states > 1);
}

TEST_CASE("lossn_manjunath: the gcd row reduction is exact") {
    // T4 of the reference test: 2 n1 + 4 n2 <= 9 and n1 + 2 n2 <= 4 have the
    // same integer solutions, so the two calls must agree exactly.
    const std::vector<double> nu{2.0, 1.0};
    const auto a = lossn_manjunath(nu, mat(1, 2, {2, 4}), std::vector<double>{9.0});
    const auto b = lossn_manjunath(nu, mat(1, 2, {1, 2}), std::vector<double>{4.0});
    for (std::size_t k = 0; k < 2; ++k) {
        INFO("route ", k);
        CHECK(a.Loss[k] == doctest::Approx(b.Loss[k]).epsilon(1e-13));
    }
    CHECK(a.lG == doctest::Approx(b.lG).epsilon(1e-13));
}

TEST_CASE("lossn_manjunath: a route in no constraint never blocks and factors out") {
    // T5 of the reference test.
    const std::vector<double> nu{2.0, 1.0, 1.5}, C{6.0};
    const Matrix<double> A = mat(1, 3, {1, 1, 0});
    const auto r = lossn_manjunath(nu, A, C);

    CHECK(r.Loss[2] == 0.0);
    CHECK(r.QLen[2] == doctest::Approx(nu[2]).epsilon(1e-13));

    // Dropping the free route perturbs neither of the others.
    const auto two =
        lossn_manjunath(std::vector<double>{nu[0], nu[1]}, mat(1, 2, {1, 1}), C);
    for (std::size_t k = 0; k < 2; ++k) {
        INFO("route ", k);
        CHECK(r.Loss[k] == doctest::Approx(two.Loss[k]).epsilon(1e-13));
    }
    // and contributes exactly the factor exp(nu_3).
    const EnumResult e = enumerate({nu[0], nu[1]}, mat(1, 2, {1, 1}), C);
    CHECK(r.lG == doctest::Approx(std::log(e.g) + nu[2]).epsilon(1e-12));
}

TEST_CASE("lossn_manjunath: a region with no bounded row at all is pure Poisson") {
    // The limit of the case above with EVERY route free: there is no admission
    // rule, g(C) = prod exp(nu_r), and no series is built.
    const std::vector<double> nu{2.0, 1.0};
    const auto r = lossn_manjunath(nu, mat(1, 2, {0, 0}), std::vector<double>{5.0});
    CHECK(r.Loss[0] == 0.0);
    CHECK(r.Loss[1] == 0.0);
    CHECK(r.QLen[0] == doctest::Approx(2.0));
    CHECK(r.QLen[1] == doctest::Approx(1.0));
    CHECK(r.lG == doctest::Approx(3.0).epsilon(1e-13));
}

TEST_CASE("lossn_manjunath: the exact transform lies inside the Monte Carlo interval") {
    // T6 of the reference test. The sampler is an independent algorithm, so this
    // is a cross-check of two ports against each other, not against a formula.
    const std::vector<double> nu{3.0, 1.5}, C{4.0, 5.0};
    const Matrix<double> A = mat(2, 2, {1, 1, 1, 2});
    const auto exact = lossn_manjunath(nu, A, C);

    line::lossn::LossnMciOptions<double> mo;
    mo.samples = 300000;
    const auto mc = line::lossn::lossn_mci<double>(nu, A, C, mo, 5);
    for (std::size_t k = 0; k < 2; ++k) {
        INFO("route ", k, ": exact ", exact.Loss[k], " CI [", mc.lossCI(k, 0), ", ",
             mc.lossCI(k, 1), "]");
        CHECK(mc.lossCI(k, 0) <= exact.Loss[k]);
        CHECK(exact.Loss[k] <= mc.lossCI(k, 1));
    }
}

TEST_CASE("lossn_manjunath: rare blocking is exact where a sampler cannot resolve it") {
    // T7 of the reference test: two routes on one link with unit demand is an
    // Erlang system at the aggregate load, so blocking is the truncated Poisson
    // tail at the capacity.
    const std::vector<double> nu{1.0, 0.5};
    const auto r = lossn_manjunath(nu, mat(1, 2, {1, 1}), std::vector<double>{8.0});

    double s = 0.0, t = 1.0;
    std::vector<double> p(9, 0.0);
    for (int k = 0; k <= 8; ++k) {
        if (k > 0) t *= 1.5 / k;
        p[k] = t;
        s += t;
    }
    CHECK(r.Loss[0] == doctest::Approx(p[8] / s).epsilon(1e-14));
    CHECK(r.Loss[1] == doctest::Approx(p[8] / s).epsilon(1e-14));
    CHECK(r.Loss[0] < 1e-3);  // the regime where a sampled throughput fails
}

TEST_CASE("lossn_manjunath: a fractional rule is refused, never rounded") {
    // T8 of the reference test.
    CHECK_THROWS_AS(lossn_manjunath(std::vector<double>{1.0, 1.0}, mat(1, 2, {1, 1.5}),
                             std::vector<double>{5.0}),
                    line::InputError);
    CHECK_THROWS_AS(lossn_manjunath(std::vector<double>{1.0}, mat(1, 1, {1.0}),
                             std::vector<double>{5.5}),
                    line::InputError);
    // Negative entries are refused for the same reason.
    CHECK_THROWS_AS(lossn_manjunath(std::vector<double>{1.0}, mat(1, 1, {-1.0}),
                             std::vector<double>{5.0}),
                    line::InputError);
    CHECK_THROWS_AS(lossn_manjunath(std::vector<double>{-1.0}, mat(1, 1, {1.0}),
                             std::vector<double>{5.0}),
                    line::InputError);
}

TEST_CASE("lossn_manjunath: a route that cannot fit even once is blocked with certainty") {
    // A capacity of 1 against a demand of 2: no state admits the route, not even
    // the empty one, so the ratio g(C - A e_r)/g(C) is not merely small but
    // undefined, and the branch that detects it must report loss 1 exactly.
    const auto r = lossn_manjunath(std::vector<double>{1.0, 1.0}, mat(1, 2, {1, 2}),
                            std::vector<double>{1.0});
    CHECK(r.Loss[1] == 1.0);
    CHECK(r.QLen[1] == 0.0);
    CHECK(r.Loss[0] > 0.0);
    CHECK(r.Loss[0] < 1.0);
}

TEST_CASE("lossn_manjunath: zero offered load neither blocks nor loads") {
    const std::vector<double> nu{4.0, 0.0}, C{5.0};
    const auto r = lossn_manjunath(nu, mat(1, 2, {1, 1}), C);
    CHECK(r.QLen[1] == doctest::Approx(0.0));
    // The idle route sees the same congestion as the busy one: on one link with
    // unit demands, blocking is a property of the state, not of the route.
    CHECK(r.Loss[1] == doctest::Approx(r.Loss[0]).epsilon(1e-13));
    CHECK(r.Loss[0] == doctest::Approx(erlangB(4.0, 5)).epsilon(1e-13));
}

TEST_CASE("lossn_manjunath: blocking rises with load and falls with capacity") {
    double prev = -1.0;
    for (double nu : {0.5, 2.0, 6.0, 20.0}) {
        const auto r = lossn_manjunath(std::vector<double>{nu}, Matrix<double>(1, 1, 1.0),
                                std::vector<double>{4.0});
        CHECK(r.Loss[0] > prev);
        prev = r.Loss[0];
    }
    double last = 2.0;
    for (double C : {2.0, 4.0, 8.0, 16.0}) {
        const auto r = lossn_manjunath(std::vector<double>{4.0}, Matrix<double>(1, 1, 1.0),
                                std::vector<double>{C});
        CHECK(r.Loss[0] < last);
        last = r.Loss[0];
    }
}

TEST_CASE("lossn_manjunath: the live-coefficient cap refuses by name") {
    LossnManjunathOptions tiny;
    tiny.max_live_states = 4;
    CHECK_THROWS_AS(lossn_manjunath(std::vector<double>{1.0}, Matrix<double>(1, 1, 1.0),
                             std::vector<double>{100.0}, tiny),
                    line::UnsupportedError);
    std::string what;
    try {
        lossn_manjunath(std::vector<double>{1.0}, Matrix<double>(1, 1, 1.0),
                 std::vector<double>{100.0}, tiny);
    } catch (const line::UnsupportedError& e) {
        what = e.what();
    }
    CHECK(what.find("lossn_mci") != std::string::npos);
}

TEST_CASE("lossn_manjunath: the transform runs at exact arithmetic") {
    // The distinguishing property of this analyzer among the three: every
    // operation on the series is rational and the reported metrics are RATIOS,
    // so Loss is a rational number with no rounding anywhere in its derivation.
    // The oracle is Erlang B as an exact rational, built by the same recursion
    // the double test uses.
    const int C = 6;
    const Rational nu = Rational(7, 2);
    Matrix<Rational> A(1, 1, Rational(1));
    const auto r = lossn_manjunath<Rational>({nu}, A, {Rational(C)});

    Rational b(1);
    for (int k = 1; k <= C; ++k) b = nu * b / (Rational(k) + nu * b);
    CHECK(r.Loss[0] == b);  // exact equality, not a tolerance
    CHECK(r.QLen[0] == nu * (Rational(1) - b));

    // g(C) = sum nu^n/n! as an exact rational, checked through its logarithm
    // since that is what the result carries.
    Rational g(0), t(1);
    for (int n = 0; n <= C; ++n) {
        if (n > 0) t = t * nu / Rational(n);
        g += t;
    }
    CHECK(r.lG == doctest::Approx(line::num_traits<Rational>::log_as_double(g)).epsilon(1e-12));

    // And it agrees with the double instantiation to the precision of double.
    const auto d = lossn_manjunath(std::vector<double>{3.5}, Matrix<double>(1, 1, 1.0),
                            std::vector<double>{double(C)});
    CHECK(d.Loss[0] == doctest::Approx(static_cast<double>(r.Loss[0])).epsilon(1e-13));
}

TEST_CASE("lossn_manjunath: exact and double agree on a coupled multirate rule") {
    // The single-link case above cannot exercise the shift-and-accumulate
    // convolution across two live links; this one does, and exact arithmetic
    // makes the double result's error visible rather than assumed.
    const Matrix<double> Ad = mat(2, 2, {1, 1, 1, 2});
    const std::vector<double> nud{3.0, 3.0 / 2.0}, Cd{4.0, 5.0};
    const auto d = lossn_manjunath(nud, Ad, Cd);

    Matrix<Rational> Ae(2, 2, Rational(0));
    Ae(0, 0) = Rational(1);
    Ae(0, 1) = Rational(1);
    Ae(1, 0) = Rational(1);
    Ae(1, 1) = Rational(2);
    const auto e = lossn_manjunath<Rational>({Rational(3), Rational(3, 2)}, Ae,
                                      {Rational(4), Rational(5)});
    for (std::size_t k = 0; k < 2; ++k) {
        INFO("route ", k);
        CHECK(d.Loss[k] == doctest::Approx(static_cast<double>(e.Loss[k])).epsilon(1e-12));
        CHECK(d.QLen[k] == doctest::Approx(static_cast<double>(e.QLen[k])).epsilon(1e-12));
    }
    CHECK(d.lG == doctest::Approx(e.lG).epsilon(1e-12));
}

TEST_CASE("lossn_manjunath: a large load does not overflow the series terms") {
    // nu^n/n! peaks near e^nu, so a load of 900 overflows a double if the terms
    // are formed directly. The reference forms them in log space and rescales;
    // the port does the same, and the test is that the answer stays finite and
    // sane rather than NaN. Blocking at C = nu is close to a half by the
    // square-root staffing heuristic; the check is only that it is in range.
    const auto r = lossn_manjunath(std::vector<double>{900.0}, Matrix<double>(1, 1, 1.0),
                            std::vector<double>{900.0});
    CHECK(std::isfinite(r.Loss[0]));
    CHECK(r.Loss[0] > 0.0);
    CHECK(r.Loss[0] < 0.1);
    CHECK(std::isfinite(r.lG));
    CHECK(r.QLen[0] == doctest::Approx(900.0 * (1.0 - r.Loss[0])).epsilon(1e-10));
    // Erlang B at C = nu, an independent recursion in the same regime.
    CHECK(r.Loss[0] == doctest::Approx(erlangB(900.0, 900)).epsilon(1e-9));
}
