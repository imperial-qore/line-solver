/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Gerasimov's residue closed form for the normalizing constant, generalized to R
 * classes. A. I. Gerasimov, "On Normalizing Constants in Multiclass Queueing
 * Networks", Operations Research 43(4):704-711, 1995.
 *
 * The paper's own worked example (Figure 1, Table I, and the two closed forms in
 * the Appendix) pins the R = 2 algorithm rather than this port. One superscript is
 * lost in the scan of the second Appendix form: its leading term prints as
 * x21/(x21-x11)^3 where the series it sums requires x21^(N1+4)/(x21-x11)^3, which
 * is what is asserted here. Everything else is checked against pfqn_ca, which is
 * exact and independent. Values agree across MATLAB, JAR, Python native and C++.
 *
 * The exact backend matters here more than anywhere else in the pfqn family: the
 * residue expansion is an alternating sum, and in Rational it carries no
 * cancellation at all, so pfqn_gerasimov and pfqn_ca must agree BIT FOR BIT.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_gerasimov.h"

using line::Matrix;
using line::Rational;
using line::num_traits;
using namespace line::pfqn;

namespace {

const double Y = 0.05202;    // Table I: x12 = x32
const double X11 = 0.06627;  // Table I: x11

template <class T>
Matrix<T> mk(const std::vector<std::vector<double>>& a) {
    Matrix<T> m(a.size(), a[0].size());
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < a[0].size(); ++j) m(i, j) = num_traits<T>::from_double(a[i][j]);
    return m;
}

Matrix<double> row(const std::vector<double>& z) {
    Matrix<double> m(1, z.size());
    for (std::size_t r = 0; r < z.size(); ++r) m(0, r) = z[r];
    return m;
}

void check(const char* label, const std::vector<std::vector<double>>& L,
           const std::vector<int>& N, const std::vector<double>& Z) {
    INFO(label);
    const NcResult<double> ge = pfqn_gerasimov(mk<double>(L), N, row(Z));
    const NcResult<double> ca = pfqn_ca(mk<double>(L), N, row(Z));
    CHECK(ge.G == doctest::Approx(ca.G).epsilon(1e-12));
}

}  // namespace

TEST_CASE("pfqn_gerasimov reproduces the paper's Appendix form (i), x11 = x21 = c") {
    const double c = 0.08;
    for (int n1 = 1; n1 <= 6; ++n1) {
        const NcResult<double> g =
            pfqn_gerasimov(mk<double>({{c, Y}, {c, 0}, {0, Y}}), std::vector<int>{n1, 2});
        const double gf = (Y * Y / 6) * std::pow(c, n1) *
                          (1.0 * n1 * n1 * n1 + 9.0 * n1 * n1 + 26.0 * n1 + 18.0);
        CHECK(g.G == doctest::Approx(gf).epsilon(1e-12));
    }
}

TEST_CASE("pfqn_gerasimov reproduces the paper's Appendix form (ii), x11 != x21") {
    const double x21 = 5.0;
    for (int n1 = 1; n1 <= 6; ++n1) {
        const NcResult<double> g =
            pfqn_gerasimov(mk<double>({{X11, Y}, {x21, 0}, {0, Y}}), std::vector<int>{n1, 2});
        const double gf = (Y * Y / X11) *
            (std::pow(x21, n1 + 4) / std::pow(x21 - X11, 3) +
             (1.0 * n1 * n1 + 7.0 * n1 + 12.0) * std::pow(X11, n1 + 2) / (2 * (X11 - x21)) -
             (n1 + 4) * std::pow(X11, n1 + 3) / std::pow(X11 - x21, 2) +
             std::pow(X11, n1 + 4) / std::pow(x21 - X11, 3) - std::pow(x21, n1 + 1));
        CHECK(g.G == doctest::Approx(gf).epsilon(1e-7));
    }
}

TEST_CASE("pfqn_gerasimov matches the convolution algorithm") {
    // The degeneracies the paper's hypotheses exclude are ordinary cases here:
    // coincident poles, tied class-2 demands, stations unvisited by a class and
    // identical station rows.
    check("fig1", {{X11, Y}, {5.0, 0}, {0, Y}}, {3, 2}, {0, 0});
    check("coincident poles x11=x21", {{0.08, Y}, {0.08, 0}, {0, Y}}, {4, 2}, {0, 0});
    check("tied class-2 demands", {{1, 2}, {3, 2}, {5, 7}}, {3, 2}, {0, 0});
    check("identical station rows", {{1, 2}, {1, 2}, {5, 7}}, {3, 3}, {0, 0});
    check("class-2 unvisited stations", {{1, 0}, {3, 0}, {5, 7}}, {4, 2}, {0, 0});
    check("think time", {{1, 2}, {3, 4}}, {3, 2}, {0.5, 1.5});
    check("think time one class", {{1, 2}, {3, 4}, {2, 1}}, {2, 3}, {1, 0});
    check("single class", {{1}, {2}, {3}}, {6}, {0});
    check("single class with delay", {{1}, {2}, {3}}, {6}, {2});
    check("empty class", {{1, 2}, {3, 4}}, {4, 0}, {0, 0});
    check("three classes", {{1, 2, 3}, {3, 4, 1}, {2, 1, 2}}, {2, 2, 2}, {0, 0, 0});
    check("three classes with delay", {{1, 2, 3}, {3, 4, 1}, {2, 1, 2}}, {3, 2, 1}, {0.7, 0, 0.3});
    check("four classes", {{1, 2, 3, 1}, {3, 4, 1, 2}, {2, 1, 2, 3}}, {2, 1, 2, 1}, {0, 0, 0, 0});
    check("four classes with delay", {{1, 2, 3, 1}, {3, 4, 1, 2}, {2, 1, 2, 3}}, {2, 1, 2, 1},
          {0.3, 0.2, 0, 0.1});
}

TEST_CASE("pfqn_gerasimov leaves the eliminated population free") {
    // A population removed by residues enters only as a pole ORDER, so it costs
    // nothing: the answer must stay exact as it grows by three orders of magnitude.
    const std::vector<std::vector<double>> L = {{1, 2}, {3, 1}, {2, 4}, {0.5, 0.7}};
    for (int n2 : {10, 100, 2000, 20000}) {
        const NcResult<double> ge = pfqn_gerasimov(mk<double>(L), std::vector<int>{6, n2});
        const NcResult<double> ca = pfqn_ca(mk<double>(L), std::vector<int>{6, n2}, Matrix<double>());
        CHECK(ge.lG == doctest::Approx(ca.lG).epsilon(1e-12));
    }
}

TEST_CASE("pfqn_gerasimov is exact in the exact backend") {
    // Every operation is a product, a quotient or a binomial coefficient, so the
    // alternating residue sum carries no cancellation at all in Rational.
    const std::vector<std::vector<double>> L = {{1, 2, 3}, {3, 4, 1}, {2, 1, 2}};
    Matrix<Rational> Ze(1, 3);
    Ze(0, 0) = num_traits<Rational>::from_rational(7, 10);
    Ze(0, 1) = num_traits<Rational>::from_int(0);
    Ze(0, 2) = num_traits<Rational>::from_rational(3, 10);
    const NcResult<Rational> ge = pfqn_gerasimov(mk<Rational>(L), std::vector<int>{3, 2, 1}, Ze);
    const NcResult<Rational> ca = pfqn_ca(mk<Rational>(L), std::vector<int>{3, 2, 1}, Ze);
    CHECK(ge.G == ca.G);

    const NcResult<Rational> ge2 = pfqn_gerasimov(mk<Rational>({{1, 2}, {1, 2}, {5, 7}}),
                                                 std::vector<int>{3, 3}, Matrix<Rational>());
    const NcResult<Rational> ca2 = pfqn_ca(mk<Rational>({{1, 2}, {1, 2}, {5, 7}}),
                                           std::vector<int>{3, 3}, Matrix<Rational>());
    CHECK(ge2.G == ca2.G);
}

TEST_CASE("pfqn_gerasimov refuses rather than truncates") {
    // A truncated residue sum is not a bound on G, it is a wrong number.
    const std::vector<std::vector<double>> L = {{1.0, 2.0, 3.0}, {3.0, 1.0, 2.0}, {2.0, 3.0, 1.0},
                                                {1.5, 2.5, 0.5}, {0.5, 1.5, 2.5}, {2.5, 0.5, 1.5}};
    CHECK_THROWS(pfqn_gerasimov(mk<double>(L), std::vector<int>{4, 8, 8}, Matrix<double>(), 1e-12, 50));
}
