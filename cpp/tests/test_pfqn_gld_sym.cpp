/**
 * @file The SYMBOLIC arm of the gld family, against the numeric one.
 *
 * The whole point of this file is that THERE IS NO SEPARATE SYMBOLIC ARM.
 * `pfqn_gld` and `pfqn_gldsingle` are already `template <class T>` over
 * `Matrix<T>`, and `num_traits<Sym>` is a fifth numeric type beside `double`,
 * `Rational`, `Real<D>` and `Complex`; instantiating the SAME algorithm at
 * `Sym` is what makes it symbolic. So these cases are as much a test of that
 * plugin point as of the arithmetic.
 *
 * EVERY ASSERTION PINS THE SUBSTITUTED VALUE against the double instantiation,
 * as the JAR's Pfqn_gldSymTest and the native python test_pfqn_symbolic.py do,
 * and for the same reason: a wrong expression must not pass by being merely
 * well-formed.
 *
 * The whole file is behind LINE_MP_USE_SYMENGINE, which is OFF by default, so
 * a host without SymEngine compiles this to nothing rather than failing.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
#include "doctest.h"

#ifdef LINE_MP_USE_SYMENGINE

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <symengine/subs.h>

#include "line/api/pfqn/pfqn_gld.h"
#include "line/api/pfqn/pfqn_gldsingle.h"
#include "line/num/sym_number.h"

using namespace line;

namespace {

/** Substitute every symbol and read the value as a double. */
double at(const Sym& e, const std::map<std::string, double>& vals) {
    SymEngine::map_basic_basic sub;
    for (std::map<std::string, double>::const_iterator it = vals.begin(); it != vals.end(); ++it)
        sub[SymEngine::symbol(it->first)] = SymEngine::real_double(it->second);
    return SymEngine::eval_double(*SymEngine::xreplace(e.expr().get_basic(), sub));
}

Matrix<double> md(const std::vector<std::vector<double> >& v) {
    Matrix<double> m(v.size(), v.empty() ? 0 : v[0].size(), 0.0);
    for (std::size_t i = 0; i < v.size(); ++i)
        for (std::size_t j = 0; j < v[i].size(); ++j) m(i, j) = v[i][j];
    return m;
}

}  // namespace

TEST_CASE("pfqn_gldsingle carries a symbolic RATE row with no threshold") {
    // mu = (m1, m2, c): a multiserver row nobody has to locate the settling
    // index of, which is the property the gld family has and lld does not.
    Matrix<Sym> L(2, 1, Sym(0));
    L(0, 0) = Sym(2);
    L(1, 0) = Sym(3);
    Matrix<Sym> mu(2, 3, Sym(1));
    mu(0, 0) = Sym::var("m1");
    mu(0, 1) = Sym::var("m2");
    mu(0, 2) = Sym::var("c");

    const std::vector<int> N(1, 3);
    const Sym G = pfqn::pfqn_gldsingle(L, 3, mu).G;

    std::map<std::string, double> ones;
    ones["m1"] = 1.0; ones["m2"] = 1.0; ones["c"] = 1.0;
    const double ref1 = pfqn::pfqn_gldsingle(md({{2.0}, {3.0}}), 3,
                                             md({{1, 1, 1}, {1, 1, 1}})).G;
    CHECK(at(G, ones) == doctest::Approx(ref1).epsilon(1e-12));

    std::map<std::string, double> ms;
    ms["m1"] = 1.0; ms["m2"] = 2.0; ms["c"] = 2.0;
    const double ref2 = pfqn::pfqn_gldsingle(md({{2.0}, {3.0}}), 3,
                                             md({{1, 2, 2}, {1, 1, 1}})).G;
    CHECK(at(G, ms) == doctest::Approx(ref2).epsilon(1e-12));
}

TEST_CASE("pfqn_gld carries symbolic DEMANDS through the multiclass recursion") {
    Matrix<Sym> L(2, 2, Sym(0));
    L(0, 0) = Sym::var("L11"); L(0, 1) = Sym::var("L12");
    L(1, 0) = Sym::var("L21"); L(1, 1) = Sym::var("L22");
    Matrix<Sym> mu(2, 4, Sym(1));

    std::vector<int> N(2);
    N[0] = 2; N[1] = 2;
    const Sym G = pfqn::pfqn_gld(L, N, mu).G;

    std::map<std::string, double> vals;
    vals["L11"] = 1.0; vals["L12"] = 0.6; vals["L21"] = 0.5; vals["L22"] = 1.1;
    const double ref = pfqn::pfqn_gld(md({{1.0, 0.6}, {0.5, 1.1}}), N,
                                      md({{1, 1, 1, 1}, {1, 1, 1, 1}})).G;
    CHECK(at(G, vals) == doctest::Approx(ref).epsilon(1e-12));
}

TEST_CASE("num_traits<Sym> declares the type the gld family expects") {
    CHECK(num_traits<Sym>::is_exact);
    // FALSE steers the family off its log-domain arm, as it does for Rational:
    // that arm exists to stop a double underflowing and an exact expression
    // cannot underflow.
    CHECK_FALSE(num_traits<Sym>::has_transcendental);
    CHECK(std::string(num_traits<Sym>::name()) == "symbolic");
    CHECK(num_traits<Sym>::from_int(7) == Sym(7));
    CHECK(num_traits<Sym>::to_double(Sym(7)) == doctest::Approx(7.0));
}

#endif  // LINE_MP_USE_SYMENGINE
