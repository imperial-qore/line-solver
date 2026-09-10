/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Primal simplex with Bland's rule. Oracles are hand-solved linear programs:
 * a two-variable problem whose optimal vertex is known in closed form, the
 * Beale instance that cycles forever under Dantzig's most-negative-reduced-cost
 * rule, an inconsistent system, and a ray. Everything is checked at
 * line::Rational as well as at double, because exactness at Rational is the
 * reason this solver exists: the optimum of an LP with rational data is
 * rational, and the test asserts it as an exact equality of rationals, not to
 * a tolerance.
 */
#include <string>
#include <vector>

#include "doctest.h"
#include "line/num/number.h"
#include "line/util/simplex.h"

using line::Rational;
using line::lp::LpModel;
using line::lp::LpSense;
using line::lp::LpStatus;
using line::lp::simplex_solve;

namespace {

Rational rat(long p, long q) { return Rational(p, q); }

}  // namespace

TEST_CASE("two-variable LP with a known vertex optimum") {
    // max 3x + 5y  s.t.  x <= 4, 2y <= 12, 3x + 2y <= 18, x,y >= 0.
    // Textbook Wyndor Glass problem: optimum at (2,6) with value 36.
    SUBCASE("double") {
        LpModel<double> m(2);
        m.set_cost(0, 3.0);
        m.set_cost(1, 5.0);
        m.row_clear();
        m.row_add(0, 1.0);
        m.emit_le(4.0);
        m.row_add(1, 2.0);
        m.emit_le(12.0);
        m.row_add(0, 3.0);
        m.row_add(1, 2.0);
        m.emit_le(18.0);
        line::lp::LpSolution<double> s = simplex_solve(m);
        REQUIRE(s.status == LpStatus::Optimal);
        CHECK(s.objective == doctest::Approx(36.0).epsilon(1e-12));
        CHECK(s.x[0] == doctest::Approx(2.0).epsilon(1e-12));
        CHECK(s.x[1] == doctest::Approx(6.0).epsilon(1e-12));
    }
    SUBCASE("exact") {
        LpModel<Rational> m(2);
        m.set_cost(0, Rational(3));
        m.set_cost(1, Rational(5));
        m.row_add(0, Rational(1));
        m.emit_le(Rational(4));
        m.row_add(1, Rational(2));
        m.emit_le(Rational(12));
        m.row_add(0, Rational(3));
        m.row_add(1, Rational(2));
        m.emit_le(Rational(18));
        line::lp::LpSolution<Rational> s = simplex_solve(m);
        REQUIRE(s.status == LpStatus::Optimal);
        CHECK(s.objective == Rational(36));  // exact, not approximate
        CHECK(s.x[0] == Rational(2));
        CHECK(s.x[1] == Rational(6));
    }
}

TEST_CASE("fractional vertex is reported as an exact rational") {
    // max x + y  s.t.  3x + 2y <= 7, x + 4y <= 6, x,y >= 0.
    // Vertex of the two active rows: x = 16/10 = 8/5, y = 11/10, value 27/10.
    LpModel<Rational> m(2);
    m.set_cost(0, Rational(1));
    m.set_cost(1, Rational(1));
    m.row_add(0, Rational(3));
    m.row_add(1, Rational(2));
    m.emit_le(Rational(7));
    m.row_add(0, Rational(1));
    m.row_add(1, Rational(4));
    m.emit_le(Rational(6));
    line::lp::LpSolution<Rational> s = simplex_solve(m);
    REQUIRE(s.status == LpStatus::Optimal);
    CHECK(s.x[0] == rat(8, 5));
    CHECK(s.x[1] == rat(11, 10));
    CHECK(s.objective == rat(27, 10));
    // The value is genuinely non-terminating in binary, so an exact solver has
    // to carry it as a fraction rather than round it.
    CHECK(line::num_traits<Rational>::denominator_str(s.objective) == std::string("10"));
}

TEST_CASE("Beale's cycling instance terminates under Bland's rule") {
    // Beale (1955). Dantzig's rule cycles through six bases forever here.
    //   max  (3/4) x1 - 150 x2 + (1/50) x3 - 6 x4
    //   s.t. (1/4) x1 -  60 x2 - (1/25) x3 + 9 x4 <= 0
    //        (1/2) x1 -  90 x2 - (1/50) x3 + 3 x4 <= 0
    //                                        x3   <= 1
    //   x >= 0.  Optimum is 1/20 at x = (1/25, 0, 1, 0).
    LpModel<Rational> m(4);
    m.set_cost(0, rat(3, 4));
    m.set_cost(1, Rational(-150));
    m.set_cost(2, rat(1, 50));
    m.set_cost(3, Rational(-6));
    m.row_add(0, rat(1, 4));
    m.row_add(1, Rational(-60));
    m.row_add(2, rat(-1, 25));
    m.row_add(3, Rational(9));
    m.emit_le(Rational(0));
    m.row_add(0, rat(1, 2));
    m.row_add(1, Rational(-90));
    m.row_add(2, rat(-1, 50));
    m.row_add(3, Rational(3));
    m.emit_le(Rational(0));
    m.row_add(2, Rational(1));
    m.emit_le(Rational(1));
    line::lp::LpSolution<Rational> s = simplex_solve(m);
    REQUIRE(s.status == LpStatus::Optimal);
    CHECK(s.objective == rat(1, 20));
    CHECK(s.x[0] == rat(1, 25));
    CHECK(s.x[1] == Rational(0));
    CHECK(s.x[2] == Rational(1));
    CHECK(s.x[3] == Rational(0));
    // Termination is the assertion: the iteration cap is never reached.
    CHECK(s.status != LpStatus::IterationLimit);
}

TEST_CASE("infeasible system is reported, not solved") {
    // x + y <= 1 and x + y >= 3 with x,y >= 0.
    LpModel<Rational> m(2);
    m.set_cost(0, Rational(1));
    m.row_add(0, Rational(1));
    m.row_add(1, Rational(1));
    m.emit_le(Rational(1));
    m.row_add(0, Rational(1));
    m.row_add(1, Rational(1));
    m.emit_ge(Rational(3));
    line::lp::LpSolution<Rational> s = simplex_solve(m);
    CHECK(s.status == LpStatus::Infeasible);
    CHECK(std::string(line::lp::lp_status_name(s.status)) == "Infeasible");
}

TEST_CASE("infeasible equality system is reported") {
    // x + y = 1 and x + y = 2.
    LpModel<double> m(2);
    m.set_cost(0, 1.0);
    m.row_add(0, 1.0);
    m.row_add(1, 1.0);
    m.emit_eq(1.0);
    m.row_add(0, 1.0);
    m.row_add(1, 1.0);
    m.emit_eq(2.0);
    CHECK(simplex_solve(m).status == LpStatus::Infeasible);
}

TEST_CASE("unbounded direction is reported") {
    // max x + y  s.t.  x - y <= 1, x,y >= 0: the ray (t,t) is feasible forever.
    LpModel<Rational> m(2);
    m.set_cost(0, Rational(1));
    m.set_cost(1, Rational(1));
    m.row_add(0, Rational(1));
    m.row_add(1, Rational(-1));
    m.emit_le(Rational(1));
    line::lp::LpSolution<Rational> s = simplex_solve(m);
    CHECK(s.status == LpStatus::Unbounded);
    // Boxing the variables makes the same model bounded.
    m.set_upper(0, Rational(5));
    m.set_upper(1, Rational(5));
    line::lp::LpSolution<Rational> s2 = simplex_solve(m);
    REQUIRE(s2.status == LpStatus::Optimal);
    CHECK(s2.objective == Rational(10));
}

TEST_CASE("equality rows with a negative right-hand side need phase 1") {
    // min x + y s.t. -x - 2y = -4, x - y <= 1, 0 <= x,y <= 3.
    // Feasible set is the segment x = 4 - 2y; minimizing x+y drives y up to
    // its largest feasible value y = 2 (x = 0), giving objective 2.
    LpModel<Rational> m(2);
    m.set_maximize(false);
    m.set_cost(0, Rational(1));
    m.set_cost(1, Rational(1));
    m.set_bounds(0, Rational(0), Rational(3));
    m.set_bounds(1, Rational(0), Rational(3));
    m.row_add(0, Rational(-1));
    m.row_add(1, Rational(-2));
    m.emit_eq(Rational(-4));
    m.row_add(0, Rational(1));
    m.row_add(1, Rational(-1));
    m.emit_le(Rational(1));
    line::lp::LpSolution<Rational> s = simplex_solve(m);
    REQUIRE(s.status == LpStatus::Optimal);
    CHECK(s.objective == Rational(2));
    CHECK(s.x[0] == Rational(0));
    CHECK(s.x[1] == Rational(2));
}

TEST_CASE("variable bounds are enforced without caller-supplied rows") {
    // max x  s.t. x <= 10 by bound only, plus a redundant row x + y <= 100.
    LpModel<Rational> m(2);
    m.set_cost(0, Rational(1));
    m.set_bounds(0, Rational(2), Rational(10));
    m.set_bounds(1, Rational(0), Rational(1));
    m.row_add(0, Rational(1));
    m.row_add(1, Rational(1));
    m.emit_le(Rational(100));
    line::lp::LpSolution<Rational> s = simplex_solve(m);
    REQUIRE(s.status == LpStatus::Optimal);
    CHECK(s.objective == Rational(10));
    CHECK(s.x[0] == Rational(10));
    // Minimizing hits the lower bound, which is not the origin.
    m.set_maximize(false);
    line::lp::LpSolution<Rational> s2 = simplex_solve(m);
    REQUIRE(s2.status == LpStatus::Optimal);
    CHECK(s2.objective == Rational(2));
}

TEST_CASE("fixed variables are substituted out") {
    // max x + y with y fixed at 3/2 and x <= 1.
    LpModel<Rational> m(2);
    m.set_cost(0, Rational(1));
    m.set_cost(1, Rational(1));
    m.set_bounds(0, Rational(0), Rational(1));
    m.fix(1, rat(3, 2));
    m.row_add(0, Rational(1));
    m.row_add(1, Rational(1));
    m.emit_le(Rational(10));
    line::lp::LpSolution<Rational> s = simplex_solve(m);
    REQUIRE(s.status == LpStatus::Optimal);
    CHECK(s.x[1] == rat(3, 2));
    CHECK(s.objective == rat(5, 2));
}

TEST_CASE("free variables are split into a difference of nonnegatives") {
    // min x s.t. x >= -4, expressed with a free variable and a >= row.
    LpModel<Rational> m(1);
    m.set_maximize(false);
    m.set_cost(0, Rational(1));
    m.set_free(0);
    m.row_add(0, Rational(1));
    m.emit_ge(Rational(-4));
    line::lp::LpSolution<Rational> s = simplex_solve(m);
    REQUIRE(s.status == LpStatus::Optimal);
    CHECK(s.x[0] == Rational(-4));
    CHECK(s.objective == Rational(-4));
}

TEST_CASE("row accumulation sums repeated entries") {
    // Emitting 2x through two add() calls must be the same LP as one call.
    LpModel<Rational> m(1);
    m.set_cost(0, Rational(1));
    m.row_add(0, Rational(1));
    m.row_add(0, Rational(1));
    m.emit_le(Rational(6));
    line::lp::LpSolution<Rational> s = simplex_solve(m);
    REQUIRE(s.status == LpStatus::Optimal);
    CHECK(s.objective == Rational(3));
    CHECK(m.num_nonzeros() == 1u);
}

TEST_CASE("degenerate transportation-style LP with redundant equalities") {
    // Supplies 1,1 and demands 1,1 with the redundant total-balance row added
    // on purpose: the extra equality is dependent and must be dropped as a
    // redundant phase-1 row rather than reported as infeasible.
    //   min 4 x11 + 1 x12 + 2 x21 + 3 x22
    // Optimal assignment x12 = x21 = 1, cost 3.
    LpModel<Rational> m(4);
    m.set_maximize(false);
    const long c[4] = {4, 1, 2, 3};
    for (int j = 0; j < 4; ++j) m.set_cost(j, line::num_traits<Rational>::from_int(c[j]));
    m.row_add(0, Rational(1));
    m.row_add(1, Rational(1));
    m.emit_eq(Rational(1));  // supply 1
    m.row_add(2, Rational(1));
    m.row_add(3, Rational(1));
    m.emit_eq(Rational(1));  // supply 2
    m.row_add(0, Rational(1));
    m.row_add(2, Rational(1));
    m.emit_eq(Rational(1));  // demand 1
    m.row_add(1, Rational(1));
    m.row_add(3, Rational(1));
    m.emit_eq(Rational(1));  // demand 2
    for (int j = 0; j < 4; ++j) m.row_add(j, Rational(1));
    m.emit_eq(Rational(2));  // redundant
    line::lp::LpSolution<Rational> s = simplex_solve(m);
    REQUIRE(s.status == LpStatus::Optimal);
    CHECK(s.objective == Rational(3));
    CHECK(s.x[1] == Rational(1));
    CHECK(s.x[2] == Rational(1));
}

TEST_CASE("exact and double solvers agree on a randomly generated feasible LP") {
    // Same data at both types; the double optimum must match the rational one.
    const int n = 6, mrows = 4;
    LpModel<Rational> mq(n);
    LpModel<double> md(n);
    unsigned seed = 12345u;
    for (int j = 0; j < n; ++j) {
        seed = seed * 1103515245u + 12345u;
        const long cj = static_cast<long>((seed >> 16) % 9) + 1;
        mq.set_cost(j, line::num_traits<Rational>::from_int(cj));
        md.set_cost(j, static_cast<double>(cj));
        mq.set_bounds(j, Rational(0), Rational(1));
        md.set_bounds(j, 0.0, 1.0);
    }
    for (int i = 0; i < mrows; ++i) {
        for (int j = 0; j < n; ++j) {
            seed = seed * 1103515245u + 12345u;
            const long a = static_cast<long>((seed >> 16) % 5) + 1;
            mq.row_add(j, line::num_traits<Rational>::from_int(a));
            md.row_add(j, static_cast<double>(a));
        }
        mq.emit_le(line::num_traits<Rational>::from_int(7));
        md.emit_le(7.0);
    }
    line::lp::LpSolution<Rational> sq = simplex_solve(mq);
    line::lp::LpSolution<double> sd = simplex_solve(md);
    REQUIRE(sq.status == LpStatus::Optimal);
    REQUIRE(sd.status == LpStatus::Optimal);
    CHECK(sd.objective ==
          doctest::Approx(line::num_traits<Rational>::to_double(sq.objective)).epsilon(1e-12));
}
