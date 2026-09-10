/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Boundary marshalling for host bindings. The round trips are the point: a
 * MATLAB gateway hands over a column-major buffer and expects one back, and a
 * transposed round trip is the classic way to lose a non-square model silently.
 */
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/io/marshal.h"

using line::Matrix;
using line::Rational;
using namespace line::io;

TEST_CASE("column-major ingest matches the MATLAB layout on a non-square matrix") {
    // MATLAB [1 2 3; 4 5 6] is stored as 1 4 2 5 3 6.
    const double buf[6] = {1, 4, 2, 5, 3, 6};
    Matrix<double> m = from_column_major<double>(buf, 2, 3);
    CHECK(m.rows() == 2);
    CHECK(m.cols() == 3);
    CHECK(m(0, 0) == 1.0);
    CHECK(m(0, 2) == 3.0);
    CHECK(m(1, 0) == 4.0);
    CHECK(m(1, 2) == 6.0);

    double out[6] = {0};
    to_column_major(m, out);
    for (int i = 0; i < 6; ++i) CHECK(out[i] == buf[i]);
}

TEST_CASE("row-major ingest is the numpy default order and differs from column-major") {
    const double buf[6] = {1, 2, 3, 4, 5, 6};
    Matrix<double> r = from_row_major<double>(buf, 2, 3);
    Matrix<double> c = from_column_major<double>(buf, 2, 3);
    CHECK(r(0, 1) == 2.0);
    CHECK(c(0, 1) == 3.0);  // the two layouts genuinely disagree
}

TEST_CASE("ingest converts into the exact instantiation without rounding") {
    const double buf[2] = {0.5, 0.25};  // both dyadic, so exactly representable
    Matrix<Rational> m = from_column_major<Rational>(buf, 2, 1);
    CHECK(m(0, 0) == Rational(1, 2));
    CHECK(m(1, 0) == Rational(1, 4));
}

TEST_CASE("an exact result crosses as numerator, denominator and an approximation") {
    Matrix<Rational> L(1, 1);
    L(0, 0) = Rational(1, 2);
    auto r = line::pfqn::pfqn_ca(L, std::vector<int>{3}, Matrix<Rational>());
    ExactValue e = marshal_scalar(r.G);
    CHECK(e.is_exact);
    CHECK(e.numerator == "1");
    CHECK(e.denominator == "8");
    CHECK(e.approx == doctest::Approx(0.125).epsilon(1e-15));

    // The double instantiation reports itself as inexact and carries no parts.
    ExactValue d = marshal_scalar(0.125);
    CHECK_FALSE(d.is_exact);
    CHECK(d.numerator.empty());
    CHECK(d.approx == 0.125);
}

TEST_CASE("error identifiers are stable and distinguish the three failure kinds") {
    CHECK(std::string(error_id(line::InputError("x"))) == "line:input");
    CHECK(std::string(error_id(line::NumericError("x"))) == "line:numeric");
    CHECK(std::string(error_id(line::UnsupportedError("x"))) == "line:unsupported");
    CHECK(std::string(error_id(line::Error("x"))) == "line:error");
}

TEST_CASE("null buffers are refused rather than dereferenced") {
    CHECK_THROWS_AS(from_column_major<double>(static_cast<const double*>(nullptr), 2, 2),
                    line::InputError);
    Matrix<double> m(2, 2, 1.0);
    CHECK_THROWS_AS(to_column_major(m, nullptr), line::InputError);
}
