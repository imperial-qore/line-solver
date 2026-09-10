/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Joint moment conversions and the house-of-moments tables. The oracles are the
 * algebraic identities: every edge pair is a mutual inverse, the separable
 * edges agree with the univariate table applied fibre by fibre, and the
 * cumulant edges invert each other. The exact backend makes the round trips the
 * IDENTITY rather than merely close, which is the point of running them here.
 *
 * The numeric rows tagged MATLAB are values printed by the reference
 * implementation in matlab/src/api/moment.
 */
#include <vector>

#include "doctest.h"
#include "line/api/moment/moment_cumulant.h"
#include "line/api/moment/moment_housematrix.h"
#include "line/api/moment/moment_joint.h"
#include "line/api/moment/moment_tail.h"
#include "line/api/moment/moment_tensor.h"

using line::Matrix;
using line::Rational;
using namespace line::moment;

namespace {

Rational rq(long n, long d = 1) { return Rational(n) / Rational(d); }

/** A nontrivial 3x4 joint array with no symmetry, entries 1,...,12. */
MomentTensor<Rational> ramp34() {
    std::vector<std::size_t> sz;
    sz.push_back(3);
    sz.push_back(4);
    MomentTensor<Rational> A(sz);
    for (std::size_t i = 0; i < 12; ++i) A.data[i] = rq(static_cast<long>(i) + 1);
    return A;
}

bool same(const MomentTensor<Rational>& a, const MomentTensor<Rational>& b) {
    if (a.sz != b.sz) return false;
    for (std::size_t i = 0; i < a.data.size(); ++i)
        if (a.data[i] != b.data[i]) return false;
    return true;
}

}  // namespace

TEST_CASE("moment_housematrix known tables") {
    // factorial_from_raw is the signed Stirling table: rows 0..3 of s(i,k).
    Matrix<Rational> s = moment_housematrix<Rational>(MomentEdge::FactorialFromRaw, 3);
    CHECK(s(3, 1) == rq(2));
    CHECK(s(3, 2) == rq(-3));
    CHECK(s(3, 3) == rq(1));
    // binomial_from_factorial divides by i!.
    Matrix<Rational> b = moment_housematrix<Rational>(MomentEdge::BinomialFromFactorial, 4);
    CHECK(b(4, 4) == rq(1, 24));
    // the tail edges are upper triangular, with T(1,k) = C(k-1,0) = 1.
    Matrix<Rational> t = moment_housematrix<Rational>(MomentEdge::BinomialFromTail, 4);
    CHECK(t(1, 4) == rq(1));
    CHECK(t(2, 4) == rq(3));
    CHECK(t(4, 1) == rq(0));
}

TEST_CASE("moment_housematrix edge pairs are mutual inverses") {
    const int n = 6;
    const MomentEdge fwd[7] = {MomentEdge::FactorialFromRaw,
                               MomentEdge::UpfactorialFromRaw,
                               MomentEdge::BinomialFromFactorial,
                               MomentEdge::NegbinomialFromUpfactorial,
                               MomentEdge::UpfactorialFromFactorial,
                               MomentEdge::NegbinomialFromBinomial,
                               MomentEdge::BinomialFromTail};
    const MomentEdge bwd[7] = {MomentEdge::RawFromFactorial,
                               MomentEdge::RawFromUpfactorial,
                               MomentEdge::FactorialFromBinomial,
                               MomentEdge::UpfactorialFromNegbinomial,
                               MomentEdge::FactorialFromUpfactorial,
                               MomentEdge::BinomialFromNegbinomial,
                               MomentEdge::TailFromBinomial};
    for (int e = 0; e < 7; ++e) {
        Matrix<Rational> A = moment_housematrix<Rational>(fwd[e], n);
        Matrix<Rational> B = moment_housematrix<Rational>(bwd[e], n);
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n; ++j) {
                Rational acc = rq(0);
                for (int k = 0; k <= n; ++k) acc += B(i, k) * A(k, j);
                CHECK(acc == rq(i == j ? 1 : 0));
            }
    }
}

TEST_CASE("moment_cumulant_from_raw matches the known central relations") {
    // For a law with raw moments 1, m1, m2, m3 the cumulants are m1, m2-m1^2 and
    // m3 - 3 m1 m2 + 2 m1^3.
    std::vector<Rational> m;
    m.push_back(rq(1));
    m.push_back(rq(2));
    m.push_back(rq(7));
    m.push_back(rq(35));
    std::vector<Rational> k = moment_cumulant_from_raw<Rational>(m);
    CHECK(k[0] == rq(0));
    CHECK(k[1] == rq(2));
    CHECK(k[2] == rq(3));
    CHECK(k[3] == rq(9));
    std::vector<Rational> back = moment_raw_from_cumulant<Rational>(k);
    for (std::size_t i = 0; i < m.size(); ++i) CHECK(back[i] == m[i]);
}

TEST_CASE("moment_tail_from_binomial inverts moment_binomial_from_tail") {
    std::vector<Rational> t;
    for (long i = 0; i < 6; ++i) t.push_back(rq(6 - i, i + 1));
    t[0] = rq(1);
    std::vector<Rational> b = moment_binomial_from_tail<Rational>(t);
    std::vector<Rational> back = moment_tail_from_binomial<Rational>(b);
    for (std::size_t i = 0; i < t.size(); ++i) CHECK(back[i] == t[i]);
}

TEST_CASE("moment_tensortrans is the mode product") {
    MomentTensor<Rational> A = ramp34();
    // Transforming by the identity leaves A alone, on either mode.
    Matrix<Rational> I3(3, 3, rq(0));
    for (std::size_t i = 0; i < 3; ++i) I3(i, i) = rq(1);
    CHECK(same(moment_tensortrans<Rational>(A, I3, 0), A));
    // A fibre along mode 0 at column j is A(:,j); scaling row 2 by 10 must hit
    // exactly the entries with first order 2.
    Matrix<Rational> S = I3;
    S(2, 2) = rq(10);
    MomentTensor<Rational> B = moment_tensortrans<Rational>(A, S, 0);
    for (std::size_t j = 0; j < 4; ++j) {
        std::vector<std::size_t> o;
        o.push_back(2);
        o.push_back(j);
        CHECK(B.at(o) == A.at(o) * rq(10));
    }
}

TEST_CASE("separable joint edges invert each other") {
    MomentTensor<Rational> A = ramp34();
    CHECK(same(moment_joint_raw_from_factorial<Rational>(
                   moment_joint_factorial_from_raw<Rational>(A)),
               A));
    CHECK(same(moment_joint_factorial_from_binomial<Rational>(
                   moment_joint_binomial_from_factorial<Rational>(A)),
               A));
    CHECK(same(moment_joint_raw_from_upfactorial<Rational>(
                   moment_joint_upfactorial_from_raw<Rational>(A)),
               A));
    CHECK(same(moment_joint_upfactorial_from_negbinomial<Rational>(
                   moment_joint_negbinomial_from_upfactorial<Rational>(A)),
               A));
    CHECK(same(moment_joint_factorial_from_upfactorial<Rational>(
                   moment_joint_upfactorial_from_factorial<Rational>(A)),
               A));
    CHECK(same(moment_joint_binomial_from_negbinomial<Rational>(
                   moment_joint_negbinomial_from_binomial<Rational>(A)),
               A));
    CHECK(same(moment_joint_binomial_from_tail<Rational>(
                   moment_joint_tail_from_binomial<Rational>(A)),
               A));
}

TEST_CASE("joint cumulant edges invert each other") {
    MomentTensor<Rational> A = ramp34();
    A.data[0] = rq(1);  // m_0 must be 1 for a moment array
    MomentTensor<Rational> k = moment_joint_cumulant_from_raw<Rational>(A);
    CHECK(k.data[0] == rq(0));
    CHECK(same(moment_joint_raw_from_cumulant<Rational>(k), A));
    // The factorial-cumulant pair is the same recurrence on the factorial sequence.
    CHECK(same(moment_joint_factcumulant_from_factorial<Rational>(A), k));
    CHECK(same(moment_joint_factorial_from_factcumulant<Rational>(k), A));
}

TEST_CASE("joint cumulants reduce to the univariate ones on a degenerate class") {
    // A d = 2 array whose second class is degenerate carries the univariate
    // cumulants in its first column.
    std::vector<Rational> m;
    m.push_back(rq(1));
    m.push_back(rq(2));
    m.push_back(rq(7));
    m.push_back(rq(35));
    std::vector<std::size_t> sz;
    sz.push_back(4);
    sz.push_back(1);
    MomentTensor<Rational> A(sz);
    for (std::size_t i = 0; i < 4; ++i) A.data[i] = m[i];
    MomentTensor<Rational> k = moment_joint_cumulant_from_raw<Rational>(A);
    std::vector<Rational> ku = moment_cumulant_from_raw<Rational>(m);
    for (std::size_t i = 0; i < 4; ++i) CHECK(k.data[i] == ku[i]);
}

TEST_CASE("joint central moments about the mean") {
    MomentTensor<Rational> A = ramp34();
    A.data[0] = rq(1);
    MomentTensor<Rational> mc = moment_joint_central_from_raw<Rational>(A);
    // The first-order central moments vanish in every class.
    std::vector<std::size_t> e0, e1;
    e0.push_back(1);
    e0.push_back(0);
    e1.push_back(0);
    e1.push_back(1);
    CHECK(mc.at(e0) == rq(0));
    CHECK(mc.at(e1) == rq(0));
    // The (1,1) entry is the covariance m11 - m10 m01.
    std::vector<std::size_t> e11;
    e11.push_back(1);
    e11.push_back(1);
    CHECK(mc.at(e11) == A.at(e11) - A.at(e0) * A.at(e1));
}

TEST_CASE("moment_joint_marking and moment_joint_aggregate are inverse directions") {
    // Binomial marking of a total count: p = [3/10, 7/10], orders [2,2].
    std::vector<Rational> f;
    f.push_back(rq(1));
    f.push_back(rq(4));
    f.push_back(rq(18));
    f.push_back(rq(96));
    f.push_back(rq(600));
    std::vector<Rational> p;
    p.push_back(rq(3, 10));
    p.push_back(rq(7, 10));
    std::vector<std::size_t> dims;
    dims.push_back(2);
    dims.push_back(2);
    MomentTensor<Rational> F = moment_joint_marking<Rational>(f, p, dims);
    std::vector<std::size_t> o;
    o.push_back(1);
    o.push_back(1);
    CHECK(F.at(o) == rq(3, 10) * rq(7, 10) * f[2]);
    // Aggregating recovers f up to the smallest per-class order.
    std::vector<Rational> fa = moment_joint_aggregate<Rational>(F);
    CHECK(fa.size() == 3u);
    for (std::size_t i = 0; i < fa.size(); ++i) CHECK(fa[i] == f[i]);
}

TEST_CASE("joint conversions match the MATLAB reference on a 3x4 ramp") {
    // A = reshape(1:12,3,4) in matlab/src/api/moment; the rows below are the
    // mat2str output of the reference implementation, read column major.
    std::vector<std::size_t> sz;
    sz.push_back(3);
    sz.push_back(4);
    MomentTensor<double> A(sz);
    for (std::size_t i = 0; i < 12; ++i) A.data[i] = static_cast<double>(i + 1);
    const double cum[12] = {0, 2, -1, 4, -3, 6, -9, 18, -54, 54, -162, 648};
    const double cen[12] = {1, 0, -1, 0, -3, 6, -9, 18, -27, 54, -81, 108};
    const double upf[12] = {1, 2, 5, 4, 5, 11, 11, 13, 28, 39, 45, 96};
    const double bft[12] = {1, 5, 3, 21, 51, 27, 27, 63, 33, 10, 23, 12};
    const double cft[12] = {1, 0, -14, 0, -54, 414, -366, 2070, -10746, 14040, -59616, 260496};
    MomentTensor<double> kc = moment_joint_cumulant_from_raw<double>(A);
    MomentTensor<double> mc = moment_joint_central_from_raw<double>(A);
    MomentTensor<double> fp = moment_joint_upfactorial_from_raw<double>(A);
    MomentTensor<double> bt = moment_joint_binomial_from_tail<double>(A);
    MomentTensor<double> ct = moment_joint_central_from_tail<double>(A);
    for (std::size_t i = 0; i < 12; ++i) {
        CHECK(kc.data[i] == doctest::Approx(cum[i]));
        CHECK(mc.data[i] == doctest::Approx(cen[i]));
        CHECK(fp.data[i] == doctest::Approx(upf[i]));
        CHECK(bt.data[i] == doctest::Approx(bft[i]));
        CHECK(ct.data[i] == doctest::Approx(cft[i]));
    }
    // moment_joint_marking(f, [0.3 0.7], [2 2]) and its aggregate.
    std::vector<double> f;
    f.push_back(1);
    f.push_back(4);
    f.push_back(18);
    f.push_back(96);
    f.push_back(600);
    std::vector<double> p;
    p.push_back(0.3);
    p.push_back(0.7);
    std::vector<std::size_t> dims;
    dims.push_back(2);
    dims.push_back(2);
    MomentTensor<double> F = moment_joint_marking<double>(f, p, dims);
    const double mrk[9] = {1, 1.2, 1.62, 2.8, 3.78, 6.048, 8.82, 14.112, 26.46};
    for (std::size_t i = 0; i < 9; ++i) CHECK(F.data[i] == doctest::Approx(mrk[i]));
    std::vector<double> ag = moment_joint_aggregate<double>(F);
    CHECK(ag.size() == 3u);
    CHECK(ag[1] == doctest::Approx(4.0));
    CHECK(ag[2] == doctest::Approx(18.0));
    // moment_tail_from_binomial([1 2 7 35]) = [1 30 -63 35].
    std::vector<double> m;
    m.push_back(1);
    m.push_back(2);
    m.push_back(7);
    m.push_back(35);
    std::vector<double> tb = moment_tail_from_binomial<double>(m);
    CHECK(tb[1] == doctest::Approx(30.0));
    CHECK(tb[2] == doctest::Approx(-63.0));
    CHECK(tb[3] == doctest::Approx(35.0));
}

TEST_CASE("moment_joint_central_from_tail agrees with the composed route") {
    MomentTensor<Rational> t = ramp34();
    t.data[0] = rq(1);
    MomentTensor<Rational> b = moment_joint_binomial_from_tail<Rational>(t);
    MomentTensor<Rational> f = moment_joint_factorial_from_binomial<Rational>(b);
    MomentTensor<Rational> expected =
        moment_joint_central_from_raw<Rational>(moment_joint_raw_from_factorial<Rational>(f));
    CHECK(same(moment_joint_central_from_tail<Rational>(t), expected));
}
