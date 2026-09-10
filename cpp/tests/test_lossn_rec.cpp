/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * `lossn_rec`: the exact normalising constant of a loss network by MDD-rec.
 *
 * WHY THIS METHOD EXISTS. The admissible set of a Kelly loss network is
 * {n >= 0 : A n <= C} and the stationary law is independent Poisson counts
 * truncated to it, so the normalising constant is a sum of a product form over a
 * set that is finite and bounded per coordinate -- exactly what a decision
 * diagram holds and `mdd_rec` walks. It complements the Manjunath-Sikdar residue
 * transform, which is equally exact but whose residue argument COUNTS WHOLE
 * UNITS and so needs an integral A and C, and this port's `lossn_erlangfp`,
 * which raises (1-E_i) to an unsigned integer power and needs integrality too.
 * A FRACTIONAL region therefore had no exact route here at all before MDD-rec,
 * only the Monte Carlo `lossn_mci` whose answer is a random variable.
 *
 * THREE ORACLES: a brute-force sum over the admissible set that shares no code
 * path with the walk; `lossn_manjunath`, which must agree on an INTEGRAL region;
 * and the other codebases, pinned at 12 decimals.
 *
 * ARITHMETIC. Everything but the reported `lG` is a sum, a product and one
 * factorial, so the method runs under exact arithmetic and the last case checks
 * that the rational G is the rational sum, exactly.
 */

#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/lossn/lossn_manjunath.h"
#include "line/api/lossn/lossn_rec.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

using namespace line;
namespace lossn = line::lossn;

namespace {

/** 2 links, 3 routes: routes 0 and 1 take one link each, route 2 takes both. */
Matrix<double> a_int() {
    Matrix<double> A(2, 3, 0.0);
    A(0, 0) = 1.0;
    A(0, 2) = 1.0;
    A(1, 1) = 1.0;
    A(1, 2) = 1.0;
    return A;
}

/** One link with FRACTIONAL class sizes and capacity. */
Matrix<double> a_frac() {
    Matrix<double> A(1, 3, 0.0);
    A(0, 0) = 1.5;
    A(0, 1) = 0.75;
    A(0, 2) = 2.25;
    return A;
}

/** MATLAB, native python and the JAR, at %.12f. */
const double INT_QLEN[3] = {2.282838019229, 1.619031944010, 0.994469884949};
const double INT_LOSS[3] = {0.086864792309, 0.100537808883, 0.171275095875};
const double INT_LG = 5.342866909910;
const double FRAC_QLEN[3] = {1.383844183191, 2.545477546770, 0.537948865597};
const double FRAC_LOSS[3] = {0.308077908404, 0.151507484410, 0.462051134403};
const double FRAC_LG = 5.451139581448;

bool feasible(const std::vector<int>& n, const Matrix<double>& A, const std::vector<double>& C,
              int plus) {
    for (std::size_t j = 0; j < C.size(); ++j) {
        double s = 0;
        for (std::size_t r = 0; r < n.size(); ++r)
            s += A(j, r) * (n[r] + (static_cast<int>(r) == plus ? 1 : 0));
        if (s > C[j] + 1e-12) return false;
    }
    return true;
}

/** The carried load, the blocking and log G, summed state by state. */
void brute(const std::vector<double>& nu, const Matrix<double>& A, const std::vector<double>& C,
           int cap, std::vector<double>& qlen, std::vector<double>& loss, double& lg) {
    const std::size_t K = nu.size();
    std::vector<int> n(K, 0);
    double tot = 0;
    std::vector<double> num(K, 0.0), acc(K, 0.0);
    long total = 1;
    for (std::size_t r = 0; r < K; ++r) total *= cap;
    for (long code = 0; code < total; ++code) {
        long c = code;
        for (std::size_t r = 0; r < K; ++r) {
            n[r] = static_cast<int>(c % cap);
            c /= cap;
        }
        if (!feasible(n, A, C, -1)) continue;
        double w = 1;
        for (std::size_t r = 0; r < K; ++r) {
            double f = 1;
            for (int k = 2; k <= n[r]; ++k) f *= k;
            w *= std::pow(nu[r], n[r]) / f;
        }
        tot += w;
        for (std::size_t r = 0; r < K; ++r) {
            num[r] += w * n[r];
            if (feasible(n, A, C, static_cast<int>(r))) acc[r] += w;
        }
    }
    qlen.assign(K, 0.0);
    loss.assign(K, 0.0);
    for (std::size_t r = 0; r < K; ++r) {
        qlen[r] = num[r] / tot;
        loss[r] = 1 - acc[r] / tot;
    }
    lg = std::log(tot);
}

}  // namespace

TEST_CASE("MDD-rec matches the brute-force sum on an integral region") {
    const std::vector<double> nu = {2.5, 1.8, 1.2};
    const std::vector<double> C = {6.0, 5.0};
    const lossn::LossnRecResult<double> r = lossn::lossn_rec<double>(nu, a_int(), C);
    std::vector<double> bq, bl;
    double blg = 0;
    brute(nu, a_int(), C, 8, bq, bl, blg);
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(r.QLen[k] == doctest::Approx(bq[k]).epsilon(1e-11));
        CHECK(r.Loss[k] == doctest::Approx(bl[k]).epsilon(1e-11));
    }
    CHECK(r.lG == doctest::Approx(blg).epsilon(1e-11));
    // one walk for G, one per class for the blocking ratios
    CHECK(r.iterations == 4);
}

TEST_CASE("MDD-rec agrees with the residue transform where both apply") {
    const std::vector<double> nu = {2.5, 1.8, 1.2};
    const std::vector<double> C = {6.0, 5.0};
    const lossn::LossnRecResult<double> r = lossn::lossn_rec<double>(nu, a_int(), C);
    const lossn::LossnManjunathResult<double> m =
        lossn::lossn_manjunath<double>(nu, a_int(), C, lossn::LossnManjunathOptions());
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(r.QLen[k] == doctest::Approx(m.QLen[k]).epsilon(1e-10));
        CHECK(r.Loss[k] == doctest::Approx(m.Loss[k]).epsilon(1e-10));
    }
    CHECK(r.lG == doctest::Approx(m.lG).epsilon(1e-10));
}

TEST_CASE("MDD-rec is exact where the residue transform cannot count") {
    const std::vector<double> nu = {2.0, 3.0, 1.0};
    const std::vector<double> C = {7.5};
    const lossn::LossnRecResult<double> r = lossn::lossn_rec<double>(nu, a_frac(), C);
    std::vector<double> bq, bl;
    double blg = 0;
    brute(nu, a_frac(), C, 12, bq, bl, blg);   // class 1 admits 10 calls: a cap of 10 truncates the set
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(r.QLen[k] == doctest::Approx(bq[k]).epsilon(1e-11));
        CHECK(r.Loss[k] == doctest::Approx(bl[k]).epsilon(1e-11));
    }
    CHECK(r.lG == doctest::Approx(blg).epsilon(1e-11));
}

TEST_CASE("the loss-network values are pinned across the codebases") {
    {
        const std::vector<double> nu = {2.5, 1.8, 1.2};
        const std::vector<double> C = {6.0, 5.0};
        const lossn::LossnRecResult<double> r = lossn::lossn_rec<double>(nu, a_int(), C);
        for (std::size_t k = 0; k < 3; ++k) {
            CHECK(r.QLen[k] == doctest::Approx(INT_QLEN[k]).epsilon(1e-10));
            CHECK(r.Loss[k] == doctest::Approx(INT_LOSS[k]).epsilon(1e-10));
        }
        CHECK(r.lG == doctest::Approx(INT_LG).epsilon(1e-10));
    }
    {
        const std::vector<double> nu = {2.0, 3.0, 1.0};
        const std::vector<double> C = {7.5};
        const lossn::LossnRecResult<double> r = lossn::lossn_rec<double>(nu, a_frac(), C);
        for (std::size_t k = 0; k < 3; ++k) {
            CHECK(r.QLen[k] == doctest::Approx(FRAC_QLEN[k]).epsilon(1e-10));
            CHECK(r.Loss[k] == doctest::Approx(FRAC_LOSS[k]).epsilon(1e-10));
        }
        CHECK(r.lG == doctest::Approx(FRAC_LG).epsilon(1e-10));
    }
}

TEST_CASE("a class consuming no resource is refused by name") {
    const std::vector<double> nu = {1.0, 1.0};
    Matrix<double> A(1, 2, 0.0);
    A(0, 0) = 1.0;
    const std::vector<double> C = {3.0};
    CHECK_THROWS_AS(lossn::lossn_rec<double>(nu, A, C), InputError);
}

TEST_CASE("MDD-rec runs in exact arithmetic on a loss network") {
    // The whole walk is a sum, a product and one factorial, so G is a RATIONAL
    // and must equal the rational sum over the admissible set exactly, with no
    // tolerance at all.
    const std::vector<Rational> nu = {Rational(5, 2), Rational(9, 5)};
    Matrix<Rational> A(1, 2, Rational(0));
    A(0, 0) = Rational(1);
    A(0, 1) = Rational(2);
    const std::vector<Rational> C = {Rational(4)};
    const lossn::LossnRecResult<Rational> r = lossn::lossn_rec<Rational>(nu, A, C);
    Rational G(0);
    for (int n0 = 0; n0 <= 4; ++n0)
        for (int n1 = 0; n1 <= 2; ++n1) {
            if (n0 + 2 * n1 > 4) continue;
            Rational w(1);
            for (int k = 1; k <= n0; ++k) w = w * nu[0] / Rational(k);
            for (int k = 1; k <= n1; ++k) w = w * nu[1] / Rational(k);
            G = G + w;
        }
    CHECK(r.G == G);
}
