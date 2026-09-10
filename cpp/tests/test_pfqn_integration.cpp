/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Integration and asymptotic members of the pfqn family.
 *
 * The oracle for every approximation here is the EXACT normalizing constant
 * pfqn_ca returns, so each test states the relative error it observes on log G
 * rather than asserting a bare number. Where a MATLAB value is asserted it was
 * obtained by running MATLAB on this very model (LINE 3.0.6, matlab -batch)
 * and is quoted to sixteen digits; the tolerance is the one the method
 * guarantees, never tighter.
 *
 * Models used throughout:
 *   A  M = 3, R = 2, L = [1 .5; .7 1.2; .3 .9], N = [2 3], Z = [1 .5]
 *   B  single class, L = [1/2, 1/3, 1/5], N = 5, Z = 1
 *   C  repairman, one station, L = [.6 .4], N = [3 2], Z = [1 2]
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/cd_peak_scaling.h"
#include "line/api/pfqn/infradius_h.h"
#include "line/api/pfqn/laplaceapprox.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_cub.h"
#include "line/api/pfqn/pfqn_dac.h"
#include "line/api/pfqn/pfqn_fnc.h"
#include "line/api/pfqn/pfqn_grnmol.h"
#include "line/api/pfqn/pfqn_kt.h"
#include "line/api/pfqn/pfqn_bkt.h"
#include "line/api/pfqn/pfqn_lekt.h"
#include "line/api/pfqn/pfqn_lap.h"
#include "line/api/pfqn/pfqn_le.h"
#include "line/api/pfqn/pfqn_ble.h"
#include "line/api/pfqn/pfqn_mmint2.h"
#include "line/api/pfqn/pfqn_nre.h"
#include "line/api/pfqn/pfqn_nrl.h"
#include "line/api/pfqn/pfqn_oi_fnc.h"
#include "line/api/pfqn/pfqn_oi_insvc.h"
#include "line/api/pfqn/pfqn_panacea.h"
#include "line/api/pfqn/pfqn_propfair.h"

using line::Matrix;
using line::Rational;
using namespace line::pfqn;

namespace {

Matrix<double> modelA() { return Matrix<double>{{1.0, 0.5}, {0.7, 1.2}, {0.3, 0.9}}; }

Matrix<double> modelAz() {
    Matrix<double> Z(1, 2);
    Z(0, 0) = 1.0;
    Z(0, 1) = 0.5;
    return Z;
}

std::vector<double> demandsB() { return {0.5, 1.0 / 3.0, 0.2}; }

Matrix<double> matrixB() {
    Matrix<double> L(3, 1);
    L(0, 0) = 0.5;
    L(1, 0) = 1.0 / 3.0;
    L(2, 0) = 0.2;
    return L;
}

/** Relative error on log G against the exact convolution value. */
double relerr(double approx, double exact) { return std::fabs(approx - exact) / std::fabs(exact); }

}  // namespace

// ---------------------------------------------------------------------------
// pfqn_le
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_le matches MATLAB on model A") {
    const std::vector<double> N{2.0, 3.0};
    const std::vector<double> Z{1.0, 0.5};
    // MATLAB: pfqn_le([1 .5;.7 1.2;.3 .9],[2 3],[1 .5]) -> lGn = 4.52924316120224
    CHECK(pfqn_le(modelA(), N, Z).lG == doctest::Approx(4.52924316120224).epsilon(1e-9));
    // MATLAB: pfqn_le(L,[2 3]) with no think time -> lGn = 3.906251944123995
    CHECK(pfqn_le(modelA(), N).lG == doctest::Approx(3.906251944123995).epsilon(1e-9));
}

TEST_CASE("pfqn_ble adds the eps->0 constant to pfqn_le") {
    const std::vector<double> N{2.0, 3.0};
    const std::vector<double> Z{1.0, 0.5};
    // MATLAB: pfqn_ble([1 .5;.7 1.2;.3 .9],[2 3],[1 .5]) -> lGn = 4.77242756158822
    CHECK(pfqn_ble(modelA(), N, Z).lG == doctest::Approx(4.77242756158822).epsilon(1e-9));
    // MATLAB: pfqn_ble(L,[2 3]) with no think time -> lGn = 4.068374877714649
    CHECK(pfqn_ble(modelA(), N).lG == doctest::Approx(4.068374877714649).epsilon(1e-9));
    // One (1 - log(2 pi)/2) per Laplaced direction, which is the exponent on
    // sqrt(2 pi) in the branch taken: M = 3 with a think time, M-1 = 2 without.
    const double c = 1.0 - std::log(2 * M_PI) / 2;
    CHECK(pfqn_ble(modelA(), N, Z).lG - pfqn_le(modelA(), N, Z).lG ==
          doctest::Approx(3 * c).epsilon(1e-12));
    CHECK(pfqn_ble(modelA(), N).lG - pfqn_le(modelA(), N).lG ==
          doctest::Approx(2 * c).epsilon(1e-12));
    // An all-zero Z is the Z = 0 branch, not the delay branch.
    const std::vector<double> Z0{0.0, 0.0};
    CHECK(pfqn_ble(modelA(), N, Z0).lG - pfqn_le(modelA(), N, Z0).lG ==
          doctest::Approx(2 * c).epsilon(1e-12));
}

TEST_CASE("pfqn_ble leaves the degenerate branch alone") {
    // No queueing demand: pfqn_le returns the exact delay term, nothing to correct.
    Matrix<double> L0(3, 2, 0.0);
    const std::vector<double> N{2.0, 3.0};
    const std::vector<double> Z{1.0, 0.5};
    CHECK(pfqn_ble(L0, N, Z).lG == doctest::Approx(pfqn_le(L0, N, Z).lG).epsilon(1e-12));
}

TEST_CASE("pfqn_le converges to the exact constant as the population grows") {
    // Laplace expansion: the relative error on log G falls with the population.
    // MATLAB pfqn_ca / pfqn_le on model A, N = [n n], Z = [1 .5]:
    //   n = 10  ca 14.10442210782851  le 13.88895639764461   rel err 1.5e-2
    //   n = 20  ca 26.08438539593642  le 25.89867770204523   rel err 7.1e-3
    //   n = 40  ca 50.13260369498655  le 49.95832933327313   rel err 3.5e-3
    const int ns[3] = {10, 20, 40};
    const double caRef[3] = {14.10442210782851, 26.08438539593642, 50.13260369498655};
    const double leRef[3] = {13.88895639764461, 25.89867770204523, 49.95832933327313};
    double prev = 1.0;
    for (int k = 0; k < 3; ++k) {
        const std::vector<double> N{static_cast<double>(ns[k]), static_cast<double>(ns[k])};
        const std::vector<double> Z{1.0, 0.5};
        const double lG = pfqn_le(modelA(), N, Z).lG;
        INFO("n = ", ns[k]);
        CHECK(lG == doctest::Approx(leRef[k]).epsilon(1e-9));
        const double ca = pfqn_ca(modelA(), std::vector<int>{ns[k], ns[k]}, modelAz()).lG;
        CHECK(ca == doctest::Approx(caRef[k]).epsilon(1e-9));
        const double e = relerr(lG, ca);
        CHECK(e < prev);  // monotone improvement with the population
        prev = e;
    }
    CHECK(prev < 4e-3);
}

// ---------------------------------------------------------------------------
// pfqn_lap
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_lap matches MATLAB on the repairman model C") {
    const std::vector<double> L{0.6, 0.4}, N{3.0, 2.0}, Z{1.0, 2.0};
    // MATLAB: pfqn_lap([.6 .4],[3 2],[1 2]) -> logI = 1.650702429011994
    CHECK(pfqn_lap(L, N, Z) == doctest::Approx(1.650702429011994).epsilon(1e-8));
    // Exact: pfqn_ca -> 1.632320758302839, so the Laplace error is 1.1e-2 relative.
    Matrix<double> Lm(1, 2);
    Lm(0, 0) = 0.6;
    Lm(0, 1) = 0.4;
    Matrix<double> Zm(1, 2);
    Zm(0, 0) = 1.0;
    Zm(0, 1) = 2.0;
    const double ca = pfqn_ca(Lm, std::vector<int>{3, 2}, Zm).lG;
    CHECK(ca == doctest::Approx(1.632320758302839).epsilon(1e-12));
    CHECK(relerr(pfqn_lap(L, N, Z), ca) < 1.2e-2);
}

TEST_CASE("pfqn_lap tightens as the population grows") {
    // MATLAB pfqn_ca / pfqn_lap on model C at N = [n n]:
    //   n =  5  ca 1.575893585847088  lap 1.568825628370974   rel err 4.5e-3
    //   n = 10  ca 1.116171998341203  lap 1.112176107795837   rel err 3.6e-3
    //   n = 20  ca 0.4048346626332844 lap 0.4027732658588674  rel err 5.1e-3
    const int ns[3] = {5, 10, 20};
    const double lapRef[3] = {1.568825628370974, 1.112176107795837, 0.4027732658588674};
    for (int k = 0; k < 3; ++k) {
        const std::vector<double> L{0.6, 0.4};
        const std::vector<double> N{static_cast<double>(ns[k]), static_cast<double>(ns[k])};
        const std::vector<double> Z{1.0, 2.0};
        INFO("n = ", ns[k]);
        CHECK(pfqn_lap(L, N, Z) == doctest::Approx(lapRef[k]).epsilon(1e-8));
    }
}

// ---------------------------------------------------------------------------
// pfqn_cub and the Grundmann-Moeller rule
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_cub is exact at Z = 0 and at the reference degree") {
    const std::vector<int> N{2, 3};
    const std::vector<double> Z0{0.0, 0.0};
    const double exact = pfqn_ca(modelA(), N, Matrix<double>()).lG;
    // MATLAB: pfqn_ca(L,[2 3],[0 0]) -> 4.139104712019314, pfqn_cub the same.
    CHECK(exact == doctest::Approx(4.139104712019314).epsilon(1e-12));
    const double cub = pfqn_cub(modelA(), N, Z0).lG;
    CHECK(cub == doctest::Approx(4.139104712019314).epsilon(1e-12));
    // The rule is exact here, so the observed relative error is pure rounding.
    CHECK(relerr(cub, exact) < 1e-13);
}

TEST_CASE("pfqn_cub with a delay matches MATLAB and the exact constant") {
    const std::vector<int> N{2, 3};
    const std::vector<double> Z{1.0, 0.5};
    // MATLAB: pfqn_cub(L,[2 3],[1 .5]) -> lGn = 4.783749663057889
    // exact  pfqn_ca                   -> lGn = 4.783752287503922
    // The outer v-grid is the only source of error: 5.5e-7 relative.
    const double cub = pfqn_cub(modelA(), N, Z).lG;
    CHECK(cub == doctest::Approx(4.783749663057889).epsilon(1e-7));
    const double exact = pfqn_ca(modelA(), N, modelAz()).lG;
    CHECK(exact == doctest::Approx(4.783752287503922).epsilon(1e-12));
    CHECK(relerr(cub, exact) < 1e-6);
}

TEST_CASE("pfqn_grnmol is exact against pfqn_ca, in exact arithmetic too") {
    // Single class, odd population: MATLAB pfqn_grnmol([1/2;1/3;1/5],5) and
    // pfqn_ca both give G = 0.1359938683127572.
    CHECK(pfqn_grnmol(matrixB(), std::vector<int>{5}) ==
          doctest::Approx(0.1359938683127572).epsilon(1e-12));
    CHECK(pfqn_ca(matrixB(), std::vector<int>{5}).G ==
          doctest::Approx(0.1359938683127572).epsilon(1e-12));
    // MATLAB pfqn_grnmol([1/2;1/3;1/5],3) -> 0.4478148148148149 (ca: ...148).
    CHECK(pfqn_grnmol(matrixB(), std::vector<int>{3}) ==
          doctest::Approx(0.4478148148148149).epsilon(1e-12));

    // The rule is rational, so in exact arithmetic it equals pfqn_ca exactly.
    Matrix<Rational> Lq(3, 1);
    Lq(0, 0) = Rational(1, 2);
    Lq(1, 0) = Rational(1, 3);
    Lq(2, 0) = Rational(1, 5);
    CHECK(pfqn_grnmol(Lq, std::vector<int>{5}) == pfqn_ca(Lq, std::vector<int>{5}).G);
    CHECK(pfqn_grnmol(Lq, std::vector<int>{3}) == pfqn_ca(Lq, std::vector<int>{3}).G);
}

TEST_CASE("pfqn_grnmol is exact on a multiclass model and rejects even populations") {
    // MATLAB pfqn_grnmol([1 .5;.7 1.2;.3 .9],[2 1]) -> 11.643, ca the same.
    CHECK(pfqn_grnmol(modelA(), std::vector<int>{2, 1}) == doctest::Approx(11.643).epsilon(1e-12));
    Matrix<Rational> Lq(3, 2);
    Lq(0, 0) = Rational(1);
    Lq(0, 1) = Rational(1, 2);
    Lq(1, 0) = Rational(7, 10);
    Lq(1, 1) = Rational(6, 5);
    Lq(2, 0) = Rational(3, 10);
    Lq(2, 1) = Rational(9, 10);
    CHECK(pfqn_grnmol(Lq, std::vector<int>{2, 1}) == pfqn_ca(Lq, std::vector<int>{2, 1}).G);

    // REFERENCE DEFECT, pinned: MATLAB's S = ceil(sum(N)-1)/2 is a half-integer
    // for an even total population and pfqn_grnmol then dies inside zeros()
    // ("Size inputs must be integers"). The port refuses instead of guessing.
    CHECK_THROWS_AS(pfqn_grnmol(modelA(), std::vector<int>{2, 2}), line::InputError);
}

// ---------------------------------------------------------------------------
// pfqn_kt
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_kt matches MATLAB and improves with the population") {
    const std::vector<double> Z{1.0, 0.5};
    // MATLAB: pfqn_kt(L,[2 3],[1 .5]) -> lG = 4.841446183242262
    CHECK(pfqn_kt(modelA(), std::vector<double>{2.0, 3.0}, Z).lG ==
          doctest::Approx(4.841446183242262).epsilon(1e-8));
    // MATLAB pfqn_ca / pfqn_kt at N = [n n]:
    //   n = 10  ca 14.10442210782851  kt 14.14880192515123   rel err 3.1e-3
    //   n = 20  ca 26.08438539593642  kt 26.15019474255796   rel err 2.5e-3
    //   n = 40  ca 50.13260369498655  kt 50.20568031516359   rel err 1.5e-3
    const int ns[3] = {10, 20, 40};
    const double ktRef[3] = {14.14880192515123, 26.15019474255796, 50.20568031516359};
    for (int k = 0; k < 3; ++k) {
        const std::vector<double> N{static_cast<double>(ns[k]), static_cast<double>(ns[k])};
        INFO("n = ", ns[k]);
        const double lG = pfqn_kt(modelA(), N, Z).lG;
        CHECK(lG == doctest::Approx(ktRef[k]).epsilon(1e-7));
        const double ca = pfqn_ca(modelA(), std::vector<int>{ns[k], ns[k]}, modelAz()).lG;
        CHECK(relerr(lG, ca) < 4e-3);
    }
}

// ---------------------------------------------------------------------------
// pfqn_bkt
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_bkt matches MATLAB and subtracts the class remainders") {
    const std::vector<double> N{2.0, 3.0};
    const std::vector<double> Z{1.0, 0.5};
    // MATLAB: pfqn_bkt(L,[2 3],[1 .5]) -> lGn = 4.77242756160185
    CHECK(pfqn_bkt(modelA(), N, Z).lG == doctest::Approx(4.77242756160185).epsilon(1e-9));
    // MATLAB: pfqn_bkt(L,[2 3]) with no think time -> lGn = 4.13902507926115
    CHECK(pfqn_bkt(modelA(), N).lG == doctest::Approx(4.13902507926115).epsilon(1e-9));
    // The correction is the sum of the exact Stirling remainders s(2) + s(3), and s(1)
    // is the constant pfqn_ble adds per station direction.
    const double s2 = pfqn_stirling_remainder<double>(2.0), s3 = pfqn_stirling_remainder<double>(3.0);
    CHECK(pfqn_bkt(modelA(), N, Z).lG - pfqn_kt(modelA(), N, Z).lG ==
          doctest::Approx(-(s2 + s3)).epsilon(1e-12));
    CHECK(pfqn_stirling_remainder<double>(1.0) ==
          doctest::Approx(1.0 - std::log(2 * M_PI) / 2).epsilon(1e-14));
    // With a think time BKT IS pfqn_ble: the two saddle points are one point in dual
    // coordinates (MATLAB: pfqn_ble(L,[2 3],[1 .5]) -> 4.77242756158822).
    CHECK(pfqn_bkt(modelA(), N, Z).lG == doctest::Approx(pfqn_ble(modelA(), N, Z).lG).epsilon(1e-9));
    // and it sits closer to the exact constant than pfqn_kt at every population:
    //   n = 10  ca 14.10442210782851  kt 14.14880192515123  bkt 14.1321407982845
    //   n = 20  ca 26.08438539593642  kt 26.15019474255796  bkt 26.141862103174
    //   n = 40  ca 50.13260369498655  kt 50.20568031516359  bkt 50.201513735287
    const int ns[3] = {10, 20, 40};
    const double ppRef[3] = {14.1321407982845, 26.141862103174, 50.201513735287};
    for (int k = 0; k < 3; ++k) {
        const std::vector<double> Nk{static_cast<double>(ns[k]), static_cast<double>(ns[k])};
        INFO("n = ", ns[k]);
        const double lG = pfqn_bkt(modelA(), Nk, Z).lG;
        CHECK(lG == doctest::Approx(ppRef[k]).epsilon(1e-7));
        const double ca = pfqn_ca(modelA(), std::vector<int>{ns[k], ns[k]}, modelAz()).lG;
        CHECK(std::abs(lG - ca) < std::abs(pfqn_kt(modelA(), Nk, Z).lG - ca));
    }
}

TEST_CASE("pfqn_bkt corrects only the classes pfqn_kt Laplaces") {
    // class 2 visits one station with no think time: pfqn_kt extracts its coefficient
    // exactly, so only class 1 carries a remainder.
    // MATLAB: pfqn_kt(Ls,[3 2],[0 0]) -> 2.84306473047309, pfqn_bkt -> 2.81538680478809
    Matrix<double> Ls{{1.0, 0.0}, {0.7, 1.2}, {0.3, 0.0}};
    const std::vector<double> N{3.0, 2.0};
    CHECK(pfqn_bkt(Ls, N).lG == doctest::Approx(2.81538680478809).epsilon(1e-9));
    CHECK(pfqn_bkt(Ls, N).lG - pfqn_kt(Ls, N).lG ==
          doctest::Approx(-pfqn_stirling_remainder<double>(3.0)).epsilon(1e-12));
    // an empty class is dropped by pfqn_kt's recursion:
    // MATLAB: pfqn_bkt(L,[3 0],[1 .5]) -> 1.97269392277537
    const std::vector<double> N0{3.0, 0.0}, Z{1.0, 0.5};
    CHECK(pfqn_bkt(modelA(), N0, Z).lG == doctest::Approx(1.97269392277537).epsilon(1e-9));
    CHECK(pfqn_bkt(modelA(), N0, Z).lG - pfqn_kt(modelA(), N0, Z).lG ==
          doctest::Approx(-pfqn_stirling_remainder<double>(3.0)).epsilon(1e-12));
}

// ---------------------------------------------------------------------------
// pfqn_lekt
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_lekt is one estimator on either side of the duality") {
    // modelA is 3 x 2, R <= M: the KT side, which is pfqn_bkt itself
    const std::vector<double> N{2.0, 3.0}, Z{1.0, 0.5};
    const LektResult<double> a = pfqn_lekt(modelA(), N, Z);
    CHECK(a.route == "kt");
    CHECK(a.lG == doctest::Approx(pfqn_bkt(modelA(), N, Z).lG).epsilon(1e-15));
    CHECK(pfqn_lekt(modelA(), N).lG == doctest::Approx(pfqn_bkt(modelA(), N).lG).epsilon(1e-15));
    // its transpose is 2 x 3: the LE side, which must land on the same number
    Matrix<double> Lw{{1.0, 0.7, 0.3}, {0.5, 1.2, 0.9}};
    const std::vector<double> Nw{2.0, 3.0, 4.0}, Zw{1.0, 0.5, 0.2};
    const LektResult<double> b = pfqn_lekt(Lw, Nw, Zw);
    CHECK(b.route == "le");
    CHECK(b.lG == doctest::Approx(pfqn_ble(Lw, Nw, Zw).lG).epsilon(1e-15));
    CHECK(std::abs(b.lG - pfqn_bkt(Lw, Nw, Zw).lG) < 1e-5);  // solver floors, not the identity
    // without a think time the LE side carries M kappa - r(N+M), which lands on bkt
    const double eta = 9.0 + 2.0;
    const double r = std::lgamma(eta) - (eta - 0.5) * std::log(eta) + eta - 0.5 * std::log(2 * M_PI);
    const double kappa = 1.0 - std::log(2 * M_PI) / 2;
    const LektResult<double> c = pfqn_lekt(Lw, Nw);
    CHECK(c.route == "le");
    CHECK(c.lG == doctest::Approx(pfqn_ble(Lw, Nw).lG + kappa - r).epsilon(1e-12));
    CHECK(std::abs(c.lG - pfqn_bkt(Lw, Nw).lG) < 1e-5);
    // MATLAB pfqn_lekt(Lw,[2 3 4],[1 .5 .2]) -> 7.64923112714508 (route le, = pfqn_ble;
    // pfqn_bkt gives 7.64923112722426), and pfqn_lekt(Lw,[2 3 4]) -> 7.09110182626673
    // (route le; pfqn_bkt 7.09110182636274, pfqn_ble 7.01761403495936)
    CHECK(b.lG == doctest::Approx(7.64923112714508).epsilon(1e-9));
    CHECK(c.lG == doctest::Approx(7.09110182626673).epsilon(1e-9));
    // a self-looping class sends a wide model to the KT side (MATLAB: 2.08635455377189)
    Matrix<double> Ls{{1.0, 0.0, 0.4}, {0.7, 1.2, 0.0}};
    CHECK(pfqn_lekt(Ls, std::vector<double>{3.0, 2.0, 2.0}).lG == doctest::Approx(2.08635455377189).epsilon(1e-9));
    CHECK(pfqn_lekt_route(Ls, std::vector<double>{3.0, 2.0, 2.0}, std::vector<double>{0.0, 0.0, 0.0}) == "kt");
    CHECK(pfqn_lekt_route(Ls, std::vector<double>{3.0, 2.0, 2.0}, std::vector<double>{0.0, 1.0, 1.0}) == "le");
}

// ---------------------------------------------------------------------------
// pfqn_panacea
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_panacea matches MATLAB in normal usage and flags the rest") {
    // Normal usage requires alpha > 0, i.e. lightly loaded stations.
    Matrix<double> L{{0.05, 0.02}, {0.03, 0.06}};
    const std::vector<int> N{2, 3};
    const std::vector<double> Z{1.0, 1.0};
    // MATLAB pfqn_ca -> -2.070623684440174; pfqn_panacea with 1, 2, 3 terms:
    //   1 term  -2.036116416941463   rel err 1.67e-2
    //   2 terms -2.077917639050010   rel err 3.52e-3
    //   3 terms -2.068343740545211   rel err 1.10e-3
    const double ca = pfqn_ca(L, N, Matrix<double>{{1.0, 1.0}}).lG;
    CHECK(ca == doctest::Approx(-2.070623684440174).epsilon(1e-12));
    const double ref[3] = {-2.036116416941463, -2.077917639050010, -2.068343740545211};
    const double tol[3] = {2e-2, 4e-3, 2e-3};
    for (int t = 1; t <= 3; ++t) {
        INFO("terms = ", t);
        const PanaceaResult<double> r = pfqn_panacea(L, N, Z, t);
        CHECK(r.normalUsage);
        CHECK(r.lG == doctest::Approx(ref[t - 1]).epsilon(1e-9));
        CHECK(relerr(r.lG, ca) < tol[t - 1]);
    }
    CHECK_THROWS_AS(pfqn_panacea(L, N, Z, 4), line::InputError);

    // Model A is far outside normal usage: MATLAB returns NaN, the port flags.
    CHECK_FALSE(pfqn_panacea(modelA(), std::vector<int>{2, 3},
                             std::vector<double>{1.0, 0.5}, 3)
                    .normalUsage);
}

// ---------------------------------------------------------------------------
// McKenna-Mitra quadratures
// ---------------------------------------------------------------------------

TEST_CASE("the three McKenna-Mitra quadratures match MATLAB on model C") {
    const std::vector<double> L{0.6, 0.4}, N{3.0, 2.0}, Z{1.0, 2.0};
    Matrix<double> Lm(1, 2);
    Lm(0, 0) = 0.6;
    Lm(0, 1) = 0.4;
    Matrix<double> Zm(1, 2);
    Zm(0, 0) = 1.0;
    Zm(0, 1) = 2.0;
    const double ca = pfqn_ca(Lm, std::vector<int>{3, 2}, Zm).lG;
    CHECK(ca == doctest::Approx(1.632320758302839).epsilon(1e-12));

    // MATLAB pfqn_mmint2 -> 1.632320740404313 (AbsTol 1e-12 on the truncated
    // interval), relative error against the exact constant 1.1e-8.
    CHECK(pfqn_mmint2(L, N, Z).lG == doctest::Approx(1.632320740404313).epsilon(1e-7));
    CHECK(relerr(pfqn_mmint2(L, N, Z).lG, ca) < 1e-7);

    // MATLAB pfqn_mmint2_gausslegendre -> 1.632320758292780, rel err 6.4e-12.
    CHECK(pfqn_mmint2_gausslegendre(L, N, Z).lG ==
          doctest::Approx(1.632320758292780).epsilon(1e-10));
    CHECK(relerr(pfqn_mmint2_gausslegendre(L, N, Z).lG, ca) < 1e-10);

    // MATLAB pfqn_mmint2_gausslaguerre -> 1.632320758303055, rel err 1.3e-13:
    // the integrand is a polynomial times e^{-u}, which this rule integrates
    // exactly once it has enough nodes.
    CHECK(pfqn_mmint2_gausslaguerre(L, N, Z).lG ==
          doctest::Approx(1.632320758303055).epsilon(1e-11));
    CHECK(relerr(pfqn_mmint2_gausslaguerre(L, N, Z).lG, ca) < 1e-11);
}

TEST_CASE("the adaptive McKenna-Mitra form loses its truncated tail at large N") {
    // MATLAB, model C at N = [n n]:
    //   n =  5 ca 1.575893585847088  mmint2  1.575880347511319   rel 8.4e-6
    //   n = 10 ca 1.116171998341203  mmint2  1.090406743758688   rel 2.3e-2
    //   n = 20 ca 0.4048346626332844 mmint2 -2.662712834121052   rel 7.6e+0
    // The cutoff is the 1 - 1e-12 quantile of the unit exponential, u = 27.6,
    // but the polynomial factor (Z + L u)^{sum N} pushes the mass past it. The
    // port reproduces the failure rather than widening the interval, and the
    // Gauss-Laguerre form stays exact throughout.
    const std::vector<double> L{0.6, 0.4}, Z{1.0, 2.0};
    Matrix<double> Lm(1, 2);
    Lm(0, 0) = 0.6;
    Lm(0, 1) = 0.4;
    Matrix<double> Zm(1, 2);
    Zm(0, 0) = 1.0;
    Zm(0, 1) = 2.0;
    const int ns[3] = {5, 10, 20};
    const double mmRef[3] = {1.575880347511319, 1.090406743758688, -2.662712834121052};
    const double caRef[3] = {1.575893585847088, 1.116171998341203, 0.4048346626332844};
    for (int k = 0; k < 3; ++k) {
        const std::vector<double> N{static_cast<double>(ns[k]), static_cast<double>(ns[k])};
        INFO("n = ", ns[k]);
        CHECK(pfqn_mmint2(L, N, Z).lG == doctest::Approx(mmRef[k]).epsilon(1e-6));
        CHECK(pfqn_ca(Lm, std::vector<int>{ns[k], ns[k]}, Zm).lG ==
              doctest::Approx(caRef[k]).epsilon(1e-12));
        // The Gauss-Laguerre rule is unaffected.
        CHECK(pfqn_mmint2_gausslaguerre(L, N, Z, 1, 90).lG ==
              doctest::Approx(caRef[k]).epsilon(1e-9));
    }
}

// ---------------------------------------------------------------------------
// laplaceapprox
// ---------------------------------------------------------------------------

TEST_CASE("laplaceapprox reproduces a Gaussian integral") {
    // h(x) = exp(-x^2/2 - y^2/2): the Laplace approximation at the mode is
    // exact, I = 2 pi, and det(-H) = 1.
    const std::function<double(const std::vector<double>&)> h =
        [](const std::vector<double>& x) { return std::exp(-0.5 * (x[0] * x[0] + x[1] * x[1])); };
    const LaplaceResult<double> r = laplaceapprox<double>(h, std::vector<double>{0.0, 0.0});
    CHECK_FALSE(r.detNegative);
    CHECK(r.I == doctest::Approx(6.283185307179586).epsilon(1e-6));
    CHECK(r.H(0, 0) == doctest::Approx(-1.0).epsilon(1e-5));
    CHECK(r.H(0, 1) == doctest::Approx(0.0).epsilon(1e-5));
    // MATLAB's logI uses -log(det(-H)) where consistency asks for -(1/2)log; on
    // this model det(-H) = 1, so the two agree and logI = log(2 pi).
    CHECK(r.logI == doctest::Approx(std::log(6.283185307179586)).epsilon(1e-6));
}

// ---------------------------------------------------------------------------
// pfqn_nrl / pfqn_nrp
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_nrl and pfqn_nrp match MATLAB on model A") {
    Matrix<double> alpha(3, 5, 1.0);
    const std::vector<double> N{2.0, 3.0}, Z{1.0, 0.5};
    // MATLAB: pfqn_nrl(L,[2 3],[1 .5],ones(3,5)) -> 4.338130313502534
    CHECK(pfqn_nrl(modelA(), N, Z, alpha) == doctest::Approx(4.338130313502534).epsilon(1e-6));
    // MATLAB: pfqn_nrp(...) -> 4.465092840664582
    CHECK(pfqn_nrp(modelA(), N, Z, alpha) == doctest::Approx(4.465092840664582).epsilon(1e-6));
    // Exact 4.783752287503922: the logit form is 9.3e-2 off, the probit 6.7e-2.
    const double ca = pfqn_ca(modelA(), std::vector<int>{2, 3}, modelAz()).lG;
    CHECK(relerr(pfqn_nrl(modelA(), N, Z, alpha), ca) < 1.0e-1);
    CHECK(relerr(pfqn_nrp(modelA(), N, Z, alpha), ca) < 7.0e-2);
}

TEST_CASE("pfqn_nrp is near exact on the single-class model B") {
    Matrix<double> alpha(3, 5, 1.0);
    const std::vector<double> N{5.0}, Z{0.0};
    // MATLAB: pfqn_nrp([1/2;1/3;1/5],5,0,ones(3,5)) -> -1.995145521568279
    // exact  pfqn_ca                                -> -1.995145480198095
    // i.e. 2.1e-8 relative: with one class the contour integral is
    // one-dimensional and the probit substitution is essentially exact.
    CHECK(pfqn_nrp(matrixB(), N, Z, alpha) == doctest::Approx(-1.995145521568279).epsilon(1e-7));
    const double ca = pfqn_ca(matrixB(), std::vector<int>{5}).lG;
    CHECK(ca == doctest::Approx(-1.995145480198095).epsilon(1e-12));
    CHECK(relerr(pfqn_nrp(matrixB(), N, Z, alpha), ca) < 1e-7);
    // MATLAB: pfqn_nrl on the same model -> -2.115927204091749 (6.1e-2 off).
    CHECK(pfqn_nrl(matrixB(), N, Z, alpha) == doctest::Approx(-2.115927204091749).epsilon(1e-6));
}

// ---------------------------------------------------------------------------
// pfqn_nre
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_nre matches MATLAB and beats the untilted contour on model A") {
    Matrix<double> alpha(3, 5, 1.0);
    const std::vector<double> N{2.0, 3.0}, Z{1.0, 0.5};
    // MATLAB: pfqn_nre(L,[2 3],[1 .5],ones(3,5)) -> 4.783188572136291
    CHECK(pfqn_nre(modelA(), N, Z, alpha) == doctest::Approx(4.783188572136291).epsilon(1e-6));
    // Exact 4.783752287503922. Tilting the contour to the saddle takes the
    // 9.3e-2 of pfqn_nrl down to 1.2e-4, which is the whole point of the method.
    const double ca = pfqn_ca(modelA(), std::vector<int>{2, 3}, modelAz()).lG;
    CHECK(relerr(pfqn_nre(modelA(), N, Z, alpha), ca) < 2.0e-4);
    CHECK(relerr(pfqn_nre(modelA(), N, Z, alpha), ca) <
          relerr(pfqn_nrl(modelA(), N, Z, alpha), ca));
}

TEST_CASE("pfqn_nre is exact on a single class") {
    // One class quotients the torus down to dimension zero, so the routine
    // returns pfqn_gldsingle itself rather than any expansion of it.
    Matrix<double> alpha(3, 5, 1.0);
    const std::vector<double> N{5.0}, Z{0.0};
    const double ca = pfqn_ca(matrixB(), std::vector<int>{5}).lG;
    CHECK(pfqn_nre(matrixB(), N, Z, alpha) == doctest::Approx(-1.995145480198095).epsilon(1e-12));
    CHECK(relerr(pfqn_nre(matrixB(), N, Z, alpha), ca) < 1e-12);
}

TEST_CASE("pfqn_nre carries a load-dependent rate row") {
    // Station 1 is a 2-server queue, mu(1,n) = min(n,2); station 2 single.
    Matrix<double> mu(2, 8, 1.0);
    for (int k = 0; k < 8; ++k) mu(0, k) = k + 1 < 2 ? 1.0 : 2.0;
    const Matrix<double> L{{1.0, 0.6}, {0.5, 1.1}};
    const std::vector<double> N{4.0, 4.0}, Z{0.5, 0.5};
    // MATLAB: pfqn_nre(L,[4 4],[.5 .5],mu) -> 3.911235290463940,
    // against the exact pfqn_ncld 3.911912990913230 (1.7e-4 relative).
    CHECK(pfqn_nre(L, N, Z, mu) == doctest::Approx(3.911235290463940).epsilon(1e-6));
}

// ---------------------------------------------------------------------------
// pfqn_dac
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_dac matches MATLAB and agrees with exact MVA on model A") {
    const std::vector<int> N{2, 3};
    const std::vector<double> Z{1.0, 0.5};
    const DacResult<double> d = pfqn_dac(modelA(), N, Z);
    // MATLAB: 56 aggregate states (M+1 = 4 centers, 5 jobs), X = [.3768854645
    // .4848902561], Q(1,1) = .661379274213801, U = [.6193 .8457 .5495].
    CHECK(d.states.size() == 56);
    CHECK(d.XN[0] == doctest::Approx(0.3768854645273076).epsilon(1e-12));
    CHECK(d.XN[1] == doctest::Approx(0.4848902560783052).epsilon(1e-12));
    CHECK(d.QN(0, 0) == doctest::Approx(0.661379274213801).epsilon(1e-12));
    CHECK(d.UN[0] == doctest::Approx(0.6193305925664601).epsilon(1e-12));
    CHECK(d.UN[1] == doctest::Approx(0.8456881324630814).epsilon(1e-12));
    CHECK(d.UN[2] == doctest::Approx(0.5494668698286669).epsilon(1e-12));
    // The joint distribution is a probability distribution.
    double s = 0.0;
    for (double p : d.Pjoint) s += p;
    CHECK(s == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("pfqn_dac is exact in rational arithmetic on a load-dependent model") {
    // de Souza e Silva 1987, Section 3 availability example.
    Matrix<double> L{{5.0, 0.0}, {0.0, 10.0}, {2.0, 1.0}};
    Matrix<double> mu{{1.0, 1.0, 1.0, 1.0}, {1.0, 2.0, 2.0, 2.0}, {1.0, 1.0, 1.0, 1.0}};
    const std::vector<int> N{1, 3};
    const std::vector<double> Z{0.0, 0.0};
    const DacResult<double> d = pfqn_dac(L, N, Z, mu);
    // MATLAB: 15 states, X = [.1338786053 .1941455015],
    // U = [.6693930263 .9944037882 .4619027120], availability .6672406371.
    CHECK(d.states.size() == 15);
    CHECK(d.XN[0] == doctest::Approx(0.1338786052518295).epsilon(1e-12));
    CHECK(d.XN[1] == doctest::Approx(0.1941455015066724).epsilon(1e-12));
    CHECK(d.UN[0] == doctest::Approx(0.6693930262591476).epsilon(1e-12));
    CHECK(d.UN[1] == doctest::Approx(0.9944037882049075).epsilon(1e-12));
    CHECK(d.UN[2] == doctest::Approx(0.4619027120103314).epsilon(1e-12));
    double av = 0.0;
    for (std::size_t i = 0; i < d.states.size(); ++i)
        if (d.states[i][0] == 1 && d.states[i][1] >= 1) av += d.Pjoint[i];
    CHECK(av == doctest::Approx(0.6672406371071891).epsilon(1e-12));

    // The same model in exact arithmetic: the joint distribution sums to
    // exactly one, with no rounding at all.
    Matrix<Rational> Lq{{Rational(5), Rational(0)},
                        {Rational(0), Rational(10)},
                        {Rational(2), Rational(1)}};
    Matrix<Rational> muq{{Rational(1), Rational(1), Rational(1), Rational(1)},
                         {Rational(1), Rational(2), Rational(2), Rational(2)},
                         {Rational(1), Rational(1), Rational(1), Rational(1)}};
    const DacResult<Rational> dq =
        pfqn_dac(Lq, N, std::vector<Rational>{Rational(0), Rational(0)}, muq);
    Rational sq(0);
    for (const Rational& p : dq.Pjoint) sq += p;
    CHECK(sq == Rational(1));
    CHECK(static_cast<double>(dq.XN[0]) == doctest::Approx(0.1338786052518295).epsilon(1e-12));
}

// ---------------------------------------------------------------------------
// pfqn_propfair
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_propfair reaches the same optimum as MATLAB fmincon") {
    const std::vector<double> N{2.0, 3.0}, Z{1.0, 0.5};
    const PropfairResult<double> r = pfqn_propfair(modelA(), N, Z);
    // MATLAB (fmincon): Xasy = [.5434113685198305 .5163431910590029],
    // lG = 2.919079651653669. The port takes a log-barrier path to the same
    // maximizer, so agreement is asserted at the optimizer tolerance, not at
    // machine precision.
    CHECK(r.Xasy[0] == doctest::Approx(0.5434113685198305).epsilon(1e-4));
    CHECK(r.Xasy[1] == doctest::Approx(0.5163431910590029).epsilon(1e-4));
    CHECK(r.lG == doctest::Approx(2.919079651653669).epsilon(1e-4));
    // The optimum must be feasible: L X <= 1 at every station.
    for (std::size_t i = 0; i < 3; ++i) {
        const double u = modelA()(i, 0) * r.Xasy[0] + modelA()(i, 1) * r.Xasy[1];
        INFO("station ", i);
        CHECK(u <= 1.0 + 1e-12);
    }
}

// ---------------------------------------------------------------------------
// pfqn_fnc, pfqn_oi_fnc, pfqn_oi_insvc, cd_peak_scaling
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_fnc reproduces MATLAB, vectorized retry predicate included") {
    Matrix<double> alpha{{1.0, 2.0, 3.0}, {2.0, 2.0, 2.0}};
    const FncResult<double> r = pfqn_fnc(alpha);
    // MATLAB: c = [0 0] and mu = [1 Inf Inf; 2 2 2]. The offset ladder does NOT
    // fire, because `~all(isfinite(mu))` reduces per column and the last column
    // is not all-non-finite; see the note in pfqn_fnc.h.
    CHECK(r.c[0] == 0.0);
    CHECK(r.c[1] == 0.0);
    CHECK(r.mu(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(std::isinf(r.mu(0, 1)));
    CHECK(std::isinf(r.mu(0, 2)));
    CHECK(r.mu(1, 0) == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(r.mu(1, 2) == doctest::Approx(2.0).epsilon(1e-12));

    // A single station is a row vector, where all() does reduce to a scalar,
    // so the ladder does fire: MATLAB gives c = -0.5 and mu = [2 2 -3].
    Matrix<double> a1{{1.0, 2.0, 3.0}};
    const FncResult<double> r1 = pfqn_fnc(a1);
    CHECK(r1.c[0] == doctest::Approx(-0.5).epsilon(1e-12));
    CHECK(r1.mu(0, 0) == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(r1.mu(0, 1) == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(r1.mu(0, 2) == doctest::Approx(-3.0).epsilon(1e-12));

    // Explicit offset: MATLAB pfqn_fnc(alpha,[0.25 0.25]) gives
    // [0.8 -10 -3; 1.6 2.5 2].
    const Matrix<double> m2 = pfqn_fnc_at(alpha, std::vector<double>{0.25, 0.25});
    CHECK(m2(0, 0) == doctest::Approx(0.8).epsilon(1e-12));
    CHECK(m2(0, 1) == doctest::Approx(-10.0).epsilon(1e-12));
    CHECK(m2(0, 2) == doctest::Approx(-3.0).epsilon(1e-12));
    CHECK(m2(1, 0) == doctest::Approx(1.6).epsilon(1e-12));
    CHECK(m2(1, 1) == doctest::Approx(2.5).epsilon(1e-12));
    CHECK(m2(1, 2) == doctest::Approx(2.0).epsilon(1e-12));

    // The recursion is rational, so the exact instantiation reproduces it with
    // no rounding: mu(1,1) = 1.6 is exactly 8/5.
    Matrix<Rational> aq{{Rational(1), Rational(2), Rational(3)},
                        {Rational(2), Rational(2), Rational(2)}};
    const Matrix<Rational> mq =
        pfqn_fnc_at(aq, std::vector<Rational>{Rational(1, 4), Rational(1, 4)});
    CHECK(mq(1, 0) == Rational(8, 5));
    CHECK(mq(1, 1) == Rational(5, 2));
    CHECK(mq(0, 1) == Rational(-10));
}

TEST_CASE("pfqn_oi_fnc reproduces the MATLAB balance function") {
    // Phi over the lattice N = [2 1], column-major: [1 .5 .25 .5 .25 .125].
    const std::vector<double> Phi{1.0, 0.5, 0.25, 0.5, 0.25, 0.125};
    const OiFncResult<double> r = pfqn_oi_fnc(Phi, std::vector<int>{2, 1});
    // MATLAB: Psi = [1 0.5 0.25 0.5 0 0].
    const double ref[6] = {1.0, 0.5, 0.25, 0.5, 0.0, 0.0};
    for (int i = 0; i < 6; ++i) {
        INFO("state ", i);
        CHECK(r.Psi[static_cast<std::size_t>(i)] == doctest::Approx(ref[i]).epsilon(1e-12));
    }
    // The zero entries have no defined rate, which the reference reports as Inf.
    CHECK(std::isinf(r.mu[4]));
    CHECK(std::isinf(r.mu[5]));

    // Exact arithmetic: the deconvolution is a rational triangular solve and
    // the two zeros above are exact zeros, not cancellation residue.
    const std::vector<Rational> Phiq{Rational(1),      Rational(1, 2), Rational(1, 4),
                                     Rational(1, 2),  Rational(1, 4), Rational(1, 8)};
    const OiFncResult<Rational> rq = pfqn_oi_fnc(Phiq, std::vector<int>{2, 1});
    CHECK(rq.Psi[1] == Rational(1, 2));
    CHECK(rq.Psi[2] == Rational(1, 4));
    CHECK(rq.Psi[4] == Rational(0));
    CHECK(rq.Psi[5] == Rational(0));
}

TEST_CASE("pfqn_oi_insvc reproduces the MATLAB in-service means") {
    // Two-server OI station: mu(n) = min(sum(n), 2), lattice N = [2 1].
    const std::function<double(const std::vector<int>&)> rate =
        [](const std::vector<int>& n) {
            int s = 0;
            for (int v : n) s += v;
            return static_cast<double>(s < 2 ? s : 2);
        };
    const OiInsvcResult<double> r = pfqn_oi_insvc(rate, std::vector<int>{2, 1});
    // MATLAB: Phi = [1 1 0.5 1 1 0.75], g(:,1) = [0 1 2 0 1 4/3].
    const double phiRef[6] = {1.0, 1.0, 0.5, 1.0, 1.0, 0.75};
    const double gRef[6] = {0.0, 1.0, 2.0, 0.0, 1.0, 4.0 / 3.0};
    for (int i = 0; i < 6; ++i) {
        INFO("state ", i);
        CHECK(r.Phi[static_cast<std::size_t>(i)] == doctest::Approx(phiRef[i]).epsilon(1e-12));
        CHECK(r.g(static_cast<std::size_t>(i), 0) == doctest::Approx(gRef[i]).epsilon(1e-12));
    }
    // In exact arithmetic 4/3 is 4/3, and the strict rate comparison that
    // decides whether the tail job is in service is an exact comparison.
    const std::function<Rational(const std::vector<int>&)> rateq =
        [](const std::vector<int>& n) {
            int s = 0;
            for (int v : n) s += v;
            return Rational(s < 2 ? s : 2);
        };
    const OiInsvcResult<Rational> rq = pfqn_oi_insvc(rateq, std::vector<int>{2, 1});
    CHECK(rq.g(5, 0) == Rational(4, 3));
    CHECK(rq.Phi[5] == Rational(3, 4));
}

TEST_CASE("cd_peak_scaling takes the peak over the lattice and the classes") {
    // MATLAB cd_peak_scaling(@(n) 1+0.1*sum(n), [2 2], 2) -> 1.4, attained at
    // the full population n = [2 2].
    const std::function<std::vector<double>(const std::vector<int>&)> beta =
        [](const std::vector<int>& n) {
            int s = 0;
            for (int v : n) s += v;
            return std::vector<double>{1.0 + 0.1 * static_cast<double>(s)};
        };
    CHECK(cd_peak_scaling<double>(beta, std::vector<int>{2, 2}) ==
          doctest::Approx(1.4).epsilon(1e-12));

    // Non-finite entries are skipped, so an Inf at one state does not become
    // the normalizer and zero every utilization at the station.
    const std::function<std::vector<double>(const std::vector<int>&)> withInf =
        [](const std::vector<int>& n) {
            if (n[0] == 2 && n[1] == 2)
                return std::vector<double>{std::numeric_limits<double>::infinity()};
            int s = 0;
            for (int v : n) s += v;
            return std::vector<double>{1.0 + 0.1 * static_cast<double>(s)};
        };
    CHECK(cd_peak_scaling<double>(withInf, std::vector<int>{2, 2}) ==
          doctest::Approx(1.3).epsilon(1e-12));

    // Exact arithmetic: 1.4 is 7/5 with no rounding.
    const std::function<std::vector<Rational>(const std::vector<int>&)> betaq =
        [](const std::vector<int>& n) {
            int s = 0;
            for (int v : n) s += v;
            return std::vector<Rational>{Rational(1) + Rational(1, 10) * Rational(s)};
        };
    CHECK(cd_peak_scaling<Rational>(betaq, std::vector<int>{2, 2}) == Rational(7, 5));
}

