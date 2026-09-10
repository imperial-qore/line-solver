/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
/**
 * Covers api/mam/dmap.h, map_algebra.h, map_max.h, map_rand.h, mmap_stats.h and
 * map_dist.h. The literals are MATLAB output on the same inputs; the identities
 * around them are what pins the parts a single reference vector cannot reach.
 */
#include <random>
#include <vector>

#include "doctest.h"
#include "line/api/mam/dmap.h"
#include "line/api/mam/map_algebra.h"
#include "line/api/mam/map_dist.h"
#include "line/api/mam/map_gamma.h"
#include "line/api/mam/map_max.h"
#include "line/api/mam/map_rand.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_stats.h"

using line::Matrix;
using namespace line::mam;

namespace {

Matrix<double> mk(std::size_t r, std::size_t c, const std::vector<double>& v) {
    Matrix<double> A(r, c);
    std::size_t k = 0;
    for (std::size_t i = 0; i < r; ++i)
        for (std::size_t j = 0; j < c; ++j) A(i, j) = v[k++];
    return A;
}

/** The three-phase two-class MMAP used throughout, and its MATLAB reference. */
Mmap<double> ref_mmap() {
    Mmap<double> mm;
    mm.D0 = mk(3, 3, {-6, 1, 0.5, 0.2, -3, 0.3, 0.1, 0.4, -2});
    const Matrix<double> Da = mk(3, 3, {1, 0.5, 0, 0.5, 0.5, 0.2, 0.3, 0.2, 0.1});
    const Matrix<double> Db = mk(3, 3, {1.5, 0.5, 1, 0.8, 0.4, 0.1, 0.6, 0.2, 0.1});
    mm.D1 = Matrix<double>(3, 3, 0.0);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) mm.D1(i, j) = Da(i, j) + Db(i, j);
    mm.Dc = {Da, Db};
    return mm;
}

Map<double> ref_map2() { return Map<double>{mk(2, 2, {-3, 1, 0.5, -4}), mk(2, 2, {1, 1, 1.5, 2})}; }

Map<double> ref_map3() {
    return Map<double>{mk(3, 3, {-6, 1, 0.5, 0.2, -3, 0.3, 0.1, 0.4, -2}),
                       mk(3, 3, {2.5, 1, 1, 1.3, 0.9, 0.3, 0.9, 0.4, 0.2})};
}

Dmap<double> ref_dmap_a() {
    return Dmap<double>{mk(2, 2, {0.3, 0.1, 0.2, 0.4}), mk(2, 2, {0.4, 0.2, 0.1, 0.3})};
}

Dmap<double> ref_dmap_b() {
    return Dmap<double>{mk(2, 2, {0.5, 0.1, 0.05, 0.6}), mk(2, 2, {0.3, 0.1, 0.15, 0.2})};
}

}  // namespace

TEST_CASE("mmap_pie and mmap_sigma are the class-transition laws") {
    const Mmap<double> mm = ref_mmap();
    const Matrix<double> pie = mmap_pie(mm);
    // One row per class, each a distribution over the phases seen at an arrival
    REQUIRE(pie.rows() == 2);
    REQUIRE(pie.cols() == 3);
    for (std::size_t c = 0; c < 2; ++c) {
        double s = 0;
        for (std::size_t j = 0; j < 3; ++j) s += pie(c, j);
        CHECK(s == doctest::Approx(1.0).epsilon(1e-12));
    }
    CHECK(pie(0, 0) == doctest::Approx(0.5259009009009).epsilon(1e-11));
    CHECK(pie(0, 1) == doctest::Approx(0.36936936936937).epsilon(1e-11));
    CHECK(pie(0, 2) == doctest::Approx(0.10472972972973).epsilon(1e-11));
    CHECK(pie(1, 0) == doctest::Approx(0.56843679880329).epsilon(1e-11));
    CHECK(pie(1, 1) == doctest::Approx(0.22139117427076).epsilon(1e-11));
    CHECK(pie(1, 2) == doctest::Approx(0.21017202692595).epsilon(1e-11));

    const Matrix<double> sig = mmap_sigma(mm);
    double tot = 0;
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) tot += sig(i, j);
    CHECK(tot == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(sig(0, 0) == doctest::Approx(0.16170477900893).epsilon(1e-11));
    CHECK(sig(0, 1) == doctest::Approx(0.23739634458658).epsilon(1e-11));
    CHECK(sig(1, 1) == doctest::Approx(0.36350253181792).epsilon(1e-11));

    const std::vector<std::vector<std::vector<double>>> s2 = mmap_sigma2(mm);
    double t2 = 0;
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            for (std::size_t h = 0; h < 2; ++h) t2 += s2[i][j][h];
    CHECK(t2 == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(s2[0][0][0] == doctest::Approx(0.065637020845023).epsilon(1e-11));
    CHECK(s2[1][1][1] == doctest::Approx(0.21997195843953).epsilon(1e-11));
    // Marginalising the third step must give back the two-step law
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) {
            double m = 0;
            for (std::size_t h = 0; h < 2; ++h) m += s2[i][j][h];
            CHECK(m == doctest::Approx(sig(i, j)).epsilon(1e-11));
        }
}

TEST_CASE("mmap counting statistics agree with the MATLAB reference") {
    const Mmap<double> mm = ref_mmap();
    const std::vector<double> cm = mmap_count_mean(mm, 2.0);
    CHECK(cm[0] == doctest::Approx(2.1397590361446).epsilon(1e-11));
    CHECK(cm[1] == doctest::Approx(3.221686746988).epsilon(1e-11));
    // The counting mean is linear in t: lambda_c * t
    const std::vector<double> cm1 = mmap_count_mean(mm, 1.0);
    CHECK(cm[0] == doctest::Approx(2.0 * cm1[0]).epsilon(1e-12));
    CHECK(cm[1] == doctest::Approx(2.0 * cm1[1]).epsilon(1e-12));

    const std::vector<double> cidc = mmap_count_idc(mm, 2.0);
    CHECK(cidc[0] == doctest::Approx(1.1104809064634).epsilon(1e-11));
    CHECK(cidc[1] == doctest::Approx(1.2181321823748).epsilon(1e-11));

    const Matrix<double> mcov = mmap_count_mcov(mm, 2.0);
    CHECK(mcov(0, 0) == doctest::Approx(2.3761615540711).epsilon(1e-11));
    CHECK(mcov(0, 1) == doctest::Approx(0.36036160954962).epsilon(1e-11));
    CHECK(mcov(1, 0) == doctest::Approx(mcov(0, 1)).epsilon(1e-12));
    CHECK(mcov(1, 1) == doctest::Approx(3.9244403080364).epsilon(1e-11));
    // The per-class index of dispersion is the diagonal of the covariance over the mean
    CHECK(mcov(0, 0) / cm[0] == doctest::Approx(cidc[0]).epsilon(1e-11));
    CHECK(mcov(1, 1) / cm[1] == doctest::Approx(cidc[1]).epsilon(1e-11));

    const std::vector<double> idc = mmap_idc(mm);
    CHECK(idc[0] == doctest::Approx(1.1313501012471).epsilon(1e-11));
    CHECK(idc[1] == doctest::Approx(1.2450365867246).epsilon(1e-11));

    const Matrix<double> cross = mmap_cross_moment(mm, 2u);
    CHECK(cross(0, 0) == doctest::Approx(0.33270581383162).epsilon(1e-11));
    CHECK(cross(0, 1) == doctest::Approx(0.30081211306525).epsilon(1e-11));
    CHECK(cross(1, 0) == doctest::Approx(0.36703533471097).epsilon(1e-11));
    CHECK(cross(1, 1) == doctest::Approx(0.32480194111176).epsilon(1e-11));

    const Matrix<double> fwd = mmap_forward_moment(mm, std::vector<unsigned>{1u, 3u});
    CHECK(fwd(0, 0) == doctest::Approx(0.36457894709482).epsilon(1e-11));
    CHECK(fwd(0, 1) == doctest::Approx(0.44641883465439).epsilon(1e-11));
    CHECK(fwd(1, 0) == doctest::Approx(0.37864913611055).epsilon(1e-11));
    CHECK(fwd(1, 1) == doctest::Approx(0.50521367645536).epsilon(1e-11));
}

TEST_CASE("mmap_timereverse is an involution and preserves the arrival rates") {
    const Mmap<double> mm = ref_mmap();
    const Mmap<double> r = mmap_timereverse(mm);
    CHECK(r.D0(0, 1) == doctest::Approx(0.29090909090909).epsilon(1e-11));
    CHECK(r.Dc[0](0, 1) == doctest::Approx(0.72727272727273).epsilon(1e-11));
    const std::vector<double> l0 = mmap_count_mean(mm, 1.0), l1 = mmap_count_mean(r, 1.0);
    CHECK(l1[0] == doctest::Approx(l0[0]).epsilon(1e-11));
    CHECK(l1[1] == doctest::Approx(l0[1]).epsilon(1e-11));
    const Mmap<double> rr = mmap_timereverse(r);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            CHECK(rr.D0(i, j) == doctest::Approx(mm.D0(i, j)).epsilon(1e-10));
            CHECK(rr.Dc[0](i, j) == doctest::Approx(mm.Dc[0](i, j)).epsilon(1e-10));
        }
}

TEST_CASE("mmap_maps, mmap_embedded and mmap_sum keep the reference layout") {
    const Mmap<double> mm = ref_mmap();
    const std::vector<Map<double>> ms = mmap_maps(mm);
    REQUIRE(ms.size() == 2);
    // The per-class MAP carries only that class's arrivals; the others fold into D0
    for (std::size_t c = 0; c < 2; ++c)
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) {
                CHECK(ms[c].D1(i, j) == doctest::Approx(mm.Dc[c](i, j)).epsilon(1e-12));
                CHECK(ms[c].D0(i, j) + ms[c].D1(i, j) ==
                      doctest::Approx(mm.D0(i, j) + mm.D1(i, j)).epsilon(1e-12));
            }

    const std::vector<Matrix<double>> emb = mmap_embedded(mm);
    REQUIRE(emb.size() == 2);
    CHECK(emb[0](0, 0) == doctest::Approx(0.21696480092325).epsilon(1e-11));
    CHECK(emb[1](0, 0) == doctest::Approx(0.33641084824005).epsilon(1e-11));
    // The two per-class embedded matrices together form a stochastic matrix
    for (std::size_t i = 0; i < 3; ++i) {
        double s = 0;
        for (std::size_t c = 0; c < 2; ++c)
            for (std::size_t j = 0; j < 3; ++j) s += emb[c](i, j);
        CHECK(s == doctest::Approx(1.0).epsilon(1e-11));
    }

    const Mmap<double> sm = mmap_sum(mm, 2u);
    REQUIRE(sm.order() == 6);
    CHECK(sm.D0(0, 3) == doctest::Approx(2.5).epsilon(1e-12));
    CHECK(sm.Dc[0](3, 0) == doctest::Approx(1.0).epsilon(1e-12));
    // Summing pairs of interarrival times halves the arrival rate
    const std::vector<double> lam = mmap_count_mean(mm, 1.0), lam2 = mmap_count_mean(sm, 1.0);
    CHECK(lam2[0] + lam2[1] == doctest::Approx(0.5 * (lam[0] + lam[1])).epsilon(1e-10));
}

TEST_CASE("map_timereverse and map_kpc match the reference construction") {
    const Map<double> a = ref_map2();
    const Map<double> r = map_timereverse(a);
    CHECK(r.D0(0, 0) == doctest::Approx(-3.0).epsilon(1e-12));
    CHECK(r.D0(0, 1) == doctest::Approx(0.5).epsilon(1e-12));
    CHECK(r.D1(0, 1) == doctest::Approx(1.5).epsilon(1e-12));
    // Time reversal preserves the interarrival moments and is an involution
    CHECK(map_mean(r) == doctest::Approx(map_mean(a)).epsilon(1e-12));
    CHECK(map_moment(r, 2u) == doctest::Approx(map_moment(a, 2u)).epsilon(1e-11));
    const Map<double> rr = map_timereverse(r);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            CHECK(rr.D1(i, j) == doctest::Approx(a.D1(i, j)).epsilon(1e-10));

    const Map<double> k = map_kpc(a, map_exponential(0.5));
    // Composing with a one-phase MAP scales both blocks, D0 with the sign flip
    CHECK(k.D0(0, 0) == doctest::Approx(-1.5).epsilon(1e-12));
    CHECK(k.D0(1, 0) == doctest::Approx(0.25).epsilon(1e-12));
    CHECK(k.D1(0, 0) == doctest::Approx(0.5).epsilon(1e-12));
    CHECK(k.D1(1, 1) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(map_largemap() == 100u);
}

TEST_CASE("map_gamma2 is the subdominant eigenvalue of the embedded chain") {
    const std::complex<double> g = map_gamma2(ref_map2());
    CHECK(g.real() == doctest::Approx(0.043478260869565).epsilon(1e-11));
    CHECK(g.imag() == doctest::Approx(0.0).epsilon(1e-14));
    // On an order-two MAP the acf really is geometric with that ratio
    const std::vector<double> acf = map_acf(ref_map2(), std::vector<unsigned>{1u, 2u});
    CHECK(acf[1] / acf[0] == doctest::Approx(g.real()).epsilon(1e-9));
}

TEST_CASE("map_max of two exponentials is the maximum of two clocks") {
    const double la = 0.7, lb = 1.3;
    const Map<double> m = map_max(map_exponential(la), map_exponential(lb));
    // E[max] = 1/la + 1/lb - 1/(la+lb)
    CHECK(map_mean(m) == doctest::Approx(1.0 / la + 1.0 / lb - 1.0 / (la + lb)).epsilon(1e-10));
    // The maximum is never faster than either component
    const Map<double> a = ref_map2();
    CHECK(map_mean(map_max(a, a)) >= map_mean(a));
}

TEST_CASE("mmap_max preserves the class count and slows the stream") {
    const Mmap<double> mm = ref_mmap();
    const Mmap<double> x = mmap_max(mm, mm, 2u);
    REQUIRE(x.classes() == 2);
    const std::vector<double> l0 = mmap_count_mean(mm, 1.0), lx = mmap_count_mean(x, 1.0);
    CHECK(lx[0] + lx[1] < l0[0] + l0[1]);
}

TEST_CASE("the random MAP generators produce feasible processes") {
    std::mt19937_64 gen(20260801u);
    for (int trial = 0; trial < 8; ++trial) {
        const Map<double> m = map_rand<double>(3, gen);
        for (std::size_t i = 0; i < 3; ++i) {
            double row = 0;
            for (std::size_t j = 0; j < 3; ++j) {
                row += m.D0(i, j) + m.D1(i, j);
                if (i != j) CHECK(m.D0(i, j) >= 0.0);
                CHECK(m.D1(i, j) >= 0.0);
            }
            CHECK(m.D0(i, i) < 0.0);
            CHECK(row == doctest::Approx(0.0).epsilon(1e-12).scale(1.0));
        }
        CHECK(map_mean(m) > 0.0);
        // An MMPP has arrivals only on the diagonal
        const Map<double> p = mmpp_rand<double>(3, gen);
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j)
                if (i != j) CHECK(p.D1(i, j) == doctest::Approx(0.0).epsilon(1e-14));
        // A hyperexponential read back through ph2hyper reproduces its own mean
        const Map<double> h = hyper_rand<double>(3, gen);
        const HyperParams<double> hp = ph2hyper(h);
        double mean = 0;
        for (std::size_t i = 0; i < 3; ++i) mean += hp.prob[i] / hp.lambda[i];
        CHECK(mean == doctest::Approx(map_mean(h)).epsilon(1e-9));
        // An acyclic PH is upper triangular in D0
        const Map<double> aph = aph_rand<double>(4, gen);
        for (std::size_t i = 0; i < 4; ++i)
            for (std::size_t j = 0; j < i; ++j)
                CHECK(aph.D0(i, j) == doctest::Approx(0.0).epsilon(1e-14));
        const Mmap<double> mr = mmap_rand<double>(3, 2, gen);
        REQUIRE(mr.classes() == 2);
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j)
                CHECK(mr.Dc[0](i, j) + mr.Dc[1](i, j) == doctest::Approx(mr.D1(i, j)).epsilon(1e-12));
    }
}

TEST_CASE("the discrete-time MAP statistics agree with the MATLAB reference") {
    const Dmap<double> a = ref_dmap_a(), b = ref_dmap_b();
    const std::vector<double> pie = dmap_pie(a);
    CHECK(pie[0] == doctest::Approx(0.5).epsilon(1e-12));
    CHECK(pie[1] == doctest::Approx(0.5).epsilon(1e-12));
    const std::vector<double> mom = dmap_moment(a, std::vector<unsigned>{1u, 2u, 3u});
    CHECK(mom[0] == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(mom[1] == doctest::Approx(6.0).epsilon(1e-12));
    CHECK(mom[2] == doctest::Approx(26.0).epsilon(1e-12));
    CHECK(dmap_isfeasible(a));
    CHECK(dmap_isfeasible(b));
    // A D-MAP with a row sum away from one is not a D-MAP
    Dmap<double> bad = a;
    bad.D1(0, 0) += 0.5;
    CHECK_FALSE(dmap_isfeasible(bad));

    CHECK(dmap_exp_mul_int(a, b, 3u) == doctest::Approx(0.020761089927891).epsilon(1e-11));
    CHECK(dmap_dist(a, b, 3u) == doctest::Approx(0.0038864114933015).epsilon(1e-10));
    CHECK(dmap_geo_mul_sum(a, b) == doctest::Approx(28.717948717949).epsilon(1e-11));
    CHECK(dmap_dist_acf(a, b) == doctest::Approx(0.13723685676097).epsilon(1e-10));
    CHECK(dmap_dist_lag1(a, b) == doctest::Approx(0.014521770615296).epsilon(1e-10));
    // Every distance vanishes between a process and itself
    CHECK(dmap_dist(a, a, 3u) == doctest::Approx(0.0).epsilon(1e-12).scale(1.0));
    CHECK(dmap_dist_acf(a, a) == doctest::Approx(0.0).epsilon(1e-12).scale(1.0));
    CHECK(dmap_dist_lag1(a, a) == doctest::Approx(0.0).epsilon(1e-12).scale(1.0));
    CHECK(dmap_dist(a, b, 3u) == doctest::Approx(dmap_dist(b, a, 3u)).epsilon(1e-11));
}

TEST_CASE("dmap_sample follows the slot lattice and reproduces the mean") {
    std::mt19937_64 gen(4242u);
    const Dmap<double> a = ref_dmap_a();
    const std::vector<unsigned> s = dmap_sample(a, 200000, gen);
    REQUIRE(s.size() == 200000);
    double sum = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        CHECK(s[i] >= 1u);  // an interarrival is at least one slot
        sum += s[i];
    }
    CHECK(sum / static_cast<double>(s.size()) == doctest::Approx(2.0).epsilon(0.02));
}

TEST_CASE("the continuous MAP distances vanish on a repeated argument") {
    const Map<double> a = ref_map2(), b = ref_map3();
    CHECK(map_dist(a, a, 3u) == doctest::Approx(0.0).epsilon(1e-12).scale(1.0));
    CHECK(map_dist_acf(a, a) == doctest::Approx(0.0).epsilon(1e-12).scale(1.0));
    CHECK(map_dist_lag1(a, a) == doctest::Approx(0.0).epsilon(1e-12).scale(1.0));
    CHECK(map_dist(a, b, 3u) == doctest::Approx(map_dist(b, a, 3u)).epsilon(1e-11));
    CHECK(map_dist_acf(a, b) == doctest::Approx(map_dist_acf(b, a)).epsilon(1e-10));
    // map_dist at lag one integrates the same joint density as map_dist_lag1
    CHECK(map_dist(a, b, 1u) == doctest::Approx(map_dist_lag1(a, b)).epsilon(1e-10));
    CHECK(map_feastol() == 8);
}

TEST_CASE("the continuous MAP distances agree with the MATLAB reference") {
    const Map<double> a = ref_map2(), b = ref_map3();
    CHECK(map_exp_mul_int(a, a, 1u) == doctest::Approx(1.3887249114522).epsilon(1e-11));
    CHECK(map_exp_mul_int(a, a, 3u) == doctest::Approx(2.6811012024799).epsilon(1e-11));
    CHECK(map_exp_mul_int(a, b, 3u) == doctest::Approx(2.9699500254906).epsilon(1e-11));
    CHECK(map_geo_mul_sum(a, a) == doctest::Approx(0.018059969121313).epsilon(1e-10));
    CHECK(map_geo_mul_sum(a, b) == doctest::Approx(0.022201392972951).epsilon(1e-10));
    CHECK(map_dist(a, b, 3u) == doctest::Approx(0.30170162728094).epsilon(1e-10));
    CHECK(map_dist_acf(a, b) == doctest::Approx(6.8724997101852e-05).epsilon(1e-8));
    CHECK(map_dist_lag1(a, b) == doctest::Approx(0.060843661413429).epsilon(1e-10));
}

TEST_CASE("map_gamma is exact below order three and fitted above it") {
    // Poisson has no correlation at all
    CHECK(map_gamma(map_exponential(2.0)) == doctest::Approx(0.0).epsilon(1e-14));
    // An order-two MAP has a geometric acf, so the rate is the exact ratio
    const Map<double> a = ref_map2();
    CHECK(map_gamma(a) == doctest::Approx(0.043478260869565).epsilon(1e-9));
    CHECK(map_gamma(a) == doctest::Approx(map_gamma2(a).real()).epsilon(1e-9));
    // A renewal PH has no lag-one correlation and so a zero decay rate
    const Map<double> ph{mk(2, 2, {-3.0, 3.0, 0.0, -1.0}), mk(2, 2, {0.0, 0.0, 1.0, 0.0})};
    CHECK(map_gamma(ph) == doctest::Approx(0.0).epsilon(1e-14));
    // Above order two rho0 is set by the scv and the rate is a genuine fit
    const MapGammaResult<double> r = map_gamma_full(ref_map3());
    const double scv = map_scv(ref_map3());
    CHECK(r.rho0 == doctest::Approx(0.5 * (1.0 - 1.0 / scv)).epsilon(1e-12));
    CHECK(r.gamma < 1.0);
}
