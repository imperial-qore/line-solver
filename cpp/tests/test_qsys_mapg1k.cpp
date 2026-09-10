/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * MAP/G/1/K with tail drop and the two per-class layers over it.
 *
 * Oracles, in the order the task prescribes.
 *  (a) Closed forms the models collapse to. M/M/1/K has the exact geometric
 *      law p_l = (1-r) r^l/(1-r^(K+1)); M/M/1/1 has Erlang's loss formula
 *      B(a,1) = a/(1+a).
 *  (b) Invariants, asserted on every instance: sum_l plevel(l) = 1;
 *      p0 = 1 - T S (Markov renewal reward over a cycle); lossProbability
 *      = 1 - T/lambda; E[N] = sum_l l plevel(l). For POISSON arrivals PASTA
 *      additionally forces pK = lossProbability, which the port satisfies
 *      without ever invoking PASTA -- it is a consequence of the construction
 *      here, not an assumption in it. For qsys_mmapg1k, splitting one MAP
 *      proportionally across classes must give IDENTICAL per-class loss
 *      ratios, since numerator and denominator both scale by the split weight;
 *      that is the check that the phase resolution is doing real work rather
 *      than reproducing the aggregate. For qsys_mapg1k_perflow at N = 1 the
 *      Poisson background is empty and the result must reduce exactly to
 *      qsys_mapg1k.
 *  (c) MATLAB, to 1e-12 relative. All four service-law paths (gamma, det, PH,
 *      density) are also cross-checked against each other on laws they can
 *      both express, which exercises four independent c_n generators against
 *      one answer.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/qsys/qsys_mapg1k.h"
#include "line/api/qsys/qsys_mapg1k_perflow.h"
#include "line/api/qsys/qsys_mmapg1k.h"

using line::Matrix;
using line::Real50;
using line::mam::Map;
using line::qsys::ServiceLaw;
using line::qsys::qsys_mapg1k;
using line::qsys::qsys_mapg1k_perflow;
using line::qsys::qsys_mmapg1k;

namespace {

template <class T>
Matrix<T> mat(const std::vector<std::vector<double>>& a) {
    Matrix<T> m(a.size(), a[0].size());
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < a[0].size(); ++j)
            m(i, j) = line::num_traits<T>::from_double(a[i][j]);
    return m;
}

template <class T>
Map<T> mkmap(const std::vector<std::vector<double>>& d0,
             const std::vector<std::vector<double>>& d1) {
    Map<T> m;
    m.D0 = mat<T>(d0);
    m.D1 = mat<T>(d1);
    return m;
}

template <class T>
Map<T> poisson(double lambda) {
    return mkmap<T>({{-lambda}}, {{lambda}});
}

/** Correlated MMPP2, lambda = 7/6. */
template <class T>
Map<T> corr_map() {
    return mkmap<T>({{-2.5, 0.2}, {0.1, -0.7}}, {{2.3, 0.0}, {0.0, 0.6}});
}

/** Every identity the construction guarantees, on any instance. */
template <class T>
void check_identities(const line::qsys::MapG1kResult<T>& r, bool poisson_arrivals) {
    T mass = line::num_traits<T>::from_int(0);
    for (const T& v : r.plevel) {
        CHECK(static_cast<double>(v) >= -1e-14);
        mass += v;
    }
    CHECK(static_cast<double>(mass) == doctest::Approx(1.0).epsilon(1e-11));
    // p0 = 1 - T S: the empty fraction of an inter-departure cycle.
    CHECK(static_cast<double>(r.p0) ==
          doctest::Approx(1.0 - static_cast<double>(r.throughput * r.meanServiceTime))
              .epsilon(1e-11));
    CHECK(static_cast<double>(r.lossProbability) ==
          doctest::Approx(1.0 - static_cast<double>(r.throughput / r.lambda)).epsilon(1e-11));
    CHECK(static_cast<double>(r.utilization) ==
          doctest::Approx(1.0 - static_cast<double>(r.p0)).epsilon(1e-14));
    double meanQ = 0.0, pKsum = 0.0;
    for (std::size_t l = 0; l < r.plevel.size(); ++l)
        meanQ += static_cast<double>(l) * static_cast<double>(r.plevel[l]);
    for (const T& v : r.pKvec) pKsum += static_cast<double>(v);
    CHECK(static_cast<double>(r.meanQueueLength) == doctest::Approx(meanQ).epsilon(1e-12));
    CHECK(pKsum == doctest::Approx(static_cast<double>(r.pK)).epsilon(1e-12));
    CHECK(static_cast<double>(r.plevel.back()) ==
          doctest::Approx(static_cast<double>(r.pK)).epsilon(1e-12));
    if (poisson_arrivals) {
        // PASTA: an arrival of a Poisson stream sees the time-stationary law,
        // so the fraction it finds full is the loss ratio. Never assumed by
        // the algorithm, so this is a real check.
        CHECK(static_cast<double>(r.pK) ==
              doctest::Approx(static_cast<double>(r.lossProbability)).epsilon(1e-11));
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// (a) closed-form collapses
// ---------------------------------------------------------------------------

TEST_CASE("M/M/1/K reproduces the truncated geometric law") {
    // lambda = 2, mu = 3, r = 2/3, K = 5.
    const double rr = 2.0 / 3.0;
    const std::size_t K = 5;
    const auto r = qsys_mapg1k(poisson<double>(2.0), ServiceLaw<double>::gamma(1.0, 1.0 / 3.0), K);
    double den = 0.0;
    for (std::size_t l = 0; l <= K; ++l) den += std::pow(rr, static_cast<double>(l));
    for (std::size_t l = 0; l <= K; ++l)
        CHECK(r.plevel[l] == doctest::Approx(std::pow(rr, static_cast<double>(l)) / den)
                                 .epsilon(1e-12));
    CHECK(r.pK == doctest::Approx(std::pow(rr, 5.0) / den).epsilon(1e-12));
    // The reference's own documented test value.
    CHECK(r.pK == doctest::Approx(0.04812030).epsilon(1e-7));
    check_identities(r, true);
}

TEST_CASE("M/M/1/1 is Erlang's loss formula") {
    // K = 1: the buffer holds only the packet in service, so the loss ratio is
    // B(a, 1) = a/(1+a) with a = lambda S = 2/3, i.e. exactly 0.4.
    const auto r = qsys_mapg1k(poisson<double>(2.0), ServiceLaw<double>::gamma(1.0, 1.0 / 3.0),
                               static_cast<std::size_t>(1));
    CHECK(r.lossProbability == doctest::Approx(0.4).epsilon(1e-13));
    CHECK(r.p0 == doctest::Approx(0.6).epsilon(1e-13));
    CHECK(r.pK == doctest::Approx(0.4).epsilon(1e-13));
    CHECK(r.meanQueueLength == doctest::Approx(0.4).epsilon(1e-13));
    check_identities(r, true);
}

// ---------------------------------------------------------------------------
// (b) the four service-law paths against each other
// ---------------------------------------------------------------------------

TEST_CASE("the four service-law generators agree where they overlap") {
    const Map<double> P = poisson<double>(2.0);
    const std::size_t K = 5;
    SUBCASE("Erlang-2 of mean 1/3, as a gamma and as a phase type") {
        const auto g = qsys_mapg1k(P, ServiceLaw<double>::gamma(2.0, 1.0 / 6.0), K);
        Matrix<double> Tm = mat<double>({{-6.0, 6.0}, {0.0, -6.0}});
        const auto p = qsys_mapg1k(P, ServiceLaw<double>::phase_type({1.0, 0.0}, Tm), K);
        CHECK(g.pK == doctest::Approx(p.pK).epsilon(1e-13));
        CHECK(g.meanQueueLength == doctest::Approx(p.meanQueueLength).epsilon(1e-13));
        CHECK(g.p0 == doctest::Approx(p.p0).epsilon(1e-13));
        check_identities(g, true);
        check_identities(p, true);
    }
    SUBCASE("exponential of mean 1/3, as a gamma and as a density") {
        const auto g = qsys_mapg1k(P, ServiceLaw<double>::gamma(1.0, 1.0 / 3.0), K);
        const auto d = qsys_mapg1k(
            P, ServiceLaw<double>::density([](const double& x) { return 3.0 * std::exp(-3.0 * x); }),
            K);
        CHECK(d.meanServiceTime == doctest::Approx(1.0 / 3.0).epsilon(1e-12));
        CHECK(d.pK == doctest::Approx(g.pK).epsilon(1e-11));
        CHECK(d.meanQueueLength == doctest::Approx(g.meanQueueLength).epsilon(1e-11));
        check_identities(d, true);
    }
    SUBCASE("deterministic service queues less than exponential at the same mean") {
        const auto e = qsys_mapg1k(P, ServiceLaw<double>::gamma(1.0, 1.0 / 3.0), K);
        const auto d = qsys_mapg1k(P, ServiceLaw<double>::deterministic(1.0 / 3.0), K);
        const auto e2 = qsys_mapg1k(P, ServiceLaw<double>::gamma(2.0, 1.0 / 6.0), K);
        // Loss falls as the service becomes less variable: Exp > Erlang-2 > Det.
        CHECK(d.pK < e2.pK);
        CHECK(e2.pK < e.pK);
        CHECK(d.meanQueueLength < e2.meanQueueLength);
        CHECK(e2.meanQueueLength < e.meanQueueLength);
        check_identities(d, true);
    }
}

// ---------------------------------------------------------------------------
// (c) MATLAB
// ---------------------------------------------------------------------------

TEST_CASE("MAP/G/1/K agrees with MATLAB") {
    SUBCASE("M/M/1/5") {
        const auto r =
            qsys_mapg1k(poisson<double>(2.0), ServiceLaw<double>::gamma(1.0, 1.0 / 3.0),
                        static_cast<std::size_t>(5));
        CHECK(r.p0 == doctest::Approx(0.365413533834586).epsilon(1e-13));
        CHECK(r.pK == doctest::Approx(0.0481203007518798).epsilon(1e-13));
        CHECK(r.throughput == doctest::Approx(1.90375939849624).epsilon(1e-13));
        CHECK(r.meanQueueLength == doctest::Approx(1.42255639097744).epsilon(1e-13));
        CHECK(r.nmax == 43u);
    }
    SUBCASE("M/D/1/5") {
        const auto r = qsys_mapg1k(poisson<double>(2.0), ServiceLaw<double>::deterministic(1.0 / 3.0),
                                   static_cast<std::size_t>(5));
        CHECK(r.p0 == doctest::Approx(0.341703649752645).epsilon(1e-13));
        CHECK(r.pK == doctest::Approx(0.0125554746289671).epsilon(1e-13));
        CHECK(r.meanQueueLength == doctest::Approx(1.20833898140121).epsilon(1e-13));
        const std::vector<double> plevel = {0.341703649752645, 0.323844180823206,
                                            0.18706378126959,  0.0917417892083818,
                                            0.0430911243172099, 0.0125554746289671};
        for (std::size_t l = 0; l < plevel.size(); ++l)
            CHECK(r.plevel[l] == doctest::Approx(plevel[l]).epsilon(1e-12));
        check_identities(r, true);
    }
    SUBCASE("M/E2/1/5 through the phase-type path") {
        Matrix<double> Tm = mat<double>({{-6.0, 6.0}, {0.0, -6.0}});
        const auto r = qsys_mapg1k(poisson<double>(2.0),
                                   ServiceLaw<double>::phase_type({1.0, 0.0}, Tm),
                                   static_cast<std::size_t>(5));
        CHECK(r.p0 == doctest::Approx(0.353280086152742).epsilon(1e-13));
        CHECK(r.pK == doctest::Approx(0.0299201292291125).epsilon(1e-13));
        CHECK(r.meanQueueLength == doctest::Approx(1.33840078973347).epsilon(1e-13));
    }
    SUBCASE("MMPP2/gamma/1/8, a shape below one so the density is singular at 0") {
        const auto r = qsys_mapg1k(corr_map<double>(), ServiceLaw<double>::gamma(0.25, 2.4),
                                   static_cast<std::size_t>(8));
        CHECK(r.p0 == doctest::Approx(0.406115208568551).epsilon(1e-12));
        CHECK(r.pK == doctest::Approx(0.103591669981134).epsilon(1e-12));
        CHECK(r.throughput == doctest::Approx(0.989807985719082).epsilon(1e-12));
        CHECK(r.lossProbability == doctest::Approx(0.151593155097929).epsilon(1e-12));
        CHECK(r.meanQueueLength == doctest::Approx(2.45026046661011).epsilon(1e-12));
        CHECK(r.rho == doctest::Approx(0.7).epsilon(1e-13));
        CHECK(r.nmax == 174u);
        CHECK(r.pKvec[0] == doctest::Approx(0.067472752326508).epsilon(1e-11));
        CHECK(r.pKvec[1] == doctest::Approx(0.0361189176546256).epsilon(1e-11));
        // Correlated arrivals: pK and the loss ratio are NOT equal, which is
        // precisely why the phase-resolved pKvec exists.
        CHECK(std::fabs(static_cast<double>(r.pK - r.lossProbability)) > 0.04);
        check_identities(r, false);
    }
}

TEST_CASE("MAP/G/1/K rejects a malformed instance") {
    CHECK_THROWS_AS(qsys_mapg1k(poisson<double>(2.0), ServiceLaw<double>::gamma(1.0, 1.0 / 3.0),
                                static_cast<std::size_t>(0)),
                    line::InputError);
    CHECK_THROWS_AS(qsys_mapg1k(poisson<double>(2.0), ServiceLaw<double>::gamma(-1.0, 1.0),
                                static_cast<std::size_t>(3)),
                    line::InputError);
    CHECK_THROWS_AS(qsys_mapg1k(poisson<double>(2.0), ServiceLaw<double>::deterministic(-1.0),
                                static_cast<std::size_t>(3)),
                    line::InputError);
}

TEST_CASE("MAP/G/1/K at Real50 reproduces its own double result") {
    const auto rd = qsys_mapg1k(poisson<double>(2.0), ServiceLaw<double>::gamma(1.0, 1.0 / 3.0),
                                static_cast<std::size_t>(5));
    const auto rq = qsys_mapg1k(poisson<Real50>(2.0),
                                ServiceLaw<Real50>::gamma(Real50(1), Real50(1) / Real50(3)),
                                static_cast<std::size_t>(5));
    CHECK(static_cast<double>(rq.pK) == doctest::Approx(rd.pK).epsilon(1e-13));
    CHECK(static_cast<double>(rq.meanQueueLength) ==
          doctest::Approx(rd.meanQueueLength).epsilon(1e-13));
}

// ---------------------------------------------------------------------------
// qsys_mmapg1k
// ---------------------------------------------------------------------------

TEST_CASE("MMAP/G/1/K tells two classes of the same rate apart") {
    const Matrix<double> D0 = mat<double>({{-2.5, 0.2}, {0.1, -0.7}});
    const std::vector<Matrix<double>> D1c = {mat<double>({{2.3, 0.0}, {0.0, 0.0}}),
                                             mat<double>({{0.0, 0.0}, {0.0, 0.6}})};
    const auto r = qsys_mmapg1k(D0, D1c, ServiceLaw<double>::gamma(0.25, 2.4),
                                static_cast<std::size_t>(8));
    CHECK(r.lambda[0] == doctest::Approx(0.766666666666667).epsilon(1e-13));
    CHECK(r.lambda[1] == doctest::Approx(0.4).epsilon(1e-13));
    CHECK(r.throughput[0] == doctest::Approx(0.611479336315698).epsilon(1e-12));
    CHECK(r.throughput[1] == doctest::Approx(0.378328649407225).epsilon(1e-12));
    CHECK(r.lossRatio[0] == doctest::Approx(0.202418256979524).epsilon(1e-12));
    CHECK(r.lossRatio[1] == doctest::Approx(0.0541783764819384).epsilon(1e-12));
    CHECK(r.lossAggregate == doctest::Approx(0.151593155097929).epsilon(1e-12));
    // The two classes differ by a factor of 3.7 in loss ratio: an
    // aggregate-only analysis would report 0.1516 for both.
    CHECK(r.lossRatio[0] / r.lossRatio[1] > 3.0);
    CHECK(r.lambdaAggregate == doctest::Approx(7.0 / 6.0).epsilon(1e-13));
}

TEST_CASE("a proportional class split must give identical loss ratios") {
    // D1c_k = w_k D1 scales numerator and denominator of L_k by the same w_k,
    // so the phase resolution cannot manufacture a difference where the
    // classes are statistically identical. This is the control for the test
    // above.
    const Matrix<double> D0 = mat<double>({{-2.5, 0.2}, {0.1, -0.7}});
    const Matrix<double> D1 = mat<double>({{2.3, 0.0}, {0.0, 0.6}});
    std::vector<Matrix<double>> D1c;
    const std::vector<double> w = {0.25, 0.35, 0.4};
    for (double wk : w) {
        Matrix<double> Dk = D1;
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 2; ++j) Dk(i, j) *= wk;
        D1c.push_back(Dk);
    }
    const auto r = qsys_mmapg1k(D0, D1c, ServiceLaw<double>::gamma(0.25, 2.4),
                                static_cast<std::size_t>(8));
    for (std::size_t k = 1; k < w.size(); ++k)
        CHECK(r.lossRatio[k] == doctest::Approx(r.lossRatio[0]).epsilon(1e-13));
    for (std::size_t k = 0; k < w.size(); ++k) {
        CHECK(r.lambda[k] == doctest::Approx(w[k] * 7.0 / 6.0).epsilon(1e-13));
        CHECK(r.throughput[k] ==
              doctest::Approx(w[k] * r.throughputAggregate).epsilon(1e-12));
    }
}

// ---------------------------------------------------------------------------
// qsys_mapg1k_perflow
// ---------------------------------------------------------------------------

TEST_CASE("per-flow analysis reduces exactly to the single-MAP model at N = 1") {
    // With one flow the Poisson background has rate zero, so T_1 = (1-p0)/S is
    // the aggregate throughput and L_1 is the aggregate loss ratio.
    const std::vector<Map<double>> one = {corr_map<double>()};
    const auto p = qsys_mapg1k_perflow(one, ServiceLaw<double>::gamma(2.0, 0.1),
                                       static_cast<std::size_t>(10));
    const auto r = qsys_mapg1k(corr_map<double>(), ServiceLaw<double>::gamma(2.0, 0.1),
                               static_cast<std::size_t>(10));
    CHECK(p.throughput[0] == doctest::Approx(r.throughput).epsilon(1e-12));
    CHECK(p.lossRatio[0] == doctest::Approx(r.lossProbability).epsilon(1e-11));
    CHECK(p.p0[0] == doctest::Approx(r.p0).epsilon(1e-13));
    CHECK(p.pK[0] == doctest::Approx(r.pK).epsilon(1e-13));
}

TEST_CASE("per-flow analysis agrees with MATLAB on three dissimilar flows") {
    // Poisson of rate 1, Erlang-2 of rate 1.5, and a correlated MMPP2 of rate
    // 7/6, sharing an Erlang-2 service of mean 0.2 into a buffer of 10.
    const std::vector<Map<double>> flows = {poisson<double>(1.0),
                                            mkmap<double>({{-3.0, 3.0}, {0.0, -3.0}},
                                                          {{0.0, 0.0}, {3.0, 0.0}}),
                                            corr_map<double>()};
    const auto p = qsys_mapg1k_perflow(flows, ServiceLaw<double>::gamma(2.0, 0.1),
                                       static_cast<std::size_t>(10));
    CHECK(p.lambda[0] == doctest::Approx(1.0).epsilon(1e-14));
    CHECK(p.lambda[1] == doctest::Approx(1.5).epsilon(1e-14));
    CHECK(p.lambda[2] == doctest::Approx(7.0 / 6.0).epsilon(1e-13));
    CHECK(p.throughput[0] == doctest::Approx(0.994689717721069).epsilon(1e-12));
    CHECK(p.throughput[1] == doctest::Approx(1.49628566968804).epsilon(1e-12));
    CHECK(p.throughput[2] == doctest::Approx(1.14152406707092).epsilon(1e-12));
    CHECK(p.lossRatio[0] == doctest::Approx(0.00531028227893149).epsilon(1e-10));
    CHECK(p.lossRatio[1] == doctest::Approx(0.00247622020797245).epsilon(1e-10));
    CHECK(p.lossRatio[2] == doctest::Approx(0.021550799653497).epsilon(1e-10));
    CHECK(p.rho == doctest::Approx(0.733333333333333).epsilon(1e-13));
    CHECK(p.lossAggregate == doctest::Approx(0.00931833059635547).epsilon(1e-10));
    // The correlated flow loses about four times as much as the smooth
    // Erlang-2 flow at the same buffer: the effect the method exists to show.
    CHECK(p.lossRatio[2] / p.lossRatio[1] > 4.0);
    // The Poisson flow's own model reproduces PASTA: its loss ratio is the
    // full-buffer probability of the model it was measured in.
    CHECK(p.lossRatio[0] == doctest::Approx(p.pK[0]).epsilon(1e-10));
}
