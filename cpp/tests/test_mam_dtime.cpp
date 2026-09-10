/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Discrete-time (slotted) matrix-analytic queues. Oracles: the Geo/Geo/1
 * closed form under LAS_DA, the MATLAB SolverMAM numbers recorded in
 * _kb/06-solver-catalog.md (which LDES slotted independently corroborated),
 * and the agreement of the two solution routes on the same queue, since the
 * QBD path uses logarithmic reduction while the batch path truncates the
 * M/G/1-type level space.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/mam/dtime.h"

using line::Matrix;
using line::lang::ProcessType;
using line::mam::DBatch;
using line::mam::Dmap;
using line::mam::Dph;
using line::mam::dmap_lambda_batch;
using line::mam::dmap_super;
using line::mam::dmap_thin;
using line::mam::dph_from_dist;
using line::mam::mg1_dt_queue;
using line::mam::q_dt_map_map_1;
using line::mam::q_dt_ph_ph_1;

namespace {

Matrix<double> scalar(double v) { return Matrix<double>(1, 1, v); }

double mean_of(const std::vector<double>& ql) {
    double q = 0;
    for (std::size_t i = 0; i < ql.size(); ++i) q += static_cast<double>(i) * ql[i];
    return q;
}

}  // namespace

TEST_CASE("Geo/Geo/1 reproduces the LAS-DA closed form") {
    const double a = 0.2, s = 0.5;
    Dmap<double> arv{scalar(1 - a), scalar(a)};
    Dmap<double> svc{scalar(1 - s), scalar(s)};
    std::vector<double> ql = q_dt_map_map_1(arv, svc);

    // 1e-7, not tighter: the level sum stops once the geometric tail carries
    // less than 1e-10, which is worth ~3e-9 of relative error on the mean
    CHECK(mean_of(ql) == doctest::Approx(a * (1 - a) / (s - a)).epsilon(1e-7));
    CHECK(ql[0] == doctest::Approx(1 - a / s).epsilon(1e-7));
}

TEST_CASE("Det and DiscreteUniform service stay on the slot lattice") {
    // MATLAB SolverMAM: 0.466667 and 0.465000; LDES slotted measured 0.466615
    Dph<double> geo = dph_from_dist<double>(ProcessType::GEOMETRIC, 1.0 / 0.2, 0.8);
    Dph<double> det = dph_from_dist<double>(ProcessType::DET, 2.0, 0.0);
    CHECK(mean_of(q_dt_ph_ph_1(geo, det)) == doctest::Approx(0.466667).epsilon(1e-5));

    Dph<double> geo2 = dph_from_dist<double>(ProcessType::GEOMETRIC, 1.0 / 0.15, 0.85);
    Dph<double> du = dph_from_dist<double>(ProcessType::DUNIFORM, 2.5, 0.2);
    CHECK(mean_of(q_dt_ph_ph_1(geo2, du)) == doctest::Approx(0.465000).epsilon(1e-5));
}

TEST_CASE("DMAP arrivals are solved on the slotted time scale") {
    // MATLAB SolverMAM 0.700000, LDES slotted 0.698332
    Matrix<double> D0(2, 2), D1(2, 2);
    D0(0, 0) = 0.5; D0(0, 1) = 0.2; D0(1, 0) = 0.1; D0(1, 1) = 0.6;
    D1(0, 0) = 0.25; D1(0, 1) = 0.05; D1(1, 0) = 0.1; D1(1, 1) = 0.2;
    Dmap<double> arv{D0, D1};
    Dmap<double> svc{scalar(0.4), scalar(0.6)};

    std::vector<double> ql = q_dt_map_map_1(arv, svc);
    CHECK(mean_of(ql) == doctest::Approx(0.700000).epsilon(1e-5));
    CHECK(1 - ql[0] == doctest::Approx(0.5).epsilon(1e-7));
}

TEST_CASE("the truncated M/G/1-type route agrees with the QBD route") {
    const double a = 0.2, s = 0.5;
    DBatch<double> batch;
    batch.push_back(scalar(1 - a));
    batch.push_back(scalar(a));
    Dmap<double> svc{scalar(1 - s), scalar(s)};

    auto r = mg1_dt_queue(batch, svc, 200, true);
    CHECK(r.QN == doctest::Approx(a * (1 - a) / (s - a)).epsilon(1e-7));
    CHECK(r.UN == doctest::Approx(a / s).epsilon(1e-7));
    CHECK(r.TN == doctest::Approx(a).epsilon(1e-12));
    CHECK(r.dep.D0.rows() == r.dep.D1.rows());
}

TEST_CASE("superposing two slotted streams produces batches, thinning halves the rate") {
    // Two Bernoulli(0.1) streams merge into batch masses 0.81 / 0.18 / 0.01:
    // the middle term is what a DMAP-only representation would have to discard.
    DBatch<double> b;
    b.push_back(scalar(0.9));
    b.push_back(scalar(0.1));
    DBatch<double> merged = dmap_super(b, b);

    REQUIRE(merged.size() == 3);
    CHECK(merged[0](0, 0) == doctest::Approx(0.81).epsilon(1e-12));
    CHECK(merged[1](0, 0) == doctest::Approx(0.18).epsilon(1e-12));
    CHECK(merged[2](0, 0) == doctest::Approx(0.01).epsilon(1e-12));
    CHECK(dmap_lambda_batch(merged) == doctest::Approx(0.2).epsilon(1e-12));
    CHECK(dmap_lambda_batch(dmap_thin(merged, 0.5)) == doctest::Approx(0.1).epsilon(1e-12));
}

TEST_CASE("an overloaded slotted station is refused rather than reported") {
    Dmap<double> arv{scalar(0.4), scalar(0.6)};
    Dmap<double> svc{scalar(0.7), scalar(0.3)};
    CHECK_THROWS(q_dt_map_map_1(arv, svc));
}
