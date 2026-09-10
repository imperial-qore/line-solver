/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The LQN parameter-identification driver stack of api/infer: the parameter
 * injection and read-back of infer_lqn_setparams and the Extended Kalman
 * Filter driver of infer_lqn.
 *
 * ORACLES.
 *  (a) Round trip: what infer_lqn_setparams writes, infer_lqn_getparams must
 *      read back exactly, and the injected law must be the exponential of the
 *      requested mean, since that is what the CASCON 2005 method assumes.
 *  (b) Identifiability on a KNOWN observation model. The driver is exercised
 *      with an analytic h(a) rather than a solver, so the estimate can be
 *      checked against the parameter that generated the measurements: a
 *      filter that tracks must drive the estimate towards the truth and the
 *      prediction error towards zero, and neither is true of a filter whose
 *      Q, R or P0 construction is wrong.
 *  (c) The covariance defaults of the paper, equations (9a) and (9b),
 *      reproduced entry by entry from the option factors.
 */
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/infer/infer_lqn.h"
#include "line/api/infer/infer_lqn_setparams.h"
#include "line/lang/lqn/lqn_struct.h"

namespace infer = line::infer;
namespace lqn = line::lqn;
using line::Matrix;
using line::lang::LqnElement;

namespace {

/**
 * A two-task, two-activity layered struct: one reference task with a think
 * time and one server task whose two activities carry host demands. Only the
 * fields the parameter stack reads are populated, which is what a struct
 * handed to a getter needs.
 */
lqn::LqnStruct<double> small_lsn() {
    lqn::LqnStruct<double> s;
    s.nhosts = 1;
    s.ntasks = 2;
    s.nentries = 0;
    s.nacts = 2;
    s.nidx = 5;
    s.hshift = 0;
    s.tshift = 1;
    s.names.resize(s.nidx + 1);
    s.hashnames.resize(s.nidx + 1);
    s.type.resize(s.nidx + 1, LqnElement::HOST);
    s.parent.assign(s.nidx + 1, 0);
    s.hostdem.assign(s.nidx + 1, line::lang::Distrib<double>::exp_mean(0.0));
    s.think.assign(s.nidx + 1, line::lang::Distrib<double>::exp_mean(0.0));
    s.actthink.assign(s.nidx + 1, line::lang::Distrib<double>::exp_mean(0.0));

    s.names[1] = "P1";
    s.type[1] = LqnElement::HOST;
    s.names[2] = "Client";
    s.type[2] = LqnElement::TASK;
    s.parent[2] = 1;
    s.think[2] = line::lang::Distrib<double>::exp_mean(1.5);
    s.names[3] = "Server";
    s.type[3] = LqnElement::TASK;
    s.parent[3] = 1;
    s.names[4] = "actA";
    s.type[4] = LqnElement::ACTIVITY;
    s.parent[4] = 3;
    s.hostdem[4] = line::lang::Distrib<double>::exp_mean(0.4);
    s.names[5] = "actB";
    s.type[5] = LqnElement::ACTIVITY;
    s.parent[5] = 3;
    s.hostdem[5] = line::lang::Distrib<double>::exp_mean(0.9);
    return s;
}

}  // namespace

TEST_CASE("infer_lqn_setparams injects and reads back the named parameters") {
    lqn::LqnStruct<double> s = small_lsn();
    std::vector<infer::LqnParamSpec> spec(3);
    spec[0].type = infer::LqnParamType::HOSTDEM;
    spec[0].name = "actA";
    spec[1].type = infer::LqnParamType::HOSTDEM;
    spec[1].name = "actB";
    spec[2].type = infer::LqnParamType::THINK;
    spec[2].name = "Client";

    // Oracle (a): read back the values the struct was built with
    const std::vector<double> a0 = infer::infer_lqn_getparams(s, spec);
    REQUIRE(a0.size() == 3);
    CHECK(a0[0] == doctest::Approx(0.4));
    CHECK(a0[1] == doctest::Approx(0.9));
    CHECK(a0[2] == doctest::Approx(1.5));

    std::vector<double> a(3);
    a[0] = 0.75;
    a[1] = 0.25;
    a[2] = 3.0;
    infer::infer_lqn_setparams(s, spec, a);
    const std::vector<double> back = infer::infer_lqn_getparams(s, spec);
    for (std::size_t i = 0; i < 3; ++i) CHECK(back[i] == doctest::Approx(a[i]).epsilon(1e-14));

    // Oracle (a): the injected law is the EXPONENTIAL of that mean, SCV 1
    CHECK(s.hostdem[4].type == line::lang::ProcessType::EXP);
    CHECK(s.hostdem[4].scv == doctest::Approx(1.0));
    CHECK(s.hostdem[4].rate() == doctest::Approx(1.0 / 0.75).epsilon(1e-12));
    CHECK(s.think[2].type == line::lang::ProcessType::EXP);
    CHECK(s.think[2].mean == doctest::Approx(3.0));
    // the untouched entries are left alone
    CHECK(s.think[3].mean == doctest::Approx(0.0));

    // the kind resolves a name that could denote either: "Server" is a TASK,
    // so a HOSTDEM request for it must fail rather than silently hit the task
    std::vector<infer::LqnParamSpec> bad(1);
    bad[0].type = infer::LqnParamType::HOSTDEM;
    bad[0].name = "Server";
    CHECK_THROWS(infer::infer_lqn_getparams(s, bad));
    bad[0].type = infer::LqnParamType::THINK;
    bad[0].name = "actA";
    CHECK_THROWS(infer::infer_lqn_getparams(s, bad));
    // a wrong-length parameter vector is refused
    CHECK_THROWS(infer::infer_lqn_setparams(s, spec, std::vector<double>(2, 1.0)));
}

TEST_CASE("infer_lqn tracks a known observation model") {
    lqn::LqnStruct<double> s = small_lsn();
    std::vector<infer::LqnParamSpec> spec(2);
    spec[0].type = infer::LqnParamType::HOSTDEM;
    spec[0].name = "actA";
    spec[1].type = infer::LqnParamType::HOSTDEM;
    spec[1].name = "actB";

    std::vector<infer::LqnObsSpec> obs(2);
    obs[0].metric = infer::LqnMetric::RespT;
    obs[0].name = "actA";
    obs[1].metric = infer::LqnMetric::RespT;
    obs[1].name = "actB";

    // a smooth, invertible surrogate for the layered solve: each activity's
    // response time inflates its own demand and is perturbed by the other, so
    // the sensitivity matrix is genuinely non-diagonal and a filter that only
    // ever moves one coordinate at a time cannot pass this test
    std::function<infer::LqnMetrics<double>(const lqn::LqnStruct<double>&)> solve =
        [](const lqn::LqnStruct<double>& m) {
            infer::LqnMetrics<double> out;
            out.RespT.assign(m.names.size(), 0.0);
            out.QLen.assign(m.names.size(), 0.0);
            out.Util.assign(m.names.size(), 0.0);
            out.Tput.assign(m.names.size(), 0.0);
            const double da = m.hostdem[4].mean, db = m.hostdem[5].mean;
            out.RespT[4] = da * (1.0 + 0.5 * db);
            out.RespT[5] = db * (1.0 + 0.8 * da);
            return out;
        };

    // measurements generated by the TRUE parameters, held constant over 12 steps
    const double trueA = 0.62, trueB = 0.31;
    const double zA = trueA * (1.0 + 0.5 * trueB);
    const double zB = trueB * (1.0 + 0.8 * trueA);
    const std::size_t nsteps = 12;
    Matrix<double> Z(2, nsteps, 0.0);
    for (std::size_t k = 0; k < nsteps; ++k) {
        Z(0, k) = zA;
        Z(1, k) = zB;
    }

    infer::InferLqnOptions<double> opt;
    opt.a_true.push_back(trueA);
    opt.a_true.push_back(trueB);
    const infer::InferLqnResult<double> r = infer::infer_lqn(s, spec, obs, Z, solve, opt);

    // Oracle (c): the defaults are the paper's, entry by entry
    REQUIRE(r.a0.size() == 2);
    CHECK(r.a0[0] == doctest::Approx(0.4));
    CHECK(r.a0[1] == doctest::Approx(0.9));
    for (std::size_t i = 0; i < 2; ++i) {
        const double q = 0.1 * 1.0 * r.a0[i];
        CHECK(r.Q(i, i) == doctest::Approx(q * q).epsilon(1e-14));
        const double p = 0.5 * r.a0[i];
        CHECK(r.P0(i, i) == doctest::Approx(p * p).epsilon(1e-14));
    }
    const double rA = (0.2 * zA) / 1.96, rB = (0.2 * zB) / 1.96;
    CHECK(r.R(0, 0) == doctest::Approx(rA * rA).epsilon(1e-14));
    CHECK(r.R(1, 1) == doctest::Approx(rB * rB).epsilon(1e-14));
    // the off-diagonals are zero: the construction is diagonal by hypothesis
    CHECK(r.Q(0, 1) == 0.0);
    CHECK(r.R(0, 1) == 0.0);
    CHECK(r.P0(1, 0) == 0.0);

    // Oracle (b): the filter must converge on the parameters that generated Z
    REQUIRE(r.ekf.ahat.cols() == nsteps);
    const double eA = r.ekf.ahat(0, nsteps - 1), eB = r.ekf.ahat(1, nsteps - 1);
    CHECK(eA == doctest::Approx(trueA).epsilon(1e-3));
    CHECK(eB == doctest::Approx(trueB).epsilon(1e-3));
    // and the estimate must improve monotonically enough that the last step is
    // far closer than the first: a filter with a broken gain drifts instead
    const double err0 = std::abs(r.ekf.ahat(0, 0) - trueA) + std::abs(r.ekf.ahat(1, 0) - trueB);
    const double errN = std::abs(eA - trueA) + std::abs(eB - trueB);
    CHECK(errN < 0.05 * err0);
    // the prediction error collapses with it, but only to the floor the drift
    // covariance sets: Q is re-injected every step, so the residual settles at
    // the process noise rather than at zero
    CHECK(std::abs(r.ekf.e(0, nsteps - 1)) < 1e-2 * std::abs(r.ekf.e(0, 0)));
    CHECK(r.ekf.has_Ea);
    CHECK(r.ekf.Er >= 0.0);

    // Oracle (b): the final estimate is written back to the model
    const std::vector<double> applied = infer::infer_lqn_getparams(s, spec);
    CHECK(applied[0] == doctest::Approx(eA).epsilon(1e-14));
    CHECK(applied[1] == doctest::Approx(eB).epsilon(1e-14));

    // an explicit covariance overrides the constructed one verbatim
    lqn::LqnStruct<double> s2 = small_lsn();
    infer::InferLqnOptions<double> opt2;
    opt2.Q = Matrix<double>(2, 2, 0.0);
    opt2.Q(0, 0) = 7.0;
    opt2.Q(1, 1) = 11.0;
    const infer::InferLqnResult<double> r2 = infer::infer_lqn(s2, spec, obs, Z, solve, opt2);
    CHECK(r2.Q(0, 0) == 7.0);
    CHECK(r2.Q(1, 1) == 11.0);

    // gammaT scales R down: a longer measurement interval is more informative
    lqn::LqnStruct<double> s3 = small_lsn();
    infer::InferLqnOptions<double> opt3;
    opt3.has_T = true;
    opt3.T_interval = 4.0;
    opt3.has_Tstar = true;
    opt3.Tstar = 1.0;
    const infer::InferLqnResult<double> r3 = infer::infer_lqn(s3, spec, obs, Z, solve, opt3);
    CHECK(r3.R(0, 0) == doctest::Approx(rA * rA / 4.0).epsilon(1e-14));

    // a Z with the wrong number of rows is refused
    Matrix<double> Zbad(3, nsteps, 1.0);
    lqn::LqnStruct<double> s4 = small_lsn();
    CHECK_THROWS(infer::infer_lqn(s4, spec, obs, Zbad, solve, opt));
}
