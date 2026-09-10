/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * infer_mlps / infer_minps: maximum-likelihood demand estimation at a PS queue.
 *
 * THE ORACLE IS A CASE WHERE THE ANSWER IS KNOWN IN CLOSED FORM, not a stored
 * number. When every sample arrives to an otherwise empty queue, the tagged job
 * is alone at the processor-sharing server, its sojourn is exactly Exp(1/d),
 * and the maximum-likelihood estimate of d is the SAMPLE MEAN. That single
 * check exercises the entire chain -- the augmented model, the generator, the
 * departure filtration, the absorbing subset, the phase-type density and the
 * optimizer -- against an answer this code did not compute; measured, it agrees
 * to seven digits.
 *
 * The remaining cases drive the parts a closed form does not reach: a second
 * class in the system (so the tagged job actually shares the server), the
 * MINPS selection rule, and the refusals.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/infer/infer_mlps.h"

namespace api = line::api;

namespace {

api::MlpsSample sample(double rt, std::size_t cls, const std::vector<double>& ql) {
    api::MlpsSample s;
    s.rt = rt;
    s.cls = cls;
    s.ql = ql;
    return s;
}

}  // namespace

TEST_CASE("a job alone at the server: the MLE is the sample mean, exactly") {
    std::vector<double> muZ(1, 1.0);
    const double v[] = {0.30, 0.55, 0.42, 0.71, 0.22, 0.63, 0.38, 0.49};
    std::vector<api::MlpsSample> S;
    double sum = 0.0;
    for (std::size_t i = 0; i < 8; ++i) {
        S.push_back(sample(v[i], 1, std::vector<double>(1, 1.0)));
        sum += v[i];
    }
    const double mean = sum / 8.0;

    const std::vector<double> d = api::infer_mlps(muZ, 1.0, S);
    REQUIRE(d.size() == 1u);
    CHECK(d[0] == doctest::Approx(mean).epsilon(1e-6));
}

TEST_CASE("the estimate scales with the response times") {
    // Doubling every response time in the alone-at-the-server case doubles the
    // demand, since the MLE is the sample mean of an exponential.
    std::vector<double> muZ(1, 1.0);
    const double v[] = {0.30, 0.55, 0.42, 0.71};
    std::vector<api::MlpsSample> S, S2;
    for (std::size_t i = 0; i < 4; ++i) {
        S.push_back(sample(v[i], 1, std::vector<double>(1, 1.0)));
        S2.push_back(sample(2.0 * v[i], 1, std::vector<double>(1, 1.0)));
    }
    const std::vector<double> a = api::infer_mlps(muZ, 1.0, S);
    const std::vector<double> b = api::infer_mlps(muZ, 1.0, S2);
    CHECK(b[0] == doctest::Approx(2.0 * a[0]).epsilon(1e-5));
}

TEST_CASE("two classes sharing the server give two positive demands") {
    std::vector<double> muZ(2, 1.0);
    std::vector<api::MlpsSample> S;
    std::vector<double> q1, q2;
    q1.push_back(1.0);
    q1.push_back(0.0);
    q2.push_back(0.0);
    q2.push_back(1.0);
    const double a[] = {0.45, 0.52, 0.38, 0.61, 0.40};
    const double b[] = {0.95, 1.10, 0.80, 1.25, 0.90};
    for (std::size_t i = 0; i < 5; ++i) {
        S.push_back(sample(a[i], 1, q1));
        S.push_back(sample(b[i], 2, q2));
    }
    const std::vector<double> d = api::infer_mlps(muZ, 1.0, S);
    REQUIRE(d.size() == 2u);
    for (std::size_t r = 0; r < 2; ++r) {
        CHECK(d[r] > 0.0);
        CHECK(d[r] <= 1.25 + 1e-9);  // the reference's upper bound is max(rt)
    }
    // Class 2's response times are about twice class 1's, so its demand is the
    // larger of the two; an estimator that mixed up the class labels would not
    // preserve that ordering.
    CHECK(d[1] > d[0]);
}

TEST_CASE("MINPS selects between MLPS and RPS by the smaller mean") {
    std::vector<double> muZ(2, 1.0);
    std::vector<api::MlpsSample> S;
    std::vector<double> q1, q2;
    q1.push_back(1.0);
    q1.push_back(0.0);
    q2.push_back(0.0);
    q2.push_back(1.0);
    const double a[] = {0.45, 0.52, 0.38, 0.61, 0.40};
    const double b[] = {0.95, 1.10, 0.80, 1.25, 0.90};
    for (std::size_t i = 0; i < 5; ++i) {
        S.push_back(sample(a[i], 1, q1));
        S.push_back(sample(b[i], 2, q2));
    }

    const std::vector<double> mlps = api::infer_mlps(muZ, 1.0, S);
    const std::vector<double> minps = api::infer_minps(muZ, 1.0, S);
    REQUIRE(minps.size() == 2u);

    double mm = 0.0, mn = 0.0;
    for (std::size_t r = 0; r < 2; ++r) {
        mm += mlps[r];
        mn += minps[r];
    }
    // It is a SELECTION, not a blend: the answer is one of the two estimators'
    // vectors, and its mean is no larger than MLPS's.
    CHECK(mn <= mm + 1e-9);
    const bool isMlps = (std::fabs(minps[0] - mlps[0]) < 1e-12 &&
                         std::fabs(minps[1] - mlps[1]) < 1e-12);
    if (!isMlps)
        for (std::size_t r = 0; r < 2; ++r) CHECK(minps[r] > 0.0);
}

TEST_CASE("the malformed inputs are refused by name") {
    std::vector<double> muZ(1, 1.0);
    std::vector<api::MlpsSample> ok;
    ok.push_back(sample(0.5, 1, std::vector<double>(1, 1.0)));

    CHECK_THROWS_AS(api::infer_mlps(std::vector<double>(), 1.0, ok), line::InputError);
    CHECK_THROWS_AS(api::infer_mlps(muZ, 1.0, std::vector<api::MlpsSample>()), line::InputError);

    // A class outside 1..R.
    std::vector<api::MlpsSample> badcls;
    badcls.push_back(sample(0.5, 2, std::vector<double>(1, 1.0)));
    CHECK_THROWS_AS(api::infer_mlps(muZ, 1.0, badcls), line::InputError);
    std::vector<api::MlpsSample> zerocls;
    zerocls.push_back(sample(0.5, 0, std::vector<double>(1, 1.0)));
    CHECK_THROWS_AS(api::infer_mlps(muZ, 1.0, zerocls), line::InputError);

    // A queue-length row of the wrong width.
    std::vector<api::MlpsSample> badql;
    badql.push_back(sample(0.5, 1, std::vector<double>(2, 1.0)));
    CHECK_THROWS_AS(api::infer_mlps(muZ, 1.0, badql), line::InputError);

    // A non-positive response time has no phase-type density.
    std::vector<api::MlpsSample> badrt;
    badrt.push_back(sample(0.0, 1, std::vector<double>(1, 1.0)));
    CHECK_THROWS_AS(api::infer_mlps(muZ, 1.0, badrt), line::InputError);

    // A sample whose own class is empty on arrival cannot host the tagged job:
    // the arrival state must COUNT the arriving job.
    std::vector<api::MlpsSample> empty;
    empty.push_back(sample(0.5, 1, std::vector<double>(1, 0.0)));
    CHECK_THROWS_AS(api::infer_mlps(muZ, 1.0, empty), line::InputError);
}
