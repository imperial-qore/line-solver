/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * infer_minps_setup: turning a raw per-class trace into a MINPS sample set.
 *
 * THIS IS ALL PREPARATION, AND EVERY STEP OF IT CHANGES WHAT THE ESTIMATOR
 * SEES, so each is driven separately here rather than checked through the
 * estimate at the end:
 *
 *  - an empty class is DROPPED and the survivors are RELABELLED, so the caller
 *    is handed the map back; a test that only looked at the sample count would
 *    not notice the labels shifting under it;
 *  - the merged stream is sorted by ARRIVAL TIME, not grouped by class,
 *    because the estimator conditions on the state a job found;
 *  - the window is CONTIGUOUS, so a request that runs off the end is refused
 *    rather than shortened;
 *  - a non-positive response time is dropped, since it has no density.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/infer/infer_minps_setup.h"

namespace api = line::api;

namespace {

api::MinpsClassTrace trace(const std::vector<double>& at_ms, const std::vector<double>& rt) {
    api::MinpsClassTrace t;
    t.arrival_ms = at_ms;
    t.rt = rt;
    return t;
}

}  // namespace

TEST_CASE("the merged stream is sorted by arrival time, not grouped by class") {
    std::vector<api::MinpsClassTrace> tr;
    // Class 1 arrives at 0 and 200 ms, class 2 at 100 and 300: interleaved.
    tr.push_back(trace({0.0, 200.0}, {0.05, 0.06}));
    tr.push_back(trace({100.0, 300.0}, {0.07, 0.08}));

    const api::MinpsSetup st = api::infer_minps_setup(tr, 1, 0);
    REQUIRE(st.samples.size() == 4u);
    CHECK(st.samples[0].cls == 1u);
    CHECK(st.samples[1].cls == 2u);
    CHECK(st.samples[2].cls == 1u);
    CHECK(st.samples[3].cls == 2u);
    // Two surviving classes, in their original order.
    REQUIRE(st.classMap.size() == 2u);
    CHECK(st.classMap[0] == 1u);
    CHECK(st.classMap[1] == 2u);
    for (std::size_t i = 0; i < st.samples.size(); ++i) CHECK(st.samples[i].ql.size() == 2u);
}

TEST_CASE("an empty class is dropped and the survivors are relabelled") {
    std::vector<api::MinpsClassTrace> tr;
    tr.push_back(trace({0.0}, {0.05}));              // class 1
    tr.push_back(trace({}, {}));                      // class 2, never seen
    tr.push_back(trace({100.0, 200.0}, {0.07, 0.08}));  // class 3

    const api::MinpsSetup st = api::infer_minps_setup(tr, 1, 0);
    REQUIRE(st.classMap.size() == 2u);
    // The MAP is what says class 2 of the sample set is class 3 of the trace.
    CHECK(st.classMap[0] == 1u);
    CHECK(st.classMap[1] == 3u);
    REQUIRE(st.samples.size() == 3u);
    // The labels are the NEW ones, dense from 1.
    CHECK(st.samples[0].cls == 1u);
    CHECK(st.samples[1].cls == 2u);
    CHECK(st.samples[2].cls == 2u);
    // The queue-length rows are as wide as the SURVIVING class count.
    for (std::size_t i = 0; i < st.samples.size(); ++i) CHECK(st.samples[i].ql.size() == 2u);
    CHECK(st.lambda.size() == 2u);
}

TEST_CASE("the window is contiguous and is applied to the sorted stream") {
    std::vector<api::MinpsClassTrace> tr;
    tr.push_back(trace({0.0, 200.0, 400.0}, {0.05, 0.06, 0.07}));
    tr.push_back(trace({100.0, 300.0, 500.0}, {0.11, 0.12, 0.13}));

    const api::MinpsSetup all = api::infer_minps_setup(tr, 1, 0);
    CHECK(all.samples.size() == 6u);

    const api::MinpsSetup mid = api::infer_minps_setup(tr, 3, 2);
    REQUIRE(mid.samples.size() == 2u);
    // Samples 3 and 4 of the sorted stream are class 1 at 200 ms and class 2
    // at 300 ms.
    CHECK(mid.samples[0].cls == 1u);
    CHECK(mid.samples[0].rt == doctest::Approx(0.06));
    CHECK(mid.samples[1].cls == 2u);
    CHECK(mid.samples[1].rt == doctest::Approx(0.12));

    // Running past the end is refused rather than shortened: a truncated window
    // would silently answer a different question.
    CHECK_THROWS_AS(api::infer_minps_setup(tr, 5, 4), line::InputError);
    CHECK_THROWS_AS(api::infer_minps_setup(tr, 7, 1), line::InputError);
}

TEST_CASE("a non-positive response time is dropped from the sample set") {
    std::vector<api::MinpsClassTrace> tr;
    tr.push_back(trace({0.0, 100.0, 200.0}, {0.05, 0.0, 0.07}));

    const api::MinpsSetup st = api::infer_minps_setup(tr, 1, 0);
    // Three in the window, two usable.
    REQUIRE(st.samples.size() == 2u);
    CHECK(st.samples[0].rt == doctest::Approx(0.05));
    CHECK(st.samples[1].rt == doctest::Approx(0.07));
}

TEST_CASE("the think-time rates are finite, non-negative and capped") {
    std::vector<api::MinpsClassTrace> tr;
    tr.push_back(trace({0.0, 200.0, 400.0}, {0.05, 0.06, 0.07}));
    tr.push_back(trace({100.0, 300.0, 500.0}, {0.11, 0.12, 0.13}));

    const api::MinpsSetup st = api::infer_minps_setup(tr, 1, 0);
    REQUIRE(st.lambda.size() == 2u);
    for (std::size_t k = 0; k < 2; ++k) {
        CHECK(std::isfinite(st.lambda[k]));
        CHECK(st.lambda[k] >= 0.0);
        CHECK(st.lambda[k] <= 1e6);
    }
    // Wexp is the largest total queue length any job found, so it is at least
    // one: every job sees itself.
    CHECK(st.threads >= 1.0);
}

TEST_CASE("the malformed traces are refused by name") {
    std::vector<api::MinpsClassTrace> none;
    CHECK_THROWS_AS(api::infer_minps_setup(none, 1, 0), line::InputError);

    std::vector<api::MinpsClassTrace> empty;
    empty.push_back(trace({}, {}));
    CHECK_THROWS_AS(api::infer_minps_setup(empty, 1, 0), line::InputError);

    std::vector<api::MinpsClassTrace> ragged;
    ragged.push_back(trace({0.0, 1.0}, {0.05}));
    CHECK_THROWS_AS(api::infer_minps_setup(ragged, 1, 0), line::InputError);

    // The window index is 1-based, as the reference's is.
    std::vector<api::MinpsClassTrace> ok;
    ok.push_back(trace({0.0, 100.0}, {0.05, 0.06}));
    CHECK_THROWS_AS(api::infer_minps_setup(ok, 0, 1), line::InputError);
}

TEST_CASE("the prepared trace drives MINPS end to end") {
    std::vector<api::MinpsClassTrace> tr;
    // Every job arrives long after the previous one finished, so each finds an
    // empty queue and the MINPS estimate is well posed.
    tr.push_back(trace({0.0, 1000.0, 2000.0, 3000.0}, {0.30, 0.55, 0.42, 0.71}));

    const std::vector<double> d = api::infer_minps_from_trace(tr, 1, 0, 1.0);
    REQUIRE(d.size() == 1u);
    CHECK(d[0] > 0.0);
    // The demands cannot exceed the largest response time observed.
    CHECK(d[0] <= 0.71 + 1e-9);
}
