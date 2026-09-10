/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * SolverLDES, the C++ client of the SSJ discrete-event engine.
 *
 * WHAT IS ASSERTED, in three kinds:
 *
 *  1. THE COMMAND LINE. The client's only real decisions are which flags to
 *     emit and which image to run, and both are silent failures when wrong: a
 *     flag omitted because its value equals a default THIS port invented would
 *     make the engine run something else, and a cost-capped cache sent to the
 *     prebuilt AOT image would simulate the UNCAPPED cache and report a number
 *     that looks fine. So the flag vector is checked value by value against
 *     `_build_flag_args` and `solveCli.m`, and the runner order against the
 *     rule those two share.
 *
 *  2. THE RESULT DOCUMENT. `parse_ldes_result` is checked against a document
 *     written by hand in the shape `jline/io/LDESResultIO.java` emits, including
 *     the parts the MATLAB and Python clients do not read (sampleCounts,
 *     impatience, totalSimulatedEvents), the `null`-is-NaN convention, and the
 *     two documents that must NOT be read as results: an `error` object and a
 *     foreign document with no `ldes-result` format tag.
 *
 *  3. THE ENGINE ITSELF, when this machine has one. Those cases are guarded by
 *     `ldes_is_available()` and skipped otherwise -- LINE ships no engine in the
 *     source tree, so a hard failure here would be a failure of the checkout and
 *     not of the port. What they assert is what only a real run can: that a
 *     closed model conserves its population, that Little's law holds on the
 *     returned table, and that one seed gives one answer twice.
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/wrappers/ldes/solver_ldes.h"

using namespace line;

namespace {

/** The model.json of a closed two-station network, N = 2. */
const char* kCqnDoc = R"({
  "format": "line-model", "version": "1.0",
  "model": {"type": "Network", "name": "ldes_cqn",
    "nodes": [
      {"name": "Delay", "type": "Delay",
       "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 1.0}}}},
      {"name": "Queue", "type": "Queue", "scheduling": "FCFS",
       "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 0.5}}}}],
    "classes": [{"name": "Class1", "type": "Closed", "population": 2, "refNode": "Delay"}],
    "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
      "Delay": {"Queue": 1.0}, "Queue": {"Delay": 1.0}}}}}})";

/** True when `flags` carries `name` followed by `value`. */
bool has_flag(const std::vector<std::string>& flags, const std::string& name,
              const std::string& value) {
    for (std::size_t i = 0; i + 1 < flags.size(); ++i)
        if (flags[i] == name && flags[i + 1] == value) return true;
    return false;
}

bool has_switch(const std::vector<std::string>& flags, const std::string& name) {
    for (std::size_t i = 0; i < flags.size(); ++i)
        if (flags[i] == name) return true;
    return false;
}

}  // namespace

TEST_CASE("ldes_flags emits the budget and the seed always, and nothing else by default") {
    const ldes::LdesOptions o;
    const std::vector<std::string> f = ldes::ldes_flags(o, std::vector<std::string>());
    // EXACTLY four tokens. Every other knob is at the engine's own default and a
    // default run has to produce the minimal command line an older AOT native
    // image still parses.
    REQUIRE(f.size() == 4);
    CHECK(f[0] == "-s");
    CHECK(f[1] == "200000");
    CHECK(f[2] == "--seed");
    CHECK(f[3] == "23000");
}

TEST_CASE("the seed is emitted even at its nominal value") {
    // The CLI's default is -1, a RANDOM stream, not 23000. Omitting --seed here
    // because "it is the default" would make the native binary pick a random
    // seed and the run irreproducible, and would silently disagree with the
    // MATLAB and Python clients, which both emit it unconditionally.
    ldes::LdesOptions o;
    o.seed = 23000;
    CHECK(has_flag(ldes::ldes_flags(o, std::vector<std::string>()), "--seed", "23000"));
    o.seed = -1;
    CHECK(has_flag(ldes::ldes_flags(o, std::vector<std::string>()), "--seed", "-1"));
}

TEST_CASE("options.events overrides options.samples as the event budget") {
    ldes::LdesOptions o;
    o.samples = 5000;
    o.events = 12345;
    CHECK(has_flag(ldes::ldes_flags(o, std::vector<std::string>()), "-s", "12345"));
}

TEST_CASE("every non-default knob reaches the command line under its engine name") {
    ldes::LdesOptions o;
    o.samples = 1000;
    o.method = "default";
    o.tranfilter = "fixed";
    o.warmupfrac = 0.0;
    o.mserbatch = 7;
    o.cimethod = "spectral";
    o.obmoverlap = 0.25;
    o.ciminbatch = 20;
    o.ciminobs = 50;
    o.spectral_low_freq_frac = 0.5;
    o.cnvgon = true;
    o.cnvgtol = 0.01;
    o.cnvgbatch = 30;
    o.cnvgchk = 500;
    o.slotted = true;
    o.slot_length = 0.5;
    o.replications = 4;
    o.numthreads = 2;
    o.has_timespan = true;
    o.t0 = 0.0;
    o.t1 = 100.0;
    o.timeout = 60.0;
    o.init_sol.push_back(2.0);
    o.init_sol.push_back(0.0);
    const std::vector<std::string> f = ldes::ldes_flags(o, std::vector<std::string>());
    CHECK(has_flag(f, "-s", "1000"));
    CHECK(has_flag(f, "--tranfilter", "fixed"));
    CHECK(has_flag(f, "--warmupfrac", "0"));
    CHECK(has_flag(f, "--mserbatch", "7"));
    CHECK(has_flag(f, "--cimethod", "spectral"));
    CHECK(has_flag(f, "--obmoverlap", "0.25"));
    CHECK(has_flag(f, "--ciminbatch", "20"));
    CHECK(has_flag(f, "--ciminobs", "50"));
    CHECK(has_flag(f, "--spectrallowfreqfrac", "0.5"));
    CHECK(has_switch(f, "--cnvgon"));
    CHECK(has_flag(f, "--cnvgtol", "0.01"));
    CHECK(has_flag(f, "--cnvgbatch", "30"));
    CHECK(has_flag(f, "--cnvgchk", "500"));
    CHECK(has_switch(f, "--slotted"));
    CHECK(has_flag(f, "--slotlength", "0.5"));
    CHECK(has_flag(f, "--replications", "4"));
    CHECK(has_flag(f, "--numthreads", "2"));
    CHECK(has_flag(f, "--timespan", "0,100"));
    CHECK(has_flag(f, "--maxtime", "60"));
    CHECK(has_flag(f, "--initsol", "2,0"));
    // `--method default` is NOT emitted: it is the engine's own, and an older
    // image that does not know the flag would refuse a run that asked for
    // nothing.
    CHECK(!has_switch(f, "--method"));
}

TEST_CASE("the convergence knobs stay silent while convergence stopping is off") {
    // They are gated on --cnvgon in both other clients: a tolerance the engine
    // never consults is a command line that misstates what the run did.
    ldes::LdesOptions o;
    o.cnvgtol = 0.01;
    o.cnvgbatch = 30;
    const std::vector<std::string> f = ldes::ldes_flags(o, std::vector<std::string>());
    CHECK(!has_switch(f, "--cnvgtol"));
    CHECK(!has_switch(f, "--cnvgbatch"));
}

TEST_CASE("an infinite timeout emits no budget at all") {
    const ldes::LdesOptions o;
    CHECK(!has_switch(ldes::ldes_flags(o, std::vector<std::string>()), "--maxtime"));
}

TEST_CASE("extra flags are appended verbatim, after the resolved knobs") {
    ldes::LdesOptions o;
    std::vector<std::string> extra;
    extra.push_back("--export-histogram");
    extra.push_back("--trajectory");
    const std::vector<std::string> f = ldes::ldes_flags(o, extra);
    REQUIRE(f.size() == 6);
    CHECK(f[4] == "--export-histogram");
    CHECK(f[5] == "--trajectory");
}

TEST_CASE("the flag formatter prints the shortest decimal that reads back") {
    CHECK(ldes::detail::shortest(0.05) == "0.05");
    CHECK(ldes::detail::shortest(0.2) == "0.2");
    CHECK(ldes::detail::shortest(1.0) == "1");
    CHECK(ldes::detail::shortest(1e-7) == "1e-07");
    // The point of the round trip: a fixed %.10g would round this to something
    // the engine would then parse as a DIFFERENT tolerance.
    const double x = 0.1234567890123456789;
    CHECK(std::strtod(ldes::detail::shortest(x).c_str(), nullptr) == x);
}

TEST_CASE("a document that postdates the prebuilt image is detected on the wire") {
    // The two other clients walk the live model object for these; the client
    // here has only the bytes, which is the same information at the point where
    // it crosses.
    CHECK(ldes::detail::document_postdates_native("{\"costCaps\": 3}"));
    CHECK(ldes::detail::document_postdates_native("{\"markedClasses\": [\"A\"]}"));
    CHECK(!ldes::detail::document_postdates_native(kCqnDoc));
}

TEST_CASE("parse_ldes_result reads every block of the ldes-result document") {
    const std::string doc = R"({
      "format": "ldes-result", "version": "1.0",
      "solver": "SolverLDES", "method": "default", "runtime": 0.25,
      "converged": true, "stoppingReason": "convergence", "convergenceBatches": 7,
      "totalSimulatedEvents": 20000,
      "dimensions": {"nstations": 2, "nclasses": 1, "nchains": 1,
                     "stationNames": ["Delay", "Queue"], "classNames": ["Class1"]},
      "metrics": {"QN": [[1.2], [0.8]], "UN": [[1.19], [0.6]], "RN": [[1.0], [0.66]],
                  "TN": [[1.19], [1.19]], "AN": [[1.19], [1.19]], "WN": [[1.0], [0.66]],
                  "CN": [[1.0]], "XN": [[1.19]], "DropRateJoin": null},
      "sampleCounts": {"QNSamples": [[995.0], [995.0]], "UNSamples": null,
                       "RNSamples": null, "TNSamples": null},
      "confidenceIntervals": {"QNCI": [[0.011], [0.011]], "UNCI": null, "RNCI": null,
                              "TNCI": null, "ANCI": null, "WNCI": null},
      "relativePrecision": {"QNRelPrec": [[0.5]], "UNRelPrec": null, "RNRelPrec": null,
                            "TNRelPrec": null},
      "fcr": {"nregions": 1, "QNfcr": [[0.4]], "UNfcr": [[null]], "RNfcr": [[0.9]],
              "TNfcr": [[0.44]], "ANfcr": [[null]], "WNfcr": [[0.9]],
              "WeightNfcr": [[0.8]], "MemOccNfcr": [[1.6]], "DropRateNfcr": [[0.0]]},
      "impatience": {"renegedCustomers": [[3.0], [0.0]], "avgRenegingWaitTime": null,
                     "renegingRate": [[0.1], [0.0]], "balkedCustomers": null,
                     "balkingProbability": null, "retriedCustomers": null,
                     "retrialDropped": null, "avgOrbitSize": null},
      "cacheMetrics": {"MyCache": {"hit": [[0.3, 0.7]], "delayed": [[0.0, 0.0]],
                                   "miss": [[0.7, 0.3]], "latency": [[1.5, 2.5]],
                                   "hitList": [[0.3], [0.7]], "itemProb": [[0.5, 0.5]],
                                   "listCost": [[2.0]]}},
      "stateHistogram": {"space": [[0, 2], [1, 1], [2, 0]], "time": [[10.0], [20.0], [70.0]],
                         "trajSpace": [[0, 2], [2, 0]], "trajTime": [[0.0], [1.0]]},
      "respTimeSamples": [[[0.5, 1.5]], [[0.25]]],
      "transient": {"t": [0.0, 1.0, 2.0],
                    "QNt": [[[[0.0, 0.0], [1.0, 1.0], [1.2, 2.0]]], [[[2.0, 0.0], [1.0, 1.0], [0.8, 2.0]]]],
                    "UNt": [[[[0.0, 0.0], [0.9, 1.0], [1.0, 2.0]]], [[[1.0, 0.0], [0.7, 1.0], [0.6, 2.0]]]],
                    "TNt": [[[[0.0, 0.0], [1.1, 1.0], [1.19, 2.0]]], [[[0.0, 0.0], [1.1, 1.0], [1.19, 2.0]]]]}
    })";
    const ldes::LdesResult r = ldes::parse_ldes_result(ldes::detail::Json::parse(doc));

    CHECK(r.method == "default");
    CHECK(r.runtime == doctest::Approx(0.25));
    CHECK(r.converged);
    CHECK(r.stopping_reason == "convergence");
    CHECK(r.convergence_batches == 7);
    CHECK(r.total_simulated_events == 20000);

    REQUIRE(r.nstations == 2);
    REQUIRE(r.nclasses == 1);
    CHECK(r.nchains == 1);
    REQUIRE(r.station_names.size() == 2);
    CHECK(r.station_names[0] == "Delay");
    CHECK(r.station_names[1] == "Queue");
    REQUIRE(r.class_names.size() == 1);
    CHECK(r.class_names[0] == "Class1");

    REQUIRE(r.QN.rows() == 2);
    CHECK(r.QN(0, 0) == doctest::Approx(1.2));
    CHECK(r.TN(1, 0) == doctest::Approx(1.19));
    CHECK(r.CN(0, 0) == doctest::Approx(1.0));
    // A `null` metric is ABSENT, not zero: the two mean different things and the
    // tables downstream distinguish them.
    CHECK(r.DropRateJoin.empty());
    CHECK(r.UNCI.empty());
    CHECK(r.QNCI(0, 0) == doctest::Approx(0.011));
    CHECK(r.QNSamples(0, 0) == doctest::Approx(995.0));

    CHECK(r.nregions == 1);
    CHECK(r.QNfcr(0, 0) == doctest::Approx(0.4));
    CHECK(r.WeightNfcr(0, 0) == doctest::Approx(0.8));
    CHECK(r.MemOccNfcr(0, 0) == doctest::Approx(1.6));
    // A `null` INSIDE a matrix is NaN: the region has no arrival rate, and a
    // zero there would read as "measured, and zero".
    REQUIRE(!r.UNfcr.empty());
    CHECK(std::isnan(r.UNfcr(0, 0)));

    CHECK(r.renegedCustomers(0, 0) == doctest::Approx(3.0));
    CHECK(r.renegingRate(0, 0) == doctest::Approx(0.1));
    CHECK(r.avgRenegingWaitTime.empty());

    REQUIRE(r.cache_metrics.count("MyCache") == 1);
    const ldes::LdesCacheMetrics& cm = r.cache_metrics.find("MyCache")->second;
    CHECK(cm.hit(0, 1) == doctest::Approx(0.7));
    CHECK(cm.miss(0, 0) == doctest::Approx(0.7));
    CHECK(cm.latency(0, 1) == doctest::Approx(2.5));
    CHECK(cm.hitList.rows() == 2);
    CHECK(cm.listCost(0, 0) == doctest::Approx(2.0));

    REQUIRE(r.histogram_space.rows() == 3);
    CHECK(r.histogram_space.cols() == 2);
    CHECK(r.histogram_time(2, 0) == doctest::Approx(70.0));
    CHECK(r.traj_space.rows() == 2);
    CHECK(r.traj_time(1, 0) == doctest::Approx(1.0));

    REQUIRE(r.t.size() == 3);
    CHECK(r.t[2] == doctest::Approx(2.0));
    REQUIRE(r.QNt.size() == 2);
    REQUIRE(r.QNt[0].size() == 1);
    CHECK(r.QNt[0][0].rows() == 3);
    // Column 0 is the VALUE and column 1 the time, the layout LDESResultIO
    // writes; reading them the other way round would silently plot time.
    CHECK(r.QNt[0][0](2, 0) == doctest::Approx(1.2));
    CHECK(r.QNt[0][0](2, 1) == doctest::Approx(2.0));
    CHECK(r.UNt[1][0](0, 0) == doctest::Approx(1.0));
    CHECK(r.TNt[1][0](2, 0) == doctest::Approx(1.19));

    REQUIRE(r.respTimeSamples.size() == 2);
    REQUIRE(r.respTimeSamples[0][0].size() == 2);
    CHECK(r.respTimeSamples[0][0][1] == doctest::Approx(1.5));
    CHECK(r.respTimeSamples[1][0][0] == doctest::Approx(0.25));
}

TEST_CASE("an engine error document is an error, not an empty result") {
    const std::string doc = R"({"error": "LDES currently supports Network models only"})";
    CHECK_THROWS_AS(ldes::parse_ldes_result(ldes::detail::Json::parse(doc)), line::NumericError);
}

TEST_CASE("a foreign document is refused rather than read as a result of NaNs") {
    // The runner falls through on failure; a runner that wrote some other JSON
    // to the output path would otherwise be read as a result with every metric
    // absent, which is indistinguishable from a model that measured nothing.
    const std::string doc = R"({"format": "line-model", "version": "1.0"})";
    CHECK_THROWS_AS(ldes::parse_ldes_result(ldes::detail::Json::parse(doc)), line::NumericError);
}

TEST_CASE("the runner list puts the native image first and the jar behind it") {
    if (!ldes::ldes_is_available()) return;  // no engine in this checkout
    const std::vector<ldes::LdesRunner> plain =
        ldes::ldes_runners(kCqnDoc, ldes::ldes_flags(ldes::LdesOptions(), std::vector<std::string>()));
    REQUIRE(!plain.empty());
    if (plain.size() < 2) return;  // only one image present, nothing to order
    CHECK(plain[0].engine == "native");
    CHECK(plain[1].engine == "jar");

    // A flag the prebuilt image predates FLIPS the order rather than extending
    // it: on that image --respt-samples is accepted and ignored, so the arm
    // would refuse for want of samples nobody asked for.
    std::vector<std::string> extra;
    extra.push_back("--respt-samples");
    const std::vector<ldes::LdesRunner> flipped =
        ldes::ldes_runners(kCqnDoc, ldes::ldes_flags(ldes::LdesOptions(), extra));
    REQUIRE(flipped.size() == 2);
    CHECK(flipped[0].engine == "jar");

    // A cost-capped cache does the same, and for the worse reason: the image
    // would simulate the UNCAPPED cache and report a number that looks fine.
    const std::vector<ldes::LdesRunner> capped = ldes::ldes_runners(
        "{\"costCaps\": [3]}", ldes::ldes_flags(ldes::LdesOptions(), std::vector<std::string>()));
    REQUIRE(capped.size() == 2);
    CHECK(capped[0].engine == "jar");
}

TEST_CASE("a closed model solved by the engine conserves its population") {
    if (!ldes::ldes_is_available()) return;
    ldes::LdesOptions o;
    o.samples = 20000;
    o.seed = 23000;
    const ldes::LdesResult r = ldes::solver_ldes_text(kCqnDoc, o, std::vector<std::string>());

    REQUIRE(r.nstations == 2);
    REQUIRE(r.nclasses == 1);
    CHECK(r.station_names[0] == "Delay");
    CHECK((r.engine == "native" || r.engine == "jar"));
    CHECK(r.total_simulated_events > 0);

    // The population is an INVARIANT of the sample path, not an estimate: every
    // job is at one of the two stations at every instant, so the two mean queue
    // lengths sum to N however short the run was.
    const double pop = r.QN(0, 0) + r.QN(1, 0);
    CHECK(pop == doctest::Approx(2.0).epsilon(1e-9));
    // Little's law on the returned table, station by station. It holds
    // ASYMPTOTICALLY and not exactly, which is a property of the estimators and
    // not a defect: Q is a time average over the post-warmup window, R the mean
    // of the per-job observations and T a completion count, and the three are
    // formed over different observation counts (the `sampleCounts` block reports
    // 995 for Q against 10000 for R on this model). So the residual is checked
    // at 1% here and then checked to FALL with the run length -- 1.9e-3 relative
    // at 2e4 events, 2.0e-4 at 2e5, 5.7e-5 at 2e6 -- which is what says it is
    // Monte Carlo noise converging to zero rather than a mis-indexed matrix,
    // which would not move at all.
    for (std::size_t i = 0; i < 2; ++i)
        CHECK(r.QN(i, 0) == doctest::Approx(r.TN(i, 0) * r.RN(i, 0)).epsilon(0.01));
    ldes::LdesOptions longer = o;
    longer.samples = 200000;
    const ldes::LdesResult rl =
        ldes::solver_ldes_text(kCqnDoc, longer, std::vector<std::string>());
    for (std::size_t i = 0; i < 2; ++i) {
        const double short_gap = std::fabs(r.QN(i, 0) - r.TN(i, 0) * r.RN(i, 0));
        const double long_gap = std::fabs(rl.QN(i, 0) - rl.TN(i, 0) * rl.RN(i, 0));
        CHECK(long_gap < short_gap);
    }
    // The Delay is an infinite server with mean 1, so its response time is that
    // mean whatever the load.
    CHECK(r.RN(0, 0) == doctest::Approx(1.0).epsilon(0.05));
    CHECK(r.UN(1, 0) > 0.0);
    CHECK(r.UN(1, 0) < 1.0);
}

TEST_CASE("one seed gives one answer twice, and a different seed does not") {
    if (!ldes::ldes_is_available()) return;
    ldes::LdesOptions o;
    o.samples = 5000;
    o.seed = 4242;
    const ldes::LdesResult a = ldes::solver_ldes_text(kCqnDoc, o, std::vector<std::string>());
    const ldes::LdesResult b = ldes::solver_ldes_text(kCqnDoc, o, std::vector<std::string>());
    CHECK(a.QN(1, 0) == doctest::Approx(b.QN(1, 0)).epsilon(1e-12));
    CHECK(a.TN(1, 0) == doctest::Approx(b.TN(1, 0)).epsilon(1e-12));

    // Reproducibility is only meaningful if the seed is what produced it: a
    // client that dropped --seed would pass the check above by accident on a
    // deterministic default and fail to be a simulation at all.
    o.seed = 777;
    const ldes::LdesResult c = ldes::solver_ldes_text(kCqnDoc, o, std::vector<std::string>());
    CHECK(a.QN(1, 0) != doctest::Approx(c.QN(1, 0)).epsilon(1e-12));
}

TEST_CASE("the transient run returns the trajectory and the steady one does not") {
    if (!ldes::ldes_is_available()) return;
    ldes::LdesOptions o;
    o.samples = 5000;
    o.seed = 23000;
    const ldes::LdesResult steady = ldes::solver_ldes_text(kCqnDoc, o, std::vector<std::string>());
    CHECK(steady.t.empty());
    CHECK(steady.QNt.empty());

    o.has_timespan = true;
    o.t0 = 0.0;
    o.t1 = 50.0;
    std::vector<std::string> extra;
    extra.push_back("--trajectory");
    const ldes::LdesResult tran = ldes::solver_ldes_text(kCqnDoc, o, extra);
    REQUIRE(!tran.t.empty());
    REQUIRE(tran.QNt.size() == 2);
    // The series are STATION-major, as `result.QNt = new Matrix[numStations]
    // [numClasses]` writes them.
    REQUIRE(!tran.QNt[0].empty());
    CHECK(tran.QNt[0][0].cols() == 2);
    // The trajectory stays inside the horizon it was asked for.
    CHECK(tran.t.back() <= 50.0 + 1e-9);
}

TEST_CASE("--respt-samples returns the per-job response times an ecdf is built from") {
    if (!ldes::ldes_is_available()) return;
    ldes::LdesOptions o;
    o.samples = 5000;
    o.seed = 23000;
    std::vector<std::string> extra;
    extra.push_back("--respt-samples");
    const ldes::LdesResult r = ldes::solver_ldes_text(kCqnDoc, o, extra);
    REQUIRE(!r.respTimeSamples.empty());
    std::size_t observed = 0;
    for (std::size_t i = 0; i < r.respTimeSamples.size(); ++i)
        for (std::size_t c = 0; c < r.respTimeSamples[i].size(); ++c)
            observed += r.respTimeSamples[i][c].size();
    CHECK(observed > 0);
}

TEST_CASE("--export-histogram returns a residence-time law over the joint state") {
    if (!ldes::ldes_is_available()) return;
    ldes::LdesOptions o;
    o.samples = 5000;
    o.seed = 23000;
    std::vector<std::string> extra;
    extra.push_back("--export-histogram");
    const ldes::LdesResult r = ldes::solver_ldes_text(kCqnDoc, o, extra);
    REQUIRE(r.histogram_space.rows() > 0);
    REQUIRE(r.histogram_time.rows() == r.histogram_space.rows());
    // Two stations and one class: the aggregate row is (station, class) wide,
    // and each visited state holds the whole population.
    CHECK(r.histogram_space.cols() == 2);
    double total = 0.0;
    for (std::size_t s = 0; s < r.histogram_time.rows(); ++s) {
        total += r.histogram_time(s, 0);
        const double held = r.histogram_space(s, 0) + r.histogram_space(s, 1);
        CHECK(held == doctest::Approx(2.0));
    }
    CHECK(total > 0.0);

    // E[n] over the histogram is the queue length the means report: the two are
    // the same measurement, so a histogram indexed the wrong way round would
    // show up here and nowhere else.
    double en = 0.0;
    for (std::size_t s = 0; s < r.histogram_time.rows(); ++s)
        en += (r.histogram_time(s, 0) / total) * r.histogram_space(s, 1);
    CHECK(en == doctest::Approx(r.QN(1, 0)).epsilon(0.02));
}

TEST_CASE("getProbAggr weighs residence times and recovers the product-form marginal") {
    if (!ldes::ldes_is_available()) return;
    // kCqnDoc is a machine repairman: a Delay of mean 1 and an FCFS Queue of
    // mean 0.5, two jobs of one class. Its queue-length law is product form,
    //   P(n) = D^n Z^(N-n)/(N-n)! / G,  D = 0.5, Z = 1, G = 1.25,
    // so P(0) = 0.4, P(1) = 0.4, P(2) = 0.2 EXACTLY. This is the model on which
    // the trajectory-based estimator of the three other clients reported 0.0101
    // for the state whose true probability is 0.4 (BUG-96): the transient QNt
    // series holds interval MEANS of the queue length, so it lands on the
    // integer state 1 only by coincidence.
    ldes::LdesOptions o;
    o.samples = 200000;
    o.seed = 23000;
    const std::size_t queue = 2;  // Delay is station 1, Queue station 2
    const double p0 = ldes::ldes_prob_aggr(kCqnDoc, o, queue, std::vector<double>(1, 0.0), 1);
    const double p1 = ldes::ldes_prob_aggr(kCqnDoc, o, queue, std::vector<double>(1, 1.0), 1);
    const double p2 = ldes::ldes_prob_aggr(kCqnDoc, o, queue, std::vector<double>(1, 2.0), 1);

    CHECK(p0 == doctest::Approx(0.4).epsilon(0.05));
    CHECK(p1 == doctest::Approx(0.4).epsilon(0.05));
    CHECK(p2 == doctest::Approx(0.2).epsilon(0.05));
    // The three states are the whole reachable space, so they exhaust the mass.
    CHECK(p0 + p1 + p2 == doctest::Approx(1.0).epsilon(1e-6));
    // A state outside the space is never visited, and that is a zero rather
    // than an error.
    CHECK(ldes::ldes_prob_aggr(kCqnDoc, o, queue, std::vector<double>(1, 3.0), 1) == 0.0);
}

TEST_CASE("getProbSysAggr constrains every station at once") {
    if (!ldes::ldes_is_available()) return;
    ldes::LdesOptions o;
    o.samples = 200000;
    o.seed = 23000;
    // The population is conserved, so (1 at the Delay, 1 at the Queue) is the
    // joint state whose marginal at the Queue is P(1) = 0.4; constraining both
    // stations cannot select less than constraining one of them.
    Matrix<double> target(2, 1, 0.0);
    target(0, 0) = 1.0;
    target(1, 0) = 1.0;
    CHECK(ldes::ldes_prob_sys_aggr(kCqnDoc, o, target) == doctest::Approx(0.4).epsilon(0.05));
    // A joint state that violates the conservation law is unreachable.
    target(0, 0) = 2.0;
    target(1, 0) = 2.0;
    CHECK(ldes::ldes_prob_sys_aggr(kCqnDoc, o, target) == 0.0);
}

TEST_CASE("a state probability without a histogram is refused, not answered as zero") {
    ldes::LdesResult empty;
    std::vector<ldes::LdesStateQuery> q(1);
    q[0].station = 1;
    q[0].counts = std::vector<double>(1, 0.0);
    // An absent histogram means the run did not carry --export-histogram. Zero
    // is a legitimate probability, so it must not be how that case reads.
    CHECK_THROWS_AS(ldes::ldes_prob_from_histogram(empty, 1, q), line::Error);
}

TEST_CASE("a model.json the engine refuses is reported with the engine's own words") {
    if (!ldes::ldes_is_available()) return;
    const ldes::LdesOptions o;
    CHECK_THROWS_AS(ldes::solver_ldes_text("{\"format\": \"line-model\", \"model\": {}}", o,
                                           std::vector<std::string>()),
                    line::Error);
}

TEST_CASE("a model built through the C++ API reaches the engine through the writer") {
    if (!ldes::ldes_is_available()) return;
    // The document entry is the primary one -- a caller holding a model.json
    // should pass the BYTES -- but a caller who built the model in C++ has only
    // a struct, and `solver_ldes` serializes it with the writer whose output the
    // reference readers consume. The two must agree on a model the writer can
    // express in full, which is what this checks.
    qn::Network<double> m("ldes_cqn");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("Class1", 2, d);
    m.set_service(d, c, lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    ldes::LdesOptions o;
    o.samples = 20000;
    o.seed = 23000;
    const ldes::LdesResult viaStruct = ldes::solver_ldes(m.get_struct(), o);
    const ldes::LdesResult viaDoc = ldes::solver_ldes_text(kCqnDoc, o, std::vector<std::string>());
    REQUIRE(viaStruct.nstations == 2);
    // Same model, same seed, same engine: the same numbers, not merely close
    // ones. A writer that renamed a station or reordered the classes would move
    // them.
    CHECK(viaStruct.QN(1, 0) == doctest::Approx(viaDoc.QN(1, 0)).epsilon(1e-12));
    CHECK(viaStruct.TN(1, 0) == doctest::Approx(viaDoc.TN(1, 0)).epsilon(1e-12));
}

TEST_CASE("$LINE_LDES_DIR is authoritative, not merely first in the search") {
    // A caller who named a directory and named a wrong one must be told there is
    // no engine, not silently served the one beside the binary: that is the
    // difference between running the engine you asked for and running another.
    const std::size_t searched = ldes::detail::engine_dirs().size();
    ::setenv("LINE_LDES_DIR", "/nonexistent-ldes-dir", 1);
    const std::vector<std::string> only = ldes::detail::engine_dirs();
    const std::string found = ldes::detail::find_engine_dir_uncached();
    ::unsetenv("LINE_LDES_DIR");
    // The forced directory REPLACES the search rather than heading it, so the
    // ancestors of the executable are not consulted at all...
    CHECK(searched > 1);
    REQUIRE(only.size() == 1);
    CHECK(only[0] == "/nonexistent-ldes-dir");
    // ... and a directory with no engine in it therefore yields no engine.
    CHECK(found.empty());
}

TEST_CASE("ldes_cdf_respt builds the empirical CDF in [F, t] order") {
    // The column order is the one every getCdfRespT follows and the REVERSE of
    // getTranCdfRespT's; getting it backwards produces a CDF that reads as a
    // time axis with no error anywhere, so it is asserted directly.
    ldes::LdesResult r;
    r.respTimeSamples.resize(1);
    // A tie at 2.0 must collapse to ONE row carrying the LARGEST CDF value.
    std::vector<double> s;
    s.push_back(3.0);
    s.push_back(1.0);
    s.push_back(2.0);
    s.push_back(2.0);
    r.respTimeSamples[0].push_back(s);

    const Matrix<double> cdf = ldes::ldes_cdf_respt(r, 0, 0);
    REQUIRE(cdf.rows() == 3);
    REQUIRE(cdf.cols() == 2);
    // Column 1 is time, strictly increasing after the collapse.
    CHECK(cdf(0, 1) == doctest::Approx(1.0));
    CHECK(cdf(1, 1) == doctest::Approx(2.0));
    CHECK(cdf(2, 1) == doctest::Approx(3.0));
    // Column 0 is the CDF, reaching exactly 1, with the tie keeping 3/4.
    CHECK(cdf(0, 0) == doctest::Approx(0.25));
    CHECK(cdf(1, 0) == doctest::Approx(0.75));
    CHECK(cdf(2, 0) == doctest::Approx(1.0));

    // A pair with no observation yields an EMPTY matrix, not a fabricated law.
    CHECK(ldes::ldes_cdf_respt(r, 0, 5).rows() == 0);
    CHECK(ldes::ldes_cdf_respt(r, 9, 0).rows() == 0);
}
