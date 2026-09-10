/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `Mrg32k3a` and `JavaRandom` against the generators they reproduce.
 *
 * THE ORACLE IS THE JAVA ITSELF, not a published table: the expected values
 * below were printed at 17 significant digits by a program run on
 * `common/jline.jar`, which bundles SSJ, using the same six-long seed idiom
 * `Solver_ssj` uses (`seed + offset .. seed + offset + 5`, here 23000 + 1000)
 * and `new java.util.Random(seed + offset)`. If these ever drift, the C++ LDES
 * engine no longer walks the Java engine's sample path and every seeded golden
 * is answering about a different simulation.
 *
 * The comparison is EXACT (bit for bit) and not to a tolerance. Both generators
 * are deterministic integer/double recurrences with no transcendental step, so
 * anything short of equality is a transcription error, and a tolerance would
 * hide exactly the drift this test exists to catch. `nextGaussian` is the one
 * value with a `sqrt` and a `log` in it; it is still required to agree to the
 * last bit, since both sides evaluate the same two library functions on the
 * same arguments.
 */

#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/util/rng_ssj.h"

using namespace line;

TEST_CASE("rng ssj: MRG32k3a reproduces SSJ's stream at the engine's seed idiom") {
    rng::Mrg32k3a s;
    s.set_seed_offset(23000, 1000);

    // umontreal.ssj.rng.MRG32k3a.nextDouble(), first eight draws.
    const std::vector<double> expected = {
        0.024016438050991653, 0.12009737989405521, 0.16343270870717325,
        0.98398385259058360,  0.062844381451520920, 0.23155767055321383,
        0.71371115731352960,  0.52251325936120900};
    for (std::size_t i = 0; i < expected.size(); ++i) {
        const double got = s.next_double();
        CHECK(got == expected[i]);
    }

    // nextInt(0, 9) continues the same stream.
    const std::vector<int> expected_int = {9, 8, 6, 8, 8};
    for (std::size_t i = 0; i < expected_int.size(); ++i) CHECK(s.next_int(0, 9) == expected_int[i]);
}

TEST_CASE("rng ssj: MRG32k3a stays in (0,1) and moves its whole state") {
    rng::Mrg32k3a s;
    s.set_seed_offset(1, 0);
    double lo = 1.0, hi = 0.0;
    for (int i = 0; i < 20000; ++i) {
        const double u = s.next_double();
        REQUIRE(u > 0.0);
        REQUIRE(u < 1.0);
        lo = std::min(lo, u);
        hi = std::max(hi, u);
    }
    // 20000 draws must have visited both tails; a stuck component would not.
    CHECK(lo < 0.001);
    CHECK(hi > 0.999);
}

TEST_CASE("rng ssj: an invalid seed is refused rather than silently repaired") {
    rng::Mrg32k3a s;
    long long zeros[6] = {0, 0, 0, 1, 2, 3};
    CHECK_THROWS_AS(s.set_seed(zeros), InputError);
    long long huge[6] = {4294967087LL, 1, 2, 3, 4, 5};  // == m1, out of range
    CHECK_THROWS_AS(s.set_seed(huge), InputError);
}

TEST_CASE("rng ssj: JavaRandom reproduces java.util.Random") {
    rng::JavaRandom r(23000 + 1000);
    const std::vector<double> expected = {
        0.75058965563818650, 0.73440928437475580, 0.97239811151721100,
        0.82317860781782020, 0.18260639907486530, 0.53753197812021680};
    for (std::size_t i = 0; i < expected.size(); ++i) CHECK(r.next_double() == expected[i]);

    // nextInt(7) exercises the rejection loop for a non-power-of-two bound.
    const std::vector<int> expected_int = {4, 3, 4, 4, 0, 0};
    for (std::size_t i = 0; i < expected_int.size(); ++i) CHECK(r.next_int(7) == expected_int[i]);

    // nextGaussian caches its second value; both must match, in order.
    const std::vector<double> expected_gauss = {-0.41551487558510697, 0.90504256441645510,
                                                0.29925929799642260, 0.55204324572790770};
    for (std::size_t i = 0; i < expected_gauss.size(); ++i)
        CHECK(r.next_gaussian() == expected_gauss[i]);

    const std::vector<long long> expected_long = {-6236566358722511993LL, 3148705076811121167LL,
                                                  6733975611032363342LL};
    for (std::size_t i = 0; i < expected_long.size(); ++i) CHECK(r.next_long() == expected_long[i]);
}

TEST_CASE("rng ssj: a power-of-two bound takes the other branch and stays in range") {
    rng::JavaRandom r(7);
    for (int i = 0; i < 1000; ++i) {
        const int v = r.next_int(8);
        REQUIRE(v >= 0);
        REQUIRE(v < 8);
    }
    CHECK_THROWS_AS(r.next_int(0), InputError);
}
