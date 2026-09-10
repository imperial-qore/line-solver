/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Maximum sustainable multiplicity of a layered software network. Oracles:
 *   1. The propagation law itself: no element ever passes on more concurrency
 *      than it received or than it can hold, except a setup task, which by
 *      design passes on its declared multiplicity. Asserted EXACTLY, since the
 *      whole computation is comparisons and additions of integers.
 *   2. A fan-in conserves what its callers hand it: the inflow of a node is
 *      the sum of the outflows of its predecessors, so the outflow of a node
 *      with unbounded multiplicity is exactly that sum.
 *   3. MATLAB reference values from matlab/src/api/lsn/lsn_max_multiplicity.m.
 */
#include <vector>

#include "doctest.h"
#include "line/api/lqn/lsn_max_multiplicity.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::lsn::LsnElementType;
using line::lsn::LsnInput;
using line::lsn::lsn_max_multiplicity;
using line::lsn::Multiplicity;

namespace {

/** The MATLAB reference chain: ref task -> entry -> task(2) -> entry. */
template <class T>
LsnInput<T> chain() {
    LsnInput<T> lsn;
    lsn.dag = Matrix<T>(4, 4, line::num_traits<T>::from_int(0));
    lsn.dag(0, 1) = line::num_traits<T>::from_int(1);
    lsn.dag(1, 2) = line::num_traits<T>::from_int(1);
    lsn.dag(2, 3) = line::num_traits<T>::from_int(1);
    lsn.mult = {Multiplicity<T>::of(5), Multiplicity<T>::inf(), Multiplicity<T>::of(2),
                Multiplicity<T>::inf()};
    lsn.type = {LsnElementType::TASK, LsnElementType::ENTRY, LsnElementType::TASK,
                LsnElementType::ENTRY};
    lsn.isref = {true, false, false, false};
    return lsn;
}

}  // namespace

TEST_CASE("lsn_max_multiplicity matches MATLAB on the reference chain") {
    // MATLAB LSN1 = [5 5 2 2]: the reference task seeds 5, the entry passes it
    // through unbounded, the downstream task clamps to its 2 threads.
    const auto out = lsn_max_multiplicity(chain<Rational>());
    REQUIRE(out.size() == 4);
    const long want[4] = {5, 5, 2, 2};
    for (std::size_t i = 0; i < 4; ++i) {
        CHECK_FALSE(out[i].infinite);
        CHECK(out[i].value == Rational(want[i]));
    }

    // the same at double and at high precision: the algorithm is integral
    const auto od = lsn_max_multiplicity(chain<double>());
    const auto or_ = lsn_max_multiplicity(chain<Real50>());
    for (std::size_t i = 0; i < 4; ++i) {
        CHECK(od[i].value == static_cast<double>(want[i]));
        CHECK(or_[i].value == Real50(want[i]));
    }
}

TEST_CASE("lsn_max_multiplicity exempts a setup task from the caller bound") {
    // MATLAB LSN2 = [5 5 7 7]: node 3 declares 7 threads and is a function
    // task, so the caller's 5 does not clamp it and the 7 propagates on.
    LsnInput<Rational> lsn = chain<Rational>();
    lsn.mult[2] = Multiplicity<Rational>::of(7);
    lsn.hassetup = {false, false, true, false};
    const auto out = lsn_max_multiplicity(lsn);
    const long want[4] = {5, 5, 7, 7};
    for (std::size_t i = 0; i < 4; ++i) {
        CHECK_FALSE(out[i].infinite);
        CHECK(out[i].value == Rational(want[i]));
    }

    // A setup task that nothing reaches is still bounded by its inflow of
    // zero: the exemption applies only once inflow is positive.
    LsnInput<Rational> iso = lsn;
    iso.dag(1, 2) = Rational(0);
    const auto out2 = lsn_max_multiplicity(iso);
    CHECK(out2[2].value == Rational(0));
}

TEST_CASE("lsn_max_multiplicity leaves an unbounded non-reference task unbounded") {
    // MATLAB LSN3 = [3 3 Inf]: the final pass overrides a non-reference TASK
    // of infinite multiplicity, whatever the caller bound was.
    LsnInput<Rational> lsn;
    lsn.dag = Matrix<Rational>(3, 3, Rational(0));
    lsn.dag(0, 1) = Rational(1);
    lsn.dag(1, 2) = Rational(1);
    lsn.mult = {Multiplicity<Rational>::of(3), Multiplicity<Rational>::inf(),
                Multiplicity<Rational>::inf()};
    lsn.type = {LsnElementType::TASK, LsnElementType::ENTRY, LsnElementType::TASK};
    lsn.isref = {true, false, false};

    const auto out = lsn_max_multiplicity(lsn);
    CHECK(out[0].value == Rational(3));
    CHECK_FALSE(out[0].infinite);
    CHECK(out[1].value == Rational(3));
    CHECK_FALSE(out[1].infinite);
    CHECK(out[2].infinite);
}

TEST_CASE("lsn_max_multiplicity seeds an entry that has an open arrival") {
    // MATLAB LSN_ARR = [1 1 1] against LSN_NOARR = [0 0 0]: an entry with an
    // open arrival is worth one thread even with no reference task anywhere.
    // The JAR omits this branch entirely and would return the second answer in
    // both cases; the port follows MATLAB, the reference implementation.
    LsnInput<Rational> lsn;
    lsn.dag = Matrix<Rational>(3, 3, Rational(0));
    lsn.dag(0, 1) = Rational(1);
    lsn.dag(1, 2) = Rational(1);
    lsn.mult = {Multiplicity<Rational>::inf(), Multiplicity<Rational>::of(4),
                Multiplicity<Rational>::inf()};
    lsn.type = {LsnElementType::ENTRY, LsnElementType::TASK, LsnElementType::ENTRY};
    lsn.isref = {false, false, false};

    LsnInput<Rational> noarr = lsn;
    const auto out0 = lsn_max_multiplicity(noarr);
    for (std::size_t i = 0; i < 3; ++i) CHECK(out0[i].value == Rational(0));

    lsn.entry_has_arrival = {true, false, false};
    const auto out1 = lsn_max_multiplicity(lsn);
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK_FALSE(out1[i].infinite);
        CHECK(out1[i].value == Rational(1));
    }
}

TEST_CASE("lsn_max_multiplicity conserves concurrency across a fan-in") {
    // Two reference tasks of 3 and 4 threads call the same entry, which has no
    // bound of its own: it must pass on exactly 7, the sum of what reaches it.
    // The entry then feeds a task of 5, which clamps to 5.
    LsnInput<Rational> lsn;
    lsn.dag = Matrix<Rational>(4, 4, Rational(0));
    lsn.dag(0, 2) = Rational(1);
    lsn.dag(1, 2) = Rational(1);
    lsn.dag(2, 3) = Rational(1);
    lsn.mult = {Multiplicity<Rational>::of(3), Multiplicity<Rational>::of(4),
                Multiplicity<Rational>::inf(), Multiplicity<Rational>::of(5)};
    lsn.type = {LsnElementType::TASK, LsnElementType::TASK, LsnElementType::ENTRY,
                LsnElementType::TASK};
    lsn.isref = {true, true, false, false};

    const auto out = lsn_max_multiplicity(lsn);
    CHECK(out[0].value == Rational(3));
    CHECK(out[1].value == Rational(4));
    CHECK(out[2].value == Rational(7));  // exactly the sum of the two callers
    CHECK(out[3].value == Rational(5));  // clamped by its own multiplicity

    // no element passes on more than it can hold
    CHECK(out[3].value <= Rational(5));
    CHECK(out[0].value <= Rational(3));
}

TEST_CASE("lsn_max_multiplicity pads a short multiplicity vector with infinity") {
    LsnInput<Rational> lsn = chain<Rational>();
    lsn.mult.resize(2);  // only the first two are declared
    const auto out = lsn_max_multiplicity(lsn);
    CHECK(out[0].value == Rational(5));
    CHECK(out[1].value == Rational(5));
    // node 2 is a non-reference TASK with an (implied) infinite multiplicity
    CHECK(out[2].infinite);
    CHECK(out[3].value == Rational(5));
}

TEST_CASE("lsn_max_multiplicity rejects a cyclic call graph") {
    // Kahn cannot order it; MATLAB then indexes with a zero and the JAR runs
    // off the end of its order list, so neither defines an answer.
    LsnInput<Rational> lsn;
    lsn.dag = Matrix<Rational>(2, 2, Rational(0));
    lsn.dag(0, 1) = Rational(1);
    lsn.dag(1, 0) = Rational(1);
    lsn.mult = {Multiplicity<Rational>::of(1), Multiplicity<Rational>::of(1)};
    lsn.type = {LsnElementType::TASK, LsnElementType::TASK};
    lsn.isref = {true, false};
    CHECK_THROWS_AS(lsn_max_multiplicity(lsn), line::InputError);

    LsnInput<Rational> bad = chain<Rational>();
    bad.isref.pop_back();
    CHECK_THROWS_AS(lsn_max_multiplicity(bad), line::InputError);
}
