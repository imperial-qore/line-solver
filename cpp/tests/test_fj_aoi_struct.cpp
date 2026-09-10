/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Structural and descriptor-building helpers of the fork-join and
 * age-of-information families. Oracles, in order of strength:
 *   1. Structure: the branches of an AND-join must PARTITION the activities
 *      strictly between the fork and the join, disjointly and exhaustively.
 *      That is a property of the answer, not a recorded number, and it is
 *      checked in exact arithmetic.
 *   2. Closed form: the mean of a PH built from a MAP must be the mean of the
 *      MAP. This is what exposes the aoi_dist2ph defect pinned below.
 *   3. MATLAB, for map_lambda, map_pie and the alpha the reference actually
 *      returns.
 */
#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

#include "doctest.h"
#include "line/api/aoi/aoi_dist2ph.h"
#include "line/api/fj/fj_branch_members.h"
#include "line/api/fj/fj_dist2fj.h"
#include "line/api/mam/map_moment.h"
#include "line/util/linalg.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::aoi::AoiPh;
using line::aoi::aoi_dist2ph;
using line::fj::FjDistKind;
using line::fj::FjProcType;
using line::fj::LqnBranchView;
using line::fj::fj_branch_members;
using line::fj::fj_dist2fj;
using line::mam::Map;

namespace {

constexpr double TOL = 1e-12;

/**
 * A two-branch AND-fork/join. Absolute indices: 1 and 2 are a task and an
 * entry, activities are 3..9 (ashift = 2, nacts = 7).
 *
 *   3 -> 4 -> 5 ------> 9
 *     -> 6 -> 7 -> 8 -> 9
 *
 * 4 and 6 are the spawned branch heads and carry POST_AND, which is how
 * lsn.actposttype records an AND-fork: getStruct writes the post type on the
 * POST activity, i.e. on the activity that was spawned, not on the forking
 * one. That is what makes the branches disjoint.
 */
template <class T>
LqnBranchView<T> forkJoinLqn() {
    using nt = line::num_traits<T>;
    LqnBranchView<T> lqn;
    lqn.ashift = 2;
    lqn.nacts = 7;
    lqn.graph = Matrix<T>(9, 9, nt::from_int(0));
    const int edges[7][2] = {{3, 4}, {3, 6}, {4, 5}, {5, 9}, {6, 7}, {7, 8}, {8, 9}};
    for (const auto& e : edges) lqn.graph(e[0] - 1, e[1] - 1) = nt::from_int(1);
    lqn.actposttype.assign(9, 0);
    lqn.actposttype[2] = line::fj::APC_POST_SEQ;  // activity 3
    lqn.actposttype[3] = line::fj::APC_POST_AND;  // activity 4, branch head
    lqn.actposttype[4] = line::fj::APC_POST_SEQ;
    lqn.actposttype[5] = line::fj::APC_POST_AND;  // activity 6, branch head
    lqn.actposttype[6] = line::fj::APC_POST_SEQ;
    lqn.actposttype[7] = line::fj::APC_POST_SEQ;
    lqn.actposttype[8] = line::fj::APC_POST_SEQ;  // activity 9, the join
    return lqn;
}

/** Erlang(2) of mean 1, written as a MAP. */
template <class T>
Map<T> erlang2() {
    using nt = line::num_traits<T>;
    Map<T> m;
    m.D0 = Matrix<T>{{nt::from_int(-2), nt::from_int(2)}, {nt::from_int(0), nt::from_int(-2)}};
    m.D1 = Matrix<T>{{nt::from_int(0), nt::from_int(0)}, {nt::from_int(2), nt::from_int(0)}};
    return m;
}

/** The two-phase MAP of the MATLAB reference run for fj_dist2fj. */
template <class T>
Map<T> map2() {
    using nt = line::num_traits<T>;
    Map<T> m;
    m.D0 = Matrix<T>{{nt::from_int(-2), nt::from_int(0)}, {nt::from_int(0), nt::from_int(-5)}};
    m.D1 = Matrix<T>{{nt::from_rational(6, 5), nt::from_rational(4, 5)},
                     {nt::from_int(2), nt::from_int(3)}};
    return m;
}

/** alpha (-T)^-1 e, the mean of the PH the conversion returns. */
template <class T>
T ph_mean(const AoiPh<T>& ph) {
    Matrix<T> negT = ph.Tmat;
    for (std::size_t i = 0; i < negT.rows(); ++i)
        for (std::size_t j = 0; j < negT.cols(); ++j) negT(i, j) = -negT(i, j);
    const std::vector<T> v = line::vecmul(ph.alpha, line::inverse(negT));
    T s = line::num_traits<T>::from_int(0);
    for (const T& x : v) s += x;
    return s;
}

}  // namespace

TEST_CASE("fj_branch_members partitions the activities between fork and join") {
    const LqnBranchView<Rational> lqn = forkJoinLqn<Rational>();
    const std::vector<std::vector<std::size_t>> br = fj_branch_members(lqn, 9);
    REQUIRE(br.size() == 2);
    // MATLAB: {[5 4], [8 7 6]}, tail first and branch head last
    CHECK(br[0] == std::vector<std::size_t>{5, 4});
    CHECK(br[1] == std::vector<std::size_t>{8, 7, 6});

    // the invariant: the branches are disjoint and cover exactly the
    // activities strictly between the fork and the join
    std::multiset<std::size_t> all;
    for (const std::vector<std::size_t>& b : br) all.insert(b.begin(), b.end());
    const std::set<std::size_t> uniq(all.begin(), all.end());
    CHECK(all.size() == uniq.size());  // disjoint
    CHECK(uniq == std::set<std::size_t>{4, 5, 6, 7, 8});
    // and every branch ends on an activity marked POST_AND, the spawned head
    for (const std::vector<std::size_t>& b : br)
        CHECK(lqn.actposttype[b.back() - 1] == line::fj::APC_POST_AND);
}

TEST_CASE("fj_branch_members ignores non-activity predecessors and stops at a merge") {
    LqnBranchView<double> lqn = forkJoinLqn<double>();
    // an entry (index 2, below ashift) also points at the join
    lqn.graph(1, 8) = 1.0;
    const std::vector<std::vector<std::size_t>> br = fj_branch_members(lqn, 9);
    CHECK(br.size() == 2);  // the entry is not an activity and is skipped

    // give activity 7 a second activity predecessor: the walk must stop there
    LqnBranchView<double> merged = forkJoinLqn<double>();
    merged.graph(4, 6) = 1.0;  // 5 -> 7, so 7 now has predecessors 6 and 5
    const std::vector<std::vector<std::size_t>> br2 = fj_branch_members(merged, 9);
    REQUIRE(br2.size() == 2);
    CHECK(br2[1] == std::vector<std::size_t>{8, 7});  // stopped before the merge
}

TEST_CASE("fj_branch_members rejects a malformed graph") {
    LqnBranchView<double> lqn = forkJoinLqn<double>();
    CHECK_THROWS_AS(fj_branch_members(lqn, 0), line::InputError);
    CHECK_THROWS_AS(fj_branch_members(lqn, 99), line::InputError);
    lqn.actposttype.resize(3);
    CHECK_THROWS_AS(fj_branch_members(lqn, 9), line::InputError);
}

TEST_CASE("fj_dist2fj builds the arrival and service descriptors") {
    const Map<double> m = map2<double>();
    const line::fj::FjDist<double> arr = fj_dist2fj(m, FjDistKind::Arrival, FjProcType::HyperExp);
    // MATLAB map_lambda and map_pie on the same MAP
    CHECK(arr.lambda == doctest::Approx(2.85714285714286).epsilon(1e-12));
    CHECK(arr.ma == 2);
    CHECK(arr.choice == 2);
    CHECK(arr.Ia(0, 0) == doctest::Approx(1.0));
    CHECK(arr.Ia(0, 1) == doctest::Approx(0.0));
    CHECK(arr.lambda1(0, 0) == doctest::Approx(1.2).epsilon(TOL));

    const line::fj::FjDist<double> svc = fj_dist2fj(m, FjDistKind::Service, FjProcType::Erlang);
    CHECK(svc.mu == doctest::Approx(2.85714285714286).epsilon(1e-12));
    CHECK(svc.choice == 3);
    CHECK(svc.tau_st[0] == doctest::Approx(0.5).epsilon(TOL));
    CHECK(svc.tau_st[1] == doctest::Approx(0.5).epsilon(TOL));
    // the exit vector is -D0 e
    CHECK(svc.St[0] == doctest::Approx(2.0).epsilon(TOL));
    CHECK(svc.St[1] == doctest::Approx(5.0).epsilon(TOL));
}

TEST_CASE("fj_dist2fj is exact and enforces the algorithm's own restrictions") {
    const Map<Rational> m = map2<Rational>();
    const line::fj::FjDist<Rational> arr = fj_dist2fj(m, FjDistKind::Arrival, FjProcType::Map);
    CHECK(arr.lambda == Rational(20, 7));  // exactly 20/7
    CHECK(arr.choice == 4);
    // a two-phase MAP is a legal arrival process but not a legal service one
    CHECK_THROWS_AS(fj_dist2fj(m, FjDistKind::Service, FjProcType::Map), line::InputError);
    // three phases are out of scope for the algorithm
    Map<Rational> big;
    big.D0 = Matrix<Rational>(3, 3, Rational(0));
    big.D1 = Matrix<Rational>(3, 3, Rational(0));
    CHECK_THROWS_AS(fj_dist2fj(big, FjDistKind::Arrival, FjProcType::Exp), line::InputError);

    // exponential, one phase: pie is trivially one and the exit rate is the rate
    Map<Rational> e;
    e.D0 = Matrix<Rational>{{Rational(-3)}};
    e.D1 = Matrix<Rational>{{Rational(3)}};
    const line::fj::FjDist<Rational> se = fj_dist2fj(e, FjDistKind::Service, FjProcType::Exp);
    CHECK(se.mu == Rational(3));
    CHECK(se.St[0] == Rational(3));
    CHECK(se.tau_st[0] == Rational(1));
    CHECK(se.choice == 1);
}

TEST_CASE("aoi_dist2ph matches MATLAB on a general MAP") {
    // The flow A of the traffic-merge tests, aggregated: D0 = [-3 1; 1 -4].
    Map<Rational> m;
    m.D0 = Matrix<Rational>{{Rational(-3), Rational(1)}, {Rational(1), Rational(-4)}};
    m.D1 = Matrix<Rational>{{Rational(3, 2), Rational(1, 2)}, {Rational(3, 2), Rational(3, 2)}};
    const AoiPh<Rational> ph = aoi_dist2ph(m);
    // MATLAB: alpha = [0.526315789473684 0.473684210526316] = [10/19 9/19]
    CHECK(ph.alpha[0] == Rational(10, 19));
    CHECK(ph.alpha[1] == Rational(9, 19));
    CHECK(ph.Tmat(0, 0) == Rational(-3));
    CHECK(ph.Tmat(1, 1) == Rational(-4));
    // and it is NOT map_pie, which MATLAB gives as [12/19 7/19]
    const std::vector<Rational> pie = line::mam::map_pie(m);
    CHECK(pie[0] == Rational(12, 19));
    CHECK(ph.alpha[0] != pie[0]);
}

TEST_CASE("aoi_dist2ph reproduces a reference defect: the entry vector is the EXIT phase") {
    // REFERENCE DEFECT, matlab/src/api/aoi/aoi_dist2ph.m lines 66-76.
    // alpha is built as theta .* (D1 e), the phase distribution where the
    // completion OCCURRED, and then used as the PH ENTRY vector. The entry
    // vector of a MAP viewed as a renewal process is map_pie = theta D1
    // normalized, the phase the process RESTARTS in. The two differ whenever
    // D1 moves probability between phases, and the error is not second order:
    // on Erlang(2) of mean 1 the reference returns alpha = [0 1], which is
    // Exp(2), so the mean of the PH it hands to the AoI solvers is 0.5 instead
    // of 1 and the SCV is 1 instead of 0.5. MATLAB reproduction:
    //   E = {[-2 2; 0 -2], [0 0; 2 0]};
    //   map_mean(E)                       -> 1
    //   [al, T] = aoi_dist2ph(E); al      -> [0 1]
    //   al * inv(-T) * ones(2,1)          -> 0.5
    //   map_pie(E) * inv(-T) * ones(2,1)  -> 1
    // The port reproduces the reference rather than silently correcting it,
    // because the aoi_* solvers are calibrated against it; the correct value
    // is asserted alongside so the gap cannot be lost.
    const Map<Rational> e = erlang2<Rational>();
    const AoiPh<Rational> ph = aoi_dist2ph(e);
    CHECK(ph.alpha[0] == Rational(0));
    CHECK(ph.alpha[1] == Rational(1));
    CHECK(ph_mean(ph) == Rational(1, 2));            // what the reference gives
    CHECK(line::mam::map_mean(e) == Rational(1));    // what the MAP actually is

    // With map_pie as the entry vector the mean comes out right, which is the
    // evidence that the defect is in alpha and not in T.
    AoiPh<Rational> fixed = ph;
    fixed.alpha = line::mam::map_pie(e);
    CHECK(ph_mean(fixed) == Rational(1));
}

TEST_CASE("aoi_dist2ph is a proper PH and agrees across arithmetics") {
    // Hyperexponential: here the exit phase and the entry phase coincide, so
    // the reference is right and the mean is preserved (MATLAB: 0.55).
    Map<double> h;
    h.D0 = Matrix<double>{{-1.0, 0.0}, {0.0, -10.0}};
    h.D1 = Matrix<double>{{0.5, 0.5}, {5.0, 5.0}};
    const AoiPh<double> ph = aoi_dist2ph(h);
    CHECK(ph.alpha[0] == doctest::Approx(0.5).epsilon(TOL));
    CHECK(ph.alpha[1] == doctest::Approx(0.5).epsilon(TOL));
    CHECK(ph_mean(ph) == doctest::Approx(0.55).epsilon(1e-12));
    CHECK(line::mam::map_mean(h) == doctest::Approx(0.55).epsilon(1e-12));
    double s = 0.0;
    for (double a : ph.alpha) s += a;
    CHECK(s == doctest::Approx(1.0).epsilon(1e-14));

    const AoiPh<Real50> phr = aoi_dist2ph(Map<Real50>{
        Matrix<Real50>{{Real50(-1), Real50(0)}, {Real50(0), Real50(-10)}},
        Matrix<Real50>{{Real50(0.5), Real50(0.5)}, {Real50(5), Real50(5)}}});
    CHECK(static_cast<double>(phr.alpha[0]) == doctest::Approx(0.5).epsilon(1e-14));

    // a positive diagonal is not a sub-generator
    Map<double> bad;
    bad.D0 = Matrix<double>{{1.0, 0.0}, {0.0, -1.0}};
    bad.D1 = Matrix<double>{{0.0, 0.0}, {1.0, 0.0}};
    CHECK_THROWS_AS(aoi_dist2ph(bad), line::InputError);
    CHECK_THROWS_AS(aoi_dist2ph(Matrix<double>(2, 2, 0.0), Matrix<double>(3, 3, 0.0)),
                    line::InputError);
}
