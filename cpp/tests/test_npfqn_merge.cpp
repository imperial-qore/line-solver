/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Superposition of marked arrival flows, with and without class switching.
 * Oracles, in order of strength:
 *   1. The conservation identity the superposition must satisfy exactly: the
 *      merged per-class rate is the sum of the operands' per-class rates, and
 *      under class switching it is the switching matrix applied to them. Both
 *      are checked for EXACT equality in the rational instantiation, where a
 *      passing test means the identity holds and not that it holds to within a
 *      tolerance.
 *   2. MATLAB, run on the same two flows (mmap_super with 'match',
 *      npfqn_traffic_merge, npfqn_traffic_merge_cs, mmap_compress): every
 *      descriptor entry and every derived rate below is a MATLAB number.
 *   3. What compression claims to preserve, checked against what it destroys:
 *      the first three aggregate moments and the class probabilities survive,
 *      the order does not.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_compress.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/npfqn/npfqn_traffic_merge.h"
#include "line/api/npfqn/npfqn_traffic_merge_cs.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::mam::Mmap;
using line::npfqn::Compress;
using line::npfqn::Merge;
using line::npfqn::MergeConfig;
using line::npfqn::mmap_mark_types;
using line::npfqn::mmap_super_match;
using line::npfqn::npfqn_traffic_merge;
using line::npfqn::npfqn_traffic_merge_cs;

namespace {

/** Flow A: two phases, two classes. Every entry is dyadic, hence exact. */
template <class T>
Mmap<T> flowA() {
    using nt = line::num_traits<T>;
    Mmap<T> m;
    m.D0 = Matrix<T>{{nt::from_double(-3), nt::from_double(1)},
                     {nt::from_double(1), nt::from_double(-4)}};
    Matrix<T> D1a{{nt::from_double(1), nt::from_double(0.5)},
                  {nt::from_double(1), nt::from_double(1)}};
    Matrix<T> D1b{{nt::from_double(0.5), nt::from_double(0)},
                  {nt::from_double(0.5), nt::from_double(0.5)}};
    m.D1 = Matrix<T>(2, 2, nt::from_int(0));
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) m.D1(i, j) = D1a(i, j) + D1b(i, j);
    m.Dc.push_back(D1a);
    m.Dc.push_back(D1b);
    return m;
}

/** Flow B: two phases, two classes. */
template <class T>
Mmap<T> flowB() {
    using nt = line::num_traits<T>;
    Mmap<T> m;
    m.D0 = Matrix<T>{{nt::from_double(-2), nt::from_double(0.5)},
                     {nt::from_double(0.25), nt::from_double(-1.25)}};
    Matrix<T> D1a{{nt::from_double(1), nt::from_double(0)},
                  {nt::from_double(0.5), nt::from_double(0.25)}};
    Matrix<T> D1b{{nt::from_double(0.25), nt::from_double(0.25)},
                  {nt::from_double(0.25), nt::from_double(0)}};
    m.D1 = Matrix<T>(2, 2, nt::from_int(0));
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) m.D1(i, j) = D1a(i, j) + D1b(i, j);
    m.Dc.push_back(D1a);
    m.Dc.push_back(D1b);
    return m;
}

/** The (4 x 2) class-switching matrix of the MATLAB reference run. */
template <class T>
Matrix<T> switchProb() {
    using nt = line::num_traits<T>;
    return Matrix<T>{{nt::from_rational(7, 10), nt::from_rational(3, 10)},
                     {nt::from_rational(2, 10), nt::from_rational(8, 10)},
                     {nt::from_rational(5, 10), nt::from_rational(5, 10)},
                     {nt::from_rational(1, 10), nt::from_rational(9, 10)}};
}

constexpr double TOL = 1e-11;

}  // namespace

TEST_CASE("mmap_super_match adds the per-class rates exactly") {
    const Mmap<Rational> A = flowA<Rational>(), B = flowB<Rational>();
    const Mmap<Rational> S = mmap_super_match(A, B);
    REQUIRE(S.order() == 4);
    REQUIRE(S.classes() == 2);
    const std::vector<Rational> la = line::mam::mmap_lambda(A);
    const std::vector<Rational> lb = line::mam::mmap_lambda(B);
    const std::vector<Rational> ls = line::mam::mmap_lambda(S);
    for (std::size_t c = 0; c < 2; ++c) CHECK(ls[c] == Rational(la[c] + lb[c]));

    // and the descriptor stays a proper MMAP: the class matrices partition D1
    CHECK(line::mam::mmap_isfeasible(S));
}

TEST_CASE("npfqn_traffic_merge without compression matches MATLAB") {
    const Mmap<double> S =
        npfqn_traffic_merge<double>({flowA<double>(), flowB<double>()},
                                    MergeConfig{Merge::Default, Compress::None});
    REQUIRE(S.order() == 4);
    // MATLAB: mmap_lambda(S)
    const std::vector<double> l = line::mam::mmap_lambda(S);
    CHECK(l[0] == doctest::Approx(2.58035714285714).epsilon(TOL));
    CHECK(l[1] == doctest::Approx(1.08035714285714).epsilon(TOL));
    const line::mam::Map<double> m = S.map();
    CHECK(line::mam::map_lambda(m) == doctest::Approx(3.66071428571429).epsilon(TOL));
    CHECK(line::mam::map_mean(m) == doctest::Approx(0.273170731707317).epsilon(TOL));
    CHECK(line::mam::map_scv(m) == doctest::Approx(1.01762864070959).epsilon(TOL));
    // D0 = krons(D0a, D0b), first row by hand
    CHECK(S.D0(0, 0) == doctest::Approx(-5.0).epsilon(TOL));
    CHECK(S.D0(0, 1) == doctest::Approx(0.5).epsilon(TOL));
    CHECK(S.D0(0, 2) == doctest::Approx(1.0).epsilon(TOL));
    CHECK(S.D0(0, 3) == doctest::Approx(0.0).epsilon(TOL));
}

TEST_CASE("npfqn_traffic_merge single flow is the identity") {
    const Mmap<Rational> A = flowA<Rational>();
    const Mmap<Rational> S = npfqn_traffic_merge<Rational>({A}, MergeConfig{});
    REQUIRE(S.order() == A.order());
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) CHECK(S.D0(i, j) == A.D0(i, j));
    // NOTE this is the case the reference's missing-compress defect does not
    // reach, because MATLAB returns before the switch.
}

TEST_CASE("npfqn_traffic_merge compressed matches MATLAB and keeps three moments") {
    const Mmap<double> S =
        npfqn_traffic_merge<double>({flowA<double>(), flowB<double>()},
                                    MergeConfig{Merge::Default, Compress::None});
    const Mmap<double> C =
        npfqn_traffic_merge<double>({flowA<double>(), flowB<double>()},
                                    MergeConfig{Merge::Default, Compress::Default});
    REQUIRE(C.order() == 4);  // one APH(2) per class

    // MATLAB npfqn_traffic_merge with compress = 'default'
    const std::vector<double> l = line::mam::mmap_lambda(C);
    CHECK(l[0] == doctest::Approx(2.58035714285714).epsilon(1e-9));
    CHECK(l[1] == doctest::Approx(1.08035714285714).epsilon(1e-9));
    const std::vector<double> pc = line::mam::mmap_pc(C);
    CHECK(pc[0] == doctest::Approx(0.704878048780488).epsilon(1e-9));
    CHECK(pc[1] == doctest::Approx(0.295121951219512).epsilon(1e-9));

    // the compressed descriptor itself, entry by entry against MATLAB
    CHECK(C.D0(0, 0) == doctest::Approx(-4.43892271651221).epsilon(1e-9));
    CHECK(C.D0(0, 1) == doctest::Approx(0.75176856161451).epsilon(1e-9));
    CHECK(C.D0(1, 1) == doctest::Approx(-3.48540668074399).epsilon(1e-9));
    CHECK(C.D0(2, 2) == doctest::Approx(-4.42110236195367).epsilon(1e-9));
    CHECK(C.D0(2, 3) == doctest::Approx(0.697514436200139).epsilon(1e-9));
    CHECK(C.D0(3, 3) == doctest::Approx(-3.48190163670019).epsilon(1e-9));
    CHECK(C.Dc[0](0, 0) == doctest::Approx(2.59899402625716).epsilon(1e-9));
    CHECK(C.Dc[0](0, 2) == doctest::Approx(1.08816012864054).epsilon(1e-9));
    CHECK(C.Dc[0](1, 0) == doctest::Approx(2.4567866603293).epsilon(1e-9));
    // class 1 can only be marked on a departure from its own component block
    CHECK(C.Dc[0](2, 0) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(C.Dc[0](3, 3) == doctest::Approx(0.0).epsilon(1e-12));

    // what compression preserves: the first three aggregate moments
    const std::vector<unsigned> ord{1u, 2u, 3u};
    for (unsigned k : ord) {
        const double before = line::mam::map_moment(S.map(), k);
        const double after = line::mam::map_moment(C.map(), k);
        INFO("moment ", k);
        CHECK(after == doctest::Approx(before).epsilon(1e-9));
    }
    // and what it destroys: the order collapses from 4 phases of a product
    // chain to 4 phases of a renewal mixture, so the autocorrelation goes
    // map_acf now carries map_acf.m's own (x - 1)/scv normalization, so the
    // renewal value is 0 and the hand-applied correction that used to stand
    // here is gone.
    const std::vector<unsigned> lag1{1u};
    const double acfC = line::mam::map_acf(C.map(), lag1)[0];
    const double acfS = line::mam::map_acf(S.map(), lag1)[0];
    CHECK(std::abs(acfC) < 1e-12);
    CHECK(std::abs(acfS) > 1e-6);
}

TEST_CASE("npfqn_traffic_merge refuses what is not ported and what is malformed") {
    const std::vector<Mmap<double>> flows{flowA<double>(), flowB<double>()};
    CHECK_THROWS_AS(npfqn_traffic_merge<double>(flows, MergeConfig{Merge::Mixture, Compress::None}),
                    line::UnsupportedError);
    // 'interpos' IS ported now (m3pp2m_fitc_theoretical + m3pp2m_interleave),
    // so what it refuses here is the INPUT and not the method: flow A fits to
    // an M3PP of order 1 -- MATLAB's m3pp2m_fitc_theoretical returns order 1 on
    // this flow too -- and the interleaving is defined for order-2 components
    // only. The reference dies with "Index in position 2 exceeds array bounds"
    // on the same pair; refusing by name is the port's improvement on it, and
    // an InputError is the right kind because the domain is the caller's.
    CHECK_THROWS_AS(npfqn_traffic_merge<double>(flows, MergeConfig{Merge::Interpos, Compress::None}),
                    line::InputError);
    // At exact arithmetic the counting-process fitter is unavailable, and THAT
    // refusal is about the method rather than the input.
    CHECK_THROWS_AS(npfqn_traffic_merge<Rational>({flowA<Rational>(), flowB<Rational>()},
                                                  MergeConfig{Merge::Interpos, Compress::None}),
                    line::UnsupportedError);
    // compression is unavailable in exact arithmetic, and says so
    CHECK_THROWS_AS(npfqn_traffic_merge<Rational>({flowA<Rational>(), flowB<Rational>()},
                                                  MergeConfig{Merge::Default, Compress::Default}),
                    line::UnsupportedError);
    CHECK_THROWS_AS(npfqn_traffic_merge<double>({}, MergeConfig{}), line::InputError);
    // a flow with a different class count cannot be matched
    Mmap<double> odd = flowA<double>();
    odd.Dc.pop_back();
    CHECK_THROWS_AS(mmap_super_match(flowA<double>(), odd), line::InputError);
}

TEST_CASE("npfqn_traffic_merge_cs carries the rates through the switching matrix") {
    const Mmap<Rational> A = flowA<Rational>(), B = flowB<Rational>();
    const Matrix<Rational> P = switchProb<Rational>();
    const Mmap<Rational> S = npfqn_traffic_merge_cs<Rational>({A, B}, P);
    REQUIRE(S.order() == 4);
    REQUIRE(S.classes() == 2);

    const std::vector<Rational> la = line::mam::mmap_lambda(A);
    const std::vector<Rational> lb = line::mam::mmap_lambda(B);
    const std::vector<Rational> ls = line::mam::mmap_lambda(S);
    for (std::size_t s = 0; s < 2; ++s) {
        Rational want = Rational(0);
        for (std::size_t r = 0; r < 2; ++r) want += la[r] * P(r, s) + lb[r] * P(2 + r, s);
        CHECK(ls[s] == want);  // exact, not to a tolerance
    }
    // a stochastic switching matrix preserves the aggregate rate
    Rational total = Rational(0), operands = Rational(0);
    for (std::size_t s = 0; s < 2; ++s) total += ls[s];
    for (std::size_t r = 0; r < 2; ++r) operands += la[r] + lb[r];
    CHECK(total == operands);
}

TEST_CASE("npfqn_traffic_merge_cs matches MATLAB entry by entry") {
    const Mmap<double> S =
        npfqn_traffic_merge_cs<double>({flowA<double>(), flowB<double>()}, switchProb<double>());
    const double D0[4][4] = {{-5, 0.5, 1, 0}, {0.25, -4.25, 0, 1}, {1, 0, -6, 0.5}, {0, 1, 0.25, -5.25}};
    const double D1a[4][4] = {{1.325, 0.025, 0.35, 0},
                              {0.275, 0.925, 0, 0.35},
                              {0.8, 0, 1.325, 0.025},
                              {0, 0.8, 0.275, 0.925}};
    const double D1b[4][4] = {{1.425, 0.225, 0.15, 0},
                              {0.475, 0.825, 0, 0.15},
                              {0.7, 0, 1.425, 0.225},
                              {0, 0.7, 0.475, 0.825}};
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) {
            INFO("entry ", i, ",", j);
            CHECK(S.D0(i, j) == doctest::Approx(D0[i][j]).epsilon(TOL));
            CHECK(S.Dc[0](i, j) == doctest::Approx(D1a[i][j]).epsilon(TOL));
            CHECK(S.Dc[1](i, j) == doctest::Approx(D1b[i][j]).epsilon(TOL));
        }
    const std::vector<double> l = line::mam::mmap_lambda(S);
    CHECK(l[0] == doctest::Approx(1.80446428571429).epsilon(TOL));
    CHECK(l[1] == doctest::Approx(1.85625).epsilon(TOL));
}

TEST_CASE("npfqn_traffic_merge_cs single flow is a pure re-marking") {
    const Mmap<Rational> A = flowA<Rational>();
    Matrix<Rational> P(2, 2);
    P(0, 0) = Rational(1, 4);
    P(0, 1) = Rational(3, 4);
    P(1, 0) = Rational(1, 2);
    P(1, 1) = Rational(1, 2);
    const Mmap<Rational> S = npfqn_traffic_merge_cs<Rational>({A}, P);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) {
            CHECK(S.D0(i, j) == A.D0(i, j));  // the hidden process is untouched
            CHECK(Rational(S.Dc[0](i, j) + S.Dc[1](i, j)) == A.D1(i, j));
        }
    CHECK_THROWS_AS(npfqn_traffic_merge_cs<Rational>({A, A}, P), line::InputError);
}

TEST_CASE("npfqn_traffic_merge exact and double agree to rounding") {
    const Mmap<Rational> Se = npfqn_traffic_merge<Rational>(
        {flowA<Rational>(), flowB<Rational>()}, MergeConfig{Merge::Super, Compress::None});
    const Mmap<double> Sd = npfqn_traffic_merge<double>(
        {flowA<double>(), flowB<double>()}, MergeConfig{Merge::Super, Compress::None});
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j)
            CHECK(Sd.D0(i, j) == doctest::Approx(static_cast<double>(Se.D0(i, j))).epsilon(1e-14));
    const Mmap<Real50> Sr = npfqn_traffic_merge<Real50>(
        {flowA<Real50>(), flowB<Real50>()}, MergeConfig{Merge::Super, Compress::None});
    const std::vector<Real50> lr = line::mam::mmap_lambda(Sr);
    const std::vector<Rational> le = line::mam::mmap_lambda(Se);
    for (std::size_t c = 0; c < 2; ++c)
        CHECK(static_cast<double>(lr[c]) ==
              doctest::Approx(static_cast<double>(le[c])).epsilon(1e-14));
}
