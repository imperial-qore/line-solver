/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * MMAP algebra. The oracles are the identities that define a marked MAP: the
 * per-class matrices partition D1, class rates sum to the aggregate rate,
 * class probabilities sum to one, and the superposition of two Poisson streams
 * is Poisson at the summed rate. All are asserted exactly in the rational
 * instantiation, which is what makes them meaningful: a partition that fails
 * only in the 16th digit is still a broken partition.
 */
#include <vector>

#include "doctest.h"
#include "line/api/mam/mmap_lambda.h"

using line::Matrix;
using line::Rational;
using line::mam::Map;
using line::mam::Mmap;
using namespace line::mam;

namespace {

/** Poisson stream of rate lambda, split into two classes with weights w, 1-w. */
template <class T>
Mmap<T> marked_poisson(const T& lambda, const T& w) {
    Map<T> base = map_exponential(lambda);
    Matrix<T> weights(1, 2);
    weights(0, 0) = w;
    weights(0, 1) = line::num_traits<T>::from_int(1) - w;
    return mmap_mark(base, weights);
}

}  // namespace

TEST_CASE("a marked Poisson stream partitions its rate exactly") {
    const Rational lam(3, 2), w(1, 4);
    Mmap<Rational> m = marked_poisson(lam, w);
    CHECK(mmap_isfeasible(m));

    const std::vector<Rational> lk = mmap_count_lambda(m);
    CHECK(lk[0] == lam * w);
    CHECK(lk[1] == lam * (Rational(1) - w));
    CHECK(lk[0] + lk[1] == map_lambda(m.map()));

    const std::vector<Rational> pc = mmap_pc(m);
    CHECK(pc[0] == w);
    CHECK(pc[1] == Rational(1) - w);
    CHECK(pc[0] + pc[1] == Rational(1));
}

TEST_CASE("superposition of two Poisson streams is Poisson at the summed rate") {
    Mmap<Rational> a = marked_poisson(Rational(2), Rational(1, 2));
    Mmap<Rational> b = marked_poisson(Rational(3), Rational(1, 3));
    Mmap<Rational> s = mmap_super(a, b);

    CHECK(s.classes() == 4);
    CHECK(mmap_isfeasible(s));
    CHECK(map_lambda(s.map()) == Rational(5));
    // The superposition of two Poisson processes is Poisson: SCV is 1 exactly.
    CHECK(map_scv(s.map()) == Rational(1));

    // Class rates are preserved by the superposition.
    const std::vector<Rational> lk = mmap_count_lambda(s);
    CHECK(lk[0] == Rational(1));      // 2 * 1/2
    CHECK(lk[1] == Rational(1));      // 2 * 1/2
    CHECK(lk[2] == Rational(1));      // 3 * 1/3
    CHECK(lk[3] == Rational(2));      // 3 * 2/3
}

TEST_CASE("mmap_normalize rebuilds D1 and the diagonal from the class matrices") {
    Mmap<Rational> m = marked_poisson(Rational(1), Rational(1, 2));
    // Perturb one class matrix negatively: normalization must clamp it and
    // restore the generator property.
    m.Dc[0](0, 0) = Rational(-1, 10);
    Mmap<Rational> n = mmap_normalize(m);
    CHECK(n.Dc[0](0, 0) == Rational(0));
    CHECK(mmap_isfeasible(n));
}

TEST_CASE("mmap_scale sets the mean exactly and preserves the class split") {
    Mmap<Rational> m = marked_poisson(Rational(4), Rational(1, 4));
    Mmap<Rational> s = mmap_scale(m, Rational(2));
    CHECK(map_mean(s.map()) == Rational(2));
    CHECK(mmap_isfeasible(s));
    // Scaling time does not change which class an arrival belongs to.
    const std::vector<Rational> pc0 = mmap_pc(m);
    const std::vector<Rational> pc1 = mmap_pc(s);
    CHECK(pc0[0] == pc1[0]);
    CHECK(pc0[1] == pc1[1]);
}

TEST_CASE("MMAP quantities agree between double and exact arithmetic") {
    Mmap<Rational> mq = marked_poisson(Rational(3, 2), Rational(1, 4));
    Mmap<double> md = marked_poisson(1.5, 0.25);
    CHECK(static_cast<double>(mmap_count_lambda(mq)[0]) ==
          doctest::Approx(mmap_count_lambda(md)[0]).epsilon(1e-12));
    CHECK(static_cast<double>(mmap_pc(mq)[1]) == doctest::Approx(mmap_pc(md)[1]).epsilon(1e-12));
}
