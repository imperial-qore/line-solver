/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Maximum-entropy (Kouvatsos) network algorithms. The oracles here are the
 * queues whose exact solution the ME building blocks are known to reproduce:
 * a GE/GE/1 station with Ca = Cs = 1 IS an M/M/1 and its departure scv must
 * come back as 1 (Burke), the same station with Cs = 0 must reproduce the
 * Pollaczek-Khinchine mean of an M/D/1, an insensitive station must return
 * the product-form mean queue length whatever the service scv, and the
 * GE/GE/c block of eq. (3.9) at Ca = Cs = 1 must reproduce M/M/c exactly.
 * For the closed algorithm the invariant is the population constraint plus
 * Little's law, and the reference values were cross-checked against the
 * MATLAB me_cqn.
 */
#include <vector>

#include "doctest.h"
#include "line/api/me/me_cqn.h"
#include "line/api/me/me_mqn.h"
#include "line/api/me/me_oqn.h"
#include "line/util/error.h"

using line::Matrix;
using namespace line::me;

namespace {

/** A single station fed by a Poisson stream, no internal routing. */
MeResult<double> single_station(double lambda, double mu, double Cs, long servers, bool insens) {
    Matrix<double> lambda0(1, 1, lambda), Ca0(1, 1, 1.0);
    Matrix<double> M(1, 1, mu), C(1, 1, Cs);
    std::vector<Matrix<double>> P(1, Matrix<double>(1, 1, 0.0));
    std::vector<long> c(1, servers);
    std::vector<char> ins(1, insens ? 1 : 0);
    return me_oqn(1u, 1u, lambda0, Ca0, M, C, P, c, ins);
}

}  // namespace

TEST_CASE("the GE/GE/1 block reproduces M/M/1 and its Burke departure process") {
    const MeResult<double> r = single_station(0.5, 1.0, 1.0, 1, false);
    CHECK(r.converged);
    CHECK(r.rho(0, 0) == doctest::Approx(0.5));
    CHECK(r.L(0, 0) == doctest::Approx(1.0).epsilon(1e-9));   // rho/(1-rho)
    CHECK(r.W(0, 0) == doctest::Approx(2.0).epsilon(1e-9));   // Little
    CHECK(r.Ca(0, 0) == doctest::Approx(1.0).epsilon(1e-9));  // Poisson input
    CHECK(r.Cd(0, 0) == doctest::Approx(1.0).epsilon(1e-9));  // Burke: Poisson output
    CHECK(r.X[0] == doctest::Approx(0.5));
}

TEST_CASE("the GE/GE/1 block reproduces the Pollaczek-Khinchine mean") {
    // M/D/1 at rho = 1/2: L = rho + rho^2 (1 + Cs) / (2 (1-rho)) = 3/4
    const MeResult<double> d = single_station(0.5, 1.0, 0.0, 1, false);
    CHECK(d.L(0, 0) == doctest::Approx(0.75).epsilon(1e-9));
    // M/H2/1 with Cs = 3 at the same load: L = 0.5 + 0.25*4/1 = 1.5
    const MeResult<double> h = single_station(0.5, 1.0, 3.0, 1, false);
    CHECK(h.L(0, 0) == doctest::Approx(1.5).epsilon(1e-9));
    // The mean grows with the service variability, as it must
    CHECK(d.L(0, 0) < h.L(0, 0));
}

TEST_CASE("an insensitive station returns the product-form mean queue length") {
    // M/G/1-PS is insensitive: L = rho/(1-rho) for every service scv
    const MeResult<double> a = single_station(0.5, 1.0, 0.0, 1, true);
    const MeResult<double> b = single_station(0.5, 1.0, 9.0, 1, true);
    CHECK(a.L(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(b.L(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("an infinite server station holds lambda/mu jobs") {
    const MeResult<double> r = single_station(2.0, 4.0, 5.0, 0, false);
    CHECK(r.L(0, 0) == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(r.W(0, 0) == doctest::Approx(0.25).epsilon(1e-9));  // the service time
    CHECK(r.Cd(0, 0) == doctest::Approx(r.Ca(0, 0)));         // no queueing distortion
}

TEST_CASE("the GE/GE/c block of eq. (3.9) reproduces M/M/2") {
    // a = lambda/mu = 1, rho = 1/2: Erlang C gives Lq = 1/3 and L = 4/3
    const MeResult<double> r = single_station(1.0, 1.0, 1.0, 2, false);
    CHECK(r.rho(0, 0) == doctest::Approx(0.5));
    CHECK(r.L(0, 0) == doctest::Approx(4.0 / 3.0).epsilon(1e-9));
    // and M/M/3 at the same offered load a = 1: P0 = 4/11, Lq = 1/22
    Matrix<double> lambda0(1, 1, 1.0), Ca0(1, 1, 1.0), mu(1, 1, 1.0), Cs(1, 1, 1.0);
    std::vector<Matrix<double>> P(1, Matrix<double>(1, 1, 0.0));
    std::vector<long> c(1, 3);
    std::vector<char> ins(1, 0);
    const MeResult<double> r3 = me_oqn(1u, 1u, lambda0, Ca0, mu, Cs, P, c, ins);
    CHECK(r3.L(0, 0) == doctest::Approx(1.0 + 1.0 / 22.0).epsilon(1e-9));
}

TEST_CASE("the feedback correction turns a self-loop into a slower server") {
    // A station with p_ii = 1/2 and mu = 2 serves a job an expected two times,
    // so the composite service is exponential with rate 1: the station is the
    // M/M/1 of the first test again, at the same effective load.
    Matrix<double> lambda0(1, 1, 0.5), Ca0(1, 1, 1.0), mu(1, 1, 2.0), Cs(1, 1, 1.0);
    std::vector<Matrix<double>> P(1, Matrix<double>(1, 1, 0.5));
    std::vector<long> c(1, 1);
    std::vector<char> ins(1, 0);
    const MeResult<double> r = me_oqn(1u, 1u, lambda0, Ca0, mu, Cs, P, c, ins);
    CHECK(r.lambda(0, 0) == doctest::Approx(1.0));  // visit inclusive
    CHECK(r.rho(0, 0) == doctest::Approx(0.5));
    CHECK(r.L(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("an open tandem conserves flow and rejects an unstable load") {
    // Two exponential stations in series, the second twice as fast
    Matrix<double> lambda0(2, 2, 0.0), Ca0(2, 2, 1.0), mu(2, 2, 0.0), Cs(2, 2, 1.0);
    lambda0(0, 0) = 0.4;
    lambda0(0, 1) = 0.2;
    mu(0, 0) = 1.0;
    mu(0, 1) = 1.0;
    mu(1, 0) = 2.0;
    mu(1, 1) = 2.0;
    std::vector<Matrix<double>> P(2, Matrix<double>(2, 2, 0.0));
    P[0](0, 1) = 1.0;
    P[1](0, 1) = 1.0;
    std::vector<long> c(2, 1);
    std::vector<char> ins(2, 0);
    const MeResult<double> r = me_oqn(2u, 2u, lambda0, Ca0, mu, Cs, P, c, ins);
    CHECK(r.converged);
    for (std::size_t k = 0; k < 2; ++k) {
        CHECK(r.lambda(0, k) == doctest::Approx(lambda0(0, k)));
        CHECK(r.lambda(1, k) == doctest::Approx(lambda0(0, k)));  // flow conservation
        CHECK(r.W(0, k) == doctest::Approx(r.L(0, k) / r.lambda(0, k)));
    }
    // The first station carries the whole exponential load, so it is an M/M/1
    // with rho = 0.6 and the two classes share L = rho/(1-rho) = 1.5
    CHECK(r.L(0, 0) + r.L(0, 1) == doctest::Approx(1.5).epsilon(1e-6));

    // An overloaded station is refused rather than reported as Inf
    Matrix<double> over = lambda0;
    over(0, 0) = 1.5;
    CHECK_THROWS_AS(me_oqn(2u, 2u, over, Ca0, mu, Cs, P, c, ins), line::NumericError);
}

TEST_CASE("the closed algorithm conserves the population and Little's law") {
    // One queue and one delay, exponential, N = 3, D = Z = 1. The exact
    // product-form solution is X = 0.9375, Lqueue = 2.0625, Ldelay = 0.9375.
    const std::size_t M = 2, R = 1;
    Matrix<double> mu(M, R, 1.0), Cs(M, R, 1.0);
    std::vector<Matrix<double>> P(R, Matrix<double>(M, M, 0.0));
    P[0](0, 1) = 1.0;
    P[0](1, 0) = 1.0;
    std::vector<long> c(M, 1);
    c[1] = 0;  // the delay
    std::vector<char> ins(M, 0);
    std::vector<long> N(R, 3), ref(R, 0);
    const MeResult<double> r = me_cqn(M, R, N, mu, Cs, P, c, ref, ins);

    // The population constraint is what the two stages exist to enforce
    CHECK(r.L(0, 0) + r.L(1, 0) == doctest::Approx(3.0).epsilon(1e-6));
    // Every queue length is nonnegative and the delay holds X*Z jobs
    CHECK(r.L(0, 0) > 0.0);
    CHECK(r.L(1, 0) == doctest::Approx(r.X[0]).epsilon(1e-6));
    // Little's law on the visit-inclusive throughputs
    for (std::size_t i = 0; i < M; ++i)
        CHECK(r.W(i, 0) == doctest::Approx(r.L(i, 0) / r.lambda(i, 0)).epsilon(1e-9));
    // The utilization of the queue is X * D
    CHECK(r.rho(0, 0) == doctest::Approx(r.X[0]).epsilon(1e-6));
    // and the throughput is within 10 percent of the exact 0.9375: the ME
    // solution is an approximation, but a GE-type one that is exact for the
    // building blocks, so it must not be far
    CHECK(r.X[0] > 0.85);
    CHECK(r.X[0] < 1.0);
}

TEST_CASE("a closed network of two queues splits its population sensibly") {
    // Two single-server exponential queues in a cycle, one twice as fast.
    const std::size_t M = 2, R = 1;
    Matrix<double> mu(M, R, 1.0), Cs(M, R, 1.0);
    mu(1, 0) = 2.0;
    std::vector<Matrix<double>> P(R, Matrix<double>(M, M, 0.0));
    P[0](0, 1) = 1.0;
    P[0](1, 0) = 1.0;
    std::vector<long> c(M, 1);
    std::vector<char> ins(M, 0);
    std::vector<long> N(R, 4), ref(R, 0);
    const MeResult<double> r = me_cqn(M, R, N, mu, Cs, P, c, ref, ins);
    CHECK(r.L(0, 0) + r.L(1, 0) == doctest::Approx(4.0).epsilon(1e-6));
    // the slow station is the bottleneck and holds the larger share
    CHECK(r.L(0, 0) > r.L(1, 0));
    // utilizations are ordered the same way and stay below one
    CHECK(r.rho(0, 0) > r.rho(1, 0));
    CHECK(r.rho(0, 0) <= 1.0);
}

TEST_CASE("the mixed algorithm degenerates to its two components") {
    // With no closed class me_mqn must return exactly the open solution
    Matrix<double> lambda0(1, 1, 0.5), Ca0(1, 1, 1.0), mu(1, 1, 1.0), Cs(1, 1, 1.0);
    std::vector<Matrix<double>> P(1, Matrix<double>(1, 1, 0.0));
    std::vector<long> c(1, 1), N(1, 0), ref(1, -1);
    std::vector<char> ins(1, 0), open(1, 1);
    const MeResult<double> mo = me_mqn(1u, 1u, open, lambda0, Ca0, N, mu, Cs, P, c, ref, ins);
    const MeResult<double> oo = me_oqn(1u, 1u, lambda0, Ca0, mu, Cs, P, c, ins);
    CHECK(mo.L(0, 0) == doctest::Approx(oo.L(0, 0)).epsilon(1e-12));
    CHECK(mo.Cd(0, 0) == doctest::Approx(oo.Cd(0, 0)).epsilon(1e-12));
    CHECK(mo.X[0] == doctest::Approx(0.5));

    // With no open class it must return exactly the closed solution
    const std::size_t M = 2;
    Matrix<double> mu2(M, 1, 1.0), Cs2(M, 1, 1.0), lam2(M, 1, 0.0), Ca2(M, 1, 1.0);
    std::vector<Matrix<double>> P2(1, Matrix<double>(M, M, 0.0));
    P2[0](0, 1) = 1.0;
    P2[0](1, 0) = 1.0;
    std::vector<long> c2(M, 1), N2(1, 2), ref2(1, 0);
    c2[1] = 0;
    std::vector<char> ins2(M, 0), open2(1, 0);
    const MeResult<double> mc = me_mqn(M, 1u, open2, lam2, Ca2, N2, mu2, Cs2, P2, c2, ref2, ins2);
    const MeResult<double> cc = me_cqn(M, 1u, N2, mu2, Cs2, P2, c2, ref2, ins2);
    CHECK(mc.L(0, 0) == doctest::Approx(cc.L(0, 0)).epsilon(1e-12));
    CHECK(mc.X[0] == doctest::Approx(cc.X[0]).epsilon(1e-12));
}

TEST_CASE("a mixed network loads the open class on top of the closed one") {
    // One queue plus one delay, one closed class of population 2 and one open
    // class arriving at the queue.
    const std::size_t M = 2, R = 2;
    Matrix<double> lambda0(M, R, 0.0), Ca0(M, R, 1.0), mu(M, R, 1.0), Cs(M, R, 1.0);
    lambda0(0, 0) = 0.2;  // class 0 is open
    std::vector<Matrix<double>> P(R, Matrix<double>(M, M, 0.0));
    P[1](0, 1) = 1.0;     // class 1 cycles queue -> delay -> queue
    P[1](1, 0) = 1.0;
    std::vector<long> c(M, 1), N(R, 0), ref(R, -1);
    c[1] = 0;
    N[1] = 2;
    ref[1] = 0;
    std::vector<char> ins(M, 0), open(R, 0);
    open[0] = 1;
    const MeResult<double> r = me_mqn(M, R, open, lambda0, Ca0, N, mu, Cs, P, c, ref, ins);
    // the closed class still conserves its population
    CHECK(r.L(0, 1) + r.L(1, 1) == doctest::Approx(2.0).epsilon(1e-6));
    // the open class is present and inflated by the closed occupancy
    CHECK(r.L(0, 0) > 0.0);
    CHECK(r.lambda(0, 0) == doctest::Approx(0.2));
    CHECK(r.X[0] == doctest::Approx(0.2));
    // its response time follows Little's law on the open throughput
    CHECK(r.W(0, 0) == doctest::Approx(r.L(0, 0) / 0.2).epsilon(1e-9));
    // and the closed class is slowed down by the open load: with the server
    // capacity reduced by rho_open the closed throughput drops below the
    // pure closed value
    const MeResult<double> pure = me_cqn(M, 1u, std::vector<long>(1, 2), Matrix<double>(M, 1, 1.0),
                                         Matrix<double>(M, 1, 1.0),
                                         std::vector<Matrix<double>>(1, P[1]), c,
                                         std::vector<long>(1, 0), ins);
    CHECK(r.X[1] < pure.X[0]);
}

TEST_CASE("me rejects malformed input") {
    Matrix<double> lambda0(1, 1, 0.5), Ca0(1, 1, 1.0), mu(1, 1, 1.0), Cs(1, 1, 1.0);
    std::vector<Matrix<double>> P(1, Matrix<double>(1, 1, 0.0));
    std::vector<long> c(1, 1);
    std::vector<char> ins(1, 0);
    // wrong number of routing matrices
    std::vector<Matrix<double>> Pbad;
    CHECK_THROWS_AS(me_oqn(1u, 1u, lambda0, Ca0, mu, Cs, Pbad, c, ins), line::InputError);
    // a multiserver station in the closed algorithm
    std::vector<long> c2(1, 2);
    CHECK_THROWS_AS(me_cqn(1u, 1u, std::vector<long>(1, 1), mu, Cs, P, c2,
                           std::vector<long>(1, 0), ins),
                    line::InputError);
    // a negative population
    CHECK_THROWS_AS(me_cqn(1u, 1u, std::vector<long>(1, -1), mu, Cs, P, c,
                           std::vector<long>(1, 0), ins),
                    line::InputError);
}
