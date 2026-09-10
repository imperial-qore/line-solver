/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * npfqn_sqd, the Smith Queue Decomposition for Blocking-After-Service networks.
 *
 * MATLAB's entry point takes an sn, so the reference values were produced by
 * building the two networks below in MATLAB, calling npfqn_sqd(sn, ...) and
 * dumping BOTH the outputs and the six arrays the routine unpacks from sn, so
 * the port is fed byte-identical inputs rather than inputs reconstructed by
 * hand. Reproduce with:
 *
 *   model = Network('A');
 *   d = Delay(model,'Think');
 *   q1 = Queue(model,'Q1',SchedStrategy.FCFS);   % likewise q2, q3
 *   cl = ClosedClass(model,'C1',6,d);
 *   d.setService(cl,Exp(1/2.0)); q1.setService(cl,Exp(1/1.5));
 *   q2.setService(cl,Exp(1/0.8)); q3.setService(cl,Exp(1/1.2));
 *   q1.setCapacity(4); q2.setCapacity(3); q3.setCapacity(5);
 *   model.link(Network.serialRouting(d,q1,q2,q3));
 *   [X,Q,U,R] = npfqn_sqd(model.getStruct(), 6, 0, true, 'downstream', 'compound');
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/npfqn/npfqn_sqd.h"

using line::Matrix;
using line::Real50;
using line::num_traits;
using line::npfqn::SqdOptions;
using line::npfqn::SqdResult;
using line::npfqn::npfqn_sqd;

namespace {

double relerr(double a, double b) {
    if (b == 0.0) return std::fabs(a);
    return std::fabs(a - b) / std::fabs(b);
}

/** Model A: Delay + three finite-capacity FCFS queues in a cycle, one class. */
template <class T>
struct ModelA {
    std::vector<T> ST, V, cap;
    std::vector<bool> isDelay;
    Matrix<T> rt;
    std::vector<std::size_t> s2sf;
    ModelA() : ST(4), V(4), cap(4), isDelay(4, false), rt(4, 4, num_traits<T>::from_int(0)), s2sf(4) {
        const double st[4] = {2.0, 1.5, 0.8, 1.2};
        const double cp[4] = {0.0, 4.0, 3.0, 5.0};
        for (std::size_t i = 0; i < 4; ++i) {
            ST[i] = num_traits<T>::from_double(st[i]);
            V[i] = num_traits<T>::from_int(1);
            cap[i] = (i == 0) ? T(std::numeric_limits<T>::infinity())
                              : num_traits<T>::from_double(cp[i]);
            s2sf[i] = i + 1;  // one-based, as sn.stationToStateful stores it
        }
        isDelay[0] = true;
        // serial routing d -> q1 -> q2 -> q3 -> d, class count 1 so rt is 4 x 4
        rt(0, 1) = num_traits<T>::from_int(1);
        rt(1, 2) = num_traits<T>::from_int(1);
        rt(2, 3) = num_traits<T>::from_int(1);
        rt(3, 0) = num_traits<T>::from_int(1);
    }
    SqdResult<T> run(int N, const SqdOptions<T>& o) const {
        return npfqn_sqd(ST, V, cap, isDelay, rt, s2sf, 1, N, o);
    }
};

/** Model B: two finite-capacity FCFS queues, no delay station at all. */
template <class T>
struct ModelB {
    std::vector<T> ST, V, cap;
    std::vector<bool> isDelay;
    Matrix<T> rt;
    std::vector<std::size_t> s2sf;
    ModelB() : ST(2), V(2), cap(2), isDelay(2, false), rt(2, 2, num_traits<T>::from_int(0)), s2sf(2) {
        ST[0] = num_traits<T>::from_int(1);
        ST[1] = num_traits<T>::from_int(2);
        V[0] = num_traits<T>::from_int(1);
        V[1] = num_traits<T>::from_int(1);
        cap[0] = num_traits<T>::from_int(3);
        cap[1] = num_traits<T>::from_int(4);
        s2sf[0] = 1;
        s2sf[1] = 2;
        rt(0, 1) = num_traits<T>::from_int(1);
        rt(1, 0) = num_traits<T>::from_int(1);
    }
    SqdResult<T> run(int N, const SqdOptions<T>& o) const {
        return npfqn_sqd(ST, V, cap, isDelay, rt, s2sf, 1, N, o);
    }
};

void check4(const SqdResult<double>& r, const double* X, const double* Q, const double* U,
            const double* R, std::size_t M) {
    for (std::size_t i = 0; i < M; ++i) {
        CHECK(relerr(r.X[i], X[i]) < 1e-13);
        CHECK(relerr(r.Q[i], Q[i]) < 1e-13);
        CHECK(relerr(r.U[i], U[i]) < 1e-13);
        CHECK(relerr(r.R[i], R[i]) < 1e-13);
    }
}

}  // namespace

TEST_CASE("npfqn_sqd reproduces MATLAB across every option switch") {
    const ModelA<double> A;

    SUBCASE("calibration mode 0, downstream, compound, server blocking time") {
        const double X[4] = {0.56111624400402627, 0.56111624400402627, 0.56111624400402627,
                             0.56111624400402627};
        const double Q[4] = {1.1222324880080525, 2.4866424905250337, 0.76665619268716534,
                             1.6244688287797482};
        const double U[4] = {1, 0.8416743660060394, 0.44889299520322101, 0.67333949280483152};
        const double R[4] = {2, 4.431599543047966, 1.3663054685717924, 2.8950664788953997};
        SqdOptions<double> o;
        check4(A.run(6, o), X, Q, U, R, 4);
    }

    SUBCASE("calibration mode 1, the fixed heuristic") {
        const double X[4] = {0.56104367021844215, 0.56104367021844215, 0.56104367021844215,
                             0.56104367021844215};
        const double Q[4] = {1.1220873404368843, 2.4862671157254237, 0.76674483943090865,
                             1.6249007044067825};
        const double U[4] = {1, 0.84156550532766317, 0.44883493617475373, 0.67325240426213051};
        const double R[4] = {2, 4.4315037272542375, 1.3666402102573885, 2.8962107419804379};
        SqdOptions<double> o;
        o.calibrationMode = 1;
        check4(A.run(6, o), X, Q, U, R, 4);
    }

    SUBCASE("calibration mode 2, blocking aware") {
        const double X[4] = {0.56115213695832689, 0.56115213695832689, 0.56115213695832689,
                             0.56115213695832689};
        const double Q[4] = {1.1223042739166538, 2.486849362054032, 0.76654659720834895,
                             1.6242997668209647};
        const double U[4] = {1, 0.84172820543749038, 0.44892170956666155, 0.67338256434999222};
        const double R[4] = {2, 4.4316847397815646, 1.3660227712280375, 2.8945800253480827};
        SqdOptions<double> o;
        o.calibrationMode = 2;
        check4(A.run(6, o), X, Q, U, R, 4);
    }

    SUBCASE("server blocking time switched off") {
        const double X[4] = {0.56891595426884867, 0.56891595426884867, 0.56891595426884867,
                             0.56891595426884867};
        const double Q[4] = {1.1378319085376973, 2.5379889175335744, 0.75966725857739215,
                             1.5645119153513363};
        const double U[4] = {1, 0.853373931403273, 0.45513276341507897, 0.68269914512261842};
        const double R[4] = {2, 4.4610964035897203, 1.3352890754376725, 2.749987768161632};
        SqdOptions<double> o;
        o.serverBlockingTime = false;
        check4(A.run(6, o), X, Q, U, R, 4);
    }

    SUBCASE("ownserver neighbour and fresh V1") {
        const double X[4] = {0.56107864558970677, 0.56107864558970677, 0.56107864558970677,
                             0.56107864558970677};
        const double Q[4] = {1.1221572911794135, 2.486853222542825, 0.76699090047315865,
                             1.6239985858046035};
        const double U[4] = {1, 0.84161796838456016, 0.44886291647176546, 0.67329437470764808};
        const double R[4] = {2, 4.4322720924962029, 1.3669935694434303, 2.8944223747773239};
        SqdOptions<double> o;
        o.calibrationMode = 1;
        o.ownServerNeighbor = true;
        o.freshV1 = true;
        check4(A.run(6, o), X, Q, U, R, 4);
    }

    SUBCASE("a shorter population sweep, mode 2, ownserver, no blocking time") {
        const double X[4] = {0.41608668818081007, 0.41608668818081007, 0.41608668818081007,
                             0.41608668818081007};
        const double Q[4] = {0.83217337636162014, 1.0031998594170748, 0.43220831972232482,
                             0.73241844449898053};
        const double U[4] = {0.83217337636162014, 0.62413003227121511, 0.33286935054464806,
                             0.49930402581697209};
        const double R[4] = {2, 2.4110356998038238, 1.0387458479193381, 1.7602544501032169};
        SqdOptions<double> o;
        o.calibrationMode = 2;
        o.serverBlockingTime = false;
        o.ownServerNeighbor = true;
        check4(A.run(3, o), X, Q, U, R, 4);
    }
}

TEST_CASE("npfqn_sqd on a network with no delay station") {
    const ModelB<double> B;
    SUBCASE("calibration mode 1") {
        const double X[2] = {0.4741605210165688, 0.4741605210165688};
        const double Q[2] = {0.9678481694136446, 3.032151830586356};
        const double U[2] = {0.4741605210165688, 0.94832104203313761};
        const double R[2] = {2.0411825247252602, 6.3947791859297416};
        SqdOptions<double> o;
        o.calibrationMode = 1;
        check4(B.run(4, o), X, Q, U, R, 2);
    }
    SUBCASE("calibration mode 2 with fresh V1") {
        const double X[2] = {0.47467870674634011, 0.47467870674634011};
        const double Q[2] = {0.96424032794984116, 3.0357596720501587};
        const double U[2] = {0.47467870674634011, 0.94935741349268021};
        const double R[2] = {2.0313536593186896, 6.3953988854031634};
        SqdOptions<double> o;
        o.calibrationMode = 2;
        o.freshV1 = true;
        check4(B.run(4, o), X, Q, U, R, 2);
    }
}

TEST_CASE("npfqn_sqd structural invariants") {
    const ModelA<double> A;
    SqdOptions<double> o;

    SUBCASE("the cycle carries one throughput and Little's law holds per station") {
        const SqdResult<double> r = A.run(6, o);
        for (std::size_t i = 1; i < 4; ++i) CHECK(r.X[i] == doctest::Approx(r.X[0]).epsilon(1e-14));
        for (std::size_t i = 0; i < 4; ++i)
            CHECK(r.Q[i] == doctest::Approx(r.X[i] * r.R[i]).epsilon(1e-12));
    }

    SUBCASE("the delay station is a pure delay: R equals its service time") {
        const SqdResult<double> r = A.run(6, o);
        CHECK(r.R[0] == doctest::Approx(2.0).epsilon(1e-14));
    }

    SUBCASE("utilization is capped at one, as the reference caps it") {
        // Station 0 is the delay and its X*ST exceeds one here, so the min(1,.)
        // in the reference is load bearing rather than cosmetic.
        const SqdResult<double> r = A.run(6, o);
        CHECK(r.U[0] == 1.0);
        CHECK(r.X[0] * 2.0 > 1.0);
    }

    SUBCASE("a zero population returns zeros") {
        const SqdResult<double> r = A.run(0, o);
        for (std::size_t i = 0; i < 4; ++i) {
            CHECK(r.X[i] == 0.0);
            CHECK(r.Q[i] == 0.0);
            CHECK(r.U[i] == 0.0);
        }
    }

    SUBCASE("throughput rises with the population") {
        double prev = 0.0;
        for (int n = 1; n <= 6; ++n) {
            const double x = A.run(n, o).X[0];
            CHECK(x > prev);
            prev = x;
        }
    }

    SUBCASE("INITIAL_V1 is load bearing, not inert") {
        // 692.192 arrives from the original contribution with no derivation. It
        // is NOT washed out by the V1 <- V1 (1 - pBlock) deflation on a sweep
        // this short, so the result depends on it: overriding it moves the
        // answer. Pinned so that a future change to the constant cannot pass
        // unnoticed.
        SqdOptions<double> o2;
        o2.initialV1.assign(4, 100.0);
        const double base = A.run(6, o).X[0];
        const double alt = A.run(6, o2).X[0];
        CHECK(relerr(alt, base) > 1e-6);
    }
}

TEST_CASE("npfqn_sqd instantiates at Real50 and tracks the double result") {
    const ModelA<Real50> A;
    SqdOptions<Real50> o;
    o.calibrationMode = 2;
    const SqdResult<Real50> r = A.run(6, o);
    CHECK(relerr(static_cast<double>(r.X[0]), 0.56115213695832689) < 1e-13);
    CHECK(relerr(static_cast<double>(r.Q[1]), 2.486849362054032) < 1e-13);
}

TEST_CASE("npfqn_sqd rejects malformed input") {
    const ModelA<double> A;
    SqdOptions<double> o;
    std::vector<double> shortST(3, 1.0);
    CHECK_THROWS_AS(npfqn_sqd(shortST, A.V, A.cap, A.isDelay, A.rt, A.s2sf, 1, 6, o),
                    line::InputError);
    CHECK_THROWS_AS(npfqn_sqd(A.ST, A.V, A.cap, A.isDelay, A.rt, A.s2sf, 1, -1, o),
                    line::InputError);
    SqdOptions<double> bad;
    bad.initialV1.assign(2, 1.0);
    CHECK_THROWS_AS(npfqn_sqd(A.ST, A.V, A.cap, A.isDelay, A.rt, A.s2sf, 1, 6, bad),
                    line::InputError);
}
