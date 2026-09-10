/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Per-type waiting times of the MMAP[K]/G[K]/1 queue.
 *
 * THE ORACLES ARE INDEPENDENT OF THE IMPLEMENTATION. With one type, Poisson
 * arrivals and exponential service the queue is M/M/1, whose waiting time is an
 * atom plus an exponential tail in closed form; with deterministic service it is
 * M/D/1 and Pollaczek-Khinchine applies. The MULTICLASS numbers are the MATLAB
 * reference, itself validated against BuTools MMAPPH1FCFS on phase-type service
 * (agreement 4e-16) and against JMT at 8e6 samples on general service (5e-4 per
 * type).
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/qsys/qsys_mmapgk1.h"

namespace qsys = line::qsys;
namespace lang = line::lang;
using line::Matrix;
using Dist = line::lang::Distrib<double>;

namespace {

Matrix<double> scalar(double v) {
    Matrix<double> m(1, 1, v);
    return m;
}

}  // namespace

TEST_CASE("qsys_mmapgk1: one type, Poisson and exponential, is M/M/1") {
    const double lam = 0.6, mu = 1.0, rho = lam / mu;
    std::vector<Matrix<double> > MM;
    MM.push_back(scalar(-lam));
    MM.push_back(scalar(lam));
    MM.push_back(scalar(lam));
    std::vector<Dist> svc;
    svc.push_back(Dist::exp_rate(mu));
    std::vector<double> pts;
    pts.push_back(0.5);
    pts.push_back(2.0);
    const qsys::MmapGk1Result<double> r = qsys::qsys_mmapgk1(MM, svc, pts, 3, 1e-12, 10000);
    CHECK(r.meanWaitingTime[0] == doctest::Approx(rho / (mu - lam)).epsilon(1e-9));
    CHECK(r.waitMoments[0][1] == doctest::Approx(2 * rho / std::pow(mu - lam, 2.0)).epsilon(1e-9));
    CHECK(r.idleVector[0] == doctest::Approx(1 - rho).epsilon(1e-12));
    CHECK(r.waitCDF[0][0] == doctest::Approx(1 - rho * std::exp(-(mu - lam) * 0.5)).epsilon(1e-6));
    CHECK(r.waitCDF[0][1] == doctest::Approx(1 - rho * std::exp(-(mu - lam) * 2.0)).epsilon(1e-6));
}

TEST_CASE("qsys_mmapgk1: deterministic service reproduces Pollaczek-Khinchine") {
    const double d = 0.8, lam = 0.9;
    std::vector<Matrix<double> > MM;
    MM.push_back(scalar(-lam));
    MM.push_back(scalar(lam));
    MM.push_back(scalar(lam));
    std::vector<Dist> svc;
    svc.push_back(Dist::det(d));
    const qsys::MmapGk1Result<double> r = qsys::qsys_mmapgk1(MM, svc);
    CHECK(r.meanWaitingTime[0] == doctest::Approx(lam * d * d / (2 * (1 - lam * d))).epsilon(1e-9));
}

TEST_CASE("qsys_mmapgk1: two types with general, type-dependent service") {
    // MMPP2(1.2, 0.3, 0.2, 0.4) marked by phase, Det for type 1, Uniform for type 2
    Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0), D2(2, 2, 0.0), Dsum(2, 2, 0.0);
    D0(0, 0) = -(1.2 + 0.2);
    D0(0, 1) = 0.2;
    D0(1, 0) = 0.4;
    D0(1, 1) = -(0.3 + 0.4);
    D1(0, 0) = 1.2;
    D2(1, 1) = 0.3;
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) Dsum(i, j) = D1(i, j) + D2(i, j);
    std::vector<Matrix<double> > MM;
    MM.push_back(D0);
    MM.push_back(Dsum);
    MM.push_back(D1);
    MM.push_back(D2);
    std::vector<Dist> svc;
    svc.push_back(Dist::det(0.5));
    svc.push_back(Dist::uniform(0.1, 0.9));
    const qsys::MmapGk1Result<double> r = qsys::qsys_mmapgk1(MM, svc);
    CHECK(r.meanWaitingTime[0] == doctest::Approx(0.319282).epsilon(1e-4));
    CHECK(r.meanWaitingTime[1] == doctest::Approx(0.105970).epsilon(1e-4));
    // the two types see DIFFERENT delays: that is the capability, not noise
    CHECK(r.meanWaitingTime[0] > 2.0 * r.meanWaitingTime[1]);
}
