/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * SolverMAM: the api pieces the analyzer needs, the analyzer itself, the
 * dispatch ladder and every refusal it makes, asserted BY NAME.
 */

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#include "doctest.h"
#include "line/api/mam/mmap_assemble.h"
#include "line/api/mam/mmapph1fcfs.h"
#include "line/util/sylvester.h"

using namespace line;

namespace {

template <class T>
Matrix<T> mat(std::size_t r, std::size_t c, const std::vector<double>& v) {
    Matrix<T> m(r, c, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < r; ++i)
        for (std::size_t j = 0; j < c; ++j) m(i, j) = num_traits<T>::from_double(v[i * c + j]);
    return m;
}

double rel(double a, double b) { return std::fabs(a - b) / std::max(1.0, std::fabs(b)); }

}  // namespace

TEST_CASE("sylvester_solve reproduces the defining equation") {
    // A random-looking but well-conditioned pair; the check is the residual,
    // not a reference value, so it is independent of the implementation.
    const Matrix<double> A = mat<double>(3, 3, {-4, 1, 0.5, 0.2, -3, 1, 0, 0.7, -2});
    const Matrix<double> B = mat<double>(2, 2, {-1.5, 0.3, 0.9, -2.5});
    const Matrix<double> C = mat<double>(3, 2, {1, 2, 3, 4, 5, 6});
    const Matrix<double> X = sylvester_solve(A, B, C);
    const Matrix<double> R = [&] {
        Matrix<double> M = matmul(A, X);
        const Matrix<double> XB = matmul(X, B);
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 2; ++j) M(i, j) += XB(i, j) - C(i, j);
        return M;
    }();
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 2; ++j) CHECK(std::fabs(R(i, j)) < 1e-13);

    // lyap(A,B,C) solves A X + X B + C = 0, the sign-flipped problem.
    const Matrix<double> Y = lyap_solve(A, B, C);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 2; ++j) CHECK(std::fabs(Y(i, j) + X(i, j)) < 1e-13);
}

TEST_CASE("sylvester_solve is exact in Rational arithmetic") {
    const Matrix<Rational> A = mat<Rational>(2, 2, {-2, 1, 0, -3});
    const Matrix<Rational> B = mat<Rational>(2, 2, {-1, 0, 0.5, -4});
    const Matrix<Rational> C = mat<Rational>(2, 2, {1, 0, 0, 1});
    const Matrix<Rational> X = sylvester_solve(A, B, C);
    Matrix<Rational> R = matmul(A, X);
    const Matrix<Rational> XB = matmul(X, B);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) R(i, j) += XB(i, j) - C(i, j);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            CHECK(R(i, j) == num_traits<Rational>::from_int(0));
}

TEST_CASE("mmapph1fcfs collapses to M/M/1") {
    // lambda = 2, mu = 3: E[N] = rho/(1-rho) = 2, P(N=n) = (1-rho) rho^n.
    mam::Mmap<double> arv;
    arv.D0 = mat<double>(1, 1, {-2.0});
    arv.D1 = mat<double>(1, 1, {2.0});
    arv.Dc.push_back(arv.D1);
    std::vector<mam::PhService<double>> svc(1);
    svc[0].sigma = {1.0};
    svc[0].S = mat<double>(1, 1, {-3.0});

    const std::vector<double> m = mam::mmapph1fcfs_ncmean(arv, svc);
    REQUIRE(m.size() == 1);
    CHECK(rel(m[0], 2.0) < 1e-10);

    const std::vector<std::vector<double>> d = mam::mmapph1fcfs_ncdistr(arv, svc, 5);
    REQUIRE(d.size() == 1);
    const double rho = 2.0 / 3.0;
    for (std::size_t n = 0; n < 5; ++n) CHECK(rel(d[0][n], (1 - rho) * std::pow(rho, n)) < 1e-9);
}

TEST_CASE("mmapph1fcfs on a two-class marked Poisson stream") {
    // Two Poisson streams, rates 0.4 and 0.6, exponential service of rates 2
    // and 1 respectively. The aggregate is an M/G/1 with rho = 0.8; the
    // reference value comes from MATLAB's MMAPPH1FCFS (see test_mam matlab
    // oracle in the port notes).
    mam::Mmap<double> arv;
    arv.D0 = mat<double>(1, 1, {-1.0});
    arv.D1 = mat<double>(1, 1, {1.0});
    arv.Dc.push_back(mat<double>(1, 1, {0.4}));
    arv.Dc.push_back(mat<double>(1, 1, {0.6}));
    std::vector<mam::PhService<double>> svc(2);
    svc[0].sigma = {1.0};
    svc[0].S = mat<double>(1, 1, {-2.0});
    svc[1].sigma = {1.0};
    svc[1].S = mat<double>(1, 1, {-1.0});

    const std::vector<double> m = mam::mmapph1fcfs_ncmean(arv, svc);
    REQUIRE(m.size() == 2);
    // Independent check: the aggregate number in system of an M/G/1 whose
    // service is the arrival-weighted mixture. rho = 0.4/2 + 0.6/1 = 0.8,
    // E[S] = 0.4*0.5 + 0.6*1 = 0.8, E[S^2] = 0.4*2*0.25 + 0.6*2*1 = 1.4,
    // E[N] = rho + lambda^2 E[S^2] / (2 (1-rho)) = 0.8 + 1*1.4/0.4 = 4.3.
    CHECK(rel(m[0] + m[1], 4.3) < 1e-9);
    // Little's law per class against the FCFS common waiting time:
    // Wq = lambda E[S^2] / (2(1-rho)) = 1.4/0.4 = 3.5, so
    // N_1 = 0.4*(3.5+0.5) = 1.6 and N_2 = 0.6*(3.5+1) = 2.7.
    CHECK(rel(m[0], 1.6) < 1e-9);
    CHECK(rel(m[1], 2.7) < 1e-9);
}

TEST_CASE("mmapph1fcfs with a correlated MMPP arrival and Erlang service") {
    // No closed form here; the check is that the ncDistr is a probability
    // distribution and that its mean agrees with ncMoms, two quantities the
    // implementation computes by different recursions.
    mam::Mmap<double> arv;
    arv.D0 = mat<double>(2, 2, {-2.5, 0.2, 0.1, -0.7});
    arv.Dc.push_back(mat<double>(2, 2, {2.3, 0, 0, 0.6}));
    arv.D1 = arv.Dc[0];
    std::vector<mam::PhService<double>> svc(1);
    svc[0].sigma = {1.0, 0.0};
    svc[0].S = mat<double>(2, 2, {-6, 6, 0, -6});

    const std::vector<double> m = mam::mmapph1fcfs_ncmean(arv, svc);
    const std::vector<std::vector<double>> d = mam::mmapph1fcfs_ncdistr(arv, svc, 200);
    double sum = 0.0, mean = 0.0;
    for (std::size_t n = 0; n < d[0].size(); ++n) {
        sum += d[0][n];
        mean += static_cast<double>(n) * d[0][n];
    }
    CHECK(rel(sum, 1.0) < 1e-8);
    CHECK(rel(mean, m[0]) < 1e-6);
}

TEST_CASE("mmap primitives for exponential marking and per-class scaling") {
    const std::vector<double> lam = {0.4, 0.6};
    const mam::Mmap<double> p = mam::mmap_exponential_vec(lam, 1);
    const std::vector<double> l = mam::mmap_count_lambda(p);
    CHECK(rel(l[0], 0.4) < 1e-12);
    CHECK(rel(l[1], 0.6) < 1e-12);

    // Re-marking a single-type stream by [0.25 0.75] splits its rate.
    mam::Mmap<double> one;
    one.D0 = mat<double>(1, 1, {-1.0});
    one.D1 = mat<double>(1, 1, {1.0});
    one.Dc.push_back(one.D1);
    const mam::Mmap<double> marked = mam::mmap_mark_probs(one, mat<double>(1, 2, {0.25, 0.75}));
    const std::vector<double> lm = mam::mmap_count_lambda(marked);
    CHECK(rel(lm[0], 0.25) < 1e-12);
    CHECK(rel(lm[1], 0.75) < 1e-12);

    // Retargeting the per-class means hits them exactly on a Poisson stream.
    const mam::Mmap<double> sc = mam::mmap_scale_perclass(marked, {4.0, 2.0});
    const std::vector<double> ls = mam::mmap_count_lambda(sc);
    CHECK(rel(ls[0], 0.25) < 1e-12);
    CHECK(rel(ls[1], 0.5) < 1e-12);
}

TEST_CASE("mmap_super_safe superposes and caps the order") {
    mam::Mmap<double> a = mam::mmap_exponential_vec(std::vector<double>{1.0}, 1);
    mam::Mmap<double> b;
    b.D0 = mat<double>(2, 2, {-4, 4, 0, -4});
    b.Dc.push_back(mat<double>(2, 2, {0, 0, 4, 0}));
    b.D1 = b.Dc[0];

    // THE CLASS ORDER IS THE INPUT ORDER. The fold still runs low-SCV first, so
    // the Erlang-2 (SCV 1/2, rate 2) is superposed before the Poisson (SCV 1,
    // rate 1), but the marks are permuted back afterwards: a caller reads mark k
    // as its own k-th class, so letting the sort decide it renamed the classes.
    const mam::Mmap<double> s = mam::mmap_super_safe(std::vector<mam::Mmap<double>>{a, b}, 128);
    CHECK(s.order() == 2);
    CHECK(s.classes() == 2);
    const std::vector<double> ls = mam::mmap_count_lambda(s);
    CHECK(rel(ls[0], 1.0) < 1e-10);
    CHECK(rel(ls[1], 2.0) < 1e-10);

    // maxorder 1 forces the Poisson collapse of every component; the input order
    // survives it, because the permutation happens after the fold.
    const mam::Mmap<double> p = mam::mmap_super_safe(std::vector<mam::Mmap<double>>{a, b}, 1);
    CHECK(p.order() == 1);
    const std::vector<double> lp = mam::mmap_count_lambda(p);
    CHECK(rel(lp[0], 1.0) < 1e-10);
    CHECK(rel(lp[1], 2.0) < 1e-10);
}

TEST_CASE("mmap_super_safe refuses the unported MAP(2) compression by name") {
    mam::Mmap<double> big;
    const std::size_t n = 5;
    big.D0 = Matrix<double>(n, n, 0.0);
    Matrix<double> D1(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        big.D0(i, i) = -5.0;
        D1(i, (i + 1) % n) = 5.0;
    }
    big.D1 = D1;
    big.Dc.push_back(D1);
    try {
        mam::mmap_super_safe(std::vector<mam::Mmap<double>>{big}, 4);
        FAIL("expected a refusal");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("mamap2m_fit_gamma_fb_mmap") != std::string::npos);
    }
}

// ---------------------------------------------------------------------------
// SolverMAM: the dispatch ladder and the dec.source analyzer.
//
// The NUMBERS come from MATLAB `SolverMAM(model, method).getAvgTable` (the
// oracle script is reproduced in _kb/14-cpp-multiprecision.md). The
// `actualmethod` strings pin the ROUTING: the ladder's branches are not
// disjoint, so a model reaching the right numbers through the wrong branch is a
// latent defect that only surfaces on the next model.
// ---------------------------------------------------------------------------

#include "line/lang/qn/network_builder.h"
#include "line/solvers/mam/solver_mam_runner.h"

using lang::Distrib;
using lang::SchedStrategy;

namespace {

using Dd = Distrib<double>;

/** Source -> Queue -> Sink with one open class. */
qn::Network<double> sqs(const std::string& name, SchedStrategy sched, const Dd& arrival,
                        const Dd& service, double servers = 1.0, double cap = -1.0) {
    qn::Network<double> m(name);
    const std::size_t s = m.add_source("Src");
    const std::size_t q = m.add_queue("Q1", sched);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("C1");
    m.set_arrival(s, o, arrival);
    m.set_service(q, o, service);
    if (servers != 1.0) m.set_number_of_servers(q, servers);
    if (cap > 0.0) m.set_capacity(q, cap);
    qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/** Source -> Queue -> Sink with two open classes, no class switching. */
qn::Network<double> sqs2(const std::string& name, SchedStrategy sched, const Dd& a1, const Dd& a2,
                         const Dd& s1, const Dd& s2) {
    qn::Network<double> m(name);
    const std::size_t s = m.add_source("Src");
    const std::size_t q = m.add_queue("Q1", sched);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o1 = m.add_open_class("C1");
    const std::size_t o2 = m.add_open_class("C2");
    m.set_arrival(s, o1, a1);
    m.set_arrival(s, o2, a2);
    m.set_service(q, o1, s1);
    m.set_service(q, o2, s2);
    qn::RoutingMatrix<double> P;
    for (std::size_t r : {o1, o2}) {
        P.set(r, r, s, q, 1.0);
        P.set(r, r, q, k, 1.0);
    }
    m.link(P);
    return m;
}

mva::AvgResult<double> run(qn::Network<double>& m, const std::string& method = "default") {
    mam::MamOptions o;
    o.method = method;
    return mam::solver_mam_run_analyzer(m.get_struct(), o);
}

mam::MamSolution<double> solve(qn::Network<double>& m, const std::string& method = "default") {
    mam::MamOptions o;
    o.method = method;
    return mam::solver_mam_solve(m.get_struct(), o);
}

const double kEps = 1e-9;

}  // namespace

TEST_CASE("dec.source on a single-class open M/M/1") {
    qn::Network<double> m = sqs("A", SchedStrategy::FCFS, Dd::exp_rate(1.0), Dd::exp_rate(2.0));
    const mam::MamSolution<double> d = solve(m);
    CHECK(d.actualmethod == "dec.source");
    const mva::AvgResult<double> r = run(m);
    CHECK(r.actualmethod == "default/dec.source");
    CHECK(r.QN(1, 0) == doctest::Approx(1.0).epsilon(kEps));
    CHECK(r.UN(1, 0) == doctest::Approx(0.5).epsilon(kEps));
    CHECK(r.RN(1, 0) == doctest::Approx(1.0).epsilon(kEps));
    CHECK(r.TN(1, 0) == doctest::Approx(1.0).epsilon(kEps));
    CHECK(r.TN(0, 0) == doctest::Approx(1.0).epsilon(kEps));
    // 'dec.source' asked for by name reports itself, without the default/ prefix.
    const mva::AvgResult<double> r2 = run(m, "dec.source");
    CHECK(r2.actualmethod == "dec.source");
    CHECK(r2.QN(1, 0) == doctest::Approx(1.0).epsilon(kEps));
}

TEST_CASE("dec.source on an M/E3/1 through MMAP[K]/PH[K]/1") {
    qn::Network<double> m =
        sqs("B", SchedStrategy::FCFS, Dd::exp_rate(1.0), Dd::erlang(6.0, 3));
    const mva::AvgResult<double> r = run(m);
    CHECK(r.QN(1, 0) == doctest::Approx(0.833333333333332).epsilon(kEps));
    CHECK(r.UN(1, 0) == doctest::Approx(0.5).epsilon(kEps));
    CHECK(r.RN(1, 0) == doctest::Approx(0.833333333333332).epsilon(kEps));
}

TEST_CASE("a correlated MMPP2 arrival takes the exact MAP/MAP/1 fast path") {
    // MMPP2(0.5, 2.0, 0.2, 0.3): D1 = diag(0.5, 2), D0 = [-0.7 0.2; 0.3 -2.3].
    Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D0(0, 0) = -0.7;
    D0(0, 1) = 0.2;
    D0(1, 0) = 0.3;
    D0(1, 1) = -2.3;
    D1(0, 0) = 0.5;
    D1(1, 1) = 2.0;
    qn::Network<double> m =
        sqs("C", SchedStrategy::FCFS, Dd::map_dist(D0, D1, lang::ProcessType::MMPP2),
            Dd::exp_rate(3.0));
    const mam::MamSolution<double> d = solve(m);
    CHECK(d.actualmethod == "exact.mapmap1");
    const mva::AvgResult<double> r = run(m);
    CHECK(r.actualmethod == "default/exact.mapmap1");
    // 4.6e-9 against MATLAB, and the gap is the ALGORITHM, not an error: the
    // port's qbd_mapmap1 evaluates the closed form pi_0 R (I-R)^-2 e while
    // Q_CT_MAP_MAP_1 sums a level distribution truncated at MaxNumComp.
    CHECK(r.QN(1, 0) == doctest::Approx(0.755144133935373).epsilon(1e-8));
    CHECK(r.UN(1, 0) == doctest::Approx(0.366666666666667).epsilon(kEps));
    CHECK(r.RN(1, 0) == doctest::Approx(0.686494667213975).epsilon(1e-8));
    CHECK(r.TN(1, 0) == doctest::Approx(1.1).epsilon(kEps));
}

TEST_CASE("dec.source on a two-class open M/M/1") {
    qn::Network<double> m = sqs2("D", SchedStrategy::FCFS, Dd::exp_rate(0.4), Dd::exp_rate(0.6),
                                 Dd::exp_rate(2.0), Dd::exp_rate(1.0));
    const mva::AvgResult<double> r = run(m);
    CHECK(r.QN(1, 0) == doctest::Approx(1.6).epsilon(kEps));
    CHECK(r.QN(1, 1) == doctest::Approx(2.7).epsilon(kEps));
    CHECK(r.UN(1, 0) == doctest::Approx(0.2).epsilon(kEps));
    CHECK(r.UN(1, 1) == doctest::Approx(0.6).epsilon(kEps));
    CHECK(r.RN(1, 0) == doctest::Approx(4.00000000000001).epsilon(kEps));
    CHECK(r.RN(1, 1) == doctest::Approx(4.50000000000001).epsilon(kEps));
}

TEST_CASE("an open M/M/2 takes the exact PH/M/c branch") {
    qn::Network<double> m =
        sqs("E", SchedStrategy::FCFS, Dd::exp_rate(1.2), Dd::exp_rate(1.0), 2.0);
    const mva::AvgResult<double> r = run(m);
    // The Erlang-C value, not the single-fast-server surrogate.
    CHECK(r.QN(1, 0) == doctest::Approx(1.875).epsilon(1e-12));
    CHECK(r.UN(1, 0) == doctest::Approx(0.6).epsilon(kEps));
    CHECK(r.RN(1, 0) == doctest::Approx(1.5625).epsilon(1e-12));
}

TEST_CASE("an open M/D/1 takes the exact MAP/D/c branch") {
    qn::Network<double> m = sqs("F", SchedStrategy::FCFS, Dd::exp_rate(1.0), Dd::det(0.5));
    const mva::AvgResult<double> r = run(m);
    // 8.4e-10 against MATLAB: both truncate the MAP/D/c lattice chain, at
    // different points.
    CHECK(r.QN(1, 0) == doctest::Approx(0.749999999372292).epsilon(1e-8));
    CHECK(r.UN(1, 0) == doctest::Approx(0.5).epsilon(kEps));
    CHECK(r.RN(1, 0) == doctest::Approx(0.749999999372292).epsilon(1e-8));
}

TEST_CASE("a finite buffer takes the exact M/M/c/K branch") {
    qn::Network<double> m =
        sqs("G", SchedStrategy::FCFS, Dd::exp_rate(1.0), Dd::exp_rate(1.5), 1.0, 5.0);
    const mva::AvgResult<double> r = run(m);
    CHECK(r.QN(1, 0) == doctest::Approx(1.42255639097744).epsilon(kEps));
    CHECK(r.UN(1, 0) == doctest::Approx(0.634586466165413).epsilon(kEps));
    CHECK(r.RN(1, 0) == doctest::Approx(1.49447077409163).epsilon(kEps));
    // The throughput is the ADMITTED rate, below the offered 1.
    CHECK(r.TN(1, 0) == doctest::Approx(0.95187969924812).epsilon(kEps));
}

TEST_CASE("an Erlang arrival takes the exact PH/M/1 branch") {
    qn::Network<double> m =
        sqs("H", SchedStrategy::FCFS, Dd::erlang(2.0, 2), Dd::exp_rate(2.0));
    const mva::AvgResult<double> r = run(m);
    CHECK(r.QN(1, 0) == doctest::Approx(0.809016994374939).epsilon(1e-12));
    CHECK(r.UN(1, 0) == doctest::Approx(0.5).epsilon(kEps));
    CHECK(r.RN(1, 0) == doctest::Approx(0.809016994374941).epsilon(1e-12));
}

TEST_CASE("a Delay station reports U = QLen = T S") {
    qn::Network<double> m("I");
    const std::size_t s = m.add_source("Src");
    const std::size_t dl = m.add_delay("D");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("C1");
    m.set_arrival(s, o, Dd::exp_rate(1.0));
    m.set_service(dl, o, Dd::exp_rate(2.0));
    m.set_service(q, o, Dd::exp_rate(1.5));
    qn::RoutingMatrix<double> P;
    P.set(s, dl, 1.0);
    P.set(dl, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    const mva::AvgResult<double> r = run(m);
    CHECK(r.QN(1, 0) == doctest::Approx(0.5).epsilon(kEps));
    CHECK(r.UN(1, 0) == doctest::Approx(0.5).epsilon(kEps));
    CHECK(r.RN(1, 0) == doctest::Approx(0.5).epsilon(kEps));
    CHECK(r.QN(2, 0) == doctest::Approx(2.0).epsilon(kEps));
    CHECK(r.UN(2, 0) == doctest::Approx(0.666666666666667).epsilon(kEps));
    CHECK(r.RN(2, 0) == doctest::Approx(2.0).epsilon(kEps));
}

TEST_CASE("a PS station takes the U/(1-Utot) form") {
    qn::Network<double> m = sqs2("J", SchedStrategy::PS, Dd::exp_rate(0.5), Dd::exp_rate(0.3),
                                 Dd::exp_rate(2.0), Dd::exp_rate(1.5));
    const mva::AvgResult<double> r = run(m);
    CHECK(r.QN(1, 0) == doctest::Approx(0.454545454545455).epsilon(kEps));
    CHECK(r.QN(1, 1) == doctest::Approx(0.363636363636364).epsilon(kEps));
    CHECK(r.UN(1, 0) == doctest::Approx(0.25).epsilon(kEps));
    CHECK(r.UN(1, 1) == doctest::Approx(0.2).epsilon(kEps));
    CHECK(r.RN(1, 0) == doctest::Approx(0.909090909090909).epsilon(kEps));
    CHECK(r.RN(1, 1) == doctest::Approx(1.21212121212121).epsilon(kEps));
}

TEST_CASE("a closed cyclic model runs the throughput fixed point") {
    qn::Network<double> m("K");
    const std::size_t dl = m.add_delay("D");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 4.0, dl);
    m.set_service(dl, c, Dd::exp_rate(1.0));
    m.set_service(q1, c, Dd::exp_rate(2.0));
    m.set_service(q2, c, Dd::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(dl, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, dl, 1.0);
    m.link(P);
    // An EXPLICIT dec.source still runs the throughput fixed point. Its answer
    // does not conserve the closed population -- 2 + 1.70 + 0.68 = 4.38 against
    // a declared N = 4 -- which is why `default` no longer routes here.
    const mva::AvgResult<double> rd = run(m, "dec.source");
    CHECK(rd.actualmethod == "dec.source");
    CHECK(rd.QN(0, 0) == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(rd.UN(0, 0) == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(rd.RN(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(rd.TN(0, 0) == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(rd.QN(1, 0) == doctest::Approx(1.7007874015748).epsilon(1e-12));
    CHECK(rd.UN(1, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(rd.RN(1, 0) == doctest::Approx(0.850393700787402).epsilon(1e-12));
    CHECK(rd.QN(2, 0) == doctest::Approx(0.682414698162729).epsilon(1e-12));
    CHECK(rd.UN(2, 0) == doctest::Approx(0.666666666666667).epsilon(1e-12));
    CHECK(rd.RN(2, 0) == doctest::Approx(0.341207349081365).epsilon(1e-12));

    // `default` routed a closed model to mna from 2026-08-14 and to BGCHAIN since
    // 2026-08-16: with no open class the background chain is the exact closed
    // CTMC at chain granularity, and this model -- single class, exponential
    // everywhere -- clears bgchain_closed_exact. The values below are therefore
    // no longer mna's fixed point but the EXACT answer, confirmed by two
    // independent oracles in MATLAB (matlab/scratch/bgchain_cyclic_oracle.m):
    // exact MVA and SolverCTMC both give 1.55480033984707, 1.61087510620221,
    // 0.834324553950722, which bgchain reproduces to 2.8e-16. The mna values this
    // case pinned before (1.39857625922123, 1.76934887863302, 0.832074862145753)
    // summed to N but were 10% out at the Delay and 9.8% out at Q1.
    const mva::AvgResult<double> rm = run(m, "default");
    CHECK(rm.actualmethod == "default/bgchain");
    CHECK(rm.QN(0, 0) == doctest::Approx(1.55480033984707).epsilon(1e-9));
    CHECK(rm.UN(0, 0) == doctest::Approx(1.55480033984707).epsilon(1e-9));
    CHECK(rm.RN(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(rm.TN(0, 0) == doctest::Approx(1.55480033984707).epsilon(1e-9));
    CHECK(rm.QN(1, 0) == doctest::Approx(1.61087510620221).epsilon(1e-9));
    CHECK(rm.UN(1, 0) == doctest::Approx(0.777400169923535).epsilon(1e-9));
    CHECK(rm.RN(1, 0) == doctest::Approx(1.03606557377049).epsilon(1e-9));
    CHECK(rm.QN(2, 0) == doctest::Approx(0.834324553950722).epsilon(1e-9));
    CHECK(rm.UN(2, 0) == doctest::Approx(0.518266779949023).epsilon(1e-9));
    CHECK(rm.RN(2, 0) == doctest::Approx(0.536612021857924).epsilon(1e-9));
    // The population is conserved, as it was under mna.
    CHECK(rm.QN(0, 0) + rm.QN(1, 0) + rm.QN(2, 0) == doctest::Approx(4.0).epsilon(1e-9));
}

TEST_CASE("dec.poisson collapses the arrival stream to marked Poisson") {
    qn::Network<double> m = sqs2("L", SchedStrategy::FCFS, Dd::erlang(0.8, 2), Dd::exp_rate(0.6),
                                 Dd::exp_rate(2.0), Dd::exp_rate(1.0));
    const mva::AvgResult<double> rd = run(m, "default");
    CHECK(rd.QN(1, 0) == doctest::Approx(1.51569143421631).epsilon(1e-8));
    CHECK(rd.QN(1, 1) == doctest::Approx(2.63676857566223).epsilon(1e-8));
    CHECK(rd.RN(1, 0) == doctest::Approx(3.78922858554077).epsilon(1e-8));
    CHECK(rd.RN(1, 1) == doctest::Approx(4.39461429277037).epsilon(1e-8));
    // Poisson-ized, the Erlang arrival's low variability is discarded and the
    // model becomes the two-class M/M/1 of case D.
    const mva::AvgResult<double> rp = run(m, "dec.poisson");
    CHECK(rp.actualmethod == "dec.poisson");
    CHECK(rp.QN(1, 0) == doctest::Approx(1.6).epsilon(kEps));
    CHECK(rp.QN(1, 1) == doctest::Approx(2.7).epsilon(kEps));
}

TEST_CASE("an unlisted MAM method refuses by name") {
    qn::Network<double> m = sqs("R", SchedStrategy::FCFS, Dd::exp_rate(1.0), Dd::exp_rate(2.0));
    // APPROVED CHANGE (2026-07-31): the LAST two rows of this case's table --
    // {"dec.mmap", "solver_mam.m"} and {"dec.source.mmap",
    // "solver_mam_basic_mmap"} -- were dropped when both analyzers landed
    // (solver_mam_decmmap.h and solver_mam_basic_mmap.h), for the reason the
    // three notes below give: a row that pins the ABSENCE of a method asserts,
    // once the method exists, that working code must fail. The table is now
    // empty and the case is what remains of it, the whitelist check; the
    // positive routing of the two is covered by "dec.mmap and dec.source.mmap
    // reproduce M/M/1" below. The case is renamed accordingly -- there is no
    // longer an unported method to refuse, only an unknown one.
    //
    // APPROVED CHANGE (2026-07-25): the {"ldqbd", "solver_mam_ldqbd"} row was
    // dropped when the analyzer landed. It had pinned the ABSENCE of ldqbd, and
    // once the capability existed the row asserted that a working method must
    // fail. The positive routing is covered by "the OPEN LD-QBD regime stays
    // opt-in" below. The rows kept here still name genuinely absent
    // analyzers; if one of them is ported, delete its row rather than relaxing
    // the assertion.
    //
    // APPROVED CHANGE (2026-07-28): the {"mna", "solver_mna_open"} row was
    // dropped for the same reason, when mam_dispatch was wired to the ported
    // solver_mna_open / solver_mna_closed. It had pinned the ABSENCE of the mna
    // analyzers on a model -- lambda = 1, mu = 2 -- that tests/test_mam_mna.cpp
    // solves exactly against a closed form, so once they landed the row
    // asserted that a working method must fail. Coverage RISES: test_mam_mna.cpp
    // carries fourteen cases for the pair, against this one refusal.
    //
    // APPROVED CHANGE (2026-07-31): the {"inap"/"inapplus"/"inapinf",
    // "solver_mam_ag"} rows were dropped. solver_mam_ag.h was already a
    // complete, unit-tested port (test_mam_ag_autocat.cpp) with an entry point
    // signature-compatible with mam_dispatch; mam_dispatch.h simply never
    // called it and threw a false "not ported" refusal instead. The rows had
    // pinned that false refusal, so once the dispatch was wired they asserted
    // that ported, working methods must fail. The positive routing is covered
    // by "inap/inapplus/inapinf solve a single-station open M/M/1" below.
    // An unlisted method is refused by the whitelist, before any analyzer.
    try {
        run(m, "amva");
        FAIL("expected a refusal");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("unsupported by this solver") != std::string::npos);
    }
}

TEST_CASE("dec.mmap and dec.source.mmap reproduce M/M/1") {
    // The wiring regression for the two analyzers ported on 2026-07-31. On a
    // single Source -> FCFS Queue -> Sink model both decompositions are EXACT,
    // and for the same reason: the source's departure process is the Poisson
    // arrival itself, so the arrival stream the traffic step delivers to the
    // queue is the model's own and MMAP[K]/PH[K]/1 answers the M/M/1 exactly.
    // The exact MAP/MAP/1 fast path does NOT claim this model -- both processes
    // are renewal, so it stands down (solver_mam_mapmap1_exact.h) -- which is
    // what lets the method dispatch be observed at all.
    qn::Network<double> m = sqs("R", SchedStrategy::FCFS, Dd::exp_rate(1.0), Dd::exp_rate(2.0));
    for (const char* method : {"dec.mmap", "dec.source.mmap"}) {
        const mva::AvgResult<double> r = run(m, method);
        CHECK(r.actualmethod == std::string(method));
        CHECK(r.QN(1, 0) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(r.UN(1, 0) == doctest::Approx(0.5).epsilon(1e-6));
        CHECK(r.RN(1, 0) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(r.TN(1, 0) == doctest::Approx(1.0).epsilon(1e-6));
        // finish_dispatch pins the Source throughput to sn.rates, as every
        // other branch of the ladder does.
        CHECK(r.TN(0, 0) == doctest::Approx(1.0).epsilon(1e-6));
    }
}

TEST_CASE("dec.mmap refuses what solver_mam.m answers with empty metrics") {
    // Both refusals are the reference's two WARNING-and-return paths, which
    // return `[]` and all-zero matrices respectively. A Delay is the sharper of
    // the two: INF is absent from solver_mam.m's discipline list, so a model
    // that every other MAM method solves gets no answer from this one.
    //
    // THE GATE NOW ANSWERS FIRST, and the message is its, not the analyzer's.
    // `mam_feature_set("dec.mmap")` withdraws SchedStrategy_INF -- the same
    // delta SolverMAM.m carries, for the same reason -- and MATLAB's
    // runAnalyzerChecks gates on the PER-METHOD set exactly as this port does,
    // so neither codebase reaches the analyzer's own "no branch for" sentence
    // any more. What the refusal must still do is name the construct, so a
    // caller knows to remove the Delay rather than the method.
    qn::Network<double> md("DECMMAP-DELAY");
    const std::size_t s = md.add_source("Src");
    const std::size_t d = md.add_delay("Think");
    const std::size_t k = md.add_sink("Sink");
    const std::size_t o = md.add_open_class("C1");
    md.set_arrival(s, o, Dd::exp_rate(1.0));
    md.set_service(d, o, Dd::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(s, d, 1.0);
    P.set(d, k, 1.0);
    md.link(P);
    try {
        run(md, "dec.mmap");
        FAIL("expected a refusal for a Delay under dec.mmap");
    } catch (const UnsupportedError& e) {
        const std::string msg(e.what());
        CHECK_MESSAGE(msg.find("SolverMAM") != std::string::npos, msg);
        CHECK_MESSAGE(msg.find("INF") != std::string::npos, msg);
    }
}

TEST_CASE("the RCAT methods are SolverAG's now, and MAM says so by name") {
    // They used to be solved here. `solver_mam_runner.h::list_valid_methods` no
    // longer advertises them and `mam_dispatch` has no branch for them, matching
    // the reference: SolverMAM.m lists them nowhere and supportsModelMethod
    // returns false, while SolverAG.m owns {default,inap,inapplus,inapinf,exact}.
    //
    // WHAT THIS CASE IS FOR is the redirect, not the refusal. A caller carrying an
    // old options.method must be told WHERE the method went; being told only that
    // it is "unsupported by this solver" sends them looking for a regression in
    // MAM. Until 2026-08-19 that was exactly what happened -- check_method threw
    // the generic message first, so the redirect in check_model_method could
    // never fire and three separate comments promised behaviour the code did not
    // have.
    //
    // The NUMERICS moved with the methods and are not duplicated here: the
    // truncated M/M/1, the exact-vs-inap alias equivalence, the M/PH/1
    // Pollaczek-Khinchine mean and the matrix-exponential refusal are all in
    // tests/test_mam_ag_autocat.cpp, asserted against the solver that owns them.
    const std::vector<std::string> valid = mam::list_valid_methods();
    for (const char* gone : {"inap", "inapplus", "inapinf", "exact"}) {
        CHECK(std::find(valid.begin(), valid.end(), gone) == valid.end());
    }
    // Deliberately NOT asserting anything about 'retrial' here. Whether MAM
    // advertises it is being changed concurrently (it is absent from
    // list_valid_methods at the commit this case was written against, and a
    // working tree was adding it), so an assertion either way would pin someone
    // else's in-flight decision rather than this case's subject. The point below
    // is the REDIRECT, which is independent of the rest of the list.

    qn::Network<double> m = sqs("AG1", SchedStrategy::FCFS, Dd::exp_rate(1.0), Dd::exp_rate(2.0));
    for (const char* gone : {"inap", "inapplus", "inapinf", "exact"}) {
        try {
            run(m, gone);
            FAIL("expected a refusal for the moved method");
        } catch (const UnsupportedError& e) {
            const std::string what = e.what();
            CHECK(what.find("moved to SolverAG") != std::string::npos);
            CHECK(what.find("-s ag") != std::string::npos);
            // The generic wording would leave the caller nowhere to go.
            CHECK(what.find("unsupported by this solver") == std::string::npos);
        }
    }
}


TEST_CASE("load dependence outside the LD-QBD shape is refused, not silently dropped") {
    // BUG A: mam_feature_set used to declare LoadDependence unconditionally,
    // and check_model_method never looked at lldscaling at all, so a
    // load-dependent model on any shape other than the single-class closed
    // Delay+Queue silently fell through to solver_mam_basic / solver_mna,
    // neither of which reads st.lldscaling, and was solved at nominal rates.
    // An open Source -> Queue -> Sink is exactly such a shape: it fails
    // is_closed_delay_queue on the population test alone.
    // Built directly, not through sqs(), to keep the Queue's node handle for
    // set_load_dependence (sqs() does not hand its stations back to the caller).
    qn::Network<double> m("LLD1");
    const std::size_t s = m.add_source("Src");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(s, c, Dd::exp_rate(1.0));
    m.set_service(q1, c, Dd::exp_rate(2.0));
    m.set_load_dependence(q1, std::vector<double>{1.0, 0.8});
    qn::RoutingMatrix<double> P;
    P.set(s, q1, 1.0);
    P.set(q1, k, 1.0);
    m.link(P);
    for (const char* method : {"default", "dec.source", "dec.poisson", "mna"}) {
        try {
            run(m, method);
            FAIL("expected a refusal for method ", method);
        } catch (const UnsupportedError& e) {
            const std::string what = e.what();
            CHECK((what.find("load-dependent") != std::string::npos));
        }
    }
}

TEST_CASE("load dependence on the LD-QBD shape still solves under 'default'") {
    // The positive counterpart of the row above: the shape 'default' actually
    // routes to solver_mam_ldqbd (the same model as "'default' routes a
    // single-class closed Delay+Queue to the LD-QBD" below) must NOT be
    // refused by the new structural check. lldscaling = {1,1,1} is a no-op --
    // solver_mam_ldqbd.h only sets hasLLD when some entry differs from one --
    // so the answer is byte-identical to the undeclared-lld baseline, which
    // is what lets this test pin an exact golden without re-deriving one.
    qn::Network<double> m("LLD2");
    const std::size_t dl = m.add_delay("D");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, dl);
    m.set_service(dl, c, Dd::exp_rate(1.0));
    m.set_service(q1, c, Dd::exp_rate(2.0));
    m.set_load_dependence(q1, std::vector<double>{1.0, 1.0, 1.0});
    qn::RoutingMatrix<double> P;
    P.set(dl, q1, 1.0);
    P.set(q1, dl, 1.0);
    m.link(P);

    const mva::AvgResult<double> r = run(m, "default");
    CHECK(r.actualmethod == "default/ldqbd");
    CHECK(r.QN(0, 0) == doctest::Approx(1.57894736842105).epsilon(1e-12));
    CHECK(r.QN(1, 0) == doctest::Approx(1.42105263157895).epsilon(1e-12));
    CHECK(r.UN(1, 0) == doctest::Approx(0.789473684210526).epsilon(1e-12));
    CHECK(r.RN(1, 0) == doctest::Approx(0.9).epsilon(1e-12));

    const mva::AvgResult<double> rl = run(m, "ldqbd");
    CHECK(rl.actualmethod == "ldqbd");
    CHECK(rl.QN(0, 0) == doctest::Approx(1.57894736842105).epsilon(1e-12));
}

TEST_CASE("FCFSPRPRIO refuses at the feature gate, not deep in the station ladder") {
    // BUG C: mam_feature_set used to declare SchedStrategy_FCFSPRPRIO even
    // though no MAM analyzer serves it (solve_fcfs_station dispatches only
    // FCFS and HOL), so the model passed the gate and only failed later, at
    // solver_mam_basic.h's station-ladder throw, with a less specific message.
    qn::Network<double> m =
        sqs("FP", SchedStrategy::FCFSPRPRIO, Dd::exp_rate(1.0), Dd::exp_rate(2.0));
    try {
        run(m, "default");
        FAIL("expected a refusal for an FCFSPRPRIO station");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("FCFSPRPRIO scheduling discipline") != std::string::npos);
    }
}

TEST_CASE("'default' routes a single-class closed Delay+Queue to the LD-QBD") {
    // APPROVED CHANGE (2026-07-25): this case was
    // "the ldqbd branch of 'default' refuses by name" and asserted an
    // UnsupportedError naming solver_mam_ldqbd. It pinned the ABSENCE of the
    // analyzer; when the analyzer landed the assertion became a claim that a
    // working exact method must fail, so it was replaced by the routing it
    // should have been asserting all along. If you are reading this because a
    // metric moved, the numbers below are MATLAB's, not a tolerance.
    qn::Network<double> m("LD");
    const std::size_t dl = m.add_delay("D");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, dl);
    m.set_service(dl, c, Dd::exp_rate(1.0));
    m.set_service(q1, c, Dd::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(dl, q1, 1.0);
    P.set(q1, dl, 1.0);
    m.link(P);

    const mva::AvgResult<double> r = run(m, "default");
    CHECK(r.actualmethod == "default/ldqbd");
    CHECK(r.QN(0, 0) == doctest::Approx(1.57894736842105).epsilon(1e-12));
    CHECK(r.QN(1, 0) == doctest::Approx(1.42105263157895).epsilon(1e-12));
    CHECK(r.UN(1, 0) == doctest::Approx(0.789473684210526).epsilon(1e-12));
    CHECK(r.RN(1, 0) == doctest::Approx(0.9).epsilon(1e-12));

    // KEPT DELIBERATELY: dec.source asked for BY NAME still runs on this same
    // model and does NOT take the ldqbd branch. That is what proves the branch
    // is gated on `method == "default"` rather than on the topology alone --
    // the part a future reader is most likely to break by "simplifying" the
    // dispatch to route on shape. Its answer differs visibly (Delay QLen 2
    // against the true 1.5789), which is why the reference prefers the LD-QBD.
    const mva::AvgResult<double> rd = run(m, "dec.source");
    CHECK(rd.actualmethod == "dec.source");
    CHECK(rd.QN(0, 0) == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(rd.QN(0, 0) > 0.0);
}

TEST_CASE("a non-Markovian service is converted and solved, not refused") {
    // Pareto(2.5, 0.3): mean 0.5, SCV 1/(alpha(alpha-2)) = 0.8. Against Exp(1)
    // arrivals that is an M/G/1 at rho = 0.5, whose Pollaczek-Khinchine queue
    // length depends on the first TWO moments only:
    //   Lq = rho^2 (1 + scv) / (2 (1 - rho)) = 0.45,  Q = Lq + rho = 0.95.
    // sn_nonmarkov_toph's default fit matches both moments exactly, so the QBD
    // solved on the surrogate must reproduce the P-K value, not merely approach
    // it. Before 2026-08-01 this model was refused by name for want of a fit.
    qn::Network<double> m =
        sqs("NM", SchedStrategy::FCFS, Dd::exp_rate(1.0), Dd::pareto(2.5, 0.3));
    const mva::AvgResult<double> r = run(m, "default");
    CHECK(r.UN(1, 0) == doctest::Approx(0.5).epsilon(1e-6));
    CHECK(r.QN(1, 0) == doctest::Approx(0.95).epsilon(1e-4));
    CHECK(r.TN(1, 0) == doctest::Approx(1.0).epsilon(1e-6));
}


TEST_CASE("ldqbd refuses a multiclass model by name") {
    qn::Network<double> m2 = sqs2("MX", SchedStrategy::FCFS, Dd::exp_rate(0.4), Dd::exp_rate(0.6),
                                  Dd::exp_rate(2.0), Dd::exp_rate(1.0));
    try {
        run(m2, "ldqbd");
        FAIL("expected a refusal");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("single-class model") != std::string::npos);
    }
}

TEST_CASE("SolverMAM refuses exact arithmetic by name") {
    qn::Network<Rational> m("EX");
    const std::size_t s = m.add_source("Src");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("C1");
    m.set_arrival(s, o, Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(1)));
    m.set_service(q, o, Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(2)));
    qn::RoutingMatrix<Rational> P;
    P.set(s, q, num_traits<Rational>::from_int(1));
    P.set(q, k, num_traits<Rational>::from_int(1));
    m.link(P);
    mam::MamOptions o2;
    try {
        mam::solver_mam_run_analyzer(m.get_struct(), o2);
        FAIL("expected a refusal");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("transcendental arithmetic") != std::string::npos);
    }
}

// ---------------------------------------------------------------------------
// The @SolverMAM class surface: getMAMResult, getProb, getProbMarg,
// getCdfRespT / getSjrnT, getPerctRespT, and the getTranAvg refusal.
//
// Reference values from MATLAB (`mam_oracle3.m`). Where an independent closed
// form exists it is asserted INSTEAD of the reference value, because an oracle
// whose derivation does not resemble the implementation is worth more: the
// M/M/1 queue-length distribution is geometric, and its sojourn time is
// Exp(mu - lambda), neither of which the age process ever computes.
// ---------------------------------------------------------------------------

#include "line/api/qsys/qsys_bmapm1.h"

namespace {

qn::Network<double> sqs_ps(const std::string& name, const Dd& arrival, const Dd& service) {
    return sqs(name, SchedStrategy::PS, arrival, service);
}

}  // namespace

TEST_CASE("qsys_bmapm1 on the Bolch batch-arrival fixture") {
    // Example 6.4 of Bolch et al., the fixture the reference's own docstring
    // carries. D0 + D1 + D2 has zero row sums by construction.
    Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0), D2(2, 2, 0.0);
    D0(0, 0) = -2.0;      D0(0, 1) = 0.5;
    D0(1, 0) = 1.0 / 3.0; D0(1, 1) = -3.0;
    D1(0, 0) = 0.25;      D1(0, 1) = 0.5;
    D1(1, 0) = 1.0 / 3.0; D1(1, 1) = 1.0;
    D2(0, 0) = 0.25;      D2(0, 1) = 0.5;
    D2(1, 0) = 1.0;       D2(1, 1) = 1.0 / 3.0;
    const qsys::BmapM1Result<double> r =
        qsys::qsys_bmapm1(std::vector<Matrix<double>>{D0, D1, D2}, 11.0);

    CHECK(r.theta[0] == doctest::Approx(0.526315789473684).epsilon(1e-13));
    CHECK(r.theta[1] == doctest::Approx(0.473684210526316).epsilon(1e-13));
    CHECK(r.lambda == doctest::Approx(3.07894736842105).epsilon(1e-13));
    CHECK(r.rho == doctest::Approx(0.279904306220096).epsilon(1e-13));
    CHECK(r.q == doctest::Approx(14.0).epsilon(1e-14));
    CHECK(r.alpha[0] == doctest::Approx(0.526315789473684).epsilon(1e-13));
    CHECK(r.G(0, 0) == doctest::Approx(0.899832866488352).epsilon(1e-12));
    CHECK(r.G(0, 1) == doctest::Approx(0.100167133510669).epsilon(1e-12));
    CHECK(r.G(1, 0) == doctest::Approx(0.11680663737588).epsilon(1e-12));
    CHECK(r.G(1, 1) == doctest::Approx(0.883193362622851).epsilon(1e-12));
    CHECK(r.drift == doctest::Approx(-0.56578947368421).epsilon(1e-13));
    CHECK(r.pi0 == doctest::Approx(0.720095693779897).epsilon(1e-12));
    CHECK(r.meanQueueLength == doctest::Approx(0.520031857079449).epsilon(1e-11));
    CHECK(r.truncLevel == 50);
    // The decay rate is read at the REFERENCE'S 1-based index; sampling one
    // level higher gives 0.413243742742746, a 2.1e-4 shift.
    CHECK(r.decayRate == doctest::Approx(0.413155865252739).epsilon(1e-9));

    // G is stochastic for a stable queue, an invariant the iteration never uses.
    for (std::size_t i = 0; i < 2; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 2; ++j) s += r.G(i, j);
        CHECK(rel(s, 1.0) < 1e-10);
    }
    CHECK(r.drift < 0.0);
}

TEST_CASE("qsys_bmapm1 collapses to M/M/1") {
    // lambda = 1, mu = 2: pi0 = 1 - rho = 1/2, E[N] = 1, decay = rho = 1/2.
    const qsys::BmapM1Result<double> r = qsys::qsys_bmapm1(
        std::vector<Matrix<double>>{Matrix<double>(1, 1, -1.0), Matrix<double>(1, 1, 1.0)}, 2.0);
    CHECK(rel(r.lambda, 1.0) < 1e-13);
    CHECK(rel(r.rho, 0.5) < 1e-13);
    CHECK(rel(r.pi0, 0.5) < 1e-12);
    CHECK(rel(r.meanQueueLength, 1.0) < 1e-10);
    CHECK(rel(r.decayRate, 0.5) < 1e-10);
    // pi0 = 1 - rho holds EXACTLY for M/M/1, and the level solve never uses it.
    CHECK(rel(r.pi0, 1.0 - num_traits<double>::to_double(r.rho)) < 1e-12);
}

TEST_CASE("qsys_bmapm1 refuses an inconsistent BMAP and a weak uniformization by name") {
    Matrix<double> D0(1, 1, -1.0), D1(1, 1, 2.0);  // row sums are not zero
    try {
        qsys::qsys_bmapm1(std::vector<Matrix<double>>{D0, D1}, 2.0);
        FAIL("expected a refusal");
    } catch (const InputError& e) {
        CHECK(std::string(e.what()).find("zero row sums") != std::string::npos);
    }
    try {
        qsys::qsys_bmapm1(std::vector<Matrix<double>>{Matrix<double>(1, 1, -1.0),
                                                      Matrix<double>(1, 1, 1.0)},
                          2.0, 0.5, static_cast<std::size_t>(0), 10000u, 1e-12, 1e-10);
        FAIL("expected a refusal");
    } catch (const InputError& e) {
        CHECK(std::string(e.what()).find("dominate the total outflow") != std::string::npos);
    }
}

TEST_CASE("getMAMResult exposes the M/G/1-type internals of an M/M/1") {
    qn::Network<double> m = sqs("MR", SchedStrategy::FCFS, Dd::exp_rate(1.0), Dd::exp_rate(2.0));
    const qn::NetworkStruct<double>& L = m.get_struct();
    const qsys::BmapM1Result<double> r = mam::solver_mam_get_mam_result(L);
    CHECK(rel(r.lambda, 1.0) < 1e-13);
    CHECK(rel(r.rho, 0.5) < 1e-13);
    CHECK(rel(r.pi0, 0.5) < 1e-12);
    CHECK(rel(r.meanQueueLength, 1.0) < 1e-10);
    CHECK(rel(r.decayRate, 0.5) < 1e-10);
    CHECK(rel(r.q, 3.0) < 1e-13);  // max(-D0) + mu = 1 + 2
}

TEST_CASE("getMAMResult refuses a multi-queue, multiclass or multiserver model by name") {
    qn::Network<double> two("MR2");
    const std::size_t s = two.add_source("Src");
    const std::size_t qa = two.add_queue("Qa", SchedStrategy::FCFS);
    const std::size_t qb = two.add_queue("Qb", SchedStrategy::FCFS);
    const std::size_t k = two.add_sink("Sink");
    const std::size_t o = two.add_open_class("C1");
    two.set_arrival(s, o, Dd::exp_rate(0.5));
    two.set_service(qa, o, Dd::exp_rate(2.0));
    two.set_service(qb, o, Dd::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(s, qa, 1.0);
    P.set(qa, qb, 1.0);
    P.set(qb, k, 1.0);
    two.link(P);
    try {
        mam::solver_mam_get_mam_result(two.get_struct());
        FAIL("expected a refusal");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("single-queue model") != std::string::npos);
    }

    qn::Network<double> mc = sqs2("MR3", SchedStrategy::FCFS, Dd::exp_rate(0.4), Dd::exp_rate(0.6),
                                  Dd::exp_rate(2.0), Dd::exp_rate(1.0));
    try {
        mam::solver_mam_get_mam_result(mc.get_struct());
        FAIL("expected a refusal");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("single-class model") != std::string::npos);
    }

    qn::Network<double> ms =
        sqs("MR4", SchedStrategy::FCFS, Dd::exp_rate(1.0), Dd::exp_rate(1.0), 2.0);
    try {
        mam::solver_mam_get_mam_result(ms.get_struct());
        FAIL("expected a refusal");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("single-server queue") != std::string::npos);
    }

    // A multi-phase service has no M/G/1-type reduction here.
    qn::Network<double> ph =
        sqs("MR5", SchedStrategy::FCFS, Dd::exp_rate(1.0), Dd::erlang(6.0, 3));
    try {
        mam::solver_mam_get_mam_result(ph.get_struct());
        FAIL("expected a refusal");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("exponential service only") != std::string::npos);
    }
}

TEST_CASE("getProb and getProbMarg on an open M/M/1 are the geometric law") {
    qn::Network<double> m = sqs("PB", SchedStrategy::FCFS, Dd::exp_rate(1.0), Dd::exp_rate(2.0));
    const qn::NetworkStruct<double>& L = m.get_struct();
    mam::MamOptions o;
    const mva::AvgResult<double> avg = mam::solver_mam_run_analyzer(L, o);

    const mam::ProbTable<double> pt = mam::solver_mam_get_prob(L, o, 2, avg);
    CHECK(pt.P.rows() == 100);   // the reference's default cutoff
    CHECK(pt.P.cols() == 1);     // one service phase
    const std::vector<double> pm = mam::solver_mam_get_prob_marg(L, o, 2, 1, avg);
    REQUIRE(pm.size() == 100);
    // P(N = n) = (1-rho) rho^n with rho = 1/2, which the age process never
    // computes and which the joint table must reproduce column by column.
    for (std::size_t n = 0; n < 10; ++n) {
        const double want = 0.5 * std::pow(0.5, static_cast<double>(n));
        CHECK(rel(pm[n], want) < 1e-12);
        CHECK(rel(pt.P(n, 0), want) < 1e-12);
    }
    double mass = 0.0;
    for (double v : pm) mass += v;
    CHECK(rel(mass, 1.0) < 1e-12);
}

TEST_CASE("getProbMarg on a closed model truncates at the class population") {
    qn::Network<double> m("PC");
    const std::size_t dl = m.add_delay("D");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, dl);
    m.set_service(dl, c, Dd::exp_rate(1.0));
    m.set_service(q1, c, Dd::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(dl, q1, 1.0);
    P.set(q1, dl, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& L = m.get_struct();
    mam::MamOptions o;
    o.method = "dec.source";
    const mva::AvgResult<double> avg = mam::solver_mam_run_analyzer(L, o);

    const std::vector<double> pm = mam::solver_mam_get_prob_marg(L, o, 2, 1, avg);
    REQUIRE(pm.size() == 4);  // N + 1
    CHECK(pm[0] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(pm[1] == doctest::Approx(0.333333353498688).epsilon(1e-11));
    CHECK(pm[2] == doctest::Approx(0.333333333333333).epsilon(1e-11));
    CHECK(pm[3] == doctest::Approx(0.333333313167979).epsilon(1e-11));
    double mass = 0.0;
    for (double v : pm) mass += v;
    CHECK(rel(mass, 1.0) < 1e-12);
}

TEST_CASE("getProb refuses a multi-queue model by name") {
    qn::Network<double> m("PD");
    const std::size_t s = m.add_source("Src");
    const std::size_t qa = m.add_queue("Qa", SchedStrategy::FCFS);
    const std::size_t qb = m.add_queue("Qb", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("C1");
    m.set_arrival(s, o, Dd::exp_rate(0.5));
    m.set_service(qa, o, Dd::exp_rate(2.0));
    m.set_service(qb, o, Dd::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(s, qa, 1.0);
    P.set(qa, qb, 1.0);
    P.set(qb, k, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& L = m.get_struct();
    mam::MamOptions o2;
    const mva::AvgResult<double> avg = mam::solver_mam_run_analyzer(L, o2);
    for (int which = 0; which < 2; ++which) {
        try {
            if (which == 0) mam::solver_mam_get_prob(L, o2, 2, avg);
            else mam::solver_mam_get_prob_marg(L, o2, 2, 1, avg);
            FAIL("expected a refusal");
        } catch (const UnsupportedError& e) {
            const std::string msg = e.what();
            CHECK(msg.find("multiple queues") != std::string::npos);
            CHECK(msg.find("SolverCTMC or SolverSSA") != std::string::npos);
        }
    }
}

TEST_CASE("getCdfRespT on an M/M/1 is the exponential sojourn law") {
    qn::Network<double> m = sqs("CD", SchedStrategy::FCFS, Dd::exp_rate(1.0), Dd::exp_rate(2.0));
    const qn::NetworkStruct<double>& L = m.get_struct();
    mam::MamOptions o;
    const std::vector<mam::RespTCdf<double>> rd = mam::solver_mam_get_cdf_respt(L, o);
    REQUIRE(rd.size() == 1);
    REQUIRE(rd[0].X.size() == 200);  // the GLOBAL num_cdf_pts, not the dead 100 fallback
    CHECK(rd[0].X.back() == doctest::Approx(19.0).epsilon(1e-12));
    // The M/M/1 sojourn time is Exp(mu - lambda) = Exp(1); the age process never
    // uses that, so it is a genuinely independent check.
    for (std::size_t i = 0; i < rd[0].X.size(); ++i) {
        const double want = 1.0 - std::exp(-rd[0].X[i]);
        CHECK(std::fabs(rd[0].F[i] - want) < 1e-12);
    }
    // Monotone, starting at zero and reaching one.
    CHECK(rd[0].F.front() == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(rd[0].F.back() > 1.0 - 1e-8);
    for (std::size_t i = 1; i < rd[0].F.size(); ++i) CHECK(rd[0].F[i] >= rd[0].F[i - 1]);

    // getSjrnT is the alias, and must agree point for point.
    const std::vector<mam::RespTCdf<double>> sj = mam::solver_mam_get_sjrn_t(L, o);
    REQUIRE(sj.size() == rd.size());
    for (std::size_t i = 0; i < sj[0].F.size(); ++i) CHECK(sj[0].F[i] == rd[0].F[i]);
}

TEST_CASE("getCdfRespT recovers the mean response time getAvg reports") {
    // The strongest check available here: integrate the CDF and compare with a
    // number produced by a completely different code path.
    // E[R] = int_0^inf (1 - F(t)) dt, by the trapezoid rule on the same grid.
    qn::Network<double> m = sqs("CM", SchedStrategy::FCFS, Dd::exp_rate(1.0), Dd::erlang(5.0, 2));
    const qn::NetworkStruct<double>& L = m.get_struct();
    mam::MamOptions o;
    const mva::AvgResult<double> avg = mam::solver_mam_run_analyzer(L, o);
    const std::vector<mam::RespTCdf<double>> rd = mam::solver_mam_get_cdf_respt(L, o);
    double mean = 0.0;
    for (std::size_t i = 1; i < rd[0].X.size(); ++i) {
        const double dx = rd[0].X[i] - rd[0].X[i - 1];
        mean += 0.5 * dx * ((1.0 - rd[0].F[i - 1]) + (1.0 - rd[0].F[i]));
    }
    CHECK(rel(mean, avg.RN(1, 0)) < 2e-4);  // trapezoid error on 200 points
}

TEST_CASE("getCdfRespT on a PS queue takes the MAP/M/1-PS law") {
    qn::Network<double> m = sqs_ps("CP", Dd::exp_rate(1.0), Dd::exp_rate(2.0));
    const qn::NetworkStruct<double>& L = m.get_struct();
    mam::MamOptions o;
    const std::vector<mam::RespTCdf<double>> rd = mam::solver_mam_get_cdf_respt(L, o);
    REQUIRE(rd[0].X.size() == 200);
    // The grid is 10x the M/M/1-PS mean 1/(mu(1-rho)) = 1.
    CHECK(rd[0].X.back() == doctest::Approx(10.0).epsilon(1e-12));
    CHECK(rd[0].F[0] == doctest::Approx(1.81898940354586e-12).epsilon(1e-6).scale(0.0));
    CHECK(rd[0].F[1] == doctest::Approx(0.066491271546077).epsilon(1e-10));
    CHECK(rd[0].F[24] == doctest::Approx(0.734774471486255).epsilon(1e-10));
    CHECK(rd[0].F[49] == doctest::Approx(0.903240160728594).epsilon(1e-10));
    CHECK(rd[0].F[99] == doctest::Approx(0.980267969534225).epsilon(1e-10));
    for (std::size_t i = 1; i < rd[0].F.size(); ++i) CHECK(rd[0].F[i] >= rd[0].F[i - 1]);
}

TEST_CASE("getPerctRespT inverts the response-time CDF") {
    qn::Network<double> m = sqs("PR", SchedStrategy::FCFS, Dd::exp_rate(1.0), Dd::exp_rate(2.0));
    const qn::NetworkStruct<double>& L = m.get_struct();
    mam::MamOptions o;
    const std::vector<std::vector<double>> pc =
        mam::solver_mam_get_perct_respt(L, o, std::vector<double>{0.5, 0.9, 0.95, 0.99});
    REQUIRE(pc.size() == 1);
    REQUIRE(pc[0].size() == 4);
    // MATLAB's values, to every digit it prints.
    CHECK(pc[0][0] == doctest::Approx(0.694030276551248).epsilon(1e-12));
    CHECK(pc[0][1] == doctest::Approx(2.30306009538054).epsilon(1e-12));
    CHECK(pc[0][2] == doctest::Approx(2.99680609911599).epsilon(1e-12));
    CHECK(pc[0][3] == doctest::Approx(4.60599179266703).epsilon(1e-12));
    // They bracket the exact quantiles -ln(1-p) of Exp(1); the residual is the
    // grid interpolation, not the law.
    const double exact[4] = {0.6931471805599453, 2.302585092994046,
                             2.995732273553991, 4.605170185988091};
    for (std::size_t i = 0; i < 4; ++i) CHECK(rel(pc[0][i], exact[i]) < 2e-3);
    // Percentages are accepted as well as fractions, as the reference does.
    const std::vector<std::vector<double>> pcp =
        mam::solver_mam_get_perct_respt(L, o, std::vector<double>{50.0, 90.0, 95.0, 99.0});
    for (std::size_t i = 0; i < 4; ++i) CHECK(pcp[0][i] == doctest::Approx(pc[0][i]));
}

TEST_CASE("getCdfRespT refuses everything outside the single-open-queue regime") {
    // Three stations: the reference warns and returns NO result at all.
    qn::Network<double> m("CE");
    const std::size_t s = m.add_source("Src");
    const std::size_t dl = m.add_delay("D");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("C1");
    m.set_arrival(s, o, Dd::exp_rate(1.0));
    m.set_service(dl, o, Dd::exp_rate(2.0));
    m.set_service(q1, o, Dd::exp_rate(1.5));
    qn::RoutingMatrix<double> P;
    P.set(s, dl, 1.0);
    P.set(dl, q1, 1.0);
    P.set(q1, k, 1.0);
    m.link(P);
    mam::MamOptions opt;
    try {
        mam::solver_mam_get_cdf_respt(m.get_struct(), opt);
        FAIL("expected a refusal");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("single open queue only") != std::string::npos);
    }

    // A PS queue with phase-type service has no MAP/M/1-PS law.
    qn::Network<double> ps = sqs_ps("CF", Dd::exp_rate(1.0), Dd::erlang(6.0, 3));
    try {
        mam::solver_mam_get_cdf_respt(ps.get_struct(), opt);
        FAIL("expected a refusal");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("exponential (order-1) service") != std::string::npos);
    }
}

TEST_CASE("the class surface refuses exact arithmetic by name") {
    qn::Network<Rational> m("XA");
    const std::size_t s = m.add_source("Src");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("C1");
    m.set_arrival(s, o, Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(1)));
    m.set_service(q, o, Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(2)));
    qn::RoutingMatrix<Rational> P;
    P.set(s, q, num_traits<Rational>::from_int(1));
    P.set(q, k, num_traits<Rational>::from_int(1));
    m.link(P);
    const qn::NetworkStruct<Rational>& L = m.get_struct();
    mam::MamOptions o2;
    mva::AvgResult<Rational> avg;
    for (int which = 0; which < 4; ++which) {
        try {
            if (which == 0) mam::solver_mam_get_mam_result(L);
            else if (which == 1) mam::solver_mam_get_prob(L, o2, 2, avg);
            else if (which == 2) mam::solver_mam_get_prob_marg(L, o2, 2, 1, avg);
            else mam::solver_mam_get_cdf_respt(L, o2);
            FAIL("expected a refusal for path ", which);
        } catch (const UnsupportedError& e) {
            CHECK(std::string(e.what()).find("--arith double or --arith real") !=
                  std::string::npos);
        }
    }
}

// ---------------------------------------------------------------------------
// The level-dependent QBD: the exact analyzer `default` prefers over
// dec.source on a single-class closed Delay+Queue, and the transient path
// behind getTranAvg.
// ---------------------------------------------------------------------------

#include "line/solvers/mam/solver_mam_ldqbd.h"

namespace {

/** Closed Delay + Queue, the shape `default` routes to the LD-QBD. */
qn::Network<double> closed_dq(const std::string& name, double N, const Dd& ds, const Dd& qs,
                              double servers = 1.0) {
    qn::Network<double> m(name);
    const std::size_t dl = m.add_delay("D");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
    if (servers != 1.0) m.set_number_of_servers(q, servers);
    const std::size_t c = m.add_closed_class("C1", N, dl);
    m.set_service(dl, c, ds);
    m.set_service(q, c, qs);
    qn::RoutingMatrix<double> P;
    P.set(dl, q, 1.0);
    P.set(q, dl, 1.0);
    m.link(P);
    return m;
}

qn::Network<double> open_sq_cap(const std::string& name, const Dd& a, const Dd& s, double servers,
                                double cap) {
    qn::Network<double> m(name);
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    if (servers != 1.0) m.set_number_of_servers(q, servers);
    if (cap > 0.0) m.set_capacity(q, cap);
    const std::size_t o = m.add_open_class("C1");
    m.set_arrival(src, o, a);
    m.set_service(q, o, s);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("a single-class closed Delay+Queue is EXACT under the LD-QBD") {
    // MATLAB SolverMAM('default') routes here, and its answer equals
    // SolverCTMC's to every printed digit -- the level-dependent arrival rate
    // (N-n) lambda IS the population constraint, not an approximation of it.
    qn::Network<double> m = closed_dq("LD1", 3.0, Dd::exp_rate(1.0), Dd::exp_rate(2.0));
    const mam::MamSolution<double> d = solve(m, "default");
    CHECK(d.actualmethod == "ldqbd");
    const mva::AvgResult<double> r = run(m, "default");
    CHECK(r.actualmethod == "default/ldqbd");
    CHECK(r.QN(0, 0) == doctest::Approx(1.57894736842105).epsilon(1e-12));
    CHECK(r.UN(0, 0) == doctest::Approx(1.57894736842105).epsilon(1e-12));
    CHECK(r.RN(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(r.TN(0, 0) == doctest::Approx(1.57894736842105).epsilon(1e-12));
    CHECK(r.QN(1, 0) == doctest::Approx(1.42105263157895).epsilon(1e-12));
    CHECK(r.UN(1, 0) == doctest::Approx(0.789473684210526).epsilon(1e-12));
    CHECK(r.RN(1, 0) == doctest::Approx(0.9).epsilon(1e-12));
    // Population is conserved: an invariant the LD-QBD never imposes, it falls
    // out of the level distribution summing to one.
    CHECK(rel(r.QN(0, 0) + r.QN(1, 0), 3.0) < 1e-12);
    // Little's law at the queue, likewise never used in the construction.
    CHECK(rel(r.QN(1, 0), r.TN(1, 0) * r.RN(1, 0)) < 1e-12);

    // Asked for by name it reports itself, without the default/ prefix.
    const mva::AvgResult<double> rn = run(m, "ldqbd");
    CHECK(rn.actualmethod == "ldqbd");
    CHECK(rn.QN(1, 0) == doctest::Approx(r.QN(1, 0)));

    // dec.source on the SAME model is visibly different, which is why the
    // reference prefers the LD-QBD here: it reports a Delay queue length of 2
    // against the true 1.579, and a response time of 0.711 against 0.9.
    const mva::AvgResult<double> rd = run(m, "dec.source");
    CHECK(rd.actualmethod == "dec.source");
    CHECK(rd.QN(0, 0) == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(rd.RN(1, 0) == doctest::Approx(0.710526315789474).epsilon(1e-12));
}

TEST_CASE("the LD-QBD closed branch on a multiserver and on PH service") {
    // c = 2, N = 5: MATLAB and SolverCTMC agree to 1e-14 here.
    qn::Network<double> m2 = closed_dq("LD2", 5.0, Dd::exp_rate(1.0), Dd::exp_rate(1.5), 2.0);
    const mva::AvgResult<double> r2 = run(m2, "default");
    CHECK(r2.actualmethod == "default/ldqbd");
    CHECK(r2.QN(0, 0) == doctest::Approx(2.53414809489576).epsilon(1e-12));
    CHECK(r2.QN(1, 0) == doctest::Approx(2.46585190510424).epsilon(1e-12));
    CHECK(r2.UN(1, 0) == doctest::Approx(0.844716031631919).epsilon(1e-11));
    CHECK(rel(r2.QN(0, 0) + r2.QN(1, 0), 5.0) < 1e-12);

    // Erlang-2 service at a SINGLE server, the regime the reference documents
    // as exact for PH; equals SolverCTMC to every printed digit.
    qn::Network<double> m3 = closed_dq("LD3", 4.0, Dd::exp_rate(1.0), Dd::erlang(4.0, 2));
    const mva::AvgResult<double> r3 = run(m3, "default");
    CHECK(r3.actualmethod == "default/ldqbd");
    CHECK(r3.QN(0, 0) == doctest::Approx(1.85178752352005).epsilon(1e-12));
    CHECK(r3.QN(1, 0) == doctest::Approx(2.14821247647995).epsilon(1e-12));
    CHECK(r3.UN(1, 0) == doctest::Approx(0.925893761760023).epsilon(1e-12));
    CHECK(r3.RN(1, 0) == doctest::Approx(1.16007503517274).epsilon(1e-12));
    CHECK(rel(r3.QN(0, 0) + r3.QN(1, 0), 4.0) < 1e-12);
}

TEST_CASE("the OPEN LD-QBD regime stays opt-in, as the reference gates it") {
    // The open Source+Queue regime truncates at options.cutoff, so the
    // reference deliberately does NOT let `default` claim it; `default` keeps
    // dec.source and only `ldqbd` by name reaches the truncated answer.
    qn::Network<double> m = open_sq_cap("LD4", Dd::exp_rate(1.0), Dd::exp_rate(2.0), 1.0, -1.0);
    const mva::AvgResult<double> rd = run(m, "default");
    CHECK(rd.actualmethod == "default/dec.source");
    CHECK(rd.QN(1, 0) == doctest::Approx(1.0).epsilon(1e-12));

    const mva::AvgResult<double> rl = run(m, "ldqbd");
    CHECK(rl.actualmethod == "ldqbd");
    // The truncation deficit against the exact M/M/1 value of 1.
    CHECK(rl.QN(1, 0) == doctest::Approx(0.999999999476131).epsilon(1e-11));
    CHECK(rl.QN(1, 0) < 1.0);

    // Multiserver open: the Erlang-C queue length 1.875, less the truncation.
    qn::Network<double> m5 = open_sq_cap("LD5", Dd::exp_rate(1.2), Dd::exp_rate(1.0), 2.0, -1.0);
    const mva::AvgResult<double> r5 = run(m5, "ldqbd");
    CHECK(r5.QN(1, 0) == doctest::Approx(1.87499999918119).epsilon(1e-11));
    CHECK(rel(r5.QN(1, 0), 1.875) < 1e-8);
}

TEST_CASE("ldqbd still refuses a multiclass model by name") {
    qn::Network<double> m = sqs2("LD6", SchedStrategy::FCFS, Dd::exp_rate(0.4), Dd::exp_rate(0.6),
                                 Dd::exp_rate(2.0), Dd::exp_rate(1.0));
    try {
        run(m, "ldqbd");
        FAIL("expected a refusal");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("single-class model") != std::string::npos);
    }
}

TEST_CASE("the LD-QBD refuses a shape outside its two regimes by name") {
    // Three stations: neither the closed Delay+Queue nor the open Source+Queue.
    qn::Network<double> m("LD7");
    const std::size_t dl = m.add_delay("D");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, dl);
    m.set_service(dl, c, Dd::exp_rate(1.0));
    m.set_service(q1, c, Dd::exp_rate(2.0));
    m.set_service(q2, c, Dd::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(dl, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, dl, 1.0);
    m.link(P);
    try {
        run(m, "ldqbd");
        FAIL("expected a refusal");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("one Delay and one Queue") != std::string::npos);
    }
    // `default` on that model is unaffected by the ldqbd refusal: it is a closed
    // model whose background chain fits and whose service is exponential, so the
    // chooser takes bgchain since 2026-08-16, as MATLAB does (verified: "MAM
    // analysis [method: default/bgchain]" on the same three stations at N = 3,
    // where it reproduces exact MVA to 1.7e-16 against mna's fixed point).
    CHECK(run(m, "default").actualmethod == "default/bgchain");
}

TEST_CASE("getTranAvg on a finite buffer approaches the steady state") {
    // M/M/1/5, lambda = 1, mu = 2. The exact stationary law is
    // p_n = rho^n (1-rho)/(1-rho^6), so E[N] = 0.904761904761905 and
    // U = 1 - p_0 = 0.492063492063492 -- neither of which the transient
    // construction uses, so they are a genuine oracle for its limit.
    qn::Network<double> m = open_sq_cap("TR1", Dd::exp_rate(1.0), Dd::exp_rate(2.0), 1.0, 5.0);
    mam::MamOptions o;
    o.timespan_start = 0.0;
    o.timespan_end = 10.0;
    const mam::TranResult<double> r = mam::solver_mam_get_tran_avg(m.get_struct(), o);
    const mam::TranCurve<double>& Q = r.Qt[1][0];
    const mam::TranCurve<double>& U = r.Ut[1][0];
    const mam::TranCurve<double>& Tp = r.Tt[1][0];
    // The reference's grid: min(101, max(11, round(10 T))) points.
    REQUIRE(Q.times.size() == 100);
    CHECK(Q.times.front() == doctest::Approx(0.0));
    CHECK(Q.times.back() == doctest::Approx(10.0));
    // An empty queue at t = 0, by construction of pi(0).
    CHECK(Q.values.front() == doctest::Approx(0.0).epsilon(1e-14));
    CHECK(U.values.front() == doctest::Approx(0.0).epsilon(1e-14));
    // Monotone approach, and within 1% of the stationary values at t = 10.
    for (std::size_t i = 1; i < Q.values.size(); ++i) CHECK(Q.values[i] >= Q.values[i - 1] - 1e-12);
    CHECK(rel(Q.values.back(), 0.904761904761905) < 3e-3);
    CHECK(rel(U.values.back(), 0.492063492063492) < 3e-3);
    // Throughput is mu times the utilization at every point, an identity the
    // metric loop computes by two different sums.
    for (std::size_t i = 0; i < Q.values.size(); ++i)
        CHECK(std::fabs(Tp.values[i] - 2.0 * U.values[i]) < 1e-12);
}

TEST_CASE("getTranAvg on an infinite buffer takes the adaptive Taylor series") {
    // M/M/1, lambda = 1, mu = 2: E[N] -> rho/(1-rho) = 1 and U -> 1/2.
    qn::Network<double> m = open_sq_cap("TR2", Dd::exp_rate(1.0), Dd::exp_rate(2.0), 1.0, -1.0);
    mam::MamOptions o;
    o.timespan_start = 0.0;
    o.timespan_end = 10.0;
    o.tol = 1e-6;
    const mam::TranResult<double> r = mam::solver_mam_get_tran_avg(m.get_struct(), o);
    const mam::TranCurve<double>& Q = r.Qt[1][0];
    const mam::TranCurve<double>& U = r.Ut[1][0];
    // libQBD's own grid, one point per 1/|min diagonal| = 1/3.
    REQUIRE(Q.times.size() >= 30);
    CHECK(Q.times[1] - Q.times[0] == doctest::Approx(1.0 / 3.0).epsilon(1e-12));
    CHECK(Q.values.front() == doctest::Approx(0.0).epsilon(1e-14));
    for (std::size_t i = 1; i < Q.values.size(); ++i) CHECK(Q.values[i] >= Q.values[i - 1] - 1e-12);
    CHECK(rel(Q.values.back(), 1.0) < 3e-2);
    CHECK(rel(U.values.back(), 0.5) < 1e-2);
}

TEST_CASE("getTranAvg with PH service, finite and infinite") {
    // M/E2/1, service mean 0.4, SCV 1/2: the M/G/1 stationary queue length is
    // rho + rho^2 (1 + cs2) / (2 (1 - rho)) = 0.6, computed here from the
    // Pollaczek-Khinchine formula the transient path never touches.
    qn::Network<double> mi = open_sq_cap("TR3", Dd::exp_rate(1.0), Dd::erlang(5.0, 2), 1.0, -1.0);
    mam::MamOptions o;
    o.timespan_start = 0.0;
    o.timespan_end = 8.0;
    o.tol = 1e-6;
    const mam::TranResult<double> ri = mam::solver_mam_get_tran_avg(mi.get_struct(), o);
    CHECK(rel(ri.Qt[1][0].values.back(), 0.6) < 5e-3);
    CHECK(rel(ri.Ut[1][0].values.back(), 0.4) < 5e-3);

    // The same service behind a finite buffer takes the expm branch instead.
    qn::Network<double> mf = open_sq_cap("TR4", Dd::exp_rate(1.0), Dd::erlang(5.0, 2), 1.0, 4.0);
    const mam::TranResult<double> rf = mam::solver_mam_get_tran_avg(mf.get_struct(), o);
    REQUIRE(rf.Qt[1][0].times.size() == 80);
    CHECK(rf.Qt[1][0].values.front() == doctest::Approx(0.0).epsilon(1e-14));
    // A finite buffer holds strictly less than the infinite one.
    CHECK(rf.Qt[1][0].values.back() < ri.Qt[1][0].values.back());
}

TEST_CASE("getTranAvg on a multiserver finite buffer conserves T = c mu U") {
    qn::Network<double> m = open_sq_cap("TR5", Dd::exp_rate(1.2), Dd::exp_rate(1.0), 2.0, 6.0);
    mam::MamOptions o;
    o.timespan_start = 0.0;
    o.timespan_end = 12.0;
    const mam::TranResult<double> r = mam::solver_mam_get_tran_avg(m.get_struct(), o);
    const mam::TranCurve<double>& U = r.Ut[1][0];
    const mam::TranCurve<double>& Tp = r.Tt[1][0];
    REQUIRE(U.values.size() == 101);
    for (std::size_t i = 0; i < U.values.size(); ++i)
        CHECK(std::fabs(Tp.values[i] - 2.0 * 1.0 * U.values[i]) < 1e-12);
}

TEST_CASE("getTranAvg refuses PH multiserver by name") {
    // HISTORY, twice over. This case began as "getTranAvg refuses by name,
    // naming the unported transient stack", pinning the ABSENCE of the whole
    // transient stack on an ordinary open M/M/1; that half died on 2026-07-25
    // when solver_mam_ldqbd_transient landed. Its replacement half asserted
    // that a correlated MMPP2 arrival refused because the Laplace route was
    // unported, and THAT half died the same day when solver_mam_transient_qbd
    // landed -- the model now runs, and its successor is "a correlated MMPP2
    // arrival now RUNS the Laplace transient". What remains is the one
    // construct that is genuinely still absent.
    //
    // A refusal test pins an absence, so it has a shelf life: when the absent
    // thing arrives, the test must be replaced by one that measures it, not
    // relaxed.
    mam::MamOptions o;
    o.timespan_start = 0.0;
    o.timespan_end = 10.0;

    // PH service with c > 1 has no level representation here.
    qn::Network<double> mp = open_sq_cap("TR7", Dd::exp_rate(1.0), Dd::erlang(6.0, 3), 2.0, 5.0);
    CHECK(!mam::mam_transient_qbd_applicable(mp.get_struct()));
    try {
        mam::solver_mam_get_tran_avg(mp.get_struct(), o);
        FAIL("expected a refusal");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("single-server") != std::string::npos);
    }
}

TEST_CASE("the LD-QBD and its transient refuse exact arithmetic by name") {
    qn::Network<Rational> m("LDX");
    const std::size_t dl = m.add_delay("D");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, dl);
    m.set_service(dl, c, Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(1)));
    m.set_service(q, c, Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(2)));
    qn::RoutingMatrix<Rational> P;
    P.set(dl, q, num_traits<Rational>::from_int(1));
    P.set(q, dl, num_traits<Rational>::from_int(1));
    m.link(P);
    mam::MamOptions o;
    o.method = "ldqbd";
    try {
        mam::solver_mam_solve(m.get_struct(), o);
        FAIL("expected a refusal");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("--arith double or --arith real") != std::string::npos);
    }
}

TEST_CASE("gammainc_lower reproduces the regularized incomplete gamma") {
    // P(1, x) = 1 - exp(-x) exactly, and P(a, x) -> 1 as x grows; both are
    // independent of the continued-fraction/series split the routine uses.
    for (double x : {0.25, 1.0, 2.0, 5.0, 20.0})
        CHECK(rel(mam::gammainc_lower(1.0, x), 1.0 - std::exp(-x)) < 1e-12);
    // P(2, x) = 1 - (1+x) exp(-x).
    for (double x : {0.5, 2.0, 7.0})
        CHECK(rel(mam::gammainc_lower(2.0, x), 1.0 - (1.0 + x) * std::exp(-x)) < 1e-12);
    CHECK(mam::gammainc_lower(5.0, 0.0) == doctest::Approx(0.0));
    CHECK(mam::gammainc_lower(5.0, 200.0) == doctest::Approx(1.0).epsilon(1e-12));
    // Monotone in x at fixed a.
    double prev = 0.0;
    for (int i = 1; i <= 50; ++i) {
        const double v = mam::gammainc_lower(3.0, 0.2 * i);
        CHECK(v >= prev - 1e-15);
        prev = v;
    }
}

TEST_CASE("the transient curves match MATLAB point for point") {
    // MATLAB `solver_mam_ldqbd_transient` on the worktree tree (which() asserted),
    // sampled at the same grid indices. Both engines are pinned: the expm branch
    // on a finite buffer and the libQBD adaptive Taylor series on an infinite
    // one. The grid SIZES are pinned too, because a transient read off a
    // different grid cannot be compared point for point at all.
    mam::MamOptions o;
    o.tol = 1e-6;

    SUBCASE("M/M/1/5, finite buffer, expm branch") {
        qn::Network<double> m = open_sq_cap("MT1", Dd::exp_rate(1.0), Dd::exp_rate(2.0), 1.0, 5.0);
        o.timespan_start = 0.0;
        o.timespan_end = 10.0;
        const mam::TranResult<double> r = mam::solver_mam_get_tran_avg(m.get_struct(), o);
        const mam::TranCurve<double>& Q = r.Qt[1][0];
        const mam::TranCurve<double>& U = r.Ut[1][0];
        const mam::TranCurve<double>& Tp = r.Tt[1][0];
        REQUIRE(Q.times.size() == 100);
        CHECK(Q.values[1] == doctest::Approx(0.091749584145099).epsilon(1e-12));
        CHECK(U.values[1] == doctest::Approx(0.0874158172647847).epsilon(1e-12));
        CHECK(Tp.values[1] == doctest::Approx(0.174831634529569).epsilon(1e-12));
        CHECK(Q.values[2] == doctest::Approx(0.168145589402765).epsilon(1e-12));
        CHECK(Q.values[99] == doctest::Approx(0.902246321661568).epsilon(1e-12));
        CHECK(U.values[99] == doctest::Approx(0.491447986334659).epsilon(1e-12));
        CHECK(Tp.values[99] == doctest::Approx(0.982895972669318).epsilon(1e-12));
    }
    SUBCASE("M/E2/1/4, finite buffer, PH levels") {
        qn::Network<double> m = open_sq_cap("MT2", Dd::exp_rate(1.0), Dd::erlang(5.0, 2), 1.0, 4.0);
        o.timespan_start = 0.0;
        o.timespan_end = 8.0;
        const mam::TranResult<double> r = mam::solver_mam_get_tran_avg(m.get_struct(), o);
        REQUIRE(r.Qt[1][0].times.size() == 80);
        CHECK(r.Qt[1][0].values[1] == doctest::Approx(0.0979711522342536).epsilon(1e-12));
        CHECK(r.Qt[1][0].values[2] == doctest::Approx(0.182184625995815).epsilon(1e-12));
        CHECK(r.Qt[1][0].values[79] == doctest::Approx(0.574680871281024).epsilon(1e-12));
        CHECK(r.Ut[1][0].values[79] == doctest::Approx(0.396644723264007).epsilon(1e-12));
        CHECK(r.Tt[1][0].values[79] == doctest::Approx(0.991607633727084).epsilon(1e-12));
    }
    SUBCASE("M/M/2/6, finite buffer, multiserver") {
        qn::Network<double> m = open_sq_cap("MT3", Dd::exp_rate(1.2), Dd::exp_rate(1.0), 2.0, 6.0);
        o.timespan_start = 0.0;
        o.timespan_end = 12.0;
        const mam::TranResult<double> r = mam::solver_mam_get_tran_avg(m.get_struct(), o);
        REQUIRE(r.Qt[1][0].times.size() == 101);
        CHECK(r.Qt[1][0].values[1] == doctest::Approx(0.13570773117266).epsilon(1e-12));
        CHECK(r.Qt[1][0].values[100] == doctest::Approx(1.62884432751765).epsilon(1e-12));
        CHECK(r.Ut[1][0].values[100] == doctest::Approx(0.58439040345618).epsilon(1e-12));
    }
    SUBCASE("M/M/1, infinite buffer, libQBD adaptive Taylor") {
        qn::Network<double> m = open_sq_cap("MT4", Dd::exp_rate(1.0), Dd::exp_rate(2.0), 1.0, -1.0);
        o.timespan_start = 0.0;
        o.timespan_end = 10.0;
        const mam::TranResult<double> r = mam::solver_mam_get_tran_avg(m.get_struct(), o);
        // libQBD's own grid: 31 points at 1/3 apart, which the port reproduces
        // exactly because it takes the same uniformization step.
        REQUIRE(r.Qt[1][0].times.size() == 31);
        CHECK(r.Qt[1][0].values[1] == doctest::Approx(0.250414662438599).epsilon(1e-11));
        CHECK(r.Qt[1][0].values[2] == doctest::Approx(0.403092662228143).epsilon(1e-11));
        CHECK(r.Qt[1][0].values[15] == doctest::Approx(0.895942508744442).epsilon(1e-11));
        CHECK(r.Qt[1][0].values[30] == doctest::Approx(0.974845769761501).epsilon(1e-11));
        CHECK(r.Ut[1][0].values[30] == doctest::Approx(0.49671100637567).epsilon(1e-11));
    }
    SUBCASE("M/E2/1, infinite buffer, PH levels under libQBD") {
        qn::Network<double> m = open_sq_cap("MT5", Dd::exp_rate(1.0), Dd::erlang(5.0, 2), 1.0, -1.0);
        o.timespan_start = 0.0;
        o.timespan_end = 8.0;
        const mam::TranResult<double> r = mam::solver_mam_get_tran_avg(m.get_struct(), o);
        REQUIRE(r.Qt[1][0].times.size() == 49);
        // The horizon carries MATLAB's own accumulation artifact, 8.00000000000001
        // rather than 8, because the grid is built by repeated addition of
        // 1/min_elem. Reproducing it bit for bit is evidence the step sequence
        // is identical, not a rounding coincidence.
        CHECK(r.Qt[1][0].times.back() == doctest::Approx(8.00000000000001).epsilon(1e-14));
        CHECK(r.Qt[1][0].values[1] == doctest::Approx(0.154261387093885).epsilon(1e-11));
        CHECK(r.Qt[1][0].values[2] == doctest::Approx(0.26690389098165).epsilon(1e-11));
        CHECK(r.Qt[1][0].values[48] == doctest::Approx(0.599174135787991).epsilon(1e-11));
        CHECK(r.Ut[1][0].values[48] == doctest::Approx(0.399810621634186).epsilon(1e-11));
    }
}

// ---------------------------------------------------------------------------
// The vendored ILT-CME table and the numerical inverse Laplace transform.
// ---------------------------------------------------------------------------

#include "line/api/mam/matlab_ilt.h"

TEST_CASE("the vendored ILT-CME table carries every reachable entry") {
    REQUIRE(mam::iltcme::kTableSize > 0);
    // Every entry either consumer can select, and no others: matlab_ilt's scan
    // over the legal maxFnEvals range (127 entries; it seeds at entry 0 before
    // testing the bound, so entry 0 is reachable whatever the budget), UNION the
    // smallest-cv2 entry at each distinct n, which is what cme_table_entry
    // selects. The union was widened from 127 to 171 on 2026-08-01 when
    // api/mam/cme.h was added: dist_fit_me walks the orders upwards, and the
    // low-n tail matlab_ilt can never reach is exactly where it starts.
    CHECK(mam::iltcme::kTableSize == 171);
    for (std::size_t i = 0; i < mam::iltcme::kTableSize; ++i) {
        const mam::iltcme::CmeEntry& e = mam::iltcme::kTable[i];
        CHECK(e.n >= 0);
        CHECK(e.len == static_cast<std::size_t>(e.n));
        CHECK(e.mu1 > 0.0);
        CHECK(e.cv2 > 0.0);
        CHECK(e.a != nullptr);
        CHECK(e.b != nullptr);
    }
    // The budget genuinely selects different entries, which is the property the
    // whole table exists for; a table that always picked one entry would not
    // need vendoring.
    auto pick = [](std::size_t budget) {
        const mam::iltcme::CmeEntry* best = &mam::iltcme::kTable[0];
        for (std::size_t i = 1; i < mam::iltcme::kTableSize; ++i) {
            const mam::iltcme::CmeEntry& e = mam::iltcme::kTable[i];
            if (e.cv2 < best->cv2 && static_cast<std::size_t>(e.n) + 1 <= budget) best = &e;
        }
        return best->n;
    };
    // Pinned against the same scan run over the FULL json table in python.
    CHECK(pick(50) == 49);
    CHECK(pick(100) == 74);
    CHECK(pick(200) == 190);
    CHECK(pick(500) == 480);
    CHECK(pick(1000) == 980);
    // Steeper budgets never select a coarser entry.
    for (std::size_t b = 12; b <= 1000; b += 37) CHECK(pick(b) <= static_cast<int>(b));
}

TEST_CASE("matlab_ilt reproduces MATLAB's inverse transform") {
    // TWO ORACLES, answering different questions. The MATLAB values pin the
    // PORT: they are `matlab_ilt(@(s) 1./(s+a), ts, 100)` run against the
    // now-vendored table, and the port must reproduce them. The analytic
    // comparison pins the METHOD: 1/(s+a) inverts to exp(-a t), and the gap
    // between the two is the CME quadrature's own truncation at this budget,
    // about 2e-5, NOT port error. MATLAB shows the identical gap digit for
    // digit, which is what distinguishes the two. Do not "tighten" the analytic
    // bound; raise maxFnEvals, which is what it is for.
    const std::vector<double> ts = {0.1, 0.5, 1.0, 2.0, 5.0};

    SUBCASE("1/(s+a), a = 1.7, maxFnEvals = 100") {
        const double a = 1.7;
        auto F = [a](const std::complex<double>& s) { return 1.0 / (s + a); };
        const std::vector<double> got = mam::matlab_ilt(F, ts, 100);
        const double matlab[5] = {0.843665658592, 0.427425613266314, 0.182701879777171,
                                  0.0333870251050902, 0.000204497341050402};
        for (std::size_t i = 0; i < ts.size(); ++i)
            CHECK(got[i] == doctest::Approx(matlab[i]).epsilon(1e-10));
        const double abserr[5] = {8.41996e-07, 1.06813e-05, 1.83557e-05, 1.37551e-05, 1.02897e-06};
        for (std::size_t i = 0; i < ts.size(); ++i)
            CHECK(std::fabs(got[i] - std::exp(-a * ts[i])) ==
                  doctest::Approx(abserr[i]).epsilon(1e-3));
    }
    SUBCASE("Erlang-2 transform, the shape the transient QBD needs") {
        const double mu = 3.0;
        auto F = [mu](const std::complex<double>& s) {
            const std::complex<double> d = s + mu;
            return (mu * mu) / (d * d);
        };
        const std::vector<double> got = mam::matlab_ilt(F, ts, 100);
        const double matlab[5] = {0.666724651480946, 1.00405901368354, 0.448128494431847,
                                  0.0446546197483045, 1.45370867009018e-05};
        for (std::size_t i = 0; i < ts.size(); ++i)
            CHECK(got[i] == doctest::Approx(matlab[i]).epsilon(1e-10));
    }
    SUBCASE("transforms the quadrature handles near-exactly") {
        // 1/s -> 1 and 1/s^2 -> t are polynomial, where the rule is far more
        // accurate, so these DO admit a tight analytic bound.
        auto Fstep = [](const std::complex<double>& s) { return 1.0 / s; };
        for (double v : mam::matlab_ilt(Fstep, ts, 100)) CHECK(std::fabs(v - 1.0) < 1e-8);
        auto Framp = [](const std::complex<double>& s) { return 1.0 / (s * s); };
        const std::vector<double> ramp = mam::matlab_ilt(Framp, ts, 200);
        for (std::size_t i = 0; i < ts.size(); ++i)
            CHECK(std::fabs(ramp[i] - ts[i]) < 1e-6 * std::max(1.0, ts[i]));
    }
    SUBCASE("a larger budget buys accuracy") {
        const double a = 1.7;
        auto F = [a](const std::complex<double>& s) { return 1.0 / (s + a); };
        const std::vector<double> g100 = mam::matlab_ilt(F, ts, 100);
        const std::vector<double> g1000 = mam::matlab_ilt(F, ts, 1000);
        double w100 = 0.0, w1000 = 0.0;
        for (std::size_t i = 0; i < ts.size(); ++i) {
            w100 = std::max(w100, std::fabs(g100[i] - std::exp(-a * ts[i])));
            w1000 = std::max(w1000, std::fabs(g1000[i] - std::exp(-a * ts[i])));
        }
        CHECK(w1000 < w100);
    }
}

TEST_CASE("matlab_ilt refuses the unported weight families by name") {
    // euler and gaver are unreachable from this tree -- solver_mam_transient_qbd
    // is the only caller and takes the default -- and they are DIFFERENT
    // quadratures, not substitutes for cme. Refusing beats shipping an
    // untested branch that looks like a fallback.
    auto F = [](const std::complex<double>& s) { return 1.0 / (s + 1.0); };
    for (mam::IltMethod m : {mam::IltMethod::Euler, mam::IltMethod::Gaver}) {
        try {
            mam::matlab_ilt(F, std::vector<double>{1.0}, 100, m);
            FAIL("expected a refusal");
        } catch (const UnsupportedError& e) {
            CHECK(std::string(e.what()).find("only the CME weights are ported") !=
                  std::string::npos);
        }
    }
}

TEST_CASE("matlab_ilt refuses a non-positive evaluation time by name") {
    auto F = [](const std::complex<double>& s) { return 1.0 / s; };
    try {
        mam::matlab_ilt(F, std::vector<double>{0.0}, 100);
        FAIL("expected a refusal");
    } catch (const InputError& e) {
        CHECK(std::string(e.what()).find("strictly positive") != std::string::npos);
    }
}

// ---------------------------------------------------------------------------
// Laplace-domain transient QBD (solver_mam_transient_qbd)
// ---------------------------------------------------------------------------

TEST_CASE("the Laplace transient QBD reproduces the expm transient on M/M/1/N") {
    // TWO ENGINES, ONE MODEL. On a Poisson arrival with exponential service the
    // finite-buffer generator is small and solver_mam_ldqbd_transient computes
    // its transient by a matrix exponential, which is exact to machine
    // precision. The Laplace route reaches the same numbers through a wholly
    // different construction -- V(s,0,m) from mam_transient2, summed over
    // levels, inverted by the CME quadrature -- so the gap between them is the
    // INVERSION's error floor and nothing else. This is the oracle that
    // characterises the METHOD; the MATLAB comparison below pins the PORT.
    qn::Network<double> m = open_sq_cap("TQ1", Dd::exp_rate(1.0), Dd::exp_rate(2.0), 1.0, 5.0);
    mam::MamOptions o;
    o.timespan_start = 0.0;
    o.timespan_end = 10.0;
    // Called directly: the dispatch sends Poisson-plus-exponential to the fast
    // path, which is exactly why that model can serve as a cross-check here.
    CHECK(!mam::mam_transient_qbd_applicable(m.get_struct()));
    const mam::TranResult<double> lap = mam::solver_mam_transient_qbd(m.get_struct(), o);
    const mam::TranResult<double> ref = mam::solver_mam_ldqbd_transient(m.get_struct(), o);

    const mam::TranCurve<double>& Ql = lap.Qt[1][0];
    const mam::TranCurve<double>& Qr = ref.Qt[1][0];
    const mam::TranCurve<double>& Ul = lap.Ut[1][0];
    const mam::TranCurve<double>& Ur = ref.Ut[1][0];
    const mam::TranCurve<double>& Tl = lap.Tt[1][0];
    const mam::TranCurve<double>& Tr = ref.Tt[1][0];
    REQUIRE(Ql.times.size() == Qr.times.size());
    for (std::size_t i = 0; i < Ql.times.size(); ++i)
        CHECK(Ql.times[i] == doctest::Approx(Qr.times[i]).epsilon(1e-14));

    // t = 0 is not inverted at all: the system starts empty, so the initial
    // condition is written exactly rather than through a singular quadrature.
    CHECK(Ql.values.front() == 0.0);
    CHECK(Ul.values.front() == 0.0);
    CHECK(Tl.values.front() == 0.0);

    double worstQ = 0.0, worstU = 0.0, worstT = 0.0;
    for (std::size_t i = 1; i < Ql.values.size(); ++i) {
        worstQ = std::max(worstQ, rel(Ql.values[i], Qr.values[i]));
        worstU = std::max(worstU, rel(Ul.values[i], Ur.values[i]));
        worstT = std::max(worstT, rel(Tl.values[i], Tr.values[i]));
    }
    CHECK(worstQ < 1e-4);
    CHECK(worstU < 1e-4);
    CHECK(worstT < 1e-4);

    // The exact stationary law of M/M/1/5 at rho = 1/2 is an oracle for the
    // limit that neither engine uses: E[N] = 0.904761904761905,
    // U = 1 - p_0 = 0.492063492063492.
    CHECK(rel(Ql.values.back(), 0.904761904761905) < 3e-3);
    CHECK(rel(Ul.values.back(), 0.492063492063492) < 3e-3);
    // Exponential service makes the departure rate mu P(busy) at every instant,
    // and the two are summed from different weight vectors.
    for (std::size_t i = 0; i < Ql.values.size(); ++i)
        CHECK(std::fabs(Tl.values[i] - 2.0 * Ul.values[i]) < 1e-9);
}

TEST_CASE("with PH service the two engines differ ONLY in the initial phase") {
    // NOT a point-for-point oracle, and the reason is a genuine modelling
    // difference rather than an error in either engine.
    //
    // The Laplace route treats service as a MAP: the phase lives at every
    // level, including level 0, so pi(0) starts it in the TIME-STATIONARY
    // phase distribution map_prob({Ds0,Ds1}) -- (0.5, 0.5) for Erlang-2, i.e.
    // the first service is on average already half elapsed. The expm route
    // treats it as a renewal PH: level 0 is a single empty state and every
    // service begins at the entry vector alpha = (1, 0). The two conventions
    // coincide exactly for exponential service (memorylessness), which is why
    // the M/M/1/N case above agrees to 1e-4 -- and they cannot be made to
    // coincide for any genuine 2-phase PH, because the time-stationary phase
    // equals alpha only when all phase rates are equal.
    //
    // So the assertion is the one that IS true: the gap decays to zero. The
    // direction is predictable too -- a partially elapsed first service
    // completes sooner, so the Laplace curve sits BELOW the expm one early.
    //
    // No dispatched model can see this: mam_transient_qbd_applicable sends a
    // Poisson arrival with renewal PH service to the expm path, and the
    // Laplace path only ever runs on arrivals or services the expm path cannot
    // represent. The two conventions never answer the same question.
    qn::Network<double> m = open_sq_cap("TQ2", Dd::exp_rate(1.0), Dd::erlang(4.0, 2), 1.0, 4.0);
    mam::MamOptions o;
    o.timespan_start = 0.0;
    o.timespan_end = 8.0;
    const mam::TranResult<double> lap = mam::solver_mam_transient_qbd(m.get_struct(), o);
    const mam::TranResult<double> ref = mam::solver_mam_ldqbd_transient(m.get_struct(), o);
    const mam::TranCurve<double>& Ql = lap.Qt[1][0];
    const mam::TranCurve<double>& Qr = ref.Qt[1][0];
    REQUIRE(Ql.values.size() == Qr.values.size());

    // Below early, by the head start the elapsed service gives it.
    CHECK(Ql.values[8] < Qr.values[8]);
    CHECK(rel(Ql.values[8], Qr.values[8]) > 1e-2);
    // Monotone decay of the gap over the second half of the horizon, where the
    // initial condition has washed out and nothing else can be shrinking.
    for (std::size_t i = Ql.values.size() / 2; i + 1 < Ql.values.size(); ++i)
        CHECK(std::fabs(Ql.values[i + 1] - Qr.values[i + 1]) <=
              std::fabs(Ql.values[i] - Qr.values[i]) + 1e-12);
    // And agreement in the limit: same chain, same stationary law.
    const std::size_t n = Ql.values.size() - 1;
    CHECK(rel(Ql.values[n], Qr.values[n]) < 1e-3);
    CHECK(rel(lap.Ut[1][0].values[n], ref.Ut[1][0].values[n]) < 1e-3);
}

TEST_CASE("the infinite-buffer branch closes the level sums without truncating") {
    // No buffer bound at all: E[N](s) = pi0 V(s,0,1) (I-R)^-2 e is a CLOSED
    // form for the whole geometric tail, so the answer carries no truncation
    // parameter. The oracle is the exact M/M/1 stationary limit, rho/(1-rho)
    // and rho, which the construction never uses.
    qn::Network<double> m = open_sq_cap("TQ3", Dd::exp_rate(1.0), Dd::exp_rate(2.0), 1.0, -1.0);
    mam::MamOptions o;
    o.timespan_start = 0.0;
    o.timespan_end = 30.0;
    const mam::TranResult<double> r = mam::solver_mam_transient_qbd(m.get_struct(), o);
    const mam::TranCurve<double>& Q = r.Qt[1][0];
    const mam::TranCurve<double>& U = r.Ut[1][0];
    const mam::TranCurve<double>& Tp = r.Tt[1][0];
    CHECK(Q.values.front() == 0.0);
    CHECK(rel(Q.values.back(), 1.0) < 5e-3);
    CHECK(rel(U.values.back(), 0.5) < 5e-3);
    // Little's law is not imposed anywhere in the transform algebra.
    CHECK(rel(Tp.values.back(), 1.0) < 5e-3);
    for (std::size_t i = 0; i < Q.values.size(); ++i)
        CHECK(std::fabs(Tp.values[i] - 2.0 * U.values[i]) < 1e-8);
}

TEST_CASE("a correlated MMPP2 arrival now RUNS the Laplace transient") {
    // SUCCESSOR to the refusal half that used to live in "getTranAvg refuses
    // the Laplace regime and PH multiserver by name": this model is the one
    // mam_transient_qbd_applicable selects, and it is now answered rather than
    // refused. The oracle is the EXACT stationary MAP/MAP/1 solution of the
    // same queue (exact.mapmap1, itself pinned against MATLAB at 1e-8), which
    // shares no code with the transient transform.
    Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D0(0, 0) = -0.7; D0(0, 1) = 0.2; D0(1, 0) = 0.3; D0(1, 1) = -2.3;
    D1(0, 0) = 0.5;  D1(1, 1) = 2.0;
    qn::Network<double> m = open_sq_cap(
        "TQ4", Dd::map_dist(D0, D1, lang::ProcessType::MMPP2), Dd::exp_rate(3.0), 1.0, -1.0);
    REQUIRE(mam::mam_transient_qbd_applicable(m.get_struct()));
    mam::MamOptions o;
    o.timespan_start = 0.0;
    o.timespan_end = 60.0;
    const mam::TranResult<double> r = mam::solver_mam_get_tran_avg(m.get_struct(), o);
    const mam::TranCurve<double>& Q = r.Qt[1][0];
    const mam::TranCurve<double>& U = r.Ut[1][0];
    const mam::TranCurve<double>& Tp = r.Tt[1][0];
    CHECK(Q.values.front() == 0.0);
    CHECK(rel(Q.values.back(), 0.755144133935373) < 1e-2);
    CHECK(rel(U.values.back(), 0.366666666666667) < 1e-2);
    CHECK(rel(Tp.values.back(), 1.1) < 1e-2);
    // Exponential service again ties throughput to the busy probability.
    for (std::size_t i = 0; i < Q.values.size(); ++i)
        CHECK(std::fabs(Tp.values[i] - 3.0 * U.values[i]) < 1e-7);
}

TEST_CASE("the Laplace transient QBD refuses outside its regime by name") {
    mam::MamOptions o;
    o.timespan_start = 0.0;
    o.timespan_end = 5.0;
    // Multiserver: the level phase would have to carry a multiset of phases.
    qn::Network<double> mc = open_sq_cap("TQ5", Dd::exp_rate(1.0), Dd::exp_rate(3.0), 2.0, 5.0);
    try {
        mam::solver_mam_transient_qbd(mc.get_struct(), o);
        FAIL("expected a refusal");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("single-server") != std::string::npos);
    }
    // A closed model has no Source at all.
    qn::Network<double> mcl = closed_dq("TQ6", 3.0, Dd::exp_rate(1.0), Dd::exp_rate(2.0));
    try {
        mam::solver_mam_transient_qbd(mcl.get_struct(), o);
        FAIL("expected a refusal");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("open model") != std::string::npos);
    }
}


TEST_CASE("the Laplace transient curves match MATLAB point for point") {
    // ORACLE ROLES. The cases above measure the METHOD -- the exact expm
    // transient bounds the inversion's error, and exact.mapmap1 bounds the
    // limit. This one measures the PORT: the values are MATLAB
    // `solver_mam_transient_qbd` on the worktree tree, run on the same three
    // models, so any gap here is a transcription difference and nothing else.
    struct Pt { std::size_t i; double t, q, u, tp; };

    SUBCASE("finite buffer, exponential service") {
        qn::Network<double> m = open_sq_cap("MO1", Dd::exp_rate(1.0), Dd::exp_rate(2.0), 1.0, 5.0);
        mam::MamOptions o;
        o.timespan_start = 0.0;
        o.timespan_end = 10.0;
        const mam::TranResult<double> r = mam::solver_mam_transient_qbd(m.get_struct(), o);
        const std::vector<Pt> pts = {
            {0,  0.0,            0.0,             0.0,             0.0},
            {24, 2.4242424242,   0.738085810584,  0.448754105470,  0.897508210961},
            {49, 4.9494949495,   0.864120349655,  0.482074523523,  0.964149047046},
            {74, 7.4747474747,   0.894652938905,  0.489589205215,  0.979178410427},
            {99, 10.0,           0.902243241473,  0.491447156323,  0.982894312648}};
        REQUIRE(r.Qt[1][0].values.size() == 100);
        for (const Pt& p : pts) {
            CHECK(r.Qt[1][0].times[p.i] == doctest::Approx(p.t).epsilon(1e-9));
            CHECK(r.Qt[1][0].values[p.i] == doctest::Approx(p.q).epsilon(1e-9));
            CHECK(r.Ut[1][0].values[p.i] == doctest::Approx(p.u).epsilon(1e-9));
            CHECK(r.Tt[1][0].values[p.i] == doctest::Approx(p.tp).epsilon(1e-9));
        }
    }

    SUBCASE("finite buffer, Erlang-2 service") {
        qn::Network<double> m = open_sq_cap("MO2", Dd::exp_rate(1.0), Dd::erlang(4.0, 2), 1.0, 4.0);
        mam::MamOptions o;
        o.timespan_start = 0.0;
        o.timespan_end = 8.0;
        const mam::TranResult<double> r = mam::solver_mam_transient_qbd(m.get_struct(), o);
        const std::vector<Pt> pts = {
            {0,  0.0,            0.0,             0.0,             0.0},
            {19, 1.9240506329,   0.621077022849,  0.422092897962,  0.871156722810},
            {39, 3.9493670886,   0.750028423738,  0.476162786818,  0.954766478040},
            {59, 5.9746835443,   0.779623414754,  0.487463397414,  0.975022320496},
            {79, 8.0,            0.785927636686,  0.489759465053,  0.979488578066}};
        REQUIRE(r.Qt[1][0].values.size() == 80);
        for (const Pt& p : pts) {
            CHECK(r.Qt[1][0].times[p.i] == doctest::Approx(p.t).epsilon(1e-9));
            CHECK(r.Qt[1][0].values[p.i] == doctest::Approx(p.q).epsilon(1e-9));
            CHECK(r.Ut[1][0].values[p.i] == doctest::Approx(p.u).epsilon(1e-9));
            CHECK(r.Tt[1][0].values[p.i] == doctest::Approx(p.tp).epsilon(1e-9));
        }
    }

    SUBCASE("infinite buffer, correlated MMPP2 arrival") {
        Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
        D0(0, 0) = -0.7; D0(0, 1) = 0.2; D0(1, 0) = 0.3; D0(1, 1) = -2.3;
        D1(0, 0) = 0.5;  D1(1, 1) = 2.0;
        qn::Network<double> m = open_sq_cap(
            "MO3", Dd::map_dist(D0, D1, lang::ProcessType::MMPP2), Dd::exp_rate(3.0), 1.0, -1.0);
        mam::MamOptions o;
        o.timespan_start = 0.0;
        o.timespan_end = 60.0;
        const mam::TranResult<double> r = mam::solver_mam_transient_qbd(m.get_struct(), o);
        const std::vector<Pt> pts = {
            {0,   0.0,   0.0,             0.0,             0.0},
            {24,  14.4,  0.753834314377,  0.366518042012,  1.099554126023},
            {50,  30.0,  0.755135135424,  0.366665699544,  1.099997098644},
            {75,  45.0,  0.755143874234,  0.366666606398,  1.099999819205},
            {100, 60.0,  0.755143995878,  0.366666623797,  1.099999871343}};
        REQUIRE(r.Qt[1][0].values.size() == 101);
        for (const Pt& p : pts) {
            CHECK(r.Qt[1][0].times[p.i] == doctest::Approx(p.t).epsilon(1e-9));
            CHECK(r.Qt[1][0].values[p.i] == doctest::Approx(p.q).epsilon(1e-9));
            CHECK(r.Ut[1][0].values[p.i] == doctest::Approx(p.u).epsilon(1e-9));
            CHECK(r.Tt[1][0].values[p.i] == doctest::Approx(p.tp).epsilon(1e-9));
        }
    }
}
