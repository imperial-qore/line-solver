/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The exact multiserver phase-type chain inside the level-dependent QBD.
 *
 * `solver_mam_ldqbd` used to collapse the c parallel PH servers into ONE PH
 * process run at min(n,c) times its speed. That keeps the aggregate service
 * rate right and gets the exponential case exactly, but it forgets which phase
 * each busy server is in, which is not a rounding detail: it makes c servers
 * behave like one fast server whose remaining work is a single phase-type
 * variable. The reference measured ~1e-2 relative against SolverCTMC there.
 *
 * Since 2026-08-18 the level carries the MULTISET of the phases the min(n,c)
 * busy servers sit in (`ldqbd_mphc`), so the chain is exact and the bar in this
 * file is machine precision rather than a percentage: an exact chain against an
 * exact chain has nothing left to differ by.
 *
 * The c == 1 and exponential cases are pinned alongside, because the multiset
 * construction has to REPRODUCE them -- at one server the multiset is just the
 * phase, and at one phase there is no configuration coordinate at all. A
 * generator check on the blocks themselves catches a miscounted transition that
 * a mean-value comparison could absorb.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/mam/ldqbd_mphc.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/mam/solver_mam_ldqbd.h"

namespace qn = line::qn;
namespace mam = line::mam;
namespace ctmc = line::ctmc;
using line::lang::SchedStrategy;
using D = line::lang::Distrib<double>;

namespace {

/** Delay(Exp) -> Queue(FCFS, c servers), one closed class of N jobs. */
qn::Network<double> delay_queue(int N, double lamDelay, const D& service, double servers,
                                const std::vector<double>& lld) {
    qn::Network<double> m("ldqbd_mphc");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("Class1", static_cast<double>(N), d);
    m.set_service(d, c, D::exp_rate(lamDelay));
    m.set_service(q, c, service);
    m.set_number_of_servers(q, servers);
    if (!lld.empty()) m.set_load_dependence(q, lld);
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    return m;
}

void check_matches_ctmc(int N, double lamDelay, const D& service, double servers) {
    qn::Network<double> mm = delay_queue(N, lamDelay, service, servers, std::vector<double>());
    mam::MamOptions mo;
    mo.method = "ldqbd";
    const mam::LdqbdSolution<double> s = mam::solver_mam_ldqbd(mm.get_struct(), mo);

    qn::Network<double> mc = delay_queue(N, lamDelay, service, servers, std::vector<double>());
    const line::mva::AvgResult<double> e = ctmc::solver_ctmc_run_analyzer(mc.get_struct(), ctmc::CtmcOptions());

    for (std::size_t i = 0; i < s.sol.Q.rows(); ++i) {
        CAPTURE(i);
        CHECK(s.sol.Q(i, 0) == doctest::Approx(e.QN(i, 0)).epsilon(1e-9));
        CHECK(s.sol.U(i, 0) == doctest::Approx(e.UN(i, 0)).epsilon(1e-9));
        CHECK(s.sol.R(i, 0) == doctest::Approx(e.RN(i, 0)).epsilon(1e-9));
        CHECK(s.sol.Tp(i, 0) == doctest::Approx(e.TN(i, 0)).epsilon(1e-9));
    }
}

}  // namespace

TEST_CASE("ldqbd: Erlang service is exact at every server count") {
    for (double c : {1.0, 2.0, 3.0}) {
        CAPTURE(c);
        check_matches_ctmc(5, 0.8, D::erlang(2.0, 2), c);   // mean 1, order 2
    }
}

TEST_CASE("ldqbd: a three-phase Erlang is exact at every server count") {
    for (double c : {1.0, 2.0, 3.0}) {
        CAPTURE(c);
        check_matches_ctmc(4, 0.9, D::erlang(3.0, 3), c);   // mean 1, order 3
    }
}

TEST_CASE("ldqbd: HyperExp service is exact at every server count") {
    for (double c : {1.0, 2.0, 3.0}) {
        CAPTURE(c);
        check_matches_ctmc(5, 0.7, D::hyperexp(0.6, 2.0, 0.5), c);
    }
}

TEST_CASE("ldqbd: the exponential path is untouched") {
    // it never enters the multiset builder and must keep matching M/M/c exactly
    for (double c : {1.0, 2.0, 4.0}) {
        CAPTURE(c);
        check_matches_ctmc(6, 0.5, D::exp_rate(1.0), c);
    }
}

TEST_CASE("ldqbd: load dependence with phase-type service") {
    // the aggregate factor is shared over the busy servers, so at c = 1 it is
    // the single server running sf(n) times faster.
    //
    // UTILIZATION IS COMPARED HERE TOO, and that is the point of the case. MAM
    // used to report P(busy) = 1 - p(0) under load dependence, which reads a
    // server running alpha(n) times faster as no busier than one at its nominal
    // rate: 0.9587 against the CTMC's 0.6612 on exactly this model. It now
    // reports the work-based sum_n p(n)*sf(n)/max(c, max(alpha)), the CTMC's own
    // convention.
    const std::vector<double> lld{1.0, 1.5, 2.0, 2.5};
    qn::Network<double> mm = delay_queue(4, 1.0, D::erlang(2.0, 2), 1.0, lld);
    mam::MamOptions mo;
    mo.method = "ldqbd";
    const mam::LdqbdSolution<double> s = mam::solver_mam_ldqbd(mm.get_struct(), mo);

    qn::Network<double> mc = delay_queue(4, 1.0, D::erlang(2.0, 2), 1.0, lld);
    const line::mva::AvgResult<double> e = ctmc::solver_ctmc_run_analyzer(mc.get_struct(), ctmc::CtmcOptions());
    for (std::size_t i = 0; i < s.sol.Q.rows(); ++i) {
        CAPTURE(i);
        CHECK(s.sol.Q(i, 0) == doctest::Approx(e.QN(i, 0)).epsilon(1e-9));
        CHECK(s.sol.U(i, 0) == doctest::Approx(e.UN(i, 0)).epsilon(1e-9));
        CHECK(s.sol.Tp(i, 0) == doctest::Approx(e.TN(i, 0)).epsilon(1e-9));
    }
}

TEST_CASE("ph_multisets: the count and the order both matter") {
    // nchoosek(k+p-1, p-1) configurations, and k = 1 must give the identity rows
    // in phase order -- that ordering is what makes c = 1 coincide with plain
    // phase indexing, so the single-server blocks come out unchanged.
    CHECK(mam::ph_multisets(3, 0) == std::vector<std::vector<int> >{{0, 0, 0}});
    CHECK(mam::ph_multisets(3, 1) ==
          std::vector<std::vector<int> >{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
    CHECK(mam::ph_multisets(2, 2) == std::vector<std::vector<int> >{{2, 0}, {1, 1}, {0, 2}});
    CHECK(mam::ph_multisets(4, 3).size() == 20u);   // comb(6,3)
    CHECK(mam::ph_multisets(3, 4).size() == 15u);   // comb(6,2)
}

TEST_CASE("ldqbd_mphc: the blocks of every level sum to zero rows") {
    // that is what makes the block-tridiagonal matrix a generator, and it
    // catches a miscounted transition anywhere in the construction
    line::Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D0(0, 0) = -4.0; D0(0, 1) = 4.0; D0(1, 1) = -4.0;
    D1(1, 0) = 4.0;
    const std::vector<double> alpha{1.0, 0.0};
    const std::size_t Nlev = 6;
    std::vector<double> arr(Nlev + 1);
    for (std::size_t n = 0; n <= Nlev; ++n) arr[n] = (static_cast<double>(Nlev - n)) * 0.7;

    for (double c : {1.0, 2.0, 3.0}) {
        CAPTURE(c);
        const mam::LdqbdMphcBlocks<double> b =
            mam::ldqbd_mphc(D0, D1, alpha, c, arr, std::vector<double>());
        REQUIRE(b.Q0.size() == Nlev);
        REQUIRE(b.Q1.size() == Nlev + 1);
        REQUIRE(b.Q2.size() == Nlev + 1);
        for (std::size_t n = 0; n <= Nlev; ++n) {
            CAPTURE(n);
            // level size is comb(min(n,c)+p-1, p-1) with p = 2, i.e. min(n,c)+1
            const std::size_t want = std::min(n, static_cast<std::size_t>(c)) + 1;
            CHECK(b.Q1[n].rows() == want);
            for (std::size_t i = 0; i < b.Q1[n].rows(); ++i) {
                double rs = 0.0;
                for (std::size_t j = 0; j < b.Q1[n].cols(); ++j) rs += b.Q1[n](i, j);
                if (n < Nlev)
                    for (std::size_t j = 0; j < b.Q0[n].cols(); ++j) rs += b.Q0[n](i, j);
                if (n >= 1)
                    for (std::size_t j = 0; j < b.Q2[n].cols(); ++j) rs += b.Q2[n](i, j);
                CHECK(std::abs(rs) < 1e-12);
            }
        }
    }
}

TEST_CASE("ldqbd_mphc: sf(n) = min(n,c) reproduces the unscaled blocks exactly") {
    // this is what keeps the no-load-dependence path free of a stray division
    line::Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D0(0, 0) = -3.0; D0(0, 1) = 3.0; D0(1, 1) = -3.0;
    D1(1, 0) = 3.0;
    const std::vector<double> alpha{1.0, 0.0};
    const std::size_t Nlev = 5;
    const std::size_t c = 2;
    std::vector<double> arr(Nlev + 1), sf(Nlev);
    for (std::size_t n = 0; n <= Nlev; ++n) arr[n] = (static_cast<double>(Nlev - n)) * 0.9;
    for (std::size_t n = 1; n <= Nlev; ++n) sf[n - 1] = static_cast<double>(std::min(n, c));

    const mam::LdqbdMphcBlocks<double> plain =
        mam::ldqbd_mphc(D0, D1, alpha, static_cast<double>(c), arr, std::vector<double>());
    const mam::LdqbdMphcBlocks<double> scaled =
        mam::ldqbd_mphc(D0, D1, alpha, static_cast<double>(c), arr, sf);
    for (std::size_t n = 0; n <= Nlev; ++n) {
        CAPTURE(n);
        for (std::size_t i = 0; i < plain.Q1[n].rows(); ++i)
            for (std::size_t j = 0; j < plain.Q1[n].cols(); ++j)
                CHECK(plain.Q1[n](i, j) == scaled.Q1[n](i, j));
        for (std::size_t i = 0; i < plain.Q2[n].rows(); ++i)
            for (std::size_t j = 0; j < plain.Q2[n].cols(); ++j)
                CHECK(plain.Q2[n](i, j) == scaled.Q2[n](i, j));
    }
}
