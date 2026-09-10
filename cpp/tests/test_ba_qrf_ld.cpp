/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The load-dependent QRF arms on delay, multiserver and load-dependent models.
 *
 * WHAT CHANGED AND WHY IT IS NOT AN APPROXIMATION. `qrf.mmi.ld` and
 * `qrf.mmi.linear` carry a scaling alpha(i,n) that multiplies every rate out of
 * station i at population n, completions and background phase changes alike.
 * That is exactly the rate law of an infinite server (alpha = n), of a c-server
 * station (alpha = min(n,c)) and of limited load dependence, so deriving alpha
 * from the model makes the relaxed chain the model's OWN chain rather than an
 * approximation of it. The two arms therefore serve models the rest of the QRF
 * family still refuses.
 *
 * THE ORACLE IS EXACTNESS, NOT CLOSENESS. At M = 2 the pairwise joint of a
 * closed chain is fully determined by the marginal, so the QRF polytope is
 * tight and the answer must EQUAL the exact CTMC one. Every two-station fixture
 * here is asserted at 1e-9 against SolverCTMC, which pins the alpha derivation,
 * the BN readout and the utilization normalizer at once; a looser tolerance
 * would let a wrong alpha through.
 */

#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/sn/sn_to_qrf_alpha.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ba/solver_ba_runner.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;
const double kTol = 1e-9;

/** Two stations in a cycle; station 0 is a Delay when delay0. */
qn::Network<double> cqn(double N, double c0, bool delay0,
                        const std::vector<double>& lld0 = std::vector<double>(),
                        int phases0 = 1) {
    qn::Network<double> m("qrfLd");
    const std::size_t s0 =
        delay0 ? m.add_delay("D0") : m.add_queue("Q0", SchedStrategy::FCFS);
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C", N, s0);
    if (!delay0 && c0 > 1.0) m.set_number_of_servers(s0, c0);
    if (!lld0.empty()) m.set_load_dependence(s0, lld0);
    if (phases0 > 1)
        m.set_service(s0, c, D::erlang(static_cast<double>(phases0), phases0));
    else
        m.set_service(s0, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(1.5));
    qn::RoutingMatrix<double> P;
    P.set(s0, q1, 1.0);
    P.set(q1, s0, 1.0);
    m.link(P);
    return m;
}

mva::AvgResult<double> run_ba(qn::Network<double>& m, const std::string& method) {
    ba::BaOptions opt;
    opt.method = method;
    return ba::solver_ba_run_analyzer(m.get_struct(), opt);
}

mva::AvgResult<double> run_ctmc(qn::Network<double>& m) {
    return ctmc::solver_ctmc_run_analyzer(m.get_struct(), ctmc::CtmcOptions());
}

void check_same_as_ctmc(qn::Network<double>& m, const std::string& method) {
    const mva::AvgResult<double> b = run_ba(m, method);
    const mva::AvgResult<double> e = run_ctmc(m);
    for (int i = 0; i < b.QN.rows(); ++i) {
        CHECK(b.QN(i, 0) == doctest::Approx(e.QN(i, 0)).epsilon(kTol));
        CHECK(b.UN(i, 0) == doctest::Approx(e.UN(i, 0)).epsilon(kTol));
        CHECK(b.TN(i, 0) == doctest::Approx(e.TN(i, 0)).epsilon(kTol));
    }
}

}  // namespace

TEST_CASE("qrf alpha: min(n,c) at a multiserver station") {
    qn::Network<double> m = cqn(4, 3, false);
    const sn::QrfAlpha a = sn::sn_to_qrf_alpha(m.get_struct());
    CHECK(a.msg.empty());
    CHECK(a.ld);
    CHECK(a.alpha[0] == std::vector<double>({1, 2, 3, 3}));
    CHECK(a.alpha[1] == std::vector<double>({1, 1, 1, 1}));
    CHECK(a.peak[0] == doctest::Approx(3.0));
    CHECK(a.peak[1] == doctest::Approx(1.0));
}

TEST_CASE("qrf alpha: the peak is the DECLARED servers, not the reachable maximum") {
    // c = 3 with N = 2: alpha reaches 2, but the station still has three servers
    // and LINE reports U = T*S/c. Normalizing by max(alpha) overstates U by 3/2.
    qn::Network<double> m = cqn(2, 3, false);
    const sn::QrfAlpha a = sn::sn_to_qrf_alpha(m.get_struct());
    CHECK(a.msg.empty());
    CHECK(a.alpha[0] == std::vector<double>({1, 2}));
    CHECK(a.peak[0] == doctest::Approx(3.0));
}

TEST_CASE("qrf alpha: n at a delay") {
    qn::Network<double> m = cqn(3, 1, true);
    const sn::QrfAlpha a = sn::sn_to_qrf_alpha(m.get_struct());
    CHECK(a.msg.empty());
    CHECK(a.ld);
    CHECK(a.alpha[0] == std::vector<double>({1, 2, 3}));
    CHECK(std::isinf(a.peak[0]));
    CHECK(a.peak[1] == doctest::Approx(1.0));
}

TEST_CASE("qrf alpha: lld composes with the server count") {
    qn::Network<double> m = cqn(3, 2, false, std::vector<double>({1.0, 1.5, 2.0}));
    const sn::QrfAlpha a = sn::sn_to_qrf_alpha(m.get_struct());
    CHECK(a.msg.empty());
    CHECK(a.ld);
    CHECK(a.alpha[0] == std::vector<double>({1.0, 3.0, 4.0}));
    // The peak is the PRODUCT of the server count and the reachable lld peak.
    CHECK(a.peak[0] == doctest::Approx(4.0));
}

TEST_CASE("qrf alpha: all ones on a single-server model") {
    // The load-independent model must reach the ld arms unchanged, which is what
    // keeps this change a no-op there.
    qn::Network<double> m = cqn(3, 1, false);
    const sn::QrfAlpha a = sn::sn_to_qrf_alpha(m.get_struct());
    CHECK(a.msg.empty());
    CHECK_FALSE(a.ld);
    for (std::size_t i = 0; i < a.alpha.size(); ++i)
        for (std::size_t n = 0; n < a.alpha[i].size(); ++n) CHECK(a.alpha[i][n] == 1.0);
    CHECK(a.peak[0] == doctest::Approx(1.0));
}

TEST_CASE("qrf alpha: PH where several jobs are served at once is refused") {
    // The QRF local state carries ONE phase per station, which describes one job
    // in service and no more, so a PH multiserver would answer a different chain
    // and the number would bound nothing.
    for (int which = 0; which < 2; ++which) {
        qn::Network<double> m =
            cqn(2, which == 0 ? 2 : 1, which == 1, std::vector<double>(), 2);
        const sn::QrfAlpha a = sn::sn_to_qrf_alpha(m.get_struct());
        CHECK(a.msg.find("one phase per station") != std::string::npos);
        CHECK(a.ld);  // ld must stay set through the refusal
        CHECK_THROWS(run_ba(m, "qrf.mmi.ld"));
    }
}

TEST_CASE("qrf: the alpha-free arms refuse and name the arms that serve") {
    const char* methods[] = {"qr", "qrf.mmi", "qrf.mem", "qrf.bethe"};
    for (int mi = 0; mi < 4; ++mi) {
        qn::Network<double> a = cqn(2, 2, false);
        qn::Network<double> b = cqn(2, 1, true);
        qn::Network<double> c = cqn(2, 1, false, std::vector<double>({1.0, 2.0}));
        CHECK_THROWS(run_ba(a, methods[mi]));
        CHECK_THROWS(run_ba(b, methods[mi]));
        CHECK_THROWS(run_ba(c, methods[mi]));
    }
}

TEST_CASE("qrf: list_valid_methods keeps the ld arms on a delay model") {
    // A caller enumerating the list must still see the only two bound methods
    // the model has; dropping them with the rest would hide them.
    qn::Network<double> m = cqn(2, 1, true);
    const std::vector<std::string> methods = ba::list_valid_methods(m.get_struct());
    bool ld = false, linear = false, mmi = false, bas = false;
    for (std::size_t i = 0; i < methods.size(); ++i) {
        if (methods[i] == "qrf.mmi.ld") ld = true;
        if (methods[i] == "qrf.mmi.linear") linear = true;
        if (methods[i] == "qrf.mmi") mmi = true;
        if (methods[i] == "qrf.bas") bas = true;
    }
    CHECK(ld);
    CHECK(linear);
    CHECK_FALSE(mmi);
    CHECK_FALSE(bas);
}

TEST_CASE("qrf ld: the two-station answers are exact") {
    // M = 2 makes the polytope tight, so these are equalities, not bounds.
    const char* methods[] = {"qrf.mmi.ld", "qrf.mmi.linear"};
    for (int mi = 0; mi < 2; ++mi) {
        qn::Network<double> a = cqn(3, 1, true);
        check_same_as_ctmc(a, methods[mi]);
        qn::Network<double> b = cqn(3, 2, false);
        check_same_as_ctmc(b, methods[mi]);
        qn::Network<double> c = cqn(4, 3, false);
        check_same_as_ctmc(c, methods[mi]);
        qn::Network<double> d = cqn(3, 1, false, std::vector<double>({1.0, 1.5, 2.0}));
        check_same_as_ctmc(d, methods[mi]);
        qn::Network<double> e = cqn(3, 1, false);
        check_same_as_ctmc(e, methods[mi]);
    }
}

TEST_CASE("qrf ld: a multiserver at population one is the single-server model") {
    // min(1,c) = 1, so the two chains are identical and only the utilization
    // normalizer differs. No CTMC needed: the oracle is the other run.
    const char* methods[] = {"qrf.mmi.ld", "qrf.mmi.linear"};
    for (int mi = 0; mi < 2; ++mi) {
        qn::Network<double> one = cqn(1, 1, false);
        qn::Network<double> three = cqn(1, 3, false);
        const mva::AvgResult<double> r1 = run_ba(one, methods[mi]);
        const mva::AvgResult<double> r3 = run_ba(three, methods[mi]);
        for (int i = 0; i < r1.QN.rows(); ++i) {
            CHECK(r1.QN(i, 0) == doctest::Approx(r3.QN(i, 0)).epsilon(kTol));
            CHECK(r1.TN(i, 0) == doctest::Approx(r3.TN(i, 0)).epsilon(kTol));
        }
        CHECK(r1.UN(0, 0) == doctest::Approx(r3.UN(0, 0) * 3.0).epsilon(kTol));
    }
}
