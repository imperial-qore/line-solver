/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The AMVA method matrix of solver_amva.m, on two closed models.
 *
 * Every number here was produced by MATLAB (`SolverMVA(model, opts).getAvg()`
 * with `opts.method` and `opts.config.multiserver` set as the case says) and is
 * asserted to the digits MATLAB prints. That is the right oracle for this
 * layer: the methods are APPROXIMATIONS, so there is no closed form to check
 * them against, and what parity means is agreeing with the reference's
 * approximation rather than with the exact answer.
 *
 * Two of the values look wrong and are the reference's:
 *   sqni    saturates on model A (U = 1, X at the bottleneck bound), because
 *           N/(D+Z) and 1/D coincide there at 2
 *   krzesinski reports U = 1.45 at a two-server station -- above one -- since
 *           the Linearizer-MS branch it selects returns the utilization per
 *           STATION and the reference passes it through undivided
 * Both are asserted as they stand: a port that "fixed" either would disagree
 * with MATLAB, the JAR and Python at once.
 */

#include <algorithm>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/** Delay(Z=1) + FCFS Queue(D=0.5), one closed class of 3. */
qn::Network<double> model_a() {
    qn::Network<double> m("A");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** The same with two servers at the queue and a population of 5. */
qn::Network<double> model_b() {
    qn::Network<double> m("B");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 5.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    m.set_number_of_servers(q, 2.0);
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

struct Metrics {
    double Qd, Qq, Uq, Rq;
};

Metrics solve(qn::Network<double>& m, const std::string& method,
              const std::string& multiserver = "default") {
    mva::MvaOptions opt;
    opt.method = method;
    opt.multiserver = multiserver;
    Matrix<double> init;
    const mva::MvaSolution<double> r = mva::solver_mva_analyzer(m.get_struct(), opt, init);
    Metrics out;
    out.Qd = r.Q(0, 0);
    out.Qq = r.Q(1, 0);
    out.Uq = r.U(1, 0);
    out.Rq = r.R(1, 0);
    return out;
}

void check(const Metrics& got, double Qd, double Qq, double Uq, double Rq) {
    CHECK(got.Qd == doctest::Approx(Qd).epsilon(1e-9));
    CHECK(got.Qq == doctest::Approx(Qq).epsilon(1e-9));
    CHECK(got.Uq == doctest::Approx(Uq).epsilon(1e-9));
    CHECK(got.Rq == doctest::Approx(Rq).epsilon(1e-9));
}

TEST_CASE("every AMVA method agrees with MATLAB on a single-server closed model") {
    qn::Network<double> m = model_a();
    check(solve(m, "exact"), 1.578947368, 1.421052632, 0.789473684, 0.900000000);
    check(solve(m, "default"), 1.578947368, 1.421052632, 0.789473684, 0.900000000);
    check(solve(m, "bs"), 1.499982834, 1.500017166, 0.749991417, 1.000022889);
    check(solve(m, "aql"), 1.580544267, 1.419554928, 0.790272134, 0.898143100);
    check(solve(m, "tay"), 1.585456789, 1.414543211, 0.792728394, 0.892199157);
    // SCAT is Linearizer with ONE Delta refresh instead of three, so it must sit
    // between bs (no refresh) and lin (three) and equal neither. It does: 9
    // inner iterations against lin's 27 on this model.
    check(solve(m, "scat"), 1.593017294, 1.406982706, 0.796508647, 0.883218727);
    check(solve(m, "lin"), 1.580792642, 1.419207358, 0.790396321, 0.897782113);
    check(solve(m, "gflin"), 1.595965230, 1.404034770, 0.797982615, 0.879740199);
    check(solve(m, "egflin"), 1.583664096, 1.416335904, 0.791832048, 0.894341109);
    check(solve(m, "qd"), 1.499999844, 1.499999667, 0.750000044, 0.999999719);
    check(solve(m, "ab"), 1.592853511, 1.407146489, 0.796426756, 0.883412366);
    check(solve(m, "schmidt"), 1.578947368, 1.421052632, 0.789473684, 0.900000000);
    check(solve(m, "schmidt-ext"), 1.578947368, 1.421052632, 0.789473684, 0.900000000);
    // sqni saturates here, exactly as the reference does
    check(solve(m, "sqni"), 2.000000000, 1.000000000, 1.000000000, 0.500000000);
}

TEST_CASE("the amva.* spellings resolve to the same methods") {
    qn::Network<double> m = model_a();
    check(solve(m, "amva.bs"), 1.499982834, 1.500017166, 0.749991417, 1.000022889);
    check(solve(m, "amva.aql"), 1.580544267, 1.419554928, 0.790272134, 0.898143100);
    check(solve(m, "amva.lin"), 1.580792642, 1.419207358, 0.790396321, 0.897782113);
    // `amva.tay` reached this line through a REFERENCE DEFECT that the port
    // found and that has since been fixed in MATLAB (2026-07-29). Recorded
    // because the evidence is what justified the fix, and because a silent
    // alias is exactly the kind of thing that regresses unnoticed.
    // SolverMVA.m:63 advertised 'amva.tay' in listValidMethods while
    // solver_amva.m carried ten `case 'amva.X'` arms and none for tay, so the
    // string survived the switch and fell through to a qd-equivalent path: an
    // ADVERTISED method silently computing a different one. Measured on this
    // model, amva.tay gave 1.499999844 / 1.499999667 / 0.750000044 /
    // 0.999999719 in 31 iterations -- the `qd` row asserted above -- against
    // tay's 1.585456789 in 8. MATLAB now has the missing arm and both
    // codebases alias amva.tay -> tay, so the row below is the SAME as `tay`
    // by construction, which is the property this line exists to hold.
    check(solve(m, "amva.tay"), 1.585456789, 1.414543211, 0.792728394, 0.892199157);
    check(solve(m, "amva.scat"), 1.593017294, 1.406982706, 0.796508647, 0.883218727);
}

TEST_CASE("the multiserver rules agree with MATLAB on a two-server closed model") {
    qn::Network<double> m = model_b();
    check(solve(m, "exact"), 3.005181347, 1.994818653, 0.751295337, 0.663793103);
    check(solve(m, "default"), 3.005181347, 1.994818653, 0.751295337, 0.663793103);
    check(solve(m, "lin"), 2.865913405, 2.134086262, 0.716478396, 0.744644316);
    check(solve(m, "ab"), 3.233444813, 1.766555187, 0.808361203, 0.546338438);
    check(solve(m, "schmidt"), 2.727272727, 2.272727273, 0.681818182, 0.833333333);
    check(solve(m, "schmidt-ext"), 2.998018587, 2.001981413, 0.749504647, 0.667768179);
    // Qd is the post-fix reference: it used to read 3.454861111, the delay
    // charged the Seidmann think time while the queue was also given the folded
    // population back, so sum(Q) came to 5.691 against N = 5
    check(solve(m, "bs"), 2.763888889, 2.236111111, 0.690972222, 0.809045226);
    // Unlike aql/qsa/tay, SCAT takes a multiserver model: it reaches the
    // algorithm Seidmann-scaled, exactly as bs does, so it is not refused here.
    check(solve(m, "scat"), 2.883835883, 2.116164117, 0.720958971, 0.733801854);
    // RE-RECORDED 2026-08-04, from MATLAB on this very model: the Conway rule
    // moved when conwayms stopped converging on the queue lengths alone and
    // started watching the marginal probabilities too. The old row
    // (2.787478236, 2.212521764, 0.696869559, 0.793735978) is that reference's
    // pre-fix output.
    check(solve(m, "lin", "conway"), 2.78751949, 2.21248051, 0.6968798725, 0.7937094316);
    // the reference's undivided utilization, above one; see the header note.
    // RE-RECORDED 2026-07-31: krzesinski is the rule that reaches
    // pfqn_linearizerms, so this row moved with that reference's own fix (the
    // partially-idle-server term, the Estimate marginals and the closed-form
    // marginal fixed point). Values from MATLAB on this very model, not from
    // this port's output.
    check(solve(m, "lin", "krzesinski"), 3.019885265, 1.980114735, 1.509942633, 0.655692041);
}

TEST_CASE("a method outside its domain is refused by name, not approximated") {
    SUBCASE("AQL has no multiserver correction") {
        qn::Network<double> m = model_b();
        CHECK_THROWS_AS(solve(m, "aql"), UnsupportedError);
    }
    SUBCASE("Tay's elasticity equations are derived for single servers") {
        // solver_amva.m:147-149 raises rather than approximating, and so does
        // this port: there is no multiserver correction to fall back on.
        qn::Network<double> m = model_b();
        CHECK_THROWS_AS(solve(m, "tay"), UnsupportedError);
        CHECK_THROWS_AS(solve(m, "amva.tay"), UnsupportedError);
    }
    SUBCASE("sqni needs one queue and one delay") {
        qn::Network<double> m("three");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
        const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
        const std::size_t c = m.add_closed_class("C1", 3.0, d);
        m.set_service(d, c, D::exp_rate(1.0));
        m.set_service(q1, c, D::exp_rate(2.0));
        m.set_service(q2, c, D::exp_rate(2.0));
        qn::RoutingMatrix<double> P;
        P.set(d, q1, 0.5);
        P.set(d, q2, 0.5);
        P.set(q1, d, 1.0);
        P.set(q2, d, 1.0);
        m.link(P);
        CHECK_THROWS_AS(solve(m, "sqni"), UnsupportedError);
    }
}

// ---------------------------------------------------------------------------
// The rest of solver_amvald_forward: the multiserver rules the Linearizer path
// does not cover, the high-variance correction, the two disciplines with a
// parameter, and the load- and class-dependent scalings.
//
// Every number below is MATLAB's, from `SolverMVA(model, opts).getAvg()` with
// `opts.method`, `opts.config.multiserver`, `opts.config.highvar` and
// `opts.config.np_priority` set as the case says.
// ---------------------------------------------------------------------------

/** Q of every (station, class) pair, which is what these cases compare. */
std::vector<double> qcol(qn::Network<double>& m, const std::string& method,
                         const std::string& multiserver = "default",
                         const std::string& highvar = "default",
                         const std::string& np_priority = "default") {
    mva::MvaOptions opt;
    opt.method = method;
    opt.multiserver = multiserver;
    opt.highvar = highvar;
    opt.np_priority = np_priority;
    Matrix<double> init;
    const mva::MvaSolution<double> r = mva::solver_mva_analyzer(m.get_struct(), opt, init);
    std::vector<double> out;
    for (std::size_t i = 0; i < r.Q.rows(); ++i)
        for (std::size_t c = 0; c < r.Q.cols(); ++c) out.push_back(r.Q(i, c));
    return out;
}

void check_all(const std::vector<double>& got, const std::vector<double>& want) {
    REQUIRE(got.size() == want.size());
    for (std::size_t i = 0; i < got.size(); ++i)
        CHECK(got[i] == doctest::Approx(want[i]).epsilon(1e-9));
}

TEST_CASE("the softmin and suri multiserver rules agree with MATLAB") {
    qn::Network<double> m = model_b();
    check_all(qcol(m, "qd", "softmin"), {3.409434599, 1.590564785});
    check_all(qcol(m, "qd", "suri"), {2.970076107, 2.029924452});
}

TEST_CASE("a DPS station charges each class by the weight ratio, as MATLAB does") {
    qn::Network<double> m("dps");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::DPS);
    const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 1.0, d);
    m.set_service(d, c1, D::exp_rate(3.0));
    m.set_service(d, c2, D::exp_rate(0.5));
    m.set_service(q, c1, D::exp_rate(1.0));
    m.set_service(q, c2, D::exp_rate(2.0));
    m.set_sched_param(q, c1, 1.0);
    m.set_sched_param(q, c2, 5.0);
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c1, c1, q, d, 1.0);
    P.set(c2, c2, d, q, 1.0);
    P.set(c2, c2, q, d, 1.0);
    m.link(P);
    check_all(qcol(m, "qd"), {0.190021753, 0.745990858, 1.809977400, 0.254009087});
}

TEST_CASE("a HOL station applies the priority scaling under both np_priority rules") {
    qn::Network<double> m("hol");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::HOL);
    const std::size_t hi = m.add_closed_class("Hi", 2.0, d, 0);
    const std::size_t lo = m.add_closed_class("Lo", 2.0, d, 1);
    m.set_service(d, hi, D::exp_rate(1.0));
    m.set_service(d, lo, D::exp_rate(1.0));
    m.set_service(q, hi, D::exp_rate(2.0));
    m.set_service(q, lo, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(hi, hi, d, q, 1.0);
    P.set(hi, hi, q, d, 1.0);
    P.set(lo, lo, d, q, 1.0);
    P.set(lo, lo, q, d, 1.0);
    m.link(P);
    // Regenerated from MATLAB after the HOL wait was split into the queued equal-or-higher
    // priority backlog and the overtaking term (CTMC: 1.014493, 0.795031, 0.985507, 1.204969)
    const std::vector<double> want{1.058118652, 0.618554316, 0.941881880, 1.381444888};
    check_all(qcol(m, "qd"), want);
    // Chandy-Lakshmi and the shadow server coincide wherever tau is zero, which
    // is every method that does not run the Linearizer recursion.
    check_all(qcol(m, "qd", "default", "default", "shadow"), want);
}

TEST_CASE("load dependence scales the demand and highvar changes the answer") {
    qn::Network<double> m("ld");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 4.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::erlang_fit(0.5, 0.25));  // SCV 1/4, so hvmva bites
    m.set_load_dependence(q, std::vector<double>{1.0, 1.5, 1.8, 2.0});
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    check_all(qcol(m, "qd"), {2.333333134, 1.666666826});
    check_all(qcol(m, "qd", "default", "hvmva"), {2.519129309, 1.480871171});
}

TEST_CASE("class dependence reaches the AMVA through pfqn_cdfun") {
    const double N = 16.0, c = 2.0;
    qn::Network<double> m("cd");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t j1 = m.add_closed_class("Class1", N, d);
    const std::size_t j2 = m.add_closed_class("Class2", N / 2.0, d);
    m.set_service(d, j1, D::exp_mean(1.0));
    m.set_service(d, j2, D::exp_mean(2.0));
    m.set_service(q, j1, D::exp_mean(1.5));
    m.set_service(q, j2, D::exp_mean(2.5));
    // multi-server for class-1 jobs only, as the MATLAB example declares it.
    // The peak is max_n beta(n) = c, and it is MANDATORY: without it SolverMVA
    // has nothing to normalize the utilization column by and refuses the model.
    m.set_class_dependence(q,
                           [c](const std::vector<double>& ni) {
                               return std::vector<double>{std::min(ni[0], c)};
                           },
                           std::vector<double>(1, c));
    qn::RoutingMatrix<double> P;
    P.set(j1, j1, d, q, 1.0);
    P.set(j1, j1, q, d, 1.0);
    P.set(j2, j2, d, q, 1.0);
    P.set(j2, j2, q, d, 1.0);
    m.link(P);
    // re-recorded 2026-08-14 against MATLAB ld_joint_dependence, after 58eb739f1
    check_all(qcol(m, "qd"), {0.889959830, 0.527875366, 15.110040622, 7.472124857});
}

TEST_CASE("the Wang-Sevcik queue-line and fraction-line estimators run") {
    qn::Network<double> m = model_a();
    check_all(qcol(m, "qli"), {1.499999844, 1.499999667});
    check_all(qcol(m, "fli"), {1.499999844, 1.499999667});
}

}  // namespace
