/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * SolverBA, the bound-analysis solver.
 *
 * Every expected number below is MATLAB's
 * `SolverBA(model,'method',<name>).getAvg` on the same model, so these are
 * parity assertions and not a record of what this port happens to produce. The
 * bracket checks against exact MVA are the part the reference values cannot
 * fake: a transcription error that agreed with MATLAB would still have to
 * bracket the exact answer, and a bound that does not is wrong whatever it
 * agrees with.
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/api/pfqn/pfqn_scb.h"
#include "line/solvers/ba/solver_ba_runner.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;
const double kEps = 1e-9;

mva::AvgResult<double> run_ba(qn::Network<double>& m, const std::string& method, int level = 2) {
    ba::BaOptions opt;
    opt.method = method;
    opt.level = level;
    return ba::solver_ba_run_analyzer(m.get_struct(), opt);
}

mva::AvgResult<double> run_mva(qn::Network<double>& m) {
    mva::MvaOptions opt;
    opt.method = "exact";
    return mva::solver_mva_run_analyzer(m.get_struct(), opt, Matrix<double>());
}

/** Delay(1) -> FCFS Queue1(2) -> FCFS Queue2(3), N = 4. */
qn::Network<double> model_a() {
    qn::Network<double> m("cqnA");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 4.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, d, 1.0);
    m.link(P);
    return m;
}

/** FCFS Queue1(2) -> Queue2(3) -> Queue3(5), N = 3, no delay station. */
qn::Network<double> model_b() {
    qn::Network<double> m("cqnB");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
    const std::size_t q3 = m.add_queue("Queue3", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, q1);
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(3.0));
    m.set_service(q3, c, D::exp_rate(5.0));
    qn::RoutingMatrix<double> P;
    P.set(q1, q2, 1.0);
    P.set(q2, q3, 1.0);
    P.set(q3, q1, 1.0);
    m.link(P);
    return m;
}

/** Delay + PS Queue1 + PS Queue2, two closed classes of 2 and 3 jobs. */
qn::Network<double> model_c() {
    qn::Network<double> m("cqnC");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::PS);
    const std::size_t k1 = m.add_closed_class("C1", 2.0, d);
    const std::size_t k2 = m.add_closed_class("C2", 3.0, d);
    m.set_service(d, k1, D::exp_rate(1.0));
    m.set_service(d, k2, D::exp_rate(2.0));
    m.set_service(q1, k1, D::exp_rate(4.0));
    m.set_service(q1, k2, D::exp_rate(3.0));
    m.set_service(q2, k1, D::exp_rate(5.0));
    m.set_service(q2, k2, D::exp_rate(6.0));
    qn::RoutingMatrix<double> P;
    for (std::size_t k : {k1, k2}) {
        P.set(k, k, d, q1, 1.0);
        P.set(k, k, q1, q2, 1.0);
        P.set(k, k, q2, d, 1.0);
    }
    m.link(P);
    return m;
}

/** Model A with two servers at Queue1, which only the ssd family accepts. */
qn::Network<double> model_d() {
    qn::Network<double> m("cqnD");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
    m.set_number_of_servers(q1, 2.0);
    const std::size_t c = m.add_closed_class("C1", 4.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, d, 1.0);
    m.link(P);
    return m;
}

TEST_CASE("ba: the asymptotic bounds on a closed single-class network with a delay") {
    qn::Network<double> m = model_a();
    SUBCASE("aba.upper") {
        const mva::AvgResult<double> r = run_ba(m, "aba.upper");
        CHECK(r.QN(0, 0) == doctest::Approx(2.0).epsilon(kEps));
        CHECK(r.QN(1, 0) == doctest::Approx(4.0).epsilon(kEps));
        CHECK(r.QN(2, 0) == doctest::Approx(2.66666666667).epsilon(kEps));
        CHECK(r.UN(0, 0) == doctest::Approx(2.0).epsilon(kEps));
        CHECK(r.UN(1, 0) == doctest::Approx(1.0).epsilon(kEps));
        CHECK(r.UN(2, 0) == doctest::Approx(0.666666666667).epsilon(kEps));
        CHECK(r.RN(0, 0) == doctest::Approx(1.0).epsilon(kEps));
        CHECK(r.RN(1, 0) == doctest::Approx(2.0).epsilon(kEps));
        CHECK(r.RN(2, 0) == doctest::Approx(1.33333333333).epsilon(kEps));
        CHECK(r.TN(0, 0) == doctest::Approx(2.0).epsilon(kEps));
        CHECK(r.TN(2, 0) == doctest::Approx(2.0).epsilon(kEps));
        // the arrival rate the runner derives from sn.rt, not from a zero matrix
        CHECK(r.AN(1, 0) == doctest::Approx(2.0).epsilon(kEps));
        CHECK(r.CN[0] == doctest::Approx(4.33333333333).epsilon(kEps));
        CHECK(r.XN[0] == doctest::Approx(2.0).epsilon(kEps));
        // runAnalyzer.m writes a zero residence time: the families define none
        CHECK(r.WN(1, 0) == doctest::Approx(0.0));
    }
    SUBCASE("aba.lower") {
        const mva::AvgResult<double> r = run_ba(m, "aba.lower");
        CHECK(r.QN(0, 0) == doctest::Approx(0.923076923077).epsilon(kEps));
        CHECK(r.QN(1, 0) == doctest::Approx(0.461538461538).epsilon(kEps));
        CHECK(r.QN(2, 0) == doctest::Approx(0.307692307692).epsilon(kEps));
        CHECK(r.RN(1, 0) == doctest::Approx(0.5).epsilon(kEps));
        CHECK(r.TN(0, 0) == doctest::Approx(0.923076923077).epsilon(kEps));
        CHECK(r.CN[0] == doctest::Approx(1.83333333333).epsilon(kEps));
    }
    SUBCASE("bjb, the balanced job bounds") {
        const mva::AvgResult<double> u = run_ba(m, "bjb.upper");
        CHECK(u.QN(0, 0) == doctest::Approx(1.66561514196).epsilon(kEps));
        CHECK(u.QN(1, 0) == doctest::Approx(3.33123028391).epsilon(kEps));
        CHECK(u.QN(2, 0) == doctest::Approx(2.22082018927).epsilon(kEps));
        CHECK(u.TN(0, 0) == doctest::Approx(1.66561514196).epsilon(kEps));
        CHECK(u.CN[0] == doctest::Approx(2.90476190476).epsilon(kEps));
        const mva::AvgResult<double> l = run_ba(m, "bjb.lower");
        CHECK(l.QN(0, 0) == doctest::Approx(1.37704918033).epsilon(kEps));
        CHECK(l.QN(1, 0) == doctest::Approx(0.688524590164).epsilon(kEps));
        CHECK(l.TN(0, 0) == doctest::Approx(1.37704918033).epsilon(kEps));
        CHECK(l.CN[0] == doctest::Approx(2.40151515152).epsilon(kEps));
    }
    SUBCASE("pb, the proportional bounds") {
        const mva::AvgResult<double> u = run_ba(m, "pb.upper");
        CHECK(u.QN(0, 0) == doctest::Approx(1.65).epsilon(kEps));
        CHECK(u.QN(1, 0) == doctest::Approx(3.3).epsilon(kEps));
        CHECK(u.QN(2, 0) == doctest::Approx(2.2).epsilon(kEps));
        CHECK(u.TN(0, 0) == doctest::Approx(1.65).epsilon(kEps));
        CHECK(u.CN[0] == doctest::Approx(2.8231292517).epsilon(kEps));
        const mva::AvgResult<double> l = run_ba(m, "pb.lower");
        CHECK(l.QN(0, 0) == doctest::Approx(1.41686746988).epsilon(kEps));
        CHECK(l.TN(0, 0) == doctest::Approx(1.41686746988).epsilon(kEps));
        CHECK(l.CN[0] == doctest::Approx(2.42424242424).epsilon(kEps));
    }
    SUBCASE("gb, the geometric bounds, which 'default' resolves to") {
        const mva::AvgResult<double> u = run_ba(m, "gb.upper");
        CHECK(u.QN(0, 0) == doctest::Approx(1.69227238434).epsilon(kEps));
        CHECK(u.QN(1, 0) == doctest::Approx(2.53039101562).epsilon(kEps));
        CHECK(u.QN(2, 0) == doctest::Approx(1.11038125).epsilon(kEps));
        CHECK(u.UN(1, 0) == doctest::Approx(0.846136192168).epsilon(kEps));
        CHECK(u.UN(2, 0) == doctest::Approx(0.564090794779).epsilon(kEps));
        // the response time follows the geometric queue length and the OPPOSITE
        // side of the throughput bracket, not the ABA residence convention
        CHECK(u.RN(0, 0) == doctest::Approx(1.0).epsilon(kEps));
        CHECK(u.RN(1, 0) == doctest::Approx(1.71460049143).epsilon(kEps));
        CHECK(u.RN(2, 0) == doctest::Approx(0.752397643355).epsilon(kEps));
        CHECK(u.TN(0, 0) == doctest::Approx(1.69227238434).epsilon(kEps));
        CHECK(u.CN[0] == doctest::Approx(2.71041191791).epsilon(kEps));
        const mva::AvgResult<double> l = run_ba(m, "gb.lower");
        CHECK(l.QN(0, 0) == doctest::Approx(1.47579044114).epsilon(kEps));
        CHECK(l.QN(1, 0) == doctest::Approx(1.0100735775).epsilon(kEps));
        CHECK(l.QN(2, 0) == doctest::Approx(0.525526995687).epsilon(kEps));
        CHECK(l.UN(1, 0) == doctest::Approx(0.737895220571).epsilon(kEps));
        CHECK(l.RN(1, 0) == doctest::Approx(0.596874112492).epsilon(kEps));
        CHECK(l.RN(2, 0) == doctest::Approx(0.310545158422).epsilon(kEps));
        CHECK(l.TN(0, 0) == doctest::Approx(1.47579044114).epsilon(kEps));
        CHECK(l.CN[0] == doctest::Approx(2.3636856791).epsilon(kEps));
        // 'default' is the alias, and reports the resolved name
        const mva::AvgResult<double> dflt = run_ba(m, "default");
        CHECK(dflt.method == "default");
        CHECK(dflt.actualmethod == "gb.upper");
        CHECK(dflt.TN(0, 0) == doctest::Approx(1.69227238434).epsilon(kEps));
    }
    SUBCASE("mwba, the Majumdar-Woodside robust box bounds") {
        // Q re-baselined 2026-08-01 onto ba_chain_qfill: neither the
        // no-contention residence nor the Theorem-1 one bounds Q on its side.
        const mva::AvgResult<double> u = run_ba(m, "mwba.upper");
        CHECK(u.QN(0, 0) == doctest::Approx(2.0).epsilon(kEps));
        CHECK(u.QN(1, 0) == doctest::Approx(4.0).epsilon(kEps));
        CHECK(u.QN(2, 0) == doctest::Approx(2.66666666667).epsilon(kEps));
        CHECK(u.RN(1, 0) == doctest::Approx(2.0).epsilon(kEps));
        CHECK(u.TN(0, 0) == doctest::Approx(2.0).epsilon(kEps));
        CHECK(u.XN[0] == doctest::Approx(2.0).epsilon(kEps));
        const mva::AvgResult<double> l = run_ba(m, "mwba.lower");
        CHECK(l.QN(0, 0) == doctest::Approx(0.923076923077).epsilon(kEps));
        CHECK(l.QN(1, 0) == doctest::Approx(0.461538461538).epsilon(kEps));
        CHECK(l.QN(2, 0) == doctest::Approx(0.307692307692).epsilon(kEps));
        CHECK(l.RN(1, 0) == doctest::Approx(0.5).epsilon(kEps));
        CHECK(l.TN(0, 0) == doctest::Approx(0.923076923077).epsilon(kEps));
    }
    SUBCASE("the hierarchical and iterative families at level 2") {
        // pbh, pbk and bjbk are the same Eager-Sevcik recursion at level k
        for (const char* fam : {"pbh", "pbk", "bjbk"}) {
            const mva::AvgResult<double> u = run_ba(m, std::string(fam) + ".upper");
            CHECK(u.TN(0, 0) == doctest::Approx(1.57130044843).epsilon(kEps));
            CHECK(u.QN(1, 0) == doctest::Approx(3.14260089686).epsilon(kEps));
            CHECK(u.XN[0] == doctest::Approx(1.57130044843).epsilon(kEps));
            CHECK(u.CN[0] == doctest::Approx(4.33333333333).epsilon(kEps));
            const mva::AvgResult<double> l = run_ba(m, std::string(fam) + ".lower");
            CHECK(l.TN(0, 0) == doctest::Approx(1.5).epsilon(kEps));
            CHECK(l.QN(1, 0) == doctest::Approx(0.75).epsilon(kEps));
            CHECK(l.XN[0] == doctest::Approx(1.5).epsilon(kEps));
            CHECK(l.CN[0] == doctest::Approx(1.83333333333).epsilon(kEps));
        }
        const mva::AvgResult<double> cu = run_ba(m, "cbh.upper");
        CHECK(cu.TN(0, 0) == doctest::Approx(1.55480033985).epsilon(kEps));
        CHECK(cu.QN(1, 0) == doctest::Approx(3.10960067969).epsilon(kEps));
        const mva::AvgResult<double> cl = run_ba(m, "cbh.lower");
        CHECK(cl.TN(0, 0) == doctest::Approx(1.55480033985).epsilon(kEps));
        CHECK(cl.QN(1, 0) == doctest::Approx(0.777400169924).epsilon(kEps));
        // 2026-07-29: ssd now ANSWERS a model with a delay station. The history
        // matters and is kept, because it is what stops the defect returning.
        // Theorem 5 is stated for Z = 0, and a bare +Z in its denominators is
        // NOT a bound: a sweep of 4000 single-server cases found 2825
        // violations, and on this very model it gave 1.2972972973, BELOW the
        // exact 1.55480033985 that an UPPER bound must dominate. Section 5.2
        // Property 4 disaggregates the delay instead, which is valid but loose
        // (TN 24/13 = 1.84615384615, QN(1,0) 48/13 = 3.69230769231 upper, 1.2
        // lower); the port refused the case outright rather than ship either.
        // What replaced both is the Lazowska et al. 1984 Table 5.2
        // terminal-workload correction, which deflates the QUEUEING term rather
        // than adding Z bare: (N-1)Yl/(1+Z/(N Rl)) lower, (N-1)Yu/(1+Z/Ru)
        // upper, the asymmetric divisors being the point. Both reduce to
        // Theorem 5 at Z = 0. MATLAB SolverBA on this model.
        const mva::AvgResult<double> su = run_ba(m, "ssd.upper");
        CHECK(su.TN(0, 0) == doctest::Approx(1.66561514195584).epsilon(kEps));
        CHECK(su.QN(1, 0) == doctest::Approx(3.33123028391167).epsilon(kEps));
        const mva::AvgResult<double> sl = run_ba(m, "ssd.lower");
        CHECK(sl.TN(0, 0) == doctest::Approx(1.33905579399142).epsilon(kEps));
        CHECK(sl.QN(1, 0) == doctest::Approx(0.669527896995708).epsilon(kEps));
        // The property the numbers exist to serve: the bracket must contain the
        // exact throughput. This is what the +Z form violated.
        CHECK(sl.TN(0, 0) <= 1.55480033985);
        CHECK(su.TN(0, 0) >= 1.55480033985);
    }
    SUBCASE("the multiclass composite bounds on a single-class model") {
        const mva::AvgResult<double> u = run_ba(m, "cub.upper");
        CHECK(u.TN(0, 0) == doctest::Approx(2.0).epsilon(kEps));
        CHECK(u.QN(1, 0) == doctest::Approx(4.0).epsilon(kEps));
        const mva::AvgResult<double> l = run_ba(m, "mbjb.lower");
        CHECK(l.TN(0, 0) == doctest::Approx(1.2).epsilon(kEps));
        CHECK(l.QN(1, 0) == doctest::Approx(0.6).epsilon(kEps));
        CHECK(l.CN[0] == doctest::Approx(3.33333333333).epsilon(kEps));
    }
}

TEST_CASE("ba: every bracket contains the exact MVA throughput") {
    qn::Network<double> m = model_a();
    const mva::AvgResult<double> ex = run_mva(m);
    const double Xex = ex.TN(0, 0);
    CHECK(Xex == doctest::Approx(1.55480033985).epsilon(kEps));
    // The bound is on the chain throughput, so the check is per family and does
    // not depend on any expected value transcribed above.
    for (const char* fam : {"aba", "bjb", "pb", "gb", "mwba", "pbh", "pbk", "bjbk", "cbh"}) {
        const double lo = run_ba(m, std::string(fam) + ".lower").TN(0, 0);
        const double hi = run_ba(m, std::string(fam) + ".upper").TN(0, 0);
        INFO("family ", fam, ": [", lo, ", ", hi, "] must contain ", Xex);
        CHECK(lo <= Xex * (1.0 + 1e-12));
        CHECK(hi >= Xex * (1.0 - 1e-12));
    }
    // cub is upper-only and mbjb is the lower bound that seeds it
    CHECK(run_ba(m, "cub.upper").TN(0, 0) >= Xex);
    CHECK(run_ba(m, "mbjb.lower").TN(0, 0) <= Xex);
    // ba_bounds returns both sides of a family in one call, at the caller's level
    const ba::BaOptions opt = [] {
        ba::BaOptions o;
        o.method = "gb.upper";
        return o;
    }();
    const ba::BaBounds<double> b = ba::ba_bounds(m.get_struct(), opt);
    CHECK(b.has_lower);
    CHECK(b.has_upper);
    CHECK(b.Tlower(0, 0) == doctest::Approx(1.47579044114).epsilon(kEps));
    CHECK(b.Tupper(0, 0) == doctest::Approx(1.69227238434).epsilon(kEps));
    CHECK(b.Qupper(1, 0) == doctest::Approx(2.53039101562).epsilon(kEps));
    // a one-sided family leaves its missing side NaN, never zero
    ba::BaOptions oc;
    oc.method = "cub.upper";
    const ba::BaBounds<double> bc = ba::ba_bounds(m.get_struct(), oc);
    CHECK(bc.has_upper);
    CHECK_FALSE(bc.has_lower);
    CHECK(std::isnan(bc.Tlower(0, 0)));
}

TEST_CASE("ba: a closed network with no delay station, where sb and sib apply") {
    qn::Network<double> m = model_b();
    SUBCASE("aba") {
        const mva::AvgResult<double> u = run_ba(m, "aba.upper");
        CHECK(u.QN(0, 0) == doctest::Approx(3.0).epsilon(kEps));
        CHECK(u.QN(1, 0) == doctest::Approx(2.0).epsilon(kEps));
        CHECK(u.QN(2, 0) == doctest::Approx(1.2).epsilon(kEps));
        CHECK(u.UN(0, 0) == doctest::Approx(1.0).epsilon(kEps));
        CHECK(u.RN(0, 0) == doctest::Approx(1.5).epsilon(kEps));
        CHECK(u.TN(0, 0) == doctest::Approx(2.0).epsilon(kEps));
        CHECK(u.CN[0] == doctest::Approx(3.1).epsilon(kEps));
        const mva::AvgResult<double> l = run_ba(m, "aba.lower");
        CHECK(l.TN(0, 0) == doctest::Approx(0.967741935484).epsilon(kEps));
        CHECK(l.CN[0] == doctest::Approx(1.03333333333).epsilon(kEps));
    }
    SUBCASE("sb, the power-sum bounds of Harel-Namn-Sturm") {
        const mva::AvgResult<double> u = run_ba(m, "sb.upper");
        CHECK(u.QN(0, 0) == doctest::Approx(0.820031428335).epsilon(kEps));
        CHECK(u.QN(1, 0) == doctest::Approx(0.54668761889).epsilon(kEps));
        CHECK(u.QN(2, 0) == doctest::Approx(0.328012571334).epsilon(kEps));
        CHECK(u.TN(0, 0) == doctest::Approx(1.64006285667).epsilon(kEps));
        CHECK(u.CN[0] == doctest::Approx(1.82919818457).epsilon(kEps));
        const mva::AvgResult<double> l = run_ba(m, "sb.lower");
        CHECK(l.QN(0, 0) == doctest::Approx(0.813169682569).epsilon(kEps));
        CHECK(l.TN(0, 0) == doctest::Approx(1.62633936514).epsilon(kEps));
        CHECK(l.CN[0] == doctest::Approx(1.8446334537).epsilon(kEps));
        // At N = 3 the power-sum cycle time A1 + 2(A1A2+A3)/(A1^2+A2) is
        // algebraically the exact convolution cycle time h3/h2, so the upper
        // bound coincides with exact MVA. Nothing transcribed here decides that.
        CHECK(u.TN(0, 0) == doctest::Approx(run_mva(m).TN(0, 0)).epsilon(1e-12));
    }
    SUBCASE("gb without a delay station") {
        const mva::AvgResult<double> u = run_ba(m, "gb.upper");
        CHECK(u.QN(0, 0) == doctest::Approx(2.08538334106).epsilon(kEps));
        CHECK(u.QN(1, 0) == doctest::Approx(1.02666611637).epsilon(kEps));
        CHECK(u.QN(2, 0) == doctest::Approx(0.477922768264).epsilon(kEps));
        CHECK(u.RN(0, 0) == doctest::Approx(1.3138930061).epsilon(kEps));
        CHECK(u.RN(2, 0) == doctest::Approx(0.301114605797).epsilon(kEps));
        CHECK(u.TN(0, 0) == doctest::Approx(1.65714824651).epsilon(kEps));
        const mva::AvgResult<double> l = run_ba(m, "gb.lower");
        CHECK(l.QN(0, 0) == doctest::Approx(1.15027928634).epsilon(kEps));
        CHECK(l.RN(0, 0) == doctest::Approx(0.694131794642).epsilon(kEps));
        CHECK(l.TN(0, 0) == doctest::Approx(1.58717896463).epsilon(kEps));
    }
    SUBCASE("sib, the successively improving bounds, which need Z = 0") {
        const mva::AvgResult<double> u = run_ba(m, "sib.upper");
        CHECK(u.TN(0, 0) == doctest::Approx(1.64113908166).epsilon(kEps));
        CHECK(u.QN(0, 0) == doctest::Approx(2.46170862249).epsilon(kEps));
        CHECK(u.CN[0] == doctest::Approx(3.1).epsilon(kEps));
        const mva::AvgResult<double> l = run_ba(m, "sib.lower");
        CHECK(l.TN(0, 0) == doctest::Approx(1.63872072724).epsilon(kEps));
        CHECK(l.QN(0, 0) == doctest::Approx(0.819360363618).epsilon(kEps));
        CHECK(l.CN[0] == doctest::Approx(1.03333333333).epsilon(kEps));
        // a delay station is refused by name rather than dropped
        qn::Network<double> a = model_a();
        CHECK_THROWS_AS(run_ba(a, "sib.upper"), UnsupportedError);
    }
    SUBCASE("harel, the sharp bounds of the same paper as sb") {
        // Not pinned to transcribed constants: the bracket around the exact
        // throughput is the property the family claims, and sb (the power-sum
        // form) is the sibling it must not be confused with.
        const mva::AvgResult<double> u = run_ba(m, "harel.upper");
        const mva::AvgResult<double> l = run_ba(m, "harel.lower");
        const double Xex = run_mva(m).TN(0, 0);
        CHECK(l.TN(0, 0) <= Xex * (1 + 1e-12));
        CHECK(u.TN(0, 0) >= Xex * (1 - 1e-12));
        // at N = 3 the extrapolation point is N itself, so UB is exact
        CHECK(u.TN(0, 0) == doctest::Approx(Xex).epsilon(1e-12));
        CHECK(u.CN[0] == doctest::Approx(3.1).epsilon(kEps));
        CHECK(l.CN[0] == doctest::Approx(1.03333333333).epsilon(kEps));
        // a delay station is refused by name rather than dropped
        qn::Network<double> a = model_a();
        CHECK_THROWS_AS(run_ba(a, "harel.upper"), UnsupportedError);
    }
    SUBCASE("ldbcmp, the closed-open equivalence lower bound") {
        const mva::AvgResult<double> l = run_ba(m, "ldbcmp.lower");
        CHECK(l.QN(0, 0) == doctest::Approx(0.253076586541).epsilon(kEps));
        CHECK(l.TN(0, 0) == doctest::Approx(0.506153173083).epsilon(kEps));
    }
    SUBCASE("the hierarchical families without a delay station") {
        const mva::AvgResult<double> pu = run_ba(m, "pbh.upper");
        CHECK(pu.TN(0, 0) == doctest::Approx(1.6577540107).epsilon(kEps));
        const mva::AvgResult<double> pl = run_ba(m, "pbh.lower");
        CHECK(pl.TN(0, 0) == doctest::Approx(1.5935334873).epsilon(kEps));
        const mva::AvgResult<double> cu = run_ba(m, "cbh.upper");
        CHECK(cu.TN(0, 0) == doctest::Approx(1.64006285667).epsilon(kEps));
        CHECK(cu.QN(0, 0) == doctest::Approx(2.46009428501).epsilon(kEps));
        const mva::AvgResult<double> cl = run_ba(m, "cbh.lower");
        CHECK(cl.TN(0, 0) == doctest::Approx(1.64006285667).epsilon(kEps));
    }
    SUBCASE("mwba without a delay station") {
        const mva::AvgResult<double> u = run_ba(m, "mwba.upper");
        CHECK(u.QN(0, 0) == doctest::Approx(3.0).epsilon(kEps));
        CHECK(u.TN(0, 0) == doctest::Approx(2.0).epsilon(kEps));
        const mva::AvgResult<double> l = run_ba(m, "mwba.lower");
        CHECK(l.QN(0, 0) == doctest::Approx(0.483870967742).epsilon(kEps));
        CHECK(l.TN(0, 0) == doctest::Approx(0.967741935484).epsilon(kEps));
    }
    SUBCASE("the brackets contain the exact answer") {
        const double Xex = run_mva(m).TN(0, 0);
        CHECK(Xex == doctest::Approx(1.64006285667).epsilon(kEps));
        // pbk and bjbk are the same PBH recursion and belong here exactly as
        // they do on the model with a delay; their earlier absence was an
        // asymmetry with no reason behind it, not an exclusion.
        for (const char* fam :
             {"aba", "bjb", "pb", "sb", "gb", "sib", "mwba", "pbh", "pbk", "bjbk", "cbh"}) {
            const double lo = run_ba(m, std::string(fam) + ".lower").TN(0, 0);
            const double hi = run_ba(m, std::string(fam) + ".upper").TN(0, 0);
            INFO("family ", fam, ": [", lo, ", ", hi, "] must contain ", Xex);
            CHECK(lo <= Xex * (1.0 + 1e-12));
            CHECK(hi >= Xex * (1.0 - 1e-12));
        }
        CHECK(run_ba(m, "ldbcmp.lower").TN(0, 0) <= Xex);
    }
}

TEST_CASE("ba: the multiclass families on a two-class closed network") {
    // Q and the R it implies re-baselined 2026-08-01 onto ba_chain_qfill.
    qn::Network<double> m = model_c();
    // station order Delay, Queue1, Queue2; class order C1, C2
    SUBCASE("mwba.upper") {
        const mva::AvgResult<double> r = run_ba(m, "mwba.upper");
        CHECK(r.QN(0, 0) == doctest::Approx(1.37931034483).epsilon(kEps));
        CHECK(r.QN(0, 1) == doctest::Approx(1.26923076923).epsilon(kEps));
        CHECK(r.QN(1, 0) == doctest::Approx(2.0).epsilon(kEps));
        CHECK(r.QN(1, 1) == doctest::Approx(3.0).epsilon(kEps));
        CHECK(r.QN(2, 0) == doctest::Approx(1.39787798408).epsilon(kEps));
        CHECK(r.QN(2, 1) == doctest::Approx(2.09681697613).epsilon(kEps));
        CHECK(r.RN(1, 0) == doctest::Approx(1.45).epsilon(kEps));
        CHECK(r.RN(1, 1) == doctest::Approx(1.18181818182).epsilon(kEps));
        CHECK(r.TN(1, 0) == doctest::Approx(1.37931034483).epsilon(kEps));
        CHECK(r.TN(1, 1) == doctest::Approx(2.53846153846).epsilon(kEps));
        CHECK(r.CN[0] == doctest::Approx(1.45).epsilon(kEps));
        CHECK(r.CN[1] == doctest::Approx(1.18181818182).epsilon(kEps));
        CHECK(r.XN[0] == doctest::Approx(1.37931034483).epsilon(kEps));
        CHECK(r.XN[1] == doctest::Approx(2.53846153846).epsilon(kEps));
    }
    SUBCASE("mwba.lower, whose residence is the Theorem 1 guarantee") {
        const mva::AvgResult<double> r = run_ba(m, "mwba.lower");
        CHECK(r.QN(0, 0) == doctest::Approx(0.615384615385).epsilon(kEps));
        CHECK(r.QN(0, 1) == doctest::Approx(0.5).epsilon(kEps));
        CHECK(r.QN(1, 0) == doctest::Approx(0.153846153846).epsilon(kEps));
        CHECK(r.QN(1, 1) == doctest::Approx(0.333333333333).epsilon(kEps));
        CHECK(r.QN(2, 0) == doctest::Approx(0.123076923077).epsilon(kEps));
        CHECK(r.QN(2, 1) == doctest::Approx(0.166666666667).epsilon(kEps));
        CHECK(r.RN(1, 0) == doctest::Approx(0.25).epsilon(kEps));
        CHECK(r.RN(1, 1) == doctest::Approx(0.333333333333).epsilon(kEps));
        CHECK(r.TN(0, 0) == doctest::Approx(0.615384615385).epsilon(kEps));
        CHECK(r.TN(0, 1) == doctest::Approx(1.0).epsilon(kEps));
        CHECK(r.CN[0] == doctest::Approx(3.25).epsilon(kEps));
        CHECK(r.CN[1] == doctest::Approx(3.0).epsilon(kEps));
    }
    SUBCASE("cub.upper and the mbjb.lower that seeds it") {
        const mva::AvgResult<double> u = run_ba(m, "cub.upper");
        CHECK(u.QN(0, 0) == doctest::Approx(2.28571428571).epsilon(kEps));
        CHECK(u.QN(0, 1) == doctest::Approx(1.19387755102).epsilon(kEps));
        CHECK(u.QN(1, 0) == doctest::Approx(2.0).epsilon(kEps));
        CHECK(u.QN(1, 1) == doctest::Approx(3.0).epsilon(kEps));
        CHECK(u.QN(2, 1) == doctest::Approx(2.56530612245).epsilon(kEps));
        CHECK(u.XN[0] == doctest::Approx(2.28571428571).epsilon(kEps));
        CHECK(u.XN[1] == doctest::Approx(2.38775510204).epsilon(kEps));
        const mva::AvgResult<double> l = run_ba(m, "mbjb.lower");
        CHECK(l.QN(0, 0) == doctest::Approx(0.816326530612).epsilon(kEps));
        CHECK(l.QN(1, 1) == doctest::Approx(0.428571428571).epsilon(kEps));
        CHECK(l.XN[0] == doctest::Approx(0.816326530612).epsilon(kEps));
        CHECK(l.XN[1] == doctest::Approx(1.28571428571).epsilon(kEps));
        CHECK(l.CN[0] == doctest::Approx(2.45).epsilon(kEps));
    }
    SUBCASE("both multiclass brackets contain the exact per-class throughput") {
        const mva::AvgResult<double> ex = run_mva(m);
        CHECK(ex.TN(1, 0) == doctest::Approx(0.98737361521).epsilon(kEps));
        CHECK(ex.TN(1, 1) == doctest::Approx(1.86355280344).epsilon(kEps));
        const mva::AvgResult<double> mu = run_ba(m, "mwba.upper");
        const mva::AvgResult<double> ml = run_ba(m, "mwba.lower");
        const mva::AvgResult<double> cu = run_ba(m, "cub.upper");
        const mva::AvgResult<double> cl = run_ba(m, "mbjb.lower");
        for (std::size_t k = 0; k < 2; ++k) {
            CHECK(ml.XN[k] <= ex.TN(1, k));
            CHECK(mu.XN[k] >= ex.TN(1, k));
            CHECK(cl.XN[k] <= ex.TN(1, k));
            CHECK(cu.XN[k] >= ex.TN(1, k));
        }
    }
    SUBCASE("a single-class family REFUSES a multiclass model") {
        // The 2026-07-25 zero-table convention was alignment with a
        // reference that has since changed: MATLAB, the JAR and python
        // all refuse by name now. Overturned by the user 2026-09-04.
        // These branches used to return zeroed metrics with an EMPTY C and X,
        // which a caller could not tell from a real answer of zero.
        for (const char* mth : {"aba.upper", "aba.lower", "bjb.upper", "pb.lower", "sb.upper",
                                "gb.upper", "default"})
            CHECK_THROWS_AS(run_ba(m, mth), UnsupportedError);
        // mwba IS multiclass, so this model is inside its domain and it answers
        CHECK(run_ba(m, "mwba.upper").XN.size() == 2);
        // the families whose gate is the reference's own line_error still refuse
        for (const char* mth : {"pbh.upper", "cbh.lower", "ssd.upper", "lr.upper"})
            CHECK_THROWS_AS(run_ba(m, mth), UnsupportedError);
    }
    SUBCASE("mwba REFUSES a MIXED model, by the same rule") {
        // mwba's premise is `nclosedjobs>0 && ~any(isinf(njobs))`, so an open
        // class present alongside a closed one puts the model outside its
        // domain. It used to answer with zeros here under the same 2026-07-25
        // convention the subcase above records; that ruling was overturned by
        // the user on 2026-09-04 and every family refuses by name now.
        qn::Network<double> mx("mx");
        const std::size_t src = mx.add_source("Src");
        const std::size_t q = mx.add_queue("Q1", SchedStrategy::PS);
        const std::size_t snk = mx.add_sink("Snk");
        const std::size_t kc = mx.add_closed_class("C1", 2.0, q);
        const std::size_t ko = mx.add_open_class("O1");
        mx.set_service(q, kc, D::exp_rate(4.0));
        mx.set_service(q, ko, D::exp_rate(5.0));
        mx.set_arrival(src, ko, D::exp_rate(0.5));
        qn::RoutingMatrix<double> P;
        P.set(kc, kc, q, q, 1.0);
        P.set(ko, ko, src, q, 1.0);
        P.set(ko, ko, q, snk, 1.0);
        mx.link(P);
        for (const char* mth : {"mwba.upper", "mwba.lower", "aba.upper", "gb.upper"})
            CHECK_THROWS_AS(run_ba(mx, mth), UnsupportedError);
        // cub's fully-closed gate IS a line_error in the reference, so it refuses
        CHECK_THROWS_AS(run_ba(mx, "cub.upper"), UnsupportedError);
        CHECK_THROWS_AS(run_ba(mx, "pbh.upper"), UnsupportedError);
    }
}

TEST_CASE("ba: a multiserver station with a delay is bounded by ssd alone") {
    qn::Network<double> m = model_d();
    // ssd is the ONLY family accepting a multiserver station, and as of
    // 2026-07-29 it accepts a delay too, so this model went from having NO
    // bound at all to having exactly one. The coverage gap recorded here
    // previously is closed. History, because it is what stops the defect
    // returning: the multiserver disaggregation (Theorem 5) and the think-time
    // treatment come from different results in the paper, and the two ways of
    // joining them that were tried both failed. Section 5.2 Property 4
    // disaggregates the delay into a (K+1)-st station: valid but loose, giving
    // upper QN(0,0) 24/11 = 2.18181818182, bound by the ABA population term
    // rather than by Theorem 5, and lower 1.41176470588. A bare +Z in the
    // Theorem 5 denominators is not a bound at all and gave 1.62711864407.
    // The Lazowska et al. 1984 Table 5.2 terminal-workload correction deflates
    // the queueing term instead and is both sound and tighter: the upper below
    // is 2.09896432681243, which beats the 24/11 = 2.18181818182 ABA
    // population term that used to bind, so Theorem 5 now sets the minimum.
    // MATLAB SolverBA on this model.
    const mva::AvgResult<double> su = run_ba(m, "ssd.upper");
    CHECK(su.TN(0, 0) == doctest::Approx(2.09896432681243).epsilon(kEps));
    CHECK(su.QN(0, 0) == doctest::Approx(2.09896432681243).epsilon(kEps));
    CHECK(su.QN(1, 0) == doctest::Approx(4.19792865362486).epsilon(kEps));
    const mva::AvgResult<double> sl = run_ba(m, "ssd.lower");
    CHECK(sl.TN(0, 0) == doctest::Approx(1.53694581280788).epsilon(kEps));
    CHECK(sl.QN(0, 0) == doctest::Approx(1.53694581280788).epsilon(kEps));
    CHECK(sl.QN(1, 0) == doctest::Approx(0.768472906403941).epsilon(kEps));
    // The bracket must contain the exact throughput, which SolverMVA puts at
    // 1.88655395797124 on this model. That is the property the numbers serve.
    CHECK(sl.TN(0, 0) <= 1.88655395797124);
    CHECK(su.TN(0, 0) >= 1.88655395797124);
    // Strictly inside the ABA population term, i.e. Theorem 5 is what binds.
    CHECK(su.TN(0, 0) < 24.0 / 11.0);
    // every family that gates on the server count refuses it and says which to use
    // cub and mbjb rejoined this list on 2026-08-01: the reference gates them on
    // the server count too (solver_ba_analyzer.m:610), and the demand-only read
    // that used to answer here made auto.upper prefer their bound over ssd's.
    for (const char* mth : {"aba.upper", "bjb.lower", "gb.upper", "pbh.upper", "cbh.lower",
                            "sib.upper", "mwba.upper", "cub.upper", "mbjb.lower"})
        CHECK_THROWS_AS(run_ba(m, mth), UnsupportedError);
    // ldbcmp does not gate on it either, and this model sits exactly at its
    // regime boundary N == Qhat = 4, where the bound degenerates to X >= 0
    CHECK(run_ba(m, "ldbcmp.lower").TN(0, 0) == doctest::Approx(0.0));
}

/** Model B over an arbitrary arithmetic: no delay station, so sb applies. */
template <class T>
mva::AvgResult<T> run_ba_arith_b(const std::string& method) {
    qn::Network<T> m("cqnB");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
    const std::size_t q3 = m.add_queue("Queue3", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, q1);
    m.set_service(q1, c, Distrib<T>::exp_rate(num_traits<T>::from_int(2)));
    m.set_service(q2, c, Distrib<T>::exp_rate(num_traits<T>::from_int(3)));
    m.set_service(q3, c, Distrib<T>::exp_rate(num_traits<T>::from_int(5)));
    qn::RoutingMatrix<T> P;
    P.set(c, c, q1, q2, num_traits<T>::from_int(1));
    P.set(c, c, q2, q3, num_traits<T>::from_int(1));
    P.set(c, c, q3, q1, num_traits<T>::from_int(1));
    m.link(P);
    ba::BaOptions opt;
    opt.method = method;
    return ba::solver_ba_run_analyzer(m.get_struct(), opt);
}

/** Model A over an arbitrary arithmetic, for the exact-arithmetic gate. */
template <class T>
mva::AvgResult<T> run_ba_arith(const std::string& method) {
    qn::Network<T> m("cqnA");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 4.0, d);
    m.set_service(d, c, Distrib<T>::exp_rate(num_traits<T>::from_int(1)));
    m.set_service(q1, c, Distrib<T>::exp_rate(num_traits<T>::from_int(2)));
    m.set_service(q2, c, Distrib<T>::exp_rate(num_traits<T>::from_int(3)));
    qn::RoutingMatrix<T> P;
    P.set(c, c, d, q1, num_traits<T>::from_int(1));
    P.set(c, c, q1, q2, num_traits<T>::from_int(1));
    P.set(c, c, q2, d, num_traits<T>::from_int(1));
    m.link(P);
    ba::BaOptions opt;
    opt.method = method;
    return ba::solver_ba_run_analyzer(m.get_struct(), opt);
}

TEST_CASE("ba: the field-arithmetic families are exact, the rest refuse by name") {
    // ABA is sums, products, minima and divisions, so under Rational it is the
    // bound with no rounding at all: X = min(1/Dmax, N/(Z+sum D)) = 2 exactly.
    const mva::AvgResult<Rational> a = run_ba_arith<Rational>("aba.upper");
    CHECK(a.TN(0, 0) == Rational(2));
    CHECK(a.QN(1, 0) == Rational(4));
    // ABA lower: X = 4/(1 + 4*5/6) = 12/13
    const mva::AvgResult<Rational> l = run_ba_arith<Rational>("aba.lower");
    CHECK(l.TN(0, 0) == Rational(12) / Rational(13));
    for (const char* mth : {"bjb.upper", "pb.lower", "mwba.upper", "pbh.upper", "cbh.lower",
                            "cub.upper", "mbjb.lower"})
        CHECK_NOTHROW(run_ba_arith<Rational>(mth));
    // ssd rejoined this list on 2026-07-29, when it regained the think time.
    // The Lazowska Table 5.2 correction is +, * and / only, so cqnA's delay
    // costs it nothing in exactness and the bound is a RATIONAL, not a rounded
    // double lifted back into one. By hand, with L = [1/2, 1/3], N = 4, Z = 1:
    //   Rl = 5/6, Yl = 1/2, Ru = 5/6, Yu = 5/12
    //   Xhi = 4 / (5/6 + 1 + 3(5/12)/(1 + 6/5)) = 4/(317/132) = 528/317
    //   Xlo = 4 / (5/6 + 1 + 3(1/2)/(1 + 3/10)) = 4/(233/78)  = 312/233
    // and the ABA terms 1/(1/2) = 2 and 24/11 both exceed 528/317, so the
    // Theorem 5 term is what the min selects. MATLAB agrees to 15 digits:
    // 1.66561514195584 and 1.33905579399142.
    const mva::AvgResult<Rational> su = run_ba_arith<Rational>("ssd.upper");
    CHECK(su.TN(0, 0) == Rational(528) / Rational(317));
    const mva::AvgResult<Rational> sl = run_ba_arith<Rational>("ssd.lower");
    CHECK(sl.TN(0, 0) == Rational(312) / Rational(233));
    // Exact arithmetic makes the bracket a THEOREM here, not a tolerance: no
    // rounding can be blamed if it ever fails.
    CHECK(sl.TN(0, 0) < su.TN(0, 0));
    CHECK_NOTHROW(run_ba_arith_b<Rational>("ssd.upper"));
    CHECK_NOTHROW(run_ba_arith_b<Rational>("ssd.lower"));
    // the geometric family and sib solve a quadratic; sb.lower takes a root
    CHECK_THROWS_AS(run_ba_arith<Rational>("gb.upper"), UnsupportedError);
    CHECK_THROWS_AS(run_ba_arith<Rational>("gb.lower"), UnsupportedError);
    CHECK_THROWS_AS(run_ba_arith<Rational>("default"), UnsupportedError);
    // sb rejects the delay station before the arithmetic, so the (N-1)-st root
    // gate has to be checked on the model without one; sb.upper is field
    // arithmetic there and stays available
    CHECK_THROWS_AS(run_ba_arith<Rational>("sb.lower"), UnsupportedError);
    CHECK_THROWS_AS(run_ba_arith_b<Rational>("sb.lower"), UnsupportedError);
    CHECK_THROWS_AS(run_ba_arith_b<Rational>("sib.upper"), UnsupportedError);
    // A1 = 31/30, A2 = 401/900, A3 = 4591/27000; the exact HNS cycle time is
    // A1 + 2(A1 A2 + A3)/(A1^2 + A2), which at N = 3 equals the exact one
    const mva::AvgResult<Rational> sbu = run_ba_arith_b<Rational>("sb.upper");
    CHECK(num_traits<Rational>::to_double(sbu.TN(0, 0)) ==
          doctest::Approx(1.64006285667).epsilon(kEps));
    // and high precision behaves as double does
    const mva::AvgResult<Real<30>> g = run_ba_arith<Real<30>>("gb.upper");
    CHECK(num_traits<Real<30>>::to_double(g.TN(0, 0)) ==
          doctest::Approx(1.69227238434).epsilon(kEps));
}

TEST_CASE("ba: the AUTO composite keeps the tightest side") {
    // MATLAB SolverBA(model_a,'auto.upper'/'auto.lower').getAvg, 2026-08-01.
    // The composite probes each family by execution, so a candidate that
    // refuses the model is skipped rather than being predicted.
    qn::Network<double> m = model_a();
    const mva::AvgResult<double> u = run_ba(m, "auto.upper");
    CHECK(u.TN(1, 0) == doctest::Approx(1.65).epsilon(kEps));
    CHECK(u.QN(1, 0) == doctest::Approx(3.3).epsilon(kEps));
    CHECK(u.QN(2, 0) == doctest::Approx(2.2).epsilon(kEps));
    CHECK(u.UN(1, 0) == doctest::Approx(0.825).epsilon(kEps));
    const mva::AvgResult<double> l = run_ba(m, "auto.lower");
    CHECK(l.TN(1, 0) == doctest::Approx(1.47579044114).epsilon(kEps));
    CHECK(l.QN(1, 0) == doctest::Approx(0.737895220571).epsilon(kEps));
    CHECK(l.QN(2, 0) == doctest::Approx(0.491930147047).epsilon(kEps));
    // it is a bracket, and bare 'auto' is its upper side
    CHECK(l.TN(1, 0) < u.TN(1, 0));
    CHECK(run_ba(m, "auto").TN(1, 0) == doctest::Approx(u.TN(1, 0)).epsilon(kEps));
    // no family under the composite reports a utilization above one, which the
    // missing nservers divisor used to do at a multiserver station
    qn::Network<double> ms = model_d();
    for (const char* mth : {"ssd.upper", "ssd.lower", "auto.upper"}) {
        const mva::AvgResult<double> r = run_ba(ms, mth);
        CHECK(r.UN(1, 0) <= 1.0);
    }
}

TEST_CASE("ba: the no-blocking QRF arms reach the optimizer through SolverBA") {
    // Two FCFS queues and no delay, which is what qrf_noblo_* admits: the
    // formulation has no infinite-server notion and models every station as one
    // server, so a delay or a c>1 station is refused by name.
    qn::Network<double> m("qrfA");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, q1);
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::erlang(6.0, 2));
    qn::RoutingMatrix<double> P;
    P.set(q1, q2, 1.0);
    P.set(q2, q1, 1.0);
    m.link(P);
    // The MMI objective is NONCONVEX and every codebase lands on the optimum its
    // start point reaches, so what is pinned here is the population constraint
    // the normalization enforces and the utilizations the ONE constraint caps,
    // not a number. The api-level agreement with the reference is
    // test_mapqn_qrf_noblo.cpp's business.
    for (const char* mth : {"qrf.mmi", "qrf.mem", "qrf.bethe", "qrf.mmi.ld",
                            "qrf.mmi.linear", "qr"}) {
        const mva::AvgResult<double> r = run_ba(m, mth);
        CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(2.0).epsilon(1e-6));
        CHECK(r.UN(0, 0) <= 1.0 + kEps);
        CHECK(r.UN(1, 0) <= 1.0 + kEps);
    }
    // a delay station and a multiserver station are refused, not approximated
    qn::Network<double> withdelay = model_a();
    CHECK_THROWS_AS(run_ba(withdelay, "qrf.mmi"), UnsupportedError);
    // qrf.bas without a hand-supplied `qrf_params` used to be REFUSED, because
    // the only alternative then was to INVENT the tables and a no-blocking
    // fallback sits ~31x farther from exact than a parameterised run. It is not
    // invented any more: `sn_to_qrf_blocking` DERIVES it from the model, and on
    // this unblocked model that derivation is exactly the trivial table (one
    // configuration, nothing blocked, F = N) a caller used to have to write out.
    // So the call succeeds, and what is checked is that it answers THIS model.
    const mva::AvgResult<double> bas = run_ba(m, "qrf.bas");
    CHECK(bas.QN(0, 0) + bas.QN(1, 0) == doctest::Approx(2.0).epsilon(1e-6));
    CHECK(bas.UN(0, 0) <= 1.0 + kEps);
    CHECK(bas.UN(1, 0) <= 1.0 + kEps);
}

TEST_CASE("ba: the LP blocking bounds solve once qrf_params is supplied") {
    // The same two-queue model, with Queue2 the finite-capacity station and one
    // blocking configuration. The tables follow the reference's 1-based queue
    // indices, which the adapter shifts to the port's 0-based ones.
    qn::Network<double> m("qrfB");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, q1);
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(q1, q2, 1.0);
    P.set(q2, q1, 1.0);
    m.link(P);

    ba::BaOptions opt;
    opt.method = "qrf.bas";
    opt.qrf_params.supplied = true;
    opt.qrf_params.f = 2;   // Queue2, 1-based as in sn_to_qrf_params
    opt.qrf_params.MR = 1;
    opt.qrf_params.BB.push_back(std::vector<int>(2, 0));
    std::vector<int> mm;
    mm.push_back(2);
    mm.push_back(2);
    opt.qrf_params.MM.push_back(mm);
    opt.qrf_params.MM1.push_back(std::vector<int>(2, 0));
    opt.qrf_params.ZZ.push_back(0);

    const mva::AvgResult<double> r = ba::solver_ba_run_analyzer(m.get_struct(), opt);
    // A BOUND, so what is pinned is the population constraint the normalization
    // enforces and the ONE constraint's cap on each utilization; the vertex
    // itself is test_mapqn_qr_bounds_bas.cpp's business.
    CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(2.0).epsilon(1e-6));
    CHECK(r.UN(0, 0) <= 1.0 + kEps);
    CHECK(r.UN(1, 0) <= 1.0 + kEps);
    CHECK(r.UN(0, 0) >= -kEps);
    CHECK(r.UN(1, 0) >= -kEps);
}

TEST_CASE("ba: qrf.bas.mem reaches the BAS entropy optimizer") {
    qn::Network<double> m("qrfBmem");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, q1);
    m.set_service(q1, c, D::exp_rate(1.0));
    m.set_service(q2, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(q1, q2, 1.0);
    P.set(q2, q1, 1.0);
    m.link(P);

    ba::BaOptions opt;
    opt.method = "qrf.bas.mem";
    // The blocking tables used to be DEMANDED here, and `F` with them: without
    // `qrf_params` this threw UnsupportedError, and with the tables but no `F`
    // it threw InputError. Both are derived now -- `sn_to_qrf_blocking` for the
    // tables, `sn_to_qrf_capacity` for `F` -- so neither refusal is reachable
    // and a bare call runs. On this UNBLOCKED model the derivation is exactly
    // the trivial table the block below spells out by hand (one configuration,
    // nothing blocked, F = N), which is what makes the agreement checked at the
    // end a real assertion rather than a restatement.
    const mva::AvgResult<double> derived = ba::solver_ba_run_analyzer(m.get_struct(), opt);

    opt.qrf_params.supplied = true;
    opt.qrf_params.f = 2;
    opt.qrf_params.MR = 1;
    opt.qrf_params.BB.push_back(std::vector<int>(2, 0));
    opt.qrf_params.MM.push_back(std::vector<int>(2, 0));
    opt.qrf_params.MM1.push_back(std::vector<int>(2, 0));
    opt.qrf_params.ZZ.push_back(0);
    opt.qrf_params.F.push_back(2);
    opt.qrf_params.F.push_back(2);
    const mva::AvgResult<double> r = ba::solver_ba_run_analyzer(m.get_struct(), opt);
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(derived.UN(i, 0) == doctest::Approx(r.UN(i, 0)).epsilon(1e-9));
        CHECK(derived.QN(i, 0) == doctest::Approx(r.QN(i, 0)).epsilon(1e-9));
    }
    CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(2.0).epsilon(1e-6));
    CHECK(r.UN(0, 0) == doctest::Approx(6.0 / 7.0).epsilon(1e-6));
    CHECK(r.UN(1, 0) == doctest::Approx(3.0 / 7.0).epsilon(1e-6));
    CHECK(r.QN(0, 0) == doctest::Approx(10.0 / 7.0).epsilon(1e-6));
    CHECK(r.QN(1, 0) == doctest::Approx(4.0 / 7.0).epsilon(1e-6));
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(std::isfinite(r.QN(i, 0)));
        CHECK(std::isfinite(r.UN(i, 0)));
        CHECK(r.QN(i, 0) >= -kEps);
        CHECK(r.UN(i, 0) >= -kEps);
        CHECK(r.UN(i, 0) <= 1.0 + kEps);
    }
    CHECK(r.method == "qrf.bas.mem");
    CHECK(r.actualmethod == "qrf.bas.mem");
}

TEST_CASE("ba: scb brackets the MULTICLASS system, not this model") {
    // Model B has no delay station, which scb requires: Theorem 3 of Dowdy et
    // al. (1992) rests on the delay-free balanced-network throughput.
    qn::Network<double> m = model_b();
    const mva::AvgResult<double> lo = run_ba(m, "scb.lower");
    const mva::AvgResult<double> hi = run_ba(m, "scb.upper");
    const mva::AvgResult<double> ex = run_mva(m);

    // UNLIKE EVERY OTHER FAMILY HERE, the lower side is not a bound on this
    // model at all -- it IS this model's exact throughput. What the pair
    // brackets is the multiclass system the single-class model aggregates.
    CHECK(lo.XN[0] == doctest::Approx(ex.XN[0]).epsilon(1e-9));
    CHECK(hi.XN[0] >= lo.XN[0]);
    // Theorem 3 gap at N = 3, K = 3: (3-1)/(3+3-1) = 2/5, capped by U <= 1.
    const double ratio = hi.XN[0] / lo.XN[0];
    CHECK(ratio <= 1.0 + pfqn::pfqn_scbgap<double>(3, 3) / (1.0 - pfqn::pfqn_scbgap<double>(3, 3)) + kEps);
    // Corollary 1: the utilization ratio is the same at every device, and the
    // busiest device is what the single-server cap binds on.
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("station ", i);
        CHECK(hi.UN(i, 0) / lo.UN(i, 0) == doctest::Approx(ratio).epsilon(1e-9));
        CHECK(hi.UN(i, 0) <= 1.0 + kEps);
    }

    // A delay station is refused by name, as it is in MATLAB.
    qn::Network<double> a = model_a();
    CHECK_THROWS_AS(run_ba(a, "scb.lower"), UnsupportedError);
    CHECK_THROWS_AS(run_ba(a, "scb.upper"), UnsupportedError);
    // A multiserver station is refused, like every family but ssd and ldbcmp.
    qn::Network<double> d = model_d();
    CHECK_THROWS_AS(run_ba(d, "scb.upper"), UnsupportedError);
    // Multiclass is refused: the aggregate is the INPUT to scb, not its output.
    qn::Network<double> c = model_c();
    CHECK_THROWS_AS(run_ba(c, "scb.upper"), UnsupportedError);

    // scb must NOT be an auto candidate: auto brackets this model's own
    // solution, and scb.lower equals it exactly, which would collapse the
    // composite lower bound onto the exact answer and misreport it as a bound.
    const std::vector<std::string> valid = ba::list_valid_methods();
    CHECK(std::find(valid.begin(), valid.end(), "scb.upper") != valid.end());
    CHECK(std::find(valid.begin(), valid.end(), "scb.lower") != valid.end());
    const mva::AvgResult<double> auto_lo = run_ba(m, "auto.lower");
    CHECK(auto_lo.XN[0] < ex.XN[0]);
}

TEST_CASE("ba: the model-aware list drops what the model cannot run") {
    // model_a carries a DELAY, which the qrf_noblo_* formulation has no notion
    // of, so the reduction bounds are refused on it. The model-independent list
    // still names them -- a direct request must reach the analyzer and get its
    // own reason -- while the model-aware overload does not, so a caller
    // enumerating the list never asks for one. The two LOAD-DEPENDENT arms read
    // the model through alpha(i,n) instead and are the exception, checked below.
    qn::Network<double> m = model_a();
    const std::vector<std::string> all = ba::list_valid_methods();
    const std::vector<std::string> fit = ba::list_valid_methods(m.get_struct());
    auto in = [](const std::vector<std::string>& v, const std::string& x) {
        return std::find(v.begin(), v.end(), x) != v.end();
    };
    for (const char* name : {"qr", "lr", "lr.upper", "lr.lower", "qrf.mmi",
                             "qrf.mem", "qrf.bethe", "qrf.bas", "qrf.bas.mmi",
                             "qrf.bas.mem", "qrf.bas.bethe", "qrf.rsrd"}) {
        CHECK(in(all, name));
        CHECK_FALSE(in(fit, name));
    }
    // 'qrf.mmi.ld' and 'qrf.mmi.linear' are the exception: alpha(i,n) = n IS a
    // delay, so they are the only two reduction arms model_a can run, and
    // dropping them would leave a caller enumerating the list with none.
    CHECK(in(fit, "qrf.mmi.ld"));
    CHECK(in(fit, "qrf.mmi.linear"));
    // 'spnlp.*' is indexed by a MARKING rather than by demands and a
    // population, so it drops on every model that is not a Petri net, exactly
    // as the demand-parameterized families drop on one.
    for (const char* name : {"spnlp.upper", "spnlp.lower", "spnlp.op.upper",
                             "spnlp.op.lower"}) {
        CHECK(in(all, name));
        CHECK_FALSE(in(fit, name));
    }
    // 'bpt.lower', 'bgt.upper' and 'snc.upper' are the mirror image: all three
    // are derived for an OPEN network, so every closed model drops them whatever
    // else it carries. 'snc.upper' joined the list after this case was written,
    // which is what made the counts below one short.
    for (const char* name : {"bpt.lower", "bgt.upper", "snc.upper"}) {
        CHECK(in(all, name));
        CHECK_FALSE(in(fit, name));
    }
    // Everything else survives: the narrowing is the reduction family less its
    // two LD arms, plus the three open-network bounds and the four Petri-net
    // ones, and nothing more.
    CHECK(in(fit, "gb.upper"));
    // ldbcmp.lower is NOT offered here: this model sits exactly at
    // N == Qhat = 4, where the bound degenerates to X >= 0 and the
    // table is all zeros. Asking by NAME still runs and still
    // publishes it; see ba::method_degenerate.
    CHECK_FALSE(in(fit, "ldbcmp.lower"));
    CHECK(in(fit, "scb.lower"));
    // This constant is the COUNT OF NAMES THIS MODEL DROPS, so it moves whenever
    // the model-independent list gains a member THIS model cannot run: 14
    // before 'qrf.bas.mmi' joined it, 15 after, 16 with 'qrf.bethe', then 21
    // once 'qrf.bas.bethe' and the four 'spnlp.*' arrived, 19 since
    // 'qrf.mmi.ld' and 'qrf.mmi.linear' started surviving on a delay model, and
    // 20 since 'ldbcmp.lower' is withheld at the N == Qhat boundary this model
    // sits on, and 22 since a79d75979 added the two 'mapamva.*' bounds, which
    // need a modulated process this model does not carry. It pins no numeric
    // answer.
    CHECK(fit.size() + 22 == all.size());

    // A single-class closed model of single-server stations keeps every bound
    // but the three open-network ones and the four Petri-net ones.
    qn::Network<double> nd("two-queue");
    const std::size_t q1 = nd.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = nd.add_queue("Q2", SchedStrategy::PS);
    const std::size_t cc = nd.add_closed_class("C1", 3, q1);
    nd.set_service(q1, cc, D::exp_rate(1.0));
    nd.set_service(q2, cc, D::exp_rate(2.0));
    qn::RoutingMatrix<double> Pn;
    Pn.set(q1, q2, 1.0);
    Pn.set(q2, q1, 1.0);
    nd.link(Pn);
    CHECK(ba::list_valid_methods(nd.get_struct()).size() + 7 == all.size());
}

TEST_CASE("ba: the method gate refuses by name") {
    qn::Network<double> m = model_a();
    CHECK_THROWS_AS(run_ba(m, "nosuchbound"), UnsupportedError);
    // the no-blocking QRF arms ARE served since 2026-08-01; they refuse on THIS
    // model because model_a has a delay station, which the qrf_noblo_*
    // formulation has no notion of, and that is the reference's own gate
    CHECK_THROWS_AS(run_ba(m, "qr"), UnsupportedError);
    CHECK_THROWS_AS(run_ba(m, "qrf.mmi"), UnsupportedError);
    // all three BAS NLP arms are served since 2026-09-03; they refuse on THIS
    // model because model_a has a delay station, the reference's own gate
    CHECK_THROWS_AS(run_ba(m, "qrf.bas.mem"), UnsupportedError);
    CHECK_THROWS_AS(run_ba(m, "qrf.bas.mmi"), UnsupportedError);
    CHECK_THROWS_AS(run_ba(m, "qrf.rsrd"), UnsupportedError);
    // lr IS ported; it refuses on THIS model only because model_a has a delay
    // station, which is the reference's own gate and not a port boundary
    CHECK_THROWS_AS(run_ba(m, "lr"), UnsupportedError);
    CHECK_THROWS_AS(run_ba(m, "lr.upper"), UnsupportedError);
    // the bound family is exactly what SolverMVA refuses and points here
    mva::MvaOptions o;
    o.method = "gb.upper";
    CHECK_THROWS_AS(mva::solver_mva_run_analyzer(m.get_struct(), o, Matrix<double>()), UnsupportedError);
    // a purely OPEN model is refused at the CLASS level (runAnalyzer's
    // nclosedjobs gate, a line_error in the reference) before any family sees
    // it -- unlike a MIXED model, which reaches the family and gets zeros
    qn::Network<double> open("mm1");
    const std::size_t s = open.add_source("Source");
    const std::size_t q = open.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t k = open.add_sink("Sink");
    const std::size_t oc = open.add_open_class("O1");
    open.set_arrival(s, oc, D::exp_rate(1.0));
    open.set_service(q, oc, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    open.link(P);
    CHECK_THROWS_AS(run_ba(open, "aba.upper"), UnsupportedError);
    // the whitelist advertises only what runs
    const std::vector<std::string> valid = ba::list_valid_methods();
    auto has = [&](const std::string& x) {
        return std::find(valid.begin(), valid.end(), x) != valid.end();
    };
    CHECK(has("gb.upper"));
    CHECK(has("ldbcmp.lower"));
    CHECK(has("qr"));
    CHECK(has("qrf.mmi"));
    CHECK(has("qrf.mem"));
    // the two LP blocking bounds are served since 2026-08-02; they read their
    // tables from options.config.qrf_params and refuse without them
    CHECK(has("qrf.bas"));
    CHECK(has("qrf.bas.mem"));
    CHECK(has("qrf.rsrd"));
    CHECK_FALSE(has("ldbcmp.upper"));
    CHECK(has("auto.upper"));
    CHECK(has("auto.lower"));
    // lr IS served: its bound is a pure LP, solved by the exact simplex
    CHECK(has("lr"));
    CHECK(has("lr.upper"));
    CHECK(has("lr.lower"));
}

/** Two FCFS queues, mu = 1 and 2, N = 2: small enough to reason about by hand. */
template <class T>
qn::Network<T> model_lr_small() {
    qn::Network<T> m("cqn2");
    const std::size_t a1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t a2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t k = m.add_closed_class("C1", 2.0, a1);
    m.set_service(a1, k, Distrib<T>::exp_rate(num_traits<T>::from_int(1)));
    m.set_service(a2, k, Distrib<T>::exp_rate(num_traits<T>::from_int(2)));
    qn::RoutingMatrix<T> P;
    const T one = num_traits<T>::from_int(1);
    P.set(k, k, a1, a2, one);
    P.set(k, k, a2, a1, one);
    m.link(P);
    return m;
}

TEST_CASE("ba: the LP linear-reduction bounds") {
    // The expected values are MATLAB's SolverBA(model,'method','lr.*').getAvg.
    // MEASURED deviation over all 14 compared quantities: max 1.02e-7 absolute,
    // 3.2e-7 relative, so 1e-6 is what the agreement actually supports and not
    // a tolerance chosen to make the test pass. The gap is linprog's, not this
    // port's -- see the identity check at the end of this case.
    const double kLr = 1e-6;
    qn::Network<double> m = model_b();
    SUBCASE("lr.upper") {
        const mva::AvgResult<double> r = run_ba(m, "lr.upper");
        CHECK(r.QN(0, 0) == doctest::Approx(2.61523045158).epsilon(kLr));
        CHECK(r.QN(1, 0) == doctest::Approx(1.74348696772).epsilon(kLr));
        CHECK(r.QN(2, 0) == doctest::Approx(1.04609218063).epsilon(kLr));
        // U is the LP objective per station, NOT the utilization law
        CHECK(r.UN(0, 0) == doctest::Approx(0.871743483859).epsilon(kLr));
        CHECK(r.UN(1, 0) == doctest::Approx(0.581162323628).epsilon(kLr));
        CHECK(r.UN(2, 0) == doctest::Approx(0.348697394367).epsilon(kLr));
        CHECK(r.RN(0, 0) == doctest::Approx(1.5).epsilon(kEps));
        CHECK(r.TN(0, 0) == doctest::Approx(1.74348696772).epsilon(kLr));
        CHECK(r.XN[0] == doctest::Approx(1.74348696772).epsilon(kLr));
        CHECK(r.CN[0] == doctest::Approx(3.1).epsilon(kEps));
    }
    SUBCASE("lr.lower") {
        const mva::AvgResult<double> r = run_ba(m, "lr.lower");
        CHECK(r.QN(0, 0) == doctest::Approx(0.793971048419).epsilon(kLr));
        CHECK(r.QN(1, 0) == doctest::Approx(0.52931403228).epsilon(kLr));
        CHECK(r.QN(2, 0) == doctest::Approx(0.317588419368).epsilon(kLr));
        CHECK(r.RN(0, 0) == doctest::Approx(0.5).epsilon(kEps));
        CHECK(r.TN(0, 0) == doctest::Approx(1.58794209684).epsilon(kLr));
        CHECK(r.CN[0] == doctest::Approx(1.03333333333).epsilon(kEps));
    }
    SUBCASE("bare 'lr' is the upper bound, as runAnalyzer resolves it") {
        const mva::AvgResult<double> r = run_ba(m, "lr");
        CHECK(r.actualmethod == "lr.upper");
        CHECK(r.TN(0, 0) == doctest::Approx(1.74348696772).epsilon(kLr));
    }
    SUBCASE("the bracket contains the exact throughput") {
        const double Xex = run_mva(m).TN(0, 0);
        CHECK(Xex == doctest::Approx(1.64006285667).epsilon(kEps));
        CHECK(run_ba(m, "lr.lower").TN(0, 0) <= Xex);
        CHECK(run_ba(m, "lr.upper").TN(0, 0) >= Xex);
    }
    SUBCASE("lr refuses a delay station and a multiserver station by name") {
        qn::Network<double> a = model_a();  // has a Delay
        CHECK_THROWS_AS(run_ba(a, "lr.upper"), UnsupportedError);
        CHECK_THROWS_AS(run_ba(a, "lr.lower"), UnsupportedError);
        qn::Network<double> d = model_d();  // 2 servers at Queue1, but also a Delay
        CHECK_THROWS_AS(run_ba(d, "lr.upper"), UnsupportedError);
    }
    SUBCASE("the LP is TIGHT on the two-station model, and exact under Rational") {
        // Both senses collapse onto the exact throughput 6/7 here, so the
        // polytope has pinned the stationary solution rather than merely
        // contained it. Under Rational the vertex is reported with no rounding
        // at all, which is what makes 6/7 checkable as an identity.
        qn::Network<double> s = model_lr_small<double>();
        const double Xex = run_mva(s).TN(0, 0);
        CHECK(Xex == doctest::Approx(6.0 / 7.0).epsilon(1e-12));
        CHECK(run_ba(s, "lr.upper").TN(0, 0) == doctest::Approx(6.0 / 7.0).epsilon(1e-12));
        CHECK(run_ba(s, "lr.lower").TN(0, 0) == doctest::Approx(6.0 / 7.0).epsilon(1e-12));

        qn::Network<Rational> sq = model_lr_small<Rational>();
        ba::BaOptions o;
        o.method = "lr.upper";
        CHECK(ba::solver_ba_run_analyzer(sq.get_struct(), o).XN[0] == Rational(6) / Rational(7));
        o.method = "lr.lower";
        CHECK(ba::solver_ba_run_analyzer(sq.get_struct(), o).XN[0] == Rational(6) / Rational(7));
    }
    SUBCASE("the port satisfies an identity MATLAB only approximates") {
        // At the optimal vertex of the lower-bound LP every station attains the
        // minimum of U_i mu_i / V_i, so that product is CONSTANT across
        // stations. The exact simplex reproduces it to machine precision;
        // MATLAB's interior-point linprog spreads it by 5.0e-7, which is larger
        // than its whole disagreement with this port (1.0e-7) and identifies
        // that disagreement as linprog's convergence gap rather than an error
        // here. Nothing transcribed above can fake this check.
        const mva::AvgResult<double> r = run_ba(m, "lr.lower");
        const double mu[3] = {2.0, 3.0, 5.0};
        const double p0 = r.UN(0, 0) * mu[0];
        for (std::size_t i = 1; i < 3; ++i)
            CHECK(r.UN(i, 0) * mu[i] == doctest::Approx(p0).epsilon(1e-12));
    }
}

/**
 * A blocking-blind family must be REFUSED on a model with a binding finite
 * buffer, not silently answered with the UNBLOCKED model's bound.
 *
 * Before the gate, `gb.upper` on `cqn_bas_blocking` returned QLen 1.28 at
 * Queue2 -- a station that can never hold more than one job -- because the cap
 * and the BAS drop rule were never read: every family here is parameterized by
 * demands and a population alone. `qrf.bas*`/`qrf.rsrd` carry the blocking
 * tables explicitly and stay available.
 */
TEST_CASE("ba: a binding finite buffer is refused, not bounded") {
    // The `cqn_bas_blocking` example: Queue1 declares BAS, Queue2 holds 1 job
    // against a population of 2, so its buffer BINDS.
    qn::Network<double> m("cqn_bas_blocking");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("Class1", 2.0, q1);
    m.set_service(q1, c, D::exp_rate(1.0));
    m.set_service(q2, c, D::exp_rate(0.8));
    m.set_capacity(q2, 1.0);
    m.set_drop_rule(q1, c, lang::DropStrategy::BAS);
    qn::RoutingMatrix<double> P;
    P.set(q1, q2, 1.0);
    P.set(q2, q1, 1.0);
    m.link(P);
    const qn::NetworkStruct<double> L = m.get_struct();
    REQUIRE(api::sn_has_blocking(L));

    SUBCASE("a blocking-blind family asked for BY NAME refuses") {
        CHECK_THROWS_AS(run_ba(m, "gb.upper"), UnsupportedError);
        CHECK_THROWS_AS(run_ba(m, "gb.lower"), UnsupportedError);
        CHECK_THROWS_AS(run_ba(m, "aba.upper"), UnsupportedError);
        // 'auto.lower' has no blocking counterpart to route to: the analyzer
        // solves qrf.bas in the 'max' direction alone.
        CHECK_THROWS_AS(run_ba(m, "auto.lower"), UnsupportedError);
    }
    SUBCASE("the upper-side aliases MEAN the QRF BAS bound here") {
        // 'default' and 'auto.upper' used to be refused with the blind
        // families, which left this model with no usable SolverBA method at
        // all: the advice was 'qrf.bas', and that in turn demanded hand-built
        // blocking tables. The tables are derived now, so the alias resolves
        // rather than refusing -- the routing SolverMVA performs to reach 'sqd'.
        const mva::AvgResult<double> named = run_ba(m, "qrf.bas");
        for (const char* alias : {"default", "auto", "auto.upper"}) {
            const mva::AvgResult<double> routed = run_ba(m, alias);
            CHECK(routed.UN(0, 0) == doctest::Approx(named.UN(0, 0)).epsilon(1e-9));
            CHECK(routed.UN(1, 0) == doctest::Approx(named.UN(1, 0)).epsilon(1e-9));
        }
    }
    SUBCASE("the narrowed list offers only the families that model the blocking") {
        const std::vector<std::string> valid = ba::list_valid_methods(L);
        CHECK_FALSE(valid.empty());
        // 'default' is offered back because it now MEANS 'qrf.bas' on this
        // model; it is the one entry whose RESOLVED form ('gb.upper') is blind,
        // which is exactly why the runner rewrites it rather than dispatching it.
        CHECK(std::find(valid.begin(), valid.end(), "default") != valid.end());
        for (std::size_t i = 0; i < valid.size(); ++i) {
            if (valid[i] == "default") continue;
            CHECK_FALSE(ba::ignores_blocking(ba::resolve_method(valid[i])));
        }
    }
    SUBCASE("getBounds refuses the routed default as one-sided, not blind") {
        // Saying "does not support blocking" here would contradict the run that
        // just succeeded; the real reason is that qrf.bas has no lower side.
        ba::BaOptions opt;
        opt.method = "default";
        CHECK_THROWS_AS(ba::ba_bounds(L, opt), UnsupportedError);
    }
    SUBCASE("an unblocked model of the same shape is untouched by the gate") {
        qn::Network<double> b = model_b();
        CHECK(run_ba(b, "gb.upper").XN[0] > 0.0);
        CHECK_FALSE(ba::list_valid_methods(b.get_struct()).empty());
    }
}

}  // namespace
