/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The conservation-law enumerator, `lqn::lqn_balance_equations`.
 *
 * Two things are checked, and they are different in kind. STRUCTURE, with no solution
 * supplied: the right relations come out with the right term sets, branches and
 * constants, which depends on the model alone. CONSISTENCY, on a converged SolverLN:
 * every relation the enumerator emits is a law the solution must satisfy, so each
 * residual must vanish to the solver's own tolerance. The second is the real content
 * -- it turns the enumerator into a conservation check on the layered solver itself.
 *
 * The model is the three-tier LQN of `lqn_basic`: T1 (reference, 50 threads, think 2)
 * calls T2 (50 threads) once, which calls T3 (25 threads) five times, on processors of
 * multiplicity 2 and 3. Twin of the MATLAB, JAR and Python tests of the same name.
 */

#include <cmath>
#include <string>

#include "doctest.h"
#include "line/api/lqn/lqn_balance_equations.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/solvers/ln/solver_ln.h"

using namespace line;
using namespace line::lang;
using D = Distrib<double>;

namespace {

const double TOL = 1e-6;

lqn::LqnStruct<double> build() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 2, SchedStrategy::PS);
    b.processor("P2", 3, SchedStrategy::PS);
    b.task("T1", 50, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(2.0));
    b.task("T2", 50, SchedStrategy::FCFS, "P1");
    b.think_time("T2", D::exp_mean(3.0));
    b.task("T3", 25, SchedStrategy::FCFS, "P2");
    b.think_time("T3", D::exp_mean(4.0));
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T3");
    b.activity("AS1", D::exp_mean(0.1), "T1");
    b.bound_to("AS1", "E1");
    b.sync_call("AS1", "E2", 1.0);
    b.activity("AS2", D::exp_mean(0.05), "T2");
    b.bound_to("AS2", "E2");
    b.sync_call("AS2", "E3", 5.0);
    b.replies_to("AS2", "E2");
    b.activity("AS3", D::exp_mean(0.02), "T3");
    b.bound_to("AS3", "E3");
    b.replies_to("AS3", "E3");
    return b.build();
}

/** The relation of a given kind anchored on the element whose hashname ENDS in NAME. */
const lqn::LqnRelation<double>& by(const lqn::LqnBalanceEquations<double>& out,
                                  const std::string& kind, const std::string& name) {
    for (const lqn::LqnRelation<double>& r : out.eqs) {
        if (r.kind != kind) continue;
        if (r.targetname == name || (r.targetname.size() >= name.size() &&
                                     r.targetname.compare(r.targetname.size() - name.size(),
                                                          name.size(), name) == 0))
            return r;
    }
    FAIL("no ", kind, " relation on ", name);
    return out.eqs[0];
}

std::size_t count(const lqn::LqnBalanceEquations<double>& out, const std::string& kind) {
    std::size_t n = 0;
    for (const lqn::LqnRelation<double>& r : out.eqs)
        if (r.kind == kind) ++n;
    return n;
}

}  // namespace

TEST_CASE("balance equations: the five families come out of the structure alone") {
    const lqn::LqnStruct<double> l = build();
    const lqn::LqnBalanceEquations<double> out = lqn::lqn_balance_equations(l);

    CHECK(count(out, "little") == 3);      // one per task
    CHECK(count(out, "callflow") == 2);    // one per call
    CHECK(count(out, "entryflow") == 2);   // E2 and E3; E1 rides the reference cycle
    CHECK(count(out, "actflow") == 3);     // one per activity
    CHECK(count(out, "hostutil") == 2);    // one per processor

    // the thread-pool branches and their populations
    const lqn::LqnRelation<double>& r1 = by(out, "little", "T1");
    CHECK(r1.branch == "ref");
    CHECK(r1.mult == 50.0);
    REQUIRE(r1.terms.size() == 1);
    CHECK(r1.termisentry[0]);              // its own entry drives its cycle
    const lqn::LqnRelation<double>& r2 = by(out, "little", "T2");
    CHECK(r2.branch == "queueing");
    CHECK(r2.mult == 50.0);
    REQUIRE(r2.terms.size() == 1);
    CHECK_FALSE(r2.termisentry[0]);        // one call class
    CHECK(r2.scaled);                      // U is normalized to [0,1] here
    CHECK(by(out, "little", "T3").mult == 25.0);

    // the host law carries the declared multiplicity and the host demands
    const lqn::LqnRelation<double>& p1 = by(out, "hostutil", "P1");
    CHECK(p1.rhsconst == 2.0);
    REQUIRE(p1.coeff.size() == 2);
    CHECK(p1.coeff[0] + p1.coeff[1] == doctest::Approx(0.15));
    const lqn::LqnRelation<double>& p2 = by(out, "hostutil", "P2");
    CHECK(p2.rhsconst == 3.0);
    REQUIRE(p2.coeff.size() == 1);
    CHECK(p2.coeff[0] == doctest::Approx(0.02));

    // the call flow carries the mean call count, and a single-activity entry one visit;
    // a call hashname prefixes BOTH sides by element kind, hence A:AS1=>E:E2
    CHECK(by(out, "callflow", "A:AS1=>E:E2").coeff[0] == doctest::Approx(1.0));
    CHECK(by(out, "callflow", "A:AS2=>E:E3").coeff[0] == doctest::Approx(5.0));
    for (const lqn::LqnRelation<double>& r : out.eqs)
        if (r.kind == "actflow") CHECK(r.coeff[0] == doctest::Approx(1.0));

    // the aggregation incidence carries the CALL classes only
    CHECK(out.A_little(2, 1) == 1.0);      // T2 is a caller class of call 1
    CHECK(out.A_little(3, 2) == 1.0);      // T3 of call 2
    double t1row = 0;
    for (std::size_t c = 0; c <= l.ncalls; ++c) t1row += out.A_little(1, c);
    CHECK(t1row == 0.0);                   // T1 is entry-driven, no call class

    // and nothing is instantiated without a solution
    CHECK_FALSE(out.has_maxresidual);
    CHECK(out.str().find("thread-pool Little") != std::string::npos);
    CHECK(out.str().find("host utilization law") != std::string::npos);
}

TEST_CASE("balance equations: every relation vanishes on a converged solve") {
    const lqn::LqnStruct<double> l = build();
    ln::SolverLN<double> s(l, ln::LnOptions());
    const ln::LnSolution<double> avg = s.get_ensemble_avg();

    lqn::LqnSolution<double> sol;
    sol.tput = s.state_tput();
    sol.util = s.state_util();
    sol.thinkt = s.state_thinkt();
    sol.servt = s.state_servt();
    sol.residt = s.state_residt();
    // The REPORTED utilization, which is not the `util` iterate: that one holds a
    // task's utilization as a server in its own task layer and is left at zero on a
    // host, and the host law is about the reported one.
    sol.un = avg.UN;

    const lqn::LqnBalanceEquations<double> out = lqn::lqn_balance_equations(l, &sol);
    for (const lqn::LqnRelation<double>& r : out.eqs) {
        if (r.degenerate || !r.has_residual) continue;
        CHECK_MESSAGE(std::fabs(r.residual) <= TOL,
                      r.kind, " ", r.targetname, ": residual ", r.residual, " (lhs ", r.lhs,
                      ", rhs ", r.rhs, ")");
    }
    REQUIRE(out.has_maxresidual);
    CHECK(out.maxresidual <= TOL);
}

TEST_CASE("balance equations: the per-class utilization reproduces the solver iterate") {
    // The per-call decomposition of the busy threads must add up to the utilization the
    // solver carries for that task, which is what makes the emitted per-class U a
    // refinement of the aggregate rather than a new quantity. The two reach it by
    // different routes through the same iterate, so they agree to the fixed point's
    // tolerance and not to machine precision.
    const lqn::LqnStruct<double> l = build();
    ln::SolverLN<double> s(l, ln::LnOptions());
    s.get_ensemble_avg();

    lqn::LqnSolution<double> sol;
    sol.tput = s.state_tput();
    sol.util = s.state_util();
    sol.thinkt = s.state_thinkt();
    sol.servt = s.state_servt();
    sol.residt = s.state_residt();   // un left empty on purpose, see below

    const lqn::LqnBalanceEquations<double> out = lqn::lqn_balance_equations(l, &sol);
    const char* names[2] = {"T2", "T3"};
    for (const char* nm : names) {
        const lqn::LqnRelation<double>& r = by(out, "little", nm);
        double sum = 0;
        for (double u : r.perclassutil) sum += u;
        CHECK(sum == doctest::Approx(s.state_util()[r.target]).epsilon(1e-6));
    }

    // Without the reported utilization the host law stays symbolic rather than report a
    // residual against a zero it never meant; every other family still instantiates.
    for (const lqn::LqnRelation<double>& r : out.eqs) {
        if (r.kind == "hostutil")
            CHECK_FALSE(r.has_residual);
        else if (!r.degenerate)
            CHECK(r.has_residual);
    }
}
