/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The derived START and PREEMPT filtrations of the CTMC generator.
 *
 * The identity that pins them is
 *
 *     startRate(i,r) == throughput(i,r) + preemptRate(i,r)
 *
 * at a lossless station with no in-service abandonment: every job starts
 * service once per entry into a server, and every preemption is followed by
 * exactly one later resume or restart. At a non-preemptive station it collapses
 * to startRate == throughput, which is an exact oracle rather than a
 * regression-recorded number.
 *
 * The tags are annotations on arcs the generator already carries, so they are
 * kept OUT of the event filtration -- adding them there would double-count the
 * off-diagonal of Q -- and every mean measure must be what it was without them.
 */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_getters.h"
#include "line/solvers/ssa/ssa_dispatch.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Source -> Queue(sched) -> Sink, one open class. */
qn::Network<double> mm1(SchedStrategy sched) {
    qn::Network<double> m("mm1");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", sched);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dist::exp_rate(0.5));
    m.set_service(q, c, Dist::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/** Source -> Queue(sched) -> Sink, an urgent and a normal open class. */
qn::Network<double> mm1_prio(SchedStrategy sched = SchedStrategy::FCFSPRPRIO) {
    qn::Network<double> m("mm1prio");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", sched);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t urgent = m.add_open_class("Urgent", 0);
    const std::size_t normal = m.add_open_class("Normal", 1);
    m.set_arrival(src, urgent, Dist::exp_rate(0.4));
    m.set_arrival(src, normal, Dist::exp_rate(0.4));
    m.set_service(q, urgent, Dist::exp_rate(1.0));
    m.set_service(q, normal, Dist::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    // PER CLASS: the two-index `set` is shorthand for class 1 alone, so routing
    // the model that way leaves the second class with no path and no traffic --
    // and then every derived rate for it is trivially zero.
    for (std::size_t r = 1; r <= 2; ++r) {
        P.set(r, r, src, q, 1.0);
        P.set(r, r, q, k, 1.0);
    }
    m.link(P);
    return m;
}

ctmc::CtmcOptions opts_with_tags(double cutoff) {
    ctmc::CtmcOptions o;
    o.cutoff = cutoff;
    // the derived filtrations ride with the event filtration
    o.keep_filtration = true;
    return o;
}

}  // namespace

TEST_CASE("ctmc start filtration: a non-preemptive station starts one service per departure") {
    // FCFS, PS and a two-server FCFS all promote without ever displacing anyone,
    // so the identity collapses to startRate == TN at each of them.
    for (int which = 0; which < 3; ++which) {
        qn::Network<double> m =
            mm1(which == 1 ? SchedStrategy::PS : SchedStrategy::FCFS);
        if (which == 2) m.set_number_of_servers(2, 2.0);
        const ctmc::CtmcSolution<double> sol =
            ctmc::solver_ctmc_analyzer(m.get_struct(), opts_with_tags(4));
        CAPTURE(which);
        CHECK(sol.avg.StartN(1, 0) == doctest::Approx(sol.avg.TN(1, 0)).epsilon(1e-9));
        CHECK(sol.avg.PreemptN(1, 0) == doctest::Approx(0.0).epsilon(1e-12));
        // a Source CREATES jobs rather than admitting them to service, so it
        // seizes nothing and its row of both derived rates is zero
        CHECK(sol.avg.StartN(0, 0) == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(sol.avg.PreemptN(0, 0) == doctest::Approx(0.0).epsilon(1e-12));
    }
}

TEST_CASE("ctmc start filtration: preempt-resume and preempt-independent displace alike") {
    // PR and PI differ in the phase the displaced job resumes in, not in how
    // often it is displaced, so the preemption rate is the SAME number.
    const ctmc::CtmcSolution<double> pr =
        ctmc::solver_ctmc_analyzer(mm1_prio(SchedStrategy::FCFSPRPRIO).get_struct(),
                                   opts_with_tags(3));
    const ctmc::CtmcSolution<double> pi =
        ctmc::solver_ctmc_analyzer(mm1_prio(SchedStrategy::FCFSPIPRIO).get_struct(),
                                   opts_with_tags(3));
    for (std::size_t r = 0; r < 2; ++r)
        CHECK(pi.avg.PreemptN(1, r) == doctest::Approx(pr.avg.PreemptN(1, r)).epsilon(1e-12));
    CHECK(pr.avg.PreemptN(1, 1) > 0.0);
}

TEST_CASE("ctmc start filtration: a closed cyclic model satisfies the identity at both stations") {
    qn::Network<double> m("cyclic");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("Jobs", 3, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    ctmc::CtmcOptions o;
    o.keep_filtration = true;
    const ctmc::CtmcSolution<double> sol = ctmc::solver_ctmc_analyzer(m.get_struct(), o);
    // both stations hold jobs here, so the identity is asserted on every row
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(sol.avg.StartN(i, 0) ==
              doctest::Approx(sol.avg.TN(i, 0) + sol.avg.PreemptN(i, 0)).epsilon(1e-9));
        CHECK(sol.avg.PreemptN(i, 0) == doctest::Approx(0.0).epsilon(1e-12));
    }
}

TEST_CASE("ctmc start filtration: the derived tags stay OUT of the event filtration") {
    // A START rides on an arc the event filtration already carries, so the two
    // must not be mixed: `filt` stays paired one-to-one with `sync` and still
    // sums to the off-diagonal of Q. Summing a filtration that had absorbed the
    // derived tags would double-count those arcs, which is the defect this
    // guards against.
    qn::Network<double> m = mm1_prio(SchedStrategy::FCFSPRPRIO);
    const qn::NetworkStruct<double> sn = m.get_struct();
    const ctmc::CtmcGenerator<double> g =
        ctmc::ctmc_get_generator(sn, ctmc::solver_ctmc_analyzer(sn, opts_with_tags(3)));
    REQUIRE(g.filt.size() == g.sync.size());
    const std::size_t n = g.Q.rows();
    double gap = 0.0;
    for (std::size_t s = 0; s < n; ++s)
        for (std::size_t ns = 0; ns < n; ++ns) {
            if (s == ns) continue;
            double acc = 0.0;
            for (std::size_t a = 0; a < g.filt.size(); ++a) acc += g.filt[a](s, ns);
            gap = std::max(gap, std::abs(acc - g.Q(s, ns)));
        }
    CHECK(gap == doctest::Approx(0.0).epsilon(1e-9));
    // and the derived filtration is present, per station and per class
    REQUIRE(g.start_filt.size() == sn.nstations);
    REQUIRE(g.preempt_filt.size() == sn.nstations);
    CHECK(g.start_filt[1].size() == sn.nclasses);
    CHECK(g.preempt_filt[1][1].rows() == n);
}

TEST_CASE("ssa start counters: BOTH engines report the same derived rates") {
    // The counters must not depend on `options.method`: the NRM engine counts
    // the events it fires and divides by the simulated time, the serial engine
    // accumulates the tag rate of the enabled transitions. At an FCFS station
    // nothing is ever displaced, so both must land on startRate == TN with no
    // preemption, which is an exact statement about the sample path rather than
    // a recorded number.
    qn::Network<double> m = mm1(SchedStrategy::FCFS);
    for (const char* method : {"nrm", "serial"}) {
        line::ssa::SsaOptions o;
        o.method = method;
        o.samples = 200000;
        o.seed = 23000;
        const line::ssa::SsaSolution sol = line::ssa::solver_ssa(m.get_struct(), o);
        CAPTURE(method);
        CHECK(sol.StartN(1, 0) == doctest::Approx(sol.TN(1, 0)).epsilon(0.05));
        CHECK(sol.StartN(1, 0) > 0.0);
        CHECK(sol.PreemptN(1, 0) == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(sol.StartN(0, 0) == doctest::Approx(0.0).epsilon(1e-12));  // the Source seizes nothing
    }
}

TEST_CASE("ctmc start filtration: fcfsprprio satisfies startRate == TN + preemptRate") {
    qn::Network<double> m = mm1_prio();
    const ctmc::CtmcSolution<double> sol =
        ctmc::solver_ctmc_analyzer(m.get_struct(), opts_with_tags(3));
    for (std::size_t r = 0; r < 2; ++r)
        CHECK(sol.avg.StartN(1, r) ==
              doctest::Approx(sol.avg.TN(1, r) + sol.avg.PreemptN(1, r)).epsilon(1e-9));
    // the urgent class is never displaced; the lower-priority one must be
    CHECK(sol.avg.PreemptN(1, 0) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(sol.avg.PreemptN(1, 1) > 0.0);
}
