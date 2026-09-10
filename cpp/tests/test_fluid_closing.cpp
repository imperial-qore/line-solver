/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * `solver_fluid_initsol` and the closing entry point built on it.
 *
 * WHY THE ORACLES HERE ARE NOT FLUID NUMBERS. The fluid limit is an
 * approximation of the queueing model, so comparing it against an exact solver
 * proves nothing unless the gap is known; every check below is instead a
 * property the fluid answer must satisfy exactly, whatever the approximation
 * error happens to be.
 *
 * For the initial condition those properties are ARITHMETIC. Phase one of a
 * block gets `nir - sum_{k>=2} kir` and the other phases get `kir`, so the
 * block sums to `nir` by cancellation, for every encoding and every service
 * process: the population the model starts with is the population the ODE
 * starts with. And on a model whose service is exponential the decode has to
 * land on the solver's own default y0, which is derived by an entirely
 * different route -- one places the population from the class table, the other
 * encodes it into a state row and decodes it back.
 *
 * For the solve they are the two places the fluid limit is EXACT rather than
 * approximate: a network of infinite servers, where the drift is linear and its
 * fixed point is the exact product-form answer, and a saturated bottleneck,
 * where the throughput is the service rate by definition. Both are stated as
 * derivations in the comments, not as tables.
 */
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/fluid/fluid_closing.h"
#include "line/solvers/fluid/solver_fluid.h"

using namespace line;
using D = lang::Distrib<double>;

namespace {

std::size_t st_index(const qn::NetworkStruct<double>& sn, const std::string& nm) {
    for (std::size_t i = 0; i < sn.nstations; ++i)
        if (sn.stations[i].name == nm) return i;
    FAIL("no station named ", nm);
    return 0;
}

/** Think(exp 1) <-> Queue(exp 2) under `sc`, one closed class of `n` jobs. */
qn::Network<double> closed_pair(double n, lang::SchedStrategy sc) {
    qn::Network<double> m("fluid_initsol_pair");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Queue", sc);
    const std::size_t c = m.add_closed_class("C1", n, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("fluid initsol: the decoded default state is the solver's own default y0") {
    // Exponential service everywhere, so every block is one phase long and the
    // two constructions must agree entry for entry. They share no code:
    // `fluid_default_initsol` reads the class table and writes the population
    // straight into the reference station's block, while `fluid_initsol` builds
    // the encoded state row through `from_marginal` and decodes it back through
    // `to_marginal`. Agreement is therefore evidence about both, and a decode
    // that lost the buffer or landed on the wrong phase shows up immediately.
    qn::Network<double> m = closed_pair(4.0, lang::SchedStrategy::PS);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const fluid::FluidLayout L = fluid::fluid_layout(sn);

    const std::vector<double> got = fluid::fluid_initsol(sn, L);
    const std::vector<double> want = fluid::detail::fluid_default_initsol(sn, L);
    REQUIRE(got.size() == L.nstates);
    REQUIRE(want.size() == got.size());
    for (std::size_t a = 0; a < got.size(); ++a)
        CHECK(got[a] == doctest::Approx(want[a]).epsilon(1e-12));

    // The population the model declares is the population the drift starts
    // with, which is the invariant the whole initial condition exists to carry.
    double total = 0;
    for (std::size_t a = 0; a < got.size(); ++a) total += got[a];
    CHECK(total == doctest::Approx(4.0).epsilon(1e-12));

    // And it starts where the class says it starts: all of it at the reference
    // station, none at the other. A start that merely conserved the mass would
    // pass the sum above and still integrate a different transient.
    const std::size_t id = st_index(sn, "Think"), iq = st_index(sn, "Queue");
    CHECK(got[L.qidx[id][0]] == doctest::Approx(4.0).epsilon(1e-12));
    CHECK(got[L.qidx[iq][0]] == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("fluid initsol: a buffered multi-phase station keeps its whole population") {
    // Erlang-2 service at an FCFS queue that is ALSO the reference station, so
    // the initial state has one job in service in some phase and the rest in
    // the buffer. The fluid state has no buffer, so those waiting jobs must
    // reappear in phase one: the block sums to nir by cancellation whichever
    // phase the encoder chose for the job in service, which is what makes this
    // a check on the rule rather than on the encoder.
    const double N = 3.0;
    qn::Network<double> m("fluid_initsol_erlang");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", N, q);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::erlang(4.0, 2));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    const fluid::FluidLayout L = fluid::fluid_layout(sn);
    const std::vector<double> y0 = fluid::fluid_initsol(sn, L);
    REQUIRE(y0.size() == L.nstates);

    const std::size_t id = st_index(sn, "Think"), iq = st_index(sn, "Queue");
    REQUIRE(L.kic[iq][0] == 2);
    double block = 0;
    for (std::size_t k = 0; k < L.kic[iq][0]; ++k) {
        // No phase may hold a negative mass: a sign error in the buffer term
        // produces exactly that, and the integrator would then clamp it away
        // and hide the mistake.
        CHECK(y0[L.qidx[iq][0] + k] >= 0.0);
        block += y0[L.qidx[iq][0] + k];
    }
    CHECK(block == doctest::Approx(N).epsilon(1e-12));
    CHECK(y0[L.qidx[id][0]] == doctest::Approx(0.0).epsilon(1e-12));

    double total = 0;
    for (std::size_t a = 0; a < y0.size(); ++a) total += y0[a];
    CHECK(total == doctest::Approx(N).epsilon(1e-12));
}

TEST_CASE("fluid initsol: a capped reference station spills the excess population") {
    // Think <-> PS Queue with the QUEUE as the reference station and a buffer of
    // two, four jobs in the class. `Network.initDefault` fills the reference
    // station up to `min(classcap, cap)` and then walks the remaining stations in
    // ascending order, so two jobs stay at the queue and two go to the delay;
    // starting all four at the queue would put the drift at a point the model
    // cannot occupy. The oracle is an OBSERVED MATLAB run, not a derivation:
    // LINE 3.0.7 solver_fluid_initsol on this model prints init_sol = [2 2] over
    // the station order (Think, Queue), with sn.classcap = [4; 2].
    const double N = 4.0, C = 2.0;
    qn::Network<double> m("fluid_initsol_capped");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", N, q);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    m.set_capacity(q, C);
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    const fluid::FluidLayout L = fluid::fluid_layout(sn);
    const std::vector<double> y0 = fluid::fluid_initsol(sn, L);
    REQUIRE(y0.size() == L.nstates);

    const std::size_t id = st_index(sn, "Think"), iq = st_index(sn, "Queue");
    // The capacity is a property of the struct and not of this test: if
    // refreshCapacity stopped deriving it the placement below would be vacuous.
    REQUIRE(sn.classcap[iq][0] == doctest::Approx(C).epsilon(1e-12));
    CHECK(y0[L.qidx[iq][0]] == doctest::Approx(C).epsilon(1e-12));
    CHECK(y0[L.qidx[id][0]] == doctest::Approx(N - C).epsilon(1e-12));

    double total = 0;
    for (std::size_t a = 0; a < y0.size(); ++a) total += y0[a];
    CHECK(total == doctest::Approx(N).epsilon(1e-12));

    // And a population that fits nowhere is refused rather than silently
    // truncated, which is the reference's own error arm in initDefault.
    qn::Network<double> tight("fluid_initsol_overfull");
    const std::size_t td = tight.add_delay("Think");
    const std::size_t tq = tight.add_queue("Queue", lang::SchedStrategy::PS);
    const std::size_t tc = tight.add_closed_class("C1", N, tq);
    tight.set_service(td, tc, D::exp_rate(1.0));
    tight.set_service(tq, tc, D::exp_rate(2.0));
    tight.set_capacity(tq, 1.0);
    tight.set_capacity(td, 1.0);
    qn::RoutingMatrix<double> TP;
    TP.set(tc, tc, td, tq, 1.0);
    TP.set(tc, tc, tq, td, 1.0);
    tight.link(TP);
    CHECK_THROWS_AS(fluid::fluid_initsol(tight.get_struct()), line::InputError);
}

TEST_CASE("fluid initsol: a Place reference station takes the whole marking and never spills") {
    // `initDefault.m:24-28` tests the reference station for NodeType.Place
    // BEFORE the placement loop and, when it is one, assigns the whole
    // population there and skips the spill entirely -- a Place is a token
    // container, so its marking is the model's own and no buffer bounds it. The
    // capacity below is what makes the branch observable: through the ordinary
    // loop two tokens would stay at P1 and two would spill to P2.
    const double N = 4.0;
    qn::Network<double> m("fluid_initsol_place");
    const std::size_t p1 = m.add_place("P1");
    const std::size_t p2 = m.add_place("P2");
    const std::size_t c = m.add_closed_class("Tok", N, p1);
    m.set_service(p1, c, D::exp_rate(1.0));
    m.set_service(p2, c, D::exp_rate(1.0));
    m.set_capacity(p1, 2.0);

    qn::TransitionParam<double> tp;
    tp.nmodes = 1;
    tp.modenames.push_back("fire");
    tp.enabling.assign(1, line::Matrix<double>(3, 1, 0.0));
    tp.inhibiting.assign(1, line::Matrix<double>(3, 1, std::numeric_limits<double>::infinity()));
    tp.firing.assign(1, line::Matrix<double>(3, 1, 0.0));
    tp.enabling[0](p1 - 1, 0) = 1.0;
    tp.firing[0](p2 - 1, 0) = 1.0;
    tp.nmodeservers.push_back(1.0);
    tp.firingphases.push_back(1);
    tp.timing.push_back(lang::TimingStrategy::TIMED);
    tp.fireweight.push_back(1.0);
    tp.firingproc.push_back(D::exp_rate(2.0));
    m.add_transition("T1", tp);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    const fluid::FluidLayout L = fluid::fluid_layout(sn);
    const std::vector<double> y0 = fluid::fluid_initsol(sn, L);
    REQUIRE(y0.size() == L.nstates);

    const std::size_t i1 = st_index(sn, "P1"), i2 = st_index(sn, "P2");
    // The buffer is real, so the branch is the only thing that can put four
    // tokens here; without it the placement would read 2 and 2.
    REQUIRE(sn.classcap[i1][0] == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(y0[L.qidx[i1][0]] == doctest::Approx(N).epsilon(1e-12));
    CHECK(y0[L.qidx[i2][0]] == doctest::Approx(0.0).epsilon(1e-12));

    double total = 0;
    for (std::size_t a = 0; a < y0.size(); ++a) total += y0[a];
    CHECK(total == doctest::Approx(N).epsilon(1e-12));
}

TEST_CASE("fluid initsol: an open model starts with the Source's unit job pool") {
    // The EXT drift keeps unit mass per class at the source and lets phase one
    // absorb whatever the other phases do not hold (`ode_rates_closing`, and
    // the EXT branch of fluid_odes.h). The initial condition has to hand it
    // that unit, so this is the drift's own precondition and not a value read
    // back out of the decoder. The queue starts empty because no job has
    // arrived yet, which is what makes an open transient meaningful.
    qn::Network<double> m("fluid_initsol_open");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, D::exp_rate(0.5));
    m.set_service(q, c, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    const fluid::FluidLayout L = fluid::fluid_layout(sn);
    const std::vector<double> y0 = fluid::fluid_initsol(sn, L);
    REQUIRE(y0.size() == L.nstates);

    const std::size_t is = st_index(sn, "Src"), iq = st_index(sn, "Queue");
    CHECK(y0[L.qidx[is][0]] == doctest::Approx(1.0).epsilon(1e-12));
    for (std::size_t k = 0; k < L.kic[iq][0]; ++k)
        CHECK(y0[L.qidx[iq][0] + k] == doctest::Approx(0.0).epsilon(1e-12));

    // The source's marginal is an infinite reservoir, so a decoder that read
    // its population instead of its in-service counts would put an infinity
    // into the state. Nothing in the initial condition may be non-finite.
    for (std::size_t a = 0; a < y0.size(); ++a) CHECK(std::isfinite(y0[a]));
}

TEST_CASE("fluid initsol: a discipline the reference cannot decode is refused by name") {
    // The reference's switch lists nine disciplines and errors on the rest. The
    // gate is checked directly for all nine so the list is asserted rather than
    // sampled, since a missing arm would not fail any model that does not use
    // it.
    const lang::SchedStrategy ok[] = {
        lang::SchedStrategy::EXT,  lang::SchedStrategy::FCFS, lang::SchedStrategy::SIRO,
        lang::SchedStrategy::PS,   lang::SchedStrategy::INF,  lang::SchedStrategy::DPS,
        lang::SchedStrategy::GPS,  lang::SchedStrategy::HOL,  lang::SchedStrategy::LCFS,
        lang::SchedStrategy::LCFSPR};
    for (std::size_t i = 0; i < sizeof(ok) / sizeof(ok[0]); ++i)
        CHECK_NOTHROW(fluid::detail::fluid_initsol_check_sched(ok[i], 1));

    // Refused for two distinct reasons: POLLING and SRPT hold a buffer this decode
    // cannot attribute to a phase, and FCFSPR records the phase of every preempted
    // job, which restarting the buffer in phase one discards.
    const lang::SchedStrategy no[] = {lang::SchedStrategy::POLLING, lang::SchedStrategy::SRPT,
                                      lang::SchedStrategy::FCFSPR, lang::SchedStrategy::SEPT};
    for (std::size_t i = 0; i < sizeof(no) / sizeof(no[0]); ++i)
        CHECK_THROWS_AS(fluid::detail::fluid_initsol_check_sched(no[i], 1), line::UnsupportedError);

    // And the refusal is reached through the ordinary entry point, not only by
    // calling the gate: a model this port cannot start must not be started.
    qn::Network<double> m = closed_pair(2.0, lang::SchedStrategy::POLLING);
    CHECK_THROWS_AS(fluid::fluid_initsol(m.get_struct()), line::UnsupportedError);

    // GPS MOVED FROM THE SECOND LIST TO THE FIRST. `solver_fluid_initsol.m:34`
    // lists it: its state encoding folds into one mass per phase exactly like PS's,
    // so the decode is fine and this port used to refuse a model the reference
    // starts. What GPS actually lacks is a first-order capacity SHARE, and that
    // refusal belongs to the per-method featset gate (only `minnormal` supplies the
    // backlog probability), not to the initial condition.
    qn::Network<double> g = closed_pair(2.0, lang::SchedStrategy::GPS);
    CHECK_NOTHROW(fluid::fluid_initsol(g.get_struct()));
}

TEST_CASE("fluid closing: a network of infinite servers hits its exact fixed point") {
    // Delay(rate 1) -> Delay(rate 3) -> Delay(rate 1), six jobs. With no
    // contention anywhere the drift is LINEAR, so the fluid limit is not an
    // approximation here: its fixed point is the exact stationary mean. Flow
    // balance q1*mu1 = q2*mu2 with q1 + q2 = N gives q1 = 4.5 and q2 = 1.5,
    // which is also the product-form answer (queue length proportional to the
    // service demand at an infinite server).
    const double N = 6.0, mu1 = 1.0, mu2 = 3.0;
    qn::Network<double> m("fluid_closing_delays");
    const std::size_t d1 = m.add_delay("D1");
    const std::size_t d2 = m.add_delay("D2");
    const std::size_t c = m.add_closed_class("C1", N, d1);
    m.set_service(d1, c, D::exp_rate(mu1));
    m.set_service(d2, c, D::exp_rate(mu2));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d1, d2, 1.0);
    P.set(c, c, d2, d1, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    const fluid::FluidSolution s = fluid::solver_fluid_closing(sn, fluid::FluidOptions());
    const std::size_t i1 = st_index(sn, "D1"), i2 = st_index(sn, "D2");

    const double q1 = N * (1.0 / mu1) / (1.0 / mu1 + 1.0 / mu2);
    CHECK(s.QN(i1, 0) == doctest::Approx(q1).epsilon(1e-4));
    CHECK(s.QN(i2, 0) == doctest::Approx(N - q1).epsilon(1e-4));
    // Conservation is exact in the drift, so it is checked far tighter than the
    // integrator tolerance the values above are compared at.
    CHECK(s.QN(i1, 0) + s.QN(i2, 0) == doctest::Approx(N).epsilon(1e-8));
    // Flow balance around the cycle: what leaves one delay enters the other.
    CHECK(s.TN(i1, 0) == doctest::Approx(s.TN(i2, 0)).epsilon(1e-6));
    // Little's law at an infinite server is an identity, not an approximation:
    // the residence time is the service time, whatever the queue length is.
    CHECK(s.RN(i1, 0) == doctest::Approx(1.0 / mu1).epsilon(1e-6));
    CHECK(s.RN(i2, 0) == doctest::Approx(1.0 / mu2).epsilon(1e-6));
}

TEST_CASE("fluid closing: the closing entry point agrees with the analyzer on the same drift") {
    // Think(rate 1) <-> PS Queue(rate 2), four jobs. The queue is the only
    // finite-server station and the population is well above one, so the fluid
    // fixed point saturates it: the throughput is then the service rate 2 by
    // definition, the delay holds X*Z = 2 jobs, and the queue holds the other
    // two. No approximation gap is being tolerated here -- these are the fluid
    // limit's own values, and they coincide with the MATLAB getAvgTable already
    // recorded for this model in test_fluid.cpp.
    qn::Network<double> m = closed_pair(4.0, lang::SchedStrategy::PS);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t id = st_index(sn, "Think"), iq = st_index(sn, "Queue");

    fluid::FluidOptions opt;
    opt.method = "closing";
    const fluid::FluidSolution s = fluid::solver_fluid_closing(sn, opt);

    CHECK(s.QN(id, 0) == doctest::Approx(2.0).epsilon(1e-4));
    CHECK(s.QN(iq, 0) == doctest::Approx(2.0).epsilon(1e-4));
    CHECK(s.TN(iq, 0) == doctest::Approx(2.0).epsilon(1e-4));
    CHECK(s.QN(id, 0) + s.QN(iq, 0) == doctest::Approx(4.0).epsilon(1e-8));

    // The analyzer reaches the same fixed point from the same start, because on
    // this model the decoded initial condition and the solver's default y0 are
    // the same vector (the first test case above). Queue lengths and
    // throughputs are therefore identical to the last digit -- the analyzer's
    // trailing correction rewrites only U and R, which is why they are not
    // compared here.
    const fluid::FluidSolution a = fluid::solver_fluid(sn, opt);
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        CHECK(s.QN(i, 0) == doctest::Approx(a.QN(i, 0)).epsilon(1e-12));
        CHECK(s.TN(i, 0) == doctest::Approx(a.TN(i, 0)).epsilon(1e-12));
    }
}

TEST_CASE("fluid closing: the initial condition is consumed, and only selects the path here") {
    qn::Network<double> m = closed_pair(4.0, lang::SchedStrategy::PS);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const fluid::FluidLayout L = fluid::fluid_layout(sn);
    const std::size_t id = st_index(sn, "Think"), iq = st_index(sn, "Queue");

    // A state of the wrong length is rejected rather than padded, which is also
    // what proves the vector is read at all: a solver that ignored `init_sol`
    // would pass every agreement check below by accident.
    fluid::FluidOptions bad;
    bad.init_sol.assign(L.nstates + 1, 0.0);
    CHECK_THROWS_AS(fluid::solver_fluid_closing(sn, bad), line::InputError);

    // Start the whole population at the queue instead of the delay. This model
    // has ONE fixed point -- a single-class network with a monotone drift -- so
    // the start can only change the transient, and the steady state must be the
    // same. The check is deliberately not made on a model with several fixed
    // points: that is precisely the case in which reproducing the reference's
    // initial condition, rather than any mass-preserving guess, is what makes
    // the answer the reference's answer.
    fluid::FluidOptions perturbed;
    perturbed.init_sol.assign(L.nstates, 0.0);
    perturbed.init_sol[L.qidx[iq][0]] = 4.0;
    const fluid::FluidSolution p = fluid::solver_fluid_closing(sn, perturbed);
    const fluid::FluidSolution d = fluid::solver_fluid_closing(sn, fluid::FluidOptions());
    CHECK(p.QN(id, 0) == doctest::Approx(d.QN(id, 0)).epsilon(1e-4));
    CHECK(p.QN(iq, 0) == doctest::Approx(d.QN(iq, 0)).epsilon(1e-4));
}

TEST_CASE("fluid closing: methods outside the closing family are refused by name") {
    // `solver_fluid_analyzer.m` routes exactly four names to
    // solver_fluid_closing. The others are different drifts or different
    // solvers, and answering for them under this name would report one method's
    // number as another's; they stay reachable through solver_fluid.
    qn::Network<double> m = closed_pair(2.0, lang::SchedStrategy::PS);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const char* elsewhere[] = {"matrix", "pnorm", "mfq", "diffusion", "rmf", "nonesuch"};
    for (std::size_t i = 0; i < sizeof(elsewhere) / sizeof(elsewhere[0]); ++i) {
        fluid::FluidOptions o;
        o.method = elsewhere[i];
        CHECK_THROWS_AS(fluid::solver_fluid_closing(sn, o), line::UnsupportedError);
    }

    // `default` at this entry point means the closing drift, not the analyzer's
    // default (which is the matrix method): the function names the family.
    fluid::FluidOptions dflt;
    const fluid::FluidSolution s = fluid::solver_fluid_closing(sn, dflt);
    CHECK(s.method == "closing");
    fluid::FluidOptions qualified;
    qualified.method = "fluid.closing";
    CHECK_NOTHROW(fluid::solver_fluid_closing(sn, qualified));
}
