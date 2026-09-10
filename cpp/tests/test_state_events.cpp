/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The event half of the `+State` port. The oracle is the encoding: each case
 * fixes an input row whose successor can be written down by hand, so a
 * disagreement is a decoding error rather than a numerical one.
 */
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/state_events.h"

namespace qn = line::qn;
using line::lang::EventType;
using line::lang::SchedStrategy;

namespace {

/** Source -> Queue -> Sink, one open class, the given discipline. */
qn::Network<double> open_queue(SchedStrategy s, double servers = 1.0) {
    qn::Network<double> m("ev");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", s);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, line::lang::Distrib<double>::exp_rate(0.5));
    m.set_service(q, c, line::lang::Distrib<double>::exp_rate(1.0));
    m.set_number_of_servers(q, servers);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("ARV at an idle FCFS server puts the job in service, not the buffer") {
    qn::Network<double> m = open_queue(SchedStrategy::FCFS);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t q = 2;  // Src, Q, Sink

    // Empty station, layout [buffer | server]: the arrival seizes the server.
    const std::vector<double> in{0.0, 0.0};
    const qn::EventOutcome<double> o = qn::after_event_station_arv(sn, q, in, 1);
    REQUIRE(o.space.size() == 1);
    CHECK(o.space[0][0] == doctest::Approx(0.0));
    CHECK(o.space[0][1] == doctest::Approx(1.0));
    // ARV is the PASSIVE half: the upstream departure sets the rate, so the
    // handler returns the -1 sentinel rather than inventing one.
    CHECK(o.rate[0] == doctest::Approx(-1.0));
    CHECK(o.prob[0] == doctest::Approx(1.0));
}

TEST_CASE("ARV at a busy FCFS server tags the waiting position with the class") {
    qn::Network<double> m = open_queue(SchedStrategy::FCFS);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // Server busy, one free buffer slot: the job waits, and the slot records
    // WHICH class waits there -- an ordered buffer, not a count.
    const std::vector<double> in{0.0, 1.0};
    const qn::EventOutcome<double> o = qn::after_event_station_arv(sn, 2, in, 1);
    REQUIRE(o.space.size() == 1);
    CHECK(o.space[0][0] == doctest::Approx(1.0));
    CHECK(o.space[0][1] == doctest::Approx(1.0));
}

TEST_CASE("ARV at a PS station never queues") {
    qn::Network<double> m = open_queue(SchedStrategy::PS);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // PS has no buffer at all: every job is in service, so a second arrival
    // increments the same phase slot rather than waiting.
    const std::vector<double> in{1.0};
    const qn::EventOutcome<double> o = qn::after_event_station_arv(sn, 2, in, 1);
    REQUIRE(o.space.size() == 1);
    CHECK(o.space[0].size() == 1);
    CHECK(o.space[0][0] == doctest::Approx(2.0));
}

TEST_CASE("ARV of an Erlang service splits over the entry phases by pie") {
    qn::Network<double> m("erl");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, line::lang::Distrib<double>::exp_rate(0.5));
    m.set_service(q, c, line::lang::Distrib<double>::erlang(0.5, 2));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // Two phases, so two candidate successors -- but an Erlang always STARTS
    // in phase 1, so pie = [1,0]. That is map_pie, the law embedded at
    // departure instants, and NOT map_prob, which for Erlang-2 is [0.5,0.5]:
    // using the latter would let a service begin mid-way through itself.
    const std::vector<double> in{0.0, 0.0};
    const qn::EventOutcome<double> o = qn::after_event_station_arv(sn, 2, in, 1);
    REQUIRE(o.space.size() == 2);
    CHECK(o.space[0][0] == doctest::Approx(1.0));
    CHECK(o.prob[0] == doctest::Approx(1.0));
    CHECK(o.space[1][1] == doctest::Approx(1.0));
    CHECK(o.prob[1] == doctest::Approx(0.0));
}

TEST_CASE("ARV at a SIRO station counts the waiting jobs per class") {
    qn::Network<double> m = open_queue(SchedStrategy::SIRO);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // SIRO stores a per-class COUNT, so the buffer column is the number of
    // class-1 jobs waiting -- there is no order to record, since the server
    // picks at random.
    const std::vector<double> in{0.0, 1.0};
    const qn::EventOutcome<double> o = qn::after_event_station_arv(sn, 2, in, 1);
    REQUIRE(o.space.size() == 1);
    CHECK(o.space[0][0] == doctest::Approx(1.0));
    CHECK(o.space[0][1] == doctest::Approx(1.0));
}

TEST_CASE("DEP at an FCFS station promotes the head of the buffer") {
    qn::Network<double> m = open_queue(SchedStrategy::FCFS);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // One in service, one waiting: the completion frees the server and the
    // waiting job takes it, so the buffer empties and the server stays busy.
    const std::vector<double> in{1.0, 1.0};
    const qn::EventOutcome<double> o = qn::after_event_station_dep(sn, 2, in, 1);
    REQUIRE(o.space.size() == 1);
    CHECK(o.space[0][0] == doctest::Approx(0.0));
    CHECK(o.space[0][1] == doctest::Approx(1.0));
    // DEP is the ACTIVE half, so it carries a real rate: D1(1,1) * kir = 1.
    CHECK(o.rate[0] == doctest::Approx(1.0));
}

TEST_CASE("DEP at an FCFS station with an empty buffer leaves the server idle") {
    qn::Network<double> m = open_queue(SchedStrategy::FCFS);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const std::vector<double> in{0.0, 1.0};
    const qn::EventOutcome<double> o = qn::after_event_station_dep(sn, 2, in, 1);
    REQUIRE(o.space.size() == 1);
    CHECK(o.space[0][1] == doctest::Approx(0.0));
    CHECK(o.rate[0] == doctest::Approx(1.0));
}

TEST_CASE("DEP at a PS station scales the rate by the server share") {
    qn::Network<double> m = open_queue(SchedStrategy::PS);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // Three jobs sharing one server: each completes at mu * (1/3) * min(3,1),
    // and there are three of them, so the total departure rate is still mu.
    // Aggregating over the class gives kir/ni * min(ni,S) = 3 * (1/3) * 1.
    const std::vector<double> in{3.0};
    const qn::EventOutcome<double> o = qn::after_event_station_dep(sn, 2, in, 1);
    REQUIRE(o.space.size() == 1);
    CHECK(o.space[0][0] == doctest::Approx(2.0));
    CHECK(o.rate[0] == doctest::Approx(1.0));
}

TEST_CASE("DEP at an INF station departs at the full phase rate") {
    qn::Network<double> m = open_queue(SchedStrategy::INF);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // No contention: every one of the 3 jobs completes at mu, so the aggregate
    // rate is 3 mu -- the difference from PS above is the whole point of INF.
    const std::vector<double> in{3.0};
    const qn::EventOutcome<double> o = qn::after_event_station_dep(sn, 2, in, 1);
    REQUIRE(o.space.size() == 1);
    CHECK(o.rate[0] == doctest::Approx(3.0));
}

TEST_CASE("DEP at a retrial station does not promote from the orbit") {
    qn::Network<double> m = open_queue(SchedStrategy::FCFS);
    const std::size_t q = 2;
    m.set_retrial(q, 1, line::lang::Distrib<double>::exp_rate(2.0), 2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // Same state as the promoting case above, but the waiting job is in an
    // ORBIT: a completion leaves the server IDLE and the orbit untouched. The
    // job re-enters only through a RETRY, which is what makes the station a
    // retrial queue rather than a queue with a differently named buffer.
    const std::vector<double> in{1.0, 1.0};
    const qn::EventOutcome<double> o = qn::after_event_station_dep(sn, q, in, 1);
    REQUIRE(o.space.size() == 1);
    CHECK(o.space[0][0] == doctest::Approx(1.0));
    CHECK(o.space[0][1] == doctest::Approx(0.0));
}

TEST_CASE("PHASE advances an Erlang service without departing") {
    qn::Network<double> m("erlph");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, line::lang::Distrib<double>::exp_rate(0.5));
    m.set_service(q, c, line::lang::Distrib<double>::erlang(2.0, 2));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // A job in phase 1 of an Erlang-2 advances to phase 2 at rate 2, with no
    // departure. The job count is unchanged: that is what separates PHASE from
    // DEP, and collapsing the two would make the service exponential.
    const std::vector<double> in{0.0, 1.0, 0.0};
    const qn::EventOutcome<double> o = qn::after_event_station_phase(sn, 2, in, 1);
    REQUIRE(o.space.size() == 1);
    CHECK(o.space[0][1] == doctest::Approx(0.0));
    CHECK(o.space[0][2] == doctest::Approx(1.0));
    CHECK(o.rate[0] == doctest::Approx(2.0));
}

TEST_CASE("PHASE under PS is slowed by the same share as a departure") {
    qn::Network<double> m("erlps");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, line::lang::Distrib<double>::exp_rate(0.5));
    m.set_service(q, c, line::lang::Distrib<double>::erlang(2.0, 2));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // Two jobs, both in phase 1, one server: each advances at 2 * (1/2) of the
    // capacity, so the aggregate is 2 * 2 * (1/2) = 2. Advancing at full speed
    // under contention would shorten the effective service time.
    const std::vector<double> in{2.0, 0.0};
    const qn::EventOutcome<double> o = qn::after_event_station_phase(sn, 2, in, 1);
    REQUIRE(o.space.size() == 1);
    CHECK(o.rate[0] == doctest::Approx(2.0));
}

TEST_CASE("RENEGE removes a waiting job, never one in service") {
    qn::Network<double> m = open_queue(SchedStrategy::FCFS);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // One in service, one waiting: only the waiting job is impatient, so the
    // rate is 1 * mu and the server is untouched.
    const std::vector<double> in{1.0, 1.0};
    const qn::EventOutcome<double> o = qn::after_event_station_renege(sn, 2, in, 1, 0.5);
    REQUIRE(o.space.size() == 1);
    CHECK(o.space[0][0] == doctest::Approx(0.0));
    CHECK(o.space[0][1] == doctest::Approx(1.0));
    CHECK(o.rate[0] == doctest::Approx(0.5));

    // With nobody waiting there is no one to renege, so the event does not fire.
    const std::vector<double> only_srv{0.0, 1.0};
    CHECK(qn::after_event_station_renege(sn, 2, only_srv, 1, 0.5).empty());
}

TEST_CASE("RETRY fires only when a server is free") {
    qn::Network<double> m = open_queue(SchedStrategy::FCFS);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // Server busy: the orbiting job retries and fails, so no transition is
    // generated at all. That refusal is what makes this an orbit rather than a
    // waiting line -- a queue would have promoted on the server's completion.
    const std::vector<double> busy{1.0, 1.0};
    CHECK(qn::after_event_station_retry(sn, 2, busy, 1, 3.0).empty());

    // Server free: the orbiting job enters service at (orbit size) * mu under
    // the default LINEAR policy, since every orbiting job carries its own timer.
    const std::vector<double> free_srv{1.0, 0.0};
    const qn::EventOutcome<double> o = qn::after_event_station_retry(sn, 2, free_srv, 1, 3.0);
    REQUIRE(o.space.size() == 1);
    CHECK(o.space[0][0] == doctest::Approx(0.0));
    CHECK(o.space[0][1] == doctest::Approx(1.0));
    CHECK(o.rate[0] == doctest::Approx(3.0));

    // CONSTANT policy: one controller retries for the whole orbit, so the rate
    // does not scale with the orbit size.
    const qn::EventOutcome<double> oc =
        qn::after_event_station_retry(sn, 2, free_srv, 1, 3.0, true);
    REQUIRE(oc.rate.size() == 1);
    CHECK(oc.rate[0] == doctest::Approx(3.0));
}

namespace {

/** Source -> Queue -> Sink with two open classes at distinct priorities. */
qn::Network<double> two_class_prio(SchedStrategy s) {
    qn::Network<double> m("prio");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", s);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("Hi", 0);   // lower value = more urgent
    const std::size_t c2 = m.add_open_class("Lo", 1);
    m.set_arrival(src, c1, line::lang::Distrib<double>::exp_rate(0.2));
    m.set_arrival(src, c2, line::lang::Distrib<double>::exp_rate(0.2));
    m.set_service(q, c1, line::lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q, c2, line::lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, k, 1.0);
    P.set(c2, c2, src, q, 1.0);
    P.set(c2, c2, q, k, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("DEP at HOL promotes by priority, not by arrival order") {
    qn::Network<double> m = two_class_prio(SchedStrategy::HOL);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // Buffer [Lo, Hi] with a class-1 job in service. Arrivals fill from the
    // right, so Lo is the OLDER job -- plain FCFS would promote it. HOL serves
    // the more urgent class first, so the Hi job enters service instead.
    const std::vector<double> in{2.0, 1.0, 1.0, 0.0};
    const qn::EventOutcome<double> o = qn::after_event_station_dep(sn, 2, in, 1);
    REQUIRE(o.space.size() == 1);
    CHECK(o.space[0][2] == doctest::Approx(1.0));  // class 1 back in service
    CHECK(o.space[0][3] == doctest::Approx(0.0));
    // The gap closes and Lo shifts right, so the buffer stays right-aligned.
    CHECK(o.space[0][0] == doctest::Approx(0.0));
    CHECK(o.space[0][1] == doctest::Approx(2.0));
    CHECK(o.rate[0] == doctest::Approx(1.0));
}

TEST_CASE("DEP at LCFS promotes the newest arrival and ignores priority") {
    qn::Network<double> m = two_class_prio(SchedStrategy::LCFS);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // Same buffer [Lo, Hi]. LCFS is NOT priority-aware: it takes the most
    // recent arrival, which is the first nonzero column -- the Lo job. Reading
    // the priorities here would silently turn LCFS into LCFSPRIO.
    const std::vector<double> in{2.0, 1.0, 1.0, 0.0};
    const qn::EventOutcome<double> o = qn::after_event_station_dep(sn, 2, in, 1);
    REQUIRE(o.space.size() == 1);
    CHECK(o.space[0][3] == doctest::Approx(1.0));  // class 2 entered service
    CHECK(o.space[0][2] == doctest::Approx(0.0));
    // LCFS clears the slot in place rather than closing the gap.
    CHECK(o.space[0][0] == doctest::Approx(0.0));
    CHECK(o.space[0][1] == doctest::Approx(1.0));
}

TEST_CASE("FAILURE and REPAIR move only the status column") {
    qn::Network<double> m = open_queue(SchedStrategy::FCFS);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // Trailing column is the server status. A job in service is NOT lost when
    // the server fails: memoryless service resumes on repair with nothing to
    // remember, so only the status flips.
    const std::vector<double> up{0.0, 1.0, 1.0};
    const qn::EventOutcome<double> f =
        qn::after_event_station_breakdown(sn, 2, up, false, 0.1);
    REQUIRE(f.space.size() == 1);
    CHECK(f.space[0][1] == doctest::Approx(1.0));  // the job is still in service
    CHECK(f.space[0][2] == doctest::Approx(0.0));
    CHECK(f.rate[0] == doctest::Approx(0.1));

    // A failure needs an UP server: firing it on a down one is not an
    // admissible transition and must generate no edge.
    const std::vector<double> down{0.0, 1.0, 0.0};
    CHECK(qn::after_event_station_breakdown(sn, 2, down, false, 0.1).empty());

    const qn::EventOutcome<double> r =
        qn::after_event_station_breakdown(sn, 2, down, true, 0.5);
    REQUIRE(r.space.size() == 1);
    CHECK(r.space[0][2] == doctest::Approx(1.0));
    CHECK(r.rate[0] == doctest::Approx(0.5));
}

TEST_CASE("after_event_station dispatches to the per-event handlers") {
    qn::Network<double> m = open_queue(SchedStrategy::FCFS);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<double> in{0.0, 1.0};

    const qn::EventOutcome<double> d =
        qn::after_event_station(sn, 2, in, EventType::DEP, 1);
    REQUIRE(d.space.size() == 1);
    CHECK(d.rate[0] == doctest::Approx(1.0));

    const qn::EventOutcome<double> a =
        qn::after_event_station(sn, 2, in, EventType::ARV, 1);
    REQUIRE(a.space.size() == 1);
    CHECK(a.rate[0] == doctest::Approx(-1.0));

    // A LOCAL event is a dummy: it moves nothing anywhere.
    CHECK(qn::after_event_station(sn, 2, in, EventType::LOCAL, 1).empty());
}

TEST_CASE("PSPRIO serves only the most urgent class present once saturated") {
    qn::Network<double> m = two_class_prio(SchedStrategy::PSPRIO);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // Two jobs, one server, so the station is saturated (ni > S). Class 1 is
    // the more urgent group and takes the whole server; class 2 gets nothing.
    const std::vector<double> in{1.0, 1.0};
    const qn::EventOutcome<double> hi = qn::after_event_station_dep(sn, 2, in, 1);
    REQUIRE(hi.space.size() == 1);
    CHECK(hi.rate[0] == doctest::Approx(1.0));

    const qn::EventOutcome<double> lo = qn::after_event_station_dep(sn, 2, in, 2);
    REQUIRE(lo.space.size() == 1);
    // The state is still emitted, with rate zero: a lower-priority class is
    // not served at all while a more urgent one is present.
    CHECK(lo.rate[0] == doctest::Approx(0.0));
}

TEST_CASE("PSPRIO below saturation shares like plain PS") {
    qn::Network<double> m = two_class_prio(SchedStrategy::PSPRIO);
    m.set_number_of_servers(2, 2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // Two jobs, two servers: nobody waits, so precedence is moot and both
    // classes are served at the full rate.
    const std::vector<double> in{1.0, 1.0};
    CHECK(qn::after_event_station_dep(sn, 2, in, 1).rate[0] == doctest::Approx(1.0));
    CHECK(qn::after_event_station_dep(sn, 2, in, 2).rate[0] == doctest::Approx(1.0));
}

TEST_CASE("refresh_sync pairs every departure with its arrival") {
    qn::Network<double> m = open_queue(SchedStrategy::FCFS);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<qn::Sync<double>> sy = qn::refresh_sync(sn);

    // Source -> Queue -> Sink, one class, exponential service: two DEP/ARV
    // pairs and no phase-change action, since a single phase has no internal
    // transition to make.
    REQUIRE(sy.size() == 2);
    for (std::size_t i = 0; i < sy.size(); ++i) {
        CHECK(sy[i].active.event == EventType::DEP);
        CHECK(sy[i].passive.event == EventType::ARV);
        CHECK(sy[i].passive.prob == doctest::Approx(1.0));
    }
    CHECK(sy[0].active.node == 1);   // Source departs
    CHECK(sy[0].passive.node == 2);  // into the Queue
    CHECK(sy[1].active.node == 2);   // Queue departs
    // ...back into the SOURCE, not into the Sink. A Sink is not stateful, so
    // the stochastic complement that builds rt eliminates it and closes the
    // open chain at the Source. That is what the EXT arrival branch means by
    // accepting a "virtual arrival from the sink": the departure has to
    // synchronize with something, and the eliminated Sink is not available.
    CHECK(sy[1].passive.node == 1);
}

TEST_CASE("refresh_sync emits a phase action only for a multi-phase service") {
    qn::Network<double> m("syncerl");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, line::lang::Distrib<double>::exp_rate(0.5));
    m.set_service(q, c, line::lang::Distrib<double>::erlang(2.0, 2));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<qn::Sync<double>> sy = qn::refresh_sync(sn);

    // The Erlang-2 queue gains a PHASE action; the exponential source does not.
    std::size_t nphase = 0;
    for (std::size_t i = 0; i < sy.size(); ++i)
        if (sy[i].active.event == EventType::PHASE) {
            ++nphase;
            CHECK(sy[i].active.node == q);
            // A phase change moves no job, so its passive half is LOCAL.
            CHECK(sy[i].passive.event == EventType::LOCAL);
        }
    CHECK(nphase == 1);
}

TEST_CASE("refresh_sync gives a retrial station a RETRY action") {
    qn::Network<double> m = open_queue(SchedStrategy::FCFS);
    m.set_retrial(2, 1, line::lang::Distrib<double>::exp_rate(2.0), 2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<qn::Sync<double>> sy = qn::refresh_sync(sn);

    std::size_t nretry = 0;
    for (std::size_t i = 0; i < sy.size(); ++i)
        if (sy[i].active.event == EventType::RETRY) {
            ++nretry;
            CHECK(sy[i].passive.event == EventType::LOCAL);
        }
    CHECK(nretry == 1);
}

namespace {

/** A two-class FCFS station whose second class is a G-network signal. */
qn::Network<double> signal_model(line::lang::SignalType ty,
                                 line::lang::RemovalPolicy pol) {
    qn::Network<double> m("gnet");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("Job");
    const std::size_t c2 = m.add_open_class("Sig");
    m.set_arrival(src, c1, line::lang::Distrib<double>::exp_rate(0.4));
    m.set_arrival(src, c2, line::lang::Distrib<double>::exp_rate(0.2));
    m.set_service(q, c1, line::lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q, c2, line::lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, k, 1.0);
    P.set(c2, c2, src, q, 1.0);
    P.set(c2, c2, q, k, 1.0);
    m.link(P);
    m.set_signal(c2, ty, pol);
    return m;
}

}  // namespace

TEST_CASE("a catastrophe signal empties the station outright") {
    qn::Network<double> m =
        signal_model(line::lang::SignalType::CATASTROPHE, line::lang::RemovalPolicy::RANDOM);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // Two class-1 jobs, one in service and one waiting: a catastrophe takes
    // every job regardless of the batch distribution, by definition.
    const std::vector<double> in{1.0, 1.0, 0.0};
    const qn::EventOutcome<double> o =
        qn::after_event_station_signal(sn, 2, in, 2);
    REQUIRE(o.space.size() == 1);
    for (std::size_t j = 0; j < o.space[0].size(); ++j)
        CHECK(o.space[0][j] == doctest::Approx(0.0));
    // The signal is passive: the source that emitted it sets the rate.
    CHECK(o.rate[0] == doctest::Approx(-1.0));
}

TEST_CASE("a negative signal with no victim vanishes without changing the state") {
    qn::Network<double> m =
        signal_model(line::lang::SignalType::NEGATIVE, line::lang::RemovalPolicy::RANDOM);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const std::vector<double> empty{0.0, 0.0, 0.0};
    const qn::EventOutcome<double> o = qn::after_event_station_signal(sn, 2, empty, 2);
    REQUIRE(o.space.size() == 1);
    CHECK(o.space[0] == empty);
    CHECK(o.prob[0] == doctest::Approx(1.0));
}

TEST_CASE("an FCFS-policy signal removes the oldest waiting job") {
    qn::Network<double> m =
        signal_model(line::lang::SignalType::NEGATIVE, line::lang::RemovalPolicy::FCFS);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // One class-1 job in service, one waiting. FCFS drains the waiting line
    // before touching a server, so the server keeps its job.
    const std::vector<double> in{1.0, 1.0, 0.0};
    const qn::EventOutcome<double> o = qn::after_event_station_signal(sn, 2, in, 2);
    REQUIRE(o.space.size() == 1);
    CHECK(o.space[0][0] == doctest::Approx(0.0));  // the waiting job is gone
    CHECK(o.space[0][1] == doctest::Approx(1.0));  // the server is untouched
}

/** As signal_model, but the job class has a two-phase Erlang service. */
static qn::Network<double> signal_model_erlang(line::lang::RemovalPolicy pol) {
    qn::Network<double> m("gnet2");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("Job");
    const std::size_t c2 = m.add_open_class("Sig");
    m.set_arrival(src, c1, line::lang::Distrib<double>::exp_rate(0.4));
    m.set_arrival(src, c2, line::lang::Distrib<double>::exp_rate(0.2));
    m.set_service(q, c1, line::lang::Distrib<double>::erlang(2.0, 2));
    m.set_service(q, c2, line::lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, k, 1.0);
    P.set(c2, c2, src, q, 1.0);
    P.set(c2, c2, q, k, 1.0);
    m.link(P);
    m.set_signal(c2, line::lang::SignalType::NEGATIVE, pol);
    return m;
}

TEST_CASE("a RANDOM-policy signal draws over waiting and in-service jobs alike") {
    qn::Network<double> m = signal_model_erlang(line::lang::RemovalPolicy::RANDOM);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // Layout [buf | job phase 1, job phase 2, signal phase]: one job waiting,
    // one job in service in PHASE 2. RANDOM is uniform over both tiers, and
    // the two outcomes are now distinguishable -- dropping the waiting job
    // leaves the server in phase 2, whereas dropping the served job promotes
    // the waiting one, which restarts in phase 1.
    const std::vector<double> in{1.0, 0.0, 1.0, 0.0};
    const qn::EventOutcome<double> o = qn::after_event_station_signal(sn, 2, in, 2);
    REQUIRE(o.space.size() == 2);
    double tot = 0;
    for (std::size_t i = 0; i < o.prob.size(); ++i) {
        tot += o.prob[i];
        CHECK(o.prob[i] == doctest::Approx(0.5));
    }
    CHECK(tot == doctest::Approx(1.0));

    // An age policy reaches the servers only once the waiting line is empty,
    // so it has a single successor here where RANDOM has two.
    qn::Network<double> mf = signal_model_erlang(line::lang::RemovalPolicy::FCFS);
    const qn::EventOutcome<double> of =
        qn::after_event_station_signal(mf.get_struct(), 2, in, 2);
    REQUIRE(of.space.size() == 1);
    CHECK(of.space[0][0] == doctest::Approx(0.0));  // the waiting job went
    CHECK(of.space[0][2] == doctest::Approx(1.0));  // the server kept phase 2
}

TEST_CASE("indistinguishable signal victims are merged into one successor") {
    qn::Network<double> m =
        signal_model(line::lang::SignalType::NEGATIVE, line::lang::RemovalPolicy::RANDOM);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // One class-1 job waiting and one class-1 job in service. Removing the
    // waiting job, or removing the served one and promoting the waiting job,
    // reach the SAME state -- jobs of a class are exchangeable. The two
    // branches must merge into a single successor of probability 1, or the
    // generator would carry the same edge twice.
    const std::vector<double> in{1.0, 1.0, 0.0};
    const qn::EventOutcome<double> o = qn::after_event_station_signal(sn, 2, in, 2);
    REQUIRE(o.space.size() == 1);
    CHECK(o.prob[0] == doctest::Approx(1.0));
    CHECK(o.space[0][0] == doctest::Approx(0.0));
    CHECK(o.space[0][1] == doctest::Approx(1.0));
}

TEST_CASE("a signal arriving is dispatched to the signal handler, not to ARV") {
    qn::Network<double> m =
        signal_model(line::lang::SignalType::CATASTROPHE, line::lang::RemovalPolicy::RANDOM);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // A signal never JOINS: routing it through the arrival branch would add a
    // job instead of removing every one of them.
    const std::vector<double> in{1.0, 1.0, 0.0};
    const qn::EventOutcome<double> o =
        qn::after_event(sn, 2, in, EventType::ARV, 2);
    REQUIRE(o.space.size() == 1);
    for (std::size_t j = 0; j < o.space[0].size(); ++j)
        CHECK(o.space[0][j] == doctest::Approx(0.0));
}

TEST_CASE("pass_and_swap ejects the end of the swap chain, not the completing job") {
    // Two classes, class 1 swappable with class 2 only.
    std::vector<std::vector<bool>> G(2, std::vector<bool>(2, false));
    G[0][1] = true;
    G[1][0] = true;

    // List (1, 2): the class-1 job at position 0 completes. It finds the
    // class-2 job swappable, takes its place, and the class-2 job -- having no
    // successor -- departs. The departing class is NOT the one that completed,
    // which is the whole point of the mechanism.
    const std::vector<std::size_t> c{1, 2};
    const std::pair<std::vector<std::size_t>, std::size_t> r =
        qn::pass_and_swap<double>(c, 0, G);
    CHECK(r.second == 2);
    REQUIRE(r.first.size() == 1);
    CHECK(r.first[0] == 1);

    // With no swappable successor the completing job departs itself.
    std::vector<std::vector<bool>> none(2, std::vector<bool>(2, false));
    const std::pair<std::vector<std::size_t>, std::size_t> r2 =
        qn::pass_and_swap<double>(c, 0, none);
    CHECK(r2.second == 1);
    REQUIRE(r2.first.size() == 1);
    CHECK(r2.first[0] == 2);
}

TEST_CASE("a PAS station serves by position increments of the rate function") {
    qn::Network<double> m("pas");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PAS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, line::lang::Distrib<double>::exp_rate(0.4));
    m.set_service(q, c, line::lang::Distrib<double>::exp_rate(1.0));
    m.set_capacity(q, 3);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    // mu(c) = |c|, so every position contributes an increment of exactly 1 and
    // the station behaves as an infinite server.
    m.set_pas(q, [](const std::vector<std::size_t>& lst) {
        return static_cast<double>(lst.size());
    });
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // The state is the ordered LIST, not a buffer/server split: two class-1
    // jobs, then a zero pad.
    const std::vector<double> in{1.0, 1.0, 0.0};
    const qn::EventOutcome<double> d = qn::after_event_station_dep(sn, 2, in, 1);
    // Both positions can complete, and with an empty swap graph each ejects
    // its own job, so there are two successors of rate 1 each.
    REQUIRE(d.space.size() == 2);
    for (std::size_t i = 0; i < d.rate.size(); ++i) CHECK(d.rate[i] == doctest::Approx(1.0));
    for (std::size_t i = 0; i < d.space.size(); ++i) {
        CHECK(d.space[i][0] == doctest::Approx(1.0));
        CHECK(d.space[i][1] == doctest::Approx(0.0));
    }

    // An arrival joins at the BACK, which is what the rate function's order
    // dependence reads.
    const qn::EventOutcome<double> a = qn::after_event_station_arv(sn, 2, in, 1);
    REQUIRE(a.space.size() == 1);
    CHECK(a.space[0][2] == doctest::Approx(1.0));
    CHECK(a.rate[0] == doctest::Approx(-1.0));

    // A full list loses the arrival rather than growing the state.
    const std::vector<double> full{1.0, 1.0, 1.0};
    CHECK(qn::after_event_station_arv(sn, 2, full, 1).empty());
}

namespace {

/** An FCFS station that holds a server across a synchronous call. */
qn::Network<double> reply_model() {
    qn::Network<double> m("sync");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t call = m.add_open_class("Call");
    const std::size_t rep = m.add_open_class("Reply");
    m.set_arrival(src, call, line::lang::Distrib<double>::exp_rate(0.3));
    m.set_arrival(src, rep, line::lang::Distrib<double>::exp_rate(0.3));
    m.set_service(q, call, line::lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q, rep, line::lang::Distrib<double>::exp_rate(1.0));
    m.set_number_of_servers(q, 1);
    qn::RoutingMatrix<double> P;
    P.set(call, call, src, q, 1.0);
    P.set(call, call, q, k, 1.0);
    P.set(rep, rep, src, q, 1.0);
    P.set(rep, rep, q, k, 1.0);
    m.link(P);
    m.set_sync_reply(q, call, rep);
    return m;
}

}  // namespace

TEST_CASE("a synchronous call keeps its server on departure") {
    qn::Network<double> m = reply_model();
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // Layout [buf | call phase, reply phase | reply-block counter]. One call
    // job in service, nothing waiting.
    const std::vector<double> in{0.0, 1.0, 0.0, 0.0};
    const qn::EventOutcome<double> d = qn::after_event_station_dep(sn, 2, in, 1);
    REQUIRE(d.space.size() == 1);
    // The job leaves for the callee but the server stays held: the counter
    // goes up, and the server is NOT free for anyone else.
    CHECK(d.space[0][1] == doctest::Approx(0.0));
    CHECK(d.space[0][3] == doctest::Approx(1.0));
}

TEST_CASE("a held server is unavailable to an arriving job") {
    qn::Network<double> m = reply_model();
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // Servers idle by raw occupancy, but one is held for a pending reply, so
    // an arriving call job must WAIT rather than seize it.
    const std::vector<double> held{0.0, 0.0, 0.0, 1.0};
    const qn::EventOutcome<double> a = qn::after_event_station_arv(sn, 2, held, 1);
    REQUIRE(a.space.size() == 1);
    CHECK(a.space[0][0] == doctest::Approx(1.0));  // queued, not served
    CHECK(a.space[0][1] == doctest::Approx(0.0));
}

TEST_CASE("a REPLY releases the held server and takes it itself") {
    qn::Network<double> m = reply_model();
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // One server held for a pending reply, one call job waiting in the buffer.
    // The reply releases the block and PASSES THROUGH into the freed server:
    // it is work the station already paid for, so it does not queue behind the
    // resident, which would misreport its residence and steal capacity.
    const std::vector<double> in{1.0, 0.0, 0.0, 1.0};
    const qn::EventOutcome<double> o = qn::after_event(sn, 2, in, EventType::ARV, 2);
    REQUIRE(o.space.size() == 1);
    CHECK(o.space[0][3] == doctest::Approx(0.0));  // the block is released
    CHECK(o.space[0][2] == doctest::Approx(1.0));  // the reply is in service
    CHECK(o.space[0][0] == doctest::Approx(1.0));  // the waiting job still waits
    CHECK(o.rate[0] == doctest::Approx(-1.0));     // passive: the callee sets it
}

namespace {

/** A two-buffer polling station with a real switchover between the buffers. */
qn::Network<double> polling_model(line::lang::PollingType ty, bool immediate_sw) {
    qn::Network<double> m("poll");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::POLLING);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("A");
    const std::size_t c2 = m.add_open_class("B");
    m.set_arrival(src, c1, line::lang::Distrib<double>::exp_rate(0.2));
    m.set_arrival(src, c2, line::lang::Distrib<double>::exp_rate(0.2));
    m.set_service(q, c1, line::lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q, c2, line::lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, k, 1.0);
    P.set(c2, c2, src, q, 1.0);
    P.set(c2, c2, q, k, 1.0);
    m.link(P);
    std::vector<line::lang::Distrib<double>> sw(2);
    if (immediate_sw) {
        sw[0] = line::lang::Distrib<double>::immediate();
        sw[1] = line::lang::Distrib<double>::immediate();
    } else {
        sw[0] = line::lang::Distrib<double>::exp_rate(4.0);
        sw[1] = line::lang::Distrib<double>::exp_rate(4.0);
    }
    m.set_polling(q, ty, sw);
    return m;
}

}  // namespace

TEST_CASE("an immediate switchover costs no controller state") {
    qn::Network<double> m = polling_model(line::lang::PollingType::EXHAUSTIVE, true);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const qn::PollingInfo<double> pi = qn::polling_info(sn, 2);
    REQUIRE(pi.valid);

    // Every switchover immediate and EXHAUSTIVE: pos and swk carry no
    // information and there is no visit budget, so the controller is
    // zero-width. Representing an Immediate walk as a state would put a ~1e8
    // rate in the generator and add a spurious state per buffer.
    CHECK(pi.width == 0);
    for (std::size_t r = 0; r < 2; ++r) CHECK(!pi.has_sw[r]);
}

TEST_CASE("a real switchover materializes the position and phase columns") {
    qn::Network<double> m = polling_model(line::lang::PollingType::EXHAUSTIVE, false);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const qn::PollingInfo<double> pi = qn::polling_info(sn, 2);
    REQUIRE(pi.valid);

    // Non-immediate walks make SWITCHING(p) a distinguishable configuration,
    // so pos and swk exist; EXHAUSTIVE still needs no budget counter.
    CHECK(pi.width == 2);
    CHECK(pi.ipos == 0);
    CHECK(pi.iswk == 1);
    CHECK(pi.ictr == static_cast<std::size_t>(-1));
    for (std::size_t r = 0; r < 2; ++r) CHECK(pi.has_sw[r]);
}

TEST_CASE("a bounded polling discipline materializes the visit budget") {
    qn::Network<double> m = polling_model(line::lang::PollingType::KLIMITED, true);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const qn::PollingInfo<double> pi = qn::polling_info(sn, 2);
    REQUIRE(pi.valid);

    // Immediate walks, so no pos/swk, but K-LIMITED bounds the visit and needs
    // the counter -- the only column that distinguishes it from EXHAUSTIVE.
    CHECK(pi.width == 1);
    CHECK(pi.ipos == static_cast<std::size_t>(-1));
    CHECK(pi.ictr == 0);

    // The budget is the discipline's own bound, not the buffer contents.
    CHECK(qn::polling_budget(pi, 7) == 1);
    qn::Network<double> mg = polling_model(line::lang::PollingType::GATED, true);
    const qn::PollingInfo<double> pg = qn::polling_info(mg.get_struct(), 2);
    CHECK(qn::polling_budget(pg, 7) == 7);   // exactly those found
    qn::Network<double> md = polling_model(line::lang::PollingType::DECREMENTING, true);
    const qn::PollingInfo<double> pd = qn::polling_info(md.get_struct(), 2);
    CHECK(qn::polling_budget(pd, 7) == 6);   // one below the level found
}

TEST_CASE("polling_next walks the cyclic order and parks when there is no work") {
    qn::Network<double> m = polling_model(line::lang::PollingType::EXHAUSTIVE, true);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const qn::PollingInfo<double> pi = qn::polling_info(sn, 2);

    // Work waiting at buffer 2: the server walks there and opens a visit.
    std::vector<long> nbuf{0, 3};
    std::size_t q = 0;
    int mode = 0;
    long budget = 0;
    qn::polling_next(pi, 1, nbuf, 2, false, q, mode, budget);
    CHECK(q == 2);
    CHECK(mode == 1);

    // Nothing anywhere: with immediate switchovers a server would otherwise
    // lap forever in zero time, so it PARKS instead.
    std::vector<long> empty{0, 0};
    qn::polling_next(pi, 1, empty, 2, false, q, mode, budget);
    CHECK(mode == 0);

    // An arrival at the buffer the server already stands at opens a visit
    // there: that switchover has already been paid for.
    std::vector<long> here{2, 0};
    qn::polling_next(pi, 1, here, 2, true, q, mode, budget);
    CHECK(q == 1);
    CHECK(mode == 1);
}

TEST_CASE("refresh_sync emits a SWITCH action only for a real switchover") {
    qn::Network<double> m = polling_model(line::lang::PollingType::EXHAUSTIVE, false);
    std::size_t nsw = 0;
    const std::vector<qn::Sync<double>> sy = qn::refresh_sync(m.get_struct());
    for (std::size_t i = 0; i < sy.size(); ++i)
        if (sy[i].active.event == EventType::SWITCH) ++nsw;
    CHECK(nsw == 2);  // one per buffer with a non-immediate entering walk

    qn::Network<double> mi = polling_model(line::lang::PollingType::EXHAUSTIVE, true);
    std::size_t nsw2 = 0;
    const std::vector<qn::Sync<double>> sy2 = qn::refresh_sync(mi.get_struct());
    for (std::size_t i = 0; i < sy2.size(); ++i)
        if (sy2[i].active.event == EventType::SWITCH) ++nsw2;
    CHECK(nsw2 == 0);  // folded walks are never a transition
}

namespace {

/** A Source -> Cache -> Sink model with one list of capacity 2 over 3 items. */
qn::Network<double> cache_model(line::lang::ReplacementStrategy rs) {
    qn::Network<double> m("cache");
    const std::size_t src = m.add_source("Src");
    const std::size_t k = m.add_sink("Sink");
    const std::size_t rd = m.add_open_class("Read");
    const std::size_t hit = m.add_open_class("Hit");
    const std::size_t mis = m.add_open_class("Miss");
    m.set_arrival(src, rd, line::lang::Distrib<double>::exp_rate(1.0));
    qn::CacheParam<double> cpar;
    cpar.nitems = 3;
    cpar.itemcap.push_back(2);
    cpar.replacestrat = rs;
    cpar.pread.assign(3, std::vector<double>());
    cpar.pread[rd - 1] = std::vector<double>{0.5, 0.3, 0.2};
    cpar.hitclass.assign(3, 0);
    cpar.missclass.assign(3, 0);
    cpar.hitclass[rd - 1] = hit;
    cpar.missclass[rd - 1] = mis;
    m.add_cache("C", cpar);
    return m;
}

}  // namespace

TEST_CASE("a cache READ on a hit moves the job to the hit class") {
    qn::Network<double> m = cache_model(line::lang::ReplacementStrategy::LRU);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t cnode = 3;  // Src, Sink, Cache

    // [read, hit, miss | slot1, slot2]: one read job present, item 1 cached in
    // slot 1 and item 2 in slot 2. Reading item 1 is a HIT.
    const std::vector<double> in{1.0, 0.0, 0.0, 1.0, 2.0};
    const qn::EventOutcome<double> o =
        qn::after_event_cache(sn, cnode, in, EventType::READ, 1);

    // One successor per item with positive read probability: items 1 and 2 hit,
    // item 3 misses.
    REQUIRE(o.space.size() == 3);
    // Item 1 hit: leaves in the hit class, already at the head so nothing moves.
    CHECK(o.space[0][0] == doctest::Approx(0.0));
    CHECK(o.space[0][1] == doctest::Approx(1.0));
    CHECK(o.space[0][3] == doctest::Approx(1.0));
    // Item 2 hit under LRU: promoted to the head, item 1 shifts down.
    CHECK(o.space[1][1] == doctest::Approx(1.0));
    CHECK(o.space[1][3] == doctest::Approx(2.0));
    CHECK(o.space[1][4] == doctest::Approx(1.0));
    // Item 3 misses: it leaves in the miss class and is admitted at the head.
    CHECK(o.space[2][2] == doctest::Approx(1.0));
    CHECK(o.space[2][3] == doctest::Approx(3.0));
}

TEST_CASE("FIFO does not reorder on a hit, LRU does") {
    const std::vector<double> in{1.0, 0.0, 0.0, 1.0, 2.0};

    qn::Network<double> mf = cache_model(line::lang::ReplacementStrategy::FIFO);
    const qn::EventOutcome<double> f =
        qn::after_event_cache(mf.get_struct(), 3, in, EventType::READ, 1);
    // A hit on item 2 leaves the list untouched: FIFO orders by INSERTION, so
    // a read carries no information about eviction order.
    REQUIRE(f.space.size() == 3);
    CHECK(f.space[1][3] == doctest::Approx(1.0));
    CHECK(f.space[1][4] == doctest::Approx(2.0));

    qn::Network<double> ml = cache_model(line::lang::ReplacementStrategy::LRU);
    const qn::EventOutcome<double> l =
        qn::after_event_cache(ml.get_struct(), 3, in, EventType::READ, 1);
    REQUIRE(l.space.size() == 3);
    CHECK(l.space[1][3] == doctest::Approx(2.0));  // promoted to the head
}

TEST_CASE("a cache read is instantaneous and weighted by the item popularity") {
    qn::Network<double> m = cache_model(line::lang::ReplacementStrategy::LRU);
    const std::vector<double> in{1.0, 0.0, 0.0, 1.0, 2.0};
    const qn::EventOutcome<double> o =
        qn::after_event_cache(m.get_struct(), 3, in, EventType::READ, 1);

    // A read is a routing decision, not a service: every branch fires at the
    // Immediate rate, scaled by the probability of reading that item.
    const double imm = line::lang::GlobalConstants::Immediate;
    REQUIRE(o.rate.size() == 3);
    CHECK(o.rate[0] == doctest::Approx(0.5 * imm));
    CHECK(o.rate[1] == doctest::Approx(0.3 * imm));
    CHECK(o.rate[2] == doctest::Approx(0.2 * imm));
}

TEST_CASE("a cache READ needs exactly one job, of the reading class") {
    qn::Network<double> m = cache_model(line::lang::ReplacementStrategy::LRU);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // No job present: nothing to read.
    const std::vector<double> none{0.0, 0.0, 0.0, 1.0, 2.0};
    CHECK(qn::after_event_cache(sn, 3, none, EventType::READ, 1).empty());

    // Two jobs present: the read is defined for a single in-flight request, so
    // the enumeration must not produce a transition here.
    const std::vector<double> two{1.0, 1.0, 0.0, 1.0, 2.0};
    CHECK(qn::after_event_cache(sn, 3, two, EventType::READ, 1).empty());
}

TEST_CASE("a cache arrival is passive and a departure is immediate") {
    qn::Network<double> m = cache_model(line::lang::ReplacementStrategy::LRU);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const std::vector<double> in{0.0, 0.0, 0.0, 1.0, 2.0};
    const qn::EventOutcome<double> a =
        qn::after_event_cache(sn, 3, in, EventType::ARV, 1);
    REQUIRE(a.space.size() == 1);
    CHECK(a.space[0][0] == doctest::Approx(1.0));
    CHECK(a.rate[0] == doctest::Approx(-1.0));  // the upstream sets the rate

    const std::vector<double> held{0.0, 1.0, 0.0, 1.0, 2.0};
    const qn::EventOutcome<double> d =
        qn::after_event_cache(sn, 3, held, EventType::DEP, 2);
    REQUIRE(d.space.size() == 1);
    CHECK(d.space[0][1] == doctest::Approx(0.0));
    CHECK(d.rate[0] == doctest::Approx(line::lang::GlobalConstants::Immediate));
}
