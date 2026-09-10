/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * `State.toMarginal` port. The oracle is the encoding itself: each case below
 * fixes a state row whose marginal counts can be read off by hand, so a
 * disagreement is a decoding error and not a numerical one.
 */
#include <cmath>
#include <limits>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/state.h"

using line::qn::Marginal;
using line::qn::to_marginal;
namespace qn = line::qn;
using line::lang::SchedStrategy;

namespace {

/** Source -> FCFS Queue -> Sink, one open class. */
qn::Network<double> open_fcfs() {
    qn::Network<double> m("st");
    const std::size_t s = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(s, c, line::lang::Distrib<double>::exp_rate(0.5));
    m.set_service(q, c, line::lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/** Source -> PS Queue -> Sink, one open class (no buffer in the encoding). */
qn::Network<double> open_ps() {
    qn::Network<double> m("stps");
    const std::size_t s = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(s, c, line::lang::Distrib<double>::exp_rate(0.5));
    m.set_service(q, c, line::lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/** Source -> FCFS Queue -> Sink, TWO open classes. */
qn::Network<double> fcfs_two_class() {
    qn::Network<double> m("st2");
    const std::size_t s = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("C1");
    const std::size_t c2 = m.add_open_class("C2");
    m.set_arrival(s, c1, line::lang::Distrib<double>::exp_rate(0.5));
    m.set_arrival(s, c2, line::lang::Distrib<double>::exp_rate(0.5));
    m.set_service(q, c1, line::lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q, c2, line::lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, s, q, 1.0); P.set(c1, c1, q, k, 1.0);
    P.set(c2, c2, s, q, 1.0); P.set(c2, c2, q, k, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("to_marginal decodes an FCFS buffer by class tag") {
    qn::Network<double> m = open_fcfs();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    // Station 2 is the FCFS queue; one class, one phase, so the row is
    // [buffer tags | one server slot]. Three tagged waiting jobs plus one in
    // service must read as nir = 4, sir = 1.
    const std::vector<double> row{1, 1, 1, 1};
    const std::vector<std::size_t> phasesz{1}, phaseshift{0};
    const Marginal<double> g = to_marginal(sn, 2, row, phasesz, phaseshift);
    CHECK(g.sir[0] == 1.0);
    CHECK(g.nir[0] == 4.0);
    CHECK(g.ni == 4.0);
    CHECK(g.kir[0][0] == 1.0);
}

TEST_CASE("to_marginal returns the EXT sentinel at a Source") {
    qn::Network<double> m = open_fcfs();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    // Station 1 is the Source. nir is +Inf BY DESIGN: a Source is an infinite
    // reservoir and the value describes the encoding. Reading it as a queue
    // length is what made MATLAB's CTMC report Q = Inf and R = Q/T = Inf at a
    // Source, so this asserts the sentinel is present rather than clamped --
    // any consumer must branch on the station being a Source instead.
    const std::vector<double> row{1};
    const std::vector<std::size_t> phasesz{1}, phaseshift{0};
    const Marginal<double> g = to_marginal(sn, 1, row, phasesz, phaseshift);
    CHECK(std::isinf(g.nir[0]));
    CHECK(std::isinf(g.ni));
    CHECK(g.sir[0] == 1.0);  // the server slot is still a finite count
}

TEST_CASE("to_marginal sums the server block over phases") {
    qn::Network<double> m = open_fcfs();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    // A two-phase class: the server block is [k=1, k=2] and sir is their sum,
    // while kir keeps them apart. Getting this wrong collapses a PH station's
    // in-service count to its first phase, which is invisible in any total.
    const std::vector<double> row{1, 0, 1};  // one buffer tag, phases (0,1)
    const std::vector<std::size_t> phasesz{2}, phaseshift{0};
    const Marginal<double> g = to_marginal(sn, 2, row, phasesz, phaseshift);
    CHECK(g.kir[0][0] == 0.0);
    CHECK(g.kir[0][1] == 1.0);
    CHECK(g.sir[0] == 1.0);
    CHECK(g.nir[0] == 2.0);  // one in service plus one tagged in the buffer
}

TEST_CASE("to_marginal refuses a row narrower than its declared server block") {
    qn::Network<double> m = open_fcfs();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    // Three phases declared, two columns supplied. Refusing beats decoding a
    // short row into a plausible marginal.
    const std::vector<double> row{1, 1};
    const std::vector<std::size_t> phasesz{3}, phaseshift{0};
    CHECK_THROWS_AS(to_marginal(sn, 2, row, phasesz, phaseshift), line::InputError);
}

TEST_CASE("space_closed_single enumerates placements of n jobs over m phases") {
    // multichoose(2,3) = 4 rows: (0,3) (1,2) (2,1) (3,0). A class with no
    // service process yields NO rows, not one empty row -- otherwise it would
    // multiply a cartesian fold by one and survive silently.
    CHECK(qn::space_closed_single<double>(2, 3).size() == 4);
    CHECK(qn::space_closed_single<double>(1, 5).size() == 1);
    CHECK(qn::space_closed_single<double>(0, 0).empty());
    CHECK(qn::space_closed_single<double>(0, 3).empty());
}

TEST_CASE("cartesian treats an empty operand as the identity") {
    const std::vector<std::vector<double>> a{{1, 2}, {3, 4}};
    const std::vector<std::vector<double>> b{{9}, {8}};
    const std::vector<std::vector<double>> e;
    CHECK(qn::cartesian(e, a).size() == 2);
    CHECK(qn::cartesian(a, e).size() == 2);
    const std::vector<std::vector<double>> ab = qn::cartesian(a, b);
    CHECK(ab.size() == 4);
    CHECK(ab[0].size() == 3);
    CHECK(ab[0][0] == 1.0);
    CHECK(ab[0][2] == 9.0);  // a's row outer, b's row inner
}

TEST_CASE("from_marginal round-trips through to_marginal on a PS station") {
    qn::Network<double> m = open_ps();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    // PS has no buffer, so each state IS the server block and to_marginal must
    // read back exactly the marginal that generated it. A decoding error shows
    // up here as a mismatch rather than as a plausible wrong number later.
    const std::vector<std::size_t> n{3}, phases{2};
    const std::vector<std::vector<double>> sp = qn::from_marginal(sn, 2, n, phases);
    CHECK(sp.size() == 4);  // 3 jobs over 2 phases
    const std::vector<std::size_t> phaseshift{0};
    for (std::size_t i = 0; i < sp.size(); ++i) {
        const qn::Marginal<double> g = qn::to_marginal(sn, 2, sp[i], phases, phaseshift);
        CHECK(g.nir[0] == 3.0);
        CHECK(g.sir[0] == 3.0);
    }
}

TEST_CASE("from_marginal enumerates the ordered FCFS buffer") {
    qn::Network<double> m = open_fcfs();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    // Single class, so the multiset {1,1} has ONE distinct permutation: with
    // one server the state is buffer [1] and one job in service. Every row must
    // decode back to the marginal that generated it.
    const std::vector<std::size_t> n{2}, phases{1};
    const std::vector<std::vector<double>> sp = qn::from_marginal(sn, 2, n, phases);
    CHECK(sp.size() == 1);
    const std::vector<std::size_t> phaseshift{0};
    for (std::size_t i = 0; i < sp.size(); ++i) {
        const qn::Marginal<double> g = qn::to_marginal(sn, 2, sp[i], phases, phaseshift);
        CHECK(g.nir[0] == 2.0);
        CHECK(g.sir[0] == 1.0);   // one in service
    }
}

TEST_CASE("from_marginal orders a two-class FCFS buffer by permutation") {
    qn::Network<double> m = fcfs_two_class();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    // One class-1 and one class-2 job over a single server: the two orderings
    // (1 waits / 2 waits) are DISTINCT states, which is exactly what a
    // count-shaped buffer would collapse into one.
    const std::vector<std::size_t> n{1, 1}, phases{1, 1};
    const std::vector<std::vector<double>> sp = qn::from_marginal(sn, 2, n, phases);
    CHECK(sp.size() == 2);
    const std::vector<std::size_t> phaseshift{0, 1};
    for (std::size_t i = 0; i < sp.size(); ++i) {
        const qn::Marginal<double> g = qn::to_marginal(sn, 2, sp[i], phases, phaseshift);
        CHECK(g.nir[0] == 1.0);
        CHECK(g.nir[1] == 1.0);
        CHECK(g.ni == 2.0);
    }
}

TEST_CASE("from_marginal yields one idle state for an empty FCFS station") {
    qn::Network<double> m = open_fcfs();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<std::size_t> n{0}, phases{1};
    const std::vector<std::vector<double>> sp = qn::from_marginal(sn, 2, n, phases);
    CHECK(sp.size() == 1);
    for (std::size_t j = 0; j < sp[0].size(); ++j) CHECK(sp[0][j] == 0.0);
}

TEST_CASE("space_closed_multi places jobs over stateful nodes") {
    // One class, 2 jobs over 2 nodes: (0,2) (1,1) (2,0).
    CHECK(qn::space_closed_multi<double>(2, std::vector<std::size_t>{2}).size() == 3);
    // Two classes fold: 3 x 3 = 9 rows, 4 columns (2 nodes x 2 classes).
    const std::vector<std::vector<double>> ss =
        qn::space_closed_multi<double>(2, std::vector<std::size_t>{2, 2});
    CHECK(ss.size() == 9);
    CHECK(ss[0].size() == 4);
}

TEST_CASE("space_closed_multi_cs enumerates the split within a chain") {
    // Two classes in ONE chain holding 2 jobs: the chain total is invariant but
    // the class split is not, so (2,0) (1,1) (0,2) must all appear. Enumerating
    // per-class populations alone would fix a split the model does not fix and
    // drop every state a class switch can reach.
    const std::vector<std::vector<bool>> one_chain{{true, true}};
    const std::vector<std::vector<double>> cs =
        qn::space_closed_multi_cs<double>(2, std::vector<std::size_t>{1, 1}, one_chain);
    // 3 splits of the chain total, each placed over 2 nodes:
    // (2,0) -> 3x1, (1,1) -> 2x2, (0,2) -> 1x3, so 3+4+3 = 10 rows.
    CHECK(cs.size() == 10);

    // The same populations with the classes in SEPARATE chains fix the split,
    // so the space is strictly smaller.
    const std::vector<std::vector<bool>> two_chains{{true, false}, {false, true}};
    const std::vector<std::vector<double>> sep =
        qn::space_closed_multi_cs<double>(2, std::vector<std::size_t>{1, 1}, two_chains);
    CHECK(sep.size() == 4);
    CHECK(sep.size() < cs.size());
}

TEST_CASE("space_generator refuses an open class with no cutoff") {
    qn::Network<double> m = open_ps();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    // An open class has unbounded population: truncating silently would return
    // a state space that looks complete and is not, so a zero cutoff errors.
    CHECK_THROWS_AS(qn::space_generator(sn, std::vector<std::size_t>{0}), line::InputError);
}

TEST_CASE("space_generator builds a closed network's states") {
    qn::Network<double> m("cq");
    const std::size_t d = m.add_delay("D");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, line::lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q, c, line::lang::Distrib<double>::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    // 2 jobs over 2 stations, one exponential phase each: the population
    // lattice has 3 rows and each fixes both local states, so 3 network states.
    const std::vector<qn::NetState<double>> ss =
        qn::space_generator(sn, std::vector<std::size_t>{0});
    CHECK(ss.size() == 3);
    for (std::size_t i = 0; i < ss.size(); ++i) CHECK(ss[i].local.size() == 2);
}

TEST_CASE("from_marginal encodes a retrial station as an (in-service, orbit) split") {
    qn::Network<double> m("retr");
    const std::size_t s = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(s, c, line::lang::Distrib<double>::exp_rate(0.5));
    m.set_service(q, c, line::lang::Distrib<double>::exp_rate(1.0));
    m.set_number_of_servers(q, 1);
    m.set_retrial(q, c, line::lang::Distrib<double>::exp_rate(2.0), 2.0);
    qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t ist = sn.nodes[q - 1].station;

    // Two jobs, one server: the admissible splits are (in service, orbit) =
    // (0,2) and (1,1). The (0,2) row is the one an ordinary queue can NEVER
    // occupy -- a job waits while the server idles -- and it exists precisely
    // because a completion does not promote from the orbit.
    const std::vector<std::vector<double>> sp =
        qn::from_marginal(sn, ist, std::vector<std::size_t>{2}, std::vector<std::size_t>{1});
    CHECK(sp.size() == 2);
    // Layout is [orbit slots | server phases], the orbit padded to n columns.
    for (std::size_t i = 0; i < sp.size(); ++i) CHECK(sp[i].size() == 3);
    // Row for csrv = 0: orbit holds both jobs, server empty.
    CHECK(sp[0][0] == doctest::Approx(1.0));
    CHECK(sp[0][1] == doctest::Approx(1.0));
    CHECK(sp[0][2] == doctest::Approx(0.0));
    // Row for csrv = 1: one in orbit, one in service.
    CHECK(sp[1][0] == doctest::Approx(0.0));
    CHECK(sp[1][1] == doctest::Approx(1.0));
    CHECK(sp[1][2] == doctest::Approx(1.0));

    // The empty station keeps one all-zero row of server width, no orbit slot.
    const std::vector<std::vector<double>> e0 =
        qn::from_marginal(sn, ist, std::vector<std::size_t>{0}, std::vector<std::size_t>{1});
    CHECK(e0.size() == 1);
    CHECK(e0[0].size() == 1);
    CHECK(e0[0][0] == doctest::Approx(0.0));
}

TEST_CASE("from_marginal_node encodes a Transition per mode, not per class") {
    qn::Network<double> m("spn");
    const std::size_t p1 = m.add_place("P1");
    const std::size_t c = m.add_closed_class("C1", 1.0, p1);
    m.set_service(p1, c, line::lang::Distrib<double>::exp_rate(1.0));
    qn::TransitionParam<double> tp;
    tp.nmodes = 2;
    tp.modenames.push_back("m1");
    tp.modenames.push_back("m2");
    tp.enabling.assign(2, line::Matrix<double>(1, 1, 1.0));
    tp.inhibiting.assign(2, line::Matrix<double>(1, 1, 0.0));
    tp.firing.assign(2, line::Matrix<double>(1, 1, 1.0));
    tp.nmodeservers.push_back(1.0);
    tp.nmodeservers.push_back(std::numeric_limits<double>::infinity());
    tp.firingphases.push_back(1);
    tp.firingphases.push_back(3);
    const std::size_t tr = m.add_transition("T1", tp);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // A Transition is stateful but is NOT a station, so it has no station index
    // and only the node-indexed entry can reach it. The row is per MODE:
    // [free servers | firing phases | fired], width nmodes + sum(firingphases)
    // + nmodes = 2 + 4 + 2 = 8, and the marginal plays no part in it.
    const std::vector<std::vector<double>> sp =
        qn::from_marginal_node(sn, tr, std::vector<std::size_t>{0}, std::vector<std::size_t>{1});
    CHECK(sp.size() == 1);
    CHECK(sp[0].size() == 8);
    CHECK(sp[0][0] == doctest::Approx(1.0));
    // An infinite server count is clamped: the row is a COUNT vector, and Inf
    // is not a count -- the same distinction that made a CTMC report Q = Inf.
    CHECK(sp[0][1] == doctest::Approx(2147483647.0));
    for (std::size_t j = 2; j < sp[0].size(); ++j) CHECK(sp[0][j] == doctest::Approx(0.0));

    // The marginal is ignored: a different one yields the same all-idle row.
    const std::vector<std::vector<double>> sp2 =
        qn::from_marginal_node(sn, tr, std::vector<std::size_t>{5}, std::vector<std::size_t>{1});
    CHECK(sp2 == sp);
}

TEST_CASE("space_generator gives a Transition its own local block") {
    qn::Network<double> m("spn2");
    const std::size_t p1 = m.add_place("P1");
    const std::size_t c = m.add_closed_class("C1", 1.0, p1);
    m.set_service(p1, c, line::lang::Distrib<double>::exp_rate(1.0));
    qn::TransitionParam<double> tp;
    tp.nmodes = 1;
    tp.modenames.push_back("m1");
    tp.enabling.assign(1, line::Matrix<double>(1, 1, 1.0));
    tp.inhibiting.assign(1, line::Matrix<double>(1, 1, 0.0));
    tp.firing.assign(1, line::Matrix<double>(1, 1, 1.0));
    tp.nmodeservers.push_back(1.0);
    tp.firingphases.push_back(1);
    m.add_transition("T1", tp);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // The lattice is over the STATEFUL NODES: iterating stations would drop the
    // Transition's block, leaving a state vector that silently omits the firing
    // state. Place and Transition are both stateful, so every state has two.
    const std::vector<qn::NetState<double>> ss =
        qn::space_generator(sn, std::vector<std::size_t>{0});
    CHECK(!ss.empty());
    for (std::size_t i = 0; i < ss.size(); ++i) {
        CHECK(ss[i].local.size() == 2);
        CHECK(ss[i].local[1].size() == 3);  // 1 server + 1 phase + 1 fired
    }
}
