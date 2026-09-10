/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `spn_lpbnd` and the SolverBA `spnlp.*` family: moment-relaxation LP bounds
 * for a stochastic Petri net.
 *
 * THE ACCEPTANCE TEST IS THE REFERENCE'S OWN TABLE 2. Liu (1998) publishes four
 * bound columns for the four-server production line of its Fig. 2b, on five
 * rate vectors, and this file asserts all four to three decimals. Those numbers
 * are what a transcription error in any constraint family would move: the
 * polytope has twelve of them and the published optimum is a function of the
 * whole set, so a bracket check alone would not catch a family that was dropped
 * or mis-signed.
 *
 * The bracket checks are the other half. A bound that agreed with the paper and
 * still failed to contain the exact answer would be wrong whatever it agreed
 * with. The exact values used here are the ones `test_spn_pf.cpp` already pins
 * for the fork-join net, which shares no code path with this LP.
 *
 * SINGLE CLASS ONLY. `spn_lpbnd` puts one level per place because
 * `NetworkStruct::transparam` carries no class dimension, the same restriction
 * `spn_mdd` states, and refuses a coloured net by name. The MATLAB, JAR and
 * python twins carry the class axis; the last case pins that refusal.
 *
 * Reference: Z. Liu, "Performance Analysis of Stochastic Timed Petri Nets Using
 * Linear Programming Approach", IEEE Trans. Software Engineering 24(11), 1998,
 * 1014-1030.
 */

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/spn/spn_lpbnd.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ba/solver_ba_runner.h"

using namespace line;
namespace qn = line::qn;
namespace spn = line::spn;
using line::lang::TimingStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** MATLAB, native python and the JAR on the fork-join net, at %.12f. */
const double FJ_TOKENS[4] = {0.710262429604, 1.864198441415, 1.864198441415, 0.425539128981};
const double FJ_MODEX = 0.607360882508;

/** Liu Table 2, p. 1023: mu[4], then l.b., u.b.2, o.l.b., o.u.b. */
struct Table2Row {
    double mu[4];
    double lb, ub2, olb, oub;
};
const Table2Row TABLE2[5] = {
    {{1.000, 1.25, 2.00, 0.50}, 1.165, 2.000, 0.930, 2.000},
    {{1.000, 1.25, 2.00, 2.50}, 1.829, 3.529, 1.481, 4.000},
    {{1.000, 1.25, 1.25, 2.50}, 1.581, 3.333, 1.333, 4.000},
    {{1.000, 1.25, 1.25, 1.00}, 1.359, 3.333, 1.111, 4.000},
    {{1.111, 1.111, 1.111, 1.111}, 1.350, 2.963, 1.111, 4.444}};

double vsum(const std::vector<double>& v) {
    double s = 0;
    for (std::size_t i = 0; i < v.size(); ++i) s += v[i];
    return s;
}

/** One timed mode consuming `from` and producing into every place in `to`. */
qn::TransitionParam<double> mode_of(std::size_t nnodes, const std::vector<std::size_t>& from,
                                    const std::vector<std::size_t>& to, double rate,
                                    double servers) {
    qn::TransitionParam<double> tp;
    tp.nmodes = 1;
    tp.modenames.push_back("fire");
    tp.enabling.assign(1, line::Matrix<double>(nnodes, 1, 0.0));
    tp.inhibiting.assign(1,
                         line::Matrix<double>(nnodes, 1, std::numeric_limits<double>::infinity()));
    tp.firing.assign(1, line::Matrix<double>(nnodes, 1, 0.0));
    for (std::size_t i = 0; i < from.size(); ++i) tp.enabling[0](from[i] - 1, 0) = 1.0;
    for (std::size_t i = 0; i < to.size(); ++i) tp.firing[0](to[i] - 1, 0) = 1.0;
    tp.nmodeservers.push_back(servers);
    tp.firingphases.push_back(1);
    tp.timing.push_back(TimingStrategy::TIMED);
    tp.fireweight.push_back(1.0);
    tp.firingproc.push_back(Dist::exp_rate(rate));
    return tp;
}

std::vector<std::size_t> one(std::size_t a) { return std::vector<std::size_t>(1, a); }

std::vector<std::size_t> two(std::size_t a, std::size_t b) {
    std::vector<std::size_t> v;
    v.push_back(a);
    v.push_back(b);
    return v;
}

/**
 * Fig. 2b: four servers, blocking before service, buffers of 3, 2 and 4.
 *
 * (p5, p2), (p4, p1) and (p3, p0) are the three buffer pairs, each conserved at
 * its capacity, so the net is a strongly connected marked graph and every
 * transition carries the same throughput. The struct carries no per-place
 * state, so the split marking is passed through `SpnLpOptions::init`.
 */
qn::Network<double> prodline(const double mu[4]) {
    qn::Network<double> m("liu98");
    const std::size_t p5 = m.add_place("p5");
    const std::size_t p4 = m.add_place("p4");
    const std::size_t p3 = m.add_place("p3");
    const std::size_t p2 = m.add_place("p2");
    const std::size_t p1 = m.add_place("p1");
    const std::size_t p0 = m.add_place("p0");
    const std::size_t c = m.add_closed_class("Class1", 9.0, p2);
    const std::size_t pl[6] = {p5, p4, p3, p2, p1, p0};
    for (int i = 0; i < 6; ++i) m.set_service(pl[i], c, Dist::exp_rate(1.0));
    const std::size_t nn = 10;
    m.add_transition("t1", mode_of(nn, one(p2), one(p5), mu[0], 1.0));
    m.add_transition("t2", mode_of(nn, two(p5, p1), two(p4, p2), mu[1], 1.0));
    m.add_transition("t3", mode_of(nn, two(p4, p0), two(p3, p1), mu[2], 1.0));
    m.add_transition("t4", mode_of(nn, one(p3), one(p0), mu[3], 1.0));
    return m;
}

/** The initial marking of `prodline`, in place-major level order. */
std::vector<double> prodline_init() {
    std::vector<double> v;
    v.push_back(0);  // p5
    v.push_back(0);  // p4
    v.push_back(0);  // p3
    v.push_back(3);  // p2
    v.push_back(2);  // p1
    v.push_back(4);  // p0
    return v;
}

/**
 * P0 -(Tf)-> P1 + P2 -(Tj)-> P3 -(Tb)-> P0.
 *
 * Tf consumes one token and produces two, so the marking is not a conserved
 * population; the two minimal-support invariants are (1,1,0,1) and (1,0,1,1) at
 * value 3. This is what exercises the WEIGHTED form of the invariant family --
 * the reference writes that family for unweighted cycles.
 */
qn::Network<double> forkjoin_spn(int ntokens) {
    qn::Network<double> m("fj");
    std::vector<std::size_t> pl;
    for (std::size_t i = 0; i < 4; ++i) pl.push_back(m.add_place("P" + std::to_string(i)));
    const std::size_t c = m.add_closed_class("C", static_cast<double>(ntokens), pl[0]);
    for (std::size_t i = 0; i < 4; ++i) m.set_service(pl[i], c, Dist::exp_rate(1.0));
    m.add_transition("Tf", mode_of(7, one(pl[0]), two(pl[1], pl[2]), 1.3, 1.0));
    m.add_transition("Tj", mode_of(7, two(pl[1], pl[2]), one(pl[3]), 0.7, 1.0));
    m.add_transition("Tb", mode_of(7, one(pl[3]), one(pl[0]), 1.9, 1.0));
    return m;
}

}  // namespace

TEST_CASE("spn_lpbnd: the reference's Table 2 upper column, to three decimals") {
    for (int c = 0; c < 5; ++c) {
        qn::Network<double> m = prodline(TABLE2[c].mu);
        spn::SpnLpOptions o;
        o.markovian = true;
        o.init = prodline_init();
        const spn::SpnLpBounds b = spn::spn_lpbnd(m.get_struct(), o);
        INFO("Table 2 case ", c + 1, " u.b.2");
        CHECK(std::abs(vsum(b.mode_tput_hi) - TABLE2[c].ub2) < 5e-4);
    }
}

TEST_CASE("spn_lpbnd: the published lower column needs the liveness rows") {
    // They are opt-in because they hold only on a live net, and this one is: a
    // strongly connected marked graph with a token on every cycle. Without them
    // the lower side falls back to the operational value, which is the second
    // check here and is why the default is not a silent loss.
    for (int c = 0; c < 5; ++c) {
        qn::Network<double> m = prodline(TABLE2[c].mu);
        spn::SpnLpOptions live;
        live.markovian = true;
        live.assumelive = true;
        live.init = prodline_init();
        const spn::SpnLpBounds bl = spn::spn_lpbnd(m.get_struct(), live);
        INFO("Table 2 case ", c + 1, " l.b.");
        CHECK(std::abs(vsum(bl.mode_tput_lo) - TABLE2[c].lb) < 5e-4);
        spn::SpnLpOptions plain;
        plain.init = prodline_init();
        const spn::SpnLpBounds bp = spn::spn_lpbnd(m.get_struct(), plain);
        CHECK(vsum(bp.mode_tput_lo) <= vsum(bl.mode_tput_lo) + 1e-9);
    }
}

TEST_CASE("spn_lpbnd: both operational columns, with no Markovian family at all") {
    for (int c = 0; c < 5; ++c) {
        qn::Network<double> m = prodline(TABLE2[c].mu);
        spn::SpnLpOptions o;
        o.markovian = false;
        o.assumelive = true;
        o.init = prodline_init();
        const spn::SpnLpBounds b = spn::spn_lpbnd(m.get_struct(), o);
        INFO("Table 2 case ", c + 1, " operational");
        CHECK(std::abs(vsum(b.mode_tput_lo) - TABLE2[c].olb) < 5e-4);
        CHECK(std::abs(vsum(b.mode_tput_hi) - TABLE2[c].oub) < 5e-4);
    }
}

TEST_CASE("spn_lpbnd: the fork-join bracket contains the exact product form") {
    for (int pass = 0; pass < 2; ++pass) {
        qn::Network<double> m = forkjoin_spn(3);
        spn::SpnLpOptions o;
        o.markovian = pass == 0;
        const spn::SpnLpBounds b = spn::spn_lpbnd(m.get_struct(), o);
        REQUIRE(b.nplacelevels == 4);
        for (std::size_t l = 0; l < 4; ++l) {
            INFO("level ", l, " markovian=", o.markovian);
            CHECK(b.tokens_lo[l] <= FJ_TOKENS[l] + 1e-6);
            CHECK(FJ_TOKENS[l] <= b.tokens_hi[l] + 1e-6);
            CHECK(b.place_tput_lo[l] <= FJ_MODEX + 1e-6);
            CHECK(FJ_MODEX <= b.place_tput_hi[l] + 1e-6);
        }
    }
}

TEST_CASE("spn_lpbnd: the level caps come off the weighted invariants") {
    // An unweighted reading would be wrong here: the two minimal-support
    // invariants are (1,1,0,1) and (1,0,1,1) at value 3, so every level is
    // capped at 3 rather than at the token count.
    qn::Network<double> m = forkjoin_spn(3);
    const spn::SpnLpBounds b = spn::spn_lpbnd(m.get_struct(), spn::SpnLpOptions());
    for (std::size_t l = 0; l < b.nplacelevels; ++l) CHECK(b.bound[l] == doctest::Approx(3.0));
}

TEST_CASE("spn_lpbnd: the bracket is never the whole box") {
    // A [0, B] bracket is a missing constraint family, not a loose bound: it
    // means the polytope does not constrain the objective at all.
    qn::Network<double> m = forkjoin_spn(3);
    const spn::SpnLpBounds b = spn::spn_lpbnd(m.get_struct(), spn::SpnLpOptions());
    for (std::size_t l = 0; l < b.nplacelevels; ++l) {
        INFO("level ", l);
        CHECK(b.tokens_hi[l] - b.tokens_lo[l] < b.bound[l] - 1e-9);
    }
}

TEST_CASE("spn_lpbnd: the Markovian polytope is the operational one plus rows") {
    qn::Network<double> m = forkjoin_spn(3);
    spn::SpnLpOptions op;
    op.markovian = false;
    const spn::SpnLpBounds bm = spn::spn_lpbnd(m.get_struct(), spn::SpnLpOptions());
    const spn::SpnLpBounds bo = spn::spn_lpbnd(m.get_struct(), op);
    for (std::size_t l = 0; l < bm.nplacelevels; ++l) {
        CHECK(bo.tokens_lo[l] <= bm.tokens_lo[l] + 1e-9);
        CHECK(bm.tokens_hi[l] <= bo.tokens_hi[l] + 1e-9);
    }
}

TEST_CASE("ba: spnlp is the only family offered on a Petri net, and only there") {
    qn::Network<double> fj = forkjoin_spn(3);
    const std::vector<std::string> onNet = ba::list_valid_methods(fj.get_struct());
    CHECK(onNet.size() == 4);
    for (std::size_t i = 0; i < onNet.size(); ++i) CHECK(ba::is_spnlp_method(onNet[i]));

    qn::Network<double> cqn("cqn");
    const std::size_t d = cqn.add_delay("Delay");
    const std::size_t q = cqn.add_queue("Q", lang::SchedStrategy::PS);
    const std::size_t c = cqn.add_closed_class("C", 4.0, d);
    cqn.set_service(d, c, Dist::exp_rate(1.0));
    cqn.set_service(q, c, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    cqn.link(P);
    const std::vector<std::string> offNet = ba::list_valid_methods(cqn.get_struct());
    for (std::size_t i = 0; i < offNet.size(); ++i) CHECK_FALSE(ba::is_spnlp_method(offNet[i]));
}

TEST_CASE("ba: the four spnlp tokens reach the analyzer and land on the place rows") {
    const char* methods[4] = {"spnlp.upper", "spnlp.lower", "spnlp.op.upper", "spnlp.op.lower"};
    for (int i = 0; i < 4; ++i) {
        qn::Network<double> m = forkjoin_spn(3);
        ba::BaOptions opt;
        opt.method = methods[i];
        const mva::AvgResult<double> r = ba::solver_ba_run_analyzer(m.get_struct(), opt);
        const bool upper = std::string(methods[i]).find("upper") != std::string::npos;
        for (std::size_t p = 0; p < 4; ++p) {
            INFO(methods[i], " place ", p);
            // U = Q at an INF station, which is what a Place is
            CHECK(r.UN(p, 0) == doctest::Approx(r.QN(p, 0)));
            if (upper) {
                CHECK(r.QN(p, 0) >= FJ_TOKENS[p] - 1e-6);
            } else {
                CHECK(r.QN(p, 0) <= FJ_TOKENS[p] + 1e-6);
            }
        }
    }
}

TEST_CASE("spn_lpbnd: an infinite-server mode is refused by name") {
    qn::Network<double> m("is");
    const std::size_t p = m.add_place("P");
    const std::size_t c = m.add_closed_class("C", 2.0, p);
    m.set_service(p, c, Dist::exp_rate(1.0));
    m.add_transition("T", mode_of(3, one(p), one(p), 1.0,
                                  std::numeric_limits<double>::infinity()));
    CHECK_THROWS_AS(spn::spn_lpbnd(m.get_struct(), spn::SpnLpOptions()), UnsupportedError);
}
