/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * PNML (ISO/IEC 15909-2) place/transition import and export.
 *
 * The oracle is not a recorded file. A round trip through pnml_save/pnml_load
 * must leave every CTMC metric of the model unchanged AND must reproduce the
 * file byte for byte, which a consistent error on both sides of the interchange
 * would not survive; and the untimed net read from another tool is checked
 * against a marginal computed by hand from the chain it defines.
 *
 * Twin of jar/src/test/java/jline/io/PnmlIOTest.java and of
 * python/tests/test_pnml_io.py.
 */
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/io/pnml.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/util/matrix.h"

using line::Matrix;

namespace {

typedef line::qn::Network<double> Net;
typedef line::qn::RoutingMatrix<double> Routing;
typedef line::lang::Distrib<double> D;
typedef line::qn::TransitionParam<double> TP;

const double kInf = std::numeric_limits<double>::infinity();

/** A temporary path inside the build's own directory. */
std::string tmp_path(const std::string& stem) {
    return "/tmp/line_pnml_test_" + stem;
}

/** One mode of a transition, with its pre, inhibitor and post arcs. */
struct ModeSpec {
    std::string name;
    D proc = D::disabled_dist();
    line::lang::TimingStrategy timing = line::lang::TimingStrategy::TIMED;
    double servers = 1.0, prio = 1.0, weight = 1.0;
    std::vector<std::pair<std::size_t, double> > enab, inhib, fire;
};

TP make_param(std::size_t nnodes, const std::vector<ModeSpec>& modes) {
    TP tp;
    tp.nmodes = modes.size();
    for (std::size_t m = 0; m < modes.size(); ++m) {
        const ModeSpec& md = modes[m];
        const bool immediate = md.timing == line::lang::TimingStrategy::IMMEDIATE;
        tp.modenames.push_back(md.name);
        tp.timing.push_back(md.timing);
        tp.firingproc.push_back(md.proc);
        tp.firingphases.push_back(immediate || md.proc.disabled
                                      ? 0
                                      : line::lang::dist_to_map(md.proc).order());
        tp.nmodeservers.push_back(md.servers);
        tp.firingprio.push_back(md.prio);
        tp.fireweight.push_back(md.weight);
        Matrix<double> en(nnodes, 1, 0.0), ih(nnodes, 1, kInf), fi(nnodes, 1, 0.0);
        for (std::size_t a = 0; a < md.enab.size(); ++a) en(md.enab[a].first - 1, 0) = md.enab[a].second;
        for (std::size_t a = 0; a < md.inhib.size(); ++a) ih(md.inhib[a].first - 1, 0) = md.inhib[a].second;
        for (std::size_t a = 0; a < md.fire.size(); ++a) fi(md.fire[a].first - 1, 0) = md.fire[a].second;
        tp.enabling.push_back(en);
        tp.inhibiting.push_back(ih);
        tp.firing.push_back(fi);
    }
    return tp;
}

/** Two places, one timed mode each, arc weight 2 on one side. */
Net two_modes() {
    Net m("twomodes");
    const std::size_t p1 = m.add_place("P1");
    const std::size_t p2 = m.add_place("P2");
    const std::size_t c1 = m.add_closed_class("Class1", 4.0, p1);

    ModeSpec a;
    a.name = "Mode1";
    a.proc = D::exp_rate(2.0);
    a.enab.push_back(std::make_pair(p1, 2.0));
    a.fire.push_back(std::make_pair(p2, 2.0));
    const std::size_t t1 = m.add_transition("T1", make_param(4, std::vector<ModeSpec>(1, a)));

    ModeSpec b;
    b.name = "Mode2";
    b.proc = D::erlang(1.5, 2);
    b.enab.push_back(std::make_pair(p2, 1.0));
    b.fire.push_back(std::make_pair(p1, 1.0));
    const std::size_t t2 = m.add_transition("T2", make_param(4, std::vector<ModeSpec>(1, b)));

    Routing R;
    R.set(c1, c1, p1, t1, 1.0);
    R.set(c1, c1, p2, t2, 1.0);
    R.set(c1, c1, t1, p2, 1.0);
    R.set(c1, c1, t2, p1, 1.0);
    m.link(R);
    m.set_initial_marking(p1, std::vector<double>(1, 4.0));
    m.set_initial_marking(p2, std::vector<double>(1, 0.0));
    return m;
}

/** One transition carrying TWO modes, plus an immediate mode with an inhibitor arc. */
Net multi_mode() {
    Net m("multimode");
    const std::size_t q1 = m.add_place("Q1");
    const std::size_t q2 = m.add_place("Q2");
    const std::size_t c1 = m.add_closed_class("Class1", 3.0, q1);

    std::vector<ModeSpec> u1modes;
    ModeSpec fast;
    fast.name = "Fast";
    fast.proc = D::exp_rate(3.0);
    fast.enab.push_back(std::make_pair(q1, 1.0));
    fast.fire.push_back(std::make_pair(q2, 1.0));
    u1modes.push_back(fast);
    ModeSpec slow;
    slow.name = "Slow";
    slow.proc = D::exp_rate(1.0);
    slow.enab.push_back(std::make_pair(q1, 2.0));
    slow.fire.push_back(std::make_pair(q2, 2.0));
    u1modes.push_back(slow);
    const std::size_t u1 = m.add_transition("U1", make_param(4, u1modes));

    ModeSpec back;
    back.name = "Back";
    back.timing = line::lang::TimingStrategy::IMMEDIATE;
    back.proc = D::immediate();
    back.weight = 2.5;
    back.enab.push_back(std::make_pair(q2, 1.0));
    back.inhib.push_back(std::make_pair(q1, 3.0));
    back.fire.push_back(std::make_pair(q1, 1.0));
    const std::size_t u2 = m.add_transition("U2", make_param(4, std::vector<ModeSpec>(1, back)));

    Routing R;
    R.set(c1, c1, q1, u1, 1.0);
    R.set(c1, c1, u1, q2, 1.0);
    R.set(c1, c1, q2, u2, 1.0);
    R.set(c1, c1, q1, u2, 1.0);
    R.set(c1, c1, u2, q1, 1.0);
    m.link(R);
    m.set_initial_marking(q1, std::vector<double>(1, 3.0));
    m.set_initial_marking(q2, std::vector<double>(1, 0.0));
    return m;
}

line::mva::AvgResult<double> solve(Net& m, double cutoff) {
    line::ctmc::CtmcOptions opt;
    opt.cutoff = cutoff;
    return line::ctmc::solver_ctmc_run_analyzer(m.get_struct(), opt);
}

std::string read_file(const std::string& path) {
    std::ifstream in(path.c_str());
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

void write_file(const std::string& path, const std::string& text) {
    std::ofstream out(path.c_str());
    out << text;
}

}  // namespace

TEST_CASE("pnml round trip leaves every metric unchanged") {
    Net m = two_modes();
    const std::string path = tmp_path("twomodes.pnml");
    line::io::pnml_save(m.get_struct(), path);
    Net back = line::io::pnml_load<double>(path);

    const line::mva::AvgResult<double> a = solve(m, 4);
    const line::mva::AvgResult<double> b = solve(back, 4);
    REQUIRE(a.QN.rows() == b.QN.rows());
    for (std::size_t i = 0; i < a.QN.rows(); ++i) {
        CHECK(a.QN(i, 0) == doctest::Approx(b.QN(i, 0)).epsilon(1e-12));
        CHECK(a.TN(i, 0) == doctest::Approx(b.TN(i, 0)).epsilon(1e-12));
        CHECK(a.UN(i, 0) == doctest::Approx(b.UN(i, 0)).epsilon(1e-12));
    }
}

TEST_CASE("a write, read and second write of a pnml net are byte identical") {
    // A field lost in the round trip changes the second file. Comparing metrics
    // sees only what the solver reads; comparing the FILE sees every attribute
    // the writer emits, and needs no solver at all.
    const char* stems[2] = {"twomodes", "multimode"};
    for (int k = 0; k < 2; ++k) {
        Net m = k == 0 ? two_modes() : multi_mode();
        const std::string first = tmp_path(std::string(stems[k]) + "-1.pnml");
        const std::string second = tmp_path(std::string(stems[k]) + "-2.pnml");
        line::io::pnml_save(m.get_struct(), first);
        Net back = line::io::pnml_load<double>(first);
        line::io::pnml_save(back.get_struct(), second);
        CHECK(read_file(first) == read_file(second));
    }
}

TEST_CASE("pnml splits a transition's modes and regroups them on the way back") {
    Net m = multi_mode();
    const std::string path = tmp_path("multimode.pnml");
    line::io::pnml_save(m.get_struct(), path);
    const std::string text = read_file(path);
    CHECK(text.find("id=\"U1.Fast\"") != std::string::npos);
    CHECK(text.find("id=\"U1.Slow\"") != std::string::npos);
    CHECK(text.find("value=\"inhibitor\"") != std::string::npos);

    Net back = line::io::pnml_load<double>(path);
    const line::qn::NetworkStruct<double>& sn = back.get_struct();
    std::size_t ntrans = 0, u1 = 0, u2 = 0, q1 = 0;
    for (std::size_t i = 0; i < sn.nodes.size(); ++i) {
        if (sn.nodes[i].nodetype == line::lang::NodeType::Transition) {
            ++ntrans;
            if (sn.nodes[i].name == "U1") u1 = i + 1;
            if (sn.nodes[i].name == "U2") u2 = i + 1;
        } else if (sn.nodes[i].name == "Q1") {
            q1 = i + 1;
        }
    }
    CHECK(ntrans == 2);
    REQUIRE(u1 != 0);
    REQUIRE(u2 != 0);
    CHECK(sn.transparam.at(u1).nmodes == 2);
    CHECK(sn.transparam.at(u2).timing[0] == line::lang::TimingStrategy::IMMEDIATE);
    CHECK(sn.transparam.at(u2).fireweight[0] == doctest::Approx(2.5));
    CHECK(sn.transparam.at(u2).inhibiting[0](q1 - 1, 0) == doctest::Approx(3.0));
}

TEST_CASE("an untimed pnml net from another tool reads and solves") {
    // No toolspecific block and no timing: every transition is read as Exp(1)
    // with one server. Two indistinguishable tokens then cycle between two
    // places, so the number in `busy` is a 3-state birth-death chain with equal
    // rates; the stationary distribution is uniform, each place holds a mean of
    // one token, and the cycle fires at 2/3.
    const std::string path = tmp_path("foreign.pnml");
    write_file(path,
               "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
               "<pnml xmlns=\"http://www.pnml.org/version-2009/grammar/pnml\">\n"
               "  <net id=\"foreign\" type=\"http://www.pnml.org/version-2009/grammar/ptnet\">\n"
               "    <page id=\"p0\">\n"
               "      <place id=\"ready\"><initialMarking><text>2</text></initialMarking></place>\n"
               "      <place id=\"busy\"><initialMarking><text>0</text></initialMarking></place>\n"
               "      <transition id=\"start\"/>\n"
               "      <transition id=\"finish\"/>\n"
               "      <arc id=\"e1\" source=\"ready\" target=\"start\"/>\n"
               "      <arc id=\"e2\" source=\"start\" target=\"busy\"/>\n"
               "      <arc id=\"e3\" source=\"busy\" target=\"finish\"/>\n"
               "      <arc id=\"e4\" source=\"finish\" target=\"ready\"/>\n"
               "    </page>\n  </net>\n</pnml>\n");

    Net m = line::io::pnml_load<double>(path);
    const line::mva::AvgResult<double> r = solve(m, 2);
    REQUIRE(r.QN.rows() == 2);
    CHECK(r.QN(0, 0) == doctest::Approx(1.0).epsilon(1e-7));
    CHECK(r.QN(1, 0) == doctest::Approx(1.0).epsilon(1e-7));
    CHECK(r.TN(1, 0) == doctest::Approx(2.0 / 3.0).epsilon(1e-7));
}

TEST_CASE("pnml refuses what the place/transition grammar cannot carry") {
    // An open class has no P/T counterpart: the grammar has no unbounded token
    // source, and writing the net without one would be a different model.
    Net m("open");
    const std::size_t src = m.add_source("Source");
    const std::size_t p1 = m.add_place("P1");
    const std::size_t sink = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("Class1", 0);
    m.set_arrival(src, c1, D::exp_rate(1.0));

    ModeSpec md;
    md.name = "Mode1";
    md.proc = D::exp_rate(2.0);
    md.enab.push_back(std::make_pair(p1, 1.0));
    md.fire.push_back(std::make_pair(sink, 1.0));
    const std::size_t t1 = m.add_transition("T1", make_param(4, std::vector<ModeSpec>(1, md)));

    Routing R;
    R.set(c1, c1, src, p1, 1.0);
    R.set(c1, c1, p1, t1, 1.0);
    R.set(c1, c1, t1, sink, 1.0);
    m.link(R);

    CHECK_THROWS_AS(line::io::pnml_save(m.get_struct(), tmp_path("bad.pnml")), line::InputError);
}
