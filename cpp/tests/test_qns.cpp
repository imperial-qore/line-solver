/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * SolverQNS: the JMVA marshalling, the qnsolver invocation and what comes back.
 *
 * THE ORACLE IS A LAW, not a recorded number. A golden would only pin whatever
 * this port produced on the day, and the point of wrapping an external tool is
 * that it is NOT this port. What is asserted instead is what must hold whichever
 * algorithm qnsolver ran: the closed population is conserved per class, Little's
 * law holds at every station, and the utilization of a single-server queue is
 * its throughput times its service time. A marshalling error -- a demand read as
 * a visit, a chain column read as a class one -- breaks all three at once.
 *
 * WHAT QNSOLVER ACTUALLY RUNS IS THE LINEARIZER, so the agreement with exact MVA
 * is only to that approximation's accuracy. The reference builds the command as
 * `qnsolver -l <model> -o <result>`, and `-l` is qnsolver's `--linearizer`: its
 * input file is POSITIONAL and it has no load flag. Measured on Delay(Z=1) ->
 * Queue(S=1), N=2, whose exact answer is X = 0.8: bare qnsolver and `-e` give
 * 0.8, `-l` gives 0.79872 and `-s` (Bard-Schweitzer) 0.763948. The port passes
 * the same `-l` as the reference does, deliberately -- the four codebases must
 * agree with each other before they agree with the exact answer -- so the
 * comparison below is at Linearizer accuracy and the tighter comparison is
 * against this port's OWN Linearizer.
 *
 * THE BINARY IS OPTIONAL. LINE ships no copy of qnsolver, so every case that
 * runs it is skipped when it is absent -- and the cases that do NOT need it (the
 * writer, the refusals, the method resolution) still run, so a machine without
 * the tool still tests everything the port owns.
 */

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/io/jmva_writer.h"
#include "line/io/qn2lqn.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/wrappers/qns/solver_qns.h"
#include "line/util/tempdir.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/** Delay -> Queue -> Delay, `nclasses` closed classes, one chain each. */
qn::Network<double> closed_cqn(std::size_t nclasses, double njobs, double servers = 1.0,
                               SchedStrategy sched = SchedStrategy::PS) {
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("Delay1");
    const std::size_t q = m.add_queue("Queue1", sched);
    if (servers != 1.0) m.set_number_of_servers(q, servers);
    qn::RoutingMatrix<double> P;
    for (std::size_t r = 0; r < nclasses; ++r) {
        const std::size_t c =
            m.add_closed_class("C" + std::to_string(r + 1), njobs, d);
        m.set_service(d, c, D::exp_rate(1.0 + 0.5 * static_cast<double>(r)));
        m.set_service(q, c, D::exp_rate(2.0 + static_cast<double>(r)));
        P.set(c, c, d, q, 1.0);
        P.set(c, c, q, d, 1.0);
    }
    m.link(P);
    return m;
}

/** Source -> Queue -> Queue -> Sink, one open class. */
qn::Network<double> open_oqn(double lambda) {
    qn::Network<double> m("oqn");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::PS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Open1");
    m.set_arrival(src, c, D::exp_rate(lambda));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, snk, 1.0);
    m.link(P);
    return m;
}

mva::AvgResult<double> run_mva(const qn::NetworkStruct<double>& sn, const std::string& method) {
    mva::MvaOptions opt;
    opt.method = method;
    Matrix<double> init;
    return mva::solver_mva_run_analyzer(sn, opt, init);
}

/**
 * The laws every product-form answer obeys, whichever algorithm produced it.
 *
 * The population invariant is the one that catches a chain column read as a
 * class column; Little's law catches a response time paired with the wrong
 * throughput; U = X*S catches a demand marshalled as a visit.
 */
void obeys_the_laws(const qn::NetworkStruct<double>& sn, const mva::AvgResult<double>& r) {
    for (std::size_t k = 0; k < sn.nclasses; ++k) {
        if (!std::isfinite(sn.classes[k].population)) continue;
        double total = 0.0;
        for (std::size_t i = 0; i < sn.nstations; ++i) total += r.QN(i, k);
        INFO("class " << k << " population");
        CHECK(std::fabs(total - sn.classes[k].population) <= 1e-4 * sn.classes[k].population);
    }
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        if (sn.stations[i].nodetype == qn::NodeType::Source) continue;
        for (std::size_t k = 0; k < sn.nclasses; ++k) {
            if (r.TN(i, k) <= 0.0) continue;
            INFO("Little at station " << i << " class " << k);
            CHECK(std::fabs(r.QN(i, k) - r.TN(i, k) * r.RN(i, k)) <=
                  1e-6 * std::max(1.0, r.QN(i, k)));
            const double c = sn.stations[i].nservers;
            if (!std::isfinite(c)) continue;
            // U = X*S/c at a queue: the service time is the reciprocal rate.
            const double S = 1.0 / sn.rates(i, k);
            INFO("U = X*S/c at station " << i << " class " << k);
            CHECK(std::fabs(r.UN(i, k) - r.TN(i, k) * S / c) <=
                  1e-5 * std::max(1.0, r.UN(i, k)));
            CHECK(r.UN(i, k) <= 1.0 + 1e-9);
        }
    }
}

/** Every finite entry of two tables agrees to `tol` in relative terms. */
void tables_agree(const Matrix<double>& a, const Matrix<double>& b, double tol,
                  const char* what) {
    REQUIRE(a.rows() == b.rows());
    REQUIRE(a.cols() == b.cols());
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < a.cols(); ++j) {
            const double x = a(i, j), y = b(i, j);
            const double scale = std::max(1.0, std::max(std::fabs(x), std::fabs(y)));
            INFO(what << " at (" << i << "," << j << "): " << x << " vs " << y);
            CHECK(std::fabs(x - y) <= tol * scale);
        }
}

std::string read_all(const std::string& path) {
    std::ifstream f(path.c_str());
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

const lqn::detail::RawActivity<double>& raw_activity(const lqn::LqnModel<double>& m,
                                                     const std::string& name) {
    for (std::size_t i = 0; i < m.acts.size(); ++i)
        if (m.acts[i].name == name) return m.acts[i];
    throw std::runtime_error("missing activity " + name);
}

const lqn::detail::RawPrecedence<double>& raw_precedence(const lqn::LqnModel<double>& m,
                                                         const std::string& pre,
                                                         lang::PrecedenceType posttype) {
    for (std::size_t t = 0; t < m.tasks.size(); ++t)
        for (std::size_t p = 0; p < m.tasks[t].precedences.size(); ++p) {
            const lqn::detail::RawPrecedence<double>& pr = m.tasks[t].precedences[p];
            if (pr.posttype == posttype && pr.preacts.size() == 1 && pr.preacts[0] == pre)
                return pr;
        }
    throw std::runtime_error("missing precedence from " + pre);
}

}  // namespace

// ---------------------------------------------------------------------------
// The JMVA document, which is what the port actually owns
// ---------------------------------------------------------------------------

TEST_CASE("writeJMVA emits the chain-level document of a closed model") {
    qn::Network<double> m = closed_cqn(1, 4.0);
    util::TempDir tmp("qns_test");
    const std::string path = tmp.file("model.jmva");
    io::write_jmva(m.get_struct(), path, "default", 10000);
    const std::string doc = read_all(path);

    CHECK(doc.find("<closedclass population=\"4\" name=\"Chain01\"/>") != std::string::npos);
    CHECK(doc.find("<delaystation name=\"Delay1\">") != std::string::npos);
    // A single-server Queue is an li station: the ld encoding is what a
    // multiserver one gets, and emitting it here would spell out a rate vector
    // for a station that has one rate.
    CHECK(doc.find("<listation name=\"Queue1\"") != std::string::npos);
    CHECK(doc.find("ldstation") == std::string::npos);
    CHECK(doc.find("refStation=\"Delay1\"") != std::string::npos);
    CHECK(doc.find("<algType name=\"MVA\"") != std::string::npos);
    // The Sink is not a station and the Source is dropped, so a closed model's
    // station count is every station it has.
    CHECK(doc.find("<stations number=\"2\">") != std::string::npos);
}

TEST_CASE("writeJMVA spells a multiserver queue out as a load-dependent station") {
    qn::Network<double> m = closed_cqn(1, 3.0, 2.0);
    util::TempDir tmp("qns_test");
    const std::string path = tmp.file("model.jmva");
    io::write_jmva(m.get_struct(), path, "default", 10000);
    const std::string doc = read_all(path);

    CHECK(doc.find("<ldstation name=\"Queue1\"") != std::string::npos);
    CHECK(doc.find("listation") == std::string::npos);
    // The vector runs to the closed population: S, S/2, S/2 for two servers and
    // three jobs, since S/min(n,c) is constant past the server count.
    const std::size_t at = doc.find("<servicetimes customerclass=\"Chain01\">");
    REQUIRE(at != std::string::npos);
    const std::string row = doc.substr(at, doc.find("</servicetimes>", at) - at);
    CHECK(std::count(row.begin(), row.end(), ';') == 2);
}

TEST_CASE("writeJMVA substitutes the reference station of an open chain") {
    qn::Network<double> m = open_oqn(0.5);
    util::TempDir tmp("qns_test");
    const std::string path = tmp.file("model.jmva");
    io::write_jmva(m.get_struct(), path, "default", 10000);
    const std::string doc = read_all(path);

    // The open chain's reference station is the Source, which the document does
    // not carry; naming it would make qnsolver reject the whole model.
    CHECK(doc.find("<openclass rate=\"0.5\" name=\"Chain01\"/>") != std::string::npos);
    CHECK(doc.find("refStation=\"Source\"") == std::string::npos);
    CHECK(doc.find("refStation=\"Queue1\"") != std::string::npos);
    CHECK(doc.find("<stations number=\"2\">") != std::string::npos);
}

// ---------------------------------------------------------------------------
// Method resolution and the refusals, none of which needs the binary
// ---------------------------------------------------------------------------

TEST_CASE("the method names resolve as the reference resolves them") {
    CHECK(qns::list_valid_methods().size() == 7);
    qns::QnsOptions o;
    // `default` is rolia, which is what runAnalyzer.m maps it to.
    CHECK(qns::resolve_multiserver(o) == "rolia");
    o.multiserver = "conway";
    CHECK(qns::resolve_multiserver(o) == "conway");
    // A named method WINS over the config, as the reference's switch overwrites
    // options.config.multiserver.
    o.method = "zhou";
    CHECK(qns::resolve_multiserver(o) == "zhou");
    CHECK(qns::is_qnsolver_multiserver("reiser"));
    CHECK_FALSE(qns::is_qnsolver_multiserver("suri"));
    CHECK_FALSE(qns::is_qnsolver_multiserver("schmidt"));
}

TEST_CASE("an unknown method is refused by name") {
    CHECK_THROWS_AS(qns::check_method("linearizer"), InputError);
    CHECK_NOTHROW(qns::check_method("conway"));
    CHECK_NOTHROW(qns::check_method("schmidt"));
}

TEST_CASE("suri and schmidt are refused rather than run as qnsolver's default") {
    qn::Network<double> m = closed_cqn(1, 4.0, 2.0);
    qns::QnsOptions o;
    o.method = "suri";
    CHECK_THROWS_AS(qns::solver_qns_run_analyzer(m.get_struct(), o), UnsupportedError);
    o.method = "schmidt";
    CHECK_THROWS_AS(qns::solver_qns_run_analyzer(m.get_struct(), o), UnsupportedError);
}

TEST_CASE("a binding station buffer is refused on both routes") {
    // NOTHING under the QNS tree reads `sn.cap` or `sn.classcap`. The JMVA
    // document is read by `qnsolver`, whose MVA-family algorithms have no
    // representation of a finite buffer, and the QN2LQN route hands the model
    // to LQNS, which has none either -- so a capped station was solved as an
    // unbounded one and the table reported the unconstrained answer under this
    // solver's name. The gate sits in `solver_qns_run_analyzer`, ahead of the split,
    // rather than in `check_supported`, which is the narrower JMVA-document
    // gate that the layered route deliberately skips.
    //
    // PS keeps the model product-form (the JMVA route); FCFS with unequal rates
    // takes it off product form and down the QN2LQN route. Both must refuse.
    for (SchedStrategy sched : {SchedStrategy::PS, SchedStrategy::FCFS}) {
        qn::Network<double> m("capped");
        const std::size_t q1 = m.add_queue("Queue1", sched);
        const std::size_t q2 = m.add_queue("Queue2", sched);
        const std::size_t k1 = m.add_closed_class("C1", 3.0, q1);
        const std::size_t k2 = m.add_closed_class("C2", 2.0, q1);
        m.set_service(q1, k1, D::exp_rate(1.0));
        m.set_service(q2, k1, D::exp_rate(0.8));
        m.set_service(q1, k2, D::exp_rate(0.5));
        m.set_service(q2, k2, D::exp_rate(1.3));
        m.set_capacity(q2, 1.0);
        qn::RoutingMatrix<double> P;
        P.set(k1, k1, q1, q2, 1.0);
        P.set(k1, k1, q2, q1, 1.0);
        P.set(k2, k2, q1, q2, 1.0);
        P.set(k2, k2, q2, q1, 1.0);
        m.link(P);

        qns::QnsOptions o;
        std::string what;
        try {
            qns::solver_qns_run_analyzer(m.get_struct(), o);
        } catch (const std::exception& e) {
            what = e.what();
        }
        CHECK(what.find("SolverQNS") != std::string::npos);
        CHECK(what.find("finite station capacity") != std::string::npos);
        CHECK(what.find("Queue2") != std::string::npos);
    }
}

TEST_CASE("a capacity that cannot bind is not refused") {
    // Only a buffer that can ACTUALLY bind is a refusal: a closed model whose
    // station capacity is at least the total population can never block a job,
    // so the declaration is a no-op and the answer stays exact. `setCap(N)` on
    // an order-independent station of an N-job model is a common idiom.
    qn::Network<double> m = closed_cqn(1, 4.0);
    qn::NetworkStruct<double> sn = m.get_struct();
    for (qn::Station<double>& st : sn.stations)
        if (st.nodetype == qn::NodeType::Queue) st.cap = 4.0;
    CHECK_NOTHROW(qns::check_method("conway"));
    std::string what;
    try {
        qns::QnsOptions o;
        qns::solver_qns_run_analyzer(sn, o);
    } catch (const std::exception& e) {
        what = e.what();
    }
    CHECK(what.find("finite station capacity") == std::string::npos);
}

TEST_CASE("QN2LQN maps a heterogeneous closed network to reference and service activities") {
    // Heterogeneous FCFS service is the canonical non-product-form branch of
    // SolverQNS.  Each class is its own chain in this fixture.
    qn::Network<double> m = closed_cqn(2, 3.0, 1.0, SchedStrategy::FCFS);
    REQUIRE_FALSE(m.get_struct().has_product_form());
    const lqn::LqnModel<double> lm = io::qn2lqn(m);
    const lqn::LqnStruct<double> ls = lqn::lqn_finalize(lm);

    CHECK(lm.procs.size() == 3);    // pseudo host, Delay1, Queue1
    CHECK(lm.tasks.size() == 4);    // two reference tasks, two service tasks
    CHECK(lm.entries.size() == 6);  // two chain entries, four service entries
    CHECK(ls.nacts == 10);          // four Q, four A and two terminal End activities

    const lqn::detail::RawActivity<double>& q11 = raw_activity(lm, "Q1_1");
    CHECK(q11.bound_to_entry == "E1_1");
    CHECK(std::fabs(q11.hostdem.mean - 1.0) < 1e-12);
    const lqn::detail::RawActivity<double>& q22 = raw_activity(lm, "Q2_2");
    CHECK(q22.bound_to_entry == "E2_2");
    CHECK(std::fabs(q22.hostdem.mean - 1.0 / 3.0) < 1e-12);

    const lqn::detail::RawActivity<double>& ref = raw_activity(lm, "A1_1");
    CHECK(ref.bound_to_entry == "Chain_1");
    REQUIRE(ref.sync_calls.size() == 1);
    CHECK(ref.sync_calls[0].dest == "E1_1");
    CHECK(ref.sync_calls[0].mean == 1.0);

    const lqn::detail::RawPrecedence<double>& cycle =
        raw_precedence(lm, "A2_1", lang::PrecedenceType::POST_OR);
    REQUIRE(cycle.postacts.size() == 1);
    CHECK(cycle.postacts[0] == "End_1_2_1");
    CHECK(cycle.postparams[0] == 1.0);
}

TEST_CASE("QN2LQN preserves Router probabilities as an OR fork") {
    qn::Network<double> m("router");
    const std::size_t d = m.add_delay("Think");
    const std::size_t r = m.add_router("R");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C", 2.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(4.0));
    qn::RoutingMatrix<double> P;
    P.set(d, r, 1.0);
    P.set(r, q1, 0.25);
    P.set(r, q2, 0.75);
    P.set(q1, d, 1.0);
    P.set(q2, d, 1.0);
    m.link(P);

    const lqn::LqnModel<double> lm = io::qn2lqn(m.get_struct());
    const lqn::detail::RawPrecedence<double>& split =
        raw_precedence(lm, "CS_1_2_1", lang::PrecedenceType::POST_OR);
    REQUIRE(split.postacts.size() == 2);
    REQUIRE(split.postparams.size() == 2);
    CHECK(split.postacts[0] == "A3_1");
    CHECK(split.postacts[1] == "A4_1");
    CHECK(split.postparams[0] == doctest::Approx(0.25));
    CHECK(split.postparams[1] == doctest::Approx(0.75));
}

TEST_CASE("QN2LQN maps Fork and Join to AND precedences") {
    qn::Network<double> m("forkjoin");
    const std::size_t d = m.add_delay("Think");
    const std::size_t f = m.add_fork("Fork");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t j = m.add_join("Join", f);
    const std::size_t c = m.add_closed_class("C", 2.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(j, d, 1.0);
    m.link(P);

    const lqn::LqnModel<double> lm = io::qn2lqn(m.get_struct());
    const lqn::detail::RawPrecedence<double>& fork =
        raw_precedence(lm, "FJ_1_2_1", lang::PrecedenceType::POST_AND);
    REQUIRE(fork.postacts.size() == 2);
    CHECK(fork.postacts[0] == "A3_1");
    CHECK(fork.postacts[1] == "A4_1");

    bool saw_join = false;
    for (std::size_t t = 0; t < lm.tasks.size(); ++t)
        for (std::size_t p = 0; p < lm.tasks[t].precedences.size(); ++p) {
            const lqn::detail::RawPrecedence<double>& pr = lm.tasks[t].precedences[p];
            if (pr.pretype != lang::PrecedenceType::PRE_AND || pr.postacts.size() != 1 ||
                pr.postacts[0] != "FJ_1_5_1")
                continue;
            REQUIRE(pr.preacts.size() == 2);
            CHECK(pr.preacts[0] == "A3_1");
            CHECK(pr.preacts[1] == "A4_1");
            saw_join = true;
        }
    CHECK(saw_join);
}

TEST_CASE("the non-product-form QNS branch now reaches SolverLQNS") {
    if (lqns::SolverLQNS<double>::is_available()) return;
    qn::Network<double> m = closed_cqn(2, 3.0, 1.0, SchedStrategy::FCFS);
    try {
        (void)qns::solver_qns_run_analyzer(m.get_struct(), qns::QnsOptions());
        FAIL("an absent LQNS binary must be diagnosed");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("SolverLQNS requires") != std::string::npos);
    }
}

TEST_CASE("the non-product-form QNS branch solves through an installed LQNS") {
    if (!lqns::SolverLQNS<double>::is_available()) return;
    qn::Network<double> m = closed_cqn(2, 3.0, 1.0, SchedStrategy::FCFS);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    REQUIRE_FALSE(sn.has_product_form());
    const mva::AvgResult<double> got = qns::solver_qns_run_analyzer(sn, qns::QnsOptions());
    CHECK(got.actualmethod == "default/rolia");
    obeys_the_laws(sn, got);
}

TEST_CASE("a station the JMVA document cannot carry is refused, not dropped") {
    qn::Network<double> m("spn");
    const std::size_t p = m.add_place("Place1");
    const std::size_t d = m.add_delay("Delay1");
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, p, 1.0);
    P.set(c, c, p, d, 1.0);
    m.link(P);
    CHECK_THROWS_AS(qns::check_supported(m.get_struct()), UnsupportedError);
}

// ---------------------------------------------------------------------------
// End to end, against exact MVA
// ---------------------------------------------------------------------------

TEST_CASE("a single-class closed solve obeys the product-form laws") {
    if (!qns::is_available()) return;  // no binary: nothing to run
    qn::Network<double> m = closed_cqn(1, 4.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const mva::AvgResult<double> got = qns::solver_qns_run_analyzer(sn, qns::QnsOptions());
    obeys_the_laws(sn, got);

    // The Linearizer is what `-l` selects, so the agreement with THIS port's
    // Linearizer is far tighter than with exact MVA. Both are asserted: the
    // first pins which algorithm ran, the second that it ran on the right model.
    tables_agree(run_mva(sn, "lin").QN, got.QN, 1e-3, "QLen vs own linearizer");
    tables_agree(run_mva(sn, "exact").QN, got.QN, 1e-2, "QLen vs exact");
    // No multiserver station, so no -m flag is passed and the label says so.
    CHECK(got.actualmethod == "default/default");
}

TEST_CASE("a two-class closed solve obeys the product-form laws") {
    if (!qns::is_available()) return;
    qn::Network<double> m = closed_cqn(2, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    REQUIRE(sn.nchains == 2);
    const mva::AvgResult<double> got = qns::solver_qns_run_analyzer(sn, qns::QnsOptions());
    obeys_the_laws(sn, got);
    tables_agree(run_mva(sn, "lin").QN, got.QN, 1e-2, "QLen vs own linearizer");
    tables_agree(run_mva(sn, "exact").QN, got.QN, 2e-2, "QLen vs exact");
}

TEST_CASE("the open-network solve recovers the chain throughput without a Source row") {
    if (!qns::is_available()) return;
    qn::Network<double> m = open_oqn(0.5);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const mva::AvgResult<double> want = run_mva(sn, "exact");
    const mva::AvgResult<double> got = qns::solver_qns_run_analyzer(sn, qns::QnsOptions());

    // An open product-form network is exact under every MVA variant, so here the
    // agreement IS tight -- capped only by the six digits qnsolver prints. The
    // Source is not in the document, so the chain throughput has to come back
    // from a station row through its visit count; getting that wrong scales
    // every metric and would show here first.
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        if (sn.stations[i].nodetype == qn::NodeType::Source) continue;
        INFO("station " << i);
        CHECK(std::fabs(want.QN(i, 0) - got.QN(i, 0)) <= 1e-5 * std::max(1.0, want.QN(i, 0)));
        CHECK(std::fabs(want.TN(i, 0) - got.TN(i, 0)) <= 1e-5 * std::max(1.0, want.TN(i, 0)));
        CHECK(std::fabs(got.TN(i, 0) - 0.5) <= 1e-5);  // every station sees lambda
    }
}

TEST_CASE("a multiserver model reports the approximation that actually ran") {
    if (!qns::is_available()) return;
    qn::Network<double> m = closed_cqn(1, 6.0, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    qns::QnsOptions o;
    const mva::AvgResult<double> def = qns::solver_qns_run_analyzer(sn, o);
    CHECK(def.actualmethod == "default/rolia");
    o.method = "conway";
    const mva::AvgResult<double> con = qns::solver_qns_run_analyzer(sn, o);
    CHECK(con.actualmethod == "conway");

    // qnsolver reports an ld station's $U as the mean number of BUSY SERVERS
    // (X*S, measured directly against the binary), so the divide by the server
    // count is what turns it into the per-server fraction LINE reports. Without
    // it a three-server station would report a utilization near 2.
    obeys_the_laws(sn, con);
    obeys_the_laws(sn, def);
}
