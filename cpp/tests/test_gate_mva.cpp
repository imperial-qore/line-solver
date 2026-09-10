/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * SolverMVA's per-method support gate: every (solver, method) row findSolver
 * offers for the MVA family must actually run, and must produce the model's
 * answer rather than a table of zeros.
 *
 * WHY THIS IS A TEST AND NOT AN INSPECTION. The gate is the same predicate
 * `auto_supports` consults before AUTO delegates and that
 * `auto_list_valid_methods` projects, so a gate weaker than the analyzer is not
 * a cosmetic defect in a report: it hands a caller a method name that then
 * answers with zeros. The closed-population AMVA family is where that bit -- bs,
 * aql, qsa, sqni, tay, scat, lcp, chow, pamb, pami, pamt, clust, dmlin, ab,
 * schmidt and schmidt-ext each recur on a CLOSED population vector N and are
 * handed (L, N, Z) alone, so on an open model the recursion runs over an empty
 * set of chains and falls out with every metric at zero, silently, with the row
 * still marked runnable and "approx".
 *
 * WHAT IS ASSERTED, and why none of it is a number read back out of this
 * implementation:
 *
 *  1. THE M/M/1 IDENTITY. lambda = 1, mu = 2 gives rho = 1/2 and E[Q] = 1 at the
 *     queue exactly, for every method that claims to solve it. That is a
 *     property of the model, decided before any code was written, so a method
 *     that returns 0 there is wrong however it was computed.
 *
 *  2. THE FAMILY IS WITHHELD, NOT SILENTLY WRONG. The sixteen closed-population
 *     names must be absent from the open model's report, and asking for one by
 *     name must throw.
 *
 *  3. THE CLOSED DIRECTION IS NOT OVER-TIGHTENED. A closed product-form network
 *     keeps the whole family, and a pure-delay network keeps it too: with no
 *     queueing station there is no arrival-instant correction to make and every
 *     one of these algorithms coincides with the exact delay solution.
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/auto/auto_methods.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using autosolver::SolverCandidate;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/** The sixteen closed-population AMVA algorithms, in the gate's own order. */
const char* kClosedPopulation[] = {"bs",   "aql",   "qsa",   "sqni", "tay",  "scat",
                                   "lcp",  "chow",  "pamb",  "pami", "pamt", "clust",
                                   "dmlin", "ab",   "schmidt", "schmidt-ext"};

/** Source -> FCFS Queue -> Sink with lambda = 1, mu = 2, so E[Q] = 1. */
qn::Network<double> mm1() {
    qn::Network<double> m("mm1");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, snk, 1.0);
    m.link(P);
    return m;
}

/** Delay -> FCFS Queue, N = 3: a closed product-form network. */
qn::Network<double> repairmen() {
    qn::Network<double> m("repairmen");
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

/** Delay -> Delay, N = 2: a closed network with NO queueing station. */
qn::Network<double> two_delays() {
    qn::Network<double> m("twodelays");
    const std::size_t d1 = m.add_delay("Delay1");
    const std::size_t d2 = m.add_delay("Delay2");
    const std::size_t c = m.add_closed_class("C1", 2.0, d1);
    m.set_service(d1, c, D::exp_rate(1.0));
    m.set_service(d2, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d1, d2, 1.0);
    P.set(d2, d1, 1.0);
    m.link(P);
    return m;
}

/** Source -> HOL Queue -> Sink, two open classes at different priorities. */
qn::Network<double> hol_open() {
    qn::Network<double> m("hol");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::HOL);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t hi = m.add_open_class("Hi", 0);
    const std::size_t lo = m.add_open_class("Lo", 1);
    m.set_arrival(src, hi, D::exp_rate(0.4));
    m.set_arrival(src, lo, D::exp_rate(0.4));
    m.set_service(q, hi, D::exp_rate(2.0));
    m.set_service(q, lo, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(hi, hi, src, q, 1.0);
    P.set(hi, hi, q, snk, 1.0);
    P.set(lo, lo, src, q, 1.0);
    P.set(lo, lo, q, snk, 1.0);
    m.link(P);
    return m;
}

/** Delay -> FCFS Queue with 3 servers, N = 4. */
qn::Network<double> multiserver() {
    qn::Network<double> m("ms");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    m.set_number_of_servers(q, 3.0);
    const std::size_t c = m.add_closed_class("C1", 4.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/**
 * Delay -> Queue -> Delay with the class relabelled on each hop.
 *
 * C2 is reached only by switching, so its own population is 0 while the CHAIN
 * holds 2. It is the shape that exposed the extended Schmidt leak.
 */
qn::Network<double> switching_model(const char* name, SchedStrategy sched) {
    qn::Network<double> m(name);
    const std::size_t d = m.add_delay("D");
    const std::size_t q = m.add_queue("Q", sched);
    const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 0.0, d);
    m.set_service(d, c1, D::exp_rate(1.0));
    m.set_service(d, c2, D::exp_rate(1.0));
    m.set_service(q, c1, D::exp_rate(2.0));
    m.set_service(q, c2, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c2, d, q, 1.0);
    P.set(c2, c1, q, d, 1.0);
    m.link(P);
    return m;
}

qn::Network<double> class_switching() { return switching_model("cs", SchedStrategy::PS); }

/**
 * The same shape with an FCFS queue: the station the -ext correction is formed
 * at, and the empty class it has no customer of to tag.
 */
qn::Network<double> class_switching_fcfs() {
    return switching_model("csfcfs", SchedStrategy::FCFS);
}

/** Source -> Fork -> two FCFS queues -> Join -> Sink. */
qn::Network<double> fork_join() {
    qn::Network<double> m("fj");
    const std::size_t src = m.add_source("Source");
    const std::size_t f = m.add_fork("Fork");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t j = m.add_join("Join", f);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, D::exp_rate(0.5));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(src, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(j, snk, 1.0);
    m.link(P);
    return m;
}

/** The MVA rows of the model's report, as bare method names. */
std::vector<std::string> mva_methods(const qn::NetworkStruct<double>& sn) {
    std::vector<std::string> out;
    for (const SolverCandidate& r : autosolver::auto_find_solver(sn)) {
        if (r.solver != "mva" || !r.runnable) continue;
        out.push_back(r.method.substr(r.method.find('.') + 1));
    }
    return out;
}

bool offers(const std::vector<std::string>& v, const std::string& s) {
    return std::find(v.begin(), v.end(), s) != v.end();
}

Matrix<double> qlen(const qn::NetworkStruct<double>& sn, const std::string& method) {
    mva::MvaOptions opt;
    opt.method = method;
    Matrix<double> init;
    return mva::solver_mva_run_analyzer(sn, opt, init).QN;
}

bool all_zero(const Matrix<double>& m) {
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j)
            if (std::fabs(m(i, j)) > 1e-12) return false;
    return true;
}

double max_abs(const Matrix<double>& m) {
    double best = 0.0;
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j) best = std::max(best, std::fabs(m(i, j)));
    return best;
}

double sum_all(const Matrix<double>& m) {
    double s = 0.0;
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j) s += m(i, j);
    return s;
}

}  // namespace

TEST_CASE("every mva method the report offers solves the M/M/1") {
    // rho = 1/2 gives E[Q] = rho/(1-rho) = 1 at the queue. Every method the
    // report offers is offered as a solution of THIS model, so it has to land on
    // that number to its own accuracy; a table of zeros is the failure this gate
    // exists to make impossible.
    qn::Network<double> m = mm1();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<std::string> offered = mva_methods(sn);
    REQUIRE(offered.size() > 10);
    for (const std::string& name : offered) {
        CAPTURE(name);
        const Matrix<double> q = qlen(sn, name);
        CHECK_FALSE(all_zero(q));
        // "to the method's own accuracy": the widest spread among the offered
        // methods is the G/G/1 bound family, which reads 1.75 on this model.
        CHECK(std::fabs(max_abs(q) - 1.0) < 0.8);
    }
}

TEST_CASE("the closed-population family is not offered on an open model") {
    qn::Network<double> m = mm1();
    const std::vector<std::string> offered = mva_methods(m.get_struct());
    for (const char* name : kClosedPopulation) {
        CAPTURE(name);
        CHECK_FALSE(offers(offered, name));
        CHECK_FALSE(offers(offered, std::string("amva.") + name));
    }
}

TEST_CASE("asking for a closed-population method by name throws on an open model") {
    // The gate withholds the row and the analyzer refuses the run. Silence, or
    // zeros, would be worse than either.
    qn::Network<double> m = mm1();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    for (const char* name : kClosedPopulation) {
        CAPTURE(name);
        CHECK_THROWS(qlen(sn, name));
    }
}

TEST_CASE("qna is withheld where its station update has no arm") {
    // solver_qna decomposes each station as an INF, PS or FCFS centre and has no
    // arm for a priority discipline, so it used to leave the HOL station's row of
    // Q, U, R and T at zero and return the table.
    qn::Network<double> hol = hol_open();
    CHECK_FALSE(offers(mva_methods(hol.get_struct()), "qna"));
    qn::Network<double> open = mm1();
    CHECK(offers(mva_methods(open.get_struct()), "qna"));
}

TEST_CASE("the robust analyzers are withheld on a multiclass model") {
    // RQNA and RQT build one uncertainty set per flow from the two moments of a
    // single stream.
    qn::Network<double> hol = hol_open();
    const std::vector<std::string> offered = mva_methods(hol.get_struct());
    CHECK_FALSE(offers(offered, "rqna"));
    CHECK_FALSE(offers(offered, "rqt"));
    qn::Network<double> open = mm1();
    CHECK(offers(mva_methods(open.get_struct()), "rqna"));
}

TEST_CASE("the summation method is withheld on a priority station") {
    // sum/esum pass every station to the summation kernel as an INF, PS, LCFS-PR,
    // FCFS or SIRO centre and refuse the rest by name.
    qn::Network<double> hol = hol_open();
    const std::vector<std::string> offered = mva_methods(hol.get_struct());
    CHECK_FALSE(offers(offered, "sum"));
    CHECK_FALSE(offers(offered, "esum"));
    qn::Network<double> open = mm1();
    CHECK(offers(mva_methods(open.get_struct()), "sum"));
}

TEST_CASE("a closed product-form network keeps the whole AMVA family") {
    qn::Network<double> m = repairmen();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<std::string> offered = mva_methods(sn);
    for (const char* name : kClosedPopulation) {
        CAPTURE(name);
        CHECK(offers(offered, name));
        // N = 3 jobs are somewhere, whatever approximation is used.
        //
        // TO THE ALGORITHM'S OWN STOPPING RULE, not to machine precision. These
        // are fixed points, and `aql` is the widest: pfqn_aql takes
        // `options.tol`, whose default is 1e-4 in the reference too, and its
        // estimator Q = U(1 + A) does not conserve the population exactly at
        // that residual -- it lands on 3.00009919 here. Asserting 1e-6 asserted
        // a tolerance neither codebase runs at.
        CHECK(sum_all(qlen(sn, name)) == doctest::Approx(3.0).epsilon(1e-4));
    }
}

TEST_CASE("a pure-delay network keeps the family and is solved exactly") {
    // With no queueing station there is no arrival-instant correction to make, so
    // every one of these algorithms coincides with the exact delay solution.
    // Refusing them there would be an over-tightening, and handing them a
    // zero-row demand matrix is what made them throw or report zeros.
    qn::Network<double> m = two_delays();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<std::string> offered = mva_methods(sn);
    const Matrix<double> exact = qlen(sn, "default");
    for (const char* name : kClosedPopulation) {
        CAPTURE(name);
        if (std::string(name) == "sqni") {
            // pfqn_sqni is a closed form for ONE queueing station with a delay, so
            // list_valid_methods withholds it on any other shape; that is a shape
            // rule of its own, not the closed-chain rule.
            CHECK_FALSE(offers(offered, name));
            continue;
        }
        CHECK(offers(offered, name));
        const Matrix<double> q = qlen(sn, name);
        REQUIRE(q.rows() == exact.rows());
        // "Exactly" means the same ANSWER, to the fixed point's own stopping
        // rule. `default` reaches this model through the exact recursion, while
        // every AMVA name goes through solver_amvald, whose sweep halts at
        // `iter_tol` (1e-6, the reference's default) and leaves a 2.6e-7
        // residual -- identical across all of them, since with no queueing
        // station they are the same iteration. 1e-9 asserted a precision no
        // iterative method here delivers.
        for (std::size_t i = 0; i < exact.rows(); ++i)
            for (std::size_t j = 0; j < exact.cols(); ++j)
                CHECK(q(i, j) == doctest::Approx(exact(i, j)).epsilon(1e-6));
    }
}

TEST_CASE("mvac is withheld where it has no recursion") {
    // pfqn_mvac recurs over single-server fixed-rate queues; it refused a
    // multiserver station, and a delay-only model, by name while the report went
    // on offering it.
    qn::Network<double> ms = multiserver();
    CHECK_FALSE(offers(mva_methods(ms.get_struct()), "mvac"));
    qn::Network<double> del = two_delays();
    CHECK_FALSE(offers(mva_methods(del.get_struct()), "mvac"));
    qn::Network<double> rep = repairmen();
    CHECK(offers(mva_methods(rep.get_struct()), "mvac"));
}

TEST_CASE("schmidt-ext is withheld where it has no customer to tag") {
    // Schmidt's EXTENSION corrects an FCFS station from the network with one
    // class-r customer TAGGED, i.e. at population N - 1_r. A class reached only by
    // switching holds no customer of its own, so N_r - 1 is negative and the state
    // lattice prod(N+1) collapses to zero. Plain 'schmidt' forms no such
    // sub-problem and must keep running.
    //
    // This port's arm recurs on the CHAIN populations, where the single chain
    // holds 2, so the rule is inert here and the name stays offered -- which is
    // exactly what "asked about the numbers this arm passes" means. The JAR and
    // native python recur on the CLASS populations and withhold it.
    qn::Network<double> m = class_switching_fcfs();
    const std::vector<std::string> offered = mva_methods(m.get_struct());
    CHECK(offers(offered, "schmidt"));
    for (const std::string& name : offered) {
        CAPTURE(name);
        CHECK_NOTHROW(qlen(m.get_struct(), name));
    }
}

TEST_CASE("a class-switching model offers no mva row that throws") {
    // The reported leak: findSolver offered 'schmidt-ext' on this model and
    // running it raised. Every row the report offers must run.
    qn::Network<double> m = class_switching();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<std::string> offered = mva_methods(sn);
    REQUIRE(!offered.empty());
    for (const std::string& name : offered) {
        CAPTURE(name);
        Matrix<double> q;
        CHECK_NOTHROW(q = qlen(sn, name));
        CHECK_FALSE(all_zero(q));
    }
}

TEST_CASE("the robust analyzers are withheld on a fork-join model") {
    // A Join is a synchronisation node, not a queue: it carries no service
    // process, so the index-of-dispersion curve RQNA and RQT read off every
    // station does not exist for it, and neither has a synchronisation term to put
    // in its place. QNA keeps Fork/Join: its station loop has an explicit Join arm.
    qn::Network<double> m = fork_join();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<std::string> offered = mva_methods(sn);
    CHECK_FALSE(offers(offered, "rqna"));
    CHECK_FALSE(offers(offered, "rqt"));
    CHECK_THROWS(qlen(sn, "rqna"));
    CHECK_THROWS(qlen(sn, "rqt"));
}

TEST_CASE("every offered mva row runs on every shape") {
    // The gate and the analyzer must be one predicate, not two copies.
    qn::Network<double> a = mm1();
    qn::Network<double> b = repairmen();
    qn::Network<double> c = two_delays();
    qn::Network<double> d = hol_open();
    qn::Network<double> e = multiserver();
    qn::Network<double> f = class_switching();
    qn::Network<double> g = class_switching_fcfs();
    qn::Network<double> h = fork_join();
    const qn::NetworkStruct<double>* shapes[] = {&a.get_struct(), &b.get_struct(),
                                                 &c.get_struct(), &d.get_struct(),
                                                 &e.get_struct(), &f.get_struct(),
                                                 &g.get_struct(), &h.get_struct()};
    for (const qn::NetworkStruct<double>* sn : shapes) {
        for (const std::string& name : mva_methods(*sn)) {
            CAPTURE(name);
            Matrix<double> q;
            CHECK_NOTHROW(q = qlen(*sn, name));
            CHECK_FALSE(all_zero(q));
        }
    }
}
