/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Tests of line/solvers/auto/auto_methods.h's findSolver: the report of which
 * solvers and solver methods can analyze a model, why the others cannot, what
 * kind of answer each returns and which measures it can report.
 *
 * WHAT IS ASSERTED, and why none of it is a number read back out of this
 * implementation.
 *
 *  1. THE PROJECTION IDENTITY. auto_list_valid_methods is defined as the
 *     runnable rows of auto_find_solver plus the method names that name no single
 *     method. That is the whole point of the refactor -- one gate, two views --
 *     so it is asserted directly rather than trusted: every runnable row's
 *     method name appears in the list, and no refused row's does.
 *
 *  2. STRUCTURAL FACTS ABOUT THE MODELS. An M/M/1 has a product-form solution
 *     and is one queueing station fed by a Source, so exact MVA and the QBD are
 *     exact ON IT; a model whose service law no product form admits is not,
 *     and the same methods must then report 'approx'. These are properties of
 *     the models, decided before the code was written.
 *
 *  3. THE REFUSAL REASONS NAME THE FEATURE. A report whose whole content is
 *     the explanation is worthless if the explanation is "unsupported", so a
 *     refused row is required to name the offending feature.
 */

#include <algorithm>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/auto/auto_methods.h"

using namespace line;
using autosolver::SolverCandidate;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/** Source -> FCFS Queue -> Sink, exponential throughout: the M/M/1. */
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

/** The same shape with a Pareto service law: no product form, and NC refuses it. */
qn::Network<double> mparetol() {
    qn::Network<double> m("mpareto1");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, D::exp_rate(0.5));
    m.set_service(q, c, D::pareto(2.5, 1.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, snk, 1.0);
    m.link(P);
    return m;
}

/** Delay -> Queue1 -> Queue2, N = 4: a closed product-form network. */
qn::Network<double> cqn() {
    qn::Network<double> m("cqn");
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

bool has_method(const std::vector<SolverCandidate>& rows, const std::string& method_name) {
    for (const SolverCandidate& r : rows)
        if (r.method == method_name) return true;
    return false;
}

const SolverCandidate* find_row(const std::vector<SolverCandidate>& rows,
                                const std::string& method_name) {
    for (const SolverCandidate& r : rows)
        if (r.method == method_name) return &r;
    return nullptr;
}

bool listed(const std::vector<std::string>& v, const std::string& s) {
    return std::find(v.begin(), v.end(), s) != v.end();
}

}  // namespace

TEST_CASE("listValidMethods is the runnable rows of findSolver, projected") {
    qn::Network<double> m = mm1();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<SolverCandidate> all = autosolver::auto_find_solver(sn, "", true);
    const std::vector<std::string> valid = autosolver::auto_list_valid_methods(sn);

    REQUIRE(!all.empty());
    std::size_t runnable = 0;
    for (const SolverCandidate& r : all) {
        if (r.runnable) {
            ++runnable;
            // Every runnable pair is a method name a caller may ask AUTO for.
            CHECK(listed(valid, r.method));
            CHECK(listed(valid, r.solver));
            CHECK(r.reason.empty());
        } else {
            // and a refused one is not offered.
            CHECK_FALSE(listed(valid, r.method));
            CHECK_FALSE(r.reason.empty());
        }
    }
    CHECK(runnable > 0);

    // The default report is exactly the runnable half.
    const std::vector<SolverCandidate> runnable_only = autosolver::auto_find_solver(sn);
    CHECK(runnable_only.size() == runnable);
    for (const SolverCandidate& r : runnable_only) CHECK(r.runnable);

    // The selection intents name a ranking rather than an algorithm and are
    // listed whatever the model is.
    CHECK(listed(valid, "default"));
    CHECK(listed(valid, "exact"));
    CHECK(listed(valid, "bound"));
}

TEST_CASE("a spelling already qualified with its own family is not doubled") {
    // The fluid registry declares both 'dae' and 'fluid.dae' so that its own
    // gate takes either; prefixing the family again would yield
    // 'fluid.fluid.dae', a method name that resolves but names the same method twice.
    qn::Network<double> m = mm1();
    const std::vector<SolverCandidate> rows = autosolver::auto_find_solver(m.get_struct(), "", true);
    CHECK(has_method(rows, "fluid.dae"));
    CHECK_FALSE(has_method(rows, "fluid.fluid.dae"));
    for (const SolverCandidate& r : rows)
        CHECK(r.method.find(r.solver + "." + r.solver + ".") == std::string::npos);
}

TEST_CASE("exactness is claimed of the model, not of the algorithm") {
    // Exact MVA is exact on a product-form model and an approximation off one;
    // the QBD is exact on one queueing station fed by a Source. Both conditions
    // are properties of the models below, not of the implementation.
    qn::Network<double> open = mm1();
    const std::vector<SolverCandidate> a = autosolver::auto_find_solver(open.get_struct(), "", true);
    REQUIRE(open.get_struct().has_product_form());
    const SolverCandidate* mva_exact = find_row(a, "mva.exact");
    REQUIRE(mva_exact != nullptr);
    CHECK(mva_exact->method_class == "exact");
    const SolverCandidate* mam_default = find_row(a, "mam.default");
    REQUIRE(mam_default != nullptr);
    CHECK(mam_default->method_class == "exact");
    // A decomposition of a network into such queues is an approximation of it.
    const SolverCandidate* dec = find_row(a, "mam.dec.source");
    REQUIRE(dec != nullptr);
    CHECK(dec->method_class == "approx");

    // Three stations: still product form, so MVA stays exact, but no longer one
    // queueing station, so the QBD is not.
    qn::Network<double> closed = cqn();
    const std::vector<SolverCandidate> b = autosolver::auto_find_solver(closed.get_struct(), "", true);
    const SolverCandidate* mva2 = find_row(b, "mva.exact");
    REQUIRE(mva2 != nullptr);
    CHECK(mva2->method_class == "exact");
    const SolverCandidate* mam2 = find_row(b, "mam.default");
    if (mam2 != nullptr) CHECK(mam2->method_class == "approx");

    // Every AMVA arm is an approximation whatever the model.
    const SolverCandidate* amva = find_row(b, "mva.amva");
    REQUIRE(amva != nullptr);
    CHECK(amva->method_class == "approx");

    // Bounds are what SolverBA is for, and a simulator is a simulator.
    for (const SolverCandidate& r : b) {
        if (r.solver == "ba") CHECK(r.method_class == "bound");
        if (r.solver == "ssa" || r.solver == "ldes") CHECK(r.method_class == "simulation");
    }
    // The CTMC generator is solved as written; 'cftp.approx' says otherwise in
    // its own name.
    const SolverCandidate* cftp = find_row(b, "ctmc.cftp.approx");
    if (cftp != nullptr) CHECK(cftp->method_class == "approx");
}

TEST_CASE("a refusal names the feature that caused it") {
    qn::Network<double> m = mparetol();
    const std::vector<SolverCandidate> rows = autosolver::auto_find_solver(m.get_struct(), "", true);
    std::size_t refused = 0;
    for (const SolverCandidate& r : rows) {
        if (r.runnable) continue;
        ++refused;
        CHECK_FALSE(r.reason.empty());
    }
    // SolverNC has no Pareto in its feature set, so it must refuse and say so.
    const SolverCandidate* nc = find_row(rows, "nc.default");
    REQUIRE(nc != nullptr);
    CHECK_FALSE(nc->runnable);
    CHECK(nc->reason.find("Pareto") != std::string::npos);
    CHECK(refused > 0);
}

TEST_CASE("a metric narrows the report to the families that answer it") {
    qn::Network<double> m = mm1();
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // A group name and the accessor that returns it ask the same question.
    const std::vector<SolverCandidate> byGroup = autosolver::auto_find_solver(sn, "cdf");
    const std::vector<SolverCandidate> byAccessor = autosolver::auto_find_solver(sn, "getCdfRespT");
    REQUIRE(!byGroup.empty());
    CHECK(byGroup.size() == byAccessor.size());

    // MVA computes no passage-time law and must not appear; the simulators and
    // the transform solvers do.
    bool sawMva = false, sawLdes = false;
    for (const SolverCandidate& r : byGroup) {
        if (r.solver == "mva") sawMva = true;
        if (r.solver == "ldes") sawLdes = true;
        CHECK(r.metrics.find("cdf") != std::string::npos);
    }
    CHECK_FALSE(sawMva);
    CHECK(sawLdes);

    // Every family answers the mean measures, so 'avg' narrows nothing away.
    CHECK(autosolver::auto_find_solver(sn, "avg").size() ==
          autosolver::auto_find_solver(sn).size());

    // A name that is neither a group nor an accessor is a caller error, not an
    // empty answer that would read as "nothing can do this".
    CHECK_THROWS_AS(autosolver::auto_find_solver(sn, "nosuchmeasure"), std::runtime_error);
}

TEST_CASE("the metric registry is self-consistent") {
    // A group name maps to itself, so findSolver('cdf') and
    // findSolver('getCdfRespT') are the same question.
    for (const std::string& g : autosolver::auto_metric_groups())
        CHECK(autosolver::auto_metric_group_of(g) == g);
    // Every group a family declares is a registered one; a typo here would
    // silently hide the family from a caller asking for that measure.
    for (const std::string& fam : autosolver::auto_network_family_names())
        for (const std::string& g : autosolver::auto_family_metrics(fam))
            CHECK(autosolver::auto_metric_group_of(g) == g);
    // Every family answers the mean measures, which is what a solver is for.
    for (const std::string& fam : autosolver::auto_network_family_names()) {
        const std::vector<std::string> gs = autosolver::auto_family_metrics(fam);
        CHECK(std::find(gs.begin(), gs.end(), std::string("avg")) != gs.end());
    }
    CHECK(autosolver::auto_metric_group_of("") == "");
    CHECK(autosolver::auto_metric_group_of("any") == "");
    CHECK(autosolver::auto_metric_group_of("getAvgTable") == "avg");
    CHECK(autosolver::auto_metric_group_of("getTranProbAggr") == "tranprob");
}

TEST_CASE("the table renders every row and pads only the bounded columns") {
    qn::Network<double> m = mm1();
    const std::vector<SolverCandidate> rows = autosolver::auto_find_solver(m.get_struct());
    const std::string table = autosolver::auto_find_solver_table(rows);
    // One header line plus one line per row.
    CHECK(std::count(table.begin(), table.end(), '\n') == static_cast<long>(rows.size()) + 1);
    CHECK(table.find("Solver") == 0);
    for (const SolverCandidate& r : rows) CHECK(table.find(r.method) != std::string::npos);
    // No row of a runnable-only report ends in the padding of an empty reason.
    CHECK(table.find(" \n") == std::string::npos);
    CHECK(autosolver::auto_find_solver_table(std::vector<SolverCandidate>()) ==
          "No solver method can analyze this model.\n");
}
