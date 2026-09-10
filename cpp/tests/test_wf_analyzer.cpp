/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * wf_analyzer: the composition of the four detectors, the collapse and the two
 * complexity reports. Ported from native Python on 2026-08-01.
 *
 * The detectors themselves are covered elsewhere; what is checked here is the
 * COMPOSITION -- that the collapse runs on what was detected, that the two
 * complexity blocks describe the before and after of the same graph, and that
 * the NetworkStruct conversion classifies the nodes and normalizes the routing
 * the way the reference's own object-model walk does.
 */
#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/wf/wf_analyzer.h"
#include "line/lang/qn/network_builder.h"

namespace wf = line::wf;
namespace qn = line::qn;
using line::Matrix;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

wf::ServiceParameters<double> expo(double lambda) {
    wf::ServiceParameters<double> p;
    p.alpha.assign(1, 1.0);
    p.T_ = Matrix<double>(1, 1, -lambda);
    return p;
}

/** 1 -> 2 -> 3 -> 4 -> 5 with 2,3,4 the service nodes. */
wf::WorkflowRepresentation<double> chain() {
    wf::WorkflowRepresentation<double> w;
    w.linkMatrix = Matrix<double>(4, 3, 0.0);
    for (std::size_t i = 0; i < 4; ++i) {
        w.linkMatrix(i, 0) = static_cast<double>(i + 1);
        w.linkMatrix(i, 1) = static_cast<double>(i + 2);
        w.linkMatrix(i, 2) = 1.0;
    }
    w.serviceNodes.push_back(2);
    w.serviceNodes.push_back(3);
    w.serviceNodes.push_back(4);
    w.serviceParameters[2] = expo(2.0);
    w.serviceParameters[3] = expo(4.0);
    w.serviceParameters[4] = expo(5.0);
    return w;
}

}  // namespace

TEST_CASE("analyze_workflow detects the chain and collapses it") {
    const wf::WorkflowRepresentation<double> w = chain();
    const wf::WorkflowAnalysis<double> a = wf::analyze_workflow(w);

    CHECK(a.detectedPatterns.sequences.size() >= 1);
    CHECK(a.detectedPatterns.parallels.empty());
    CHECK(a.detectedPatterns.loops.empty());

    // The collapse ran on what was detected: fewer links, fewer laws.
    CHECK(a.optimizedWorkflow.linkMatrix.rows() < w.linkMatrix.rows());
    CHECK(a.optimizedWorkflow.serviceParameters.size() < w.serviceParameters.size());
    CHECK(a.statistics.updateStats.originalLinks == 4);
    CHECK(a.statistics.updateStats.linksReduced > 0);
}

TEST_CASE("the two complexity blocks describe the before and the after") {
    const wf::WorkflowRepresentation<double> w = chain();
    const wf::WorkflowAnalysis<double> a = wf::analyze_workflow(w);

    const wf::WorkflowComplexity<double>& o = a.statistics.originalComplexity;
    CHECK(o.serviceNodes == 3);
    CHECK(o.controlNodes == 0);
    CHECK(o.totalNodes == 3);
    CHECK(o.totalLinks == 4);
    CHECK(o.connectedNodes == 5);  // ids 1..5 appear in the matrix
    // 4 edges contribute 8 endpoints over 5 nodes.
    CHECK(o.avgDegree == doctest::Approx(8.0 / 5.0).epsilon(1e-12));

    const wf::WorkflowComplexity<double>& p = a.statistics.optimizedComplexity;
    CHECK(p.totalNodes == a.optimizedWorkflow.serviceParameters.size());
    CHECK(p.totalLinks == a.optimizedWorkflow.linkMatrix.rows());
    CHECK(p.connectedNodes <= o.connectedNodes);
}

TEST_CASE("an empty workflow analyzes to zeros rather than throwing") {
    wf::WorkflowRepresentation<double> w;
    w.linkMatrix = Matrix<double>(0, 3, 0.0);
    const wf::WorkflowAnalysis<double> a = wf::analyze_workflow(w);
    CHECK(a.detectedPatterns.sequences.empty());
    CHECK(a.statistics.originalComplexity.totalNodes == 0);
    CHECK(a.statistics.originalComplexity.avgDegree == doctest::Approx(0.0));
    CHECK(wf::get_optimization_recommendations(a).empty());
    CHECK(wf::validate_analysis(a) == true);
}

TEST_CASE("the recommendations name each detected family") {
    const wf::WorkflowAnalysis<double> a = wf::analyze_workflow(chain());
    const std::vector<std::string> rec = wf::get_optimization_recommendations(a);
    REQUIRE(!rec.empty());
    bool sawSequence = false;
    for (std::size_t i = 0; i < rec.size(); ++i)
        if (rec[i].find("sequence patterns that can be simplified") != std::string::npos)
            sawSequence = true;
    CHECK(sawSequence);
}

TEST_CASE("wf_from_struct classifies the nodes and normalizes the routing") {
    qn::Network<double> m("wfsn");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t rt = m.add_router("R");
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, rt, 1.0);
    P.set(rt, d, 1.0);
    m.link(P);

    const wf::WorkflowRepresentation<double> w = wf::wf_from_struct(m.get_struct());
    CHECK(w.serviceNodes.size() == 2);   // the Delay and the Queue
    CHECK(w.routerNodes.size() == 1);    // the Router
    CHECK(w.forkNodes.empty());
    CHECK(w.joinNodes.empty());
    CHECK(w.serviceParameters.size() == 2);

    // Every row is a probability and every source row sums to one.
    REQUIRE(w.linkMatrix.rows() > 0);
    std::map<int, double> rowsum;
    for (std::size_t i = 0; i < w.linkMatrix.rows(); ++i) {
        CHECK(w.linkMatrix(i, 2) > 0.0);
        CHECK(w.linkMatrix(i, 2) <= 1.0 + 1e-12);
        rowsum[static_cast<int>(w.linkMatrix(i, 0))] += w.linkMatrix(i, 2);
    }
    for (std::map<int, double>::const_iterator it = rowsum.begin(); it != rowsum.end(); ++it)
        CHECK(it->second == doctest::Approx(1.0).epsilon(1e-12));

    // The declared service laws reach the representation, not a placeholder.
    CHECK(w.serviceParameters.at(static_cast<int>(q) - 1).T_(0, 0) == doctest::Approx(-2.0));
}
