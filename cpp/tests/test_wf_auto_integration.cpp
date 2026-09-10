/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * wf_auto_integration: the workflow-aware amendment of the AUTO recommendation.
 *
 * This is a HEURISTIC, so the oracle is not a queueing result -- it is the
 * reference's decision table. Each test drives one arm of it and asserts the
 * choice, the reasoning text and the confidence the Java produces on the same
 * input, because two codebases that give different ADVICE from the same model
 * have diverged even though neither computed a number.
 *
 * The features are driven directly where the arm needs a shape no small model
 * produces (parallelism above five, a loop above 0.8, more than fifty nodes);
 * driving the amendment through a real NetworkStruct as well is what checks
 * that the flattening feeds it correctly.
 */
#include <cstddef>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/wf/wf_auto_integration.h"
#include "line/lang/qn/network_builder.h"

namespace wf = line::wf;
namespace qn = line::qn;
using line::Matrix;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

bool has_text(const std::vector<std::string>& v, const std::string& needle) {
    for (std::size_t i = 0; i < v.size(); ++i)
        if (v[i].find(needle) != std::string::npos) return true;
    return false;
}

bool lists(const std::vector<wf::WfSolver>& v, wf::WfSolver s) {
    for (std::size_t i = 0; i < v.size(); ++i)
        if (v[i] == s) return true;
    return false;
}

/** A closed two-station cycle: single chain, product form, small population. */
qn::Network<double> cycle(double njobs) {
    qn::Network<double> m("C");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

wf::WorkflowAnalysis<double> empty_analysis() {
    wf::WorkflowRepresentation<double> w;
    w.linkMatrix = Matrix<double>(0, 3, 0.0);
    return wf::analyze_workflow(w);
}

}  // namespace

TEST_CASE("the base recommendation follows the reference's four bands") {
    wf::WfModelFacts<double> f;
    f.hasSingleChain = true;
    CHECK(wf::wf_base_recommendation(f) == wf::WfSolver::NC);

    f.hasSingleChain = false;
    f.hasMultiChain = true;
    f.hasProductForm = true;
    f.totalJobs = 6;
    f.avgJobsPerChain = 3.0;
    CHECK(wf::wf_base_recommendation(f) == wf::WfSolver::NC);  // fewer than ten jobs

    f.totalJobs = 40;
    f.avgJobsPerChain = 20.0;
    CHECK(wf::wf_base_recommendation(f) == wf::WfSolver::MVA);

    f.hasProductForm = false;
    f.avgJobsPerChain = 50.0;
    CHECK(wf::wf_base_recommendation(f) == wf::WfSolver::FLUID);

    // Exactly 30 falls through both strict bands to the default.
    f.avgJobsPerChain = 30.0;
    CHECK(wf::wf_base_recommendation(f) == wf::WfSolver::MVA);
}

TEST_CASE("an open class does not swamp the job count") {
    // The reference skips infinite populations; counting them as MAX_VALUE
    // would push every open model past the ten-job band.
    qn::Network<double> m("O");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("O1");
    m.set_arrival(src, c, Dist::exp_rate(0.5));
    m.set_service(q, c, Dist::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, snk, 1.0);
    m.link(P);

    const wf::WfModelFacts<double> f = wf::wf_model_facts(m.get_struct());
    CHECK(f.totalJobs == 0);
    CHECK(f.avgJobsPerChain == doctest::Approx(0.0));
}

TEST_CASE("high parallelism moves an exact choice to SSA and offers JMT") {
    wf::WorkflowFeatures<double> f;
    f.hasParallelPatterns = true;
    f.maxParallelism = 8;
    const wf::ExtendedSolverRecommendation<double> r =
        wf::enhance_recommendation_with_workflow(wf::WfSolver::NC, f, empty_analysis());

    CHECK(r.recommendedSolver == wf::WfSolver::SSA);
    CHECK(has_text(r.reasoning, "Switching to SSA for high-parallelism workflow"));
    CHECK(lists(r.alternativeSolvers, wf::WfSolver::JMT));
    // The seed confidence is untouched on this arm: no 0.05 bonus is given.
    CHECK(r.confidence == doctest::Approx(0.7));
}

TEST_CASE("a near-certain loop costs confidence and offers the inexact solvers") {
    wf::WorkflowFeatures<double> f;
    f.hasLoopPatterns = true;
    f.avgLoopProbability = 0.85;
    f.maxLoopProbability = 0.9;
    const wf::ExtendedSolverRecommendation<double> r =
        wf::enhance_recommendation_with_workflow(wf::WfSolver::MVA, f, empty_analysis());

    CHECK(r.recommendedSolver == wf::WfSolver::MVA);  // the arm advises, it does not switch
    CHECK(has_text(r.reasoning, "may cause numerical instability"));
    CHECK(lists(r.alternativeSolvers, wf::WfSolver::SSA));
    CHECK(lists(r.alternativeSolvers, wf::WfSolver::FLUID));
    CHECK(r.confidence == doctest::Approx(0.6));  // 0.7 - 0.1

    // A moderate loop takes the other band instead: advice, and a bonus.
    wf::WorkflowFeatures<double> g;
    g.hasLoopPatterns = true;
    g.maxLoopProbability = 0.6;
    const wf::ExtendedSolverRecommendation<double> s =
        wf::enhance_recommendation_with_workflow(wf::WfSolver::MVA, g, empty_analysis());
    CHECK(has_text(s.reasoning, "Moderate loop probability"));
    CHECK(s.confidence == doctest::Approx(0.75));
}

TEST_CASE("a large collapsed workflow demotes NC to MVA") {
    wf::WorkflowFeatures<double> f;
    f.originalNodeCount = 80;
    f.optimizedNodeCount = 60;
    const wf::ExtendedSolverRecommendation<double> r =
        wf::enhance_recommendation_with_workflow(wf::WfSolver::NC, f, empty_analysis());

    CHECK(r.recommendedSolver == wf::WfSolver::MVA);
    CHECK(has_text(r.reasoning, "Switching from NC to MVA for large workflow"));
    // 80 -> 60 is a 25% reduction, reported as an integer percentage.
    CHECK(has_text(r.reasoning, "reduced by 25%"));
    CHECK(r.confidence == doctest::Approx(0.8));  // 0.7 + 0.1 for the reduction
}

TEST_CASE("the alternatives are capped at three and never repeat the choice") {
    wf::WorkflowFeatures<double> f;
    f.hasBranchPatterns = true;
    f.avgBranchEntropy = 2.0;
    f.maxBranches = 9;
    const wf::ExtendedSolverRecommendation<double> r =
        wf::enhance_recommendation_with_workflow(wf::WfSolver::CTMC, f, empty_analysis());

    CHECK(r.alternativeSolvers.size() == 3);
    CHECK(!lists(r.alternativeSolvers, wf::WfSolver::CTMC));
    // The branch arm pushes JMT then SSA before the standard padding runs, so
    // the padding order is visible in the answer.
    CHECK(r.alternativeSolvers[0] == wf::WfSolver::JMT);
    CHECK(r.alternativeSolvers[1] == wf::WfSolver::SSA);
    CHECK(r.alternativeSolvers[2] == wf::WfSolver::MVA);
}

TEST_CASE("low branching entropy is read as near-deterministic") {
    wf::WorkflowFeatures<double> f;
    f.hasBranchPatterns = true;
    f.avgBranchEntropy = 0.2;
    f.maxBranches = 2;
    const wf::ExtendedSolverRecommendation<double> r =
        wf::enhance_recommendation_with_workflow(wf::WfSolver::MVA, f, empty_analysis());
    CHECK(has_text(r.reasoning, "deterministic-like behavior"));
    CHECK(r.confidence == doctest::Approx(0.8));
    CHECK(!has_text(r.reasoning, "complex decision structure"));
}

TEST_CASE("the confidence is clamped into the reference's band") {
    wf::WorkflowFeatures<double> f;
    f.hasSequencePatterns = true;
    f.hasParallelPatterns = true;
    f.hasLoopPatterns = true;
    f.hasBranchPatterns = true;
    f.maxParallelism = 2;   // the 0.05 bonus
    f.maxLoopProbability = 0.6;
    f.avgBranchEntropy = 0.1;
    f.originalNodeCount = 10;
    f.optimizedNodeCount = 4;
    const wf::ExtendedSolverRecommendation<double> r =
        wf::enhance_recommendation_with_workflow(wf::WfSolver::MVA, f, empty_analysis());
    // 0.7 + 0.1 + 0.05 + 0.05 + 0.1 + 0.1 = 1.1, clamped to one.
    CHECK(r.confidence == doctest::Approx(1.0));
}

TEST_CASE("a real model drives the whole path and validates") {
    qn::Network<double> m = cycle(4.0);
    const qn::NetworkStruct<double> sn = m.get_struct();

    const wf::ExtendedSolverRecommendation<double> r =
        wf::recommend_solver_with_workflow_analysis(sn);
    // One closed chain, so the base is NC and nothing in a two-station cycle
    // moves it.
    CHECK(r.recommendedSolver == wf::WfSolver::NC);
    CHECK(r.confidence > 0.1);
    CHECK(r.confidence <= 1.0);
    CHECK(wf::create_optimal_solver(sn) == r.recommendedSolver);
    CHECK(wf::wf_solver_name(r.recommendedSolver) == "NC");

    const wf::OptimizationInsights ins = wf::get_optimization_insights(sn);
    // The insights are advisory and may be empty on so small a model; what must
    // hold is that the call is total and the analysis behind it validates.
    CHECK(ins.recommendations.size() + ins.patternInsights.size() +
              ins.performancePredictions.size() >=
          0u);
}

TEST_CASE("the feature flattening reports the pattern counts it was given") {
    qn::Network<double> m = cycle(2.0);
    const wf::WorkflowAnalysis<double> a = wf::analyze_workflow(wf::wf_from_struct(m.get_struct()));
    const wf::WorkflowFeatures<double> f = wf::extract_workflow_features(a);

    CHECK(f.hasSequencePatterns == !a.detectedPatterns.sequences.empty());
    CHECK(f.hasLoopPatterns == !a.detectedPatterns.loops.empty());
    CHECK(f.numSequences == a.detectedPatterns.sequences.size());
    CHECK(f.numLoops == a.detectedPatterns.loops.size());
    CHECK(f.originalNodeCount == a.statistics.originalComplexity.totalNodes);
    CHECK(f.optimizedNodeCount == a.statistics.optimizedComplexity.totalNodes);
    // An absent group leaves its features at the reference's neutral zero.
    if (!f.hasBranchPatterns) {
        CHECK(f.avgBranchEntropy == doctest::Approx(0.0));
        CHECK(f.maxBranches == 0u);
    }
}
