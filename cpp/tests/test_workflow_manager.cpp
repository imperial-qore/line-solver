/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * workflow_manager: the facade over the analyzer and the AUTO integration.
 *
 * Nothing here computes a queueing result, so the oracle is the reference's
 * arithmetic and the reference's BYTES. The export layouts are asserted
 * character for character because a consumer parses them, the efficiency
 * scores are asserted against the closed forms they are, and the complexity
 * bands are driven across their three cut points.
 *
 * The one deliberate divergence from the Java is also pinned: its int
 * `factorial` overflows at 13 and goes NEGATIVE at 17, which would make the
 * complexity score DROP as a fork widens. The port computes in double, so the
 * test asserts monotonicity -- the property the Java loses.
 */
#include <cstddef>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/wf/workflow_manager.h"
#include "line/lang/qn/network_builder.h"

namespace wf = line::wf;
namespace qn = line::qn;
using line::Matrix;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

qn::Network<double> cycle() {
    qn::Network<double> m("C");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

std::vector<std::vector<int>> groups(const std::vector<std::size_t>& sizes) {
    std::vector<std::vector<int>> out;
    for (std::size_t i = 0; i < sizes.size(); ++i) {
        std::vector<int> g;
        for (std::size_t k = 0; k < sizes[i]; ++k) g.push_back(static_cast<int>(k));
        out.push_back(g);
    }
    return out;
}

}  // namespace

TEST_CASE("the efficiency scores are the reference's closed forms") {
    // Empty means perfect, by the reference's convention.
    CHECK(wf::calculate_sequence_efficiency(std::vector<std::vector<int>>()) ==
          doctest::Approx(1.0));
    CHECK(wf::calculate_parallel_efficiency(std::vector<std::vector<int>>()) ==
          doctest::Approx(1.0));
    CHECK(wf::calculate_loop_efficiency(std::vector<int>()) == doctest::Approx(1.0));

    // Sequences: avg length over 5, saturating at one.
    std::vector<std::size_t> s2;
    s2.push_back(2);
    s2.push_back(4);
    CHECK(wf::calculate_sequence_efficiency(groups(s2)) == doctest::Approx(3.0 / 5.0));
    std::vector<std::size_t> s9(1, 9);
    CHECK(wf::calculate_sequence_efficiency(groups(s9)) == doctest::Approx(1.0));

    // Parallels: 1 - (avg - 2)/10, floored at 0.1.
    std::vector<std::size_t> p4(1, 4);
    CHECK(wf::calculate_parallel_efficiency(groups(p4)) == doctest::Approx(0.8));
    std::vector<std::size_t> p20(1, 20);
    CHECK(wf::calculate_parallel_efficiency(groups(p20)) == doctest::Approx(0.1));

    // Loops: the reference's placeholder, a flat one half.
    std::vector<int> one(1, 3);
    CHECK(wf::calculate_loop_efficiency(one) == doctest::Approx(0.5));
    std::vector<int> many(7, 3);
    CHECK(wf::calculate_loop_efficiency(many) == doctest::Approx(0.5));
}

TEST_CASE("the branch efficiency is the mean normalized entropy") {
    std::vector<wf::BranchPattern<double>> b;
    wf::BranchPattern<double> even;
    even.branchNodes.push_back(1);
    even.branchNodes.push_back(2);
    even.probabilities.push_back(0.5);
    even.probabilities.push_back(0.5);
    b.push_back(even);
    // An even two-way split has maximal entropy, so the normalized value is one.
    CHECK(wf::calculate_branch_efficiency(b) == doctest::Approx(1.0).epsilon(1e-9));

    wf::BranchPattern<double> skewed;
    skewed.branchNodes.push_back(1);
    skewed.branchNodes.push_back(2);
    skewed.probabilities.push_back(0.99);
    skewed.probabilities.push_back(0.01);
    b.push_back(skewed);
    // Averaging in a near-deterministic branch must lower it.
    CHECK(wf::calculate_branch_efficiency(b) < 1.0);
    CHECK(wf::calculate_branch_efficiency(b) > 0.5);

    CHECK(wf::calculate_branch_efficiency(std::vector<wf::BranchPattern<double>>()) ==
          doctest::Approx(1.0));
}

TEST_CASE("the complexity score stays monotone where the Java overflows") {
    // The Java's int factorial goes NEGATIVE at 17, so its score DROPS as the
    // fork widens. This is the property that divergence buys.
    double prev = -1.0;
    for (std::size_t width = 2; width <= 20; ++width) {
        wf::DetectedPatterns<double> p;
        std::vector<std::size_t> one(1, width);
        p.parallels = groups(one);
        const double score = wf::calculate_complexity_score<double>(10, 12, p);
        CHECK(score > prev);
        prev = score;
    }

    // Below the overflow the two agree exactly: 5! * 0.2 = 24 on top of the
    // 10 + 12*0.5 size term.
    wf::DetectedPatterns<double> p5;
    std::vector<std::size_t> five(1, 5);
    p5.parallels = groups(five);
    CHECK(wf::calculate_complexity_score<double>(10, 12, p5) ==
          doctest::Approx(10.0 + 6.0 + 120.0 * 0.2));
}

TEST_CASE("the complexity level follows the reference's three cut points") {
    // The score is nodes + 0.5*links with no patterns, so the level is driven
    // directly through each band.
    wf::DetectedPatterns<double> none;
    CHECK(wf::calculate_complexity_score<double>(4, 4, none) == doctest::Approx(6.0));
    CHECK(wf::calculate_complexity_score<double>(40, 4, none) == doctest::Approx(42.0));
    CHECK(wf::calculate_complexity_score<double>(80, 20, none) == doctest::Approx(90.0));
    CHECK(wf::calculate_complexity_score<double>(200, 20, none) == doctest::Approx(210.0));

    // A single loop is worth ten points on its own.
    wf::DetectedPatterns<double> oneLoop;
    oneLoop.loops.push_back(3);
    CHECK(wf::calculate_complexity_score<double>(0, 0, oneLoop) == doctest::Approx(10.0));
}

TEST_CASE("the complexity report bands a real model and tallies its patterns") {
    qn::Network<double> m = cycle();
    const wf::ComplexityReport<double> r = wf::generate_complexity_report(m.get_struct());

    CHECK(r.sequences.count == 0);
    CHECK(r.overallComplexityScore >= 0.0);
    // Whatever band it lands in, the label must agree with the score.
    if (r.overallComplexityScore < 10.0) CHECK(r.complexityLevel == "LOW");
    else if (r.overallComplexityScore < 50.0) CHECK(r.complexityLevel == "MEDIUM");
    else if (r.overallComplexityScore < 100.0) CHECK(r.complexityLevel == "HIGH");
    else CHECK(r.complexityLevel == "VERY_HIGH");

    // The two complexity blocks describe the same graph before and after the
    // collapse, so the collapse never adds nodes.
    CHECK(r.optimizedMetrics.totalNodes <= r.originalMetrics.totalNodes);
}

TEST_CASE("the summary export is byte-for-byte the reference's layout") {
    wf::WorkflowAnalysisResult<double> a;
    a.solverRecommendation.recommendedSolver = wf::WfSolver::MVA;
    a.solverRecommendation.confidence = 0.85;
    a.solverRecommendation.reasoning.push_back("First reason");
    a.solverRecommendation.reasoning.push_back("Second reason");

    const std::string s = wf::export_analysis(a, wf::WfExportFormat::Summary);
    const std::string expected =
        "=== Workflow Analysis Summary ===\n"
        "\n"
        "Recommended Solver: MVA\n"
        "Confidence: 0.85\n"
        "\n"
        "Detected Patterns:\n"
        "- Sequences: 0\n"
        "- Parallels: 0\n"
        "- Loops: 0\n"
        "- Branches: 0\n"
        "\n"
        "Reasoning:\n"
        "1. First reason\n"
        "2. Second reason\n";
    CHECK(s == expected);
}

TEST_CASE("the JSON and CSV exports are the reference's layouts") {
    wf::WorkflowAnalysisResult<double> a;
    a.solverRecommendation.recommendedSolver = wf::WfSolver::SSA;
    a.solverRecommendation.confidence = 0.5;

    const std::string j = wf::export_analysis(a, wf::WfExportFormat::Json);
    const std::string expectedJson =
        "{\n"
        "  \"solver_recommendation\": \"SSA\",\n"
        "  \"confidence\": 0.5,\n"
        "  \"patterns\": {\n"
        "    \"sequences\": 0,\n"
        "    \"parallels\": 0,\n"
        "    \"loops\": 0,\n"
        "    \"branches\": 0\n"
        "  }\n"
        "}";
    CHECK(j == expectedJson);

    const std::string c = wf::export_analysis(a, wf::WfExportFormat::Csv);
    const std::string expectedCsv =
        "Metric,Value\n"
        "Recommended Solver,SSA\n"
        "Confidence,0.5\n"
        "Sequences,0\n"
        "Parallels,0\n"
        "Loops,0\n"
        "Branches,0\n";
    CHECK(c == expectedCsv);
}

TEST_CASE("the benchmark records a success, an aggregate and a failure") {
    std::vector<wf::WfSolver> which = wf::default_benchmark_solvers();
    CHECK(which.size() == 4);

    // NC throws; the rest return a 2x1 queue-length vector summing to 3.
    std::function<Matrix<double>(wf::WfSolver)> runner = [](wf::WfSolver s) {
        if (s == wf::WfSolver::NC) throw line::InputError("no product form");
        Matrix<double> QN(2, 1, 0.0);
        QN(0, 0) = 1.0;
        QN(1, 0) = 2.0;
        return QN;
    };
    const std::vector<std::pair<wf::WfSolver, wf::BenchmarkRow>> rows =
        wf::benchmark_solvers<double>(which, runner);

    CHECK(rows.size() == 4);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        if (rows[i].first == wf::WfSolver::NC) {
            CHECK(rows[i].second.success == false);
            CHECK(rows[i].second.error.find("product form") != std::string::npos);
        } else {
            CHECK(rows[i].second.success == true);
            CHECK(rows[i].second.hasResults == true);
            CHECK(rows[i].second.totalQueueLength == doctest::Approx(3.0));
            CHECK(rows[i].second.avgQueueLength == doctest::Approx(1.5));
            CHECK(rows[i].second.solveTimeMs >= 0.0);
        }
    }
}

TEST_CASE("the facade validates a real model and reports its size") {
    qn::Network<double> m = cycle();
    const qn::NetworkStruct<double> sn = m.get_struct();

    const wf::WorkflowValidation v = wf::validate_workflow(sn);
    CHECK(v.nodeCount == sn.nodes.size());
    CHECK(v.hasRouting == true);

    // THE REFERENCE CALLS A PATTERN-FREE MODEL INVALID, and this pins it. Its
    // validateWorkflowEnhancement requires a NON-EMPTY reasoning list, but every
    // reasoning line is produced by a pattern arm, so a plain two-station cycle
    // -- which is a perfectly valid workflow -- fails the check. The pattern
    // analysis itself passes; only the enhancement test objects.
    CHECK(v.isValid == false);
    CHECK(v.issues.size() == 1);
    CHECK(v.issues[0] == "Workflow-enhanced solver selection validation failed");
    CHECK(wf::validate_analysis(wf::analyze_workflow(wf::wf_from_struct(sn))) == true);

    const wf::WorkflowAnalysisResult<double> a = wf::analyze_workflow_full(sn);
    CHECK(a.performanceMetrics.solverConfidence == a.solverRecommendation.confidence);
    CHECK(a.performanceMetrics.sequenceEfficiency ==
          doctest::Approx(wf::calculate_sequence_efficiency(a.patternAnalysis.detectedPatterns.sequences)));

    // quick_analysis is the summary of exactly that analysis.
    CHECK(wf::quick_analysis(sn) == wf::export_analysis(a, wf::WfExportFormat::Summary));
    CHECK(wf::get_optimal_solver(sn) == a.solverRecommendation.recommendedSolver);
}
