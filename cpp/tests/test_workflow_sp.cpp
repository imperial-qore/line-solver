/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Workflow: the series-parallel composition and the geometric loop.
 *
 * The oracle is CLOSED FORM throughout, not a golden: a serial pair of
 * exponentials is Erlang-2, a parallel pair is the maximum of two
 * exponentials, and a geometric compound of mean COUNT has mean COUNT*m and
 * SCV cx/COUNT + 1 - 1/COUNT. Those identities also pin the twin tests in
 * jar/src/test/java/jline/lang/workflow/WorkflowSPTest.java and
 * python/tests/test_workflow_sp.py, so the three ports are checked against the
 * same numbers rather than against each other.
 *
 * The case the BLOCK heuristic cannot reduce -- an AND-fork nested inside a
 * loop -- is asserted here as well, because it is the whole reason the
 * series-parallel path exists.
 */
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/io/workflow_reader.h"
#include "line/lang/workflow/workflow.h"

namespace wfl = line::workflow;
using line::Matrix;
using Dist = line::lang::Distrib<double>;
using Wf = wfl::Workflow<double>;
using Law = wfl::PhLaw<double>;

namespace {

constexpr double kTol = 1e-9;

/** Mean and SCV of the law (alpha, S), by -alpha S^-1 e and 2 alpha S^-2 e. */
void moments(const Law& law, double& mean, double& scv) {
    const std::size_t n = law.S.rows();
    // Solve S x = -e for x = -S^-1 e, then S y = -x for y = S^-2 e
    Matrix<double> A(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) A(i, j) = law.S(i, j);
    std::vector<double> x(n, -1.0), y(n, 0.0);

    // Gaussian elimination with partial pivoting on [A | -e], reused for -x
    std::vector<std::vector<double>> M(n, std::vector<double>(n + 2, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) M[i][j] = A(i, j);
        M[i][n] = -1.0;
    }
    for (std::size_t c = 0; c < n; ++c) {
        std::size_t piv = c;
        for (std::size_t r = c + 1; r < n; ++r)
            if (std::abs(M[r][c]) > std::abs(M[piv][c])) piv = r;
        std::swap(M[c], M[piv]);
        for (std::size_t r = 0; r < n; ++r) {
            if (r == c) continue;
            const double f = M[r][c] / M[c][c];
            for (std::size_t j = c; j <= n; ++j) M[r][j] -= f * M[c][j];
        }
    }
    for (std::size_t i = 0; i < n; ++i) x[i] = M[i][n] / M[i][i];

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) M[i][j] = A(i, j);
        M[i][n] = -x[i];
    }
    for (std::size_t c = 0; c < n; ++c) {
        std::size_t piv = c;
        for (std::size_t r = c + 1; r < n; ++r)
            if (std::abs(M[r][c]) > std::abs(M[piv][c])) piv = r;
        std::swap(M[c], M[piv]);
        for (std::size_t r = 0; r < n; ++r) {
            if (r == c) continue;
            const double f = M[r][c] / M[c][c];
            for (std::size_t j = c; j <= n; ++j) M[r][j] -= f * M[c][j];
        }
    }
    for (std::size_t i = 0; i < n; ++i) y[i] = M[i][n] / M[i][i];

    double m1 = 0.0, m2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        m1 += law.alpha[i] * x[i];
        m2 += 2.0 * law.alpha[i] * y[i];
    }
    mean = m1;
    scv = (m2 - m1 * m1) / (m1 * m1);
}

void moments(const Dist& d, double& mean, double& scv) {
    Law law;
    law.S = d.D0;
    law.alpha.assign(d.params.begin(), d.params.begin() + static_cast<long>(d.D0.rows()));
    moments(law, mean, scv);
}

Law exp_law(double mean) {
    Law l;
    l.alpha.assign(1, 1.0);
    l.S = Matrix<double>(1, 1, -1.0 / mean);
    return l;
}

/** E[max(X, Y)] for independent exponentials of the given means. */
double emax(double mx, double my) {
    const double rx = 1.0 / mx, ry = 1.0 / my;
    return 1 / rx + 1 / ry - 1 / (rx + ry);
}

Wf fork_in_loop(double mean_c = 2.0) {
    Wf wf("ForkInLoop");
    wf.add_activity("A", Dist::exp_mean(0.5));
    wf.add_activity("B", Dist::exp_mean(1.0));
    wf.add_activity("C", Dist::exp_mean(mean_c));
    wf.add_activity("D", Dist::exp_mean(1.5));
    wf.add_activity("E", Dist::exp_mean(0.25));
    wf.add_activity("F", Dist::exp_mean(0.75));
    wf.add_precedence(Wf::Loop("A", {"B", "F"}, 2.0));
    wf.add_precedence(Wf::AndFork("B", {"C", "D"}));
    wf.add_precedence(Wf::AndJoin({"C", "D"}, "E"));
    return wf;
}

}  // namespace

TEST_CASE("workflow: a serial pair of exponentials is Erlang-2") {
    Wf wf("serial2");
    wf.add_activity("A", Dist::exp_mean(1.0));
    wf.add_activity("B", Dist::exp_mean(1.0));
    wf.add_precedence(Wf::Serial("A", "B"));

    double m = 0.0, c = 0.0;
    moments(wf.to_ph(), m, c);
    CHECK(m == doctest::Approx(2.0).epsilon(kTol));
    CHECK(c == doctest::Approx(0.5).epsilon(kTol));
    CHECK(wf.sp_tree() != nullptr);
}

TEST_CASE("workflow: a parallel pair is the maximum of the branches") {
    const Law par = Wf::compose_parallel(exp_law(2.0), exp_law(1.5));
    double m = 0.0, c = 0.0;
    moments(par, m, c);
    CHECK(m == doctest::Approx(emax(2.0, 1.5)).epsilon(kTol));
}

TEST_CASE("workflow: a geometric loop of an exponential is exponential") {
    const Law law = Wf::compose_loop_geometric(exp_law(2.0), 3.0);
    double m = 0.0, c = 0.0;
    moments(law, m, c);
    CHECK(m == doctest::Approx(6.0).epsilon(kTol));
    CHECK(c == doctest::Approx(1.0).epsilon(kTol));
}

TEST_CASE("workflow: a geometric loop matches the compound moments") {
    // Erlang-2 body of mean 2 and SCV 1/2
    Law body;
    body.alpha = {1.0, 0.0};
    body.S = Matrix<double>(2, 2, 0.0);
    body.S(0, 0) = -1.0;
    body.S(0, 1) = 1.0;
    body.S(1, 1) = -1.0;

    const double count = 3.0, scv_body = 0.5;
    const Law law = Wf::compose_loop_geometric(body, count);
    double m = 0.0, c = 0.0;
    moments(law, m, c);

    CHECK(m == doctest::Approx(count * 2.0).epsilon(kTol));
    CHECK(c == doctest::Approx(scv_body / count + 1 - 1 / count).epsilon(kTol));
    // the order is that of the body, and the generator is cyclic
    CHECK(law.S.rows() == 2u);
    CHECK_FALSE(Wf::is_acyclic_generator(law.S));
}

TEST_CASE("workflow: a fractional loop count runs the body with that probability") {
    const Law law = Wf::compose_loop_geometric(exp_law(2.0), 0.25);
    double m = 0.0, c = 0.0;
    moments(law, m, c);
    CHECK(m == doctest::Approx(0.5).epsilon(1e-6));
}

TEST_CASE("workflow: a loop keeps the mean and the body order") {
    Wf wf("LoopWorkflow");
    wf.add_activity("A", Dist::exp_mean(1.0));
    wf.add_activity("B", Dist::exp_mean(2.0));
    wf.add_activity("C", Dist::exp_mean(0.5));
    wf.add_precedence(Wf::Loop("A", {"B", "C"}, 3.0));

    const Dist ph = wf.to_ph();
    double m = 0.0, c = 0.0;
    moments(ph, m, c);
    CHECK(m == doctest::Approx(7.5).epsilon(kTol));
    // A, the geometric loop over B, and C
    CHECK(ph.D0.rows() == 3u);
    CHECK(c == doctest::Approx((1.0 + 36.0 + 0.25) / (7.5 * 7.5)).epsilon(kTol));
}

TEST_CASE("workflow: an OR-fork mixes the branches") {
    Wf wf("BranchingWorkflow");
    wf.add_activity("A", Dist::exp_mean(1.0));
    wf.add_activity("B", Dist::exp_mean(2.0));
    wf.add_activity("C", Dist::exp_mean(5.0));
    wf.add_activity("D", Dist::exp_mean(0.5));
    wf.add_precedence(Wf::OrFork("A", {"B", "C"}, {0.6, 0.4}));
    wf.add_precedence(Wf::OrJoin({"B", "C"}, "D"));

    double m = 0.0, c = 0.0;
    moments(wf.to_ph(), m, c);
    CHECK(m == doctest::Approx(1.0 + 0.6 * 2.0 + 0.4 * 5.0 + 0.5).epsilon(kTol));
    CHECK(wf.sp_tree() != nullptr);
}

TEST_CASE("workflow: a fork nested inside a loop is reduced exactly") {
    Wf wf = fork_in_loop();
    double m = 0.0, c = 0.0;
    moments(wf.to_ph(), m, c);
    const double body = 1.0 + emax(2.0, 1.5) + 0.25;
    CHECK(m == doctest::Approx(0.5 + 2 * body + 0.75).epsilon(kTol));
    CHECK(wf.sp_tree() != nullptr);
}

TEST_CASE("workflow: an incremental refresh matches a rebuilt workflow") {
    Wf wf = fork_in_loop();
    wf.to_ph();

    Wf ref = fork_in_loop(3.0);
    const Dist ph_ref = ref.to_ph();
    double m_ref = 0.0, c_ref = 0.0;
    moments(ph_ref, m_ref, c_ref);

    wf.set_activity_demand("C", Dist::exp_mean(3.0));
    const Dist ph_inc = wf.refresh_ph();
    double m_inc = 0.0, c_inc = 0.0;
    moments(ph_inc, m_inc, c_inc);

    CHECK(m_inc == doctest::Approx(m_ref).epsilon(kTol));
    CHECK(c_inc == doctest::Approx(c_ref).epsilon(kTol));
    CHECK(ph_inc.D0.rows() == ph_ref.D0.rows());
}

TEST_CASE("workflow: a mean-only rescale keeps the shape") {
    Wf wf("Rescale");
    wf.add_activity("A", line::lang::aph_fit_mean_scv(2.0, 0.3));
    wf.add_activity("B", Dist::exp_mean(1.0));
    wf.add_precedence(Wf::Serial("A", "B"));
    wf.to_ph();

    wf.set_activity_demand_mean("A", 5.0);
    const Dist scaled = wf.refresh_ph();
    double m = 0.0, c = 0.0;
    moments(scaled, m, c);

    Wf ref("RescaleRef");
    ref.add_activity("A", line::lang::aph_fit_mean_scv(5.0, 0.3));
    ref.add_activity("B", Dist::exp_mean(1.0));
    ref.add_precedence(Wf::Serial("A", "B"));
    const Dist ph_ref = ref.to_ph();
    double m_ref = 0.0, c_ref = 0.0;
    moments(ph_ref, m_ref, c_ref);

    CHECK(m == doctest::Approx(6.0).epsilon(kTol));
    CHECK(c == doctest::Approx(c_ref).epsilon(1e-6));
    CHECK(scaled.D0.rows() == ph_ref.D0.rows());
}

TEST_CASE("workflow: a quorum join is refused by name") {
    Wf wf("Quorum");
    wf.add_activity("A", Dist::exp_mean(1.0));
    wf.add_activity("B", Dist::exp_mean(1.0));
    wf.add_activity("C", Dist::exp_mean(1.0));
    wf.add_activity("D", Dist::exp_mean(1.0));
    wf.add_precedence(Wf::AndFork("A", {"B", "C"}));
    wf.add_precedence(Wf::AndJoin({"B", "C"}, "D", {1.0}));

    CHECK_THROWS_AS(wf.to_ph(), line::UnsupportedError);
}

TEST_CASE("workflow: a full AND-join is accepted") {
    Wf wf("FullJoin");
    wf.add_activity("A", Dist::exp_mean(1.0));
    wf.add_activity("B", Dist::exp_mean(1.0));
    wf.add_activity("C", Dist::exp_mean(1.0));
    wf.add_activity("D", Dist::exp_mean(1.0));
    wf.add_precedence(Wf::AndFork("A", {"B", "C"}));
    wf.add_precedence(Wf::AndJoin({"B", "C"}, "D"));

    double m = 0.0, c = 0.0;
    moments(wf.to_ph(), m, c);
    CHECK(m == doctest::Approx(1.0 + emax(1.0, 1.0) + 1.0).epsilon(kTol));
}

TEST_CASE("workflow: a graph that is not series-parallel falls back") {
    Wf wf("NotSP");
    wf.add_activity("A", Dist::exp_mean(1.0));
    wf.add_activity("B", Dist::exp_mean(1.0));
    wf.add_activity("C", Dist::exp_mean(1.0));
    wf.add_precedence(Wf::Serial("A", "B"));
    wf.add_precedence(Wf::Serial("A", "C"));

    double m = 0.0, c = 0.0;
    moments(wf.to_ph(), m, c);
    CHECK(m > 0.0);
    CHECK(wf.sp_tree() == nullptr);
}

TEST_CASE("workflow: execution counts weight the loop body") {
    Wf wf("Execs");
    wf.add_activity("A", Dist::exp_mean(1.0));
    wf.add_activity("B", Dist::exp_mean(1.0));
    wf.add_activity("C", Dist::exp_mean(1.0));
    wf.add_precedence(Wf::Loop("A", {"B", "C"}, 3.0));

    const wfl::SPTree<double>* tree = wf.sp_tree();
    REQUIRE(tree != nullptr);
    CHECK(tree->execs[tree->leaf_of[0]] == doctest::Approx(1.0));  // A
    CHECK(tree->execs[tree->leaf_of[1]] == doctest::Approx(3.0));  // B, the body
    CHECK(tree->execs[tree->leaf_of[2]] == doctest::Approx(1.0));  // C
}

TEST_CASE("workflow: the wf_* example suite reproduces its reference means") {
    // wf_serial: A -> B -> C
    {
        Wf wf("SerialWorkflow");
        wf.add_activity("A", Dist::exp_mean(1.0));
        wf.add_activity("B", Dist::exp_mean(2.0));
        wf.add_activity("C", Dist::exp_mean(1.5));
        wf.add_precedence(Wf::Serial("A", "B"));
        wf.add_precedence(Wf::Serial("B", "C"));
        double m = 0.0, c = 0.0;
        moments(wf.to_ph(), m, c);
        CHECK(m == doctest::Approx(4.5).epsilon(kTol));
    }
    // wf_complex: A -> B -> [C || D] -> F -> G
    {
        Wf wf("ComplexWorkflow");
        wf.add_activity("A", Dist::exp_mean(0.5));
        wf.add_activity("B", Dist::exp_mean(1.0));
        wf.add_activity("C", Dist::exp_mean(2.0));
        wf.add_activity("D", Dist::exp_mean(1.5));
        wf.add_activity("F", Dist::exp_mean(1.0));
        wf.add_activity("G", Dist::exp_mean(0.5));
        wf.add_precedence(Wf::Serial("A", "B"));
        wf.add_precedence(Wf::AndFork("B", {"C", "D"}));
        wf.add_precedence(Wf::AndJoin({"C", "D"}, "F"));
        wf.add_precedence(Wf::Serial("F", "G"));
        const Dist ph = wf.to_ph();
        double m = 0.0, c = 0.0;
        moments(ph, m, c);
        CHECK(m == doctest::Approx(0.5 + 1.0 + emax(2.0, 1.5) + 1.0 + 0.5).epsilon(kTol));
        CHECK(ph.D0.rows() == 7u);
    }
}


/*
 * The model.json wire path. The document below is what the PYTHON writer
 * emitted for `python/examples/basic/workflowModels/wf_loop.py`, copied
 * verbatim rather than hand-authored: the point of the test is that C++ reads
 * what the other codebases WRITE, and a hand-written approximation of the
 * schema would pass while the real document still failed.
 *
 * It also pins the loop arm across the wire, which is where the semantics
 * changed: `postType: "post-LOOP"` with `postParams: [3]` must compose the
 * GEOMETRIC law (mean 7.5, three phases), not the 3-fold convolution (mean
 * 7.5, five phases). The mean alone cannot tell those apart, which is exactly
 * how the original defect survived.
 */
TEST_CASE("workflow: a model.json written by the Python writer is read back") {
    const std::string doc = R"JSON({
  "format": "line-model",
  "version": "1.0",
  "model": {
    "type": "Workflow",
    "name": "LoopWorkflow",
    "activities": [
      {"name": "A", "hostDemand": {"type": "Exp", "params": {"lambda": 1.0}}},
      {"name": "B", "hostDemand": {"type": "Exp", "params": {"lambda": 0.5}}},
      {"name": "C", "hostDemand": {"type": "Exp", "params": {"lambda": 2.0}}}
    ],
    "precedences": [
      {
        "preActs": ["A"],
        "postActs": ["B", "C"],
        "preType": "pre",
        "postType": "post-LOOP",
        "postParams": [3]
      }
    ]
  }
})JSON";

    const line::io::detail::json root = line::io::detail::json::parse(doc);
    Wf wf = line::io::build_workflow_from_json<double>(root);

    CHECK(wf.name() == "LoopWorkflow");
    const Dist ph = wf.to_ph();
    double m = 0.0, c = 0.0;
    moments(ph, m, c);
    CHECK(m == doctest::Approx(7.5).epsilon(kTol));
    CHECK(ph.D0.rows() == 3u);
    // SCV of the geometric compound: cx/k + 1 - 1/k on the composed law.
    CHECK(c == doctest::Approx(0.66222222).epsilon(1e-6));
}

/*
 * The refusals. Each names the thing it will not do, because a reader that
 * DROPS one of these returns a confident law for a different workflow.
 */
TEST_CASE("workflow: the reader refuses what it does not implement") {
    // A Network envelope is not a Workflow, and says which reader takes it.
    {
        const std::string doc = R"JSON({"model": {"type": "Network", "name": "n"}})JSON";
        const line::io::detail::json root = line::io::detail::json::parse(doc);
        CHECK_THROWS_AS(line::io::build_workflow_from_json<double>(root), line::UnsupportedError);
    }
    // An unknown model-level key is refused rather than dropped.
    {
        const std::string doc =
            R"JSON({"model": {"type": "Workflow", "name": "w", "jobClasses": []}})JSON";
        const line::io::detail::json root = line::io::detail::json::parse(doc);
        CHECK_THROWS_AS(line::io::build_workflow_from_json<double>(root), line::UnsupportedError);
    }
    // An unknown precedence spelling is refused by name.
    {
        const std::string doc = R"JSON({"model": {"type": "Workflow", "name": "w",
            "activities": [{"name": "A", "hostDemand": {"type": "Exp", "params": {"lambda": 1.0}}},
                           {"name": "B", "hostDemand": {"type": "Exp", "params": {"lambda": 1.0}}}],
            "precedences": [{"preActs": ["A"], "postActs": ["B"],
                             "preType": "pre", "postType": "post-FORK"}]}})JSON";
        const line::io::detail::json root = line::io::detail::json::parse(doc);
        CHECK_THROWS_AS(line::io::build_workflow_from_json<double>(root), line::UnsupportedError);
    }
    // A precedence naming an activity that was never declared is an input error
    // raised by validate(), not silently a no-op edge.
    {
        const std::string doc = R"JSON({"model": {"type": "Workflow", "name": "w",
            "activities": [{"name": "A", "hostDemand": {"type": "Exp", "params": {"lambda": 1.0}}}],
            "precedences": [{"preActs": ["A"], "postActs": ["ghost"],
                             "preType": "pre", "postType": "post"}]}})JSON";
        const line::io::detail::json root = line::io::detail::json::parse(doc);
        Wf wf = line::io::build_workflow_from_json<double>(root);
        CHECK_THROWS(wf.to_ph());
    }
}
