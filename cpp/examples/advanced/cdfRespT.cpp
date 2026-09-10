/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/advanced/cdfRespT/`: the response-time CDF, not the mean.
 *
 * Every entry reads its law off SolverFLD's passage-time analysis
 * (`fluid::solver_fluid_cdf_respt`, reached through `cdf_respt("FLD", ...)`),
 * which is the solver the reference uses. The simulation curves the scripts put
 * beside them -- `JMT.getCdfRespT` and `JMT.getTranCdfRespT` -- are refused by
 * name: JMT has no C++ port, and a fluid curve printed under a simulation
 * caption would be a silent solver substitution.
 *
 * The mean and SCV are recomputed from the curve exactly as the reference
 * does, sum_j (F_{j+1} - F_j) t_{j+1} and the same sum against t^2, so the two
 * codebases summarize the same object the same way.
 */

#include <cstdio>
#include <string>
#include <vector>

#include "example_util.h"
#include "examples_common.h"
#include "line/lang/dist_fitters.h"

namespace line {
namespace examples {

namespace {

/** The Riemann-Stieltjes moments the reference reads off a CDF curve. */
struct CdfStats {
    double mean = 0.0;
    double scv = 0.0;
    bool valid = false;
};

/**
 * Where in an interval the mass dF is charged, and the choice is the LAW's and
 * not a taste.
 *
 * A SIMULATED CDF jumps AT its samples, so RIGHT is not an approximation there
 * at all: the right endpoint IS the observed sample and the sum is exact. A
 * FLUID CDF is a CONTINUOUS law read off an integrator grid, where the mass of
 * an interval sits somewhere inside it -- RIGHT is then first-order and carries
 * a bias that moves with whatever grid the integrator chose, while MIDPOINT is
 * second-order. On cdf_respt_closed, whose exact answers are its own service
 * laws, MIDPOINT reads SCV 1.00032 and 0.33301 against 1 and 1/3 where RIGHT
 * reads 1.01715 and 0.33409.
 *
 * The other examples in this file still read their fluid curves with RIGHT:
 * their goldens were recorded under it and re-recording them is its own change.
 */
enum CdfQuad { CDF_QUAD_RIGHT = 0, CDF_QUAD_MIDPOINT = 1 };

CdfStats cdf_stats(const CdfCurve& c, CdfQuad quad = CDF_QUAD_RIGHT) {
    CdfStats s;
    if (c.t.size() < 2) return s;
    double m1 = 0.0, m2 = 0.0;
    for (std::size_t j = 0; j + 1 < c.t.size(); ++j) {
        const double dF = c.F[j + 1] - c.F[j];
        const double t =
            quad == CDF_QUAD_MIDPOINT ? 0.5 * (c.t[j] + c.t[j + 1]) : c.t[j + 1];
        m1 += dF * t;
        m2 += dF * t * t;
    }
    s.mean = m1;
    s.scv = m1 > 0.0 ? (m2 - m1 * m1) / (m1 * m1) : 0.0;
    s.valid = true;
    return s;
}

/** One row of the per-(station, class) CDF summary the examples print. */
void print_cdf_stats(const Sn& sn, const std::vector<std::vector<CdfCurve> >& rd,
                     const std::string& caption, CdfQuad quad = CDF_QUAD_RIGHT) {
    std::printf("%s\n", caption.c_str());
    std::printf("%-16s %-14s %12s %12s %10s\n", "Station", "JobClass", "MeanFromCdf", "ScvFromCdf",
                "Points");
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        // A Source has no response time of its own; the reference skips station 0.
        if (sn.stations[i].sched == SchedStrategy::EXT) continue;
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            const CdfStats s = cdf_stats(rd[i][c], quad);
            std::printf("%-16s %-14s %12.6g %12.6g %10zu\n", sn.stations[i].name.c_str(),
                        sn.classes[c].name.c_str(), s.mean, s.scv, rd[i][c].t.size());
        }
    }
}

/**
 * The labelled matrix the reference prints, and the ONE thing the parity parser
 * reads out of this example.
 *
 * `print_cdf_stats` above is the readable per-row summary; `parse_cdf_statistics`
 * in the shared parser keys on `Average Response Time from CDF (<solver>):` (or
 * the SCV line) followed by a nested list, which is what
 * `print(avg_respt_from_cdf_sim)` writes in the Python twin and what
 * `AvgRespTfromCDFSim =` writes in the MATLAB one. `(Simulation)` names JMT and
 * `(Fluid)` names FLD; an unlabelled line is read as FLD.
 *
 * EVERY station is emitted, a Source included, because the reference loops over
 * `getNumberOfStations()` and stores a zero where the law is empty: the parser
 * drops the all-zero rows and renumbers, so a twin that skipped them itself
 * would shift `Station2` onto `Station1`.
 */
void print_cdf_matrix(const std::string& label, const Sn& sn,
                      const std::vector<std::vector<CdfCurve> >& rd, bool scv,
                      const std::string& solver, CdfQuad quad) {
    std::vector<std::vector<double> > values(sn.nstations, std::vector<double>(sn.nclasses, 0.0));
    std::printf("\n%s:\n[", label.c_str());
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        if (i) std::printf(", ");
        std::printf("[");
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            if (c) std::printf(", ");
            const CdfStats s = cdf_stats(rd[i][c], quad);
            values[i][c] = s.valid ? (scv ? s.scv : s.mean) : 0.0;
            std::printf("%.10g", values[i][c]);
        }
        std::printf("]");
    }
    std::printf("]\n\n");
    // These statistics are a QUADRATURE THE EXAMPLE PERFORMS over the law the
    // solver returned, so no result table carries them and no getter can be
    // hooked to record them: the example declares the keys itself. See
    // `record_cdf_matrix` for why the row number is not the station index.
    record_cdf_matrix(solver, scv, values);
}

/** Both statistics of one law, under the reference's two labels. */
void print_cdf_stats_pair(const Sn& sn, const std::vector<std::vector<CdfCurve> >& rd,
                          const std::string& suffix, CdfQuad quad = CDF_QUAD_RIGHT) {
    // `(Simulation)` names JMT and `(Fluid)` -- or no suffix at all -- names FLD,
    // which is the same reading the shared parser made of these labels.
    const std::string solver = suffix.find("Simulation") != std::string::npos ? "JMT" : "FLD";
    print_cdf_matrix("Average Response Time from CDF" + suffix, sn, rd, false, solver, quad);
    print_cdf_matrix("Squared Coefficient of Variation from CDF" + suffix, sn, rd, true, solver,
                     quad);
}

/** `Network.circul(M)` as a closed single-chain ring over the given nodes. */
void link_ring(Net& m, const std::vector<std::size_t>& nodes) {
    Routing P;
    const std::size_t K = m.raw_struct().classes.size();
    for (std::size_t r = 1; r <= K; ++r)
        for (std::size_t j = 0; j < nodes.size(); ++j)
            P.set(r, r, nodes[j], nodes[(j + 1) % nodes.size()], 1.0);
    m.link(P);
}

/** The `cdf_respt_populations` model at one population. */
Net populations_model(double N) {
    Net m("model");
    Delay d(m, "Delay");
    Queue q1(m, "Queue1", SchedStrategy::PS);
    Queue q2(m, "Queue2", SchedStrategy::PS);
    ClosedClass k(m, "Class1", N, d);
    m.raw_struct().classes[k - 1].completes = false;
    d.set_service(k, Exp(1.0));
    q1.set_service(k, Exp(0.5));
    q2.set_service(k, Exp(0.5));
    link_ring(m, {d, q1, q2});
    return m;
}

}  // namespace

void cdf_respt_closed() {
    note("=== CDF Response Time - Closed Network ===");

    const D serv1 = Exp(1.0 / 0.1);
    const D serv2 = erlang_fit(1.0, 1.0 / 3.0);

    Net m("model");
    Delay d(m, "Delay");
    Queue q(m, "Queue2", SchedStrategy::PS);
    ClosedClass k(m, "Class1", 1, d);
    d.set_service(k, serv1);
    q.set_service(k, serv2);
    link_serial(m, {d, q});

    note("Computing CDF with JMT simulation...");
    section("JMT");
    const std::vector<std::vector<CdfCurve> > rd_sim =
        cdf_respt("JMT", m, sim_opts(23000, 10000));
    print_cdf_stats(m.get_struct(), rd_sim, "Response-time CDF (Simulation):");

    note("Computing CDF with Fluid solver...");
    section("FLD");
    const std::vector<std::vector<CdfCurve> > rd = cdf_respt("FLD", m);
    print_cdf_stats(m.get_struct(), rd, "Response-time CDF (Fluid):", CDF_QUAD_MIDPOINT);

    print_cdf_stats_pair(m.get_struct(), rd_sim, " (Simulation)");
    print_cdf_stats_pair(m.get_struct(), rd, " (Fluid)", CDF_QUAD_MIDPOINT);

    note("\nSince there is a single job, mean and squared coefficient of variation");
    note("of response times are close, up to fluid approximation precision, to those");
    note("of the service time distribution.");
    print_vector("Average Response Time from Theory", {serv1.mean, serv2.mean});
    print_vector("Squared Coefficient of Variation from Theory", {serv1.scv, serv2.scv});
}

LINE_EXAMPLE("advanced/cdfRespT", cdf_respt_closed);

void cdf_respt_closed_threeclasses() {
    note("=== CDF Response Time - Closed Network with Three Classes ===");

    Net m("model");
    Delay d(m, "Delay");
    Queue q(m, "Queue2", SchedStrategy::PS);
    ClosedClass k1(m, "Class1", 1, d);
    ClosedClass k2(m, "Class2", 0, d);
    ClosedClass k3(m, "Class3", 0, d);
    // Class 1 never completes, so the population is carried by the switch cycle.
    m.raw_struct().classes[k1 - 1].completes = false;
    d.set_service(k1, Exp(1.0));
    d.set_service(k2, Exp(1.0));
    d.set_service(k3, Exp(1.0));
    q.set_service(k1, Exp(1.0));
    q.set_service(k2, Erlang(0.5, 2));
    q.set_service(k3, Exp(1.0 / 0.01));

    Routing P;
    P.set(k1, k1, d, q, 1.0);
    P.set(k1, k2, q, d, 1.0);
    P.set(k2, k2, d, q, 1.0);
    P.set(k2, k1, q, d, 1.0);
    P.set(k3, k3, d, q, 1.0);
    P.set(k3, k3, q, d, 1.0);
    m.link(P);

    note("Computing with Fluid solver (state-dependent method)...");
    SolverOpts o;
    o.method = "statedep";
    o.iter_max = 100;

    section("FLD");
    print_vector("Average Response Time", solve_avg("FLD", m, o).column("RespT"));
    const std::vector<std::vector<CdfCurve> > rd = cdf_respt("FLD", m, o);
    print_cdf_stats(m.get_struct(), rd, "Response-time CDF (Fluid):");
    // The reference labels neither block here, and an unlabelled one is FLD.
    print_cdf_stats_pair(m.get_struct(), rd, "");

    note("\nNote: This example demonstrates class-switching behavior:");
    note("  - Class 1 jobs switch to Class 2 at the queue");
    note("  - Class 2 jobs switch back to Class 1 at the delay");
    note("  - Class 3 jobs follow a simple circular route without switching");
}

LINE_EXAMPLE("advanced/cdfRespT", cdf_respt_closed_threeclasses);

void cdf_respt_distrib() {
    note("=== CDF Response Time - Different Service Distributions ===");

    Net m("model");
    Delay d(m, "Delay");
    Queue q(m, "Queue1", SchedStrategy::PS);
    ClosedClass k1(m, "Class1", 1, d);
    d.set_service(k1, D::exp_mean(1.0));
    q.set_service(k1, D::exp_mean(2.0));
    ClosedClass k2(m, "Class2", 3, d);
    d.set_service(k2, erlang_fit_order(4.0, 2));
    // HyperExp.fitMeanAndSCV is map_hyperexp(mean, scv, p=0.99); example_util.h's
    // `hyperexp_fit` is the BALANCED-means fit, which is a different law.
    q.set_service(k2, lang::hyperexp_fit_mean_scv<double>(5.0, 30.0));
    link_ring(m, {d, q});

    note("Computing CDF with Fluid solver (steady-state)...");
    section("FLD");
    const std::vector<std::vector<CdfCurve> > rd = cdf_respt("FLD", m);
    print_cdf_stats(m.get_struct(), rd, "Response-time CDF (Fluid):");

    section("JMT");
    const std::vector<std::vector<CdfCurve> > rd_sim =
        cdf_respt("JMT", m, sim_opts(23000, 100000));
    print_cdf_stats(m.get_struct(), rd_sim, "Response-time CDF (Simulation):");

    print_cdf_stats_pair(m.get_struct(), rd, " (Fluid)");
    print_cdf_stats_pair(m.get_struct(), rd_sim, " (Simulation)");

    // The reference guards getServiceProcess with try/except and falls back to
    // 'N/A'; a pair the class never visits is Disabled here and reads the same.
    const Sn& sn = m.get_struct();
    note("\nService Processes by Station and Class:");
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            const bool has = sn.service[i][c].type != lang::ProcessType::DISABLED;
            kv("proc[" + sn.stations[i].name + "," + sn.classes[c].name + "]",
               has ? std::string(lang::process_to_text(sn.service[i][c].type))
                   : std::string("N/A"));
        }

    note("\nNote: Class 1 uses Exponential service distributions.");
    note("      Class 2 uses Erlang (delay) and HyperExponential (queue) distributions.");
    note("      This shows how different service distributions affect response time CDFs.");
}

LINE_EXAMPLE("advanced/cdfRespT", cdf_respt_distrib);

void cdf_respt_open_twoclasses() {
    note("=== CDF Response Time - Open Network with Two Classes ===");

    Net m("myModel");
    Source src(m, "Source");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    Sink snk(m, "Sink");
    OpenClass k1(m, "Class1", 0);
    OpenClass k2(m, "Class2", 0);
    src.set_arrival(k1, D::exp_mean(4.0));
    src.set_arrival(k2, D::exp_mean(4.0));
    q1.set_service(k1, D::exp_mean(1.0));
    q1.set_service(k2, D::exp_mean(1.0));
    q2.set_service(k1, D::exp_mean(1.0));
    q2.set_service(k2, D::exp_mean(1.0));

    // Each class switches into the other on every hop, as P{r,s}(i,j) declares.
    Routing P;
    P.set(k1, k1, src, q1, 1.0);
    P.set(k1, k2, q1, q2, 1.0);
    P.set(k2, k1, q2, snk, 1.0);
    P.set(k2, k2, src, q1, 1.0);
    P.set(k2, k1, q1, q2, 1.0);
    P.set(k1, k2, q2, snk, 1.0);
    m.link(P);

    note("Computing CDF with Fluid solver...");
    SolverOpts o;
    o.iter_max = 300;
    section("FLD");
    const std::vector<std::vector<CdfCurve> > rd = cdf_respt("FLD", m, o);
    print_cdf_stats(m.get_struct(), rd, "Response-time CDF (Fluid):");

    section("JMT");
    const std::vector<std::vector<CdfCurve> > rd_sim =
        cdf_respt("JMT", m, sim_opts(23000, 10000));
    print_cdf_stats(m.get_struct(), rd_sim, "Response-time CDF (Simulation):");

    print_cdf_stats_pair(m.get_struct(), rd, " (Fluid)");
    print_cdf_stats_pair(m.get_struct(), rd_sim, " (Simulation)");

    note("\nNote: Response time tail distributions are available in the curves above");
    note("      for detailed analysis and plotting.");
}

LINE_EXAMPLE("advanced/cdfRespT", cdf_respt_open_twoclasses);

void cdf_respt_populations() {
    const double n_jobs[3] = {1, 4, 8};
    note("=== CDF Response Time - Varying Population Sizes ===");

    for (int i = 0; i < 3; ++i) {
        std::printf("\n--- Population N = %g jobs ---\n", n_jobs[i]);
        Net m = populations_model(n_jobs[i]);

        SolverOpts o;
        o.iter_max = 100;
        section("FLD");
        const AvgTable t = solve_avg("FLD", m, o);
        // The reference prints the vector twice: once plain, once tagged with N as
        // its codebase-independent parity target.
        print_vector("Average Response Time", t.column("RespT"));
        // `FLD getAvgRespT N=<k>: v1 v2 ...`, character for character as the
        // reference prints it: this is the ONE line the parity parser reads for
        // this example (parser.parse_direct_respt), and `print_vector` writes
        // `label = ...`, which that pattern does not match.
        const std::vector<double> respt = t.column("RespT");
        std::printf("FLD getAvgRespT N=%d:", static_cast<int>(n_jobs[i]));
        for (std::size_t k = 0; k < respt.size(); ++k) std::printf(" %.6f", respt[k]);
        std::printf("\n");
        // THE TWO AGGREGATE ROWS ARE WHAT THE GOLDEN ASSERTS AT THE LARGER
        // POPULATIONS, and deliberately: once every queue saturates the fluid
        // equilibrium is a SET, and which member an engine's ODE trajectory
        // lands on is not fixed. The cycle time and the queueing total are
        // invariant along that set, so those are the invariants to compare.
        const std::string pop = "N" + std::to_string(static_cast<int>(n_jobs[i]));
        double total = 0.0, queues = 0.0;
        for (std::size_t k = 0; k < respt.size(); ++k) {
            derived("FLD", "Station" + std::to_string(k + 1), pop, respt[k], "RespT");
            total += respt[k];
            if (k > 0) queues += respt[k];
        }
        if (!respt.empty()) {
            derived("FLD", "AllStations", pop, total, "RespT");
            derived("FLD", "Queues", pop, queues, "RespT");
        }
        print_cdf_stats(m.get_struct(), cdf_respt("FLD", m, o), "Response-time CDF (Fluid):");
    }

    note("\n=== Summary ===");
    note("As population increases, response times increase due to queueing effects.");
    note("The CDF shifts to the right, showing higher probability of longer response times.");
}

LINE_EXAMPLE("advanced/cdfRespT", cdf_respt_populations);

}  // namespace examples
}  // namespace line
