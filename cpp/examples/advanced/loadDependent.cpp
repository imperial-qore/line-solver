/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/advanced/loadDependent/`: load-, class- and joint-dependent
 * service rates, and the flow-equivalent-server aggregation built on them.
 *
 * THE FOUR FES ENTRIES RUN THEIR AGGREGATION BLOCK through
 * `line::fes::fes_aggregate`, the port of `ModelAdapter.aggregateFES`: stochastic
 * complement of the routing, isolated subnetwork, throughput table over the whole
 * population lattice, then a rebuilt network carrying that table as a
 * class-dependence handle. They refused it while only the pieces underneath
 * (`fes_build_isolated`, `fes_compute_throughputs`, `fes_beta_handle`) were
 * ported; the transformation itself has been here since, so the refusal outlived
 * its gap and these entries printed nothing where the reference prints a solve.
 */

#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

#include "example_util.h"
#include "examples_common.h"
#include "line/api/fes/fes_aggregate.h"
#include "line/api/sn/sn_gd_balance.h"

namespace line {
namespace examples {

namespace {

/** `[min(i, c) for i in range(1, n + 1)]`, the multi-server lld vector. */
std::vector<double> min_ramp(std::size_t n, double c) {
    std::vector<double> alpha;
    for (std::size_t i = 1; i <= n; ++i) alpha.push_back(std::min(static_cast<double>(i), c));
    return alpha;
}

/**
 * The aggregated model, with `subset` given as 1-BASED station indices in the
 * reference's own base. `fes_aggregate` returns the rebuilt network; the caller
 * solves it exactly as it solves the original, which is the point of Norton's
 * theorem: the replacement is another closed network, not a special object.
 */
Net aggregate_fes(Net& m, const std::vector<std::size_t>& subset) {
    return fes::fes_aggregate(m.get_struct(), subset).model;
}

/** The 4-station tandem the two `fes_*` scripts aggregate, at the given demands. */
Net fes_tandem(const std::string& nm, const std::string& delay_name, double n1, double n2,
               const std::vector<double>& delay_means, const std::vector<double>& q1_means,
               const std::vector<double>& q2_means, const std::vector<double>& q3_means,
               const std::vector<std::string>& queue_names) {
    Net m(nm);
    Delay d(m, delay_name);
    Queue q1(m, queue_names[0], SchedStrategy::PS);
    Queue q2(m, queue_names[1], SchedStrategy::PS);
    Queue q3(m, queue_names[2], SchedStrategy::PS);
    ClosedClass k1(m, "Class1", n1, d);
    d.set_service(k1, D::exp_mean(delay_means[0]));
    q1.set_service(k1, D::exp_mean(q1_means[0]));
    q2.set_service(k1, D::exp_mean(q2_means[0]));
    q3.set_service(k1, D::exp_mean(q3_means[0]));
    if (delay_means.size() > 1) {
        ClosedClass k2(m, "Class2", n2, d);
        d.set_service(k2, D::exp_mean(delay_means[1]));
        q1.set_service(k2, D::exp_mean(q1_means[1]));
        q2.set_service(k2, D::exp_mean(q2_means[1]));
        q3.set_service(k2, D::exp_mean(q3_means[1]));
    }
    link_serial(m, {d, q1, q2, q3});
    return m;
}

}  // namespace

// ---------------------------------------------------------------------------
// Flow-equivalent server aggregation
// ---------------------------------------------------------------------------

void fes_aggregation() {
    note("=== Flow-Equivalent Server (FES) Aggregation Example ===");
    note("Creating original 4-station network...");
    Net m = fes_tandem("OriginalModel", "ThinkTime", 3, 2, {5.0, 4.0}, {1.5, 2.0}, {1.0, 1.2},
                       {0.8, 1.0}, {"Queue1", "Queue2", "Queue3"});

    note("\n--- Solving Original Model ---");
    section("MVA");
    const AvgTable orig = solve_avg("MVA", m);
    print_avg(orig, "Original model results:");

    note("\n--- Creating FES Model ---");
    note("Aggregating Queue1 and Queue2 into a single FES...");
    Net fes = aggregate_fes(m, {2, 3});  // Queue1, Queue2
    const AvgTable agg = solve_avg("MVA", fes);
    print_avg(agg, "FES model results:");

    // Norton's theorem is exact for this network, so the two throughputs agree
    // and the relative error is the arithmetic's, not the transform's.
    note("\n--- Throughput Comparison ---");
    kv("Class1 Original Tput", orig.get("Tput", "ThinkTime", "Class1"));
    kv("Class1 FES Tput", agg.get("Tput", "ThinkTime", "Class1"));
    kv("Class2 Original Tput", orig.get("Tput", "ThinkTime", "Class2"));
    kv("Class2 FES Tput", agg.get("Tput", "ThinkTime", "Class2"));
}

LINE_EXAMPLE("advanced/loadDependent", fes_aggregation);

void fes_single_class() {
    note("=== Single-Class FES Aggregation Example ===");
    note("Creating original 4-station network...");
    Net m = fes_tandem("OriginalModel", "ThinkTime", 5, 0, {5.0}, {1.5}, {1.0}, {0.8},
                       {"Queue1", "Queue2", "Queue3"});

    note("\n--- Solving Original Model ---");
    section("MVA");
    const AvgTable orig = solve_avg("MVA", m);
    print_avg(orig, "Original model results:");

    note("\n--- Creating FES Model ---");
    note("Aggregating Queue1 and Queue2 into a single FES...");
    Net fes = aggregate_fes(m, {2, 3});  // Queue1, Queue2
    const AvgTable agg = solve_avg("MVA", fes);
    print_avg(agg, "FES model results:");

    note("\n--- Throughput Comparison ---");
    kv("Original Tput", orig.get("Tput", "ThinkTime", "Class1"));
    kv("FES Tput", agg.get("Tput", "ThinkTime", "Class1"));
}

LINE_EXAMPLE("advanced/loadDependent", fes_single_class);

void ld_fes_singleclass() {
    Net m = fes_tandem("OriginalModel", "Delay", 5, 0, {5.0}, {1.5}, {1.0}, {0.8},
                       {"Q1", "Q2", "Q3"});
    // NO `section` HERE. This block prints TWO solvers -- MVA on the original
    // model and NC on the FES aggregation -- and declaring the whole block
    // `MVA` filed the NC table under MVA's key, which is what the golden's
    // separate `NC` row went missing behind. Each `AvgTable` names its own
    // solver, and `print_avg` declares it.
    SolverOpts o;
    o.method = "exact";
    print_avg(solve_avg("MVA", m, o), "MVA (original):");
    Net fes = aggregate_fes(m, {2, 3, 4});  // Q1, Q2, Q3
    SolverOpts nco;
    nco.method = "exact";
    print_avg(solve_avg("NC", fes, nco), "NC (FES model):");
}

LINE_EXAMPLE("advanced/loadDependent", ld_fes_singleclass);

void ld_fes_multiclass() {
    Net m = fes_tandem("OriginalModel", "Delay", 3, 2, {1.0, 1.5}, {0.5, 0.8}, {0.3, 0.6},
                       {0.4, 0.7}, {"Q1", "Q2", "Q3"});
    // NO `section` HERE. This block prints TWO solvers -- MVA on the original
    // model and NC on the FES aggregation -- and declaring the whole block
    // `MVA` filed the NC table under MVA's key, which is what the golden's
    // separate `NC` row went missing behind. Each `AvgTable` names its own
    // solver, and `print_avg` declares it.
    SolverOpts o;
    o.method = "exact";
    print_avg(solve_avg("MVA", m, o), "MVA (original):");
    Net fes = aggregate_fes(m, {2, 3, 4});  // Q1, Q2, Q3
    SolverOpts nco;
    nco.method = "exact";
    print_avg(solve_avg("NC", fes, nco), "NC (FES model):");
}

LINE_EXAMPLE("advanced/loadDependent", ld_fes_multiclass);

// ---------------------------------------------------------------------------
// Class and joint dependence
// ---------------------------------------------------------------------------

void ld_class_dependence() {
    const double N = 16.0, c = 2.0;
    note("=== Load-Dependent Service - Class Dependence (product-form) ===");

    Net m("model");
    Delay d(m, "Delay");
    Queue q(m, "Queue1", SchedStrategy::PS);
    ClosedClass k1(m, "Class1", N, d);
    ClosedClass k2(m, "Class2", N / 2.0, d);
    d.set_service(k1, D::exp_mean(1.0));
    d.set_service(k2, D::exp_mean(2.0));
    q.set_service(k1, D::exp_mean(1.5));
    q.set_service(k2, D::exp_mean(2.5));
    // beta_{i,r}(n_{i,r}): class 1 scales up to c servers with its OWN count.
    m.set_class_dependence(q,
                           [c](const std::vector<double>& ni) {
                               return std::vector<double>{std::min(ni[0], c), 1.0};
                           },
                           std::vector<double>{c, 1.0});
    link_serial(m, {d, q});

    section("CTMC");
    print_avg(solve_avg("CTMC", m), "CTMC (exact):");

    section("MVA");
    SolverOpts o;
    o.method = "qd";
    print_avg(solve_avg("MVA", m, o), "MVA with QD method:");

    // The reference runs no JMT block here: the JSIM writer has no representation
    // for the class-dependence handle, so SolverJMT rejects the model.
    note("\nJMT is not solved: SolverJMT.getFeatureSet rejects the class-dependence handle.");

    note("\nNote: Class-dependent (product-form) service scales each class by its");
    note("      own per-class population, preserving the BCMP product form.");
}

LINE_EXAMPLE("advanced/loadDependent", ld_class_dependence);

void ld_joint_dependence() {
    const double N = 16.0, c = 2.0;
    note("=== Load-Dependent Service - Joint Dependence ===");

    Net m("model");
    Delay d(m, "Delay");
    Queue q(m, "Queue1", SchedStrategy::PS);
    ClosedClass k1(m, "Class1", N, d);
    ClosedClass k2(m, "Class2", N / 2.0, d);
    d.set_service(k1, D::exp_mean(1.0));
    d.set_service(k2, D::exp_mean(2.0));
    q.set_service(k1, D::exp_mean(1.5));
    q.set_service(k2, D::exp_mean(2.5));
    // Non-product-form: the rate reads the class-1 marginal and is shared.
    m.set_joint_dependence(q,
                           [c](const std::vector<double>& ni) {
                               return std::vector<double>{std::min(ni[0], c)};
                           },
                           std::vector<double>{c});
    link_serial(m, {d, q});

    section("CTMC");
    print_avg(solve_avg("CTMC", m), "CTMC (exact):");

    section("MVA");
    SolverOpts o;
    o.method = "qd";
    print_avg(solve_avg("MVA", m, o), "MVA with QD method:");

    // The reference runs no JMT block here: the JSIM writer has no representation
    // for the joint-dependence handle, so SolverJMT rejects the model.
    note("\nJMT is not solved: SolverJMT.getFeatureSet rejects the joint-dependence handle.");

    note("\nNote: Joint-dependent service is non-product-form (see setJointDependence).");
    note("      In this example, service rate scales with Class 1 population only,");
    note("      modeling c servers available exclusively for Class 1 jobs.");
    note("      Class 2 jobs do not benefit from multi-server parallelism.");
}

LINE_EXAMPLE("advanced/loadDependent", ld_joint_dependence);

// ---------------------------------------------------------------------------
// Multi-server queues expressed through load dependence
// ---------------------------------------------------------------------------

void ld_multiserver_fcfs() {
    const double N = 16.0, c = 2.0;
    note("=== Load-Dependent Service - Multi-Server FCFS Queue ===");

    note("--- Standard Multi-Server Model ---");
    Net ms("model");
    {
        Delay d(ms, "Delay");
        Queue q(ms, "Queue1", SchedStrategy::FCFS);
        ClosedClass k(ms, "Class1", N, d);
        ms.set_service(d, k, D::exp_mean(1.0));
        ms.set_service(q, k, D::exp_mean(1.5));
        ms.set_number_of_servers(q, c);
        link_serial(ms, {d, q});
    }
    section("NC");
    print_avg(solve_avg("NC", ms), "NC Solver:");

    note("\n--- Load-Dependent Model ---");
    Net ld("model");
    {
        Delay d(ld, "Delay");
        Queue q(ld, "Queue1", SchedStrategy::FCFS);
        ClosedClass k(ld, "Class1", N, d);
        ld.set_service(d, k, D::exp_mean(1.0));
        ld.set_service(q, k, D::exp_mean(1.5));
        ld.set_load_dependence(q, min_ramp(static_cast<std::size_t>(N), c));
        link_serial(ld, {d, q});
    }

    section("CTMC");
    print_avg(solve_avg("CTMC", ld), "CTMC (exact):");

    const char* nc_methods[] = {"default", "rd", "nrp", "nrl"};
    const char* nc_captions[] = {"NC (exact):", "NC with RD method:", "NC with NRP method:",
                                 "NC with NRL method:"};
    for (int i = 0; i < 4; ++i) {
        section("NC");
        SolverOpts o;
        o.method = nc_methods[i];
        print_avg(solve_avg("NC", ld, o), nc_captions[i]);
    }

    const char* mva_methods[] = {"exact", "qd"};
    const char* mva_captions[] = {"MVA (exact):", "MVA with QD method:"};
    for (int i = 0; i < 2; ++i) {
        section("MVA");
        SolverOpts o;
        o.method = mva_methods[i];
        print_avg(solve_avg("MVA", ld, o), mva_captions[i]);
    }

    section("JMT");
    print_avg(ld.get_struct(), jmt_avg(ld, sim_opts(23000)));

    note("\n--- Class-Dependent Model ---");
    Net cd("model");
    {
        Delay d(cd, "Delay");
        Queue q(cd, "Queue1", SchedStrategy::FCFS);
        ClosedClass k(cd, "Class1", N, d);
        cd.set_service(d, k, D::exp_mean(1.0));
        cd.set_service(q, k, D::exp_mean(1.5));
        cd.set_class_dependence(q,
                                [c](const std::vector<double>& ni) {
                                    double tot = 0.0;
                                    for (std::size_t r = 0; r < ni.size(); ++r) tot += ni[r];
                                    return std::vector<double>{std::min(tot, c)};
                                },
                                std::vector<double>{c});
        link_serial(cd, {d, q});
    }

    section("CTMC");
    print_avg(solve_avg("CTMC", cd), "CTMC (exact):");

    section("MVA");
    SolverOpts qd;
    qd.method = "qd";
    print_avg(solve_avg("MVA", cd, qd), "MVA with QD method:");

    // The reference runs no JMT block on the class-dependent model: the JSIM writer
    // has no representation for the handle, so SolverJMT rejects it. Only the lld
    // model above is JMT-solvable, because lldscaling maps onto the server count.
    note("\nJMT is not solved on the class-dependent model: the handle has no JSIM form.");
}

LINE_EXAMPLE("advanced/loadDependent", ld_multiserver_fcfs);

void ld_multiserver_ps() {
    const double N = 4.0, c = 3.0;
    note("=== Load-Dependent Service - Multi-Server PS Queue ===");

    note("--- Standard Multi-Server Model ---");
    Net ms("model");
    {
        Delay d(ms, "Delay");
        Queue q1(ms, "Queue1", SchedStrategy::PS);
        Queue q2(ms, "Queue2", SchedStrategy::PS);
        ClosedClass k1(ms, "Class1", N, d);
        ClosedClass k2(ms, "Class2", N / 2.0, d);
        ms.set_service(d, k1, D::exp_mean(1.0));
        ms.set_service(d, k2, D::exp_mean(2.0));
        ms.set_service(q1, k1, D::exp_mean(1.5));
        ms.set_service(q1, k2, D::exp_mean(2.5));
        ms.set_number_of_servers(q1, c);
        ms.set_service(q2, k1, D::exp_mean(3.5));
        ms.set_service(q2, k2, D::exp_mean(4.5));
        ms.set_number_of_servers(q2, c);
        link_serial(ms, {d, q1, q2});
    }
    section("MVA");
    SolverOpts exact;
    exact.method = "exact";
    print_avg(solve_avg("MVA", ms, exact), "MVA (exact):");

    note("\n--- Load-Dependent Model ---");
    Net ld("ldmodel");
    {
        Delay d(ld, "Delay");
        Queue q1(ld, "Queue1", SchedStrategy::PS);
        Queue q2(ld, "Queue2", SchedStrategy::PS);
        ClosedClass k1(ld, "Class1", N, d);
        ClosedClass k2(ld, "Class2", N / 2.0, d);
        ld.set_service(d, k1, D::exp_mean(1.0));
        ld.set_service(d, k2, D::exp_mean(2.0));
        ld.set_service(q1, k1, D::exp_mean(1.5));
        ld.set_service(q1, k2, D::exp_mean(2.5));
        ld.set_load_dependence(q1, min_ramp(static_cast<std::size_t>(N + N / 2.0), c));
        ld.set_service(q2, k1, D::exp_mean(3.5));
        ld.set_service(q2, k2, D::exp_mean(4.5));
        ld.set_load_dependence(q2, min_ramp(static_cast<std::size_t>(N + N / 2.0), c));
        link_serial(ld, {d, q1, q2});
    }

    section("CTMC");
    print_avg(solve_avg("CTMC", ld), "CTMC (exact):");

    const char* nc_methods[] = {"exact", "rd", "nrp", "nrl"};
    const char* nc_captions[] = {"NC (exact):", "NC with RD method:", "NC with NRP method:",
                                 "NC with NRL method:"};
    for (int i = 0; i < 4; ++i) {
        section("NC");
        SolverOpts o;
        o.method = nc_methods[i];
        print_avg(solve_avg("NC", ld, o), nc_captions[i]);
    }

    const char* mva_methods[] = {"exact", "qd"};
    const char* mva_captions[] = {"MVA (exact):", "MVA with QD method:"};
    for (int i = 0; i < 2; ++i) {
        section("MVA");
        SolverOpts o;
        o.method = mva_methods[i];
        print_avg(solve_avg("MVA", ld, o), mva_captions[i]);
    }

    section("JMT");
    print_avg(ld.get_struct(), jmt_avg(ld, sim_opts(23000)));

    note("\nNote: Load dependence allows modeling multi-server queues by scaling");
    note("      service rates based on the number of jobs present.");
}

LINE_EXAMPLE("advanced/loadDependent", ld_multiserver_ps);

void ld_multiserver_ps_twoclasses() {
    const double N = 4.0, c = 2.0;
    note("=== Load-Dependent Service - Multi-Server PS Queue (Two Classes) ===");

    note("--- Standard Multi-Server Model ---");
    Net ms("model");
    {
        Delay d(ms, "Delay");
        Queue q(ms, "Queue1", SchedStrategy::PS);
        ClosedClass k1(ms, "Class1", N, d);
        ClosedClass k2(ms, "Class2", N / 2.0, d);
        ms.set_service(d, k1, D::exp_mean(1.0));
        ms.set_service(d, k2, D::exp_mean(2.0));
        ms.set_service(q, k1, D::exp_mean(1.5));
        ms.set_service(q, k2, D::exp_mean(2.5));
        ms.set_number_of_servers(q, c);
        link_serial(ms, {d, q});
    }
    section("MVA");
    SolverOpts exact;
    exact.method = "exact";
    print_avg(solve_avg("MVA", ms, exact), "MVA (exact):");

    note("\n--- Load-Dependent Model ---");
    Net ld("ldmodel");
    {
        Delay d(ld, "Delay");
        Queue q(ld, "Queue1", SchedStrategy::PS);
        ClosedClass k1(ld, "Class1", N, d);
        ClosedClass k2(ld, "Class2", N / 2.0, d);
        ld.set_service(d, k1, D::exp_mean(1.0));
        ld.set_service(d, k2, D::exp_mean(2.0));
        ld.set_service(q, k1, D::exp_mean(1.5));
        ld.set_service(q, k2, D::exp_mean(2.5));
        ld.set_load_dependence(q, min_ramp(static_cast<std::size_t>(N + N / 2.0), c));
        link_serial(ld, {d, q});
    }

    section("CTMC");
    print_avg(solve_avg("CTMC", ld), "CTMC (exact):");

    const char* nc_methods[] = {"default", "rd", "nrp", "nrl"};
    const char* nc_captions[] = {"NC (exact):", "NC with RD method:", "NC with NRP method:",
                                 "NC with NRL method:"};
    for (int i = 0; i < 4; ++i) {
        section("NC");
        SolverOpts o;
        o.method = nc_methods[i];
        print_avg(solve_avg("NC", ld, o), nc_captions[i]);
    }

    const char* mva_methods[] = {"exact", "qd"};
    const char* mva_captions[] = {"MVA (exact):", "MVA with QD method:"};
    for (int i = 0; i < 2; ++i) {
        section("MVA");
        SolverOpts o;
        o.method = mva_methods[i];
        print_avg(solve_avg("MVA", ld, o), mva_captions[i]);
    }

    section("JMT");
    print_avg(ld.get_struct(), jmt_avg(ld, sim_opts(23000, 5000)));

    note("\nNote: Load dependence models multi-server behavior across multiple classes.");
    note("      Service rate scales with total queue length up to c servers.");
}

LINE_EXAMPLE("advanced/loadDependent", ld_multiserver_ps_twoclasses);

void ld_global_dependence() {
    const double N = 3.0;
    note("=== Load-Dependent Service - Global (Whittle) Dependence ===");

    // set_global_dependence declares a rate scaling phi(n) over the FULL
    // (nstations x nclasses) population matrix, not the population local to one
    // station. Here two PS stations share one unit of capacity,
    // phi_s(n) = n_s/|n|, the single-link allocation every alpha-fair rule
    // collapses to. It satisfies the Whittle balance property, so the chain is
    // reversible, has the product form pi(n) ~ Phi(n) prod rho^n, and is
    // INSENSITIVE to the service distribution beyond its mean.
    Net m("model");
    Queue q1(m, "Queue1", SchedStrategy::PS);
    Queue q2(m, "Queue2", SchedStrategy::PS);
    ClosedClass k1(m, "Class1", N, q1);
    q1.set_service(k1, D::exp_mean(1.0));
    q2.set_service(k1, D::exp_mean(0.5));
    link_serial(m, {q1, q2});

    const std::size_t M = m.get_struct().nstations;
    m.set_global_dependence(
        [M](const std::vector<double>& n) {
            std::vector<double> v(M, 1.0);
            double tot = 0;
            for (std::size_t i = 0; i < M; ++i) tot += n[i];
            if (tot > 0)
                for (std::size_t i = 0; i < M; ++i) v[i] = n[i] / tot;
            return v;
        },
        std::vector<double>{1.0});

    section("CTMC");
    print_avg(solve_avg("CTMC", m), "CTMC (exact):");

    // SolverSSA carries the SAME factorization on the sample path: phi(n) is a
    // CONSTANT within a state, so it is evaluated once per state and multiplies
    // every station service rate there. Its NRM engine cannot (its propensity
    // closures see one station's population slice), so the model is routed to
    // the serial engine.
    section("SSA");
    SolverOpts ssa;
    ssa.seed = 23000;
    ssa.samples = 200000;
    print_avg(solve_avg("SSA", m, ssa), "SSA (serial, 2e5 samples):");

    // Only SolverCTMC and SolverSSA declare GlobalDependence; every other solver
    // refuses the model rather than solving it unscaled.
    note("\nOnly SolverCTMC and SolverSSA are solved: phi(n) reads the whole state,");
    note("      which no per-station scaling can express, so the other feature sets");
    note("      refuse it.");

    // The balance property is checkable, and is what separates a Whittle network
    // from an arbitrary state-dependent rate.
    const double viol = line::sn::sn_gd_balance<double>(
        [](const std::vector<double>& n) {
            double tot = 0;
            for (std::size_t i = 0; i < n.size(); ++i) tot += n[i];
            std::vector<double> v(n.size(), 0.0);
            if (tot > 0)
                for (std::size_t i = 0; i < n.size(); ++i) v[i] = n[i] / tot;
            return v;
        },
        std::vector<std::size_t>{4, 4});
    note("\nworst relative balance violation: " + std::to_string(viol) +
         " (balanced when ~0)");
}

LINE_EXAMPLE("advanced/loadDependent", ld_global_dependence);

namespace {

/** Mixed-radix state code over the bounded lattice. */
std::size_t bf_index(const std::vector<std::size_t>& n, std::size_t base) {
    std::size_t idx = 0, mult = 1;
    for (std::size_t s = 0; s < n.size(); ++s) {
        idx += n[s] * mult;
        mult *= base;
    }
    return idx;
}

/** Phi(n) = max_l (1/C_l) sum_{s in l} Phi(n-e_s), Phi(0) = 1. */
std::vector<double> bf_balance_function(const std::vector<std::vector<int>>& A,
                                        const std::vector<double>& C, std::size_t cutoff) {
    const std::size_t S = A[0].size(), base = cutoff + 1;
    std::size_t ns = 1;
    for (std::size_t s = 0; s < S; ++s) ns *= base;
    std::vector<double> Phi(ns, 0.0);
    Phi[0] = 1.0;
    // walk in increasing total population, so every Phi(n-e_s) is already set
    std::vector<std::pair<std::size_t, std::size_t>> order;
    order.reserve(ns);
    for (std::size_t k = 0; k < ns; ++k) {
        std::size_t rem = k, tot = 0;
        for (std::size_t s = 0; s < S; ++s) {
            tot += rem % base;
            rem /= base;
        }
        order.push_back(std::make_pair(tot, k));
    }
    std::stable_sort(order.begin(), order.end());
    for (std::size_t oi = 0; oi < order.size(); ++oi) {
        const std::size_t k = order[oi].second;
        std::vector<std::size_t> n(S, 0);
        std::size_t rem = k, tot = 0;
        for (std::size_t s = 0; s < S; ++s) {
            n[s] = rem % base;
            rem /= base;
            tot += n[s];
        }
        if (tot == 0) continue;
        double best = 0;
        for (std::size_t l = 0; l < A.size(); ++l) {
            double acc = 0;
            for (std::size_t s = 0; s < S; ++s) {
                if (A[l][s] > 0 && n[s] > 0) {
                    n[s] -= 1;
                    acc += Phi[bf_index(n, base)];
                    n[s] += 1;
                }
            }
            best = std::max(best, acc / C[l]);
        }
        Phi[k] = best;
    }
    return Phi;
}

}  // namespace

void ld_whittle_bandwidth() {
    note("=== Load-Dependent Service - Open Whittle (Bandwidth Sharing) ===");

    // One route holds SEVERAL links at once, which no per-station rate scaling
    // can express: the 2-link linear network, route 1 crossing both links and
    // routes 2 and 3 one link each. Capacity is shared by BALANCED FAIRNESS,
    // x_s(n) = Phi(n-e_s)/Phi(n), which satisfies the Whittle balance property by
    // construction, so pi(n) ~ Phi(n) prod rho_s^n_s and the model is insensitive.
    const std::vector<std::vector<int>> A = {{1, 1, 0}, {1, 0, 1}};
    const std::vector<double> C = {1.0, 1.0};
    const double nu[3] = {0.20, 0.30, 0.30};
    const double mu[3] = {1.00, 1.00, 1.00};
    const std::size_t CUTOFF = 3, S = 3, base = CUTOFF + 1;
    const std::vector<double> Phi = bf_balance_function(A, C, CUTOFF);

    Net m("model");
    Source src(m, "Source");
    std::vector<std::size_t> routes(S);
    for (std::size_t s = 0; s < S; ++s)
        routes[s] = m.add_queue("Route" + std::to_string(s + 1), SchedStrategy::PS);
    Sink snk(m, "Sink");
    std::vector<std::size_t> cls(S);
    for (std::size_t s = 0; s < S; ++s) {
        cls[s] = m.add_open_class("Route" + std::to_string(s + 1) + "Flows");
        src.set_arrival(cls[s], Exp(nu[s]));
    }
    for (std::size_t s = 0; s < S; ++s)
        for (std::size_t t = 0; t < S; ++t)
            // BOTH ARMS AS `D`: Exp and Disabled are distinct subclasses of
            // lang::Distrib<double>, so an unwrapped ternary has no common type.
            m.set_service(routes[t], cls[s], s == t ? D(Exp(mu[s])) : D(Disabled()));
    qn::RoutingMatrix<double> P;
    for (std::size_t s = 0; s < S; ++s) {
        P.set(cls[s], cls[s], src, routes[s], 1.0);
        P.set(cls[s], cls[s], routes[s], snk, 1.0);
    }
    m.link(P);

    const qn::NetworkStruct<double>& sn0 = m.get_struct();
    std::vector<std::size_t> idx(S);
    for (std::size_t s = 0; s < S; ++s) idx[s] = sn0.nodes[routes[s] - 1].station - 1;
    const std::size_t M = sn0.nstations, K = sn0.nclasses;

    m.set_global_dependence(
        [Phi, idx, M, K, S, base](const std::vector<double>& n) {
            std::vector<std::size_t> npop(S, 0);
            std::size_t tot = 0;
            for (std::size_t s = 0; s < S; ++s) {
                npop[s] = static_cast<std::size_t>(n[idx[s] * K + s] + 0.5);
                tot += npop[s];
            }
            std::vector<double> v(M * K, 1.0);
            if (tot == 0) return v;
            const double den = Phi[bf_index(npop, base)];
            for (std::size_t s = 0; s < S; ++s) {
                if (npop[s] > 0) {
                    npop[s] -= 1;
                    v[idx[s] * K + s] = Phi[bf_index(npop, base)] / den;
                    npop[s] += 1;
                } else {
                    v[idx[s] * K + s] = 0.0;
                }
            }
            return v;
        },
        // The third argument is the per-slot open-class truncation used when phi
        // is materialized onto the JSON wire; it matches the cutoff below, which
        // is also the range over which the balance function Phi was built.
        std::vector<double>{1.0}, static_cast<int>(CUTOFF));

    section("CTMC");
    SolverOpts o;
    o.cutoff = static_cast<double>(CUTOFF);
    print_avg(solve_avg("CTMC", m, o), "CTMC (exact, truncated):");

    note("\nTruncating a REVERSIBLE chain preserves the conditional law, so the");
    note("      means agree with the closed-form product form rather than only");
    note("      approximating it.");
}

LINE_EXAMPLE("advanced/loadDependent", ld_whittle_bandwidth);

}  // namespace examples
}  // namespace line


