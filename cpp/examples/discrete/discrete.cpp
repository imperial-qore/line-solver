/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `matlab/examples/discrete`: models on a slot lattice rather than on the
 * continuous time axis.
 *
 * Time advances one slot at a time; in each slot a job in service completes
 * with probability p and an arrival occurs with probability b, both recorded at
 * the end of the slot with the departure resolved before the arrival (Daduna's
 * LA rule and D/A rule). Every rate below is therefore a per-slot probability
 * and every time a number of slots.
 *
 * SolverNC enters this route only when `NcSolverOptions::slotted` is set, the
 * same switch SolverLDES uses; it is NEVER inferred from the presence of a
 * Geometric distribution, because a Geometric service time is an ordinary
 * continuous-time model unless the caller declares the lattice. When the switch
 * is on and the model falls outside the discrete-time product form, the solver
 * refuses it instead of approximating.
 *
 * Reference: H. Daduna, Queueing Networks with Discrete Time Scale, LNCS 2046,
 * Springer, 2001.
 */

#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

#include "example_util.h"
#include "examples_common.h"
#include "line/api/dpfqn/dpfqn_nc.h"
#include "line/api/dqsys/dqsys_bernoulli1.h"
#include "line/solvers/nc/solver_nc_runner.h"

namespace line {
namespace examples {

namespace {

/** SolverNC on the slot lattice. */
void run_nc_slotted(const Sn& sn) {
    nc::NcSolverOptions opt;
    opt.slotted = true;
    print_avg(sn, nc::solver_nc_run_analyzer(sn, opt));
}

void print_vec(const char* label, const std::vector<double>& v, std::size_t upto) {
    std::printf("%s", label);
    for (std::size_t n = 0; n < std::min(upto, v.size()); ++n) std::printf(" %.6f", v[n]);
    std::printf("\n");
}

/** The service probabilities the cycle examples share. */
std::vector<double> cycle_rates() {
    std::vector<double> p;
    p.push_back(0.5);
    p.push_back(0.25);
    p.push_back(0.7);
    return p;
}

}  // namespace

// ---------------------------------------------------------------------------
// dt_geogeo1
// ---------------------------------------------------------------------------

/**
 * Geo/Geo/1 with an unbounded buffer.
 *
 * The closed form of theorem 2.3 collapses for constant a and s to the
 * geometric law of corollary 2.7: mean queue length a(1-a)/(s-a) and mean
 * sojourn time (1-a)/(s-a) slots.
 */
void dt_geogeo1() {
    const double a = 0.2, s = 0.5;
    Net m("GeoGeo1");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1", 0);
    source.set_arrival(c1, Geometric(a));
    queue.set_service(c1, Geometric(s));
    link_serial(m, {source, queue, sink});

    section("NC (slotted)");
    run_nc_slotted(m.get_struct());
    std::printf("closed form: E[N] = %g, E[T] = %g slots\n",
                a * (1 - a) / (s - a), (1 - a) / (s - a));
}

// ---------------------------------------------------------------------------
// dt_geogeo1_loss
// ---------------------------------------------------------------------------

/**
 * Geo/Geo/1/L, the loss system of corollary 2.8.
 *
 * An arrival that lands in a slot which already holds L jobs is lost, i.e.
 * b(n) = 0 for n >= L, which is exactly the assumption corollary 2.8 places on
 * the arrival probabilities.
 */
void dt_geogeo1_loss() {
    const double a = 0.2, s = 0.5;
    const std::size_t L = 4;
    Net m("GeoGeo1L");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1", 0);
    source.set_arrival(c1, Geometric(a));
    queue.set_service(c1, Geometric(s));
    queue.set_capacity(static_cast<double>(L));
    link_serial(m, {source, queue, sink});

    section("NC (slotted)");
    run_nc_slotted(m.get_struct());

    const dqsys::Bernoulli1Result<double> r =
            dqsys::dqsys_bernoulli1(std::vector<double>(1, a), std::vector<double>(1, s), L);
    print_vec("queue length law   :", r.pmf, L + 1);
    std::printf("loss probability   : %g\n", r.lossProb);
    std::printf("carried throughput : %g of the %g offered per slot\n", r.throughput, a);
}

// ---------------------------------------------------------------------------
// dt_bernoulli_loaddep
// ---------------------------------------------------------------------------

/**
 * Load-dependent Bernoulli server, and the discrete-time arrival theorem.
 *
 * Example 2.10 notes that a discrete-time M/M/c queue has no exactly equivalent
 * state dependent single server, but that p(n) = p min(n,c) reproduces its
 * conditional service intensity. The arrival law of theorem 2.11 is printed
 * next to the time-stationary one: discrete time has no PASTA analogue, and the
 * two differ even for a state independent Bernoulli arrival stream.
 */
void dt_bernoulli_loaddep() {
    const double a = 0.6, s = 0.3;
    const std::size_t c = 3, L = 20;
    Net m("LoadDepBernoulli");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1", 0);
    source.set_arrival(c1, Geometric(a));
    queue.set_service(c1, Geometric(s));
    queue.set_capacity(static_cast<double>(L));
    std::vector<double> alpha(L);
    for (std::size_t n = 1; n <= L; ++n) alpha[n - 1] = static_cast<double>(std::min(n, c));
    queue.set_load_dependence(alpha);
    link_serial(m, {source, queue, sink});

    section("NC (slotted)");
    run_nc_slotted(m.get_struct());

    std::vector<double> p(L);
    for (std::size_t n = 1; n <= L; ++n) p[n - 1] = s * static_cast<double>(std::min(n, c));
    const dqsys::Bernoulli1Result<double> r =
            dqsys::dqsys_bernoulli1(std::vector<double>(1, a), p, L);
    print_vec("time-stationary law pi(0..6)  :", r.pmf, 7);
    print_vec("arrival law         pi_1(0..6):", r.arrivalPmf, 7);
    std::printf("no PASTA in discrete time: the two rows above are different laws\n");
}

// ---------------------------------------------------------------------------
// dt_cycle
// ---------------------------------------------------------------------------

/**
 * Closed cycle of Bernoulli servers.
 *
 * The stationary queue length vector has the product form of corollary 3.4,
 * whose extra factor (1/q_j) on the busy nodes is what separates it from the
 * continuous-time Gordon-Newell form: a homogeneous cycle is uniform on the
 * state space in continuous time and is not here.
 */
void dt_cycle() {
    const std::vector<double> p = cycle_rates();
    const std::size_t N = 5;
    Net m("BernoulliCycle");
    std::vector<std::size_t> station;
    for (std::size_t j = 0; j < p.size(); ++j) {
        station.push_back(m.add_queue("Queue" + std::to_string(j + 1), SchedStrategy::FCFS));
    }
    ClosedClass c1(m, "Jobs", static_cast<double>(N), station[0], 0);
    for (std::size_t j = 0; j < p.size(); ++j) m.set_service(station[j], c1, Geometric(p[j]));
    link_serial(m, station);

    section("NC (slotted)");
    run_nc_slotted(m.get_struct());

    const dpfqn::DtNcResult<double> nc = dpfqn::dpfqn_nc(p, N);
    std::printf("log G(N,J)   = %g\n", nc.lG);
    std::printf("throughput   = %g jobs per slot\n", nc.throughput());
    std::printf("utilizations =");
    for (std::size_t j = 0; j < p.size(); ++j) std::printf(" %.6f", nc.throughput() / p[j]);
    std::printf("\n");
}

// ---------------------------------------------------------------------------
// dt_cycle_loaddep
// ---------------------------------------------------------------------------

/**
 * Closed cycle of state dependent Bernoulli servers.
 *
 * Station 2 runs at p_2(n) = p_2 min(n,2), the discrete-time analogue of adding
 * a second server. That is admissible only because the dependence is expressed
 * as a state dependent SINGLE server: a genuine multiserver node inside a cycle
 * of geometrical queues destroys the product form for every finite server count
 * (Pestien and Ramakrishnan, cited before example 2.10), and the analyzer
 * refuses one rather than approximating it.
 */
void dt_cycle_loaddep() {
    const std::vector<double> p = cycle_rates();
    const std::size_t N = 5;
    Net m("BernoulliCycleLD");
    std::vector<std::size_t> station;
    for (std::size_t j = 0; j < p.size(); ++j) {
        station.push_back(m.add_queue("Queue" + std::to_string(j + 1), SchedStrategy::FCFS));
    }
    ClosedClass c1(m, "Jobs", static_cast<double>(N), station[0], 0);
    for (std::size_t j = 0; j < p.size(); ++j) m.set_service(station[j], c1, Geometric(p[j]));
    std::vector<double> alpha(N);
    for (std::size_t n = 1; n <= N; ++n) alpha[n - 1] = static_cast<double>(std::min<std::size_t>(n, 2));
    m.set_load_dependence(station[1], alpha);
    link_serial(m, station);

    section("NC (slotted)");
    run_nc_slotted(m.get_struct());

    std::vector<std::vector<double> > P(p.size(), std::vector<double>(N));
    for (std::size_t j = 0; j < p.size(); ++j) {
        for (std::size_t n = 1; n <= N; ++n) {
            P[j][n - 1] = (j == 1) ? p[j] * static_cast<double>(std::min<std::size_t>(n, 2)) : p[j];
        }
    }
    const dpfqn::DtNcLdResult<double> nc = dpfqn::dpfqn_ncld(P, N);
    print_vec("P(X_2 = 0..N) =", nc.marginal(1), N + 1);
}

// ---------------------------------------------------------------------------
// dt_cycle_multiclass
// ---------------------------------------------------------------------------

/**
 * Multichain closed cycle.
 *
 * Service in the cycle is type independent and FCFS forbids overtaking, so the
 * cyclic order of the jobs is frozen for all time and the joint queue length
 * law is the unichain one at the aggregate population. Each chain then holds a
 * share of every station equal to its share of the population, which is the
 * sense in which section 3.2 calls the multichain case a direct adaptation.
 */
void dt_cycle_multiclass() {
    const std::vector<double> p = cycle_rates();
    const std::size_t N1 = 3, N2 = 2;
    Net m("BernoulliCycleMC");
    std::vector<std::size_t> station;
    for (std::size_t j = 0; j < p.size(); ++j) {
        station.push_back(m.add_queue("Queue" + std::to_string(j + 1), SchedStrategy::FCFS));
    }
    ClosedClass c1(m, "Chain1", static_cast<double>(N1), station[0], 0);
    ClosedClass c2(m, "Chain2", static_cast<double>(N2), station[0], 0);
    for (std::size_t j = 0; j < p.size(); ++j) {
        m.set_service(station[j], c1, Geometric(p[j]));
        m.set_service(station[j], c2, Geometric(p[j]));
    }
    // Both chains follow the same cycle and never switch class.
    Routing P;
    const std::size_t J = station.size();
    for (std::size_t j = 0; j < J; ++j) {
        P.set(c1, c1, station[j], station[(j + 1) % J], 1.0);
        P.set(c2, c2, station[j], station[(j + 1) % J], 1.0);
    }
    m.link(P);

    section("NC (slotted)");
    run_nc_slotted(m.get_struct());

    const dpfqn::DtNcResult<double> nc = dpfqn::dpfqn_nc(p, N1 + N2);
    const double X = nc.throughput();
    std::printf("aggregate throughput  = %g jobs per slot\n", X);
    std::printf("per-chain throughput  = %g and %g\n",
                X * N1 / (N1 + N2), X * N2 / (N1 + N2));
}

LINE_EXAMPLE("discrete", dt_geogeo1);
LINE_EXAMPLE("discrete", dt_geogeo1_loss);
LINE_EXAMPLE("discrete", dt_bernoulli_loaddep);
LINE_EXAMPLE("discrete", dt_cycle);
LINE_EXAMPLE("discrete", dt_cycle_loaddep);
LINE_EXAMPLE("discrete", dt_cycle_multiclass);

}  // namespace examples
}  // namespace line
