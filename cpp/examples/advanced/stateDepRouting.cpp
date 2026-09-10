/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/advanced/stateDepRouting/`: round-robin and join-the-
 * shortest-queue dispatching.
 *
 * WHAT THIS PORT DOES WITH A STATE-DEPENDENT STRATEGY, and why three of these
 * four examples report a refusal instead of a table.
 *
 *   RROBIN   `NetworkStruct::refresh_routing` expands it uniformly, because the
 *            routing PROBABILITIES a round-robin pointer induces are uniform
 *            and only the higher moments of the split are deterministic. A
 *            solver that does not carry the pointer in its state must therefore
 *            not declare `RoutingStrategy_RROBIN`, or it answers a random-
 *            routing model under a round-robin name -- and none of the solvers
 *            these examples call does declare it, so each refuses BY NAME
 *            through its feature set. The refusal is printed and not rethrown:
 *            it is the port's answer to the model.
 *   JSQ      refused one step earlier, when the struct is built, so no solver
 *            can be asked at all.
 *
 * JMT and LDES are refused by name as well: this port carries neither.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "example_util.h"
#include "example_node_table.h"
#include "examples_common.h"
#include "line/api/pfqn/pfqn_sdr.h"
#include "line/lang/dist_fitters.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/nc/solver_nc.h"
#include "line/solvers/ssa/ssa_dispatch.h"

namespace line {
namespace examples {

namespace {

/** The CTMC table, or the refusal that stands in for it. */
void ctmc_table(const Sn& sn, const ctmc::CtmcOptions& opt) {
    try {
        print_avg(sn, ctmc::solver_ctmc_avg_table(sn, ctmc::solver_ctmc_analyzer(sn, opt),
                                                  opt.method));
    } catch (const Error& e) {
        std::printf("N/A: %s\n", e.what());
    }
}

/** The same, as the NODE table the state-dependent-routing references print. */
void ctmc_node_table(const Sn& sn, const ctmc::CtmcOptions& opt) {
    try {
        print_avg_node(sn, ctmc::solver_ctmc_avg_table(sn, ctmc::solver_ctmc_analyzer(sn, opt),
                                                       opt.method));
    } catch (const Error& e) {
        std::printf("N/A: %s\n", e.what());
    }
}

/** The SSA table, or the refusal that stands in for it. */
void ssa_table(const Sn& sn, std::size_t samples) {
    ssa::SsaOptions opt;
    opt.samples = samples;
    opt.seed = 23000;
    try {
        print_avg_sim(sn, ssa::solver_ssa(sn, opt));
        note("Note: the serial SSA engine simulates the routing table it is given, so an RROBIN "
             "dispatcher reaches it as its uniform expansion and the split is random rather than "
             "cyclic.");
    } catch (const Error& e) {
        std::printf("N/A: %s\n", e.what());
    }
}

/** One matrix printed as a row of six-figure numbers, in the reference's layout. */
std::string row_of(const Matrix<double>& m) {
    std::string out = "[";
    char buf[64];
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j) {
            std::snprintf(buf, sizeof(buf), "%.6g", m(i, j));
            if (out.size() > 1) out += " ";
            out += buf;
        }
    return out + "]";
}

/** Largest elementwise absolute difference of two equally shaped matrices. */
double max_abs_diff(const Matrix<double>& a, const Matrix<double>& b) {
    double worst = 0.0;
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < a.cols(); ++j)
            worst = std::max(worst, std::fabs(a(i, j) - b(i, j)));
    return worst;
}

/** `MAP([[..]],[[..]])` written out, the shape sdroute_twoclasses_closed uses. */
D map2(double a00, double a01, double a10, double a11, double b00, double b01, double b10,
       double b11) {
    Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D0(0, 0) = a00;
    D0(0, 1) = a01;
    D0(1, 0) = a10;
    D0(1, 1) = a11;
    D1(0, 0) = b00;
    D1(0, 1) = b01;
    D1(1, 0) = b10;
    D1(1, 1) = b11;
    return D::map_dist(D0, D1, lang::ProcessType::MAP);
}

}  // namespace

/**
 * `sdroute_closed.py`: a closed model whose Delay dispatches round-robin over
 * itself and the two queues.
 */
void sdroute_closed() {
    Net m("model");
    Delay delay(m, "Delay");
    Queue q1(m, "Queue1", SchedStrategy::PS);
    Queue q2(m, "Queue2", SchedStrategy::PS);
    ClosedClass c1(m, "Class1", 1, delay, 0);
    delay.set_service(c1, lang::hyperexp_fit_mean_scv<double>(1.0, 25.0));
    q1.set_service(c1, Exp(1.0));
    q2.set_service(c1, Exp(2.0));

    // The links of `addLink`: the Delay reaches itself and both queues, and each
    // queue returns to it with probability one.
    Routing P;
    P.set(c1, c1, delay, delay, 1.0);
    P.set(c1, c1, delay, q1, 1.0);
    P.set(c1, c1, delay, q2, 1.0);
    P.set(c1, c1, q1, delay, 1.0);
    P.set(c1, c1, q2, delay, 1.0);
    m.set_routing(delay, c1, RoutingStrategy::RROBIN);
    m.link(P);

    const Sn& sn = m.get_struct();

    section("CTMC");
    ctmc_table(sn, ctmc::CtmcOptions());

    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 100000)));

    section("SSA");
    ssa_table(sn, 10000);
}

/**
 * `sdroute_twoclasses_closed.py`: the same dispatcher with two classes and
 * Markovian service, renewal (APH, PH) at Queue1 and non-renewal (MAP) at
 * Queue2.
 */
void sdroute_twoclasses_closed() {
    Net m("model");
    Delay delay(m, "Delay");
    Queue q1(m, "Queue1", SchedStrategy::PS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    ClosedClass c1(m, "Class1", 1, delay, 0);
    ClosedClass c2(m, "Class2", 2, delay, 0);

    Matrix<double> A21(2, 2, 0.0);
    A21(0, 0) = -2.0;
    A21(0, 1) = 2.0;
    A21(1, 1) = -0.5;
    const std::vector<double> alpha21{1.0, 0.0};
    Matrix<double> A22(1, 1, -1.0);
    const std::vector<double> alpha22{1.0};

    delay.set_service(c1, lang::hyperexp_fit_mean_scv<double>(1.0, 25.0));
    q1.set_service(c1, APH(alpha21, A21));
    q2.set_service(c1, map2(-20.0, 0.0, 0.0, -1.0, 0.0, 20.0, 0.8, 0.2));
    delay.set_service(c2, lang::hyperexp_fit_mean_scv<double>(1.0, 25.0));
    q1.set_service(c2, PH(alpha22, A22));
    q2.set_service(c2, map2(-4.0, 3.0, 4.0, -6.0, 1.0, 0.0, 0.0, 2.0));

    Routing P;
    for (std::size_t r : {c1, c2}) {
        P.set(r, r, delay, delay, 1.0);
        P.set(r, r, delay, q1, 1.0);
        P.set(r, r, delay, q2, 1.0);
        P.set(r, r, q1, delay, 1.0);
        P.set(r, r, q2, delay, 1.0);
        m.set_routing(delay, r, RoutingStrategy::RROBIN);
    }
    m.link(P);

    const Sn& sn = m.get_struct();

    section("CTMC");
    ctmc_table(sn, ctmc::CtmcOptions());

    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 100000)));

    section("SSA");
    ssa_table(sn, 5000);
}

/**
 * `sdroute_open.py`: an open model whose Router alternates between two queues.
 *
 * The reference prints `getAvgNodeTable` and so does this twin. It used to
 * print the STATION table and say the port had no node-level one, which was
 * already untrue -- `line-cli -a node` has been printing it -- and the
 * difference is the whole model: the Router is the dispatcher under test and is
 * not a station, so it was absent from every row.
 */
void sdroute_open() {
    Net m("myModel");
    Source source(m, "Source");
    Router router(m, "Router");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1");
    source.set_arrival(c1, Exp(1.0));
    q1.set_service(c1, Exp(2.0));
    q2.set_service(c1, Exp(2.0));

    Routing P;
    P.set(c1, c1, source, router, 1.0);
    P.set(c1, c1, router, q1, 1.0);
    P.set(c1, c1, router, q2, 1.0);
    P.set(c1, c1, q1, sink, 1.0);
    P.set(c1, c1, q2, sink, 1.0);
    // RAND on the linked nodes, as `addLink` leaves them: JMT draws from the
    // node's own stream for a random split even at one destination, so PROB
    // there would shift the whole sample path at the same seed.
    m.set_routing(source, c1, RoutingStrategy::RAND);
    m.set_routing(q1, c1, RoutingStrategy::RAND);
    m.set_routing(q2, c1, RoutingStrategy::RAND);
    m.set_routing(sink, c1, RoutingStrategy::RAND);
    m.set_routing(router, c1, RoutingStrategy::RROBIN);
    m.link(P);

    const Sn& sn = m.get_struct();

    section("JMT");
    print_avg_node(sn, jmt_avg(m, sim_opts(23000)));

    // The cutoff is presentation, not a second solver: the golden keys this
    // table `CTMC`, and filing it under the decorated spelling would read as a
    // solver that produced nothing. See `section(title, golden_key)`.
    section("CTMC (cutoff = 5)", "CTMC");
    ctmc::CtmcOptions copt;
    copt.cutoff = 5.0;
    ctmc_node_table(sn, copt);
}

/**
 * `sdroute_jsq.py`: join-the-shortest-queue dispatching over three queues.
 *
 * NO SOLVER RUNS. The reference solves it with JMT and LDES, neither of which
 * this port carries, and the model itself is unreachable here as well:
 * `refresh_routing` refuses JSQ by name when the struct is built, because
 * expanding it uniformly would answer a random-routing model under a
 * shortest-queue name. The model is built here so the refusal names a real
 * model rather than a hypothetical one.
 */
void sdroute_jsq() {
    const std::size_t N = 3;
    const double rho = 0.7, mu = 1.0;
    const double lam = static_cast<double>(N) * rho * mu;

    Net m("myModel");
    Source source(m, "Source");
    Router router(m, "Router");
    std::vector<std::size_t> queue;
    for (std::size_t i = 0; i < N; ++i)
        queue.push_back(m.add_queue("Queue" + std::to_string(i + 1), SchedStrategy::FCFS));
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1");
    source.set_arrival(c1, Exp(lam));
    for (std::size_t i = 0; i < N; ++i) m.set_service(queue[i], c1, Exp(mu));

    Routing P;
    P.set(c1, c1, source, router, 1.0);
    for (std::size_t i = 0; i < N; ++i) {
        P.set(c1, c1, router, queue[i], 1.0);
        P.set(c1, c1, queue[i], sink, 1.0);
    }
    // `addLink` LEAVES A NODE ON RAND, and this is not cosmetic: JMT draws from
    // that node's stream for a random split even when it has ONE destination,
    // so a model built with an explicit routing matrix (PROB) consumes the
    // stream differently and its sample path parts company with the
    // reference's at the same seed. The reference links every node this way.
    m.set_routing(source, c1, RoutingStrategy::RAND);
    for (std::size_t i = 0; i < N; ++i) m.set_routing(queue[i], c1, RoutingStrategy::RAND);
    m.set_routing(sink, c1, RoutingStrategy::RAND);
    m.set_routing(router, c1, RoutingStrategy::JSQ);
    m.link(P);

    // BOTH SIMULATORS RUN THIS MODEL. `refresh_routing` lists JSQ among the
    // strategies it EXPANDS, and the routing itself is resolved at dispatch time
    // by each engine: JMT writes a "Join the Shortest Queue (JSQ)" strategy into
    // the JSIM document and the LDES engine picks the least-loaded candidate.
    // NODE tables, which is what the reference prints: the Router is the
    // dispatcher under test and is not a station.
    section("JMT");
    print_avg_node(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    section("LDES");
    print_avg_node(m.get_struct(), ldes_avg(m, sim_opts(23000)));
}

/**
 * Product-form state-dependent routing, A. E. Krzesinski, "Multiclass Queueing
 * Networks with State-Dependent Routing", Perform. Eval. 7:125-143, 1987.
 *
 * Unlike RROBIN and JSQ above, this routing is state dependent and STILL has a
 * product form, eq. (16), which `solver_nc_sdr` evaluates exactly. Sec. 2.5's
 * central server with two peripheral centres and two levels of nesting: with
 * C = (-1,-1), d_12 = 1, d_13 = 2 and d_23 = 3 the probabilities out of the
 * central server are Table 1 of the paper, and the branch populations are
 * capped at m_2 <= 1 and m_3 <= 3 by the coefficients themselves. When both
 * branches are full the customer returns to the central server and is served
 * again before retrying, the busy form of waiting of Sec. 2.5, which is why the
 * entry node carries a self-loop.
 */
void sdroute_krzesinski() {
    const double mu[3] = {1.0, 0.8, 0.5};
    Net m("sdrCentralServer");
    Queue cpu(m, "CPU", SchedStrategy::FCFS);
    Queue d1(m, "Disk1", SchedStrategy::FCFS);
    Queue d2(m, "Disk2", SchedStrategy::FCFS);
    ClosedClass c1(m, "Class1", 3, cpu, 0);
    cpu.set_service(c1, Exp(mu[0]));
    d1.set_service(c1, Exp(mu[1]));
    d2.set_service(c1, Exp(mu[2]));

    Routing P;
    P.set(c1, c1, cpu, cpu, 1.0 / 3);  // denied entry: the busy form of waiting
    P.set(c1, c1, cpu, d1, 1.0 / 3);
    P.set(c1, c1, cpu, d2, 1.0 / 3);
    P.set(c1, c1, d1, cpu, 1.0);
    P.set(c1, c1, d2, cpu, 1.0);
    m.link(P);

    Matrix<double> dcoeff(2, 3, 0.0);
    dcoeff(0, 1) = 1.0;
    dcoeff(0, 2) = 2.0;
    dcoeff(1, 2) = 3.0;
    m.set_state_dep_routing(cpu, cpu, {{}, {d1}, {d2}}, {0, 1, 2}, {-1, -1}, dcoeff);

    const Sn& sn = m.get_struct();

    // `note`, not `section`: this block is a table of routing probabilities, not
    // a solver run, and a parity record filed under it would read as a solver
    // that produced no output. See the section(title, golden_key) doc.
    note("\nSDR routing probabilities (Table 1, eq. 10)");
    const pfqn::SdrCoeff co = pfqn::pfqn_sdrcoeff(sn.sdr);
    const int tab[6][2] = {{0, 0}, {0, 1}, {1, 0}, {0, 2}, {1, 1}, {1, 2}};
    for (int t = 0; t < 6; ++t) {
        const std::vector<double> n = {3.0 - tab[t][0] - tab[t][1], double(tab[t][0]),
                                       double(tab[t][1])};
        const std::vector<double> Pb = pfqn::pfqn_sdrprob(co, n);
        std::printf("  (m2,m3)=(%d,%d)  P12=%.4f  P13=%.4f  P11=%.4f\n", tab[t][0], tab[t][1],
                    Pb[1], Pb[2], pfqn::pfqn_sdrped(Pb));
    }

    section("NC");
    print_avg(solve_avg("NC", m));

    // The generator carries eq. (10) per state through `rt_state`, so this is
    // the same model solved twice by unrelated routes: an exact chain and a
    // product form.
    section("CTMC");
    ctmc_table(sn, ctmc::CtmcOptions());
}

/**
 * `sdroute_multibranch.m`: a branch holding SEVERAL centres.
 *
 * Eq. (16) still holds, but the coefficients xi are the branch traffic
 * equations rather than Section 3.2's shorthand xi = xi_e, which is exact only
 * when the branch departure centre is visited once. Case A is a plain series
 * branch, where the two readings coincide; case B feeds the branch departure
 * centre back onto its entry centre, where they do not -- the traffic-equation
 * xi reproduces the exact CTMC to machine precision and the literal xi = 1 is
 * out by 2.98e-1 in Q. See _kb/16-state-dependent-routing.md
 */
void sdroute_multibranch() {
    const double mu[4] = {1.0, 0.9, 0.7, 0.5};  // centres 1, 2a, 2b, 3
    const std::size_t N = 3;
    const double pbacks[2] = {0.0, 0.5};
    const char* labels[2] = {"A: branch 2 = 2a -> 2b (series)",
                             "B: branch 2 = 2a -> 2b, 2b -> 2a w.p. 0.5 (feedback onto the "
                             "branch departure)"};

    std::printf("\n==== SDR with multi-centre branches ====\n");
    for (int cs = 0; cs < 2; ++cs) {
        const double pback = pbacks[cs];
        Net m("sdr_multi");
        Queue cpu(m, "CPU", SchedStrategy::FCFS);
        Queue b2a(m, "B2a", SchedStrategy::FCFS);
        Queue b2b(m, "B2b", SchedStrategy::FCFS);
        Queue b3(m, "B3", SchedStrategy::FCFS);
        ClosedClass c1(m, "Class1", N, cpu, 0);
        const std::size_t centre[4] = {cpu, b2a, b2b, b3};
        for (int i = 0; i < 4; ++i) m.set_service(centre[i], c1, Exp(mu[i]));

        Routing P;
        P.set(c1, c1, cpu, cpu, 1.0);  // denied entry: the busy form of waiting
        P.set(c1, c1, cpu, b2a, 1.0);
        P.set(c1, c1, cpu, b3, 1.0);
        P.set(c1, c1, b2a, b2b, 1.0);
        if (pback > 0) {
            P.set(c1, c1, b2b, b2a, pback);
            P.set(c1, c1, b2b, cpu, 1.0 - pback);
        } else {
            P.set(c1, c1, b2b, cpu, 1.0);
        }
        P.set(c1, c1, b3, cpu, 1.0);
        m.link(P);

        Matrix<double> dcoeff(2, 3, 0.0);
        dcoeff(0, 1) = 2.0;
        dcoeff(0, 2) = 2.0;
        dcoeff(1, 2) = 2.0;
        m.set_state_dep_routing(cpu, cpu, {{}, {b2a, b2b}, {b3}}, {0, 1, 2}, {-1, -1}, dcoeff);

        const Sn& sn = m.get_struct();

        // The state-INDEPENDENT part of the routing: the complement M-V and the
        // arcs inside a branch. The state-dependent arcs out of the entry centre
        // are not part of it.
        Matrix<double> Psir(4, 4, 0.0);
        Psir(0, 0) = 1.0;  // complement M-V = {1}
        Psir(1, 2) = 1.0;  // inside branch 2
        if (pback > 0) Psir(2, 1) = pback;
        const Matrix<double> xi =
            pfqn::pfqn_sdrvisits<double>(sn.sdr, std::vector<Matrix<double>>(1, Psir));
        Matrix<double> S(4, 1, 0.0);
        for (int i = 0; i < 4; ++i) S(i, 0) = 1.0 / mu[i];
        const pfqn::SdrResult<double> res =
            pfqn::pfqn_sdr<double>(S, xi, std::vector<std::size_t>(1, N), sn.sdr);

        const mva::AvgResult<double> ct = ctmc::solver_ctmc_avg_table(
            sn, ctmc::solver_ctmc_analyzer(sn, ctmc::CtmcOptions()), "default");

        std::printf("\n  %s\n", labels[cs]);
        std::printf("    xi          = %s\n", row_of(xi).c_str());
        std::printf("    pfqn_sdr Q  = %s\n", row_of(res.QN).c_str());
        std::printf("    CTMC     Q  = %s\n", row_of(ct.QN).c_str());
        std::printf("    pfqn_sdr X  = %s\n", row_of(res.XN).c_str());
        std::printf("    CTMC     X  = %s\n", row_of(ct.TN).c_str());
        std::printf("    max|dQ| = %.3e   max|dX| = %.3e\n", max_abs_diff(ct.QN, res.QN),
                    max_abs_diff(ct.TN, res.XN));

        // the reading asserted verbatim by the paper for E+D: xi = xi_e everywhere
        Matrix<double> xiflat(4, 1, 1.0);
        const pfqn::SdrResult<double> flat =
            pfqn::pfqn_sdr<double>(S, xiflat, std::vector<std::size_t>(1, N), sn.sdr);
        std::printf("    xi==1 everywhere: max|dQ| = %.3e   max|dX| = %.3e\n",
                    max_abs_diff(ct.QN, flat.QN), max_abs_diff(ct.TN, flat.XN));
    }
}

LINE_EXAMPLE("advanced/stateDepRouting", sdroute_closed);
LINE_EXAMPLE("advanced/stateDepRouting", sdroute_twoclasses_closed);
LINE_EXAMPLE("advanced/stateDepRouting", sdroute_open);
LINE_EXAMPLE("advanced/stateDepRouting", sdroute_jsq);
LINE_EXAMPLE("advanced/stateDepRouting", sdroute_krzesinski);
LINE_EXAMPLE("advanced/stateDepRouting", sdroute_multibranch);

}  // namespace examples
}  // namespace line
