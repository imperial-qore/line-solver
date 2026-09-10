/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `matlab/examples/advanced/networkCalculus/`,
 * `python/examples/advanced/networkCalculus/`: stochastic network calculus.
 *
 * Every other solver in LINE answers with a MEAN. The 'snc' family of SolverBA
 * answers with a TAIL: given a violation probability eps it returns a delay d
 * for which P{D > d} <= eps holds, and the guarantee is valid for any
 * work-conserving scheduling policy at the station. That is the quantity a
 * service-level objective is written against.
 *
 * The models are an M/M/1 and a feed-forward tandem, whose exact answers are
 * known in closed form, so every number printed can be checked against them.
 */

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "example_util.h"
#include "examples_common.h"
#include "line/api/snc/snc_bound_delay.h"
#include "line/api/snc/snc_conv.h"
#include "line/api/snc/snc_env_poisson.h"
#include "line/api/snc/snc_perc_delay.h"
#include "line/api/snc/snc_srv_exp.h"
#include "line/solvers/ba/solver_ba_runner.h"
#include "line/solvers/ba/solver_ba_snc.h"

namespace line {
namespace examples {

namespace {

mva::AvgResult<double> run_snc(const Sn& sn) {
    ba::BaOptions opt;
    opt.method = "snc.upper";
    return ba::solver_ba_run_analyzer(sn, opt);
}

void snc_delay_quantile() {
    note("=== Stochastic network calculus: delay quantile ===");

    const double lambda = 0.6, mu = 1.0;
    Net m("SncQuantile");
    Source src(m, "Source");
    Queue q(m, "Queue", SchedStrategy::FCFS);
    Sink snk(m, "Sink");
    OpenClass c(m, "Class1", 0);
    src.set_arrival(c, Exp(lambda));
    q.set_service(c, Exp(mu));
    Routing P;
    P.set(src, q, 1.0);
    P.set(q, snk, 1.0);
    m.link(P);
    const Sn sn = m.get_struct();

    // The quantile is the native output of the family: the smallest d and b for
    // which P{D>d} <= eps and P{Q>b} <= eps are certified.
    section("Quantiles at eps = 1e-3 (snc.upper)");
    const ba::SncPercentiles p3 = ba::solver_ba_snc_perc(sn, 1e-3);
    std::printf("%-10s %-10s %12s %12s\n", "Station", "JobClass", "RespTPerc", "QLenPerc");
    std::printf("%-10s %-10s %12.5f %12.5f\n", "Queue", "Class1", p3.D(1, 0), p3.B(1, 0));

    // The exact M/M/1 sojourn tail is exp(-(mu-lambda)*d) and the exact
    // queue-length tail is rho^(n+1). The bound reproduces both DECAY RATES
    // exactly and pays a constant prefactor, so the ratio of the bounded
    // quantile to the exact one falls towards 1 as eps is tightened: the family
    // is at its best exactly where simulation is at its worst.
    section("How the bound tightens as the guarantee gets stricter");
    std::printf("%-8s %10s %10s %8s   %10s %10s %8s\n", "eps", "d bound", "d exact", "ratio",
                "n bound", "n exact", "ratio");
    const double epsList[5] = {1e-2, 1e-3, 1e-6, 1e-9, 1e-12};
    for (int i = 0; i < 5; ++i) {
        const ba::SncPercentiles pp = ba::solver_ba_snc_perc(sn, epsList[i]);
        const double dexact = -std::log(epsList[i]) / (mu - lambda);
        const double nexact = std::log(epsList[i]) / std::log(lambda / mu) - 1.0;
        std::printf("%-8.0e %10.4f %10.4f %8.3f   %10.4f %10.4f %8.3f\n", epsList[i], pp.D(1, 0),
                    dexact, pp.D(1, 0) / dexact, pp.B(1, 0), nexact, pp.B(1, 0) / nexact);
    }

    // The mean columns still work: the response time is the integral of the tail
    // bound, hence an upper bound on the mean. It is loose, and deliberately so
    // -- integrating over the whole axis is dominated by the prefactor rather
    // than by the decay rate the family gets right.
    section("Mean columns (snc.upper), and why they are the loose end");
    const mva::AvgResult<double> r = run_snc(sn);
    const double exactR = 1.0 / (mu - lambda);
    const double exactQ = (lambda / mu) / (1.0 - lambda / mu);
    std::printf("exact M/M/1: R = %.4f, Q = %.4f\n", exactR, exactQ);
    std::printf("snc.upper  : R = %.4f (%.1fx), Q = %.4f (%.1fx)\n", r.RN(1, 0),
                r.RN(1, 0) / exactR, r.QN(1, 0), r.QN(1, 0) / exactQ);

    // The solver is a thin wrapper over line/api/snc. Note snc_srv_exp, not
    // snc_srv_rate: the work unit here is the JOB, so the server is the counting
    // process of an Exp(mu) service, and a constant-rate element would model an
    // M/D/1 and understate the delay.
    section("The same answer from the api");
    const snc::Envelope arv = snc::snc_env_poisson_fn(lambda);
    const snc::Envelope srv = snc::snc_srv_exp_fn(mu);
    const snc::SncResult d = snc::snc_perc_delay(arv, srv, 1e-3);
    std::printf("api: d(1e-3) = %.4f at theta = %.4f\n", d.value, d.theta);
    std::printf("     the optimal theta approaches log(mu/lambda) = %.4f, which is\n",
                std::log(mu / lambda));
    std::printf("     what makes the backlog decay rate exact\n");
    const snc::SncResult back = snc::snc_bound_delay(arv, srv, d.value);
    std::printf("api: P{D > %.4f} <= %.3e (theta = %.4f), exact tail %.3e\n", d.value, back.value,
                back.theta, std::exp(-(mu - lambda) * d.value));
}

LINE_EXAMPLE("advanced/networkCalculus", snc_delay_quantile);

void snc_tandem_multiclass() {
    note("=== Stochastic network calculus: tandem and shared server ===");

    const double lambda = 0.6;
    const double rates[3] = {1.5, 1.2, 1.0};
    Net m("SncTandem");
    Source src(m, "Source");
    Queue q1(m, "Q1", SchedStrategy::FCFS);
    Queue q2(m, "Q2", SchedStrategy::FCFS);
    Queue q3(m, "Q3", SchedStrategy::FCFS);
    Sink snk(m, "Sink");
    OpenClass c(m, "Class1", 0);
    src.set_arrival(c, Exp(lambda));
    q1.set_service(c, Exp(rates[0]));
    q2.set_service(c, Exp(rates[1]));
    q3.set_service(c, Exp(rates[2]));
    Routing P;
    P.set(src, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, q3, 1.0);
    P.set(q3, snk, 1.0);
    m.link(P);

    // Each hop replaces the arrival envelope by the DEPARTURE envelope of the
    // station upstream, which carries the burst the server has added. The exact
    // answer does not degrade this way -- by Burke's theorem the departure
    // process of an M/M/1 is again Poisson -- so the ratio to the exact response
    // time grows hop by hop. The bound stays valid; it is the price of assuming
    // nothing about the departure process beyond its envelope.
    section("Tandem: what envelope propagation costs");
    const mva::AvgResult<double> t = run_snc(m.get_struct());
    std::printf("%-8s %10s %10s %8s\n", "station", "R bound", "R exact", "ratio");
    const char* names[3] = {"Q1", "Q2", "Q3"};
    for (int i = 0; i < 3; ++i) {
        const double exact = 1.0 / (rates[i] - lambda);
        std::printf("%-8s %10.4f %10.4f %8.2f\n", names[i], t.RN(i + 1, 0), exact,
                    t.RN(i + 1, 0) / exact);
    }

    // Summing the per-station bounds pays the burst term at every hop.
    // Concatenating the three service envelopes with snc_conv first and bounding
    // the composed element once pays it only once, the classical result.
    section("End-to-end: pay bursts only once");
    const snc::Envelope arv = snc::snc_env_poisson_fn(lambda);
    const snc::Envelope endToEnd = [&](double theta) {
        snc::Env acc = snc::snc_srv_exp(rates[0], theta);
        for (int i = 1; i < 3; ++i)
            acc = snc::snc_conv(acc, snc::snc_srv_exp(rates[i], theta), theta);
        return acc;
    };
    const snc::SncResult concat = snc::snc_perc_delay(arv, endToEnd, 1e-3);
    double hopByHop = 0.0;
    for (int i = 0; i < 3; ++i)
        hopByHop += snc::snc_perc_delay(arv, snc::snc_srv_exp_fn(rates[i]), 1e-3 / 3.0).value;
    std::printf("  concatenated (snc_conv) : %8.4f  at theta = %.4f\n", concat.value, concat.theta);
    std::printf("  summed per hop          : %8.4f\n", hopByHop);
    std::printf("  pay bursts once saves   : %7.1f%%\n", 100.0 * (1.0 - concat.value / hopByHop));

    // A class sharing a station sees the server minus whatever the other classes
    // take from it: snc_leftover subtracts the cross-flow arrival envelope from
    // the service envelope. The result holds for ANY work-conserving discipline
    // at that station, which is why it is well above the FCFS answer.
    section("Two classes sharing one server (blind multiplexing)");
    Net s("SncShared");
    Source src2(s, "Source");
    Queue qs(s, "Shared", SchedStrategy::FCFS);
    Sink snk2(s, "Sink");
    OpenClass ca(s, "ClassA", 0);
    OpenClass cb(s, "ClassB", 0);
    s.set_arrival(src2, ca, Exp(0.3));
    s.set_arrival(src2, cb, Exp(0.3));
    s.set_service(qs, ca, Exp(1.0));
    s.set_service(qs, cb, Exp(1.0));
    Routing Ps;
    Ps.set(ca, ca, src2, qs, 1.0);
    Ps.set(ca, ca, qs, snk2, 1.0);
    Ps.set(cb, cb, src2, qs, 1.0);
    Ps.set(cb, cb, qs, snk2, 1.0);
    s.link(Ps);
    const Sn ssn = s.get_struct();
    const mva::AvgResult<double> rs = run_snc(ssn);
    const ba::SncPercentiles ps = ba::solver_ba_snc_perc(ssn, 1e-3);
    std::printf("%-10s %10s %12s %12s\n", "JobClass", "R bound", "RespTPerc", "QLenPerc");
    std::printf("%-10s %10.4f %12.4f %12.4f\n", "ClassA", rs.RN(1, 0), ps.D(1, 0), ps.B(1, 0));
    std::printf("%-10s %10.4f %12.4f %12.4f\n", "ClassB", rs.RN(1, 1), ps.D(1, 1), ps.B(1, 1));
    std::printf("exact per-class response time (aggregate M/M/1, lambda=0.6, mu=1): %.4f\n",
                1.0 / (1.0 - 0.6));

    // The elementary envelope algebra has real limits, and the analyzer states
    // them rather than returning a number that looks plausible. Each refusal is
    // a modelling assumption of the calculus, not an implementation gap.
    section("What the family refuses, and why");
    Net u("Unequal");
    Source src3(u, "Source");
    Queue q(u, "Q", SchedStrategy::FCFS);
    Sink snk3(u, "Sink");
    OpenClass ua(u, "A", 0);
    OpenClass ub(u, "B", 0);
    u.set_arrival(src3, ua, Exp(0.2));
    u.set_arrival(src3, ub, Exp(0.2));
    u.set_service(q, ua, Exp(1.0));
    u.set_service(q, ub, Exp(2.0));
    Routing Pu;
    Pu.set(ua, ua, src3, q, 1.0);
    Pu.set(ua, ua, q, snk3, 1.0);
    Pu.set(ub, ub, src3, q, 1.0);
    Pu.set(ub, ub, q, snk3, 1.0);
    u.link(Pu);
    try {
        run_snc(u.get_struct());
        std::printf("  unequal service rates at a shared station: NO REFUSAL (unexpected)\n");
    } catch (const std::exception& e) {
        std::printf("  unequal service rates at a shared station: %s\n", e.what());
    }
}

LINE_EXAMPLE("advanced/networkCalculus", snc_tandem_multiclass);

}  // namespace
}  // namespace examples
}  // namespace line
