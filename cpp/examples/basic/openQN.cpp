/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/basic/openQN`: the open-network half of the basic gallery.
 *
 * Every entry is its reference script's model and its `__main__` block, in the
 * same order and with the same solvers. Three of them ask for something this
 * port does not carry and say so BY NAME rather than answering with a
 * neighbour: JMT and LDES, which have no C++ engine at all, `getAvgNodeTable`,
 * which has no C++ counterpart, and `oqn_cs_routing`'s per-(station, class)
 * CTMC cutoff, which `CtmcOptions` cannot express -- it carries a scalar and a
 * per-class vector, so a matrix would have to be flattened, and flattening it
 * changes which states are enumerated.
 *
 * NONE OF THESE SCRIPTS DRAWS A RANDOM PARAMETER, so nothing here is a
 * materialized sample; the two seeds that appear (`SSA`, `JMT`) are solver
 * seeds and travel with the solver call.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "example_util.h"
#include "examples_common.h"
#include "line/api/mam/map_moment.h"
#include "line/lang/dist_fitters.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/fluid/fluid_kp.h"
#include "line/solvers/fluid/fluid_runner.h"
#include "line/solvers/mam/solver_mam_runner.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/nc/solver_nc_runner.h"
#include "line/solvers/ssa/ssa_dispatch.h"

namespace line {
namespace examples {

namespace {

/** Linear interpolation of a sampled series, `numpy.interp` on a sorted grid. */
double interp_at(double x, const std::vector<double>& t, const std::vector<double>& y) {
    if (t.empty()) return 0.0;
    if (x <= t.front()) return y.front();
    if (x >= t.back()) return y.back();
    const std::size_t j =
        static_cast<std::size_t>(std::lower_bound(t.begin(), t.end(), x) - t.begin());
    const double t0 = t[j - 1], t1 = t[j];
    if (!(t1 > t0)) return y[j];
    return y[j - 1] + (y[j] - y[j - 1]) * (x - t0) / (t1 - t0);
}

/** The segment of a cyclic piecewise-constant schedule that covers `t`. */
std::size_t schedule_segment(double t, const std::vector<double>& breakpoints, bool cyclic) {
    const double t0 = breakpoints.front(), period = breakpoints.back() - breakpoints.front();
    double u = t;
    if (cyclic && period > 0.0) u = t0 + std::fmod(std::fmod(t - t0, period) + period, period);
    for (std::size_t k = 0; k + 1 < breakpoints.size(); ++k)
        if (u >= breakpoints[k] && u < breakpoints[k + 1]) return k;
    return breakpoints.size() - 2;
}

}  // namespace

// ---------------------------------------------------------------------------
// oqn_basic
// ---------------------------------------------------------------------------

/**
 * Source -> Delay -> Queue1 -> Sink, one open class, solved by every native
 * solver the reference runs.
 *
 * SSA carries the reference's own enlarged budget: at 5e3 events the standard
 * error of Util and Tput on this rho=0.1 queue is about 4%, which is the size
 * of the cross-codebase parity tolerance, and 5e4 brings it to about 1.4%.
 */
void oqn_basic() {
    Net m("model");
    Delay delay(m, "Delay");
    Queue queue(m, "Queue1", SchedStrategy::FCFS);
    Source source(m, "Source");
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1", 0);

    delay.set_service(c1, HyperExp(0.5, 3.0, 10.0));
    queue.set_service(c1, Exp(1.0));
    source.set_arrival(c1, Exp(0.1));

    Routing P;
    P.set(c1, c1, delay, queue, 1.0);
    P.set(c1, c1, queue, sink, 1.0);
    P.set(c1, c1, source, delay, 1.0);
    m.link(P);
    const Sn& sn = m.get_struct();

    section("CTMC");
    ctmc::CtmcOptions copt;
    copt.cutoff = 10.0;
    print_avg(sn, ctmc::solver_ctmc_run_analyzer(sn, copt));

    section("FLD");
    print_avg_sim(sn, fluid::solver_fluid_run_analyzer(sn, fluid::FluidOptions()));

    section("MVA");
    Matrix<double> init;
    print_avg(sn, mva::solver_mva_run_analyzer(sn, mva::MvaOptions(), init));

    section("MAM");
    print_avg(sn, mam::solver_mam_run_analyzer(sn, mam::MamOptions()));

    section("NC");
    print_avg(sn, nc::solver_nc_run_analyzer(sn, nc::NcSolverOptions()));

    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));

    section("SSA");
    ssa::SsaOptions sopt;
    sopt.seed = 23000;
    sopt.samples = 50000;
    print_avg_sim(sn, ssa::solver_ssa(sn, sopt));

    section("LDES");
    print_avg(m.get_struct(), ldes_avg(m, sim_opts(23000, 200000)));
}

// ---------------------------------------------------------------------------
// oqn_cs_routing
// ---------------------------------------------------------------------------

/**
 * Three open classes over two PS queues and a ClassSwitch node: A and B are
 * both switched to C on the way from Queue 1 to Queue 2.
 *
 * The ClassSwitch node is created AFTER the classes, unlike the reference,
 * because `add_class_switch` takes its (nclasses x nclasses) matrix at
 * construction; the node order changes, the station order does not, so the
 * printed table is the reference's.
 */
void oqn_cs_routing() {
    Net m("myModel");
    Source source(m, "Source 1");
    Queue q1(m, "Queue 1", SchedStrategy::PS);
    Sink sink(m, "Sink 1");
    Queue q2(m, "Queue 2", SchedStrategy::PS);

    OpenClass ca(m, "Class A", 0);
    OpenClass cb(m, "Class B", 0);
    OpenClass cc(m, "Class C", 0);

    Matrix<double> C(3, 3, 0.0);
    for (std::size_t i = 0; i < 3; ++i) C(i, i) = 1.0;
    const std::size_t cs = m.add_class_switch("ClassSwitch 1", C);

    source.set_arrival(ca, D::exp_mean(0.5));
    source.set_arrival(cb, D::exp_mean(1.0));
    source.set_arrival(cc, Disabled());

    q1.set_service(ca, D::exp_mean(0.2));
    q1.set_service(cb, D::exp_mean(0.3));
    q1.set_service(cc, D::exp_mean(0.333333));

    q2.set_service(ca, D::exp_mean(1.0));
    q2.set_service(cb, D::exp_mean(1.0));
    q2.set_service(cc, D::exp_mean(0.15));

    Routing P;
    P.set(ca, ca, source, q1, 1.0);
    P.set(ca, ca, q1, cs, 1.0);
    P.set(ca, cc, cs, q2, 1.0);
    P.set(ca, ca, q2, sink, 1.0);

    P.set(cb, cb, source, q1, 1.0);
    P.set(cb, cb, q1, cs, 1.0);
    P.set(cb, cc, cs, q2, 1.0);
    P.set(cb, cb, q2, sink, 1.0);

    P.set(cc, cc, source, q1, 1.0);
    P.set(cc, cc, q1, cs, 1.0);
    P.set(cc, cc, cs, q2, 1.0);
    P.set(cc, cc, q2, sink, 1.0);
    m.link(P);
    const Sn& sn = m.get_struct();

    na("CTMC",
       "the reference passes options.cutoff as a (station x class) matrix; CtmcOptions carries a "
       "scalar and a per-class vector only, and collapsing the matrix would enumerate a different "
       "state space");

    section("FLD");
    print_avg_sim(sn, fluid::solver_fluid_run_analyzer(sn, fluid::FluidOptions()));

    section("MVA");
    Matrix<double> init;
    print_avg(sn, mva::solver_mva_run_analyzer(sn, mva::MvaOptions(), init));

    section("MAM");
    print_avg(sn, mam::solver_mam_run_analyzer(sn, mam::MamOptions()));

    section("NC");
    print_avg(sn, nc::solver_nc_run_analyzer(sn, nc::NcSolverOptions()));

    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 100000)));

    section("SSA");
    ssa::SsaOptions sopt;
    sopt.seed = 23000;
    sopt.samples = 100000;
    print_avg_sim(sn, ssa::solver_ssa(sn, sopt));

    // TODO(cpp): solver = np.append(solver, LDES(model, keep=True, verbose=True, seed=23000,
    // samples=100000))
    na("LDES", "SolverLDES(model, seed=23000, samples=100000) has no C++ engine in this port");
}

// ---------------------------------------------------------------------------
// oqn_fourqueues
// ---------------------------------------------------------------------------

/** Three open classes over four queues, each fed back into Queue1 with p=1/4. */
void oqn_fourqueues() {
    Net m("MyNetwork");
    Source source(m, "Source");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    Queue q3(m, "Queue3", SchedStrategy::PS);
    Queue q4(m, "Queue4", SchedStrategy::FCFS);
    Sink sink(m, "Sink");

    const std::vector<std::size_t> cls = {m.add_open_class("Class1", 0),
                                          m.add_open_class("Class2", 0),
                                          m.add_open_class("Class3", 0)};

    const double arrival[3] = {5.0, 8.0, 7.0};
    const double svc[3][4] = {{0.3, 1.1, 2.0, 1.5}, {0.5, 1.3, 2.1, 0.9}, {0.6, 1.5, 1.9, 2.3}};
    const std::size_t queues[4] = {q1, q2, q3, q4};
    for (std::size_t r = 0; r < 3; ++r) {
        source.set_arrival(cls[r], D::exp_mean(arrival[r]));
        for (std::size_t i = 0; i < 4; ++i) m.set_service(queues[i], cls[r], D::exp_mean(svc[r][i]));
    }

    Routing P;
    for (std::size_t r = 0; r < 3; ++r) {
        P.set(cls[r], cls[r], source, q1, 1.0);
        P.set(cls[r], cls[r], q1, q2, 0.25);
        P.set(cls[r], cls[r], q1, q3, 0.25);
        P.set(cls[r], cls[r], q1, q4, 0.25);
        P.set(cls[r], cls[r], q1, sink, 0.25);
        P.set(cls[r], cls[r], q2, q1, 1.0);
        P.set(cls[r], cls[r], q3, q1, 1.0);
        P.set(cls[r], cls[r], q4, q1, 1.0);
    }
    m.link(P);
    const Sn& sn = m.get_struct();

    section("CTMC");
    ctmc::CtmcOptions copt;
    copt.cutoff = 1.0;
    print_avg(sn, ctmc::solver_ctmc_run_analyzer(sn, copt));

    section("MVA");
    Matrix<double> init;
    print_avg(sn, mva::solver_mva_run_analyzer(sn, mva::MvaOptions(), init));

    section("MAM");
    print_avg(sn, mam::solver_mam_run_analyzer(sn, mam::MamOptions()));

    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 1000000)));
    // TODO(cpp): solver = np.append(solver, LDES(model, seed=23000, samples=1000000))
    na("LDES", "SolverLDES(model, seed=23000, samples=1000000) has no C++ engine in this port");
}

// ---------------------------------------------------------------------------
// oqn_mapt
// ---------------------------------------------------------------------------

namespace {

/** The two MAPt segments of `oqn_mapt`, as the reference writes them. */
const double kMaptBreak[3] = {0.0, 1.0, 2.5};

mam::Map<double> mapt_segment(std::size_t k) {
    mam::Map<double> s;
    s.D0 = Matrix<double>(2, 2, 0.0);
    s.D1 = Matrix<double>(2, 2, 0.0);
    if (k == 0) {
        s.D0(0, 0) = -5.0;  s.D0(0, 1) = 1.0;
        s.D0(1, 0) = 2.0;   s.D0(1, 1) = -4.0;
        s.D1(0, 0) = 3.0;   s.D1(0, 1) = 1.0;
        s.D1(1, 0) = 1.0;   s.D1(1, 1) = 1.0;
    } else {
        s.D0(0, 0) = -12.0; s.D0(0, 1) = 3.0;
        s.D0(1, 0) = 5.0;   s.D0(1, 1) = -9.0;
        s.D1(0, 0) = 7.0;   s.D1(0, 1) = 2.0;
        s.D1(1, 0) = 2.0;   s.D1(1, 1) = 2.0;
    }
    return s;
}

/**
 * `MAPt.getRateAt(t)`: the stationary arrival rate of the segment covering `t`.
 *
 * A MAPt has no scalar intensity of its own -- only the pair (D0, D1) of the
 * segment does -- so the rate the reference prints is `map_lambda` of that pair.
 */
double mapt_rate_at(double t) {
    const std::vector<double> bp(kMaptBreak, kMaptBreak + 3);
    return mam::map_lambda(mapt_segment(schedule_segment(t, bp, true)));
}

/**
 * The (offset, width) of a (station, class) block in the `kp` state vector.
 *
 * `FluidKpTransient` carries the covariance and the raw state, not a per-pair
 * mean, so the mean trajectory the reference reads out of `result.QNt` is
 * summed off the state here. The layout is `solver_fluid_kp_core`'s: the
 * enabled (station, class) pairs in station-major order, each as wide as its
 * D0, and a pair that is disabled or carries no positive rate takes no room.
 */
std::pair<std::size_t, std::size_t> kp_block_of(const Sn& sn, std::size_t ist, std::size_t cls) {
    std::size_t off = 0;
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            if (sn.disabled[i][r]) continue;
            const std::size_t h = sn.service[i][r].D0.rows();
            const double rate = sn.rates(i, r);
            if (h == 0 || !std::isfinite(rate) || rate <= 0.0) continue;
            if (i == ist && r == cls) return std::make_pair(off, h);
            off += h;
        }
    return std::make_pair(off, static_cast<std::size_t>(0));
}

}  // namespace

/**
 * A MAP_t arrival feeding an infinite server, under the Ko-Pender limits.
 *
 * `SolverFLD`'s 'kp' method integrates the fluid AND diffusion limits of Ko and
 * Pender (Oper. Res. Lett. 45, 2017), so it is the only fluid method that
 * returns a second moment; at an infinite-server station the rate functions are
 * affine in the state, so the variance is exact rather than asymptotic. The
 * tolerance is the reference's 1e-9 and not the default 1e-4, at which the
 * integrators of the different codebases disagree by about 0.2%.
 */
void oqn_mapt() {
    Net m("model");
    Source source(m, "Source");
    Delay delay(m, "Delay");
    Sink sink(m, "Sink");
    OpenClass oc(m, "OpenClass", 0);

    std::vector<double> bp(kMaptBreak, kMaptBreak + 3);
    std::vector<Matrix<double>> D0segs, D1segs;
    for (std::size_t k = 0; k < 2; ++k) {
        const mam::Map<double> s = mapt_segment(k);
        D0segs.push_back(s.D0);
        D1segs.push_back(s.D1);
    }
    source.set_arrival(oc, D::mapt(bp, D0segs, D1segs, true));
    delay.set_service(oc, Exp(2.0));

    Routing P;
    serial(P, {source, delay, sink});
    m.link(P);
    const Sn& sn = m.get_struct();

    fluid::FluidOptions opt;
    opt.method = "kp";
    opt.tol = 1e-9;
    // The method is presentation: the golden keys this table plain `FLD`.
    section("FLD (kp)", "FLD");
    print_avg_sim(sn, fluid::solver_fluid_run_analyzer(sn, opt));

    // Transient mean AND variance over two periods of the cyclic schedule.
    fluid::FluidOptions topt = opt;
    topt.timespan_end = 5.0;
    const fluid::FluidKpTransient tran = fluid::solver_fluid_tran_avg_var(sn, topt);
    const std::size_t ist = m.station_index(delay) - 1;
    const std::pair<std::size_t, std::size_t> blk = kp_block_of(sn, ist, oc - 1);

    std::vector<double> mean(tran.t.size(), 0.0), var(tran.t.size(), 0.0);
    for (std::size_t n = 0; n < tran.t.size(); ++n) {
        for (std::size_t p = 0; p < blk.second; ++p) mean[n] += tran.q[n][blk.first + p];
        var[n] = tran.QVar[n](ist, oc - 1);
    }

    std::printf("\n   t     lambda(t)    QLen mean     QLen var\n");
    const double probe[7] = {0.5, 1.0, 1.5, 2.5, 3.0, 4.0, 5.0};
    for (std::size_t i = 0; i < 7; ++i)
        std::printf("%6.2f  %10.4f  %11.6f  %11.6f\n", probe[i], mapt_rate_at(probe[i]),
                    interp_at(probe[i], tran.t, mean), interp_at(probe[i], tran.t, var));
}

// ---------------------------------------------------------------------------
// oqn_multichain_cs
// ---------------------------------------------------------------------------

/**
 * Five open classes in three chains: two pairs that switch into each other at
 * the Router, and one class that does not.
 */
void oqn_multichain_cs() {
    Net m("oqn_multichain");
    Source source(m, "RequestSource");
    Queue router(m, "Router", SchedStrategy::FCFS);
    Queue cache(m, "WebCache", SchedStrategy::PS);
    Queue app(m, "AppServer", SchedStrategy::FCFS);
    Queue data(m, "DataServer", SchedStrategy::PS);
    Sink sink(m, "ResponseSink");

    const std::vector<std::size_t> cls = {
        m.add_open_class("InteractiveA", 0), m.add_open_class("InteractiveB", 0),
        m.add_open_class("BatchA", 0), m.add_open_class("BatchB", 0),
        m.add_open_class("RealTime", 0)};

    const std::size_t st[4] = {router, cache, app, data};
    const double svc[4][5] = {{0.15, 0.18, 0.2, 0.22, 0.1},
                              {0.3, 0.35, 1.0, 1.2, 0.2},
                              {0.5, 0.6, 1.5, 1.8, 0.3},
                              {0.4, 0.45, 0.8, 1.0, 0.25}};
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t r = 0; r < 5; ++r) m.set_service(st[i], cls[r], D::exp_mean(svc[i][r]));

    const double lambda[5] = {0.3, 0.1, 0.15, 0.12, 0.2};
    for (std::size_t r = 0; r < 5; ++r) source.set_arrival(cls[r], Exp(lambda[r]));

    Routing P;
    for (std::size_t r = 0; r < 5; ++r) {
        P.set(cls[r], cls[r], source, router, 1.0);
        P.set(cls[r], cls[r], cache, app, 1.0);
        P.set(cls[r], cls[r], app, data, 1.0);
        P.set(cls[r], cls[r], data, sink, 1.0);
    }
    // The class switch happens on the Router -> WebCache edge alone.
    P.set(cls[0], cls[0], router, cache, 0.6);
    P.set(cls[0], cls[1], router, cache, 0.4);
    P.set(cls[1], cls[1], router, cache, 0.7);
    P.set(cls[1], cls[0], router, cache, 0.3);
    P.set(cls[2], cls[2], router, cache, 0.5);
    P.set(cls[2], cls[3], router, cache, 0.5);
    P.set(cls[3], cls[3], router, cache, 0.6);
    P.set(cls[3], cls[2], router, cache, 0.4);
    P.set(cls[4], cls[4], router, cache, 1.0);
    m.link(P);
    const Sn& sn = m.get_struct();

    note("Open QN with 3 chains, 5 classes, and class switching:");
    note("- Chain 1: InteractiveA and InteractiveB (switch at Router)");
    note("- Chain 2: BatchA and BatchB (switch at Router)");
    note("- Chain 3: RealTime");

    section("MVA");
    Matrix<double> init;
    print_avg(sn, mva::solver_mva_run_analyzer(sn, mva::MvaOptions(), init));

    section("FLD");
    print_avg_sim(sn, fluid::solver_fluid_run_analyzer(sn, fluid::FluidOptions()));

    section("SSA");
    ssa::SsaOptions sopt;
    sopt.seed = 23000;
    print_avg_sim(sn, ssa::solver_ssa(sn, sopt));
}

// ---------------------------------------------------------------------------
// oqn_nhpp
// ---------------------------------------------------------------------------

/**
 * An M(t)/M/1 whose arrival intensity is a cyclic step schedule.
 *
 * Steady state is the time-average rate, so the interesting answer is the
 * transient: `getTranAvg` integrates the closing drift with the intensity as a
 * time-varying multiplier, and the queue throughput tracks lambda(t).
 */
void oqn_nhpp() {
    Net m("model");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass oc(m, "OpenClass", 0);

    // Rates 2, 8, 4 held for 3, 1, 2 time units, repeating cyclically.
    const std::vector<double> bp = {0.0, 3.0, 4.0, 6.0};
    const std::vector<double> rates = {2.0, 8.0, 4.0};
    source.set_arrival(oc, D::nhpp(bp, rates, true));
    queue.set_service(oc, Exp(10.0));

    Routing P;
    serial(P, {source, queue, sink});
    m.link(P);
    const Sn& sn = m.get_struct();

    // TODO(cpp): avg_table_1 = LDES(model, seed=1234, samples=100000).avg_table(); print(...)
    na("LDES", "SolverLDES(model, seed=1234, samples=100000) has no C++ engine in this port");

    fluid::FluidOptions opt;
    opt.timespan_end = 12.0;
    const std::vector<fluid::FluidTranPoint> tran = fluid::solver_fluid_tran_avg(sn, opt);
    const std::size_t ist = m.station_index(queue) - 1;
    std::vector<double> t, y;
    for (std::size_t n = 0; n < tran.size(); ++n) {
        t.push_back(tran[n].t);
        y.push_back(tran[n].TN(ist, oc - 1));
    }

    note("FLD transient:");
    const double probe[6] = {1.5, 3.5, 5.0, 7.5, 9.5, 11.0};
    for (std::size_t i = 0; i < 6; ++i)
        std::printf("  t=%5.2f lambda(t)=%6.3f queue Tput=%6.3f\n", probe[i],
                    rates[schedule_segment(probe[i], bp, true)], interp_at(probe[i], t, y));
}

// ---------------------------------------------------------------------------
// oqn_oneline
// ---------------------------------------------------------------------------

/**
 * `Network.tandemPsInf(lambda, D, Z)`, written out.
 *
 * The factory is not ported, so the model it builds is built here: a Source, one
 * Delay per row of Z, one PS queue per row of D -- all named `Station<i>` as the
 * factory names them -- and serial routing into the Sink. The arrival rate of
 * class r is `lambda(r)` and the mean service time at station i is `D(i,r)`.
 */
void oqn_oneline() {
    note("This example shows a compact solution of a tandem open queueing network.");

    const double lambda[2] = {1.0 / 50.0, 2.0 / 50.0};
    const double Z[1][2] = {{91.0, 92.0}};
    const double Dm[2][2] = {{10.0, 5.0}, {5.0, 9.0}};

    Net m("Model");
    std::vector<std::size_t> line;
    line.push_back(m.add_source("Source"));
    line.push_back(m.add_delay("Station1"));
    line.push_back(m.add_queue("Station2", SchedStrategy::PS));
    line.push_back(m.add_queue("Station3", SchedStrategy::PS));
    line.push_back(m.add_sink("Sink"));

    const std::vector<std::size_t> cls = {m.add_open_class("Class1"), m.add_open_class("Class2")};
    for (std::size_t r = 0; r < 2; ++r) {
        m.set_arrival(line[0], cls[r], D::exp_mean(1.0 / lambda[r]));
        m.set_service(line[1], cls[r], D::exp_mean(Z[0][r]));
        for (std::size_t i = 0; i < 2; ++i) m.set_service(line[i + 2], cls[r], D::exp_mean(Dm[i][r]));
    }

    Routing P;
    for (std::size_t r = 0; r < 2; ++r) serial(P, cls[r], line);
    m.link(P);
    const Sn& sn = m.get_struct();

    section("MVA");
    Matrix<double> init;
    print_avg(sn, mva::solver_mva_run_analyzer(sn, mva::MvaOptions(), init));
}

// ---------------------------------------------------------------------------
// oqn_trace_driven
// ---------------------------------------------------------------------------

/**
 * Trace-driven service: the queue replays the reference's own sample file, so
 * both codebases serve the same 10000 samples in the same order.
 *
 * BOTH SOLVERS RUN. The Replayer is built through `replayer_from_file`, which
 * keeps the trace PATH beside the samples -- JMT is handed a `ReplayerPar`
 * naming a file and cannot take a trace inline, so a Replayer built from bare
 * samples was refused by the exporter. LDES has a native C++ engine
 * (`common/ldes`), which the old refusal here predated.
 */
void oqn_trace_driven() {
    Net m("model");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass oc(m, "OpenClass", 0);

    source.set_arrival(oc, Exp(1.0));
    queue.set_service(oc,
                  replayer_from_file(std::string(LINE_EXAMPLES_REPO_ROOT) +
                                     "/matlab/examples/basic/openQN/example_trace.txt"));

    Routing P;
    serial(P, {source, queue, sink});
    m.link(P);
    m.get_struct();

    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    section("LDES");
    print_avg(m.get_struct(), ldes_avg(m, sim_opts(23000)));
}

// ---------------------------------------------------------------------------
// oqn_vsinks
// ---------------------------------------------------------------------------

/**
 * Two virtual sinks: Router nodes that split the departure stream so that each
 * exit path carries its own throughput.
 *
 * The point of the reference is the NODE table -- a Router is not a station, so
 * it has no row in the station table -- and this port has no `getAvgNodeTable`,
 * so that half is refused by name.
 */
void oqn_vsinks() {
    Net m("model");
    Source source(m, "Source");
    Queue queue(m, "Queue1", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    Router vsink1(m, "VSink1");
    Router vsink2(m, "VSink2");

    OpenClass c1(m, "Class1");
    OpenClass c2(m, "Class2");

    source.set_arrival(c1, Exp(1.0));
    queue.set_service(c1, Exp(100.0));
    source.set_arrival(c2, Exp(1.0));
    queue.set_service(c2, Exp(100.0));

    Routing P;
    P.set(c1, c1, source, queue, 1.0);
    P.set(c1, c1, queue, vsink1, 0.6);
    P.set(c1, c1, queue, vsink2, 0.4);
    P.set(c1, c1, vsink1, sink, 1.0);
    P.set(c1, c1, vsink2, sink, 1.0);

    P.set(c2, c2, source, queue, 1.0);
    P.set(c2, c2, queue, vsink1, 0.1);
    P.set(c2, c2, queue, vsink2, 0.9);
    P.set(c2, c2, vsink1, sink, 1.0);
    P.set(c2, c2, vsink2, sink, 1.0);
    m.link(P);
    const Sn& sn = m.get_struct();

    section("MVA");
    Matrix<double> init;
    print_avg(sn, mva::solver_mva_run_analyzer(sn, mva::MvaOptions(), init));
    na("MVA getAvgNodeTable",
       "the node-level table, which is where the virtual sinks appear, has no C++ counterpart");

    section("MAM");
    print_avg(sn, mam::solver_mam_run_analyzer(sn, mam::MamOptions()));
    na("MAM getAvgNodeTable", "the node-level table has no C++ counterpart");

    section("NC");
    print_avg(sn, nc::solver_nc_run_analyzer(sn, nc::NcSolverOptions()));
    na("NC getAvgNodeTable", "the node-level table has no C++ counterpart");
}

LINE_EXAMPLE("basic/openQN", oqn_basic);
LINE_EXAMPLE("basic/openQN", oqn_cs_routing);
LINE_EXAMPLE("basic/openQN", oqn_fourqueues);
LINE_EXAMPLE("basic/openQN", oqn_mapt);
LINE_EXAMPLE("basic/openQN", oqn_multichain_cs);
LINE_EXAMPLE("basic/openQN", oqn_nhpp);
LINE_EXAMPLE("basic/openQN", oqn_oneline);
LINE_EXAMPLE("basic/openQN", oqn_trace_driven);
LINE_EXAMPLE("basic/openQN", oqn_vsinks);

}  // namespace examples
}  // namespace line
