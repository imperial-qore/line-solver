/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `matlab/examples/advanced/busyPeriod/`, `python/examples/advanced/busyPeriod/`:
 * the mean busy period of order n of a subnetwork.
 *
 * The busy period of order n for a set of stations is the time from the instant
 * a job entering the set finds n-1 jobs in it up to the next instant when fewer
 * than n remain. `pfqn_busyp` evaluates it exactly from the normalizing
 * constants of the subnetwork and of its complement (H. Daduna, "Busy Periods
 * for Subnetworks in Stochastic Networks: Mean Value Analysis", J. ACM 35(3),
 * 1988).
 *
 * The reference scripts print the same table twice, once exactly and once as
 * measured by LDES along a simulated sample path, and so does this port: the
 * client emits `--busyperiod` and reads the `busyPeriods` block back.
 */

#include <algorithm>
#include <cstdio>
#include <vector>

#include "example_util.h"
#include "examples_common.h"
#include "line/api/pfqn/pfqn_busyp.h"
#include "line/api/sn/sn_rt_stations.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/solvers/wrappers/ldes/solver_ldes.h"

namespace line {
namespace examples {

namespace {

/** Station-to-station chain routing: the class blocks weighted by class visits. */
Matrix<double> chain_routing(const Sn& sn) {
    const api::SnRtStations<double> rt = api::sn_rt_stations(sn);
    const std::size_t M = sn.nstations, K = sn.nclasses;
    Matrix<double> Pst(M, M, 0.0);
    for (std::size_t i = 0; i < M; ++i) {
        double vtot = 0.0;
        for (std::size_t r = 0; r < K; ++r) vtot += rt.Vst(i, r);
        for (std::size_t j = 0; j < M; ++j) {
            double flow = 0.0;
            for (std::size_t r = 0; r < K; ++r)
                for (std::size_t s = 0; s < K; ++s)
                    flow += rt.Vst(i, r) * rt.rtst(i * K + r, j * K + s);
            Pst(i, j) = (vtot > 0) ? flow / vtot : 0.0;
        }
    }
    return Pst;
}

void busyp_row(const char* label, const std::vector<double>& b) {
    std::printf("%-14s %10.4f %10.4f %10.4f\n", label, b[0], b[1], b[2]);
}

}  // namespace

/** `busyp_subnetwork`: busy periods of orders 1, 3 and 5 of a closed network. */
void busyp_subnetwork() {
    const double N = 5.0;
    Net m("busyPeriodModel");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    Queue q3(m, "Queue3", SchedStrategy::FCFS);
    ClosedClass c1(m, "Class1", N, q1, 0);

    const double rate[3] = {1.5, 0.9, 2.0};
    q1.set_service(c1, Exp(rate[0]));
    q2.set_service(c1, Exp(rate[1]));
    q3.set_service(c1, Exp(rate[2]));

    Routing P;
    P.set(c1, c1, q1, q2, 0.6);
    P.set(c1, c1, q1, q3, 0.4);
    P.set(c1, c1, q2, q1, 0.7);
    P.set(c1, c1, q2, q3, 0.3);
    P.set(c1, c1, q3, q1, 0.5);
    P.set(c1, c1, q3, q2, 0.5);
    m.link(P);
    const Sn& sn = m.get_struct();

    // alpha are the chain visit ratios; the closed formula is invariant to their
    // scale, so the reference's unnormalized visits serve as they are
    const mva::ChainDemands<double> d = mva::sn_get_demands_chain(sn);
    std::vector<double> alpha;
    for (std::size_t i = 0; i < sn.nstations; ++i) alpha.push_back(d.Vchain(i, 0));

    // mu(j,k) is a RATE: the load-dependent scaling over the chain service time
    Matrix<double> mu(sn.nstations, static_cast<std::size_t>(N), 0.0);
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t k = 0; k < static_cast<std::size_t>(N); ++k)
            mu(i, k) = 1.0 / d.STchain(i, 0);

    const Matrix<double> Pst = chain_routing(sn);
    const std::vector<std::size_t> orders = {1, 3, 5};
    const std::vector<std::vector<std::size_t>> subnets = {{0}, {1}, {2}, {0, 1}};

    section("Mean busy period of order n (exact, pfqn_busyp)");
    std::printf("%-14s %10s %10s %10s\n", "subnetwork", "n=1", "n=3", "n=5");
    const char* labels[4] = {"[1]", "[2]", "[3]", "[1 2]"};
    for (std::size_t s = 0; s < subnets.size(); ++s) {
        const pfqn::BusyPeriodResult r =
            pfqn::pfqn_busyp(alpha, mu, Pst, N, subnets[s], orders);
        busyp_row(labels[s], r.b);
    }

    // the same quantity measured along a simulated sample path
    ldes::LdesOptions opt;
    opt.samples = 2000000;
    opt.seed = 23000;
    opt.busy_period_orders = 5;
    opt.busy_period_subnets.push_back(std::vector<std::size_t>{0, 1});
    opt.verbose = false;
    const ldes::LdesResult sim = ldes::solver_ldes(sn, opt);

    section("Mean busy period of order n (measured, LDES)");
    std::printf("%-14s %10s %10s %10s\n", "subnetwork", "n=1", "n=3", "n=5");
    for (std::size_t s = 0; s < subnets.size(); ++s) {
        std::vector<double> b(orders.size(), 0.0);
        for (std::size_t t = 0; t < sim.busy_periods.size(); ++t) {
            const ldes::LdesResult::BusyPeriodTarget& tgt = sim.busy_periods[t];
            if (tgt.job_class != -1) continue;
            std::vector<std::size_t> st = tgt.stations;
            std::sort(st.begin(), st.end());
            if (st != subnets[s]) continue;
            for (std::size_t k = 0; k < orders.size(); ++k)
                b[k] = (orders[k] <= tgt.mean.size()) ? tgt.mean[orders[k] - 1] : 0.0;
            break;
        }
        busyp_row(labels[s], b);
    }
}

LINE_EXAMPLE("advanced/busyPeriod", busyp_subnetwork);

}  // namespace examples
}  // namespace line
