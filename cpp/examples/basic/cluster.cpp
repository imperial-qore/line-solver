/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/basic/cluster/`: server farms behind a dispatcher.
 *
 * The reference scripts build their models with `Network.cluster*` and the
 * `Cluster` builder of `python/line_solver/gen/cluster.py`, neither of which
 * this port carries. The three topologies those factories produce are
 * transcribed below -- Source -> Dispatcher -> Servers -> Sink for the open
 * farm, Think -> Dispatcher -> Servers -> Think for the closed one, and both
 * arcs at once for the mixed one -- so the models are the reference's and only
 * the spelling differs. `Cluster.build()` converts service RATES to mean
 * service times, which is why every caller below passes 1/mu.
 *
 * WHAT IS REFUSED. `dispatch_mixed` cross-checks MVA with SolverLDES, which has
 * no C++ port. `cl_compare` sweeps the dispatching policy over RAND and RROBIN
 * under SolverSSA: the NRM sub-engine refuses RROBIN because it carries no
 * round-robin pointer, and the dispatch falls through to the serial simulator,
 * which resolves the pointer at firing time -- so the round-robin arm is a
 * round-robin answer and not a random-routing one under that name.
 */

#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

#include "example_util.h"
#include "examples_common.h"

namespace line {
namespace examples {

namespace {

/** A cluster's per-(server, class) mean service times. */
typedef std::vector<std::vector<double> > Demands;

/** `Demands` of M identical single-class servers, each of mean service time d. */
Demands uniform_demands(std::size_t M, double d) {
    return Demands(M, std::vector<double>(1, d));
}

/**
 * `Network.cluster`: Source -> Dispatcher -> M parallel servers -> Sink.
 *
 * The dispatcher is a Router whose strategy spreads the arrivals over the
 * servers it is linked to; `d[i][r]` is the mean service time of class r there.
 */
Net cluster_open(const std::vector<double>& lambda, const Demands& d,
                 const std::vector<SchedStrategy>& strategies, const std::vector<double>& S,
                 RoutingStrategy dispatching, double service_scv = 1.0) {
    const std::size_t M = d.size(), R = lambda.size();
    Net model("Cluster");
    Source source(model, "Source");
    Router dispatcher(model, "Dispatcher");
    std::vector<std::size_t> servers;
    for (std::size_t i = 0; i < M; ++i) {
        servers.push_back(model.add_queue("Station" + std::to_string(i + 1), strategies[i]));
        if (S[i] > 1.0) model.set_number_of_servers(servers[i], S[i]);
    }
    Sink sink(model, "Sink");

    Routing P;
    for (std::size_t r = 0; r < R; ++r) {
        OpenClass cls(model, "Class" + std::to_string(r + 1), 0);
        source.set_arrival(cls, D::exp_mean(1.0 / lambda[r]));
        // Cluster.applyDistributionScvs fits an APH of the same mean once the
        // service SCV departs from one, and leaves the exponential otherwise.
        for (std::size_t i = 0; i < M; ++i)
            model.set_service(servers[i], cls, service_scv == 1.0
                                                   ? D::exp_mean(d[i][r])
                                                   : aph_fit(d[i][r], service_scv));
        P.set(cls, cls, source, dispatcher, 1.0);
        for (std::size_t i = 0; i < M; ++i) {
            P.set(cls, cls, dispatcher, servers[i], 1.0);
            P.set(cls, cls, servers[i], sink, 1.0);
        }
    }
    model.link(P);
    for (std::size_t r = 1; r <= R; ++r) model.set_routing(dispatcher, r, dispatching);
    return model;
}

/** `Network.cluster_closed`: Think -> Dispatcher -> M servers -> Think. */
Net cluster_closed(const std::vector<double>& N, const std::vector<double>& Z, const Demands& d,
                   const std::vector<SchedStrategy>& strategies, const std::vector<double>& S,
                   RoutingStrategy dispatching) {
    const std::size_t M = d.size(), R = N.size();
    Net model("Cluster");
    Delay think(model, "Think");
    Router dispatcher(model, "Dispatcher");
    std::vector<std::size_t> servers;
    for (std::size_t i = 0; i < M; ++i) {
        servers.push_back(model.add_queue("Station" + std::to_string(i + 1), strategies[i]));
        if (S[i] > 1.0) model.set_number_of_servers(servers[i], S[i]);
    }

    Routing P;
    for (std::size_t r = 0; r < R; ++r) {
        ClosedClass cls(model, "Class" + std::to_string(r + 1), N[r], think, 0);
        think.set_service(cls, D::exp_mean(Z[r]));
        for (std::size_t i = 0; i < M; ++i) model.set_service(servers[i], cls, D::exp_mean(d[i][r]));
        P.set(cls, cls, think, dispatcher, 1.0);
        for (std::size_t i = 0; i < M; ++i) {
            P.set(cls, cls, dispatcher, servers[i], 1.0);
            P.set(cls, cls, servers[i], think, 1.0);
        }
    }
    model.link(P);
    for (std::size_t r = 1; r <= R; ++r) model.set_routing(dispatcher, r, dispatching);
    return model;
}

/**
 * `Network.cluster_mixed`: the open and the closed farm sharing one dispatcher
 * and one set of servers.
 *
 * The classes are ordered open first, so column r of `d` is the open class r
 * while r < Ro and the closed class r - Ro after that. The shared arcs carry
 * per-class probabilities: a server feeds the Sink for an open class and the
 * Think delay for a closed one.
 */
Net cluster_mixed(const std::vector<double>& lambda, const std::vector<double>& N,
                  const std::vector<double>& Z, const Demands& d,
                  const std::vector<SchedStrategy>& strategies, const std::vector<double>& S,
                  RoutingStrategy dispatching) {
    const std::size_t M = d.size(), Ro = lambda.size(), Rc = N.size();
    Net model("Cluster");
    Source source(model, "Source");
    Delay think(model, "Think");
    Router dispatcher(model, "Dispatcher");
    std::vector<std::size_t> servers;
    for (std::size_t i = 0; i < M; ++i) {
        servers.push_back(model.add_queue("Station" + std::to_string(i + 1), strategies[i]));
        if (S[i] > 1.0) model.set_number_of_servers(servers[i], S[i]);
    }
    Sink sink(model, "Sink");

    Routing P;
    std::vector<std::size_t> jobclasses;
    for (std::size_t r = 0; r < Ro; ++r) {
        OpenClass cls(model, "Class" + std::to_string(r + 1), 0);
        jobclasses.push_back(cls);
        source.set_arrival(cls, D::exp_mean(1.0 / lambda[r]));
        think.set_service(cls, Disabled());
        for (std::size_t i = 0; i < M; ++i) model.set_service(servers[i], cls, D::exp_mean(d[i][r]));
        P.set(cls, cls, source, dispatcher, 1.0);
        for (std::size_t i = 0; i < M; ++i) {
            P.set(cls, cls, dispatcher, servers[i], 1.0);
            P.set(cls, cls, servers[i], sink, 1.0);
        }
    }
    for (std::size_t c = 0; c < Rc; ++c) {
        const std::size_t r = Ro + c;
        ClosedClass cls(model, "Class" + std::to_string(r + 1), N[c], think, 0);
        jobclasses.push_back(cls);
        think.set_service(cls, D::exp_mean(Z[c]));
        for (std::size_t i = 0; i < M; ++i) model.set_service(servers[i], cls, D::exp_mean(d[i][r]));
        P.set(cls, cls, think, dispatcher, 1.0);
        for (std::size_t i = 0; i < M; ++i) {
            P.set(cls, cls, dispatcher, servers[i], 1.0);
            P.set(cls, cls, servers[i], think, 1.0);
        }
    }
    model.link(P);
    for (std::size_t r = 0; r < jobclasses.size(); ++r)
        model.set_routing(dispatcher, jobclasses[r], dispatching);
    return model;
}

/** `SchedStrategy` repeated across M servers, which is all the examples ask for. */
std::vector<SchedStrategy> same_sched(std::size_t M, SchedStrategy s) {
    return std::vector<SchedStrategy>(M, s);
}

/** Single-server multiplicity at every station, the `Cluster` default. */
std::vector<double> single_servers(std::size_t M) { return std::vector<double>(M, 1.0); }

}  // namespace

/** An open farm: Poisson arrivals at rate 0.4 spread over four PS servers. */
void cl_basic() {
    const std::vector<double> lambda(1, 0.4);
    const Demands d = uniform_demands(4, 1.0);
    Net model = cluster_open(lambda, d, same_sched(4, SchedStrategy::PS), single_servers(4),
                             RoutingStrategy::RAND);
    print_avg(solve_avg("MVA", model));
}

/** A closed farm: eight jobs cycling between a Think delay and three PS servers. */
void cl_closed() {
    const std::vector<double> N(1, 8.0), Z(1, 1.0);
    const Demands d = uniform_demands(3, 1.0);
    Net model = cluster_closed(N, Z, d, same_sched(3, SchedStrategy::PS), single_servers(3),
                               RoutingStrategy::RAND);
    print_avg(solve_avg("MVA", model));
}

/** Two open classes over the same two servers: interactive against batch. */
void cl_multiclass() {
    std::vector<double> lambda;
    lambda.push_back(0.3);
    lambda.push_back(0.2);
    Demands d(2, std::vector<double>());
    for (std::size_t i = 0; i < 2; ++i) {
        d[i].push_back(1.0);
        d[i].push_back(0.5);
    }
    Net model = cluster_open(lambda, d, same_sched(2, SchedStrategy::PS), single_servers(2),
                             RoutingStrategy::RAND);
    print_avg(solve_avg("MVA", model));
}

/**
 * The same four-server farm under two dispatching policies, simulated.
 *
 * `Cluster.set_service_rate(0.4)` is a RATE, so the mean service time is 2.5.
 */
void cl_compare() {
    const std::vector<double> lambda(1, 1.0);
    const Demands d = uniform_demands(4, 1.0 / 0.4);

    SolverOpts ssa;
    ssa.seed = 23000;
    ssa.samples = 2000;

    const RoutingStrategy policies[2] = {RoutingStrategy::RAND, RoutingStrategy::RROBIN};
    const char* policy_name[2] = {"RoutingStrategy.RAND", "RoutingStrategy.RROBIN"};
    for (std::size_t p = 0; p < 2; ++p) {
        std::printf("\n=== Dispatching: %s ===\n", policy_name[p]);
        Net model = cluster_open(lambda, d, same_sched(4, SchedStrategy::PS), single_servers(4),
                                 policies[p]);
        print_avg(solve_avg("SSA", model, ssa));
    }
}

/** The arrival-rate sweep of a two-server PS farm, up to 75% utilization. */
void cl_sweep() {
    const Demands d = uniform_demands(2, 1.0);
    const double rates[4] = {0.2, 0.5, 0.9, 1.5};
    for (std::size_t k = 0; k < 4; ++k) {
        std::printf("\n=== lambda = %g ===\n", rates[k]);
        const std::vector<double> lambda(1, rates[k]);
        Net model = cluster_open(lambda, d, same_sched(2, SchedStrategy::PS), single_servers(2),
                                 RoutingStrategy::RAND);
        print_avg(solve_avg("MVA", model));
    }
}

/** Three PS servers, ten closed jobs: MVA against simulation. */
void dispatch_closed() {
    const std::vector<double> N(1, 10.0), Z(1, 1.0);
    const Demands d = uniform_demands(3, 1.0);
    Net model = cluster_closed(N, Z, d, same_sched(3, SchedStrategy::PS), single_servers(3),
                               RoutingStrategy::RAND);
    print_avg(solve_avg("MVA", model));

    SolverOpts ssa;
    ssa.seed = 23000;
    ssa.samples = 20000;
    print_avg(solve_avg("SSA", model, ssa));
}

/** Two PS servers shared by one open and one closed class. */
void dispatch_mixed() {
    const std::vector<double> lambda(1, 0.5), N(1, 3.0), Z(1, 1.0);
    Demands d(2, std::vector<double>());
    for (std::size_t i = 0; i < 2; ++i) {
        d[i].push_back(1.0 / 2.0);   // open class, mu = 2.0
        d[i].push_back(1.0 / 1.5);   // closed class, mu = 1.5
    }
    Net model = cluster_mixed(lambda, N, Z, d, same_sched(2, SchedStrategy::PS), single_servers(2),
                              RoutingStrategy::RAND);
    print_avg(solve_avg("MVA", model));

    // TODO(cpp): LDES(cluster.build(), seed=23000).get_avg_table()
    na("LDES", "SolverLDES is the SSJ simulation engine of the JAR and has no C++ port");
}

/**
 * `Cluster.compareScheduling`: the same farm under FCFS and PS.
 *
 * With exponential service the two disciplines agree in the mean, so the
 * reference makes the service SCV 4.0 to separate them: `Cluster` then fits
 * `APH.fitMeanAndSCV(1/mu, scv)` at every server instead of an exponential,
 * which is what `aph_fit` builds here.
 */
void cl_scheduling() {
    const std::size_t M = 3;
    const double lambda = 0.9, mu = 0.5, scv = 4.0;

    SolverOpts ssa;
    ssa.seed = 23000;
    ssa.samples = 20000;

    const SchedStrategy disciplines[2] = {SchedStrategy::FCFS, SchedStrategy::PS};
    const char* discipline_name[2] = {"FCFS", "PS"};
    for (std::size_t k = 0; k < 2; ++k) {
        std::printf("\n=== Scheduling: %s ===\n", discipline_name[k]);
        Net model = cluster_open(std::vector<double>(1, lambda), uniform_demands(M, 1.0 / mu),
                                 same_sched(M, disciplines[k]), single_servers(M),
                                 RoutingStrategy::RAND, scv);
        print_avg(solve_avg("SSA", model, ssa));
    }
}

/**
 * `Cluster.sweepNumStations`: response time against the number of servers at a
 * fixed total arrival rate.
 *
 * The single-station service rate is replicated across M servers, so adding a
 * server splits the same arrival stream: per-station utilization falls and
 * response time drops towards the bare service time.
 */
void cl_stations() {
    const double lambda = 1.6, mu = 1.0;
    const std::size_t counts[4] = {2, 3, 4, 6};
    for (std::size_t k = 0; k < 4; ++k) {
        const std::size_t M = counts[k];
        Net model = cluster_open(std::vector<double>(1, lambda), uniform_demands(M, 1.0 / mu),
                                 same_sched(M, SchedStrategy::PS), single_servers(M),
                                 RoutingStrategy::RAND);
        const AvgTable t = solve_avg("MVA", model);
        double sum_respt = 0.0, max_util = 0.0;
        std::size_t n = 0;
        for (std::size_t i = 0; i < t.Station.size(); ++i) {
            if (t.Station[i] == "Source") continue;
            sum_respt += t.RespT[i];
            max_util = std::max(max_util, t.Util[i]);
            ++n;
        }
        std::printf("\n=== M = %zu  (mean RespT %.4f, max Util %.4f) ===\n", M,
                    n ? sum_respt / static_cast<double>(n) : 0.0, max_util);
        print_avg(t);
    }
}

/** Three PS servers, one open class at rate 0.4: MVA against simulation. */
void dispatch_open() {
    const std::vector<double> lambda(1, 0.4);
    const Demands d = uniform_demands(3, 1.0);
    Net model = cluster_open(lambda, d, same_sched(3, SchedStrategy::PS), single_servers(3),
                             RoutingStrategy::RAND);
    print_avg(solve_avg("MVA", model));

    SolverOpts ssa;
    ssa.seed = 23000;
    ssa.samples = 20000;
    print_avg(solve_avg("SSA", model, ssa));
}

LINE_EXAMPLE("basic/cluster", cl_basic);
LINE_EXAMPLE("basic/cluster", cl_closed);
LINE_EXAMPLE("basic/cluster", cl_compare);
LINE_EXAMPLE("basic/cluster", cl_multiclass);
LINE_EXAMPLE("basic/cluster", cl_scheduling);
LINE_EXAMPLE("basic/cluster", cl_stations);
LINE_EXAMPLE("basic/cluster", cl_sweep);
LINE_EXAMPLE("basic/cluster", dispatch_closed);
LINE_EXAMPLE("basic/cluster", dispatch_mixed);
LINE_EXAMPLE("basic/cluster", dispatch_open);

}  // namespace examples
}  // namespace line
