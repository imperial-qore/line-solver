/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/advanced/stateProbabilities/`: the four state-probability
 * queries, `getProbAggr`, `getProb`, `getProbSysAggr` and `getProbSys`.
 *
 * WHAT A `station.setState(n)` MEANS HERE. The reference scripts set a state on
 * the model and then ask a solver for its probability; this port has no state
 * on the Network object, so the queried state is built where it is asked for,
 * through `State.fromMarginal` -- the same encoding `setState` would have
 * produced. A station the reference leaves at -1 is IGNORED by the query it
 * feeds (`getProbAggr` reads one station's row, `getProbSysAggr` reads them
 * all), and is carried here at the model's default marginal so the object is a
 * well-formed encoding rather than a hole.
 *
 * THE DETAILED QUERIES TAKE THE FIRST REALIZATION. `fromMarginal` returns every
 * local row realizing a marginal -- phases, and the buffer order of an FCFS
 * station -- and `getProb`/`getProbSys` are per-ROW answers, so one row has to
 * be named. Row zero is taken, which is what `Network.initDefault` and this
 * port's `default_init_state` both take.
 *
 * JMT is refused by name: this port does not carry it.
 */

#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "examples_common.h"
#include "line/lang/qn/state.h"
#include "line/api/pfqn/pfqn_comb_common.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_prob.h"
#include "line/solvers/nc/solver_nc_prob.h"
#include "line/solvers/ssa/solver_ssa_getters.h"
#include "line/solvers/ssa/solver_ssa_serial.h"

namespace line {
namespace examples {

namespace {

/** The per-station, per-class counts a `setState` sequence names. */
using Counts = std::vector<std::vector<std::size_t>>;

/**
 * "Compute the normalizing constant here", the `nargin == 2` branch of
 * `solver_nc_margaggr.m`.
 *
 * A LOAD-DEPENDENT MODEL MAKES THIS LOAD-BEARING. `@@SolverNC/getProbAggr.m`
 * hands the analyzer's constant on only when a previous probability query
 * cached one; on a first call it passes none, and the constant is then the
 * `pfqn_ncld` one, which carries the multiserver rate vector. Feeding the
 * analyzer's `cub` constant instead reports 0.3254 where the reference and the
 * exact chain both report 0.34 on `statepr_aggr`.
 */
const double kNoLg = std::numeric_limits<double>::quiet_NaN();

/** The model's own default marginal: every closed class at its reference station. */
Counts default_counts(const Sn& sn) {
    Counts n(sn.nstations, std::vector<std::size_t>(sn.nclasses, 0));
    for (std::size_t r = 0; r < sn.nclasses; ++r) {
        const double pop = sn.classes[r].population;
        const std::size_t rs = sn.classes[r].refstat;
        if (pop > 0 && pop < 1e18 && rs >= 1 && rs <= sn.nstations)
            n[rs - 1][r] = static_cast<std::size_t>(pop);
    }
    return n;
}

/** `State.fromMarginal` over every stateful node: the queried network state. */
qn::NetState<double> state_of(const Sn& sn, const Counts& n) {
    qn::NetState<double> st;
    st.local.assign(sn.stateful_nodes.size(), std::vector<double>());
    for (std::size_t f = 0; f < sn.stateful_nodes.size(); ++f) {
        const std::size_t ind = sn.stateful_nodes[f];
        const std::size_t ist = sn.nodes[ind - 1].station;
        std::vector<std::size_t> marg(sn.nclasses, 0), ph(sn.nclasses, 1);
        if (ist != 0) {
            for (std::size_t r = 0; r < sn.nclasses; ++r) ph[r] = sn.phases_of(ist, r + 1);
            marg = n[ist - 1];
        }
        const std::vector<std::vector<double>> rows = qn::from_marginal_node(sn, ind, marg, ph);
        if (rows.empty())
            throw InputError("stateProbabilities: node '" + sn.nodes[ind - 1].name +
                             "' admits no state with the requested per-class counts");
        st.local[f] = rows[0];
    }
    return st;
}

/** The same state at the encoding width a sample path lives at, left-padded. */
qn::NetState<double> widen(const qn::NetState<double>& st, const qn::NetState<double>& ref) {
    qn::NetState<double> out = st;
    for (std::size_t f = 0; f < out.local.size() && f < ref.local.size(); ++f) {
        const std::size_t w = ref.local[f].size();
        if (out.local[f].size() < w)
            out.local[f].insert(out.local[f].begin(), w - out.local[f].size(), 0.0);
    }
    return out;
}

/** The counts as SolverNC's marginal state, a station per row. */
nc::MarginalState nc_marginal(const Sn& sn, const Counts& n) {
    nc::MarginalState nir(sn.nstations, std::vector<int>(sn.nclasses, 0));
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t r = 0; r < sn.nclasses; ++r) nir[i][r] = static_cast<int>(n[i][r]);
    return nir;
}

/** `[a, b, ...]`, the shape the reference prints a queried state in. */
std::string counts_text(const std::vector<std::size_t>& n) {
    std::string s = "[";
    for (std::size_t r = 0; r < n.size(); ++r) {
        if (r) s += ", ";
        s += std::to_string(n[r]);
    }
    return s + "]";
}

/** Delay + Queue1 + Queue2, the three-station shape five of the six share. */
Net three_station_model(SchedStrategy sched2, std::size_t& delay, std::size_t& q1,
                        std::size_t& q2) {
    Net m("model");
    delay = m.add_delay("Delay");
    q1 = m.add_queue("Queue1", SchedStrategy::PS);
    q2 = m.add_queue("Queue2", sched2);
    m.set_number_of_servers(q2, 2.0);
    return m;
}

}  // namespace

/**
 * `statepr_aggr.py`: the marginal aggregate probability of one station, by the
 * exact chain and by the normalizing constants.
 */
void statepr_aggr() {
    std::size_t delay = 0, q1 = 0, q2 = 0;
    Net m = three_station_model(SchedStrategy::PS, delay, q1, q2);
    ClosedClass c1(m, "Class1", 2, delay, 0);
    ClosedClass c2(m, "Class2", 0, delay, 0);
    m.set_service(delay, c1, Exp(1.0));
    m.set_service(delay, c2, Exp(1.0));
    m.set_service(q1, c1, Exp(3.0));
    m.set_service(q1, c2, Exp(4.0));
    m.set_service(q2, c1, Exp(1.0));
    m.set_service(q2, c2, Exp(3.0));
    Routing P;
    cyclic(P, c1, {delay, q1, q2});
    cyclic(P, c2, {delay, q1, q2});
    m.link(P);

    const Sn& sn = m.get_struct();
    const std::size_t ist = m.station_index(q2);
    Counts n = default_counts(sn);
    n[ist - 1].assign(sn.nclasses, 0);
    kv("Queried station", sn.stations[ist - 1].name);
    kv("Queried state", counts_text(n[ist - 1]));

    section("CTMC");
    const ctmc::CtmcOptions copt;
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, copt);
    const double pr_ctmc = ctmc::solver_ctmc_margaggr(sn, d, state_of(sn, n))[ist - 1];
    kv("getProbAggr", pr_ctmc);
    bare(pr_ctmc);  // the line the reference prints bare; see examples_common.h

    section("NC");
    const nc::NcSolverOptions nopt;
    const std::vector<int> row(n[ist - 1].begin(), n[ist - 1].end());
    kv("getProbAggr", nc::solver_nc_getprob_aggr(sn, nopt, ist, row, kNoLg));
}

/**
 * `statepr_aggr_large.py`: the same query on a four-class model whose class
 * switches all happen at Queue2.
 */
void statepr_aggr_large() {
    std::size_t delay = 0, q1 = 0, q2 = 0;
    Net m = three_station_model(SchedStrategy::PS, delay, q1, q2);
    ClosedClass c1(m, "Class1", 1, delay, 0);
    ClosedClass c2(m, "Class2", 0, delay, 0);
    ClosedClass c3(m, "Class3", 4, delay, 0);
    ClosedClass c4(m, "Class4", 0, delay, 0);
    m.set_service(delay, c1, D::exp_mean(1.0));
    m.set_service(delay, c2, D::exp_mean(0.5));
    m.set_service(delay, c3, D::exp_mean(1.0));
    m.set_service(delay, c4, D::exp_mean(1.0));
    m.set_service(q1, c1, D::exp_mean(1.0 / 3.0));
    m.set_service(q1, c2, D::exp_mean(0.25));
    m.set_service(q1, c3, D::exp_mean(0.2));
    m.set_service(q1, c4, D::exp_mean(1.0));
    m.set_service(q2, c1, D::exp_mean(1.0));
    m.set_service(q2, c2, D::exp_mean(1.0 / 3.0));
    m.set_service(q2, c3, D::exp_mean(0.2));
    m.set_service(q2, c4, D::exp_mean(0.5));

    Routing P;
    serial(P, c1, {delay, q1, q2});
    P.set(c1, c2, q2, delay, 1.0);
    serial(P, c2, {delay, q1, q2});
    P.set(c2, c1, q2, delay, 1.0);
    serial(P, c3, {delay, q1, q2});
    P.set(c3, c4, q2, delay, 1.0);
    P.set(c4, c3, q2, delay, 1.0);
    P.set(c4, c4, delay, q2, 1.0);
    m.link(P);

    const Sn& sn = m.get_struct();
    const std::size_t ist = m.station_index(q2);
    Counts n = default_counts(sn);
    n[ist - 1] = std::vector<std::size_t>{1, 0, 2, 1};
    kv("Queried station", sn.stations[ist - 1].name);
    kv("Queried state", counts_text(n[ist - 1]));

    section("CTMC");
    const ctmc::CtmcOptions copt;
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, copt);
    const double pr_ctmc = ctmc::solver_ctmc_margaggr(sn, d, state_of(sn, n))[ist - 1];
    kv("getProbAggr", pr_ctmc);
    bare(pr_ctmc);  // the line the reference prints bare; see examples_common.h

    section("NC");
    const nc::NcSolverOptions nopt;
    const std::vector<int> row(n[ist - 1].begin(), n[ist - 1].end());
    kv("getProbAggr", nc::solver_nc_getprob_aggr(sn, nopt, ist, row, kNoLg));
}

namespace {

/**
 * The two-class class-switching model `statepr_allprobs_*` share: Class1 is
 * served at Queue1 and becomes Class2 on its way to Queue2, which returns it as
 * Class1.
 */
Net allprobs_model(SchedStrategy sched2, double rate21, double rate22) {
    std::size_t delay = 0, q1 = 0, q2 = 0;
    Net m = three_station_model(sched2, delay, q1, q2);
    ClosedClass c1(m, "Class1", 2, delay, 0);
    ClosedClass c2(m, "Class2", 0, delay, 0);
    m.set_service(delay, c1, Exp(1.0));
    m.set_service(delay, c2, Exp(1.0));
    m.set_service(q1, c1, Exp(3.0));
    m.set_service(q1, c2, Exp(4.0));
    m.set_service(q2, c1, Exp(rate21));
    m.set_service(q2, c2, Exp(rate22));
    Routing P;
    P.set(c1, c1, delay, q1, 1.0);
    P.set(c1, c1, q2, delay, 1.0);
    P.set(c1, c2, q1, q2, 1.0);
    P.set(c2, c1, delay, q1, 1.0);
    P.set(c2, c1, q2, delay, 1.0);
    P.set(c2, c2, q1, q2, 1.0);
    m.link(P);
    return m;
}

/** The four queries of `statepr_allprobs_*`, under CTMC, NC and SSA. */
void allprobs_report(Net& m, std::size_t samples) {
    const Sn& sn = m.get_struct();
    const std::size_t ist = sn.nstations;
    const std::size_t ind = sn.node_of_station(ist);
    Counts n(sn.nstations, std::vector<std::size_t>(sn.nclasses, 0));
    n[1] = std::vector<std::size_t>{1, 0};
    n[2] = std::vector<std::size_t>{0, 1};
    const qn::NetState<double> st = state_of(sn, n);
    kv("Queried station", sn.stations[ist - 1].name);
    kv("Queried state", counts_text(n[ist - 1]));

    const ctmc::CtmcOptions copt;
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, copt);
    const nc::NcSolverOptions nopt;
    ssa::SsaSerialOptions sopt;
    sopt.method = "serial";
    sopt.samples = samples;
    sopt.seed = 23000;
    const ssa::SsaSerialSolution<double> sim = ssa::solver_ssa_serial_analyzer(sn, sopt);
    const qn::NetState<double> wide = sim.run.space.empty() ? st : widen(st, sim.run.space[0]);
    const std::vector<double> counts =
        ssa::getters_detail::node_counts(sn, ind, wide.local[sn.stateful_of_station(ist) - 1]);

    section("CTMC");
    const double pr_ctmc = ctmc::solver_ctmc_margaggr(sn, d, st)[ist - 1];
    kv("getProbAggr", pr_ctmc);
    bare(pr_ctmc);  // the line the reference prints bare; see examples_common.h
    kv("getProb", ctmc::solver_ctmc_marg(sn, d, st)[ist - 1]);
    kv("getProbSysAggr", ctmc::solver_ctmc_jointaggr(sn, d, st));
    kv("getProbSys", ctmc::solver_ctmc_joint(sn, d, st));

    section("NC");
    const std::vector<int> row(n[ist - 1].begin(), n[ist - 1].end());
    kv("getProbAggr", nc::solver_nc_getprob_aggr(sn, nopt, ist, row, kNoLg));
    kv("getProbSysAggr", nc::solver_nc_getprob_sys_aggr(sn, nopt, nc_marginal(sn, n)));

    section("SSA");
    kv("getProbAggr", ssa::ssa_prob_aggr(sn, sim.run, ind, counts).prob);
    kv("getProb", ssa::ssa_prob(sn, sim.run, ind, wide.local[sn.stateful_of_station(ist) - 1]).prob);
    kv("getProbSysAggr", ssa::ssa_prob_sys_aggr(sn, sim.run, wide).prob);
    kv("getProbSys", ssa::ssa_prob_sys(sn, sim.run, wide).prob);

    // TODO(cpp): solver_jmt = JMT(model, options); pr = solver_jmt.getProbAggr(target_station)
    // TODO(cpp): pr = solver_jmt.getProbSysAggr()
    // The JMT engine runs here (`jmt_avg`); what it has no port of is the
    // PROBABILITY getter. `getProbAggr` weighs the trajectory by the time spent
    // in the model's CURRENT state, and this port's NetworkStruct carries an
    // initial marking, a state prior and a state space -- no current state row
    // to weigh against. `line-cli -s jmt -a prob` refuses by the same name.
    na("JMT", "SolverJMT getProbAggr/getProbSysAggr weigh the trajectory against the model's "
              "CURRENT state, which this port's NetworkStruct does not carry; the mean tables of "
              "the same engine DO run here");
}

}  // namespace

/** `statepr_allprobs_fcfs.py`: all four queries, Queue2 FCFS with two servers. */
void statepr_allprobs_fcfs() {
    Net m = allprobs_model(SchedStrategy::FCFS, 3.0, 3.0);
    allprobs_report(m, 100000);
}

/** `statepr_allprobs_ps.py`: the same four, Queue2 processor sharing. */
void statepr_allprobs_ps() {
    Net m = allprobs_model(SchedStrategy::PS, 1.0, 3.0);
    allprobs_report(m, 20000);
}

/**
 * `statepr_sys_aggr.py`: the joint aggregate probability that every job sits at
 * Queue2, on a four-class model with class switching there.
 */
void statepr_sys_aggr() {
    std::size_t delay = 0, q1 = 0, q2 = 0;
    Net m = three_station_model(SchedStrategy::PS, delay, q1, q2);
    ClosedClass c1(m, "Class1", 1, delay, 0);
    ClosedClass c2(m, "Class2", 0, delay, 0);
    ClosedClass c3(m, "Class3", 3, delay, 0);
    ClosedClass c4(m, "Class4", 0, delay, 0);
    m.set_service(delay, c1, Exp(1.0));
    m.set_service(delay, c2, Exp(2.0));
    m.set_service(delay, c3, Exp(1.0));
    m.set_service(delay, c4, Exp(1.0));
    m.set_service(q1, c1, Exp(3.0));
    m.set_service(q1, c2, Exp(4.0));
    m.set_service(q1, c3, Exp(5.0));
    m.set_service(q1, c4, Exp(1.0));
    m.set_service(q2, c1, Exp(1.0));
    m.set_service(q2, c2, Exp(3.0));
    m.set_service(q2, c3, Exp(5.0));
    m.set_service(q2, c4, Exp(2.0));

    Routing P;
    serial(P, c1, {delay, q1, q2});
    P.set(c1, c2, q2, delay, 1.0);
    serial(P, c2, {delay, q1, q2});
    P.set(c2, c1, q2, delay, 1.0);
    serial(P, c3, {delay, q1, q2});
    P.set(c3, c4, q2, delay, 1.0);
    P.set(c4, c3, q2, delay, 1.0);
    P.set(c4, c4, delay, q2, 1.0);
    m.link(P);

    const Sn& sn = m.get_struct();
    Counts n(sn.nstations, std::vector<std::size_t>(sn.nclasses, 0));
    n[sn.nstations - 1] = std::vector<std::size_t>{1, 0, 3, 0};
    kv("Queried state", counts_text(n[sn.nstations - 1]) + " at " +
                            sn.stations[sn.nstations - 1].name);

    section("CTMC");
    const ctmc::CtmcOptions copt;
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, copt);
    const double pr_ctmc = ctmc::solver_ctmc_jointaggr(sn, d, state_of(sn, n));
    kv("getProbSysAggr", pr_ctmc);
    bare(pr_ctmc);  // the line the reference prints bare; see examples_common.h

    section("NC (exact)");
    nc::NcSolverOptions nopt;
    nopt.method = "exact";
    kv("getProbSysAggr", nc::solver_nc_getprob_sys_aggr(sn, nopt, nc_marginal(sn, n)));

    // TODO(cpp): solver_jmt = JMT(model, samples=int(1e4), seed=532733); pr_jmt = solver_jmt.getProbSysAggr()
    na("JMT", "SolverJMT getProbSysAggr weighs the trajectory against the model's CURRENT "
              "state, which this port's NetworkStruct does not carry; the mean tables of the "
              "same engine DO run here");
}

/**
 * `statepr_sys_aggr_large.py`: the same joint query on three queues and four
 * unit-population classes. The MATLAB reference reports about 0.000348.
 */
void statepr_sys_aggr_large() {
    Net m("model");
    Queue q1(m, "Queue1", SchedStrategy::PS);
    Queue q2(m, "Queue2", SchedStrategy::PS);
    Queue q3(m, "Queue3", SchedStrategy::PS);
    q3.set_number_of_servers(3.0);
    ClosedClass c1(m, "Class1", 1, q1, 0);
    ClosedClass c2(m, "Class2", 1, q1, 0);
    ClosedClass c3(m, "Class3", 1, q1, 0);
    ClosedClass c4(m, "Class4", 1, q1, 0);
    m.set_service(q1, c1, Exp(1.0));
    m.set_service(q1, c2, Exp(2.0));
    m.set_service(q1, c3, Exp(1.0));
    m.set_service(q1, c4, Exp(1.0));
    m.set_service(q2, c1, Exp(3.0));
    m.set_service(q2, c2, Exp(4.0));
    m.set_service(q2, c3, Exp(5.0));
    m.set_service(q2, c4, Exp(1.0));
    q3.set_service(c1, Exp(1.0));
    q3.set_service(c2, Exp(3.0));
    q3.set_service(c3, Exp(5.0));
    q3.set_service(c4, Exp(2.0));

    Routing P;
    serial(P, c1, {q1, q2, q3});
    P.set(c1, c2, q3, q1, 1.0);
    serial(P, c2, {q1, q2, q3});
    P.set(c2, c1, q3, q1, 1.0);
    serial(P, c3, {q1, q2, q3});
    P.set(c3, c4, q3, q1, 1.0);
    P.set(c4, c3, q3, q1, 1.0);
    P.set(c4, c4, q1, q3, 1.0);
    m.link(P);

    const Sn& sn = m.get_struct();
    Counts n(sn.nstations, std::vector<std::size_t>(sn.nclasses, 0));
    n[sn.nstations - 1] = std::vector<std::size_t>{1, 1, 1, 1};
    kv("Queried state", counts_text(n[sn.nstations - 1]) + " at " +
                            sn.stations[sn.nstations - 1].name);

    section("CTMC");
    const ctmc::CtmcOptions copt;
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, copt);
    const double pr_ctmc = ctmc::solver_ctmc_jointaggr(sn, d, state_of(sn, n));
    kv("getProbSysAggr", pr_ctmc);
    bare(pr_ctmc);  // the line the reference prints bare; see examples_common.h
    kv("MATLAB reference", 0.000348);

    section("NC (exact)");
    nc::NcSolverOptions nopt;
    nopt.method = "exact";
    kv("getProbSysAggr", nc::solver_nc_getprob_sys_aggr(sn, nopt, nc_marginal(sn, n)));

    // TODO(cpp): solver_jmt = JMT(model, samples=int(1e4), seed=532733); pr_jmt = solver_jmt.getProbSysAggr()
    na("JMT", "SolverJMT getProbSysAggr weighs the trajectory against the model's CURRENT "
              "state, which this port's NetworkStruct does not carry; the mean tables of the "
              "same engine DO run here");
}

/**
 * `statepr_sys_marg.m` / `.py`: the JOINT law of the per-station TOTAL queue
 * lengths, all classes summed out, from `getProbSysMarg`.
 *
 * NEITHER of the two above. `statepr_sys_aggr` fixes the PER-CLASS population of
 * every station and is a product form; each probability here is the SUM of that
 * one over every per-class table with these row sums. The fibre grows
 * combinatorially, so the quantity is evaluated as a permanent of the demand
 * matrix replicated once per job (Ryser 1963) rather than by enumerating it.
 *
 * The whole lattice of total states is swept, so the printed column sums to one
 * and the first moments of the law are the mean queue lengths -- which is what
 * makes the CTMC comparison below a check on every state and not only on the
 * normalization.
 */
void statepr_sys_marg() {
    Net m("model");
    Delay delay(m, "Delay");
    Queue q1(m, "Queue1", SchedStrategy::PS);
    Queue q2(m, "Queue2", SchedStrategy::PS);
    ClosedClass c1(m, "Class1", 2, delay, 0);
    ClosedClass c2(m, "Class2", 1, delay, 0);
    m.set_service(delay, c1, Exp(1.0 / 1.5));
    m.set_service(delay, c2, Exp(1.0 / 2.0));
    m.set_service(q1, c1, Exp(1.0 / 0.7));
    m.set_service(q1, c2, Exp(1.0 / 0.4));
    m.set_service(q2, c1, Exp(1.0 / 0.3));
    m.set_service(q2, c2, Exp(1.0 / 0.9));

    // CYCLIC, not serial: `Network.serialRouting` closes the cycle for a closed
    // model, and an unclosed chain leaves the normalizing constant zero.
    Routing P;
    cyclic(P, c1, {delay, q1, q2});
    cyclic(P, c2, {delay, q1, q2});
    m.link(P);

    const Sn& sn = m.get_struct();
    const int total = 3;
    const std::vector<std::vector<int> > states =
        pfqn::multichoose_rows(static_cast<int>(sn.nstations), total);

    section("NC (exact permanent)");
    const nc::NcSolverOptions nopt;
    std::vector<double> p(states.size(), 0.0), mean(sn.nstations, 0.0);
    double sum = 0.0;
    std::printf("  n(Delay) n(Queue1) n(Queue2)        P(n)\n");
    for (std::size_t j = 0; j < states.size(); ++j) {
        p[j] = nc::solver_nc_getprob_sys_marg(sn, nopt, states[j]);
        sum += p[j];
        for (std::size_t i = 0; i < sn.nstations; ++i) mean[i] += p[j] * states[j][i];
        std::printf("  %8d %9d %9d  %10.6f\n", states[j][0], states[j][1], states[j][2], p[j]);
    }
    std::printf("  ------------------------------------------\n");
    kv("sum", sum);

    // The law is exact, so its first moments are the queue lengths.
    const mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer(sn, ctmc::CtmcOptions());
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        double q = 0.0;
        for (std::size_t c = 0; c < r.QN.cols(); ++c) q += r.QN(i, c);
        kv("E[n] " + sn.stations[i].name + " (joint law)", mean[i]);
        kv("QLen " + sn.stations[i].name + " (CTMC)", q);
    }

    // The approximate engines trade accuracy for cost on models whose class
    // count makes the exact expansion dear. They need FULL SUPPORT and refuse a
    // structural zero rather than flooring it.
    section("NC (Bethe permanent)");
    double sb = 0.0, err = 0.0;
    std::vector<double> pb(states.size(), 0.0);
    for (std::size_t j = 0; j < states.size(); ++j) {
        pb[j] = nc::solver_nc_getprob_sys_marg(sn, nopt, states[j], "bethe");
        sb += pb[j];
    }
    for (std::size_t j = 0; j < states.size(); ++j) err += std::fabs(pb[j] / sb - p[j]) / p[j];
    kv("mean relative error (%)", 100.0 * err / static_cast<double>(states.size()));
}

LINE_EXAMPLE("advanced/stateProbabilities", statepr_aggr);
LINE_EXAMPLE("advanced/stateProbabilities", statepr_aggr_large);
LINE_EXAMPLE("advanced/stateProbabilities", statepr_allprobs_fcfs);
LINE_EXAMPLE("advanced/stateProbabilities", statepr_allprobs_ps);
LINE_EXAMPLE("advanced/stateProbabilities", statepr_sys_aggr);
LINE_EXAMPLE("advanced/stateProbabilities", statepr_sys_aggr_large);
LINE_EXAMPLE("advanced/stateProbabilities", statepr_sys_marg);

}  // namespace examples
}  // namespace line
