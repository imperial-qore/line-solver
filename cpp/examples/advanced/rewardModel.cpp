/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/advanced/rewardModel/`: `setReward` and the SolverCTMC
 * reward surface, `getAvgReward` and `getTranReward`.
 *
 * WHAT A REWARD SEES. The reference hands its callback a `RewardState` wrapping
 * the AGGREGATE state row and offers `state.at(node)` / `state.at(node, class)`
 * on it; this port hands the same row as a plain vector in `(ist-1)*K + r`
 * order, so the accessors are the two helpers below and nothing else changes.
 * The station index is the STATION one, not the node one -- a Source occupies a
 * block of its own -- which is why each example resolves it with
 * `station_index` before declaring a reward.
 *
 * THE `Reward.*` TEMPLATES ARE FUNCTIONS HERE. `Reward.queue_length`,
 * `Reward.utilization` and `Reward.blocking` are Python factory methods with no
 * C++ counterpart, so they are reproduced as the file-local factories below, at
 * their reference semantics: utilization is min(jobs, nservers) and blocking is
 * the indicator of the station sitting at its capacity.
 *
 * The example is spelled two ways and both are registered here. MATLAB, which
 * names the corpus, calls it `rewardModel_*`; the Python twin's file stem is
 * `reward_model_*`. There is NO symlink reconciling them -- an earlier version
 * of this comment said there was. The mapping is data, in
 * `goldens/pythonScriptAliases` (`corpus.json`), read by the parity suite and
 * by `doc/latex/exdoc.py`.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "examples_common.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_reward.h"

namespace line {
namespace examples {

namespace {

using RewardFn = std::function<double(const std::vector<double>&)>;

/** `state.at(station, class)`: one entry of the aggregate row. */
double at(const std::vector<double>& s, std::size_t ist, std::size_t K, std::size_t r) {
    const std::size_t j = (ist - 1) * K + (r - 1);
    return j < s.size() ? s[j] : 0.0;
}

/** `state.at(station).total()`: the station's whole block. */
double total_at(const std::vector<double>& s, std::size_t ist, std::size_t K) {
    double acc = 0.0;
    for (std::size_t r = 1; r <= K; ++r) acc += at(s, ist, K, r);
    return acc;
}

/** `Reward.queue_length(node)` and `Reward.queue_length(node, class)`. */
RewardFn reward_qlen(std::size_t ist, std::size_t K) {
    return [ist, K](const std::vector<double>& s) { return total_at(s, ist, K); };
}

RewardFn reward_qlen(std::size_t ist, std::size_t K, std::size_t r) {
    return [ist, K, r](const std::vector<double>& s) { return at(s, ist, K, r); };
}

/** `Reward.utilization(node[, class])`: min(jobs, nservers). */
RewardFn reward_util(std::size_t ist, std::size_t K, double nservers) {
    return [ist, K, nservers](const std::vector<double>& s) {
        return std::min(total_at(s, ist, K), nservers);
    };
}

RewardFn reward_util(std::size_t ist, std::size_t K, std::size_t r, double nservers) {
    return [ist, K, r, nservers](const std::vector<double>& s) {
        return std::min(at(s, ist, K, r), nservers);
    };
}

/** `Reward.blocking(node)`: the indicator of the station sitting at capacity. */
RewardFn reward_blocking(std::size_t ist, std::size_t K, double cap) {
    return [ist, K, cap](const std::vector<double>& s) {
        return total_at(s, ist, K) >= cap ? 1.0 : 0.0;
    };
}

/** The steady-state expectations, printed one per line as the reference does.
 *
 * The banner is part of the shape, not decoration: the shared parity parser
 * recognises a reward block by `=== Steady-State Expected Rewards ===` and
 * reads the `name: value` lines that follow it. Without it the values are
 * printed, parsed by nothing, and the row compares no cell. */
void print_rewards(const std::vector<std::string>& names, const std::vector<double>& R) {
    std::printf("=== Steady-State Expected Rewards ===\n");
    // `%.6f`, the reference's own width, and NOT `kv`'s `%.6g`: the golden holds
    // seven significant digits (2.123077) and the parity tolerance on an exact
    // CTMC reward is 1e-6, which a six-significant-digit print cannot meet.
    // Recorded under the solver `section()` declared, keyed as the golden keys
    // it: the NAMES are the example's labels for elements of one reward vector,
    // so the vector alone cannot supply them and nothing downstream can recover
    // them from the printed block without parsing it back.
    for (std::size_t i = 0; i < R.size(); ++i) {
        std::printf("%15s: %.6f\n", names[i].c_str(), R[i]);
        derived_here(names[i], "Reward", R[i]);
    }
}

/** One `LINE = ..., Analytical = ..., Error = ...` comparison row. */
void compare(const std::string& name, double line_value, double analytical) {
    std::printf("%-24s LINE = %12.6f  Analytical = %12.6f  Error = %.2e\n", name.c_str(),
                line_value, analytical, std::fabs(line_value - analytical));
}

/** The M/M/1/K stationary law pi(n), n = 0..K, at the reference's parameters. */
std::vector<double> mm1k_pi(double rho, std::size_t K) {
    std::vector<double> pi(K + 1, 0.0);
    if (rho == 1.0) {
        for (std::size_t n = 0; n <= K; ++n) pi[n] = 1.0 / static_cast<double>(K + 1);
        return pi;
    }
    for (std::size_t n = 0; n <= K; ++n)
        pi[n] = (1.0 - rho) * std::pow(rho, static_cast<double>(n)) /
                (1.0 - std::pow(rho, static_cast<double>(K + 1)));
    return pi;
}

/** Source -> Queue -> Sink, the shape every reward example is built on. */
Net reward_model(const std::string& nm, SchedStrategy sched, double nservers, double cap,
                 std::size_t& source, std::size_t& queue, std::size_t& sink) {
    Net m(nm);
    source = m.add_source("Source");
    queue = m.add_queue("Queue", sched);
    sink = m.add_sink("Sink");
    m.set_number_of_servers(queue, nservers);
    m.set_capacity(queue, cap);
    return m;
}

}  // namespace

/**
 * `reward_model_mm1k.py`: the four rewards of an M/M/1/K, in steady state and
 * over a finite horizon, against the analytical law.
 */
void reward_model_mm1k() {
    std::size_t source = 0, queue = 0, sink = 0;
    Net m = reward_model("RewardExample", SchedStrategy::FCFS, 1.0, 3.0, source, queue, sink);
    OpenClass oclass(m, "Class1");
    m.set_arrival(source, oclass, Exp(2.0));
    m.set_service(queue, oclass, Exp(3.0));
    Routing P;
    serial(P, oclass, {source, queue, sink});
    m.link(P);

    const std::size_t ist = m.station_index(queue);
    const std::size_t K = 1;
    m.set_reward("QueueLength", reward_qlen(ist, K, 1));
    m.set_reward("Utilization", reward_util(ist, K, 1, 1.0));
    m.set_reward("BlockingProb", reward_blocking(ist, K, 3.0));
    m.set_reward("QueueCost", [ist, K](const std::vector<double>& s) {
        const double n = at(s, ist, K, 1);
        return n * n;
    });

    const Sn& sn = m.get_struct();
    const ctmc::CtmcOptions opt;
    std::vector<std::string> names;
    const std::vector<double> R = ctmc::solver_ctmc_avg_reward(sn, opt, &names);

    section("CTMC");
    print_rewards(names, R);

    section("CTMC transient, timespan [0, 5]");
    std::vector<double> t;
    std::vector<std::string> tnames;
    const std::vector<std::vector<double>> Rt =
        ctmc::solver_ctmc_tran_reward(sn, opt, 0.0, 5.0, &t, &tnames);
    kv("Time points", static_cast<double>(t.size()));
    for (std::size_t i = 0; i < Rt.size(); ++i)
        std::printf("%-24s E[r(X(0))] = %12.6f  E[r(X(T))] = %12.6f\n", tnames[i].c_str(),
                    Rt[i].front(), Rt[i].back());

    // The analytical comparison of the reference, at ITS K: the script fixes
    // K = 10 here while the model's capacity is 3, and that is reproduced.
    section("M/M/1/K analytical");
    const double rho = 2.0 / 3.0;
    const std::size_t Kbuf = 10;
    const std::vector<double> pi = mm1k_pi(rho, Kbuf);
    double L = 0.0, cost = 0.0;
    for (std::size_t n = 0; n <= Kbuf; ++n) {
        L += static_cast<double>(n) * pi[n];
        cost += static_cast<double>(n) * static_cast<double>(n) * pi[n];
    }
    compare("QueueLength", R[0], L);
    compare("Utilization", R[1], 1.0 - pi[0]);
    compare("BlockingProb", R[2], pi[Kbuf]);
    compare("QueueCost", R[3], cost);
}

/** `reward_model_templates.py`: the three reward templates on an M/M/1/10. */
void reward_model_templates() {
    std::size_t source = 0, queue = 0, sink = 0;
    Net m = reward_model("RewardTemplatesExample", SchedStrategy::FCFS, 1.0, 10.0, source, queue,
                         sink);
    OpenClass oclass(m, "Class1");
    m.set_arrival(source, oclass, Exp(1.5));
    m.set_service(queue, oclass, Exp(2.0));
    Routing P;
    serial(P, oclass, {source, queue, sink});
    m.link(P);

    const std::size_t ist = m.station_index(queue);
    const std::size_t K = 1;
    m.set_reward("QueueLength", reward_qlen(ist, K));
    m.set_reward("QueueLength_Class1", reward_qlen(ist, K, 1));
    m.set_reward("Utilization", reward_util(ist, K, 1.0));
    m.set_reward("Utilization_Class1", reward_util(ist, K, 1, 1.0));
    m.set_reward("BlockingProb", reward_blocking(ist, K, 10.0));

    const Sn& sn = m.get_struct();
    const ctmc::CtmcOptions opt;
    std::vector<std::string> names;
    const std::vector<double> R = ctmc::solver_ctmc_avg_reward(sn, opt, &names);

    section("CTMC");
    print_rewards(names, R);

    section("M/M/1/K analytical");
    const double rho = 1.5 / 2.0;
    const std::size_t Kbuf = 10;
    const std::vector<double> pi = mm1k_pi(rho, Kbuf);
    double L = 0.0;
    for (std::size_t n = 0; n <= Kbuf; ++n) L += static_cast<double>(n) * pi[n];
    compare("QueueLength", R[0], L);
    compare("Utilization", R[2], 1.0 - pi[0]);
    compare("BlockingProb", R[4], pi[Kbuf]);
}

/** `reward_model_aggregation.py`: the aggregation operations of a state view. */
void reward_model_aggregation() {
    std::size_t source = 0, queue = 0, sink = 0;
    Net m = reward_model("RewardAggregationExample", SchedStrategy::FCFS, 1.0, 4.0, source, queue,
                         sink);
    OpenClass c1(m, "HighPriority");
    OpenClass c2(m, "LowPriority");
    m.set_arrival(source, c1, Exp(1.0));
    m.set_arrival(source, c2, Exp(0.8));
    m.set_service(queue, c1, Exp(3.0));
    m.set_service(queue, c2, Exp(3.0));
    Routing P;
    serial(P, c1, {source, queue, sink});
    serial(P, c2, {source, queue, sink});
    m.link(P);

    const std::size_t ist = m.station_index(queue);
    const std::size_t K = 2;
    m.set_reward("TotalJobs", reward_qlen(ist, K));
    m.set_reward("MaxClass", [ist, K](const std::vector<double>& s) {
        return std::max(at(s, ist, K, 1), at(s, ist, K, 2));
    });
    m.set_reward("ClassCount", [ist, K](const std::vector<double>& s) {
        return (at(s, ist, K, 1) > 0 ? 1.0 : 0.0) + (at(s, ist, K, 2) > 0 ? 1.0 : 0.0);
    });
    m.set_reward("HP_Jobs", reward_qlen(ist, K, 1));
    m.set_reward("LP_Jobs", reward_qlen(ist, K, 2));
    m.set_reward("WeightedLoad", [ist, K](const std::vector<double>& s) {
        return 2.0 * at(s, ist, K, 1) + 1.0 * at(s, ist, K, 2);
    });

    const Sn& sn = m.get_struct();
    const ctmc::CtmcOptions opt;
    std::vector<std::string> names;
    const std::vector<double> R = ctmc::solver_ctmc_avg_reward(sn, opt, &names);

    section("CTMC");
    print_rewards(names, R);
}

/** `reward_model_multiclass.py`: per-class rewards at a two-server PS station. */
void reward_model_multiclass() {
    std::size_t source = 0, queue = 0, sink = 0;
    Net m = reward_model("MultiClassRewardExample", SchedStrategy::PS, 2.0, 6.0, source, queue,
                         sink);
    OpenClass ci(m, "Interactive");
    OpenClass cb(m, "Batch");
    m.set_arrival(source, ci, Exp(2.0));
    m.set_arrival(source, cb, Exp(1.5));
    m.set_service(queue, ci, Exp(0.5));
    m.set_service(queue, cb, Exp(1.0));
    Routing P;
    serial(P, ci, {source, queue, sink});
    serial(P, cb, {source, queue, sink});
    m.link(P);

    const std::size_t ist = m.station_index(queue);
    const std::size_t K = 2;
    m.set_reward("Interactive_QLen", reward_qlen(ist, K, 1));
    m.set_reward("Interactive_Util", reward_util(ist, K, 1, 2.0));
    m.set_reward("Batch_QLen", reward_qlen(ist, K, 2));
    m.set_reward("Batch_Util", reward_util(ist, K, 2, 2.0));
    m.set_reward("Total_QLen", reward_qlen(ist, K));
    m.set_reward("Total_Util", reward_util(ist, K, 2.0));
    m.set_reward("Interactive_Ratio", [ist, K](const std::vector<double>& s) {
        return at(s, ist, K, 1) / std::max(at(s, ist, K, 2), 0.001);
    });
    m.set_reward("Weighted_Cost", [ist, K](const std::vector<double>& s) {
        return 3.0 * at(s, ist, K, 1) + 1.0 * at(s, ist, K, 2);
    });
    m.set_reward("Fairness", [ist, K](const std::vector<double>& s) {
        return (at(s, ist, K, 1) > 0 && at(s, ist, K, 2) > 0) ? 1.0 : 0.0;
    });

    const Sn& sn = m.get_struct();
    ctmc::CtmcOptions opt;
    opt.cutoff = 6.0;
    std::vector<std::string> names;
    const std::vector<double> R = ctmc::solver_ctmc_avg_reward(sn, opt, &names);

    section("CTMC");
    print_rewards(names, R);

    section("Analysis");
    const double lambda_int = 2.0, mu_int = 2.0, lambda_batch = 1.5, mu_batch = 1.0, c = 2.0;
    const double rho_total = lambda_int / (c * mu_int) + lambda_batch / (c * mu_batch);
    kv("Servers", c);
    kv("Total utilization", rho_total);
    if (R[4] > 0) {
        kv("Interactive share (%)", 100.0 * R[0] / R[4]);
        kv("Batch share (%)", 100.0 * R[2] / R[4]);
    }
    kv("Interactive resp. time", R[0] / lambda_int);
    kv("Batch resp. time", R[2] / lambda_batch);
}

// The `rewardModel_*.py` spellings are symlinks to the `reward_model_*.py`
// files above, so each delegates to the one it points at.
void rewardModel_mm1k() { reward_model_mm1k(); }
void rewardModel_templates() { reward_model_templates(); }
void rewardModel_aggregation() { reward_model_aggregation(); }
void rewardModel_multiclass() { reward_model_multiclass(); }

LINE_EXAMPLE("advanced/rewardModel", reward_model_mm1k);
LINE_EXAMPLE("advanced/rewardModel", reward_model_templates);
LINE_EXAMPLE("advanced/rewardModel", reward_model_aggregation);
LINE_EXAMPLE("advanced/rewardModel", reward_model_multiclass);
LINE_EXAMPLE("advanced/rewardModel", rewardModel_mm1k);
LINE_EXAMPLE("advanced/rewardModel", rewardModel_templates);
LINE_EXAMPLE("advanced/rewardModel", rewardModel_aggregation);
LINE_EXAMPLE("advanced/rewardModel", rewardModel_multiclass);

}  // namespace examples
}  // namespace line
