/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_GEN_LAYERED_NETWORK_GENERATOR_H
#define LINE_GEN_LAYERED_NETWORK_GENERATOR_H

/**
 * @file
 * @ingroup line_gen
 * Random layered-queueing-network generation: the C++ twin of MATLAB
 * `@LayeredNetworkGenerator`, the JAR `jline.gen.LayeredNetworkGenerator` and
 * the native Python `line_solver.gen.layered_network_generator`.
 *
 * WHAT IT PRODUCES. A layered model of `numClients` reference tasks calling
 * down through `numLevels` layers of `numTasks` server tasks, hosted on
 * `numProcessors` processors. Every client is a REF task on its own
 * infinite-server processor with a think-time activity; every server task has
 * one entry and one activity that replies to it; the layers are wired top-down
 * so a task of level l is called by a task of level l-1, and the tasks are
 * spread over the processors in contiguous blocks.
 *
 * THIS ONE IS SAMPLE-PATH IDENTICAL TO THE JAR, unlike its flat twin
 * `NetworkGenerator`. Every draw the reference makes goes through the ONE
 * `java.util.Random` that `setSeed` replaces -- there is no `randGraph` and no
 * `Collections.shuffle` on a separate source here -- and `rng::JavaRandom` is a
 * bit-exact reproduction of that generator. The draw ORDER is preserved too:
 * the builder needs a processor before the task that sits on it, and a task
 * before its entry, where the reference creates the objects first and wires
 * them afterwards, so the PARAMETERS are drawn in the reference's order into
 * vectors and the builder calls are issued from those. Nothing reorders a draw.
 *
 * RANGES. Each `*_range` is a closed interval, sampled as the reference samples
 * it: an INTEGER range rounds its bounds inwards (`ceil` of the lower, `floor`
 * of the upper) and draws uniformly over the integers between, a REAL range
 * draws `lo + (hi - lo) * U`. The defaults are all `{1, 1}` with both infinity
 * probabilities zero, which generates the deterministic skeleton -- one job per
 * client, unit think times, unit demands, single-server tasks and processors.
 *
 * MEANS BECOME DISTRIBUTIONS the way `setHostDemand` does: a mean at or below
 * `FineTol` is an `Immediate`, anything else an `Exp` of that mean.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/util/error.h"
#include "line/util/rng_ssj.h"

namespace line {
namespace gen {

/** A closed interval the generator samples from. */
struct Range {
    double lo = 1.0;
    double hi = 1.0;
    Range() {}
    Range(double a, double b) : lo(a), hi(b) {}
};

/**
 * A random `lqn::LqnBuilder<T>` source, configured once and then drawn from.
 *
 * `generate` returns the BUILDER rather than the finalised struct: a caller
 * that wants the model `SolverLN` consumes calls `build()` on it, and one that
 wants the `.lqnx` document calls `lqn::write_lqnx(g.model(), ...)` or
 * `lqn::lqnx_to_string(g.model(), ...)`. Handing back the struct alone would
 * lose the raw model the writer needs.
 */
template <class T>
class LayeredNetworkGenerator {
  public:
    LayeredNetworkGenerator() : rng_(0) {}
    explicit LayeredNetworkGenerator(long long seed) : rng_(seed) {}

    /** Reseed the single stream every draw comes from. */
    void set_seed(long long seed) { rng_.set_seed(seed); }

    /** The `name` of the generated model. `lnw` is the reference's. */
    void set_model_name(const std::string& nm) { model_name_ = nm; }
    const std::string& model_name() const { return model_name_; }

    /** Reference-task populations. Strictly positive: a client with no job runs nothing. */
    void set_population_range(const Range& r) {
        population_range_ = check_range(r, "Population range", true);
    }
    /** Client think times. Zero is allowed and means an Immediate think. */
    void set_think_time_range(const Range& r) {
        think_time_range_ = check_range(r, "Think time range", false);
    }
    /** Probability that a server task is infinite-server rather than FCFS. */
    void set_task_inf_probability(double p) {
        task_inf_probability_ = check_probability(p, "Task infinite probability");
    }
    /** Probability that a processor is infinite-server rather than PS. */
    void set_proc_inf_probability(double p) {
        proc_inf_probability_ = check_probability(p, "Processor infinite probability");
    }
    /** Task multiplicities, for the tasks that are not infinite-server. */
    void set_task_multi_range(const Range& r) {
        task_multi_range_ = check_range(r, "Task multiplicity range", true);
    }
    /** Processor multiplicities, for the processors that are not infinite-server. */
    void set_proc_multi_range(const Range& r) {
        proc_multi_range_ = check_range(r, "Processor multiplicity range", true);
    }
    /** Activity host demands. Zero is allowed and means an Immediate activity. */
    void set_host_demand_range(const Range& r) {
        host_demand_range_ = check_range(r, "Host demand range", false);
    }
    /** Mean number of synchronous calls on a call arc. Strictly positive. */
    void set_synch_call_range(const Range& r) {
        synch_call_range_ = check_range(r, "Synchronous call range", true);
    }

    const Range& population_range() const { return population_range_; }
    const Range& think_time_range() const { return think_time_range_; }
    double task_inf_probability() const { return task_inf_probability_; }
    double proc_inf_probability() const { return proc_inf_probability_; }
    const Range& task_multi_range() const { return task_multi_range_; }
    const Range& proc_multi_range() const { return proc_multi_range_; }
    const Range& host_demand_range() const { return host_demand_range_; }
    const Range& synch_call_range() const { return synch_call_range_; }

    /** How many tasks landed on each level, after the last `generate`. */
    const std::vector<int>& tasks_per_level() const { return tasks_per_level_; }
    /** How many tasks landed on each processor, after the last `generate`. */
    const std::vector<int>& tasks_per_processor() const { return tasks_per_processor_; }

    // -----------------------------------------------------------------------
    // generate
    // -----------------------------------------------------------------------

    lqn::LqnBuilder<T> generate(int num_clients, int num_levels, int num_tasks,
                                int num_processors) {
        validate_args(num_clients, num_levels, num_tasks, num_processors);

        lqn::LqnBuilder<T> b;
        create_clients(b, num_clients);

        // THE PARAMETERS ARE DRAWN HERE, in the reference's order, and the
        // builder calls are issued below: `LqnBuilder::task` names its processor
        // at creation, while the reference assigns processors only after every
        // task and processor exists. Drawing first is what keeps the stream
        // identical to the reference's despite the different construction order.
        std::vector<TaskSpec> tspec;
        tspec.reserve(static_cast<std::size_t>(num_tasks));
        for (int t = 0; t < num_tasks; ++t) {
            TaskSpec s;
            if (choose_boolean(task_inf_probability_)) {
                s.mult = std::numeric_limits<double>::infinity();
                s.sched = lang::SchedStrategy::INF;
            } else {
                s.mult = double(sample_integer(task_multi_range_));
                s.sched = lang::SchedStrategy::FCFS;
            }
            s.hostdem = as_distribution(sample_real(host_demand_range_));
            tspec.push_back(s);
        }

        std::vector<TaskSpec> pspec;
        pspec.reserve(static_cast<std::size_t>(num_processors));
        for (int p = 0; p < num_processors; ++p) {
            TaskSpec s;
            if (choose_boolean(proc_inf_probability_)) {
                s.mult = std::numeric_limits<double>::infinity();
                s.sched = lang::SchedStrategy::INF;
            } else {
                s.mult = double(sample_integer(proc_multi_range_));
                s.sched = lang::SchedStrategy::PS;
            }
            pspec.push_back(s);
        }

        tasks_per_level_ = make_integer_vector(num_levels, num_tasks);
        tasks_per_processor_ = make_integer_vector(num_processors, num_tasks);

        for (int p = 0; p < num_processors; ++p)
            b.processor(proc_name(p), pspec[std::size_t(p)].mult, pspec[std::size_t(p)].sched);

        // The tasks take the processors in contiguous blocks, which is exactly
        // what `connectTasksToProcessors` does after the fact; it makes no draw,
        // so folding it into the creation loop changes nothing but the order of
        // two assignments.
        std::vector<std::string> host(static_cast<std::size_t>(num_tasks));
        {
            int filled = 0;
            for (int p = 0; p < num_processors; ++p)
                for (int k = 0; k < tasks_per_processor_[std::size_t(p)]; ++k)
                    host[std::size_t(filled++)] = proc_name(p);
        }
        for (int t = 0; t < num_tasks; ++t) {
            const TaskSpec& s = tspec[std::size_t(t)];
            b.task(task_name(t), s.mult, s.sched, host[std::size_t(t)]);
            b.entry(entry_name(t), task_name(t));
            b.activity(act_name(t), s.hostdem, task_name(t));
            b.bound_to(act_name(t), entry_name(t));
            b.replies_to(act_name(t), entry_name(t));
        }

        connect_clients_to_tasks(b, num_clients);
        connect_tasks_to_tasks(b, num_levels);
        return b;
    }

  private:
    struct TaskSpec {
        double mult = 1.0;
        lang::SchedStrategy sched = lang::SchedStrategy::FCFS;
        lang::Distrib<T> hostdem;
    };

    // -----------------------------------------------------------------------
    // Names
    // -----------------------------------------------------------------------

    static std::string cproc_name(int c) { return "c_processor_" + std::to_string(c + 1); }
    static std::string ctask_name(int c) { return "c_task_" + std::to_string(c + 1); }
    static std::string centry_name(int c) { return "c_entry_" + std::to_string(c + 1); }
    static std::string cact_name(int c) { return "c_activity_" + std::to_string(c + 1); }
    static std::string proc_name(int p) { return "processor_" + std::to_string(p + 1); }
    static std::string task_name(int t) { return "task_" + std::to_string(t + 1); }
    static std::string entry_name(int t) { return "entry_" + std::to_string(t + 1); }
    static std::string act_name(int t) { return "activity_" + std::to_string(t + 1); }

    // -----------------------------------------------------------------------
    // Validation
    // -----------------------------------------------------------------------

    static Range check_range(const Range& r, const std::string& name, bool strictly_positive) {
        const bool lower_ok = strictly_positive ? r.lo > 0.0 : r.lo >= 0.0;
        if (!(lower_ok && r.lo <= r.hi))
            throw InputError("LayeredNetworkGenerator: " + name + " is not valid");
        return r;
    }

    static double check_probability(double v, const std::string& name) {
        if (!(v >= 0.0 && v <= 1.0))
            throw InputError("LayeredNetworkGenerator: " + name + " is not valid");
        return v;
    }

    static void validate_args(int clients, int levels, int tasks, int procs) {
        if (clients < 1) throw InputError("LayeredNetworkGenerator: the number of clients is less "
                                          "than one");
        if (levels < 1) throw InputError("LayeredNetworkGenerator: the number of levels is less "
                                         "than one");
        if (tasks < 1) throw InputError("LayeredNetworkGenerator: the number of tasks is less than "
                                        "one");
        if (procs < 1) throw InputError("LayeredNetworkGenerator: the number of processors is less "
                                        "than one");
        if (levels > tasks)
            throw InputError("LayeredNetworkGenerator: the number of levels is greater than that "
                             "of tasks");
        if (procs > tasks)
            throw InputError("LayeredNetworkGenerator: the number of processors is greater than "
                             "that of tasks");
    }

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    void create_clients(lqn::LqnBuilder<T>& b, int num_clients) {
        for (int c = 0; c < num_clients; ++c) {
            const int population = sample_integer(population_range_);
            const double think = sample_real(think_time_range_);
            b.processor(cproc_name(c), std::numeric_limits<double>::infinity(),
                        lang::SchedStrategy::INF);
            b.task(ctask_name(c), double(population), lang::SchedStrategy::REF, cproc_name(c));
            b.entry(centry_name(c), ctask_name(c));
            // THE THINK TIME IS THE CLIENT ACTIVITY'S HOST DEMAND, not
            // `Task.setThinkTime`. That is what the reference builds, and the
            // two are not the same model: a think time is served by nothing,
            // while this demand is served by the client's own INF processor.
            b.activity(cact_name(c), as_distribution(think), ctask_name(c));
            b.bound_to(cact_name(c), centry_name(c));
        }
    }

    /**
     * Every first-level task is called by some client, and every client calls
     * some first-level task. The first loop guarantees the former, the second
     * sweeps up the clients the draws happened to miss.
     */
    void connect_clients_to_tasks(lqn::LqnBuilder<T>& b, int num_clients) {
        std::vector<bool> connected(static_cast<std::size_t>(num_clients), false);
        for (int t = 0; t < tasks_per_level_[0]; ++t) {
            const double calls = sample_real(synch_call_range_);
            const int c = sample_integer(Range(1.0, double(num_clients))) - 1;
            b.sync_call(cact_name(c), entry_name(t), num_traits<T>::from_double(calls));
            connected[std::size_t(c)] = true;
        }
        for (int c = 0; c < num_clients; ++c) {
            if (connected[std::size_t(c)]) continue;
            const double calls = sample_real(synch_call_range_);
            const int t = sample_integer(Range(1.0, double(tasks_per_level_[0]))) - 1;
            b.sync_call(cact_name(c), entry_name(t), num_traits<T>::from_double(calls));
            connected[std::size_t(c)] = true;
        }
    }

    /** Each task of level l is called by one task drawn from level l-1. */
    void connect_tasks_to_tasks(lqn::LqnBuilder<T>& b, int num_levels) {
        int seen = tasks_per_level_[0];
        for (int l = 1; l < num_levels; ++l) {
            for (int t2 = seen; t2 < seen + tasks_per_level_[std::size_t(l)]; ++t2) {
                const double calls = sample_real(synch_call_range_);
                const int t1 =
                    sample_integer(Range(double(seen - tasks_per_level_[std::size_t(l - 1)] + 1),
                                         double(seen))) - 1;
                b.sync_call(act_name(t1), entry_name(t2), num_traits<T>::from_double(calls));
            }
            seen += tasks_per_level_[std::size_t(l)];
        }
    }

    // -----------------------------------------------------------------------
    // The draws
    // -----------------------------------------------------------------------

    /** `randi([ceil(lo), floor(hi)])`: the bounds are rounded INWARDS. */
    int sample_integer(const Range& r) {
        const int lower = static_cast<int>(std::ceil(r.lo));
        const int upper = static_cast<int>(std::floor(r.hi));
        if (upper < lower)
            throw InputError("LayeredNetworkGenerator: an integer range rounds inwards to an "
                             "empty interval");
        return lower + rng_.next_int(upper - lower + 1);
    }

    double sample_real(const Range& r) { return r.lo + (r.hi - r.lo) * rng_.next_double(); }

    bool choose_boolean(double probability) { return rng_.next_double() < probability; }

    /**
     * `length` strictly positive integers summing to `sum`: one each, then the
     * surplus handed out one unit at a time to a uniformly drawn element. Not
     * the same law as `randintfixedsum` (this one is multinomial about the mean,
     * that one is uniform over the compositions), and the reference uses each
     * where it uses it.
     */
    std::vector<int> make_integer_vector(int length, int sum) {
        if (length < 1) throw InputError("LayeredNetworkGenerator: makeIntegerVector needs a "
                                         "positive length");
        if (sum < length)
            throw InputError("LayeredNetworkGenerator: makeIntegerVector cannot reach a sum below "
                             "its length");
        std::vector<int> v(static_cast<std::size_t>(length), 1);
        for (int s = 0; s < sum - length; ++s) ++v[std::size_t(rng_.next_int(length))];
        return v;
    }

    /** `setHostDemand(mean)`: an Immediate below FineTol, an Exp of that mean above. */
    static lang::Distrib<T> as_distribution(double mean) {
        if (mean <= lang::GlobalConstants::FineTol) return lang::Distrib<T>::immediate();
        return lang::Distrib<T>::exp_rate(num_traits<T>::from_double(1.0 / mean));
    }

    // -----------------------------------------------------------------------
    // State
    // -----------------------------------------------------------------------

    rng::JavaRandom rng_;
    std::string model_name_ = "lnw";
    Range population_range_ = Range(1.0, 1.0);
    Range think_time_range_ = Range(1.0, 1.0);
    Range task_multi_range_ = Range(1.0, 1.0);
    Range proc_multi_range_ = Range(1.0, 1.0);
    Range host_demand_range_ = Range(1.0, 1.0);
    Range synch_call_range_ = Range(1.0, 1.0);
    double task_inf_probability_ = 0.0;
    double proc_inf_probability_ = 0.0;

    std::vector<int> tasks_per_level_;
    std::vector<int> tasks_per_processor_;
};

}  // namespace gen
}  // namespace line

#endif  // LINE_GEN_LAYERED_NETWORK_GENERATOR_H
