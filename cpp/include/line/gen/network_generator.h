/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_GEN_NETWORK_GENERATOR_H
#define LINE_GEN_NETWORK_GENERATOR_H

/**
 * @file
 * @ingroup line_gen
 * Random queueing-network generation: the C++ twin of MATLAB
 * `@NetworkGenerator`, the JAR `jline.gen.NetworkGenerator` and the native
 * Python `line_solver.gen.network_generator`.
 *
 * WHAT IT PRODUCES. A `qn::Network<T>` with `numQueues` queues, `numDelays`
 * delay stations, `numOClass` open classes and `numCClass` closed classes, wired
 * over a random strongly connected topology, with random service laws, random
 * scheduling, random routing and (optionally) random ClassSwitch nodes on the
 * links. It is the model source the test suites, the benchmark sweeps and the
 * solver-comparison harnesses draw from, so the shapes it can emit matter more
 * than any single draw: every arm of the reference is reproduced, including the
 * multi-chain class-switch masks and the load bands.
 *
 * HOW IT DIFFERS FROM THE REFERENCES, and this is deliberate:
 *
 * 1. ONE SEEDED STREAM. MATLAB, the JAR and Python all leave part of the draw
 *    on a source `setSeed` does not reach -- the JAR's `randGraph` builds its
 *    own `new Random()` and `Collections.shuffle` uses the shared static one,
 *    and MATLAB's `randGraph` draws from the global stream while the object
 *    seeds nothing at all. A seeded run there is reproducible in its service
 *    laws and not in its topology. Here EVERY draw, the spanning tree and both
 *    permutations included, comes from one `rng::JavaRandom`, so `set_seed(s)`
 *    reproduces the whole model. The consequence is that this port is NOT
 *    sample-path identical to the JAR or MATLAB at a shared seed; it is a
 *    DISTRIBUTIONAL port, and the parity claim is over the family of models it
 *    can emit, not over one draw. Do not baseline a seeded golden across the
 *    codebases from this generator.
 *
 * 2. LINKS ARE A ROUTING MATRIX, NOT `addLink`. The reference calls
 *    `model.addLink(a, b)` and then `setProbRouting`/`setRouting`; this port has
 *    no `addLink` -- `qn::Network::link(P)` derives the connection graph from
 *    the routing matrix itself (see `network_struct.h`, the RAND expansion).
 *    So a RANDOM-routed (node, class) pair is given its arcs in `P` as well as
 *    its strategy: the arcs are what make the pair connected, and `route_eff`
 *    then overwrites the probabilities with the uniform split RAND means. The
 *    two spellings describe the same model.
 *
 * 3. THE REJECTION TEST IS STRUCTURAL. The JAR resamples when `sn_refresh_visits`
 *    throws a message containing "no recurrent flow"; no such throw exists any
 *    more in any codebase, so that gate is dead code there. This port tests what
 *    the gate was meant to test: after the struct is materialised, EVERY CLOSED
 *    CHAIN must actually visit its own reference station. A closed chain
 *    stranded in a component that does not contain its reference station gets
 *    zero visits there, its visit vector is not normalisable, and the model is
 *    not a valid closed network -- such a draw is discarded and resampled, up to
 *    `max_generate_attempts()`.
 *
 * 4. `initializeStates` IS NOT A FIELD HERE. Its only effect in the JAR is to
 *    force `getStruct()` so the rejection test has something to test; this port
 *    materialises the struct unconditionally for exactly that reason. The C++
 *    model has no per-node `setState` to fill in either: an initial state is
 *    declared through `Network::set_state_prior`, over a state SPACE rather
 *    than one row, so there is nothing here to port the flag onto.
 *
 * ARITHMETIC. Every random quantity is drawn as a double (the reference draws
 * are `nextDouble`/`nextInt`) and lifted into T, so the generator instantiates
 * at exact arithmetic as well. The one transcendental step is the HyperExp
 * moment fit, which is done in double and its three parameters lifted, the same
 * bargain `sn_aggregate_chains` strikes for the same reason.
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "line/lang/dist_fitters.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/rng_ssj.h"

namespace line {
namespace gen {

/** The two topology generators the reference ships, `randGraph` and `cyclicGraph`. */
enum class TopologyKind { Rand, Cyclic };

// ---------------------------------------------------------------------------
// Topology generators (module-level in Python, static in the JAR)
// ---------------------------------------------------------------------------

/**
 * `randSpanningTree(n)`: vertex i (i >= 1) is attached to a uniformly chosen
 * earlier vertex, which is a uniform draw over the labelled rooted trees this
 * construction reaches, and always a tree rooted at 0.
 */
inline Matrix<double> rand_spanning_tree(std::size_t num_vertices, rng::JavaRandom& r) {
    Matrix<double> tree(num_vertices, num_vertices, 0.0);
    for (std::size_t i = 1; i < num_vertices; ++i)
        tree(static_cast<std::size_t>(r.next_int(static_cast<int32_t>(i))), i) = 1.0;
    return tree;
}

namespace detail {

/** The DFS bookkeeping of `strongConnect`, Tarjan's low-link walk. */
struct RandGraphState {
    int global_start_time = 0;
    std::vector<int> start_time, lowest_link, inv_start_time;
    explicit RandGraphState(std::size_t n)
        : start_time(n, -1), lowest_link(n, -1), inv_start_time(n, -1) {}
};

/**
 * The repair half of `randGraph`: a DFS that, on leaving the root of a strongly
 * connected component that is not the root of the whole graph, adds ONE back
 * edge from a random descendant to a random ancestor. That is what turns the
 * spanning tree into a strongly connected digraph without making it complete.
 */
inline void strong_connect(Matrix<double>& g, std::size_t v, RandGraphState& st,
                           rng::JavaRandom& r) {
    st.start_time[v] = st.global_start_time;
    st.inv_start_time[st.start_time[v]] = static_cast<int>(v);
    st.lowest_link[v] = st.start_time[v];
    ++st.global_start_time;

    std::vector<std::size_t> out;
    for (std::size_t w = 0; w < g.cols(); ++w)
        if (g(v, w) > 0.0) out.push_back(w);

    for (std::size_t k = 0; k < out.size(); ++k) {
        const std::size_t w = out[k];
        if (st.start_time[w] == -1) {
            strong_connect(g, w, st, r);
            st.lowest_link[v] = std::min(st.lowest_link[v], st.lowest_link[w]);
        } else {
            st.lowest_link[v] = std::min(st.lowest_link[v], st.start_time[w]);
        }
    }

    if (st.lowest_link[v] == st.start_time[v] && st.start_time[v] > 0) {
        const int descendant = r.next_int(st.global_start_time - st.start_time[v]) + st.start_time[v];
        const int ancestor = r.next_int(st.start_time[v]);
        g(static_cast<std::size_t>(st.inv_start_time[descendant]),
          static_cast<std::size_t>(st.inv_start_time[ancestor])) = 1.0;
        st.lowest_link[v] = ancestor;
    }
}

/**
 * `java.util.Collections.shuffle(list, rnd)`, the Fisher-Yates sweep it
 * specifies: for i from n-1 down to 1, swap element i with a uniform element of
 * 0..i. Reproduced rather than replaced by `std::shuffle` so the permutation is
 * the reference's function of the stream.
 */
inline void java_shuffle(std::vector<std::size_t>& v, rng::JavaRandom& r) {
    for (std::size_t i = v.size(); i > 1; --i)
        std::swap(v[i - 1], v[static_cast<std::size_t>(r.next_int(static_cast<int32_t>(i)))]);
}

}  // namespace detail

/**
 * `randGraph(n)`: a random strongly connected digraph on n vertices, as a
 * random spanning tree closed up by `strongConnect` and then relabelled by a
 * random permutation (without which vertex 0 would always be the DFS root).
 */
inline Matrix<double> rand_graph(std::size_t num_vertices, rng::JavaRandom& r) {
    if (num_vertices == 0) throw InputError("randGraph: the number of vertices must be positive");
    if (num_vertices == 1) {
        Matrix<double> adj(1, 1, 0.0);
        adj(0, 0) = 1.0;
        return adj;
    }
    Matrix<double> g = rand_spanning_tree(num_vertices, r);
    detail::RandGraphState st(num_vertices);
    detail::strong_connect(g, 0, st, r);

    std::vector<std::size_t> perm(num_vertices);
    for (std::size_t i = 0; i < num_vertices; ++i) perm[i] = i;
    detail::java_shuffle(perm, r);

    Matrix<double> out(num_vertices, num_vertices, 0.0);
    for (std::size_t i = 0; i < num_vertices; ++i)
        for (std::size_t j = 0; j < num_vertices; ++j)
            if (g(i, j) > 0.0) out(perm[i], perm[j]) = 1.0;
    return out;
}

/** `cyclicGraph(n)`: the single cycle 0 -> 1 -> ... -> n-1 -> 0. */
inline Matrix<double> cyclic_graph(std::size_t num_vertices) {
    if (num_vertices == 0) throw InputError("cyclicGraph: the number of vertices must be positive");
    Matrix<double> adj(num_vertices, num_vertices, 0.0);
    for (std::size_t i = 0; i + 1 < num_vertices; ++i) adj(i, i + 1) = 1.0;
    adj(num_vertices - 1, 0) = 1.0;
    return adj;
}

// ---------------------------------------------------------------------------
// The two fixed-sum samplers
// ---------------------------------------------------------------------------

/**
 * `randfixedsumone(n)`: n probabilities summing to exactly 1.
 *
 * The reference's own recipe, quirks kept: n uniforms are normalised, each is
 * rounded UP to a multiple of 1e-3 (so the rounded vector oversums), and the
 * LARGEST entry absorbs the whole residual. Rounding to three digits is what
 * keeps the emitted `model.json` readable; the absorption is what keeps the row
 * stochastic to the last bit after it.
 */
inline std::vector<double> randfixedsumone(std::size_t num_elems, rng::JavaRandom& r) {
    if (num_elems == 0) return std::vector<double>();
    if (num_elems == 1) return std::vector<double>(1, 1.0);

    std::vector<double> values(num_elems, 0.0);
    double sum = 0.0;
    for (std::size_t i = 0; i < num_elems; ++i) {
        values[i] = r.next_double();
        sum += values[i];
    }
    for (std::size_t i = 0; i < num_elems; ++i)
        values[i] = std::ceil(values[i] / sum * 1000.0) / 1000.0;

    std::size_t max_idx = 0;
    for (std::size_t i = 1; i < num_elems; ++i)
        if (values[i] > values[max_idx]) max_idx = i;
    double current = 0.0;
    for (std::size_t i = 0; i < num_elems; ++i) current += values[i];
    values[max_idx] -= (current - 1.0);
    return values;
}

/**
 * `randintfixedsum(s, n)`: n STRICTLY POSITIVE integers summing to s.
 *
 * The reference's recursion: draw the first part uniformly from 1..s-n (which
 * leaves at least one unit for each remaining part), recurse on the rest, and
 * shuffle so the first part carries no special distribution.
 */
inline std::vector<int> randintfixedsum(int s, int n, rng::JavaRandom& r) {
    if (n <= 0) throw InputError("randintfixedsum: the number of parts must be positive");
    if (s < n) throw InputError("randintfixedsum: the sum must be at least the number of parts");
    if (n == 1) return std::vector<int>(1, s);
    if (s == n) return std::vector<int>(static_cast<std::size_t>(n), 1);

    const int first = r.next_int(s - n) + 1;
    std::vector<int> out = randintfixedsum(s - first, n - 1, r);
    out.insert(out.begin(), first);
    std::vector<std::size_t> idx(out.size());
    for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = i;
    detail::java_shuffle(idx, r);
    std::vector<int> shuffled(out.size(), 0);
    for (std::size_t i = 0; i < idx.size(); ++i) shuffled[i] = out[idx[i]];
    return shuffled;
}

// ---------------------------------------------------------------------------
// The generator
// ---------------------------------------------------------------------------

/**
 * A random `qn::Network<T>` source, configured once and then drawn from.
 *
 * Every setter validates its argument against the same vocabulary the reference
 * accepts, so a typo is an error at configuration time rather than a silently
 * different model.
 */
template <class T>
class NetworkGenerator {
  public:
    NetworkGenerator() : rng_(0) {}
    explicit NetworkGenerator(long long seed) : rng_(seed) {}

    /** Reseed the single stream every draw comes from. */
    void set_seed(long long seed) { rng_.set_seed(seed); }

    /** The name given to the generated model. `nw` is the reference's. */
    void set_model_name(const std::string& nm) { model_name_ = nm; }
    const std::string& model_name() const { return model_name_; }

    /** fcfs, ps, inf, lcfs, lcfspr, siro, sjf, ljf, sept, lept, or randomize. */
    void set_sched_strat(const std::string& strat) {
        if (!is_one_of(strat, sched_vocabulary()))
            throw InputError("NetworkGenerator: scheduling strategy '" + strat +
                             "' does not exist or is not supported");
        sched_strat_ = strat;
    }
    const std::string& sched_strat() const { return sched_strat_; }

    /** Probabilities, Random, or randomize. */
    void set_routing_strat(const std::string& strat) {
        if (strat != "Probabilities" && strat != "Random" && strat != "randomize")
            throw InputError("NetworkGenerator: routing strategy '" + strat +
                             "' does not exist or is not supported");
        routing_strat_ = strat;
    }
    const std::string& routing_strat() const { return routing_strat_; }

    /** Exp, Erlang, HyperExp, or randomize. */
    void set_distribution(const std::string& d) {
        if (!ieq(d, "Exp") && !ieq(d, "Erlang") && !ieq(d, "HyperExp") && !ieq(d, "randomize"))
            throw InputError("NetworkGenerator: distribution '" + d +
                             "' does not exist or is not supported");
        distribution_ = d;
    }
    const std::string& distribution() const { return distribution_; }

    /** high, medium, low, or randomize: the population band of a closed class. */
    void set_cclass_job_load(const std::string& load) {
        if (!ieq(load, "high") && !ieq(load, "medium") && !ieq(load, "low") &&
            !ieq(load, "randomize"))
            throw InputError("NetworkGenerator: model load can only be high, medium, low or "
                             "randomize");
        cclass_job_load_ = load;
    }
    const std::string& cclass_job_load() const { return cclass_job_load_; }

    /** Service means spread over 2^-6 .. 2^6 instead of all being 1. */
    void set_varying_service_rates(bool v) { varying_service_rates_ = v; }
    bool varying_service_rates() const { return varying_service_rates_; }

    /** Queues get 1..40 servers instead of exactly one. */
    void set_multi_server_queues(bool v) { multi_server_queues_ = v; }
    bool multi_server_queues() const { return multi_server_queues_; }

    /** Each link gets a ClassSwitch node inserted on a coin flip. */
    void set_random_cs_nodes(bool v) { random_cs_nodes_ = v; }
    bool random_cs_nodes() const { return random_cs_nodes_; }

    /**
     * Classes are partitioned into random chains and switching is confined to a
     * chain, instead of the default "all open classes switch among themselves
     * and all closed classes among themselves".
     */
    void set_multi_chain_cs(bool v) { multi_chain_cs_ = v; }
    bool multi_chain_cs() const { return multi_chain_cs_; }

    /** Pick one of the two shipped topology generators. */
    void set_topology(TopologyKind kind) {
        if (kind == TopologyKind::Cyclic)
            topology_fcn_ = [](std::size_t n, rng::JavaRandom&) { return cyclic_graph(n); };
        else
            topology_fcn_ = [](std::size_t n, rng::JavaRandom& r) { return rand_graph(n, r); };
    }

    /**
     * A custom topology: any function from a vertex count (and the generator's
     * own stream, so a custom topology stays reproducible too) to a square 0/1
     * adjacency matrix. Validated on a two-vertex call, as the reference does.
     */
    void set_topology_fcn(const std::function<Matrix<double>(std::size_t, rng::JavaRandom&)>& fcn) {
        if (!fcn) throw InputError("NetworkGenerator: the topology function is empty");
        rng::JavaRandom probe(0);
        const Matrix<double> adj = fcn(2, probe);
        if (adj.rows() != 2 || adj.cols() != 2)
            throw InputError("NetworkGenerator: topologyFcn must take a positive integer and "
                             "return a square adjacency matrix of that order");
        topology_fcn_ = fcn;
    }

    /** Resampling budget before an unsatisfiable configuration is reported. */
    static std::size_t max_generate_attempts() { return 100; }

    // -----------------------------------------------------------------------
    // generate
    // -----------------------------------------------------------------------

    /**
     * The full form. `num_delays < 0` means "decide for me", which is the
     * reference's `null`: one delay when there is a single queue, otherwise a
     * coin flip between none and one.
     */
    qn::Network<T> generate(int num_queues, int num_delays, int num_oclass, int num_cclass) {
        if (num_delays < 0) num_delays = (num_queues > 1) ? rng_.next_int(2) : 1;
        validate_args(num_queues, num_delays, num_oclass, num_cclass);

        std::string last_reject;
        for (std::size_t attempt = 0; attempt < max_generate_attempts(); ++attempt) {
            try {
                qn::Network<T> model(model_name_);
                build(model, num_queues, num_delays, num_oclass, num_cclass);
                const qn::NetworkStruct<T>& sn = model.get_struct();
                const std::string why = why_invalid(sn);
                if (!why.empty()) {
                    last_reject = why;
                    continue;  // resample: the draw is not a valid closed network
                }
                return model;
            } catch (const Error& e) {
                last_reject = e.what();
                continue;
            }
        }
        throw InputError("NetworkGenerator: could not generate a model with recurrent closed-chain "
                         "routing in " + std::to_string(max_generate_attempts()) +
                         " attempts. The requested topology may not admit a strongly connected "
                         "closed chain. Last rejection: " + last_reject);
    }

    /** `generate(numQueues, numDelays, numOClass)`, closed-class count drawn 1..4. */
    qn::Network<T> generate(int num_queues, int num_delays, int num_oclass) {
        return generate(num_queues, num_delays, num_oclass, rng_.next_int(4) + 1);
    }
    /** `generate(numQueues, numDelays)`, no open classes. */
    qn::Network<T> generate(int num_queues, int num_delays) {
        return generate(num_queues, num_delays, 0, rng_.next_int(4) + 1);
    }
    /** `generate(numQueues)`, delays decided by the rule above. */
    qn::Network<T> generate(int num_queues) {
        return generate(num_queues, -1, 0, rng_.next_int(4) + 1);
    }
    /** Everything drawn: 1..8 queues, 1..4 closed classes, no open classes. */
    qn::Network<T> generate() {
        const int nq = rng_.next_int(8) + 1;
        return generate(nq, -1, 0, rng_.next_int(4) + 1);
    }

  private:
    // -----------------------------------------------------------------------
    // Configuration
    // -----------------------------------------------------------------------

    static const std::vector<std::string>& sched_vocabulary() {
        static const std::vector<std::string> v = {"fcfs", "ps",   "inf",  "lcfs", "lcfspr", "siro",
                                                   "sjf",  "ljf",  "sept", "lept", "randomize"};
        return v;
    }

    static bool ieq(const std::string& a, const std::string& b) {
        if (a.size() != b.size()) return false;
        for (std::size_t i = 0; i < a.size(); ++i)
            if (std::tolower(static_cast<unsigned char>(a[i])) !=
                std::tolower(static_cast<unsigned char>(b[i])))
                return false;
        return true;
    }

    static bool is_one_of(const std::string& s, const std::vector<std::string>& v) {
        return std::find(v.begin(), v.end(), s) != v.end();
    }

    static void validate_args(int nq, int nd, int no, int nc) {
        if (nq < 0 || nd < 0 || no < 0 || nc < 0)
            throw InputError("NetworkGenerator: the station and class counts must be non-negative");
        if (nq + nd <= 0 || no + nc <= 0)
            throw InputError("NetworkGenerator: at least one station and one job class are "
                             "required");
    }

    // -----------------------------------------------------------------------
    // The draws
    // -----------------------------------------------------------------------

    double choose_num_servers() {
        return multi_server_queues_ ? double(rng_.next_int(MAX_SERVERS) + 1) : 1.0;
    }

    double choose_num_jobs() {
        if (ieq(cclass_job_load_, "high")) return band(HIGH_LO, HIGH_HI);
        if (ieq(cclass_job_load_, "medium")) return band(MED_LO, MED_HI);
        if (ieq(cclass_job_load_, "low")) return band(LOW_LO, LOW_HI);
        return double(rng_.next_int(HIGH_HI) + 1);  // randomize: the whole 1..40 range
    }

    double band(int lo, int hi) { return double(rng_.next_int(hi - lo + 1) + lo); }

    lang::SchedStrategy choose_sched_strat() {
        if (ieq(sched_strat_, "randomize"))
            return (rng_.next_int(2) + 1) == 1 ? lang::SchedStrategy::FCFS
                                               : lang::SchedStrategy::PS;
        return sched_from_name(sched_strat_);
    }

    static lang::SchedStrategy sched_from_name(const std::string& s) {
        if (s == "fcfs") return lang::SchedStrategy::FCFS;
        if (s == "ps") return lang::SchedStrategy::PS;
        if (s == "inf") return lang::SchedStrategy::INF;
        if (s == "lcfs") return lang::SchedStrategy::LCFS;
        if (s == "lcfspr") return lang::SchedStrategy::LCFSPR;
        if (s == "siro") return lang::SchedStrategy::SIRO;
        if (s == "sjf") return lang::SchedStrategy::SJF;
        if (s == "ljf") return lang::SchedStrategy::LJF;
        if (s == "sept") return lang::SchedStrategy::SEPT;
        if (s == "lept") return lang::SchedStrategy::LEPT;
        throw InputError("NetworkGenerator: scheduling strategy '" + s + "' is not supported");
    }

    /** "Random" or "Probabilities"; `randomize` picks between them per class. */
    std::string choose_routing_strat() {
        if (ieq(routing_strat_, "randomize"))
            return (rng_.next_int(2) + 1) == 1 ? std::string("Random") : std::string("Probabilities");
        return routing_strat_;
    }

    /**
     * A service or arrival law. The mean is 1 unless varying rates are asked
     * for, in which case it is a power of two in 2^-6 .. 2^6; the Erlang order
     * and the HyperExp SCV are powers of two in 1 .. 64, which is the reference's
     * way of covering both sides of the exponential in a few draws.
     */
    lang::Distrib<T> choose_distribution() {
        int id;
        if (ieq(distribution_, "Exp")) id = 1;
        else if (ieq(distribution_, "Erlang")) id = 2;
        else if (ieq(distribution_, "HyperExp")) id = 3;
        else id = rng_.next_int(3) + 1;

        const double mean =
            varying_service_rates_ ? std::pow(2.0, double(rng_.next_int(13) - 6)) : 1.0;

        if (id == 1) return lang::Distrib<T>::exp_rate(num_traits<T>::from_double(1.0 / mean));
        if (id == 2) {
            const std::size_t k = std::size_t(std::pow(2.0, double(rng_.next_int(7))));
            return lang::erlang_fit_mean_order(num_traits<T>::from_double(mean), k);
        }
        const double scv = std::pow(2.0, double(rng_.next_int(7)));
        return fit_hyperexp(mean, scv);
    }

    /**
     * `HyperExp.fitMeanAndSCV`, which is `map_hyperexp` at p = 0.99 read back as
     * (p, mu1, mu2). The fit needs a square root of the moment discriminant, so
     * at exact arithmetic it is done in double and the three parameters lifted:
     * a two-moment fit carries no exactness claim, and refusing would leave the
     * HyperExp arm with no exact-arithmetic path at all.
     */
    static lang::Distrib<T> fit_hyperexp(double mean, double scv) {
        const lang::Distrib<double> dd = lang::hyperexp_fit_mean_scv<double>(mean, scv);
        return lang::Distrib<T>::hyperexp(num_traits<T>::from_double(dd.params[0]),
                                          num_traits<T>::from_double(dd.params[1]),
                                          num_traits<T>::from_double(dd.params[2]));
    }

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    void build(qn::Network<T>& model, int num_queues, int num_delays, int num_oclass,
               int num_cclass) {
        stations_.clear();
        source_ = 0;
        sink_ = 0;

        create_stations(model, num_queues, num_delays, num_oclass > 0);
        create_classes(model, num_oclass, num_cclass);
        set_service_processes(model);
        define_topology(model, num_oclass > 0);
    }

    void create_stations(qn::Network<T>& model, int num_queues, int num_delays, bool has_oclass) {
        for (int i = 0; i < num_queues; ++i) {
            const lang::SchedStrategy sched = choose_sched_strat();
            const std::size_t nd = model.add_queue("queue" + std::to_string(i + 1), sched);
            // The draw is unconditional, as in the reference: an INF-scheduled
            // queue ignores the count but must not change the stream.
            const double servers = choose_num_servers();
            model.set_number_of_servers(nd, servers);
            stations_.push_back(nd);
        }
        for (int i = 0; i < num_delays; ++i)
            stations_.push_back(model.add_delay("delay" + std::to_string(i + 1)));
        if (has_oclass) {
            source_ = model.add_source("source");
            sink_ = model.add_sink("sink");
        }
    }

    void create_classes(qn::Network<T>& model, int num_oclass, int num_cclass) {
        classes_.clear();
        for (int i = 0; i < num_oclass; ++i) {
            const std::size_t r = model.add_open_class("OClass" + std::to_string(i + 1));
            model.set_arrival(source_, r, choose_distribution());
            classes_.push_back(r);
        }
        if (num_cclass > 0) {
            // THE REFERENCE STATION IS NEVER THE SOURCE. MATLAB draws over
            // `getNumberOfStations - (numOClass > 0)`, which excludes it; the
            // JAR draws over `getStations()`, which includes it and can hand a
            // closed class a Source reference station. MATLAB is ground truth.
            const std::size_t ref =
                stations_[std::size_t(rng_.next_int(static_cast<int32_t>(stations_.size())))];
            for (int i = 0; i < num_cclass; ++i)
                classes_.push_back(model.add_closed_class("CClass" + std::to_string(i + 1),
                                                          choose_num_jobs(), ref));
        }
    }

    void set_service_processes(qn::Network<T>& model) {
        for (std::size_t c = 0; c < classes_.size(); ++c)
            for (std::size_t s = 0; s < stations_.size(); ++s)
                model.set_service(stations_[s], classes_[c], choose_distribution());
    }

    /**
     * The topology, the ClassSwitch insertions and the routing, in the
     * reference's order: the class-switch mask is drawn BEFORE the Source and
     * Sink rows are grafted onto the adjacency matrix.
     */
    void define_topology(qn::Network<T>& model, bool has_oclass) {
        Matrix<double> topo = topology_fcn_(stations_.size(), rng_);
        if (topo.rows() != stations_.size() || topo.cols() != stations_.size())
            throw InputError("NetworkGenerator: the topology function returned an adjacency "
                             "matrix of the wrong order");

        const Matrix<T> mask = gen_cs_mask(model);

        // Stations own adjacency rows 0..S-1; the Source and Sink take the two
        // rows appended here, matching the node order the builder created.
        std::vector<std::size_t> row_node = stations_;
        if (has_oclass) {
            const std::size_t S = stations_.size();
            Matrix<double> grown(S + 2, S + 2, 0.0);
            for (std::size_t i = 0; i < S; ++i)
                for (std::size_t j = 0; j < S; ++j) grown(i, j) = topo(i, j);
            grown(S, std::size_t(rng_.next_int(static_cast<int32_t>(S)))) = 1.0;
            grown(std::size_t(rng_.next_int(static_cast<int32_t>(S))), S + 1) = 1.0;
            topo = grown;
            row_node.push_back(source_);
            row_node.push_back(sink_);
        }

        // The Sink is not a station and routes nothing, so it owns no row here:
        // `apply_sink_closure` writes the Sink -> Source arc during the refresh.
        const std::size_t nrows = has_oclass ? stations_.size() + 1 : stations_.size();

        qn::RoutingMatrix<T> P;
        for (std::size_t i = 0; i < nrows; ++i) {
            std::vector<std::size_t> dest;  // 1-based node indices, ascending
            for (std::size_t j = 0; j < topo.cols(); ++j)
                if (topo(i, j) > 0.0) dest.push_back(row_node[j]);
            std::sort(dest.begin(), dest.end());
            const std::vector<std::size_t> outgoing =
                add_outgoing_links(model, P, row_node[i], dest, mask);
            set_routing_strategies(model, P, row_node[i], outgoing);
        }
        model.link(P);
    }

    /**
     * `genCSMask`: which class may switch into which at an inserted ClassSwitch.
     *
     * The default confines switching to the open block and the closed block,
     * which is what keeps a closed chain closed. `multi_chain_cs` refines that
     * into randomly sized chains inside each block, so a model can carry several
     * independent chains of the same kind.
     */
    Matrix<T> gen_cs_mask(qn::Network<T>& model) {
        const qn::NetworkStruct<T>& sn = model.raw_struct();
        std::size_t num_open = 0, num_closed = 0;
        for (std::size_t r = 0; r < sn.classes.size(); ++r) {
            if (sn.classes[r].type == lang::JobClassType::OPEN) ++num_open;
            else ++num_closed;
        }
        const std::size_t K = num_open + num_closed;
        const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
        Matrix<T> mask(K, K, zero);

        if (!multi_chain_cs_) {
            for (std::size_t i = 0; i < num_open; ++i)
                for (std::size_t j = 0; j < num_open; ++j) mask(i, j) = one;
            for (std::size_t i = num_open; i < K; ++i)
                for (std::size_t j = num_open; j < K; ++j) mask(i, j) = one;
            return mask;
        }

        std::vector<int> sizes;
        if (num_open > 0) {
            const std::vector<int> s = randintfixedsum(
                static_cast<int>(num_open), rng_.next_int(static_cast<int32_t>(num_open)) + 1, rng_);
            sizes.insert(sizes.end(), s.begin(), s.end());
        }
        if (num_closed > 0) {
            const std::vector<int> s =
                randintfixedsum(static_cast<int>(num_closed),
                                rng_.next_int(static_cast<int32_t>(num_closed)) + 1, rng_);
            sizes.insert(sizes.end(), s.begin(), s.end());
        }
        std::size_t start = 0;
        for (std::size_t c = 0; c < sizes.size(); ++c) {
            const std::size_t end = start + std::size_t(sizes[c]);
            for (std::size_t i = start; i < end; ++i)
                for (std::size_t j = start; j < end; ++j) mask(i, j) = one;
            start = end;
        }
        return mask;
    }

    /** A ClassSwitch matrix: each row is a `randfixedsumone` over its mask. */
    Matrix<T> rand_class_switch_matrix(const Matrix<T>& mask) {
        const std::size_t K = mask.rows();
        const T zero = num_traits<T>::from_int(0);
        Matrix<T> C(K, K, zero);
        for (std::size_t i = 0; i < K; ++i) {
            std::vector<std::size_t> valid;
            for (std::size_t j = 0; j < K; ++j)
                if (num_traits<T>::to_double(mask(i, j)) > 0.0) valid.push_back(j);
            if (valid.empty()) {
                C(i, i) = num_traits<T>::from_int(1);  // a class no mask reaches keeps itself
                continue;
            }
            const std::vector<double> probs = randfixedsumone(valid.size(), rng_);
            for (std::size_t k = 0; k < valid.size(); ++k)
                C(i, valid[k]) = num_traits<T>::from_double(probs[k]);
        }
        return C;
    }

    /**
     * Walk the row's destinations, optionally interposing a ClassSwitch node.
     *
     * Returns the nodes the row actually routes INTO, which is the destination
     * itself or the interposed node; the second leg cs -> dest is deterministic
     * in the arriving class and is written here.
     */
    std::vector<std::size_t> add_outgoing_links(qn::Network<T>& model, qn::RoutingMatrix<T>& P,
                                                std::size_t from,
                                                const std::vector<std::size_t>& dest,
                                                const Matrix<T>& mask) {
        const qn::NetworkStruct<T>& sn = model.raw_struct();
        const std::size_t K = sn.classes.size();
        const T one = num_traits<T>::from_int(1);
        std::vector<std::size_t> outgoing;
        for (std::size_t d = 0; d < dest.size(); ++d) {
            const std::size_t to = dest[d];
            const bool from_source = sn.nodes[from - 1].nodetype == lang::NodeType::Source;
            const bool to_sink = sn.nodes[to - 1].nodetype == lang::NodeType::Sink;
            if (random_cs_nodes_ && rng_.next_boolean() && !from_source && !to_sink) {
                const Matrix<T> C = rand_class_switch_matrix(mask);
                const std::size_t cs = model.add_class_switch(
                    "cs_" + sn.nodes[from - 1].name + "_" + sn.nodes[to - 1].name, C);
                outgoing.push_back(cs);
                for (std::size_t r = 1; r <= K; ++r) P.set(r, r, cs, to, one);
            } else {
                outgoing.push_back(to);
            }
        }
        return outgoing;
    }

    /**
     * The per-class routing of one row.
     *
     * A Source routes each OPEN class into every outgoing node with probability
     * one (its row has a single destination by construction, so the row is
     * stochastic). Every other row draws a strategy per class: `Probabilities`
     * spreads a `randfixedsumone` over the destinations, `Random` declares RAND
     * and lets `route_eff` spread uniformly over the node's connections.
     *
     * A CLOSED CLASS IS NEVER ROUTED INTO THE SINK, under either strategy: the
     * Sink would absorb a job the population has to conserve. The reference
     * zeroes the trailing Sink entry in the Probabilities arm and relies on
     * `getRoutingMatrix` to drop it in the Random arm; here it is dropped in
     * both, because the raw routing matrix is also what the model DOCUMENT
     * carries and what the reducibility test reads.
     */
    void set_routing_strategies(qn::Network<T>& model, qn::RoutingMatrix<T>& P, std::size_t from,
                                const std::vector<std::size_t>& outgoing) {
        const qn::NetworkStruct<T>& sn = model.raw_struct();
        const std::size_t K = sn.classes.size();
        const T one = num_traits<T>::from_int(1);

        if (sn.nodes[from - 1].nodetype == lang::NodeType::Source) {
            for (std::size_t r = 1; r <= K; ++r) {
                if (sn.classes[r - 1].type != lang::JobClassType::OPEN) continue;
                for (std::size_t k = 0; k < outgoing.size(); ++k) P.set(r, r, from, outgoing[k], one);
            }
            return;
        }

        for (std::size_t r = 1; r <= K; ++r) {
            const bool closed = sn.classes[r - 1].type != lang::JobClassType::OPEN;
            std::vector<std::size_t> allowed;
            for (std::size_t k = 0; k < outgoing.size(); ++k) {
                if (closed && sn.nodes[outgoing[k] - 1].nodetype == lang::NodeType::Sink) continue;
                allowed.push_back(outgoing[k]);
            }
            const std::string strat = choose_routing_strat();
            if (strat == "Random") {
                model.set_routing(from, r, lang::RoutingStrategy::RAND);
                if (allowed.empty()) continue;
                const double share = 1.0 / double(allowed.size());
                for (std::size_t k = 0; k < allowed.size(); ++k)
                    P.set(r, r, from, allowed[k], num_traits<T>::from_double(share));
            } else {
                model.set_routing(from, r, lang::RoutingStrategy::PROB);
                if (allowed.empty()) continue;
                const std::vector<double> probs = randfixedsumone(allowed.size(), rng_);
                for (std::size_t k = 0; k < allowed.size(); ++k)
                    P.set(r, r, from, allowed[k], num_traits<T>::from_double(probs[k]));
            }
        }
    }

    // -----------------------------------------------------------------------
    // The rejection test
    // -----------------------------------------------------------------------

    /**
     * Empty when the draw is a valid network, otherwise why it is not.
     *
     * A closed chain whose reference station collects no visits has no
     * normalisable visit vector: the jobs are declared at a station the chain
     * never reaches, which is exactly the stranded-chain draw the reference
     * resamples. Every other structural defect the builder can produce throws
     * out of `get_struct()` and is caught by the caller.
     */
    static std::string why_invalid(const qn::NetworkStruct<T>& sn) {
        for (std::size_t c = 0; c < sn.nchains; ++c) {
            if (c >= sn.inchain.size() || sn.inchain[c].empty()) continue;
            const std::vector<std::size_t>& ic = sn.inchain[c];
            bool closed = true;
            for (std::size_t x = 0; x < ic.size(); ++x)
                if (!std::isfinite(num_traits<T>::to_double(sn.classes[ic[x] - 1].population)))
                    closed = false;
            if (!closed) continue;
            const std::size_t rs = sn.classes[ic[0] - 1].refstat;
            if (rs == 0 || rs > sn.station_to_node.size()) continue;
            const std::size_t node = sn.station_to_node[rs - 1];
            if (c >= sn.nodevisits.size() || node == 0 || node > sn.nodevisits[c].rows()) continue;
            double total = 0.0;
            for (std::size_t x = 0; x < ic.size(); ++x)
                total += num_traits<T>::to_double(sn.nodevisits[c](node - 1, ic[x] - 1));
            if (!(total > lang::GlobalConstants::FineTol))
                return "closed chain " + std::to_string(c + 1) +
                       " never visits its reference station '" + sn.nodes[node - 1].name + "'";
        }
        return std::string();
    }

    // -----------------------------------------------------------------------
    // State
    // -----------------------------------------------------------------------

    static constexpr int MAX_SERVERS = 40;
    static constexpr int HIGH_LO = 31, HIGH_HI = 40;
    static constexpr int MED_LO = 11, MED_HI = 20;
    static constexpr int LOW_LO = 1, LOW_HI = 5;

    rng::JavaRandom rng_;
    std::string model_name_ = "nw";
    std::string sched_strat_ = "randomize";
    std::string routing_strat_ = "randomize";
    std::string distribution_ = "randomize";
    std::string cclass_job_load_ = "randomize";
    bool varying_service_rates_ = false;
    bool multi_server_queues_ = false;
    bool random_cs_nodes_ = false;
    bool multi_chain_cs_ = false;
    std::function<Matrix<double>(std::size_t, rng::JavaRandom&)> topology_fcn_ =
        [](std::size_t n, rng::JavaRandom& r) { return rand_graph(n, r); };

    std::vector<std::size_t> stations_;  // node indices of the queues and delays
    std::vector<std::size_t> classes_;   // class indices, open first
    std::size_t source_ = 0, sink_ = 0;
};

}  // namespace gen
}  // namespace line

#endif  // LINE_GEN_NETWORK_GENERATOR_H
