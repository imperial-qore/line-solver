/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_QN_NETWORK_BUILDER_H
#define LINE_LANG_QN_NETWORK_BUILDER_H

/**
 * The Network constructor API: Queue, Delay, Source, Sink, Router,
 * ClassSwitch, Cache, Fork and Join, the job classes, and `link`.
 *
 * This is the port of what a MATLAB model script writes -- `Network`,
 * `Queue(model, name, sched)`, `queue.setService(class, dist)`,
 * `model.link(P)` -- feeding the SAME refresh as every other front end, so a
 * model built here and a model read from a file produce the same
 * `NetworkStruct`. The split is deliberate and mirrors `lqn_builder.h` /
 * `lqn_reader.h`: stage one constructs, stage two (`NetworkStruct::refresh_*`)
 * derives, and only stage two is allowed to compute anything.
 *
 * INDEX SPACES. Every method takes and returns 1-based NODE indices, which is
 * what a model script deals in; the station index is an internal detail of the
 * struct and is looked up here. `add_*` returns the node index of the node it
 * created, and `add_*_class` the class index.
 *
 * WHAT IS REFUSED, and where. The builder accepts the feature set SolverMVA
 * declares (`SolverMVA.getFeatureSet`); the refresh refuses a state-dependent
 * routing strategy by name, and each solver refuses what it cannot analyse.
 * Nothing is silently mapped onto a neighbour.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "line/api/sn/sn_fj_nodevisits_mmt.h"
#include "line/lang/distribution.h"
#include "line/lang/prior.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace qn {

/**
 * The routing matrix a model script fills in, MATLAB's `P` cell array.
 *
 * `P{r,s}(i,j)` is the probability that a job leaving node i in class r enters
 * node j in class s. The single-class form `set(i, j, p)` is the shorthand a
 * one-class model uses, and is exactly `set(1, 1, i, j, p)`.
 */
template <class T>
class RoutingMatrix {
  public:
    void set(std::size_t r, std::size_t s, std::size_t i, std::size_t j, const T& p) {
        entries[Key(r, s, i, j)] = p;
    }
    void set(std::size_t i, std::size_t j, const T& p) { set(1, 1, i, j, p); }

    T get(std::size_t r, std::size_t s, std::size_t i, std::size_t j) const {
        auto it = entries.find(Key(r, s, i, j));
        return it == entries.end() ? num_traits<T>::from_int(0) : it->second;
    }

    struct Key {
        std::size_t r, s, i, j;
        Key(std::size_t r_, std::size_t s_, std::size_t i_, std::size_t j_)
            : r(r_), s(s_), i(i_), j(j_) {}
        bool operator<(const Key& o) const {
            if (r != o.r) return r < o.r;
            if (s != o.s) return s < o.s;
            if (i != o.i) return i < o.i;
            return j < o.j;
        }
    };
    std::map<Key, T> entries;
};

/**
 * A queueing network under construction.
 *
 * The object owns a `NetworkStruct` and hands out a refreshed reference to it.
 * `get_struct()` runs the whole refresh chain each time it is called, because
 * every setter can have moved a quantity the chain derives; a caller that
 * changes nothing and asks twice pays for the second refresh, which is the same
 * bargain MATLAB's `hasStruct` cache makes on the other side.
 */
template <class T>
class Network {
  public:
    explicit Network(const std::string& nm) { sn_.name = nm; }

    // -----------------------------------------------------------------------
    // Nodes
    // -----------------------------------------------------------------------

    /**
     * A queueing station. The default discipline is FCFS, as in MATLAB.
     *
     * INF SCHEDULING BUILDS A DELAY. MATLAB's Queue constructor sets
     * numberOfServers = Inf on that branch and getNodeTypes then reports the
     * station as NodeType.Delay (the JAR does both with Integer.MAX_VALUE), and
     * the analyzers partition the stations on those two fields:
     * sn_get_product_form_chain_params splits by NODETYPE, so a Queue node left
     * at one server was handed to the linearizer as a finite-server queue and
     * reported utilization 1 and a saturated response time where the reference
     * reports a delay (2.4494/19.5506 against 0.1253/21.8747 on a two-station
     * closed model with 22 jobs).
     */
    std::size_t add_queue(const std::string& nm, SchedStrategy sched = SchedStrategy::FCFS) {
        Station<T> st;
        st.name = nm;
        st.nodetype = (sched == SchedStrategy::INF) ? NodeType::Delay : NodeType::Queue;
        st.sched = sched;
        st.nservers = (sched == SchedStrategy::INF)
                          ? std::numeric_limits<double>::infinity()
                          : 1.0;
        const std::size_t ist = sn_.add_station(st);
        init_node(sn_.station_to_node[ist - 1]);
        return sn_.station_to_node[ist - 1];
    }

    /** An infinite-server station (a Delay, MATLAB's `Delay` / `DelayStation`). */
    std::size_t add_delay(const std::string& nm) {
        Station<T> st;
        st.name = nm;
        st.nodetype = NodeType::Delay;
        st.sched = SchedStrategy::INF;
        st.nservers = std::numeric_limits<double>::infinity();
        const std::size_t ist = sn_.add_station(st);
        init_node(sn_.station_to_node[ist - 1]);
        return sn_.station_to_node[ist - 1];
    }

    /**
     * The external arrival station.
     *
     * It IS a station -- its "service" process is the arrival process -- with
     * the EXT discipline and one server, exactly as MATLAB's Source is built.
     */
    std::size_t add_source(const std::string& nm) {
        if (sn_.sourceIdx != 0) throw InputError("Network: the model already has a Source");
        Station<T> st;
        st.name = nm;
        st.nodetype = NodeType::Source;
        st.sched = SchedStrategy::EXT;
        st.nservers = 1.0;
        const std::size_t ist = sn_.add_station(st);
        sn_.sourceIdx = ist;
        init_node(sn_.station_to_node[ist - 1]);
        return sn_.station_to_node[ist - 1];
    }

    /** The external departure node. It is NOT a station and holds no jobs. */
    std::size_t add_sink(const std::string& nm) {
        if (sn_.sinkNode != 0) throw InputError("Network: the model already has a Sink");
        const std::size_t nd = sn_.add_node(nm, NodeType::Sink, false);
        sn_.sinkNode = nd;
        init_node(nd);
        return nd;
    }

    /** A stateless routing node. */
    std::size_t add_router(const std::string& nm) {
        const std::size_t nd = sn_.add_node(nm, NodeType::Router, false);
        init_node(nd);
        return nd;
    }

    /**
     * A Logger node: a pass-through that records every job crossing it.
     *
     * It holds no jobs and changes no routing probability, so it is eliminated
     * by the same stochastic complement that removes a Router; what makes it a
     * distinct node type is that `used_lang_features` emits Logger/LogTunnel
     * from it, so a solver with no logging refuses the model by name rather
     * than silently dropping the trace the user asked for.
     */
    std::size_t add_logger(const std::string& nm, const std::string& log_file = std::string()) {
        const std::size_t nd = sn_.add_node(nm, NodeType::Logger, false);
        init_node(nd);
        // `Logger.m` splits the argument and keeps only the base name; the
        // directory is the model's log path, which every Logger shares.
        std::string base = log_file;
        const std::size_t slash = base.find_last_of('/');
        if (slash != std::string::npos) base = base.substr(slash + 1);
        sn_.nodes[nd - 1].logger.file_name = base;
        return nd;
    }

    /** `Network.setLogPath`: the directory every Logger of this model writes into. */
    void set_log_path(const std::string& path) { sn_.log_path = path; }

    /**
     * A ClassSwitch node carrying the (nclasses x nclasses) switching matrix.
     *
     * The classes must exist before the node, since the matrix is indexed by
     * them; that is also the order a MATLAB script writes.
     */
    std::size_t add_class_switch(const std::string& nm, const Matrix<T>& C) {
        if (C.rows() != sn_.classes.size() || C.cols() != sn_.classes.size())
            throw InputError("ClassSwitch '" + nm +
                             "': the matrix must be (nclasses x nclasses); declare the classes "
                             "before the node");
        const std::size_t nd = sn_.add_node(nm, NodeType::ClassSwitch, false);
        sn_.csmatrix[nd] = C;
        init_node(nd);
        return nd;
    }

    /**
     * A ClassSwitch node whose matrix is installed LATER, by
     * `set_class_switch_matrix`.
     *
     * FOR A FILE READER, which meets the nodes before the classes. A `.jsimg`
     * and a `.lqnx` both list their nodes first, and a closed class names its
     * reference STATION, so neither order can be satisfied without deferring
     * one of the two -- MATLAB's own `ClassSwitch` constructor stores the matrix
     * without checking it against a class list that does not exist yet, which is
     * the same deferral by another name.
     *
     * THE MATRIX IS LEFT EMPTY, not filled with an identity: an identity is a
     * valid switching matrix (every class keeps its own), so a caller who forgot
     * to install the real one would get a plausible model instead of an error.
     */
    std::size_t add_class_switch(const std::string& nm) {
        const std::size_t nd = sn_.add_node(nm, NodeType::ClassSwitch, false);
        sn_.csmatrix[nd] = Matrix<T>();
        init_node(nd);
        return nd;
    }

    /** Install the switching matrix of a ClassSwitch created without one. */
    void set_class_switch_matrix(std::size_t node, const Matrix<T>& C) {
        if (node == 0 || node > sn_.nodes.size())
            throw InputError("set_class_switch_matrix: node index is out of range");
        if (sn_.nodes[node - 1].nodetype != NodeType::ClassSwitch)
            throw InputError("set_class_switch_matrix: node '" + sn_.nodes[node - 1].name +
                             "' is not a ClassSwitch");
        if (C.rows() != sn_.classes.size() || C.cols() != sn_.classes.size())
            throw InputError("set_class_switch_matrix: the matrix of '" +
                             sn_.nodes[node - 1].name +
                             "' must be (nclasses x nclasses)");
        sn_.csmatrix[node] = C;
    }

    /** A Fork node. It holds no jobs and is removed by the stochastic complement.
     *  `tasks_per_link` (Fork.output.tasksPerLink) defaults to 1. */
    std::size_t add_fork(const std::string& nm, double tasks_per_link = 1.0) {
        const std::size_t nd = sn_.add_node(nm, NodeType::Fork, false);
        sn_.nodes[nd - 1].tasks_per_link = tasks_per_link;
        init_node(nd);
        return nd;
    }

    /**
     * Variable forking levels on an existing Fork, the twin of MATLAB
     * `Fork.setTasksPerLink(jobclass, n)`,
     * `Fork.setTasksPerLinkDistribution(jobclass, dist [, destNode])` and
     * `Fork.setBranchProbability(jobclass, destNode, p)`.
     *
     * The block is allocated lazily and only on a fork that actually declares an
     * override, so a plain fork has no `sn.forkparam` entry at all and every
     * consumer can tell the classic case from the variable one by asking
     * `sn.fork_param_of(f)` for a null.
     *
     * `dest_node` 0 means every outgoing link of that class. Classes and nodes
     * are 1-based, as everywhere in this port.
     *
     * An override is RECORDED and replayed when the routing is installed, so
     * these may be called in either order with respect to `link()` -- as the
     * MATLAB, Java and Python setters may, which store the override on the
     * Forker section and materialise it at refresh time. Recording rather than
     * writing is what buys that: `dest_node = 0` means "every link this class
     * takes", a set that does not exist until the routing does.
     */
    void set_fork_tasks_per_link(std::size_t fork_node, std::size_t jobclass,
                                 double tasks, std::size_t dest_node = 0) {
        ForkOverride ov;
        ov.kind = ForkOverride::TASKS;
        ov.fork = fork_node;
        ov.cls = jobclass;
        ov.dest = dest_node;
        ov.value = tasks;
        record_fork_override(ov);
    }

    /** A random jobs-per-link degree, redrawn per link and per forked job. */
    void set_fork_tasks_per_link_dist(std::size_t fork_node, std::size_t jobclass,
                                      const lang::Distrib<T>& dist,
                                      std::size_t dest_node = 0) {
        if (dist.type != lang::ProcessType::DISCRETESAMPLER)
            throw InputError("set_fork_tasks_per_link_dist: the jobs-per-link "
                             "distribution must be a DiscreteSampler");
        ForkOverride ov;
        ov.kind = ForkOverride::DIST;
        ov.fork = fork_node;
        ov.cls = jobclass;
        ov.dest = dest_node;
        ov.dist = dist;
        record_fork_override(ov);
    }

    /** A branch that fires only with probability `prob`. */
    void set_fork_branch_probability(std::size_t fork_node, std::size_t jobclass,
                                     std::size_t dest_node, double prob) {
        if (prob < 0.0 || prob > 1.0)
            throw InputError("set_fork_branch_probability: a branch activation "
                             "probability must lie in [0,1]");
        ForkOverride ov;
        ov.kind = ForkOverride::PROB;
        ov.fork = fork_node;
        ov.cls = jobclass;
        ov.dest = dest_node;
        ov.value = prob;
        record_fork_override(ov);
    }

    /**
     * A Join node, which IS a station: it serves at an infinite rate, and the
     * synchronisation delay is supplied by the fork-join transform.
     */
    std::size_t add_join(const std::string& nm, std::size_t fork_node) {
        const std::size_t nd = add_join_unbound(nm);
        bind_join(nd, fork_node);
        return nd;
    }

    /**
     * The Join station on its own, with the fork left to `bind_join`.
     *
     * A Join IS a station, so creating it late shifts every station index after
     * it, and the result document is indexed by station row. A reader that has
     * to see the Fork first therefore declares the Join here, in the position
     * the model gives it, and binds the pair once the Fork exists.
     */
    std::size_t add_join_unbound(const std::string& nm) {
        Station<T> st;
        st.name = nm;
        st.nodetype = NodeType::Join;
        st.sched = SchedStrategy::INF;
        st.nservers = std::numeric_limits<double>::infinity();
        const std::size_t ist = sn_.add_station(st);
        const std::size_t nd = sn_.station_to_node[ist - 1];
        init_node(nd);
        return nd;
    }

    /** Record which Fork a Join created by `add_join_unbound` closes. */
    void bind_join(std::size_t join_node, std::size_t fork_node) {
        if (join_node == 0 || join_node > sn_.nodes.size() ||
            sn_.nodes[join_node - 1].nodetype != NodeType::Join)
            throw InputError("bind_join: the node being bound is not a Join");
        if (fork_node == 0 || fork_node > sn_.nodes.size() ||
            sn_.nodes[fork_node - 1].nodetype != NodeType::Fork)
            throw InputError("Join '" + sn_.nodes[join_node - 1].name +
                             "': the node it closes is not a Fork");
        sn_.fj.emplace_back(fork_node, join_node);
    }

    /**
     * A Place: an SPN token container. Modelled as an INF-scheduled station so
     * the marginal machinery treats its tokens as "in service", which is the
     * encoding `State.toMarginal` folds the buffer slot back into.
     */
    std::size_t add_place(const std::string& nm) {
        Station<T> st;
        st.name = nm;
        st.nodetype = NodeType::Place;
        st.sched = SchedStrategy::INF;
        st.nservers = std::numeric_limits<double>::infinity();
        const std::size_t ist = sn_.add_station(st);
        init_node(sn_.station_to_node[ist - 1]);
        return sn_.station_to_node[ist - 1];
    }

    /**
     * A Transition: the firing rules of an SPN, as `Transition` in MATLAB.
     *
     * A transition has MODES, not classes: each mode has its own enabling and
     * inhibiting conditions over the places, its own firing effect, and its own
     * firing process. The parameters are transcribed from the reference's
     * `refreshPetriNetNodes.m`, so `State.fromMarginal` finds the fields it
     * reads instead of a node it cannot decode.
     */
    std::size_t add_transition(const std::string& nm, const TransitionParam<T>& par) {
        if (par.nmodes == 0)
            throw InputError("add_transition: a transition needs at least one mode");
        if (par.enabling.size() != par.nmodes || par.firing.size() != par.nmodes)
            throw InputError(
                "add_transition: enabling and firing must have one entry per mode");
        const std::size_t nd = sn_.add_node(nm, NodeType::Transition, true);
        sn_.transparam[nd] = par;
        init_node(nd);
        return nd;
    }

    /**
     * `Queue.setRetrial(...)`: a station with an ORBIT instead of a waiting
     * line. An arrival finding every server busy joins the orbit and re-attempts
     * at `rate`; a completion does NOT promote from the orbit, so the state is
     * an (in-service, orbit) split rather than an ordered buffer.
     */
    void set_retrial(std::size_t node, std::size_t cls, const Distrib<T>& proc, const T& rate,
                     int max_attempts = 0) {
        if (node == 0 || node > sn_.nodes.size())
            throw InputError("set_retrial: node index is out of range");
        const std::size_t ist = sn_.nodes[node - 1].station;
        if (ist == 0) throw InputError("set_retrial: node is not a station");
        // Size against the LIVE class list: sn_.nclasses is a finalize-time
        // field and is still zero while the model is being built, so sizing
        // against it left the vectors empty and the write below went out of
        // bounds.
        const std::size_t K = sn_.classes.size();
        if (cls == 0 || cls > K)
            throw InputError("set_retrial: class index is out of range");
        RetrialParam<T>& rp = sn_.retrialparam[ist];
        if (rp.retrial_proc.size() < K) {
            rp.retrial_proc.resize(K);
            rp.retrial_rate.resize(K, num_traits<T>::from_int(0));
            rp.max_attempts.resize(K, 0);
        }
        rp.retrial_proc[cls - 1] = proc;
        rp.retrial_rate[cls - 1] = rate;
        rp.max_attempts[cls - 1] = max_attempts;
    }

    /**
     * `Queue.setPatience(class, dist, type)`: the abandonment timer of a job
     * WAITING at the station, and which impatience rule the timer belongs to.
     */
    void set_patience(std::size_t node, std::size_t cls, const Distrib<T>& dist,
                      lang::ImpatienceType kind = lang::ImpatienceType::RENEGING) {
        Station<T>& st = station_ref(node, cls, "set_patience");
        grow_class_slot(st.patience, cls, Distrib<T>::disabled_dist());
        grow_class_slot(st.impatience, cls, lang::ImpatienceType::NONE);
        st.patience[cls - 1] = dist;
        st.impatience[cls - 1] = kind;
    }

    /** `Queue.setOrbitImpatience(class, dist)`: abandonment from the retrial orbit. */
    void set_orbit_impatience(std::size_t node, std::size_t cls, const Distrib<T>& dist) {
        Station<T>& st = station_ref(node, cls, "set_orbit_impatience");
        grow_class_slot(st.orbit_impatience, cls, Distrib<T>::disabled_dist());
        st.orbit_impatience[cls - 1] = dist;
    }

    /** `Queue.setBatchRejectProbability(class, p)`. */
    void set_batch_reject(std::size_t node, std::size_t cls, const T& p) {
        Station<T>& st = station_ref(node, cls, "set_batch_reject");
        grow_class_slot(st.batch_reject, cls, num_traits<T>::from_int(0));
        st.batch_reject[cls - 1] = p;
    }

    /**
     * `Queue.setBalking(class, strategy, thresholds)`: an arrival that refuses
     * to JOIN, on the state it finds. Distinct from reneging, which abandons a
     * job that has already joined.
     */
    void set_balking(std::size_t node, std::size_t cls, lang::BalkingStrategy strategy,
                     const std::vector<typename Station<T>::BalkingThreshold>& thresholds) {
        Station<T>& st = station_ref(node, cls, "set_balking");
        grow_class_slot(st.balking, cls, typename Station<T>::BalkingParam());
        st.balking[cls - 1].strategy = strategy;
        st.balking[cls - 1].thresholds = thresholds;
    }

    /** `Queue.addServerType(...)`: one heterogeneous server pool of the station. */
    void add_server_type(std::size_t node, const typename Station<T>::ServerType& stype) {
        const std::size_t ist = station_of(node, "add_server_type");
        sn_.stations[ist - 1].server_types.push_back(stype);
        // The pools size the station, as in MATLAB/JAR updateTotalServerCount
        double total = 0.0;
        for (const auto& pool : sn_.stations[ist - 1].server_types) total += pool.count;
        sn_.stations[ist - 1].nservers = total;
    }

    /**
     * `Queue.setServerParallelism(class, n)`: the servers a job seizes for the
     * whole of its service. The station then serves at most floor(c/n) such jobs
     * at a time.
     */
    void set_server_parallelism(std::size_t node, std::size_t cls, std::size_t n) {
        Station<T>& st = station_ref(node, cls, "set_server_parallelism");
        if (n < 1) {
            throw InputError("set_server_parallelism: parallelism must be a positive integer");
        }
        if (std::isfinite(st.nservers) && static_cast<double>(n) > st.nservers) {
            throw InputError(
                "set_server_parallelism: parallelism " + std::to_string(n) + " exceeds the " +
                std::to_string(static_cast<long long>(st.nservers)) + " servers of station '" +
                sn_.nodes[node - 1].name + "', so a job of this class could never enter service");
        }
        grow_class_slot(st.server_parallelism, cls, static_cast<std::size_t>(1));
        st.server_parallelism[cls - 1] = n;
    }

    /** `Queue.setHeteroSchedPolicy(...)`: how the server pools are picked among. */
    void set_hetero_sched_policy(std::size_t node, lang::HeteroSchedPolicy policy) {
        sn_.stations[station_of(node, "set_hetero_sched_policy") - 1].hetero_policy = policy;
    }

    /**
     * `Source.setArrivalBatch(class, dist)`: the batch-size law released at each
     * arrival epoch. The arrival process itself only spaces the epochs.
     */
    void set_arrival_batch(std::size_t node, std::size_t cls, const Distrib<T>& dist) {
        Station<T>& st = station_ref(node, cls, "set_arrival_batch");
        grow_class_slot(st.arrival_batch, cls, Distrib<T>::disabled_dist());
        st.arrival_batch[cls - 1] = dist;
    }

    /** `Source.markedClasses`: the 1-based class carried by each mark of an MMAP. */
    void set_marked_classes(std::size_t node, const std::vector<std::size_t>& classes) {
        sn_.stations[station_of(node, "set_marked_classes") - 1].marked_classes = classes;
    }

    /** `Place.setDepartureDiscipline(class, rule)`. */
    void set_departure_discipline(std::size_t node, std::size_t cls,
                                  lang::DepartureDiscipline rule) {
        Station<T>& st = station_ref(node, cls, "set_departure_discipline");
        grow_class_slot(st.departure_discipline, cls, lang::DepartureDiscipline::NORMAL);
        st.departure_discipline[cls - 1] = rule;
    }

    /** `Place.setState(marking)`: the initial token count of the place, per class. */
    void set_initial_marking(std::size_t node, const std::vector<T>& tokens) {
        if (node == 0 || node > sn_.nodes.size())
            throw InputError("set_initial_marking: node index is out of range");
        if (sn_.nodes[node - 1].nodetype != NodeType::Place)
            throw InputError("set_initial_marking: only a Place carries an initial marking");
        sn_.initmarking[node] = tokens;
    }

    /**
     * `StatefulNode.setStatePrior(space, prior)`: a distribution over the rows
     * of a DECLARED state space. The two are set together because a prior
     * indexes that space and means nothing without it.
     */
    void set_state_prior(std::size_t node, const Matrix<T>& space, const std::vector<T>& prior) {
        if (node == 0 || node > sn_.nodes.size())
            throw InputError("set_state_prior: node index is out of range");
        if (space.rows() != prior.size())
            throw InputError(
                "set_state_prior: the prior has one entry per ROW of the declared state space");
        sn_.statespace[node] = space;
        sn_.stateprior[node] = prior;
    }

    /** `Join.setStrategy(...)`: STD waits for every sibling, PARTIAL for a quorum. */
    void set_join_strategy(std::size_t node, lang::JoinStrategy strategy, double quorum = 0.0) {
        if (node == 0 || node > sn_.nodes.size())
            throw InputError("set_join_strategy: node index is out of range");
        if (sn_.nodes[node - 1].nodetype != NodeType::Join)
            throw InputError("set_join_strategy: the node is not a Join");
        typename NetworkStruct<T>::JoinDecl jd;
        jd.strategy = strategy;
        jd.quorum = quorum;
        sn_.joindecl[node] = jd;
    }

    /** The per-destination weights of a WRROBIN dispatcher, per (node, class). */
    void set_routing_weights(std::size_t node, std::size_t cls,
                             const std::map<std::size_t, double>& weights) {
        if (node == 0 || node > sn_.nodes.size())
            throw InputError("set_routing_weights: node index is out of range");
        std::vector<std::map<std::size_t, double>>& rw = sn_.nodes[node - 1].routing_weights;
        if (rw.size() < sn_.classes.size()) rw.resize(sn_.classes.size());
        if (cls == 0 || cls > rw.size())
            throw InputError("set_routing_weights: class index is out of range");
        rw[cls - 1] = weights;
    }

    /** The d of a power-of-d (SQ) dispatcher, per (node, class). */
    void set_routing_param(std::size_t node, std::size_t cls, int d) {
        if (node == 0 || node > sn_.nodes.size())
            throw InputError("set_routing_param: node index is out of range");
        std::vector<int>& rp = sn_.nodes[node - 1].routing_param;
        if (rp.size() < sn_.classes.size()) rp.resize(sn_.classes.size(), 0);
        if (cls == 0 || cls > rp.size())
            throw InputError("set_routing_param: class index is out of range");
        rp[cls - 1] = d;
    }

    /**
     * `Queue.setDelayOff(class, setupTime, delayoffTime)`: the station powers
     * down after sitting idle for the delay-off time, and the next arrival pays
     * the setup time before service. The pair is set together because a setup
     * with no delay-off never fires and a delay-off with no setup is free.
     */
    void set_setup_delayoff(std::size_t node, std::size_t cls, const Distrib<T>& setup,
                            const Distrib<T>& delayoff) {
        const std::size_t ist = station_of(node, "set_setup_delayoff");
        const std::size_t K = sn_.classes.size();
        if (cls == 0 || cls > K)
            throw InputError("set_setup_delayoff: class index is out of range");
        SetupDelayOffParam<T>& sp = sn_.setupparam[ist];
        if (sp.setup.size() < K) {
            sp.setup.resize(K, Distrib<T>::disabled_dist());
            sp.delayoff.resize(K, Distrib<T>::disabled_dist());
        }
        sp.setup[cls - 1] = setup;
        sp.delayoff[cls - 1] = delayoff;
    }

    /**
     * `Queue.setBreakdown(failure, repair, downService)`: the server alternates
     * up and down on the two clocks.
     *
     * BOTH CLOCKS ARE REQUIRED. A server that fails and is never repaired is a
     * different model -- an absorbing one -- and the reference declines to infer
     * it from a missing repair time rather than treating it as infinite.
     *
     * `down_service` is per class and OPTIONAL; an entry left disabled means the
     * class gets no service while the server is down, which is the ordinary
     * reading. A declared one must be EXPONENTIAL: a phase-type degraded service
     * would need its own phase block in the joint chain and no codebase builds
     * one, so it is refused by name rather than approximated by its mean.
     */
    void set_breakdown(std::size_t node, const Distrib<T>& failure, const Distrib<T>& repair,
                       const std::vector<Distrib<T>>& down_service = std::vector<Distrib<T>>()) {
        const std::size_t ist = station_of(node, "set_breakdown");
        if (failure.disabled || repair.disabled)
            throw InputError(
                "set_breakdown: both a failure and a repair time are required; a server that "
                "never recovers is an absorbing model, not a breakdown");
        const double fm = num_traits<T>::to_double(failure.mean);
        const double rm = num_traits<T>::to_double(repair.mean);
        if (!(fm > 0.0) || !(rm > 0.0))
            throw InputError("set_breakdown: the failure and repair times must have positive means");
        const std::size_t K = sn_.classes.size();
        BreakdownParam<T>& bp = sn_.breakdownparam[ist];
        bp.failure = failure;
        bp.repair = repair;
        bp.failure_rate = T(num_traits<T>::from_int(1) / failure.mean);
        bp.repair_rate = T(num_traits<T>::from_int(1) / repair.mean);
        bp.down_service_rates.assign(K, num_traits<T>::from_int(0));
        for (std::size_t r = 0; r < K && r < down_service.size(); ++r) {
            const Distrib<T>& d = down_service[r];
            if (d.disabled) continue;
            if (d.type != lang::ProcessType::EXP)
                throw UnsupportedError(
                    "set_breakdown: station '" + sn_.stations[ist - 1].name +
                    "': the down-server service distribution must be exponential; a phase-type "
                    "degraded service would need its own phase block in the joint chain");
            if (num_traits<T>::to_double(d.mean) > 0.0)
                bp.down_service_rates[r] = T(num_traits<T>::from_int(1) / d.mean);
        }
    }

    /** A Cache node with its item population, list capacities and popularity. */
    std::size_t add_cache(const std::string& nm, const CacheParam<T>& par) {
        const std::size_t nd = sn_.add_node(nm, NodeType::Cache, true);
        sn_.nodeparam[nd] = par;
        init_node(nd);
        return nd;
    }

    /**
     * `Cache.setItemReadClasses(readClasses, hitClasses)`: declare that
     * `read_classes[i]` is the request stream for item i at this cache. Use at the
     * cache the exogenous requests enter, where the per-item classes are the
     * caller's own; item popularity is then carried by the per-class request rates
     * rather than by a popularity the cache draws from. This is what keeps a cache
     * network free of arc-level class switching, so no class acquires a default
     * route into the cache the model never intended. `hit_classes` is either one
     * class shared by every item or one per item.
     */
    void set_item_read_classes(std::size_t cache_node,
                               const std::vector<std::size_t>& read_classes,
                               const std::vector<std::size_t>& hit_classes) {
        auto it = sn_.nodeparam.find(cache_node);
        if (it == sn_.nodeparam.end())
            throw InputError("setItemReadClasses: node is not a Cache");
        CacheParam<T>& cp = it->second;
        const std::size_t nitems = cp.nitems;
        if (read_classes.size() != nitems)
            throw InputError("setItemReadClasses: pass exactly one read class per item");
        const std::vector<std::size_t> hit =
            per_item_classes(hit_classes, nitems, "setItemReadClasses hit");
        const std::size_t K = sn_.classes.size();
        cp.pread.resize(K);
        cp.preadkind.resize(K);
        cp.hitclass.resize(K, 0);
        cp.missclass.resize(K, 0);
        cp.classitem.resize(K, 0);
        for (std::size_t i = 0; i < nitems; ++i) {
            std::vector<T> onehot(nitems, num_traits<T>::from_int(0));
            onehot[i] = num_traits<T>::from_int(1);
            cp.pread[read_classes[i] - 1] = onehot;
            cp.hitclass[read_classes[i] - 1] = hit[i];
            cp.classitem[read_classes[i] - 1] = i + 1;
        }
        cache_item_classes_[cache_node] = read_classes;
    }

    /**
     * `Cache.setMissCache(readClass, nextCache, hitClassAtNext)`: send this cache's
     * misses to `next_cache` preserving item identity, by minting one class per item
     * there and making the miss class of item i here its read class for item i.
     * Returns the minted classes. The cache-to-cache arc itself is registered here
     * and injected by `link()`, so the caller routes only its own topology.
     */
    std::vector<std::size_t> set_miss_cache(std::size_t cache_node, std::size_t next_cache,
                                            const std::vector<std::size_t>& hit_classes_at_next) {
        auto it = sn_.nodeparam.find(cache_node);
        auto itn = sn_.nodeparam.find(next_cache);
        if (it == sn_.nodeparam.end() || itn == sn_.nodeparam.end())
            throw InputError("setMissCache: both nodes must be Caches");
        if (it->second.nitems != itn->second.nitems)
            throw InputError("setMissCache: a cache network requires one common item set");
        auto self_it = cache_item_classes_.find(cache_node);
        if (self_it == cache_item_classes_.end())
            throw InputError("setMissCache: call setItemReadClasses on the source cache first");
        const std::vector<std::size_t> self_cls = self_it->second;
        const std::size_t nitems = it->second.nitems;
        const std::vector<std::size_t> hit =
            per_item_classes(hit_classes_at_next, nitems, "setMissCache hit");

        std::vector<std::size_t> minted;
        minted.reserve(nitems);
        const bool closed = sn_.classes[self_cls[0] - 1].type == JobClassType::CLOSED;
        const std::size_t refstat = sn_.classes[self_cls[0] - 1].refstat;
        for (std::size_t i = 0; i < nitems; ++i) {
            std::size_t cls;
            if (closed)
                cls = add_closed_class(sn_.nodes[next_cache - 1].name + "_item" + std::to_string(i + 1),
                                       0.0, sn_.station_to_node[refstat - 1]);
            else
                cls = add_open_class(sn_.nodes[next_cache - 1].name + "_item" + std::to_string(i + 1));
            minted.push_back(cls);
        }
        // the class count grew, so re-size both caches' per-class vectors once
        const std::size_t K = sn_.classes.size();
        for (std::size_t nd : {cache_node, next_cache}) {
            CacheParam<T>& c = sn_.nodeparam.at(nd);
            c.pread.resize(K);
            c.preadkind.resize(K);
            c.hitclass.resize(K, 0);
            c.missclass.resize(K, 0);
            c.classitem.resize(K, 0);
        }
        CacheParam<T>& src = sn_.nodeparam.at(cache_node);
        CacheParam<T>& dst = sn_.nodeparam.at(next_cache);
        for (std::size_t i = 0; i < nitems; ++i) {
            std::vector<T> onehot(nitems, num_traits<T>::from_int(0));
            onehot[i] = num_traits<T>::from_int(1);
            dst.pread[minted[i] - 1] = onehot;
            dst.hitclass[minted[i] - 1] = hit[i];
            dst.classitem[minted[i] - 1] = i + 1;
            src.missclass[self_cls[i] - 1] = minted[i];
            // the miss hop is part of the construction, not of the user's topology:
            // register it here and let link() inject it, as the MATLAB helper does
            // through retrievalRoutingEntries
            cache_miss_arcs_.push_back(CacheMissArc(minted[i], cache_node, next_cache));
        }
        cache_item_classes_[next_cache] = minted;
        return minted;
    }

    /**
     * `Cache.setItemMissClass(readClass, missClasses)`: terminate a cache network,
     * every per-item class of this cache reporting a miss as the matching entry of
     * `miss_classes`, which the caller routes onward.
     */
    void set_item_miss_class(std::size_t cache_node, const std::vector<std::size_t>& miss_classes) {
        auto it = sn_.nodeparam.find(cache_node);
        if (it == sn_.nodeparam.end())
            throw InputError("setItemMissClass: node is not a Cache");
        auto self_it = cache_item_classes_.find(cache_node);
        if (self_it == cache_item_classes_.end())
            throw InputError("setItemMissClass: the cache has no per-item classes");
        const std::vector<std::size_t>& self_cls = self_it->second;
        const std::vector<std::size_t> miss =
            per_item_classes(miss_classes, self_cls.size(), "setItemMissClass miss");
        CacheParam<T>& cp = it->second;
        cp.missclass.resize(sn_.classes.size(), 0);
        for (std::size_t i = 0; i < self_cls.size(); ++i)
            cp.missclass[self_cls[i] - 1] = miss[i];
    }

    /**
     * `Cache.setRetrievalSystem(readClass, missClass, queues)`: a delayed-hit
     * cache whose misses are fetched by circulating a per-item retrieval class
     * through `queue_nodes` and back to the cache. Creates one retrieval class
     * per item, each inheriting the read class's service at every retrieval queue
     * (call `set_service(queue, readClass, ...)` first); the routing among the
     * cache and the queues is inherited from the read class in `link()`. Records
     * the retrieval capacity (nitems - total cache capacity), the queue node
     * list, and the item -> retrieval-class map that `cache_retrieval_inputs`
     * reads. Must be called after the read/miss classes and the cache exist.
     */
    void set_retrieval_system(std::size_t cache_node, std::size_t read_class,
                              std::size_t miss_class, const std::vector<std::size_t>& queue_nodes) {
        auto it = sn_.nodeparam.find(cache_node);
        if (it == sn_.nodeparam.end())
            throw InputError("setRetrievalSystem: node is not a Cache");
        CacheParam<T>& cp = it->second;
        if (queue_nodes.empty())
            throw InputError("setRetrievalSystem: the retrieval system has no stations");
        const std::size_t nitems = cp.nitems;
        int totalcap = 0;
        for (int c : cp.itemcap) totalcap += c;
        cp.retrieval_capacity = static_cast<int>(nitems) - totalcap;
        cp.retrieval_queues[read_class - 1] = queue_nodes;  // 0-based key, 1-based nodes

        // capture the read class's service distribution at each queue up front
        std::vector<Distrib<T> > svc;
        svc.reserve(queue_nodes.size());
        for (std::size_t q : queue_nodes) {
            const std::size_t st = station_of(q, "setRetrievalSystem");
            svc.push_back(sn_.service[st - 1][read_class - 1]);
        }

        cp.retrieval_classes.assign(nitems, std::vector<std::size_t>(sn_.classes.size(), 0));
        const bool closed = sn_.classes[read_class - 1].type == JobClassType::CLOSED;
        const std::size_t refstat = sn_.classes[read_class - 1].refstat;
        for (std::size_t i = 0; i < nitems; ++i) {
            std::size_t rc;
            if (closed)
                rc = add_closed_class(sn_.classes[read_class - 1].name + "_retrievalClass_" +
                                          std::to_string(i + 1),
                                      0.0, sn_.station_to_node[refstat - 1]);
            else
                rc = add_open_class(sn_.classes[read_class - 1].name + "_retrievalClass_" +
                                    std::to_string(i + 1));
            for (std::size_t s = 0; s < queue_nodes.size(); ++s) {
                set_service(queue_nodes[s], rc, svc[s]);
                // at most one retrieval of a given item is ever in flight
                set_class_capacity(queue_nodes[s], rc, 1.0);
            }
            // grow the retrieval_classes rows to the new class count and record
            for (std::size_t k = 0; k < nitems; ++k)
                cp.retrieval_classes[k].resize(sn_.classes.size(), 0);
            cp.retrieval_classes[i][read_class - 1] = rc;
            // The returning READ of a retrieval class reads ITS OWN item, and is
            // logged as the miss that started the fetch. Without these two the
            // state-based solvers see a class that reads nothing, so the fetch
            // never completes and the retrieval sub-network is never entered.
            const std::size_t K = sn_.classes.size();
            cp.pread.resize(K);
            cp.preadkind.resize(K);
            cp.hitclass.resize(K, 0);
            cp.missclass.resize(K, 0);
            std::vector<T> onehot(nitems, num_traits<T>::from_int(0));
            onehot[i] = num_traits<T>::from_int(1);
            cp.pread[rc - 1] = onehot;
            cp.missclass[rc - 1] = miss_class;
        }
    }

    // -----------------------------------------------------------------------
    // Classes
    // -----------------------------------------------------------------------

    /** A closed class of the given population, referencing a station node. */
    std::size_t add_closed_class(const std::string& nm, double njobs, std::size_t refstat_node,
                                 int prio = 0) {
        if (!(njobs >= 0.0) || std::isinf(njobs))
            throw InputError("ClosedClass '" + nm + "': the population must be finite");
        JobClass cl;
        cl.name = nm;
        cl.type = JobClassType::CLOSED;
        cl.population = njobs;
        cl.refstat = station_of(refstat_node, "ClosedClass '" + nm + "'");
        cl.prio = prio;
        const std::size_t r = sn_.add_class(cl);
        grow_class_vectors();
        return r;
    }

    /**
     * An open class. Its reference station is the Source, which must exist:
     * an open class with no arrival station has no reference for its visits.
     */
    std::size_t add_open_class(const std::string& nm, int prio = 0) {
        if (sn_.sourceIdx == 0)
            throw InputError("OpenClass '" + nm +
                             "': the model has no Source to reference; add one first");
        JobClass cl;
        cl.name = nm;
        cl.type = JobClassType::OPEN;
        cl.population = std::numeric_limits<double>::infinity();
        cl.refstat = sn_.sourceIdx;
        cl.prio = prio;
        const std::size_t r = sn_.add_class(cl);
        grow_class_vectors();
        return r;
    }

    /**
     * `SelfLoopingClass(model, name, njobs, refstat, prio)`: a closed class
     * whose jobs perpetually cycle at their reference station.
     *
     * Built as a closed class because that is all it is -- the self-loop is in
     * the routing, and `SelfLoopingClass.m` adds no state -- with the subclass
     * recorded so the wire type survives a round trip.
     */
    std::size_t add_self_looping_class(const std::string& nm, double njobs,
                                       std::size_t refstat_node, int prio = 0) {
        const std::size_t r = add_closed_class(nm, njobs, refstat_node, prio);
        sn_.classes[r - 1].self_looping = true;
        return r;
    }

    /** `JobClass.setReferenceClass(true)`: `sn.refclass(c)` picks this class. */
    void set_reference_class(std::size_t cls) {
        class_ref(cls, "set_reference_class").is_ref_class = true;
    }

    /** `JobClass.deadline`: the soft deadline EDD, EDF and JMT's tardiness use. */
    void set_class_deadline(std::size_t cls, double due) {
        class_ref(cls, "set_class_deadline").deadline = due;
    }

    /**
     * `JobClass.spawnClass` (`sn.classspawn`): the class injected at the same
     * station on every completion of `cls`.
     */
    void set_class_spawn(std::size_t cls, std::size_t spawn_cls) {
        const std::size_t K = sn_.classes.size();
        if (spawn_cls == 0 || spawn_cls > K)
            throw InputError("set_class_spawn: the spawned class index is out of range");
        class_ref(cls, "set_class_spawn").spawn = spawn_cls;
    }

    /**
     * `JobClass.setPatience(kind, dist)`: the CLASS-WIDE abandonment law.
     *
     * The reference has no class-indexed patience in `sn`: `refreshStruct`
     * reads it through `Queue.getPatience`, which falls back to the class
     * setting wherever the station declares none, so the class-level law is
     * materialized onto every Queue and Delay here for exactly that reason.
     * A station-level `setPatience` therefore wins, as it does in MATLAB.
     */
    void set_class_patience(std::size_t cls, const Distrib<T>& dist,
                            lang::ImpatienceType kind = lang::ImpatienceType::RENEGING) {
        const std::size_t K = sn_.classes.size();
        if (cls == 0 || cls > K) throw InputError("set_class_patience: class index out of range");
        for (std::size_t ist = 0; ist < sn_.stations.size(); ++ist) {
            const std::size_t ind = sn_.station_to_node[ist];
            const lang::NodeType nt = sn_.nodes[ind - 1].nodetype;
            if (nt != lang::NodeType::Queue && nt != lang::NodeType::Delay) continue;
            Station<T>& st = sn_.stations[ist];
            grow_class_slot(st.patience, cls, Distrib<T>::disabled_dist());
            grow_class_slot(st.impatience, cls, lang::ImpatienceType::NONE);
            if (!st.patience[cls - 1].disabled) continue;  // the station setting wins
            st.patience[cls - 1] = dist;
            st.impatience[cls - 1] = kind;
        }
    }

    /**
     * `JobClass.setReplySignalClass(reply)` (`sn.syncreply`), plus the
     * `sn.replyblock` the state layer needs.
     *
     * The reference derives the block rather than being told it
     * (`refreshLocalVars.m:340-386`): a server is held at every non-Source,
     * non-INF station the REPLY class can be routed INTO, and every such
     * station must be FCFS because a held server is encoded as a per-class
     * counter. That derivation needs the routing, so it runs in `finalize()`;
     * this only records the binding.
     */
    void set_reply_signal_class(std::size_t call_cls, std::size_t reply_cls) {
        const std::size_t K = sn_.classes.size();
        if (call_cls == 0 || call_cls > K || reply_cls == 0 || reply_cls > K)
            throw InputError("set_reply_signal_class: class index is out of range");
        if (sn_.syncreply.size() < K) sn_.syncreply.assign(K, 0);
        sn_.syncreply[call_cls - 1] = reply_cls;
        set_signal(reply_cls, lang::SignalType::REPLY);
    }

    // -----------------------------------------------------------------------
    // Station parameters
    // -----------------------------------------------------------------------

    /**
     * `station.setService(class, dist)`.
     *
     * A Prior gets its MIXTURE moments here, where a Markovian family gets the
     * moments of its (D0,D1): both are the "what does the struct report before
     * anything solves" question, and leaving a Prior at mean 0 would make a
     * struct dump read as an Immediate.
     */
    void set_service(std::size_t node, std::size_t cls, const Distrib<T>& d) {
        Distrib<T> dd = d;
        if (dd.is_prior())
            lang::prior_refresh_moments(dd);
        else
            lang::dist_refresh_moments(dd);
        sn_.set_service(station_of(node, "setService"), cls, dd);
    }

    /** `source.setArrival(class, dist)`: the same table, at the Source. */
    void set_arrival(std::size_t node, std::size_t cls, const Distrib<T>& d) {
        const std::size_t ist = station_of(node, "setArrival");
        if (sn_.stations[ist - 1].nodetype != NodeType::Source)
            throw InputError("setArrival: node '" + sn_.nodes[node - 1].name + "' is not a Source");
        Distrib<T> dd = d;
        if (dd.is_prior())
            lang::prior_refresh_moments(dd);
        else
            lang::dist_refresh_moments(dd);
        sn_.set_service(ist, cls, dd);
    }

    /**
     * `queue.setNumberOfServers(n)`.
     *
     * IT IS A NO-OP ON AN INF-SCHEDULED STATION, which is what MATLAB does: the
     * method switches on the discipline and ignores the request for
     * SchedStrategy.INF. Lowering the multiplicity onto the station instead
     * looks harmless and is not -- utilization at a finite-server station is
     * divided by the server count, so an inf-scheduled station would report a
     * utilization a factor `n` too small.
     */
    void set_number_of_servers(std::size_t node, double n) {
        const std::size_t ist = station_of(node, "setNumberOfServers");
        if (sn_.stations[ist - 1].sched == SchedStrategy::INF) return;
        if (!(n >= 1.0)) throw InputError("setNumberOfServers: the server count must be >= 1");
        sn_.stations[ist - 1].nservers = n;
        // getNodeTypes reads the COUNT: a queue given infinitely many servers is
        // reported as a Delay station, which is the field the analyzers partition on
        if (std::isinf(n) && sn_.stations[ist - 1].nodetype == NodeType::Queue)
            sn_.stations[ist - 1].nodetype = NodeType::Delay;
    }

    /** `station.setCapacity(k)`, the K of Kendall's notation. */
    void set_capacity(std::size_t node, double k) {
        sn_.stations[station_of(node, "setCapacity") - 1].cap = k;
    }

    /** `station.setChainCapacity(class, k)`. */
    void set_class_capacity(std::size_t node, std::size_t cls, double k) {
        Station<T>& st = sn_.stations[station_of(node, "setChainCapacity") - 1];
        st.classcap.resize(sn_.classes.size(), std::numeric_limits<double>::infinity());
        st.classcap[cls - 1] = k;
    }

    /**
     * `queue.setImmediateFeedback(class)`: a completing job of that class is fed
     * straight back into service, HOLDING THE SERVER, rather than being routed
     * out and re-queued.
     *
     * Node-scoped. `set_class_immediate_feedback` is the class-wide spelling;
     * `sn.immfeed` is the OR of the two, as `refreshStruct` computes it.
     */
    void set_immediate_feedback(std::size_t node, std::size_t cls) {
        Station<T>& st = sn_.stations[station_of(node, "setImmediateFeedback") - 1];
        st.immfeed.resize(sn_.classes.size(), false);
        st.immfeed[cls - 1] = true;
    }

    /** `jobclass.setImmediateFeedback()`: the same property, class-wide. */
    void set_class_immediate_feedback(std::size_t cls) {
        sn_.classes[cls - 1].immfeed = true;
    }

    /** `station.setDropRule(class, rule)`. */
    void set_drop_rule(std::size_t node, std::size_t cls, DropStrategy rule) {
        Station<T>& st = sn_.stations[station_of(node, "setDropRule") - 1];
        st.droprule.resize(sn_.classes.size(), 0);
        st.droprule[cls - 1] = static_cast<int>(rule);
    }

    /**
     * `Queue.setServiceRateFunction(muFun)`: the TOTAL service rate of a PAS or
     * OI station as a function of the ordered microstate, a 1-based list of
     * class indices in queue order.
     *
     * Only PAS and OI take one, and the reference errors on any other
     * discipline. As `Queue.setServiceRateFunction` does, this ALSO installs a
     * representative per-class service distribution `Exp(mu([r]))`, so that the
     * ordinary rate/procid machinery stays consistent; the authoritative
     * description of the station remains mu(c). A class whose mu([r]) is not
     * positive and finite is disabled there, again as the reference does.
     *
     * `swap_graph` is `sn.nodeparam{ind}.swapGraph`, empty (all zero) for a
     * genuinely order-independent station.
     */
    void set_service_rate_function(
        std::size_t node, const std::function<T(const std::vector<std::size_t>&)>& muFun,
        const Matrix<T>& swap_graph = Matrix<T>()) {
        const std::size_t ist = station_of(node, "setServiceRateFunction");
        Station<T>& st = sn_.stations[ist - 1];
        if (st.sched != SchedStrategy::PAS && st.sched != SchedStrategy::OI)
            throw InputError(
                "setServiceRateFunction is only applicable to PAS (pass-and-swap) and OI "
                "(order-independent) queues");
        if (!muFun) throw InputError("setServiceRateFunction: the rate function is empty");
        st.svc_rate_fun = muFun;
        st.swap_graph = swap_graph;
        // ONE DECLARATION, TWO READERS. `solver_nc_oi` and `solver_mva_oi` read
        // the rate function off the Station; `to_marginal`, the PAS event
        // handler and the JSON writer read it off `sn.pasparam`. Filling only
        // one left the other with no rate function AND no error: a model built
        // here fell back to "every job present is in service" in every
        // state-space consumer, and a model loaded from JSON (which fills
        // `pasparam` alone, via `set_pas`) was refused by the OI solvers.
        pas_mirror(ist, muFun, swap_graph);
        for (std::size_t r = 1; r <= sn_.classes.size(); ++r) {
            const T rate_r = muFun(std::vector<std::size_t>{r});
            const double v = num_traits<T>::to_double(rate_r);
            set_service(node, r, (std::isfinite(v) && v > 0.0) ? Distrib<T>::exp_rate(rate_r)
                                                              : Distrib<T>::disabled_dist());
        }
    }

    /**
     * `Queue.setPollingType(rule, par)`: the polling discipline of a POLLING
     * station, identical across all class buffers as the reference assumes. Only
     * K-limited carries a parameter; every other rule ignores it.
     */
    void set_polling_type(std::size_t node, lang::PollingType rule, int par = 0) {
        Station<T>& st = sn_.stations[station_of(node, "setPollingType") - 1];
        if (st.sched != SchedStrategy::POLLING)
            throw InputError("setPollingType is only applicable to a POLLING station");
        if (rule == lang::PollingType::KLIMITED && par < 1)
            throw InputError("K-limited polling requires a parameter K >= 1");
        st.polling_type.assign(sn_.classes.size(), rule);
        st.polling_par = (rule == lang::PollingType::KLIMITED) ? par : 0;
    }

    /** `Queue.setSwitchover(jobclass, distrib)`: the switchover time of a class. */
    void set_switchover(std::size_t node, std::size_t cls, const Distrib<T>& so) {
        Station<T>& st = sn_.stations[station_of(node, "setSwitchover") - 1];
        if (st.sched != SchedStrategy::POLLING)
            throw InputError("setSwitchover is only applicable to a POLLING station");
        st.switchover.resize(sn_.classes.size(), Distrib<T>::immediate());
        st.switchover[cls - 1] = so;
    }

    /** The DPS / GPS weight of a class at a station. */
    void set_sched_param(std::size_t node, std::size_t cls, const T& weight) {
        Station<T>& st = sn_.stations[station_of(node, "setSchedParam") - 1];
        st.schedparam.resize(sn_.classes.size(), num_traits<T>::from_int(1));
        st.schedparam[cls - 1] = weight;
    }

    /**
     * `station.setLoadDependence(alpha)`: the rate multiplier at population
     * 1, 2, ... The vector is indexed from population one, as `sn.lldscaling`
     * is, so entry 0 is the multiplier of a station holding one job.
     */
    void set_load_dependence(std::size_t node, const std::vector<T>& alpha) {
        if (alpha.empty()) throw InputError("setLoadDependence: the scaling vector is empty");
        sn_.stations[station_of(node, "setLoadDependence") - 1].lldscaling = alpha;
    }

    /**
     * `station.setClassDependence(beta, peakRatePerClass)`.
     *
     * The peak is a scalar broadcast across the classes, or one value per class;
     * it is the declared max_n beta_r(n) that utilization is normalized by, and
     * the reference makes it mandatory because it cannot be recovered from beta
     * without sweeping the whole lattice -- a sweep that needs a bound the
     * handle does not carry, and an open class has none.
     *
     * An empty peak is accepted HERE and refused where it is READ, matching
     * `getLimitedClassDependencePeak`. That contract is only worth anything if
     * every reader honours it, and they do: SolverCTMC, `solver_mva_run_analyzer`,
     * `solver_nc_conv`, both SSA engines and the LDES engine each throw by name.
     * SolverMVA was the exception until 2026-08-19, silently writing a column of
     * ZEROS into U instead. `set_joint_dependence` below takes the peak as a
     * REQUIRED argument and refuses at declaration; both shapes end in an error.
     */
    void set_class_dependence(std::size_t node, const CdScaling<T>& fun,
                              const std::vector<T>& peak = std::vector<T>()) {
        Station<T>& st = sn_.stations[station_of(node, "setClassDependence") - 1];
        const std::size_t K = sn_.classes.size();
        st.cdscaling = fun;
        if (peak.empty())
            st.cdscalingpeak.clear();
        else if (peak.size() == 1)
            st.cdscalingpeak.assign(K, peak[0]);
        else if (peak.size() == K)
            st.cdscalingpeak = peak;
        else
            throw InputError(
                "setClassDependence: peakRatePerClass must be a scalar or a vector of length "
                "nclasses");
    }

    /**
     * `station.setJointDependence(eta, peakRatePerClass)`: MATLAB's
     * `Station.ljdScaling` / `ljdScalingPeak`.
     *
     * The peak is MANDATORY, exactly as in `Station.setJointDependence`, and for
     * the same reason as `setClassDependence`: utilization at a dependent station
     * is reported as T*S/peak, and max_n eta_i(n) is not recoverable from the
     * handle without sweeping the whole lattice.
     */
    void set_joint_dependence(std::size_t node, const CdScaling<T>& fun,
                              const std::vector<T>& peak) {
        Station<T>& st = sn_.stations[station_of(node, "setJointDependence") - 1];
        const std::size_t K = sn_.classes.size();
        if (peak.empty())
            throw InputError(
                "setJointDependence: joint dependence requires an explicit peak rate; pass a "
                "scalar (identical peak for every class) or a per-class vector");
        st.jdscaling = fun;
        if (peak.size() == 1)
            st.jdscalingpeak.assign(K, peak[0]);
        else if (peak.size() == K)
            st.jdscalingpeak = peak;
        else
            throw InputError(
                "setJointDependence: peakRatePerClass must be a scalar or a vector of length "
                "nclasses");
    }

    /**
     * `model.setGlobalDependence(phi, peak)`: MATLAB's `Network.gdScaling`.
     *
     * Declares a globally state-dependent rate scaling phi(n) whose argument is
     * the FULL (nstations x nclasses) population matrix, row-major, rather than
     * one station's slice. This is the Whittle primitive: when phi satisfies
     * phi_s(n) phi_t(n-e_s) = phi_t(n) phi_s(n-e_t) the chain is reversible with
     * pi(n) ~ Phi(n) prod rho_s^n_s and is insensitive; it also expresses
     * bandwidth sharing, where a route holds several links at once.
     *
     * phi returns one scalar (broadcast), one entry per station, or one entry per
     * (station, class) in row-major order. The peak is MANDATORY for the same
     * reason as `set_class_dependence`: utilization is reported as T*S/peak.
     * Only SolverCTMC honours the handle.
     */
    void set_global_dependence(const GdScaling<T>& fun, const std::vector<T>& peak) {
        set_global_dependence(fun, peak, 10);
    }

    /**
     * As above, with an explicit per-slot OPEN-class truncation used when phi is
     * materialized onto the JSON wire (closed classes are tabulated up to their own
     * population). It plays no part in solving, and exists because a handle cannot
     * cross a language boundary: the writer needs to know how far the lattice
     * extends. Set it to the cutoff the model is solved at.
     */
    void set_global_dependence(const GdScaling<T>& fun, const std::vector<T>& peak,
                               int wire_cutoff) {
        if (!fun)
            throw InputError("setGlobalDependence: the scaling must be a callable");
        const std::size_t M = sn_.stations.size(), K = sn_.classes.size();
        if (peak.empty())
            throw InputError(
                "setGlobalDependence: a global dependence requires an explicit peak rate; pass a "
                "scalar, one entry per station, or one entry per (station, class)");
        for (std::size_t i = 0; i < peak.size(); ++i)
            if (num_traits<T>::to_double(peak[i]) <= 0)
                throw InputError("setGlobalDependence: peak must be positive");
        // Probe now so a wrong output shape is refused at declaration time rather
        // than midway through state-space generation.
        for (int probe = 0; probe < 2; ++probe) {
            const std::vector<T> n(M * K, num_traits<T>::from_int(probe));
            const std::vector<T> v = fun(n);
            if (v.size() != 1 && v.size() != M && v.size() != M * K)
                throw InputError(
                    "setGlobalDependence: the handle must return a scalar, one entry per station, "
                    "or one entry per (station, class)");
            for (std::size_t j = 0; j < v.size(); ++j)
                if (!(num_traits<T>::to_double(v[j]) >= 0))
                    throw InputError(
                        "setGlobalDependence: the handle must return finite nonnegative scalings");
        }
        if (wire_cutoff < 1)
            throw InputError("setGlobalDependence: wireCutoff must be a positive integer");
        sn_.gdscaling = fun;
        sn_.gdscalingcutoff = wire_cutoff;
        if (peak.size() == 1)
            sn_.gdscalingpeak.assign(M * K, peak[0]);
        else if (peak.size() == M) {
            sn_.gdscalingpeak.assign(M * K, num_traits<T>::from_int(1));
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < K; ++r) sn_.gdscalingpeak[i * K + r] = peak[i];
        } else if (peak.size() == M * K)
            sn_.gdscalingpeak = peak;
        else
            throw InputError(
                "setGlobalDependence: peak must be a scalar, one entry per station, or one entry "
                "per (station, class)");
    }

    /** `node.setRouting(class, strategy)`. */
    void set_routing(std::size_t node, std::size_t cls, RoutingStrategy rs) {
        NodeDef& nd = sn_.nodes[node - 1];
        nd.routing.resize(sn_.classes.size(), RoutingStrategy::PROB);
        nd.routing[cls - 1] = rs;
    }

    /**
     * `node.setStateDepRouting(class, departure, branches, level, C, d)`.
     *
     * Declares `entry` the entry centre e of a subnetwork Q(V,V) served by the
     * product-form state-dependent routing of A. E. Krzesinski, "Multiclass
     * Queueing Networks with State-Dependent Routing", Performance Evaluation
     * 7(2):125-143, 1987.
     *
     * All node indices are 1-based, as elsewhere in the builder. `departure`
     * may equal `entry` in a central server model. `branches` follows the
     * paper's own indexing: `branches[0]` must be empty because branch index 1
     * denotes the complement M-V, and `branches[b]` lists the nodes of branch b
     * with its entry centre first and its departure centre last. `level[b]` is
     * the index t of the subnetwork with B_b in V_t - V_{t+1}, and `level[0]`
     * is ignored. `C` holds the T coefficients C_t and `d` the T by B
     * coefficients d_tb, read for 1 <= t <= level[b].
     *
     * Negative C_t and positive d_tb make the routing prefer the least
     * congested branches and impose the population bounds m_b <= d_tb/(-C_t)
     * and v_t <= D_tt/(-C_t). The residual probability returns the customer to
     * the departure centre, the busy form of waiting of Sec. 2.5, so the entry
     * node needs a self-loop when entry and departure coincide.
     */
    void set_state_dep_routing(std::size_t entry, std::size_t departure,
                               const std::vector<std::vector<std::size_t>>& branches,
                               const std::vector<std::size_t>& level,
                               const std::vector<double>& C, const Matrix<double>& d,
                               std::size_t cls = 0) {
        if (branches.size() < 2 || !branches[0].empty())
            throw InputError("set_state_dep_routing: branches[0] must be empty, branch index 1 "
                             "denotes the complement M-V");
        const std::size_t B = branches.size();
        if (level.size() != B)
            throw InputError("set_state_dep_routing: level must have one entry per branch index, "
                             "including the unused index 0");
        pfqn::SdrStruct nodesdr;
        nodesdr.entry = entry - 1;
        nodesdr.departure = departure - 1;
        nodesdr.branch.assign(B, std::vector<std::size_t>());
        nodesdr.entryOf.assign(B, 0);
        nodesdr.departureOf.assign(B, 0);
        for (std::size_t b = 1; b < B; ++b) {
            if (branches[b].empty())
                throw InputError("set_state_dep_routing: branch " + std::to_string(b + 1) +
                                 " is empty");
            for (std::size_t k = 0; k < branches[b].size(); ++k)
                nodesdr.branch[b].push_back(branches[b][k] - 1);
            nodesdr.entryOf[b] = branches[b].front() - 1;
            nodesdr.departureOf[b] = branches[b].back() - 1;
        }
        nodesdr.level = level;
        nodesdr.C = C;
        nodesdr.d = d;

        // Station-indexed twin. Every centre of an SDR network must be a
        // station: the product form is over queue lengths, and a stateless node
        // holds none.
        std::vector<std::size_t> node_to_station(sn_.nodes.size(), 0);
        std::vector<bool> is_station(sn_.nodes.size(), false);
        for (std::size_t k = 0; k < sn_.station_to_node.size(); ++k) {
            node_to_station[sn_.station_to_node[k] - 1] = k;
            is_station[sn_.station_to_node[k] - 1] = true;
        }
        struct Map {
            const std::vector<std::size_t>& n2s;
            const std::vector<bool>& isst;
            const NetworkStruct<T>& sn;
            std::size_t operator()(std::size_t nd) const {
                if (nd >= isst.size() || !isst[nd])
                    throw InputError("set_state_dep_routing: node '" + sn.nodes[nd].name +
                                     "' takes part in state-dependent routing but is not a "
                                     "station: the product form is over queue lengths, and a "
                                     "stateless node holds none");
                return n2s[nd];
            }
        } to_station{node_to_station, is_station, sn_};

        pfqn::SdrStruct stsdr = nodesdr;
        stsdr.entry = to_station(nodesdr.entry);
        stsdr.departure = to_station(nodesdr.departure);
        for (std::size_t b = 1; b < B; ++b) {
            for (std::size_t k = 0; k < nodesdr.branch[b].size(); ++k)
                stsdr.branch[b][k] = to_station(nodesdr.branch[b][k]);
            stsdr.entryOf[b] = to_station(nodesdr.entryOf[b]);
            stsdr.departureOf[b] = to_station(nodesdr.departureOf[b]);
        }
        pfqn::pfqn_sdrcoeff(stsdr);  // validates the declaration and its population bounds

        sn_.sdr_nodes = nodesdr;
        sn_.sdr = stsdr;
        const std::size_t K = sn_.classes.size();
        if (cls == 0) {
            for (std::size_t r = 1; r <= K; ++r) set_routing(entry, r, RoutingStrategy::SDR);
        } else {
            set_routing(entry, cls, RoutingStrategy::SDR);
        }
    }

    // Routing

    /** An empty routing matrix, MATLAB's `model.initRoutingMatrix`. */
    RoutingMatrix<T> init_routing_matrix() const { return RoutingMatrix<T>(); }

    /**
     * `model.link(P)`: install the routing.
     *
     * The probabilities are stored as given. A node whose strategy is RAND
     * needs only the CONNECTIONS -- any positive entry marks one -- and the
     * refresh replaces them by the uniform split.
     */
    void link(const RoutingMatrix<T>& Pm) {
        for (const auto& kv : Pm.entries) {
            const auto& k = kv.first;
            if (k.r == 0 || k.r > sn_.classes.size() || k.s == 0 || k.s > sn_.classes.size())
                throw InputError("link: the routing matrix names a class that does not exist");
            if (k.i == 0 || k.i > sn_.nodes.size() || k.j == 0 || k.j > sn_.nodes.size())
                throw InputError("link: the routing matrix names a node that does not exist");
        }
        // A CLASS SWITCH ON A LINK BECOMES A NODE, as `@MNetwork/link.m:225-329`
        // makes it. This port used to keep `P{r,s}(i,j)` with r != s as an EDGE
        // attribute and synthesize nothing, which the analytical solvers read
        // correctly through route_eff but which cost two things nothing else
        // could supply: the node table had no CS_ row to report (a node the
        // reference counts, so cache_replc_fifo showed 2 nodes against 3), and
        // the JSIM export had nowhere to put the switch at all -- `jmt_writer.h`
        // can only emit a ClassSwitch for a node whose type IS ClassSwitch, so
        // JMT silently simulated the UNSWITCHED model.
        //
        // The rewrite makes every surviving route SAME-CLASS; the switching
        // lives entirely in the inserted node's matrix.
        const std::size_t K = sn_.classes.size();
        const std::size_t I = sn_.nodes.size();  // before any insertion
        const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
        RoutingMatrix<T> P = Pm;
        // Deferred miss hops from Cache.setMissCache. They are part of the cache
        // network's own construction, and the per-item classes they carry do not
        // exist when the caller builds its routing matrix.
        for (std::size_t a = 0; a < cache_miss_arcs_.size(); ++a) {
            const CacheMissArc& mc = cache_miss_arcs_[a];
            P.set(mc.cls, mc.cls, mc.from_node, mc.to_node, one);
        }
        std::map<std::pair<std::size_t, std::size_t>, std::size_t> csid;  // (i,j) -> 1-based node
        for (std::size_t i = 1; i <= I; ++i)
            for (std::size_t j = 1; j <= I; ++j) {
                Matrix<T> C(K, K, zero);
                bool any = false;
                for (std::size_t r = 1; r <= K; ++r)
                    for (std::size_t s = 1; s <= K; ++s) {
                        const T p = P.get(r, s, i, j);
                        if (num_traits<T>::to_double(p) != 0.0) any = true;
                        C(r - 1, s - 1) = p;
                    }
                if (!any) continue;  // no link at all: the reference's identity, no node
                // CONDITIONED ON THE JOB TAKING THIS LINK, so each row is
                // renormalised; a class that never leaves i for j keeps itself.
                bool offdiag = false;
                for (std::size_t r = 1; r <= K; ++r) {
                    T S = zero;
                    for (std::size_t s = 1; s <= K; ++s) S += C(r - 1, s - 1);
                    if (num_traits<T>::to_double(S) > 0) {
                        for (std::size_t s = 1; s <= K; ++s) C(r - 1, s - 1) = T(C(r - 1, s - 1) / S);
                    } else {
                        C(r - 1, r - 1) = one;
                    }
                    for (std::size_t s = 1; s <= K; ++s)
                        if (s != r && num_traits<T>::to_double(C(r - 1, s - 1)) != 0.0) offdiag = true;
                }
                if (!offdiag) continue;  // `~isdiag`: a pure same-class link needs no node
                // THE SINK -> SOURCE ARC IS THE OPEN-NETWORK CLOSURE, NOT A LINK.
                // A job cannot be routed out of a Sink, so an arc from one into a
                // Source is never something a user asked for: it is what makes an
                // open chain irreducible for the visit computation, and the
                // reference adds it in refreshStruct AFTER link.m has run, so
                // link.m never sees it and synthesizes nothing for it. This port
                // reads a model.json where linemodel_save has already written the
                // closure in, and where it changes class -- an open model that
                // enters as InitClass and leaves as HitClass closes as
                // `P{HitClass,InitClass}(Sink,Source)` -- the class change looked
                // like a switch on a link and grew a CS_Sink_to_Source node. The
                // extra hop HALVED every downstream visit: on cache_replc_routing
                // both Delays reported throughput 0.2/0.3 against MATLAB's
                // 0.4/0.6, in every engine at once. The route is still installed
                // below; only the node is not.
                if (sn_.nodes[i - 1].nodetype == NodeType::Sink &&
                    sn_.nodes[j - 1].nodetype == NodeType::Source)
                    continue;
                csid[std::make_pair(i, j)] =
                    add_class_switch("CS_" + sn_.nodes[i - 1].name + "_to_" + sn_.nodes[j - 1].name, C);
            }
        // RE-ROUTE i -> cs -> j, the reference's own three assignments. The
        // switching probability is folded into the FIRST leg in the departing
        // class, and the second leg is deterministic in the arriving class.
        for (const auto& kv : csid) {
            const std::size_t i = kv.first.first, j = kv.first.second, c = kv.second;
            for (std::size_t r = 1; r <= K; ++r)
                for (std::size_t s = 1; s <= K; ++s) {
                    const T p = P.get(r, s, i, j);
                    if (!(num_traits<T>::to_double(p) > 0)) continue;
                    P.set(r, r, i, c, T(P.get(r, r, i, c) + p));
                    P.set(r, s, i, j, zero);
                    P.set(s, s, c, j, one);
                }
        }
        for (const auto& kv : P.entries) {
            const auto& k = kv.first;
            sn_.set_route(k.r, k.s, k.i, k.j, kv.second);
        }
        // A retrieval class inherits the read class's routing among the cache and
        // the retrieval queues: it circulates them the same way, entering at the
        // cache and returning to it. Set once P is installed (setRetrievalSystem
        // runs before link, as the reference documents).
        for (const auto& np : sn_.nodeparam) {
            const std::size_t ci = np.first;  // 1-based cache node
            const CacheParam<T>& cp = np.second;
            if (cp.retrieval_capacity <= 0) continue;
            for (const auto& rq : cp.retrieval_queues) {
                const std::size_t rd = rq.first + 1;  // 1-based read class
                std::vector<std::size_t> nodeset = rq.second;  // 1-based queue nodes
                nodeset.push_back(ci);
                bool minted = false;
                for (std::size_t i = 0; i < cp.nitems; ++i) {
                    const std::size_t rc = cp.retrieval_classes[i][rd - 1];
                    if (rc == 0) continue;
                    minted = true;
                    for (std::size_t a : nodeset)
                        for (std::size_t b : nodeset) {
                            const T p = Pm.get(rd, rd, a, b);
                            if (p > num_traits<T>::from_int(0)) sn_.set_route(rc, rc, a, b, p);
                        }
                }
                // CONSUME the read class's template edges over the queue set,
                // as `link.m:187-194` does. They were only ever a TEMPLATE from
                // which each retrieval class's circulation is copied; leaving
                // them in place makes the read class circulate the fetch
                // stations on its own, so each retrieval class becomes a
                // separate communicating class and the chain decomposition
                // reports one singleton chain per item instead of one chain
                // over all of them. Measured before the fix on a closed
                // Delay+Cache+Fetch model: 4 chains against MATLAB's 1, and
                // every visit at the fetch station zero against MATLAB's 1.5
                // per retrieval class.
                if (minted)
                    for (std::size_t a : nodeset)
                        for (std::size_t b : nodeset)
                            if (Pm.get(rd, rd, a, b) > num_traits<T>::from_int(0))
                                sn_.set_route(rd, rd, a, b, num_traits<T>::from_int(0));
            }
        }

        // A (node, class) PAIR THE CALLER NEVER ROUTED IS LEFT ON RAND, which is
        // where `addLink` leaves it in the reference: `setProbRouting` is what
        // turns a pair PROB, and it is called only for the entries the routing
        // matrix actually carries. The distinction is not bookkeeping.
        //
        // It decides the VISITS. `getRoutingMatrix` expands RAND uniformly over
        // the node's connections, so the pair gets a row; PROB with no entries
        // leaves the row EMPTY, and an empty row is an UNVISITED station. On
        // `gallery_erlerl1` -- Class1 declares Source->Queue, Class2 declares
        // Queue->Sink -- the Queue came out unvisited by Class1 and every metric
        // of the model was zero. The `served` mask in `refresh_visits` is what
        // keeps the fill honest: it drops the (station, class) pairs with no
        // service law.
        //
        // It also decides the SAMPLE PATH. `saveRoutingStrategy` writes a
        // RandomStrategy for a RAND pair and an EmpiricalStrategy for a PROB
        // one, and JMT draws from the node's stream for a random split EVEN AT
        // ONE DESTINATION -- so the two consume the stream differently and part
        // company at the same seed. On `fj_cs_prefork` that was a 14% gap on
        // Queue1 against a golden JMT itself produced.
        //
        // A SELF-LOOPING CLASS AND A SIGNAL ARE LEFT ALONE: the first never
        // leaves its station, and the second is routed explicitly by P.
        const std::size_t Kc = sn_.classes.size(), Ic = sn_.nodes.size();
        for (std::size_t i = 1; i <= Ic; ++i) {
            NodeDef& nd = sn_.nodes[i - 1];
            nd.routing.resize(Kc, RoutingStrategy::PROB);
            bool linked = false;
            for (std::size_t j = 1; j <= Ic && !linked; ++j)
                for (std::size_t a = 1; a <= Kc && !linked; ++a)
                    for (std::size_t b = 1; b <= Kc && !linked; ++b)
                        if (num_traits<T>::to_double(sn_.get_route(a, b, i, j)) > 0.0) linked = true;
            if (!linked) continue;
            for (std::size_t r = 1; r <= Kc; ++r) {
                if (nd.routing[r - 1] != RoutingStrategy::PROB) continue;
                if (sn_.classes[r - 1].self_looping) continue;
                if (r <= sn_.issignal.size() && sn_.issignal[r - 1]) continue;
                bool routed = false;
                for (std::size_t j = 1; j <= Ic && !routed; ++j)
                    for (std::size_t s = 1; s <= Kc && !routed; ++s)
                        if (num_traits<T>::to_double(sn_.get_route(r, s, i, j)) > 0.0) routed = true;
                if (!routed) nd.routing[r - 1] = RoutingStrategy::RAND;
            }
        }

        // A variable forking level declared BEFORE the routing existed. Its
        // destination set is read off the routing, so this is the first moment
        // it can be resolved; a call made after `link()` was applied where it
        // stood and replaying it here is idempotent.
        apply_fork_overrides();
    }

    // -----------------------------------------------------------------------
    // Struct
    // -----------------------------------------------------------------------

    /** The refreshed struct, MATLAB's `model.getStruct()`. */
    const NetworkStruct<T>& get_struct() {
        if (routing_installed()) apply_fork_overrides();
        validate();
        sn_.refresh_struct();
        // The fork block of refreshStruct.m, and it belongs HERE rather than in
        // refresh_struct(): it needs the MMT transformation, which is built on a
        // refreshed struct, exactly as the reference calls ModelAdapter.mmt after
        // refreshChains. fj_tag's own refresh_struct() therefore skips it, which
        // is what `~isFJAugmented` buys the reference.
        api::sn_fj_nodevisits_mmt(sn_);
        return sn_;
    }

    /** The struct WITHOUT refreshing it, for a caller that is still building. */
    NetworkStruct<T>& raw_struct() { return sn_; }

    /**
     * `Queue.setService(@(c) ...)` for a pass-and-swap / order-independent
     * station: the total service rate mu(c) of an ordered list of 1-based
     * class indices, and the swap graph saying which class may take another's
     * place. An OI station is the special case of an empty swap graph.
     */
    void set_pas(std::size_t node,
                 const std::function<T(const std::vector<std::size_t>&)>& mu,
                 const std::vector<std::vector<bool>>& swap_graph =
                     std::vector<std::vector<bool>>()) {
        const std::size_t ist = station_of(node, "set_pas");
        typename NetworkStruct<T>::PasParam pp;
        pp.svc_rate_fun = mu;
        pp.swap_graph = swap_graph;
        if (pp.swap_graph.empty())
            pp.swap_graph.assign(sn_.classes.size(), std::vector<bool>(sn_.classes.size(), false));
        sn_.pasparam[ist] = pp;
        // The other half of the same declaration; see `set_service_rate_function`.
        Station<T>& stp = sn_.stations[ist - 1];
        stp.svc_rate_fun = mu;
        const std::size_t R = pp.swap_graph.size();
        Matrix<T> G(R, R, num_traits<T>::from_int(0));
        for (std::size_t a = 0; a < R; ++a)
            for (std::size_t b = 0; b < pp.swap_graph[a].size(); ++b)
                if (pp.swap_graph[a][b]) G(a, b) = num_traits<T>::from_int(1);
        stp.swap_graph = G;
    }

    /**
     * Copy a PAS/OI declaration into `sn.pasparam`, the form the state-space
     * layer reads. The (R x R) swap adjacency crosses as a boolean matrix
     * because that is the shape `after_event`'s swap walk indexes.
     */
    void pas_mirror(std::size_t ist,
                    const std::function<T(const std::vector<std::size_t>&)>& mu,
                    const Matrix<T>& swap_graph) {
        typename NetworkStruct<T>::PasParam pp;
        pp.svc_rate_fun = mu;
        const std::size_t R = sn_.classes.size();
        pp.swap_graph.assign(R, std::vector<bool>(R, false));
        for (std::size_t a = 0; a < R && a < swap_graph.rows(); ++a)
            for (std::size_t b = 0; b < R && b < swap_graph.cols(); ++b)
                pp.swap_graph[a][b] = num_traits<T>::to_double(swap_graph(a, b)) != 0.0;
        sn_.pasparam[ist] = pp;
    }

    /**
     * `Queue.setPollingType(...)`: the polling discipline of a POLLING station
     * and the switchover walks between its buffers.
     *
     * `switchover[r-1]` is the walk the server takes when LEAVING buffer r. An
     * Immediate walk is folded rather than represented, so it costs no state.
     */
    void set_polling(std::size_t node, lang::PollingType ptype,
                     const std::vector<Distrib<T>>& switchover = std::vector<Distrib<T>>(),
                     std::size_t pk = 1) {
        const std::size_t ist = station_of(node, "set_polling");
        if (sn_.stations[ist - 1].sched != SchedStrategy::POLLING)
            throw InputError("set_polling: the station is not POLLING-scheduled");
        typename NetworkStruct<T>::PollingParam pp;
        pp.ptype = ptype;
        pp.pk = pk;
        pp.switchover = switchover;
        if (pp.switchover.size() < sn_.classes.size())
            pp.switchover.resize(sn_.classes.size(), Distrib<T>::disabled_dist());
        sn_.pollingparam[ist] = pp;
    }

    /**
     * Declare a SYNCHRONOUS call: a job of `call_cls` leaving `node` keeps its
     * server until a job of the REPLY class `reply_cls` arrives back here.
     *
     * The block is one nvars column per (node, calling class), so it stays
     * zero-width -- and every other model's state width unchanged -- unless a
     * model actually declares a reply. `node` is checked here and marked, but
     * `refresh_replyblock` re-derives the whole block from the routing on every
     * refresh, exactly as the reference does; naming a node is therefore an
     * assertion about this model, not the definition of the block.
     */
    void set_sync_reply(std::size_t node, std::size_t call_cls, std::size_t reply_cls) {
        const std::size_t K = sn_.classes.size();
        if (node == 0 || node > sn_.nodes.size())
            throw InputError("set_sync_reply: node index is out of range");
        if (call_cls == 0 || call_cls > K || reply_cls == 0 || reply_cls > K)
            throw InputError("set_sync_reply: class index is out of range");
        const std::size_t ist = sn_.nodes[node - 1].station;
        if (ist == 0) throw InputError("set_sync_reply: node is not a station");
        if (sn_.stations[ist - 1].sched != SchedStrategy::FCFS)
            throw InputError(
                "set_sync_reply: synchronous calls are supported only at FCFS stations; "
                "holding a server across a call has no representation in the state of the "
                "other disciplines");
        if (sn_.replyblock.size() < sn_.nodes.size())
            sn_.replyblock.assign(sn_.nodes.size(), std::vector<bool>(K, false));
        if (sn_.syncreply.size() < K) sn_.syncreply.assign(K, 0);
        sn_.replyblock[node - 1][call_cls - 1] = true;
        sn_.syncreply[call_cls - 1] = reply_cls;
        set_signal(reply_cls, lang::SignalType::REPLY);
    }

    /**
     * `FiniteCapacityRegion(model, nodes)`: a cap on the jobs held ACROSS a set
     * of stations.
     *
     * `class_max_jobs[r]` and `global_max_jobs` are -1 for unbounded, the
     * reference's sentinel; the per-class cap is additionally tightened by
     * `floor(class_max_memory[r] / class_size[r])` when a memory budget is set,
     * because a class whose footprint exceeds the budget cannot have as many
     * jobs resident as its job cap alone would allow.
     *
     * @param nodes 1-based node indices; each must be a station
     * @param class_max_jobs (K) per-class job cap inside the region, -1 for unbounded
     * @param global_max_jobs cap on the total jobs inside the region, -1 for unbounded
     * @param rule (K) per-class drop strategy applied when the cap is reached
     * @param class_max_memory (K) per-class memory budget, -1 for unbounded
     * @param class_size (K) per-job memory footprint of each class
     * @param global_max_memory cap on the total memory inside the region, -1 for unbounded
     */
    std::size_t add_region(const std::vector<std::size_t>& nodes,
                           const std::vector<double>& class_max_jobs,
                           double global_max_jobs = -1.0,
                           const std::vector<DropStrategy>& rule = std::vector<DropStrategy>(),
                           const std::vector<double>& class_max_memory = std::vector<double>(),
                           const std::vector<T>& class_size = std::vector<T>(),
                           double global_max_memory = -1.0,
                           const std::string& name = std::string()) {
        const std::size_t M = sn_.stations.size(), K = sn_.classes.size();
        typename NetworkStruct<T>::Region rg;
        rg.name = name;
        rg.cap.assign(M, std::vector<double>(K + 1, -1.0));
        rg.maxmem.assign(M, -1.0);
        rg.members.assign(M, false);
        // A class beyond the region's own vectors defaults to WAITQ, weight 1
        // and size 1, exactly as refreshRegions does -- NOT to DROP, which
        // would silently start losing jobs of a class added after the region.
        rg.rule.assign(K, DropStrategy::WAITQ);
        rg.weight.assign(K, num_traits<T>::from_int(1));
        rg.size.assign(K, num_traits<T>::from_int(1));
        for (std::size_t r = 0; r < K && r < rule.size(); ++r) rg.rule[r] = rule[r];
        for (std::size_t r = 0; r < K && r < class_size.size(); ++r) rg.size[r] = class_size[r];

        for (std::size_t j = 0; j < nodes.size(); ++j) {
            const std::size_t ist = station_of(nodes[j], "addRegion");
            rg.members[ist - 1] = true;
            for (std::size_t r = 0; r < K; ++r) {
                // A class beyond the region's own job-cap vector is UNBOUNDED and
                // is not clamped by the memory budget either: refreshRegions.m
                // takes the `continue` before it reaches the memory block. Applying
                // the clamp here would cap a class MATLAB leaves free, whenever
                // class_max_memory is the longer of the two vectors.
                if (r >= class_max_jobs.size()) {
                    rg.cap[ist - 1][r] = -1.0;
                    continue;
                }
                double c = class_max_jobs[r];
                if (r < class_max_memory.size() && class_max_memory[r] != -1.0) {
                    const double sz = num_traits<T>::to_double(rg.size[r]);
                    if (sz > 0) {
                        const double memjobs = std::floor(class_max_memory[r] / sz);
                        c = c == -1.0 ? memjobs : std::min(c, memjobs);
                    }
                }
                rg.cap[ist - 1][r] = c;
            }
            rg.cap[ist - 1][K] = global_max_jobs;
            rg.maxmem[ist - 1] = global_max_memory;
        }
        sn_.regions.push_back(rg);
        return sn_.regions.size();
    }

    /**
     * `FiniteCapacityRegion.setClassWeight`: the per-class weight the region's
     * global cap counts a job against, defaulting to 1. A weight of 2 makes one
     * job of the class consume two of the region's slots.
     */
    void set_region_weights(std::size_t region, const std::vector<T>& weight) {
        if (region == 0 || region > sn_.regions.size())
            throw InputError("set_region_weights: region index is out of range");
        typename NetworkStruct<T>::Region& rg = sn_.regions[region - 1];
        for (std::size_t r = 0; r < rg.weight.size() && r < weight.size(); ++r)
            rg.weight[r] = weight[r];
    }

    /** The optional linear constraint A n <= b a region may carry beyond its caps. */
    void set_region_constraint(std::size_t region, const Matrix<T>& A, const std::vector<T>& b) {
        if (region == 0 || region > sn_.regions.size())
            throw InputError("set_region_constraint: region index is out of range");
        if (A.rows() != b.size())
            throw InputError("set_region_constraint: A and b disagree on the number of rows");
        sn_.regions[region - 1].lincon_A = A;
        sn_.regions[region - 1].lincon_b = b;
    }

    /**
     * `model.setReward(name, fn)`: a named reward evaluated on the AGGREGATE
     * state row, the per-(station, class) job counts in `(ist-1)*K + k` order.
     *
     * Redeclaring a name REPLACES it rather than adding a second reward under
     * the same name, so a caller refining a definition does not end up with two
     * answers labelled identically.
     */
    void set_reward(const std::string& nm,
                    const std::function<T(const std::vector<T>&)>& fn,
                    const std::string& kind = std::string(), std::size_t node = 0,
                    std::size_t cls = 0) {
        for (std::size_t i = 0; i < sn_.reward.size(); ++i)
            if (sn_.reward[i].name == nm) {
                sn_.reward[i].fn = fn;
                sn_.reward[i].kind = kind;
                sn_.reward[i].node = node;
                sn_.reward[i].cls = cls;
                return;
            }
        typename NetworkStruct<T>::Reward rw;
        rw.name = nm;
        rw.fn = fn;
        rw.kind = kind;
        rw.node = node;
        rw.cls = cls;
        sn_.reward.push_back(rw);
    }

    /**
     * Declare a class to be a G-network SIGNAL rather than a job.
     *
     * A signal never joins a station: it removes jobs already there and is
     * annihilated. `target` is the 1-based class it may remove, or 0 for the
     * classic untargeted Gelenbe customer. `remdist` is the batch-size pmf
     * indexed by batch size 0,1,2,...; empty means "remove exactly one".
     */
    void set_signal(std::size_t cls, lang::SignalType type,
                    lang::RemovalPolicy policy = lang::RemovalPolicy::RANDOM,
                    std::size_t target = 0,
                    const std::vector<T>& remdist = std::vector<T>()) {
        const std::size_t K = sn_.classes.size();
        if (cls == 0 || cls > K) throw InputError("set_signal: class index is out of range");
        if (target > K) throw InputError("set_signal: target class index is out of range");
        if (sn_.issignal.size() < K) {
            sn_.issignal.assign(K, false);
            sn_.signaltype.assign(K, lang::SignalType::NEGATIVE);
            sn_.signaltarget.assign(K, 0);
            sn_.signalrempolicy.assign(K, lang::RemovalPolicy::RANDOM);
            sn_.signalremdist.assign(K, std::vector<T>());
        }
        sn_.issignal[cls - 1] = true;
        sn_.signaltype[cls - 1] = type;
        sn_.signaltarget[cls - 1] = target;
        sn_.signalrempolicy[cls - 1] = policy;
        sn_.signalremdist[cls - 1] = remdist;
    }

    std::size_t station_index(std::size_t node) const { return station_of(node, "station_index"); }

  private:
    NetworkStruct<T> sn_;

    /** Per-item read classes of each cache in a cache network, keyed by cache node. */
    std::map<std::size_t, std::vector<std::size_t> > cache_item_classes_;

    /** A miss hop registered by `set_miss_cache`, injected into P by `link()`. */
    struct CacheMissArc {
        std::size_t cls, from_node, to_node;
        CacheMissArc(std::size_t c, std::size_t f, std::size_t t)
            : cls(c), from_node(f), to_node(t) {}
    };
    std::vector<CacheMissArc> cache_miss_arcs_;

    /** One class shared by every item, or one per item. */
    static std::vector<std::size_t> per_item_classes(const std::vector<std::size_t>& spec,
                                                     std::size_t nitems, const char* what) {
        if (spec.size() == nitems) return spec;
        if (spec.size() == 1) return std::vector<std::size_t>(nitems, spec[0]);
        throw InputError(std::string(what) +
                         ": pass one class per item or a single class shared by all");
    }

    /** The class record, with the index checked against the live class list. */
    JobClass& class_ref(std::size_t cls, const char* what) {
        if (cls == 0 || cls > sn_.classes.size())
            throw InputError(std::string(what) + ": class index is out of range");
        return sn_.classes[cls - 1];
    }

    void init_node(std::size_t nd) {
        sn_.nodes[nd - 1].routing.assign(sn_.classes.size(), RoutingStrategy::PROB);
    }

    /** Grow the per-class vectors of every node and station after a new class. */
    void grow_class_vectors() {
        const std::size_t K = sn_.classes.size();
        for (NodeDef& nd : sn_.nodes) nd.routing.resize(K, RoutingStrategy::PROB);
        for (Station<T>& st : sn_.stations) {
            if (!st.classcap.empty())
                st.classcap.resize(K, std::numeric_limits<double>::infinity());
            if (!st.droprule.empty()) st.droprule.resize(K, 0);
            if (!st.schedparam.empty()) st.schedparam.resize(K, num_traits<T>::from_int(1));
        }
    }

    /** The station of a node, with the class index checked against the live class list. */
    Station<T>& station_ref(std::size_t node, std::size_t cls, const char* what) {
        const std::size_t ist = station_of(node, what);
        if (cls == 0 || cls > sn_.classes.size())
            throw InputError(std::string(what) + ": class index is out of range");
        return sn_.stations[ist - 1];
    }

    /**
     * Widen an OPTIONAL per-class vector to hold `cls`, filling with the
     * "not declared" value. The vectors start empty on purpose -- an empty one
     * means the station declares none of this property at all -- so they cannot
     * be sized in `grow_class_vectors` without turning every station into one
     * that declares it.
     */
    template <class V>
    void grow_class_slot(std::vector<V>& v, std::size_t cls, const V& fill) {
        const std::size_t K = sn_.classes.size();
        if (v.size() < K) v.resize(K < cls ? cls : K, fill);
    }

    std::size_t station_of(std::size_t node, const std::string& what) const {
        if (node == 0 || node > sn_.nodes.size())
            throw InputError(what + ": node index " + std::to_string(node) + " does not exist");
        const std::size_t ist = sn_.nodes[node - 1].station;
        if (ist == 0)
            throw InputError(what + ": node '" + sn_.nodes[node - 1].name +
                             "' is not a station (it serves no jobs)");
        return ist;
    }

    /**
     * The fork's node record, with its variable-forking-level matrices sized on
     * first use.
     *
     * They start EMPTY on every fork, so a consumer can tell a plain fork from
     * one with overrides without inspecting entries; the first override is what
     * allocates them, seeded from `tasks_per_link` and probability 1 on the
     * links the model actually declares.
     */
    /** One recorded `Fork.set*` call, replayed once the routing is installed. */
    struct ForkOverride {
        enum Kind { TASKS, DIST, PROB };
        Kind kind = TASKS;
        std::size_t fork = 0, cls = 0, dest = 0;
        double value = 0.0;
        lang::Distrib<T> dist;
    };
    std::vector<ForkOverride> fork_overrides_;

    /**
     * Record an override, and apply it at once when the routing already exists.
     *
     * Applying eagerly is not an optimisation: it is what makes an override
     * naming a node that the fork does not reach fail AT THE CALL, where the
     * caller can see which line is wrong, rather than at `get_struct()`.
     */
    void record_fork_override(const ForkOverride& ov) {
        if (ov.fork == 0 || ov.fork > sn_.nodes.size() ||
            sn_.nodes[ov.fork - 1].nodetype != NodeType::Fork)
            throw InputError("the node given is not a Fork of this model");
        if (ov.cls == 0 || ov.cls > sn_.classes.size())
            throw InputError("a fork override names a class that does not exist");
        fork_overrides_.push_back(ov);
        if (routing_installed()) apply_fork_override(ov);
    }

    /** True once `link()` has written a routing block a fork override can read. */
    bool routing_installed() const {
        return sn_.rtnodes.rows() >= sn_.nodes.size() * sn_.classes.size();
    }

    /** Replay every recorded override, in the order the model declared them. */
    void apply_fork_overrides() {
        for (std::size_t i = 0; i < fork_overrides_.size(); ++i)
            apply_fork_override(fork_overrides_[i]);
    }

    void apply_fork_override(const ForkOverride& ov) {
        qn::ForkParam<T>& f = fork_param(ov.fork);
        const std::vector<std::size_t> dests = fork_dests(ov.fork, ov.dest);
        for (std::size_t x = 0; x < dests.size(); ++x) {
            const std::size_t k = dests[x];
            switch (ov.kind) {
                case ForkOverride::TASKS:
                    f.fan_out_link(k - 1, ov.cls - 1) = num_traits<T>::from_double(ov.value);
                    break;
                case ForkOverride::DIST:
                    f.fan_out_dist[k - 1][ov.cls - 1] = ov.dist;
                    // the scalar slot carries the mean, so a consumer that only
                    // reads fan_out_link still sees E[tasks per link]
                    f.fan_out_link(k - 1, ov.cls - 1) = ov.dist.mean;
                    break;
                case ForkOverride::PROB:
                    f.fan_out_prob(k - 1, ov.cls - 1) = num_traits<T>::from_double(ov.value);
                    break;
            }
        }
        refresh_fork_scalar(sn_.nodes[ov.fork - 1], f);
    }

    qn::ForkParam<T>& fork_param(std::size_t fork_node) {
        if (fork_node == 0 || fork_node > sn_.nodes.size() ||
            sn_.nodes[fork_node - 1].nodetype != NodeType::Fork)
            throw InputError("the node given is not a Fork of this model");
        qn::ForkParam<T>& f = sn_.forkparam[fork_node];
        const std::size_t I = sn_.nodes.size(), K = sn_.classes.size();
        if (f.fan_out_link.rows() == I && f.fan_out_link.cols() == K) return f;
        const T zero = num_traits<T>::from_int(0);
        f.fan_out_link = Matrix<T>(I, K, zero);
        f.fan_out_prob = Matrix<T>(I, K, zero);
        f.fan_out_dist.assign(I, std::vector<lang::Distrib<T> >(K));
        const T tpl = num_traits<T>::from_double(sn_.nodes[fork_node - 1].tasks_per_link);
        const T one = num_traits<T>::from_int(1);
        for (std::size_t k = 1; k <= I; ++k)
            for (std::size_t r = 1; r <= K; ++r)
                if (fork_links_to(fork_node, k, r)) {
                    f.fan_out_link(k - 1, r - 1) = tpl;
                    f.fan_out_prob(k - 1, r - 1) = one;
                }
        return f;
    }

    /**
     * Keep the scalar `tasks_per_link` consistent with the per-link mean, so a
     * solver that has not been taught the matrices degrades to E[tasks per
     * link] and not to a value the fork never emits.
     *
     * The branch probability is folded in HERE and not into `fan_out_link`,
     * because JMT and LDES read the two separately: `fan_out_link` is the count
     * GIVEN the branch fires, `fan_out_prob` is whether it fires at all.
     */
    void refresh_fork_scalar(qn::NodeDef& nd, const qn::ForkParam<T>& f) {
        double acc = 0.0;
        std::size_t cnt = 0;
        for (std::size_t k = 0; k < f.fan_out_link.rows(); ++k)
            for (std::size_t r = 0; r < f.fan_out_link.cols(); ++r) {
                if (num_traits<T>::to_double(f.fan_out_prob(k, r)) == 0.0) continue;
                acc += num_traits<T>::to_double(f.fan_out_link(k, r)) *
                       num_traits<T>::to_double(f.fan_out_prob(k, r));
                ++cnt;
            }
        if (cnt > 0) nd.tasks_per_link = acc / static_cast<double>(cnt);
    }

    /** True when class r of `fork_node` routes to node k, read off `rtnodes`. */
    bool fork_links_to(std::size_t fork_node, std::size_t k, std::size_t r) const {
        const std::size_t K = sn_.classes.size(), I = sn_.nodes.size();
        if (sn_.rtnodes.rows() < I * K) return false;
        for (std::size_t s = 1; s <= K; ++s)
            if (num_traits<T>::to_double(
                    sn_.rtnodes((fork_node - 1) * K + r - 1, (k - 1) * K + s - 1)) != 0.0)
                return true;
        return false;
    }

    /**
     * Node indexes a fork override applies to. `dest_node` 0 means every
     * destination the fork actually links to, so the routing must be in place
     * by the time this runs -- which is what the recorded-override replay in
     * `link()` guarantees whichever order the caller used.
     */
    std::vector<std::size_t> fork_dests(std::size_t fork_node, std::size_t dest_node) const {
        std::vector<std::size_t> out;
        const std::size_t I = sn_.nodes.size(), K = sn_.classes.size();
        if (dest_node != 0) {
            if (dest_node > I)
                throw InputError("a fork override names a node that does not exist");
            out.push_back(dest_node);
            return out;
        }
        for (std::size_t k = 1; k <= I; ++k)
            for (std::size_t r = 1; r <= K; ++r)
                if (fork_links_to(fork_node, k, r)) { out.push_back(k); break; }
        if (out.empty())
            throw InputError("a fork override was set on '" + sn_.nodes[fork_node - 1].name +
                             "', which links nowhere yet: call link() before the override");
        return out;
    }

    /**
     * The checks a model must pass before its struct means anything.
     *
     * They are the ones whose absence produces a struct that solves to a
     * plausible wrong answer rather than to an error: a class with no service
     * anywhere, an open class with no arrival, a Fork with no Join.
     */
    void validate() const {
        if (sn_.classes.empty()) throw InputError("Network '" + sn_.name + "': it has no classes");
        if (sn_.stations.empty())
            throw InputError("Network '" + sn_.name + "': it has no stations");
        // WHO SWITCHES INTO WHOM. A job can enter a class by arriving in it, or
        // by CLASS SWITCHING into it from another class -- through a routing
        // block with r != s, through a ClassSwitch node's matrix, or through a
        // cache's hit/miss/retrieval switch. `switches` is that edge relation,
        // and the only test below that reads it is the "no service process
        // anywhere" one: a class reached ONLY by switching legitimately has no
        // service of its own at some stations.
        //
        // There is deliberately NO reachability closure over these edges any
        // more. It existed to refuse an open class no job can enter, which is
        // not an error -- see the note on `gallery_erlerl1` below. MATLAB and
        // python have no check of this kind at all, so there was never a
        // reference predicate to copy. `validate()` runs BEFORE
        // `refresh_routing`, so the edges are read from the raw `P` and from
        // `csmatrix`, not from `Peff`, which does not exist yet.
        const std::size_t K = sn_.classes.size();
        const T zero = num_traits<T>::from_int(0);
        std::vector<std::vector<bool>> switches(K, std::vector<bool>(K, false));
        for (std::size_t r = 0; r < K; ++r)
            for (std::size_t sc = 0; sc < K; ++sc) {
                if (r == sc) continue;
                for (std::size_t i = 1; i <= sn_.nodes.size() && !switches[r][sc]; ++i)
                    for (std::size_t j = 1; j <= sn_.nodes.size(); ++j)
                        if (sn_.get_route(r + 1, sc + 1, i, j) > zero) {
                            switches[r][sc] = true;
                            break;
                        }
            }
        for (const auto& kv : sn_.csmatrix) {
            const Matrix<T>& C = kv.second;
            for (std::size_t r = 0; r < K && r < C.rows(); ++r)
                for (std::size_t sc = 0; sc < K && sc < C.cols(); ++sc)
                    if (r != sc && C(r, sc) > zero) switches[r][sc] = true;
        }
        for (const auto& kv : sn_.nodeparam) {
            const CacheParam<T>& cp = kv.second;
            for (std::size_t r = 0; r < K; ++r) {
                if (r < cp.hitclass.size() && cp.hitclass[r] >= 1 && cp.hitclass[r] <= K)
                    switches[r][cp.hitclass[r] - 1] = true;
                if (r < cp.missclass.size() && cp.missclass[r] >= 1 && cp.missclass[r] <= K)
                    switches[r][cp.missclass[r] - 1] = true;
            }
            for (const auto& row : cp.retrieval_classes)
                for (std::size_t r = 0; r < K && r < row.size(); ++r)
                    if (row[r] >= 1 && row[r] <= K) switches[r][row[r] - 1] = true;
        }
        // SPAWN ON COMPLETION REACHES A CLASS TOO. A phase-2 continuation is
        // injected by the completion of its trigger and never arrives at a
        // Source or crosses a switch, so without this edge it reads as a class
        // no job can enter and a well-formed LQN phase-2 model is refused.
        for (std::size_t r = 0; r < K; ++r)
            if (sn_.classes[r].spawn >= 1 && sn_.classes[r].spawn <= K)
                switches[r][sn_.classes[r].spawn - 1] = true;
        for (std::size_t r = 0; r < sn_.classes.size(); ++r) {
            // A class reached only by switching legitimately has no service of
            // its own at some stations; the served test still applies to the
            // rest, so it is kept for every class that is not switched into.
            bool switched_into = false;
            for (std::size_t q = 0; q < K; ++q)
                if (switches[q][r]) switched_into = true;
            bool served = false;
            for (std::size_t i = 0; i < sn_.stations.size(); ++i)
                if (!sn_.service[i][r].disabled) served = true;
            // In an SPN the timing lives in the transition modes, not in station
            // service: a Place is a token container and only a QUEUEING place
            // carries a service process, so a token class served nowhere is a
            // well-formed net rather than an incomplete one.
            bool petri = false;
            for (std::size_t nd = 0; nd < sn_.nodes.size(); ++nd)
                if (sn_.nodes[nd].nodetype == NodeType::Transition) petri = true;
            if (!served && !switched_into && !petri)
                throw InputError("Network '" + sn_.name + "': class '" + sn_.classes[r].name +
                                 "' has no service process at any station");
            // AN UNREACHABLE OPEN CLASS IS WELL FORMED, and this used to refuse
            // it. `gallery_erlerl1` ships in all three reference suites with a
            // second open class whose arrival is Disabled and which nothing
            // routes into; MATLAB and native Python both SOLVE it and simply
            // report no row for that class, because a class no job can enter
            // carries zero of every metric and the table drops an all-zero row.
            // The refusal made this port the only one that could not read its
            // own gallery. Same argument as the fork-with-no-join case below:
            // a construct the reference suites ship and the reference solvers
            // answer is not an input error, whatever it looks like in isolation.
        }
        // A FORK WITH NO JOIN IS WELL FORMED ONLY IF ITS SIBLINGS CAN LEAVE.
        // `fj_nojoin` ships in all three reference suites as an OPEN model
        // whose fork branches each end at the Sink, and MATLAB and native
        // Python both solve it, so a blanket refusal is wrong: the
        // synchronisation point is what a Join provides, and a model that
        // never synchronises simply has none (see fj_driver.h, which drives
        // forkLambda from the fork's own firing rate in that case).
        //
        // A join-less fork whose branches RETURN INTO THE MODEL is a different
        // object. Every firing turns one job into k siblings, none of them ever
        // merges and none ever departs, so the population is not conserved and
        // grows without bound. The reference has no check for it and cannot
        // solve it either: `sortForks` calls `nestedForks(f, [])`, whose
        // `startNode == endNode` test can never hold against an empty join, and
        // on a closed model it recurses until MATLAB reports "Out of memory.
        // The likely cause is an infinite recursion" (measured 2026-07-30).
        // Refusing here names the node the caller has to close.
        for (std::size_t i = 0; i < sn_.nodes.size(); ++i) {
            if (sn_.nodes[i].nodetype != NodeType::Fork) continue;
            bool closed = false;
            for (const auto& fjp : sn_.fj)
                if (fjp.first == i + 1) closed = true;
            if (closed) continue;
            // Can a sibling ever leave? Reachability of a Sink from the Fork
            // over the raw routing graph, any class pair -- `validate()` runs
            // before refresh_routing, so `Peff` does not exist yet.
            std::vector<bool> seen(sn_.nodes.size(), false);
            std::vector<std::size_t> stack(1, i + 1);
            seen[i] = true;
            bool departs = false;
            while (!stack.empty() && !departs) {
                const std::size_t u = stack.back();
                stack.pop_back();
                for (std::size_t v = 1; v <= sn_.nodes.size() && !departs; ++v) {
                    if (seen[v - 1]) continue;
                    bool edge = false;
                    for (std::size_t r = 1; r <= K && !edge; ++r)
                        for (std::size_t s = 1; s <= K; ++s)
                            if (sn_.get_route(r, s, u, v) > zero) {
                                edge = true;
                                break;
                            }
                    if (!edge) continue;
                    if (sn_.nodes[v - 1].nodetype == NodeType::Sink) {
                        departs = true;
                        break;
                    }
                    seen[v - 1] = true;
                    stack.push_back(v);
                }
            }
            if (!departs)
                throw InputError("Network '" + sn_.name + "': the Fork '" + sn_.nodes[i].name +
                                 "' is not closed by a Join and no Sink is reachable from it, "
                                 "so its siblings can neither merge nor depart");
        }
    }
};

}  // namespace qn
}  // namespace line

#endif  // LINE_LANG_QN_NETWORK_BUILDER_H
