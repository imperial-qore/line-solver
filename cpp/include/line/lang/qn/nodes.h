/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_QN_NODES_H
#define LINE_LANG_QN_NODES_H

/**
 * @file
 * @ingroup line_lang
 * @ingroup line_public
 * The model API a user writes, spelled as its Python twin.
 *
 *     Network model("model");
 *     Delay  delay (model, "Delay");
 *     Queue  queue1(model, "Queue1", SchedStrategy::PS);
 *     Source source(model, "Source");
 *     Sink   sink  (model, "Sink");
 *     ClosedClass closed(model, "ClosedClass", 2, delay, 0);
 *     queue1.set_service(closed, Exp(1.0));
 *
 * against Python's
 *
 *     model  = Network('model')
 *     delay  = Delay(model, 'Delay')
 *     queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
 *     closed = ClosedClass(model, 'ClosedClass', 2, delay, 0)
 *     queue1.set_service(closed, Exp(1.0))
 *
 * EVERY CLASS HERE IS A HANDLE, NOT A NODE. It holds the model and the 1-based
 * index `qn::Network` already returns, and every method forwards to the builder
 * call of the same name -- `network_builder.h` remains the engine and is not
 * touched. That is what lets old and new code mix: a handle CONVERTS to its
 * index, so it drops straight into `RoutingMatrix::set` and into any call still
 * written against the index API.
 *
 * The handles are `double`-only, as Python is. The templated
 * `qn::Network<T>` stays reachable for the multiprecision paths.
 */

#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_builder.h"

namespace line {

typedef qn::Network<double> NetworkModel;
typedef qn::RoutingMatrix<double> Routing;
typedef lang::Distrib<double> Dist;

class JobClass;

/** A node of the model: the index it was given, and the model that owns it. */
class Node {
  public:
    /** `get_index()`: the 1-based node index the builder assigned. */
    std::size_t get_index() const { return idx_; }
    /** `getName()`. */
    const std::string& get_name() const { return name_; }
    /** The model this node belongs to. */
    NetworkModel& model() const { return *model_; }

    /**
     * A handle IS its index wherever one is expected.
     *
     * This is what lets `P.set(closed, closed, delay, queue1, 1.0)` take
     * handles, and what lets a half-migrated example keep compiling.
     */
    operator std::size_t() const { return idx_; }

    /** `setRouting(class, strategy)`. */
    void set_routing(std::size_t cls, lang::RoutingStrategy rs) {
        model_->set_routing(idx_, cls, rs);
    }
    /** `setRoutingWeight(class, weights)`: the WRROBIN share per destination. */
    void set_routing_weights(std::size_t cls, const std::map<std::size_t, double>& weights) {
        model_->set_routing_weights(idx_, cls, weights);
    }

  protected:
    Node(NetworkModel& m, std::size_t idx, const std::string& nm)
        : model_(&m), idx_(idx), name_(nm) {}
    NetworkModel* model_;
    std::size_t idx_;
    std::string name_;
};

/** A node that holds jobs and serves them: MATLAB's `Station`. */
class Station : public Node {
  public:
    /** `setService(class, distribution)`. */
    void set_service(std::size_t cls, const Dist& d) { model_->set_service(idx_, cls, d); }
    /** `setNumberOfServers(n)`. */
    void set_number_of_servers(double n) { model_->set_number_of_servers(idx_, n); }
    /** `setCapacity(k)`, the K of Kendall's notation. */
    void set_capacity(double k) { model_->set_capacity(idx_, k); }
    /** `setClassCapacity(class, k)`. */
    void set_class_capacity(std::size_t cls, double k) {
        model_->set_class_capacity(idx_, cls, k);
    }
    /** `setDropRule(class, rule)`. */
    void set_drop_rule(std::size_t cls, lang::DropStrategy rule) {
        model_->set_drop_rule(idx_, cls, rule);
    }
    /** `setSchedParam(class, weight)`: the DPS/GPS share. */
    void set_sched_param(std::size_t cls, double weight) {
        model_->set_sched_param(idx_, cls, weight);
    }
    /** `setLoadDependence(alpha)`. */
    void set_load_dependence(const std::vector<double>& alpha) {
        model_->set_load_dependence(idx_, alpha);
    }
    /** `setPollingType(type, k)`. */
    void set_polling_type(lang::PollingType rule, int par = 0) {
        model_->set_polling_type(idx_, rule, par);
    }
    /** `setSwitchover(class, distribution)`. */
    void set_switchover(std::size_t cls, const Dist& so) {
        model_->set_switchover(idx_, cls, so);
    }
    /** `setSwitchover(fromClass, toClass, distribution)`. */
    void set_switchover(std::size_t from_cls, std::size_t to_cls, const Dist& so) {
        model_->set_switchover(idx_, from_cls, to_cls, so);
    }
    /** `setBreakdown(failure, repair)`. */
    void set_breakdown(const Dist& failure, const Dist& repair) {
        model_->set_breakdown(idx_, failure, repair);
    }
    /** `addServerType(type)`: one pool of a heterogeneous station. */
    void add_server_type(const qn::Station<double>::ServerType& stype) {
        model_->add_server_type(idx_, stype);
    }
    /** `getStationIndex()`: the 1-based station index, distinct from the node one. */
    std::size_t get_station_index() const { return model_->station_index(idx_); }

  protected:
    Station(NetworkModel& m, std::size_t idx, const std::string& nm) : Node(m, idx, nm) {}
};

/** `Queue(model, name, strategy)`. */
class Queue : public Station {
  public:
    Queue(NetworkModel& m, const std::string& nm,
          lang::SchedStrategy sched = lang::SchedStrategy::FCFS)
        : Station(m, m.add_queue(nm, sched), nm) {}
};

/** `Delay(model, name)`: the infinite-server station. */
class Delay : public Station {
  public:
    Delay(NetworkModel& m, const std::string& nm) : Station(m, m.add_delay(nm), nm) {}
};

/** `Source(model, name)`: the external arrival station. */
class Source : public Station {
  public:
    Source(NetworkModel& m, const std::string& nm) : Station(m, m.add_source(nm), nm) {}
    /** `setArrival(class, distribution)`. */
    void set_arrival(std::size_t cls, const Dist& d) { model_->set_arrival(idx_, cls, d); }
    /** `setArrivalBatch(class, batchSizeDist)`. */
    void set_arrival_batch(std::size_t cls, const Dist& d) {
        model_->set_arrival_batch(idx_, cls, d);
    }
    /**
     * `setMarkedArrival(process, classes)`: bind a MARKED arrival process, one
     * class per mark.
     *
     * A mark is a label the process attaches to each event, and here it selects
     * the CLASS of the arriving job: mark k emits jobs of `classes[k-1]`, with
     * all marks driven by ONE shared modulating chain. That sharing is the point
     * and it is why this is not `set_arrival` per class with copies of the
     * process: K copies are K independent streams. The binding IS the order of
     * `classes`.
     *
     * Takes any marked family the engine declares: MMAP and MPH (stationary),
     * MMAPt and MPHt (schedule-bearing), or BMMAPt, whose epochs also release a
     * BATCH of jobs that all carry the mark's class.
     */
    void set_marked_arrival(const Dist& d, const std::vector<std::size_t>& classes) {
        for (std::size_t k = 0; k < classes.size(); ++k)
            model_->set_arrival(idx_, classes[k], d);
        model_->set_marked_classes(idx_, classes);
    }
};

/** `Sink(model, name)`: the external departure node, which holds no jobs. */
class Sink : public Node {
  public:
    Sink(NetworkModel& m, const std::string& nm) : Node(m, m.add_sink(nm), nm) {}
};

/** `Router(model, name)`: a stateless routing node. */
class Router : public Node {
  public:
    Router(NetworkModel& m, const std::string& nm) : Node(m, m.add_router(nm), nm) {}
};

/** `ClassSwitch(model, name, C)`. */
class ClassSwitch : public Node {
  public:
    ClassSwitch(NetworkModel& m, const std::string& nm, const Matrix<double>& C)
        : Node(m, m.add_class_switch(nm), nm) {
        m.set_class_switch_matrix(idx_, C);
    }
    /** `setClassSwitchingMatrix(C)`. */
    void set_class_switching_matrix(const Matrix<double>& C) {
        model_->set_class_switch_matrix(idx_, C);
    }
};

/** `Fork(model, name)`. */
class Fork : public Node {
  public:
    Fork(NetworkModel& m, const std::string& nm, double tasks_per_link = 1.0)
        : Node(m, m.add_fork(nm, tasks_per_link), nm) {}
};

/** `Join(model, name, fork)`. */
class Join : public Node {
  public:
    Join(NetworkModel& m, const std::string& nm, std::size_t fork_node)
        : Node(m, m.add_join(nm, fork_node), nm) {}
    /** `setStrategy(strategy, quorum)`. */
    void set_strategy(lang::JoinStrategy strategy, double quorum = 0.0) {
        model_->set_join_strategy(idx_, strategy, quorum);
    }
};

/** `Cache(model, name, params)`. */
class Cache : public Node {
  public:
    Cache(NetworkModel& m, const std::string& nm, const qn::CacheParam<double>& par)
        : Node(m, m.add_cache(nm, par), nm) {}
    /** `setItemReadClasses(read_classes, hit_classes)`, one entry per item. */
    void set_item_read_classes(const std::vector<std::size_t>& read_classes,
                               const std::vector<std::size_t>& hit_classes) {
        model_->set_item_read_classes(idx_, read_classes, hit_classes);
    }
};

/** `Place(model, name)`: a Petri-net place. */
class Place : public Station {
  public:
    Place(NetworkModel& m, const std::string& nm) : Station(m, m.add_place(nm), nm) {}
    /**
     * `Place(model, name, schedStrategy)`: a place with an EMBEDDED QUEUE.
     *
     * The place becomes a queueing place once `set_service` gives it a law, as
     * `Place.setService` does in MATLAB; until then the discipline is declared
     * and unused and the place behaves as an ordinary one.
     */
    Place(NetworkModel& m, const std::string& nm, lang::SchedStrategy sched)
        : Station(m, m.add_place(nm, sched), nm) {}
    /** `setDepartureDiscipline(class, rule)`: how a served token leaves the depository. */
    void set_departure_discipline(std::size_t cls, lang::DepartureDiscipline rule) {
        model_->set_departure_discipline(idx_, cls, rule);
    }
    /** `setInitialMarking(tokens)`. */
    void set_initial_marking(const std::vector<double>& tokens) {
        model_->set_initial_marking(idx_, tokens);
    }
};

/**
 * `Transition(model, name)`: a Petri-net transition, with its modes declared
 * after it as MATLAB, the JAR and Python declare theirs.
 *
 *     Transition t1(model, "T1");
 *     const std::size_t mode1 = t1.add_mode("Mode1");
 *     t1.set_distribution(mode1, Exp(2.0));
 *     t1.set_enabling_conditions(mode1, c1, p1, 2);
 *     t1.set_firing_outcome(mode1, c1, p2, 2);
 *
 * THE ARCS TAKE HANDLES, not indices: a `Place` where a place is meant and a
 * `JobClass` where a class is meant, so mixing the two is a compile error.
 * Both convert to `std::size_t`, which is exactly what made them
 * interchangeable at a call site written against the index API.
 */
class Transition : public Node {
  public:
    Transition(NetworkModel& m, const std::string& nm) : Node(m, m.add_transition(nm), nm) {}
    /** The all-at-once form, for a caller that already holds a finished block. */
    Transition(NetworkModel& m, const std::string& nm, const qn::TransitionParam<double>& par)
        : Node(m, m.add_transition(nm, par), nm) {}

    /** `addMode(name)`: the new mode's 1-based index. */
    std::size_t add_mode(const std::string& nm) { return model_->add_mode(idx_, nm); }
    /** `setDistribution(mode, dist)`: the mode's firing law. */
    void set_distribution(std::size_t mode, const Dist& d) {
        model_->set_mode_distribution(idx_, mode, d);
    }
    /** `setTimingStrategy(mode, strategy)`: TIMED or IMMEDIATE. */
    void set_timing_strategy(std::size_t mode, lang::TimingStrategy ts) {
        model_->set_mode_timing(idx_, mode, ts);
    }
    /** `setNumberOfServers(mode, n)`; `GlobalConstants::MaxInt` is infinite. */
    void set_number_of_servers(std::size_t mode, double n) {
        model_->set_mode_servers(idx_, mode, n);
    }
    /** `setFiringPriorities(mode, priority)`. */
    void set_firing_priorities(std::size_t mode, double prio) {
        model_->set_firing_priority(idx_, mode, prio);
    }
    /** `setFiringWeights(mode, weight)`: the share among simultaneously enabled modes. */
    void set_firing_weights(std::size_t mode, double w) {
        model_->set_firing_weight(idx_, mode, w);
    }
    /** `setFiringRateDependence(mode, g)`: g(marking) scales the firing rate. */
    void set_firing_rate_dependence(
        std::size_t mode, const std::function<double(const std::vector<double>&)>& g) {
        model_->set_mode_firing_dependence(idx_, mode, g);
    }
    /** `setEnablingConditions(mode, class, place, tokens)`. */
    void set_enabling_conditions(std::size_t mode, const JobClass& cls, const Place& place,
                                 double tokens);
    /** `setInhibitingConditions(mode, class, place, tokens)`; absent means never. */
    void set_inhibiting_conditions(std::size_t mode, const JobClass& cls, const Place& place,
                                   double tokens);
    /** `setFiringOutcome(mode, class, node, tokens)`; the node may be a Sink. */
    void set_firing_outcome(std::size_t mode, const JobClass& cls, const Node& dest,
                            double tokens);
};

// ---------------------------------------------------------------------------
// Job classes
// ---------------------------------------------------------------------------

/** A job class: the index it was given, and the model that owns it. */
class JobClass {
  public:
    /** `get_index()`: the 1-based class index. */
    std::size_t get_index() const { return idx_; }
    /** `getName()`. */
    const std::string& get_name() const { return name_; }
    /** A handle IS its index wherever one is expected. */
    operator std::size_t() const { return idx_; }

  protected:
    JobClass(NetworkModel& m, std::size_t idx, const std::string& nm)
        : model_(&m), idx_(idx), name_(nm) {}
    NetworkModel* model_;
    std::size_t idx_;
    std::string name_;
};

/** `ClosedClass(model, name, njobs, refstat, prio)`. */
class ClosedClass : public JobClass {
  public:
    ClosedClass(NetworkModel& m, const std::string& nm, double njobs, std::size_t refstat_node,
                int prio = 0)
        : JobClass(m, m.add_closed_class(nm, njobs, refstat_node, prio), nm) {}
};

/** `OpenClass(model, name, prio)`. */
class OpenClass : public JobClass {
  public:
    OpenClass(NetworkModel& m, const std::string& nm, int prio = 0)
        : JobClass(m, m.add_open_class(nm, prio), nm) {}
};

/** `SelfLoopingClass(model, name, njobs, refstat, prio)`. */
class SelfLoopingClass : public JobClass {
  public:
    SelfLoopingClass(NetworkModel& m, const std::string& nm, double njobs,
                     std::size_t refstat_node, int prio = 0)
        : JobClass(m, m.add_self_looping_class(nm, njobs, refstat_node, prio), nm) {}
};


// The three arc setters of `Transition`, defined here because they name
// `JobClass`, which is declared below the node handles.

inline void Transition::set_enabling_conditions(std::size_t mode, const JobClass& cls,
                                                const Place& place, double tokens) {
    model_->set_enabling_conditions(idx_, mode, cls.get_index(), place.get_index(), tokens);
}

inline void Transition::set_inhibiting_conditions(std::size_t mode, const JobClass& cls,
                                                  const Place& place, double tokens) {
    model_->set_inhibiting_conditions(idx_, mode, cls.get_index(), place.get_index(), tokens);
}

inline void Transition::set_firing_outcome(std::size_t mode, const JobClass& cls,
                                           const Node& dest, double tokens) {
    model_->set_firing_outcome(idx_, mode, cls.get_index(), dest.get_index(), tokens);
}

// ---------------------------------------------------------------------------
// Routing helpers, MATLAB's static Network methods
// ---------------------------------------------------------------------------

/** `Network.serialRouting(nodes)` for one class pair: 1 -> 2 -> ... -> n. */
inline void serial_routing(Routing& P, std::size_t r, std::size_t s,
                           const std::vector<std::size_t>& nodes) {
    for (std::size_t k = 0; k + 1 < nodes.size(); ++k) P.set(r, s, nodes[k], nodes[k + 1], 1.0);
}

/** `Network.serialRouting(nodes)` on one class of a model. */
inline void serial_routing(Routing& P, std::size_t r, const std::vector<std::size_t>& nodes) {
    serial_routing(P, r, r, nodes);
}

/** The same chain closed into a cycle, which is how a closed model circulates. */
inline void cyclic_routing(Routing& P, std::size_t r, const std::vector<std::size_t>& nodes) {
    serial_routing(P, r, r, nodes);
    if (nodes.size() > 1) P.set(r, r, nodes.back(), nodes.front(), 1.0);
}

}  // namespace line

#endif  // LINE_LANG_QN_NODES_H
