/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_QN_NETWORK_STRUCT_H
#define LINE_LANG_QN_NETWORK_STRUCT_H

/**
 * A queueing network and its refreshed NetworkStruct.
 *
 * SCOPE. This is the `sn` of `matlab/src/lang/@@MNetwork/refreshStruct.m`, held
 * together with the model it was refreshed from, because every consumer of a
 * struct in this port also mutates the model and re-derives it (SolverLN
 * re-parameterises service processes between outer iterations, the fork-join
 * transform rewrites routing). MATLAB does the same through the `sn` its
 * Network caches.
 *
 * It grew out of the SolverLN layer -- which is now `qn::Layer<T>`, a
 * NetworkStruct plus the LQN element annotations -- so the fields a layer never
 * carries are being filled in as the solvers that read them are ported. What is
 * ABSENT is absent by name, never silently: a model needing an unported field
 * is refused where the field would be read.
 *
 * NODES vs STATIONS vs STATEFUL NODES. There are three nested index spaces,
 * exactly as `sn` has, and they stopped coinciding the moment fork/join
 * arrived:
 *
 *   nodes     everything routing passes through: the stations, plus a Fork, its
 *             output Routers, a Join, a Sink
 *   stateful  the nodes that hold jobs: every node except the Fork
 *   stations  the nodes that serve jobs, including a Join (which serves at rate
 *             Inf) and a Source (whose service process is the arrival process)
 *
 * Routing (`P`, and the `rtnodes` derived from it) lives at NODE level; `rt` is
 * its stochastic complement over the stateful nodes, which is what eliminates
 * the Fork; `visits` is indexed by stateful node and `nodevisits` by node.
 *
 * ClassSwitch NODES ARE MATERIALISED, as `@MNetwork/link.m:225-329` does it:
 * `link()` inserts one `CS_<i>_to_<j>` per ordered pair of linked nodes whose
 * routing switches class, folds the switching probability into the first leg
 * and leaves every surviving route SAME-CLASS. The stochastic complement in
 * `rt` removes them again, since they are not stateful, so `rt` is unchanged by
 * their presence; `rtnodes` and `nnodes` carry the extra hop, which is what the
 * node table reports and what the JSIM export needs (`jmt_writer.h` can only
 * emit a ClassSwitch for a node whose type IS ClassSwitch). Until 2026-08-04
 * this port kept the switch on the EDGE and synthesized nothing; the analytical
 * solvers read that correctly through `route_eff`, but the node count was one
 * short per switch and JMT silently simulated the UNSWITCHED model. Verified
 * rather than assumed: the regression compares chains, inchain, refstat,
 * refclass, njobs, nservers, rates, scv and visits against MATLAB dumps.
 *
 * ARITHMETIC. Rates, service times, routing probabilities and visits are T.
 * Populations and server counts are `double`, since they are counts that may
 * be infinite (a delay station, an open class).
 */

#include "line/util/line_console.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "line/api/mc/dtmc_solve.h"
#include "line/api/mc/dtmc_solve_reducible.h"
#include "line/api/mc/dtmc_stochcomp.h"
#include "line/api/mc/stronglyconncomp.h"
#include "line/api/pfqn/pfqn_sdr.h"
#include "line/lang/lang_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace qn {

using lang::CdScaling;
using lang::GdScaling;
using lang::Distrib;
using lang::DropStrategy;
using lang::GlobalConstants;
using lang::JobClassType;
using lang::ReplacementStrategy;
using lang::NodeType;
using lang::ProcessType;
using lang::RoutingStrategy;
using lang::routing_to_text;
using lang::SchedStrategy;

/** One fork firing synchronization: `sn.fjsync{k}`. */
template <class T>
struct FjSync {
    std::size_t fork = 0;   ///< 1-based Fork node
    std::size_t join = 0;   ///< 1-based Join node that closes it
    std::size_t cls = 0;    ///< 1-based ORIGINAL class being forked
    std::size_t tag = 0;    ///< 1-based tag this entry allocates
    std::vector<std::size_t> branchheads;  ///< 1-based node per branch
    std::vector<std::size_t> auxclasses;   ///< the tag's auxiliary class per branch
    /** (B x T) every auxiliary class of this (fork, class), for the tag scan. */
    std::vector<std::vector<std::size_t>> auxall;
    std::size_t weight = 1;  ///< tasksPerLink: siblings emitted per branch
    /**
     * Per-branch tasksPerLink, EMPTY when every branch carries `weight`.
     *
     * A fork may send a different number of tasks down each link
     * (`Fork.setTasksPerLink(class, n, dest)`), and the count is what sizes the
     * auxiliary class capacity and the join's required count, so it cannot be
     * collapsed to the node-wide mean. It stays empty in the uniform case so
     * that a plain fork carries exactly the shape it always did; the two paths
     * emit the same interleaved order there, so this is a shape convention and
     * not a behavioural fork.
     */
    std::vector<std::size_t> weightlink;
    T prob = num_traits<T>::from_int(1);
};

/**
 * `sn.nodeparam{j}.fj`: what a Join node needs to fire on identity.
 *
 * `auxmatrix[r]` is the (B x T) matrix of auxiliary sibling classes minted for
 * original class r -- row b is branch b, column t is tag t -- and `required[r][b]`
 * is how many siblings of branch b a firing consumes (`tasksPerLink` under the
 * only join strategy the exact implementation accepts). It lives here rather
 * than in `fj_tag.h` because `NetworkStruct` carries it and the event layer reads
 * it; `fj_tag.h` is what FILLS it.
 */
struct FjJoinParam {
    std::size_t fork = 0;
    std::vector<std::size_t> origclasses;
    std::map<std::size_t, std::vector<std::vector<std::size_t>>> auxmatrix;
    std::map<std::size_t, std::vector<std::size_t>> required;
};

/**
 * Variable forking levels, the twin of MATLAB `sn.nodeparam{f}.fanOutLink` /
 * `.fanOutProb` / `.fanOutDist`.
 *
 * All three are (nnodes x nclasses) and indexed by DESTINATION NODE rather than
 * by link ordinal, because the link order is an artefact of `connmatrix`
 * traversal and renumbers whenever the model is relinked.
 *
 * `fan_out_link(k,r)`: expected tasks sent to destination k for class r, zero
 * on a link class r does not take. `fan_out_prob(k,r)`: probability the branch
 * fires at all, so the SIBLING COUNT is random even when each link carries a
 * fixed number. `fan_out_dist[k][r]`: the jobs-per-link distribution, with an
 * unset (DISABLED) entry meaning degenerate at `fan_out_link(k,r)`.
 *
 * This is a per-node BLOCK keyed by node index -- the shape `joindecl` already
 * uses, and the shape MATLAB's `nodeparam{f}` has -- rather than three fields on
 * `NodeDef`, because `NodeDef` is not templated on the arithmetic and a
 * `Matrix<T>` cannot live there. A fork with no override has no entry at all,
 * which is how a consumer tells the classic case from the variable one without
 * inspecting a single number.
 */
template <class T>
struct ForkParam {
    Matrix<T> fan_out_link;
    Matrix<T> fan_out_prob;
    std::vector<std::vector<lang::Distrib<T> > > fan_out_dist;
};

/**
 * A node of the network.
 *
 * `station` is the 1-based station index when the node serves jobs and 0 when
 * it does not (a Fork, a Router, a ClassSwitch, a Sink). `stateful` is false
 * for the nodes that hold no jobs and are eliminated by the stochastic
 * complement that builds rt.
 */
struct NodeDef {
    std::string name;
    NodeType nodetype = NodeType::Queue;
    bool stateful = true;
    std::size_t station = 0;
    /**
     * `sn.routing`, per class. PROB means the routing block P carries the
     * probabilities as given; RAND and RROBIN mean the refresh spreads them
     * uniformly over the nodes this one is connected to, which is the routing
     * MATRIX of a round-robin dispatcher (its determinism lives in the higher
     * moments, see refresh_routing). Anything else is state dependent, and the
     * refresh refuses it by name: a state-dependent strategy silently treated
     * as PROB returns a product-form answer for a model that has none.
     */
    std::vector<RoutingStrategy> routing;
    /**
     * The per-destination weights of a WRROBIN dispatcher, per class: a map
     * from 1-based destination NODE index to its weight.
     *
     * MATLAB and the JAR keep these on the node object and the native Python
     * struct calls the field `sn.routingweights`; this port carries them here
     * so a WRROBIN model round-trips through model.json without loss. The
     * refresh still refuses WRROBIN by name, because spreading the weights into
     * a routing MATRIX would answer a state-independent model instead.
     */
    std::vector<std::map<std::size_t, double>> routing_weights;
    /**
     * The scalar parameter of a parameterized dispatcher, per class: the d of
     * a power-of-d (SQ) choice. Zero means the strategy takes none. Carried for
     * the same reason as `routing_weights`.
     */
    std::vector<int> routing_param;
    /**
     * `Fork.output.tasksPerLink` == MATLAB `sn.nodeparam{f}.fanOut`: how many
     * tasks a fork emits per outgoing link. It scales both the MMT auxiliary
     * arrival rate and the synchronisation delay; the default 1 leaves a plain
     * fork unchanged.
     *
     * For a fork whose degree is RANDOM this is the MEAN, so a consumer that
     * only knows this field gets E[tasks per link] rather than a number the
     * fork never emits.
     */
    double tasks_per_link = 1.0;

    /**
     * A Logger node's trace configuration, MATLAB's `Logger` properties and
     * `sn.nodeparam{ind}` for a Logger.
     *
     * The defaults are `Logger.m`'s constructor: timestamp, job id and job
     * class are recorded, the wall-clock start time, the logger name and the
     * two inter-departure columns are not. They are carried on the node rather
     * than derived because a Logger with no file name exports to JMT as a
     * LogTunnel that writes nowhere -- the solver runs, and the trace the user
     * asked for is silently absent.
     */
    struct LoggerParam {
        std::string file_name;   ///< base name, no directory
        std::string file_path;   ///< directory, MATLAB's model.getLogPath
        bool start_time = false;
        bool logger_name = false;
        bool timestamp = true;
        bool job_id = true;
        bool job_class = true;
        bool time_same_class = false;
        bool time_any_class = false;
    };
    LoggerParam logger;
};

/**
 * The parameters of a Cache node, MATLAB's `sn.nodeparam{ind}` for a Cache.
 *
 * `itemcap` is the capacity of each of the `h` cache lists (a single-level
 * cache has one entry). `nitems` is the item population `n`.
 *
 * `pread[v]` is the read distribution of class v over the n items; an EMPTY row
 * is MATLAB's `NaN` placeholder, "class v does not read this cache", and leaves
 * that class's rates at zero. `accost[v][k]` is the ((h+1) x (h+1)) list-to-list
 * routing matrix of class v on item k; an EMPTY `accost` selects the reference
 * default, the linear cache in which an item moves from list l to l+1 on a hit.
 * These two match `da::CacheParam` exactly, so `da_cache_isolate` consumes them
 * without a conversion.
 *
 * `hitclass` and `missclass` are the 1-based classes a job of class r switches
 * into on a hit and on a miss, 0 meaning the transition is disabled.
 */
/**
 * The parameters of a Transition node: MATLAB `sn.nodeparam{ind}` for an SPN
 * transition, as `refreshPetriNetNodes.m` writes them.
 *
 * A transition has MODES, not classes, but its ARCS carry a class: MATLAB's
 * `enablingConditions{m}` is an (nnodes x nclasses) matrix, so one mode can
 * require two tokens of Class1 at a place while another requires one of Class2
 * at the same place. Each mode carries its own enabling and inhibiting
 * conditions over the input places, its firing effect, and its own firing
 * process -- which is why State.fromMarginal treats a Transition's state as
 * per-mode and cannot reuse the per-class station encoding.
 *
 * THE CLASS DIMENSION IS NOT DECORATION. Collapsing the arcs onto the place, as
 * this port did until 2026-08-12, lets a token of one class satisfy another
 * class's pre-arc -- a DIFFERENT net, not an approximation of this one. Every
 * consumer therefore reads a (node, class) pair, and the few that are genuinely
 * class-blind (the MDD level aggregation, the S-invariants) say so by refusing a
 * multiclass net rather than by summing it.
 */
template <class T>
struct TransitionParam {
    std::size_t nmodes = 0;
    std::vector<std::string> modenames;
    /** enabling[m](p,r): class-r tokens of place p (0-based node) mode m needs. */
    std::vector<Matrix<T>> enabling;
    /** inhibiting[m](p,r): class-r tokens of p that BLOCK mode m (Inf = never). */
    std::vector<Matrix<T>> inhibiting;
    /** firing[m](p,r): class-r tokens mode m moves to/from place p when it fires. */
    std::vector<Matrix<T>> firing;

    /**
     * The arcs of one mode summed over classes, for a consumer that is class
     * blind BECAUSE THE NET IS SINGLE CLASS.
     *
     * `is_multiclass()` is what makes that safe, and a caller that cannot honour
     * the class dimension must test it and refuse: the sum of a multiclass net's
     * arcs describes a net whose tokens are interchangeable, which is not this
     * one.
     *
     * A NON-FINITE ENTRY IS SKIPPED, and a row of nothing but non-finite entries
     * stays non-finite: on an enabling arc Inf is JMT's "any number of tokens"
     * and on an inhibiting one it is the ABSENCE of the arc, so adding it in
     * would turn either into a row that no marking satisfies.
     */
    static std::vector<T> arc_total(const std::vector<Matrix<T>>& a, std::size_t m) {
        std::vector<T> v;
        if (m >= a.size()) return v;
        v.assign(a[m].rows(), num_traits<T>::from_int(0));
        for (std::size_t p = 0; p < a[m].rows(); ++p) {
            double s = 0.0;
            bool any = false;
            for (std::size_t r = 0; r < a[m].cols(); ++r) {
                const double x = num_traits<T>::to_double(a[m](p, r));
                if (!std::isfinite(x)) continue;
                s += x;
                any = true;
            }
            v[p] = any || a[m].cols() == 0
                       ? num_traits<T>::from_double(s)
                       : num_traits<T>::from_double(std::numeric_limits<double>::infinity());
        }
        return v;
    }

    /**
     * The inhibiting THRESHOLD of one mode per place, class blind.
     *
     * A threshold is not a count and does not add up: the mode is blocked as
     * soon as ANY class reaches its own bound, so the class-blind reduction is
     * the smallest finite threshold, and Inf where no class declares one.
     */
    static std::vector<T> inhibit_total(const std::vector<Matrix<T>>& a, std::size_t m) {
        std::vector<T> v;
        if (m >= a.size()) return v;
        const double inf = std::numeric_limits<double>::infinity();
        v.assign(a[m].rows(), num_traits<T>::from_double(inf));
        for (std::size_t p = 0; p < a[m].rows(); ++p) {
            double best = inf;
            for (std::size_t r = 0; r < a[m].cols(); ++r) {
                const double x = num_traits<T>::to_double(a[m](p, r));
                if (std::isfinite(x) && x < best) best = x;
            }
            v[p] = num_traits<T>::from_double(best);
        }
        return v;
    }

    /** True when any mode's arcs touch more than one class. */
    bool is_multiclass() const {
        const std::vector<Matrix<T>>* all[3] = {&enabling, &inhibiting, &firing};
        for (int w = 0; w < 3; ++w)
            for (std::size_t m = 0; m < all[w]->size(); ++m) {
                std::size_t touched = 0;
                for (std::size_t r = 0; r < (*all[w])[m].cols(); ++r) {
                    bool any = false;
                    for (std::size_t p = 0; p < (*all[w])[m].rows() && !any; ++p) {
                        const double v = num_traits<T>::to_double((*all[w])[m](p, r));
                        // An inhibiting Inf is the ABSENCE of an arc, not one.
                        any = w == 1 ? (std::isfinite(v) && v > 0.0) : v > 0.0;
                    }
                    if (any) ++touched;
                }
                if (touched > 1) return true;
            }
        return false;
    }

    /** The one class every arc of every mode touches, 1-based; 0 when none does. */
    std::size_t single_class() const {
        const std::vector<Matrix<T>>* all[3] = {&enabling, &inhibiting, &firing};
        for (int w = 0; w < 3; ++w)
            for (std::size_t m = 0; m < all[w]->size(); ++m)
                for (std::size_t r = 0; r < (*all[w])[m].cols(); ++r)
                    for (std::size_t p = 0; p < (*all[w])[m].rows(); ++p) {
                        const double v = num_traits<T>::to_double((*all[w])[m](p, r));
                        if (w == 1 ? (std::isfinite(v) && v > 0.0) : v > 0.0) return r + 1;
                    }
        return 0;
    }
    std::vector<double> nmodeservers;   ///< servers per mode, may be infinite
    std::vector<double> firingprio;     ///< firing priority per mode
    std::vector<T> fireweight;          ///< weight among simultaneously enabled modes
    std::vector<lang::TimingStrategy> timing;  ///< immediate or timed
    std::vector<lang::Distrib<T>> firingproc;  ///< firing distribution per mode
    std::vector<std::size_t> firingphases;     ///< phase count per mode, 0 when non-Markovian
    /**
     * Marking-dependent firing-rate multiplier g_m(marking); an empty entry is
     * the unit multiplier. Transition.setFiringRateDependence.
     */
    std::vector<std::function<T(const std::vector<T>&)>> firingdep;
};

/**
 * The parameters of a retrial station: MATLAB `sn.retrialProc` and friends.
 *
 * A retrial station has NO waiting line. An arrival that finds every server
 * busy joins an ORBIT and re-attempts at the retrial rate, so its state is an
 * (in-service, orbit) split rather than an ordered buffer -- State.fromMarginal
 * enumerates that split, and the freed server is NOT filled on a completion.
 */
template <class T>
struct RetrialParam {
    /** retrial_proc[r] is the class-r retrial process; empty = not a retrial class. */
    std::vector<lang::Distrib<T>> retrial_proc;
    std::vector<T> retrial_rate;        ///< mu_r, the per-class orbit retry rate
    std::vector<int> max_attempts;      ///< 0 = unbounded
};

/**
 * Setup and delay-off of a station that powers down when it falls idle.
 *
 * MATLAB `Queue.setDelayOff(class, setupTime, delayoffTime)`, which
 * refreshLocalVars copies into `sn.nodeparam{node}{class}`. The server starts
 * a delay-off timer when it empties and shuts down when the timer expires; the
 * next arrival to a shut-down server pays the setup time before service.
 *
 * Stored per class because the reference stores it per class, but every
 * consumer reads ONE pair per station -- `solver_mam_basic.m` takes the LAST
 * class's, so `last()` below is what a solver should call rather than picking
 * a class itself.
 */
template <class T>
struct SetupDelayOffParam {
    std::vector<lang::Distrib<T>> setup;     ///< per class, disabled = not declared
    std::vector<lang::Distrib<T>> delayoff;  ///< per class
    /** The pair the solvers use: the last class that declares one. */
    bool last(lang::Distrib<T>& su, lang::Distrib<T>& doff) const {
        for (std::size_t r = setup.size(); r > 0; --r)
            if (!setup[r - 1].disabled) {
                su = setup[r - 1];
                doff = r - 1 < delayoff.size() ? delayoff[r - 1] : lang::Distrib<T>::disabled_dist();
                return true;
            }
        return false;
    }
};

/**
 * Server breakdown and repair of a station whose server fails and is repaired.
 *
 * MATLAB `Queue.setBreakdown(failure, repair, downService)`, which
 * `refreshStruct` spreads over `sn.hasbreakdown` (per NODE),
 * `sn.breakdownMu` / `sn.repairMu` / `sn.breakdownProc` / `sn.repairProc` (per
 * STATION) and `sn.downServiceRates` (per station and class). They are one
 * feature and are kept as one record here, keyed by station, for the same
 * reason `SetupDelayOffParam` is: presence IS the flag.
 *
 * WHAT A BREAKDOWN IS, and what it is not. The server alternates up and down on
 * independent clocks. A job IN SERVICE when the server goes down is NOT evicted
 * and does not restart: it holds its residual work across the outage, so the
 * outage is a pure interruption. `down_service_rates(r)` is a DEGRADED server,
 * not a stopped one -- a positive rate there means class r keeps being served
 * while the server is down, more slowly. Zero means no service at all, which is
 * the ordinary reading of "broken".
 *
 * THE DEGRADED SERVICE MUST BE EXPONENTIAL, and the reference refuses anything
 * else by name: a phase-type degraded service would need its own phase block in
 * the joint chain, which no codebase builds. Only the RATE is stored for that
 * reason -- there is no distribution left to carry.
 */
template <class T>
struct BreakdownParam {
    lang::Distrib<T> failure;   ///< time to failure of an up server
    lang::Distrib<T> repair;    ///< time to repair of a down server
    T failure_rate = num_traits<T>::from_int(0);  ///< `sn.breakdownMu`: 1 / mean failure time
    T repair_rate = num_traits<T>::from_int(0);   ///< `sn.repairMu`: 1 / mean repair time
    /** `sn.downServiceRates(ist, :)`: per class, 0 = no service while down. */
    std::vector<T> down_service_rates;

    T down_rate_of(std::size_t cls_1based) const {
        if (cls_1based == 0 || cls_1based > down_service_rates.size())
            return num_traits<T>::from_int(0);
        return down_service_rates[cls_1based - 1];
    }
};

template <class T>
struct CacheParam {
    std::size_t nitems = 0;
    std::vector<int> itemcap;
    /**
     * Item read by each per-item class of a cache network (MATLAB
     * `Cache.setItemReadClasses`, `sn.nodeparam{i}.classitem`), 1-based, 0 where
     * the class is not one. Stored rather than inferred from a one-hot pread,
     * which is ambiguous against a genuine single-item popularity.
     */
    std::vector<std::size_t> classitem;
    /**
     * Per-item storage cost (size) and per-list cap on the total cost of the
     * resident items (ton21cache Sec. IX). BOTH EMPTY = unconstrained, the
     * classic model. `costcapglobal` records that the caps came from a single
     * cache-wide value.
     */
    std::vector<int> itemsize;
    std::vector<int> costcap;
    bool costcapglobal = false;
    std::vector<std::vector<T>> pread;           ///< (u) x (n), empty row = NaN
    /**
     * The popularity LAW each class declared, beside the pmf it expands to.
     *
     * Every solver reads `pread`, which is why the pmf is what the reader
     * materializes. The law is kept because JMT's Cache section takes a
     * PARAMETRIC popularity -- a Zipf exponent or a uniform range -- and
     * cannot be given a pmf; without this, exporting a Zipf cache would have to
     * either guess the exponent back out of the pmf or drop the popularity.
     * `type` is NONE when the pmf was supplied directly.
     */
    struct Popularity {
        lang::ProcessType type = lang::ProcessType::NONE;
        double s = 0.0;          ///< Zipf exponent
        std::size_t n = 0;       ///< support size of a Zipf or a DiscreteSampler
    };
    std::vector<Popularity> preadkind;  ///< per class, parallel to `pread`
    lang::ReplacementStrategy replacestrat = lang::ReplacementStrategy::RR;
    std::vector<std::size_t> hitclass, missclass;
    std::vector<std::vector<Matrix<T>>> accost;  ///< (u) x (n) of (h+1)x(h+1), or empty
    /**
     * Delayed-hit retrieval system (Cache.setRetrievalSystem). Zero capacity =
     * none. `retrieval_queues[r]` are the 1-based retrieval-station node indices
     * a read of class r (0-based key) circulates on a miss; `retrieval_classes`
     * is (nitems x nclasses), item i of read class r -> the per-item retrieval
     * class (1-based), 0 where none. Built by set_retrieval_system.
     */
    /** q-LRU admission probability; 1 admits every miss (plain LRU). */
    T qlru = num_traits<T>::from_int(1);
    int retrieval_capacity = 0;
    std::map<std::size_t, std::vector<std::size_t>> retrieval_queues;  ///< read class(0-based)->nodes
    std::vector<std::vector<std::size_t>> retrieval_classes;           ///< (nitems x nclasses), 1-based
    /**
     * The DECLARED initial contents of the cache, as the reference dumps the
     * node's state row: the per-class job counts, then the list contents, then
     * the retrieval bitmap. Empty means the cache starts empty, which is what
     * `initDefault` builds; a warm cache is not derivable from anything else.
     */
    std::vector<T> initstate;
    /**
     * Truncation level of block B: how many secondary requests may be merged
     * onto the in-flight fetches of this cache at once. -1 is UNBOUNDED, which
     * is what a sample path needs and what `State.afterEventCache` calls the
     * `isSimulation` branch; a non-negative value is the enumeration bound an
     * exact solver generates its local state space under.
     */
    long max_pending_retrieval = -1;
};

/**
 * Port of `State.cacheRetrievalClassMap`: the canonical order of a cache's
 * retrieval classes, which is the column order of block B.
 *
 * Block B keys the merged requests by RETRIEVAL CLASS rather than by item so
 * that the originating class, hence its hit class, is recoverable when the
 * fetch completes and the merged requests are released as delayed hits.
 */
template <class T>
void cache_retrieval_class_map(const CacheParam<T>& cp, std::vector<std::size_t>& rc_list,
                               std::vector<std::size_t>& rc_items,
                               std::vector<std::size_t>& rc_orig) {
    rc_list.clear();
    rc_items.clear();
    rc_orig.clear();
    for (std::size_t k = 0; k < cp.retrieval_classes.size(); ++k)
        for (std::size_t c = 0; c < cp.retrieval_classes[k].size(); ++c)
            if (cp.retrieval_classes[k][c] != 0) {
                rc_list.push_back(cp.retrieval_classes[k][c]);
                rc_items.push_back(k + 1);
                rc_orig.push_back(c + 1);
            }
    // The reference sorts by the retrieval class index and carries the two
    // parallel arrays along; a plain index sort reproduces that permutation.
    std::vector<std::size_t> ord(rc_list.size());
    for (std::size_t i = 0; i < ord.size(); ++i) ord[i] = i;
    std::stable_sort(ord.begin(), ord.end(),
                     [&](std::size_t a, std::size_t b) { return rc_list[a] < rc_list[b]; });
    std::vector<std::size_t> l(rc_list.size()), it(rc_list.size()), oc(rc_list.size());
    for (std::size_t i = 0; i < ord.size(); ++i) {
        l[i] = rc_list[ord[i]];
        it[i] = rc_items[ord[i]];
        oc[i] = rc_orig[ord[i]];
    }
    rc_list.swap(l);
    rc_items.swap(it);
    rc_orig.swap(oc);
}

/** One station of the network. */
template <class T>
struct Station {
    std::string name;
    NodeType nodetype = NodeType::Queue;
    SchedStrategy sched = SchedStrategy::FCFS;
    double nservers = 1.0;  ///< may be infinite (a Delay, or an inf-scheduled task)
    bool attr_ishost = false;
    std::size_t attr_idx = 0;  ///< LQN element this station stands for

    /**
     * `sn.schedparam`, per class: the DPS / GPS weight, or the SEPT / LEPT rank.
     * Empty means the discipline takes no parameter; the refresh fills it with
     * ones for DPS and GPS, which is MATLAB's default weight.
     */
    std::vector<T> schedparam;
    /**
     * Station capacity in Kendall's K, as `setCapacity` sets it. Infinite means
     * unbounded; `sn.cap` is derived from it and from the class capacities.
     */
    double cap = std::numeric_limits<double>::infinity();
    /** Per-class buffer from `setChainCapacity`; infinite where unset. */
    std::vector<double> classcap;
    /**
     * Node-level immediate feedback, per class; empty when the station sets
     * none. Only a Queue carries it in the reference, which is also the only
     * node the JSON writer emits it for.
     */
    std::vector<bool> immfeed;
    /**
     * Per-class blocking rule as an INT, with 0 meaning "not set".
     *
     * The sentinel is MATLAB's: DropStrategy has no member with value 0, so a
     * zero entry is what an unset rule looks like there, and the refresh
     * derives those from the capacity. It matters that the two are
     * distinguishable, because an EXPLICIT WAITQ for an open class at a finite
     * buffer is rejected while the derived one is not.
     */
    std::vector<int> droprule;
    /**
     * `sn.lldscaling` for this station: the multiplier at population 1, 2, ...
     * Empty when the station is not load dependent.
     */
    std::vector<T> lldscaling;
    /** `sn.cdscaling` for this station: the class-dependence map, empty when unset. */
    CdScaling<T> cdscaling;
    /**
     * `sn.cdscalingpeak` for this station: the DECLARED peak rate scaling per
     * class, empty when the station is not class dependent.
     *
     * It is not derivable from `cdscaling`: finding max_n beta_r(n) would mean
     * sweeping the whole population lattice, and the reference does not. It is
     * what utilization at a class-dependent station is normalized by, so that
     * U = T*S/peak keeps the T*S/c convention of an ordinary multiserver
     * station; without it a beta emulating two servers reports twice the true
     * utilization. MATLAB's `setClassDependence` requires it.
     */
    std::vector<T> cdscalingpeak;
    /**
     * `sn.jdscaling` for this station: MATLAB's `Station.ljdScaling`, the JOINT
     * dependence map eta_i(n), empty when unset.
     *
     * It has the same signature as `cdscaling` and is folded into it
     * multiplicatively wherever a rate is scaled, exactly as
     * `State.afterEventInit` does. It is kept as a SEPARATE field rather than
     * pre-multiplied into `cdscaling` because the two carry different modelling
     * claims: a `cdscaling` beta_r(n) keeps the product form (it is the
     * class-dependent rate lattice `pfqn_cdfun` evaluates), while an eta_i(n)
     * does not, and `pfqn_mvajd` / `pfqn_ncjd` are selected on that distinction.
     */
    CdScaling<T> jdscaling;
    /**
     * `sn.jdscalingpeak` for this station: the declared peak joint-dependent
     * scaling per class. `setJointDependence` makes it mandatory for the same
     * reason `setClassDependence` does -- utilization is reported as T*S/peak.
     */
    std::vector<T> jdscalingpeak;

    /**
     * `sn.nodeparam{ind}.svcRateFun` for a PAS / OI station: the TOTAL service
     * rate as a function of the ordered microstate, a 1-based list of class
     * indices in queue order.
     *
     * A pass-and-swap or order-independent queue is parameterized by mu(c) as a
     * whole; there is no per-class service distribution, and MATLAB's
     * `setServiceRateFunction` rejects one. The refresh still derives a
     * representative per-class rate mu([r]) so the ordinary rate machinery
     * stays consistent, exactly as `Queue.setServiceRateFunction` does.
     */
    std::function<T(const std::vector<std::size_t>&)> svc_rate_fun;

    /**
     * Polling parameters for a POLLING station, MATLAB's `pollingType`,
     * `switchoverTime` and `pollingPar` on the Queue.
     *
     * `polling_type[r]` is the discipline of class r's buffer (the reference
     * assumes it is identical across buffers), `switchover[r]` its switchover
     * distribution, and `polling_par` the K of a K-limited discipline. Empty
     * `polling_type` means the station is not a polling station.
     */
    std::vector<lang::PollingType> polling_type;
    std::vector<Distrib<T>> switchover;
    int polling_par = 0;
    /**
     * `sn.nodeparam{ind}.swapGraph`: which class a departing job promotes the
     * jobs behind it into. All zero is the ORDER-INDEPENDENT case, where no
     * swapping happens and the station is product-form; a nonzero entry makes
     * it a genuine pass-and-swap station, which the OI analyzer refuses.
     */
    Matrix<T> swap_graph;

    // ---- impatience, balking and heterogeneous servers ---------------------
    //
    // These are declared per class on the station, as the reference declares
    // them, and they are all OPTIONAL: an empty vector means the station
    // declares none, which is not the same as declaring a disabled one.
    // `used_lang_features` emits Reneging / Balking / HeteroServers from them,
    // so a solver that does not implement one refuses the model by name
    // instead of solving a station without it.

    /**
     * `Queue.setPatience(class, dist)`: the abandonment timer of a WAITING job,
     * with `impatience[r]` naming which rule it is. A disabled entry is a class
     * that declares none.
     */
    std::vector<Distrib<T>> patience;
    std::vector<lang::ImpatienceType> impatience;
    /**
     * `Queue.setOrbitImpatience(class, dist)`: abandonment from the RETRIAL
     * ORBIT, which is a different population from the waiting line above -- a
     * job that gave up retrying never occupied a buffer slot.
     */
    std::vector<Distrib<T>> orbit_impatience;
    /** `Queue.setBatchRejectProbability`: per-class rejection of a whole batch. */
    std::vector<T> batch_reject;

    /**
     * One balking threshold: with `min_jobs <= n <= max_jobs` at the station,
     * an arriving job of the class refuses to join with `probability`.
     * `max_jobs = -1` is the wire's spelling of an unbounded upper end.
     */
    struct BalkingThreshold {
        double min_jobs = 0.0;
        double max_jobs = -1.0;
        T probability = num_traits<T>::from_int(0);
    };
    /** Per class; `strategy == NONE` is a class that declares no balking. */
    struct BalkingParam {
        lang::BalkingStrategy strategy = lang::BalkingStrategy::NONE;
        std::vector<BalkingThreshold> thresholds;
    };
    std::vector<BalkingParam> balking;

    /**
     * A heterogeneous server pool: `count` servers that serve only
     * `compatible` classes, each with its own service law.
     *
     * The station's own `service` row stays the class-level default and is what
     * every homogeneous consumer reads; `server_types` is the refinement, and
     * `hetero_policy` says how the pools are picked among.
     */
    struct ServerType {
        std::string name;
        double count = 1.0;
        std::vector<bool> compatible;      ///< per class; empty = every class
        std::vector<Distrib<T>> service;   ///< per class
    };
    std::vector<ServerType> server_types;
    lang::HeteroSchedPolicy hetero_policy = lang::HeteroSchedPolicy::ORDER;

    /**
     * `Queue.setServerParallelism(class, n)`: the servers a job seizes for the
     * whole of its service, JMT's job parallelism. Per class, empty or all ones
     * when every job seizes one server.
     */
    std::vector<std::size_t> server_parallelism;

    /**
     * `Source.setArrivalBatch(class, dist)`: the batch-size law released at
     * each arrival epoch. It does NOT space the epochs -- the arrival process
     * in `service` does -- so the two are separate and both are needed.
     */
    std::vector<Distrib<T>> arrival_batch;
    /** `Source.markedClasses`: the 1-based class of each mark of an MMAP arrival. */
    std::vector<std::size_t> marked_classes;
    /** `Place.departureDiscipline`, per class. */
    std::vector<lang::DepartureDiscipline> departure_discipline;
};

/** One job class of the network. */
struct JobClass {
    std::string name;
    JobClassType type = JobClassType::CLOSED;
    double population = 0.0;   ///< infinite for an open class
    std::size_t refstat = 1;   ///< 1-based reference station
    /**
     * Whether passage through the reference station is a COMPLETION.
     *
     * TRUE BY DEFAULT, as in `JobClass.m:34`, `JobClass.java:108` and the native
     * Python `classes.py:37` -- every other codebase constructs a class that
     * completes, and nothing on the model.json wire carries the flag, so a
     * `false` default here made the SAME model mean different things in this port
     * than in the three it is a port of. The visible consequence was that every
     * response-time law refused: `solver_ctmc_cdf_respt` needs one completing
     * class in the tagged chain to have an event to end the passage at, so
     * `getCdfRespT` on an ordinary two-station model reported that no class
     * completes. The paths that need a NON-completing class -- an auxiliary
     * fork-join sibling, an LN pseudo-class -- set it to false explicitly, and
     * did so already.
     */
    bool completes = true;
    bool is_ref_class = false; ///< marks the chain's reference class
    int attr_kind = -1;        ///< LayeredNetworkElement of the element it stands for
    std::size_t attr_idx = 0;  ///< index of that element
    int prio = 0;
    /**
     * Class-level immediate feedback, ORed with the station's own setting into
     * `sn.immfeed`. It is the class-wide spelling of the same property, and the
     * JSON wire carries it as a bare `"immediateFeedback": true` on the class
     * where the node-level form is a per-class map on the node.
     */
    bool immfeed = false;
    /**
     * `sn.classdeadline(r)`: the soft deadline EDD and EDF order by, and the
     * tardiness JMT reports. Infinite where the class declares none, which is
     * `NetworkStruct.m:17`'s "Inf = no deadline" sentinel.
     */
    double deadline = std::numeric_limits<double>::infinity();
    /**
     * `sn.classspawn(r)`: the 1-based class injected at the SAME station on
     * every completion of this class, 0 where none. MATLAB stores -1 for none;
     * the 0 here is this port's usual absent-index sentinel.
     */
    std::size_t spawn = 0;
    /**
     * A `SelfLoopingClass`: a closed class that perpetually cycles at its
     * reference station. It carries no state beyond `ClosedClass`, so it is
     * built as one; the marker exists so the writer does not silently downgrade
     * the wire type to `Closed` and so `getUsedLangFeatures` can name it.
     */
    bool self_looping = false;
};

/**
 * A network plus its refreshed NetworkStruct.
 *
 * The two halves are one object because every caller that mutates the model
 * (service processes, class populations, routing probabilities) re-derives the
 * struct, exactly as MATLAB's Network does through its cached `sn`. The refresh
 * entry points mirror the MATLAB ones and have the same granularity, which
 * matters for cost: `refresh_rates` is called on every outer iteration of
 * SolverLN for every layer, `refresh_chains` only where the routing changed.
 */
template <class T>
class NetworkStruct {
public:
    std::string name;
    /**
     * `Network.setLogPath` / `getLogPath`: the directory every Logger writes
     * into, and the `logPath` attribute of an exported JMT model.
     *
     * Model-level rather than per-Logger because that is where the reference
     * keeps it -- `Logger`'s constructor REFUSES when it is unset -- and
     * because the JSIM header carries one such path for the whole model.
     */
    std::string log_path;
    std::vector<NodeDef> nodes;        ///< every node, in creation order
    std::vector<Station<T>> stations;  ///< stations[k-1] is the k-th station
    std::vector<std::size_t> station_to_node;  ///< (nstations) 1-based node index
    std::vector<std::size_t> stateful_nodes;   ///< 1-based node indices, ascending
    std::vector<JobClass> classes;
    /** service[i][r], 0-based station and class; a disabled entry marks a pair never visited. */
    std::vector<std::vector<Distrib<T>>> service;

    /**
     * Does (station i, class r) have a service law an analyzer may convert?
     * 0-based, and the ONLY correct precondition for `dist_to_map(service[i][r])`.
     *
     * `disabled` alone is not enough, and the gap is not hypothetical. It is the
     * twin of MATLAB's NaN-in-`sn.rates` sentinel, so it answers "does the class
     * visit this station". A JOIN is visited -- `refreshRates` gives it
     * `rates = Inf`, `scv = 0` -- and yet it has NO service law:
     * `refreshProcessRepresentations` hands it a `Coxian(NaN,NaN)`, which every
     * MATLAB analyzer then skips through a SECOND and separate `any(isnan(D0))`
     * guard. Rational has no NaN, so the C++ twin of that Coxian is a DISABLED
     * `Distrib`, and a guard that tests only `disabled[i][r] || rates <= 0` sails
     * past a Join (Inf > 0) straight into `dist_to_map`, which throws
     * "the distribution is disabled" with no station named. `fj_basic_open`
     * failed exactly this way under `-s mam`.
     *
     * So: `disabled[i][r]` is about the VISIT, `service[i][r].disabled` is about
     * the LAW, and a station can be visited without having one.
     */
    bool has_service_law(std::size_t i, std::size_t r) const {
        return !disabled[i][r] && !service[i][r].disabled;
    }
    /** P[(r,s)] is an (nnodes x nnodes) block; absent means all zero. */
    std::map<std::pair<std::size_t, std::size_t>, Matrix<T>> P;
    /**
     * The routing after refresh_routing() has expanded the non-PROB strategies
     * and folded the class switches in. EMPTY when the expansion is the
     * identity, which is the case for every model whose routing is given as
     * probabilities and has no ClassSwitch node -- and then route_eff() reads P
     * directly, so the two representations never drift apart.
     */
    std::map<std::pair<std::size_t, std::size_t>, Matrix<T>> Peff;
    /**
     * Krzesinski (1987) product-form state-dependent routing, in 0-based
     * STATION indices; empty unless a node declares it with
     * `set_state_dep_routing`. Branch index 1 denotes the complement M-V and is
     * unused. `sdr_nodes` is the same structure in 0-based NODE indices.
     *
     * The routing is state dependent yet keeps a product form of its own, so
     * `sn_has_sd_routing` is true for it while the normalizing-constant solver
     * still accepts it. See _kb/16-state-dependent-routing.md
     */
    pfqn::SdrStruct sdr;
    /** Node-indexed twin of `sdr`. */
    pfqn::SdrStruct sdr_nodes;
    /** fj(f,j): the Join node j that closes the Fork node f, 1-based. */
    std::vector<std::pair<std::size_t, std::size_t>> fj;
    /**
     * `sn.isfjaugmented`: this struct came out of `fj_tag`, so its Fork nodes are
     * STATEFUL and its Join nodes carry a per-class sibling count instead of the
     * ordinary buffer/server split.
     *
     * It is a flag and not an inference from `fj` being non-empty, because the
     * un-augmented struct of the SAME model also has `fj` populated: what
     * distinguishes them is the auxiliary class block, and the event layer must
     * not take the Join branch before it exists.
     */
    bool isfjaugmented = false;
    /**
     * `sn.fjclassmap`: the ORIGINAL class of each auxiliary sibling class,
     * 0 for an original class. Empty unless `isfjaugmented`.
     */
    std::vector<std::size_t> fjclassmap;
    /**
     * `sn.nodeparam{j}.fj` for each Join node: the tag matrix and the required
     * sibling multiplicity `after_event_join` fires on. Declared as an opaque
     * map here and defined in `fj_tag.h`, which owns its layout.
     */
    std::map<std::size_t, FjJoinParam> fjjoinparam;

    // ---- the Source/Sink pair, when the model has open classes -------------
    std::size_t sourceIdx = 0;  ///< 1-based station index of the Source, 0 = none
    std::size_t sinkNode = 0;   ///< 1-based NODE index of the Sink, 0 = none (it is not a station)

    /** Cache parameters by 1-based NODE index; only Cache nodes have an entry. */
    std::map<std::size_t, CacheParam<T>> nodeparam;

    /**
     * The DECLARED initial state of a stateful node, by 1-based node index.
     *
     * `initmarking` is a Place's token count per class, which is the initial
     * marking of an SPN and is not derivable from anything else -- an SPN with
     * no tokens anywhere is a dead net, so dropping it changes the answer to
     * "nothing ever fires". `stateprior` and `statespace` are the pair MATLAB
     * writes together (`StatefulNode.statePrior` over `StatefulNode.space`):
     * the prior is a distribution over the ROWS of that space, so neither is
     * meaningful without the other and the reader refuses a lone one.
     */
    std::map<std::size_t, std::vector<T>> initmarking;
    std::map<std::size_t, std::vector<T>> stateprior;
    std::map<std::size_t, Matrix<T>> statespace;

    /**
     * The DECLARED join rule of a Join node, by 1-based node index.
     *
     * STD waits for every sibling; PARTIAL fires on `quorum` of them. The
     * quorum is what `FjJoinParam::required` holds once the fork-join tagging
     * has run, but that runs on the augmented struct: this is the declaration
     * as the model.json carries it, so the two do not replace each other.
     */
    struct JoinDecl {
        lang::JoinStrategy strategy = lang::JoinStrategy::STD;
        double quorum = 0.0;  ///< 0 = every sibling
    };
    std::map<std::size_t, JoinDecl> joindecl;

    /** Variable forking levels, by 1-based Fork node; absent on a plain fork. */
    std::map<std::size_t, ForkParam<T> > forkparam;

    /** The fork's override block, or null when it declares none. */
    const ForkParam<T>* fork_param_of(std::size_t node) const {
        const typename std::map<std::size_t, ForkParam<T> >::const_iterator it =
            forkparam.find(node);
        return it == forkparam.end() ? static_cast<const ForkParam<T>*>(0) : &it->second;
    }

    /**
     * `sn.reward`: the user-declared reward functions, MATLAB's
     * `model.setReward(name, fn)`.
     *
     * `fn` is evaluated on the AGGREGATE state row -- the per-(station, class)
     * job counts in `(ist-1)*K + k` order, which is what
     * `ctmc_state_space_aggr` builds -- and not on the detailed state. That is
     * the reference's `RewardState` contract, and it is what makes a reward
     * portable across disciplines: the caller writes `state[q]` for a queue
     * length without having to know whether the buffer stores class tags or
     * per-class counts.
     */
    struct Reward {
        std::string name;
        std::function<T(const std::vector<T>&)> fn;
        /**
         * The DECLARATIVE form the reward was built from, when it was: the
         * template name (QLen, Util, Blocking), the 1-based node it is declared
         * at, and the 1-based class it covers (0 = every class).
         *
         * A lambda cannot be written back out, so a reward built from one has
         * an empty `kind` and the writer omits it -- which is what
         * `linemodel_save` does, with a warning, rather than emitting a
         * definition that would be wrong on reload.
         */
        std::string kind;
        std::size_t node = 0;
        std::size_t cls = 0;
    };
    std::vector<Reward> reward;

    /**
     * FINITE CAPACITY REGIONS, MATLAB's `refreshRegions` output.
     *
     * A region caps the jobs (and the memory) held ACROSS a set of stations,
     * which no per-station capacity can express: three stations each able to
     * hold 5 jobs but at most 6 between them is a region, not three caps.
     *
     * `cap(i,r)` is the class-r bound at station i, `cap(i,K)` the station's
     * global bound, and -1 means unbounded -- the reference's sentinel, kept
     * rather than translated to infinity because it is compared with `~= -1`
     * to test MEMBERSHIP as well as boundedness.
     *
     * `rule(r)` decides what an arrival that would violate the region does:
     * DROP loses it, WAITQ parks it outside every station's own queue, which is
     * why the CTMC needs a separate waiting-room state per region and not just
     * a filter on the enumerated space.
     */
    struct Region {
        /** The region's declared name, as the wire carries it; a generated one otherwise. */
        std::string name;
        std::vector<std::vector<double>> cap;   ///< (nstations x nclasses+1), -1 = unbounded
        std::vector<double> maxmem;             ///< per member station, -1 = unbounded
        std::vector<bool> members;              ///< membership, independent of the caps
        std::vector<DropStrategy> rule;         ///< per class
        std::vector<T> weight, size;            ///< per class; size is the memory footprint
        Matrix<T> lincon_A;                     ///< optional linear constraint A n <= b
        std::vector<T> lincon_b;
    };
    std::vector<Region> regions;
    /** Transition (SPN) parameters, keyed by 1-based node index. */
    std::map<std::size_t, TransitionParam<T>> transparam;
    /** Retrial parameters, keyed by 1-based STATION index. */
    std::map<std::size_t, RetrialParam<T>> retrialparam;
    /**
     * Setup / delay-off, keyed by 1-based STATION index.
     *
     * Presence IS MATLAB's `sn.hassetup(ist)`: the reference sets that flag
     * from `~isempty(station.setupTime)`, so a station appears here exactly
     * when it is a setup task's server.
     */
    std::map<std::size_t, SetupDelayOffParam<T>> setupparam;
    /**
     * Server breakdown / repair, keyed by 1-based STATION index.
     *
     * Presence IS MATLAB's `sn.hasbreakdown(node)`, which the reference sets
     * from `~isempty(node.breakdownFailure) && ~isempty(node.breakdownRepair)`
     * -- BOTH, since a server that fails and is never repaired is a different
     * model and the reference declines to infer one.
     */
    std::map<std::size_t, BreakdownParam<T>> breakdownparam;

    /** `sn.hasbreakdown(ind)`: does the NODE's server break down? */
    bool has_breakdown_node(std::size_t ind) const {
        if (ind == 0 || ind > nodes.size()) return false;
        const std::size_t ist = nodes[ind - 1].station;
        return ist != 0 && breakdownparam.find(ist) != breakdownparam.end();
    }
    /** Any station at all, the guard every solver gate needs first. */
    bool has_breakdown() const { return !breakdownparam.empty(); }
    /**
     * The class-switch matrix of a ClassSwitch node, by 1-based NODE index.
     *
     * (nclasses x nclasses), row-stochastic. It is applied on the way OUT of
     * the node, exactly as MATLAB's ClassSwitch does, and the node is not
     * stateful, so the stochastic complement folds it into the edges around it.
     */
    std::map<std::size_t, Matrix<T>> csmatrix;

    // ---- refreshed struct -------------------------------------------------
    std::size_t nstations = 0, nclasses = 0, nchains = 0;
    /**
     * `sn.gdscaling`: the network-level globally state-dependent (Whittle)
     * rate scaling phi(n). Its argument is the FULL (nstations x nclasses)
     * population matrix in row-major order, NOT one station's slice, so unlike
     * `Station::jdscaling` it lives on the struct rather than on a station.
     * Empty when the model declares none. See `set_global_dependence`.
     */
    GdScaling<T> gdscaling;
    /**
     * `sn.gdscalingpeak`: the declared (nstations x nclasses) peak of
     * `gdscaling`, row-major, used to report Util = T*S/peak.
     */
    std::vector<T> gdscalingpeak;
    /**
     * `sn.gdscalingcutoff`: the per-slot OPEN-class truncation used when
     * `gdscaling` is materialized onto the JSON wire (closed classes are
     * tabulated up to their own population). Solving ignores it entirely.
     */
    int gdscalingcutoff = 10;
    /**
     * (nstations x nclasses) service rates and SCVs, with a PARALLEL disabled
     * flag instead of MATLAB's NaN sentinel.
     *
     * MATLAB writes NaN into sn.rates for a (station, class) pair the class
     * never visits, and every consumer tests isnan. Rational has no NaN -- it
     * is a field, not a floating-point format -- so the marker has to be
     * carried out of band or the exact instantiation could not represent a
     * disabled pair at all. The flag is the marker; rates and scv hold zero
     * there, and no consumer may read them without consulting `disabled`.
     */
    Matrix<T> rates, scv;
    std::vector<std::vector<bool>> disabled;

    /**
     * `sn.immfeed`: (nstations x nclasses) IMMEDIATE FEEDBACK, the reference's
     * `refreshStruct` field.
     *
     * A job of class r completing at station i is fed straight back into
     * service, HOLDING THE SERVER, instead of being routed out and re-queued.
     * The flag is the OR of the station's own setting (Queue only) and the
     * class-level one, exactly as MATLAB `@@MNetwork/refreshStruct.m` computes
     * it, so either spelling reaches the same matrix.
     *
     * IT IS NOT A `Feature`, deliberately. `feature_name` is byte-for-byte the
     * MATLAB registry field name, and the reference has no such entry: MVA and
     * NC WARN that they approximate it as class-switching with re-queueing,
     * while CTMC and SSA implement the sample-path semantics. Making it a
     * feature would invent a name no other codebase carries and would turn the
     * reference's warning into a refusal.
     */
    std::vector<std::vector<bool>> immfeed;
    std::vector<std::vector<bool>> chains;        ///< (nchains x nclasses)
    std::vector<std::vector<std::size_t>> inchain;///< 1-based class indices per chain
    std::vector<std::size_t> refclass;            ///< (nchains) 1-based class, 0 = none
    std::vector<Matrix<T>> visits;      ///< (nchains) each (nstateful x nclasses)
    std::vector<Matrix<T>> nodevisits;  ///< (nchains) each (nnodes x nclasses)

    /**
     * `sn.cap` and `sn.classcap`: the total and per-class buffers.
     *
     * Derived by refresh_capacity() from the station capacities and the chain
     * populations, exactly as MATLAB's refreshCapacity does, so a station with
     * no explicit capacity still carries the population bound of the chains
     * that reach it.
     */
    std::vector<double> cap;
    std::vector<std::vector<double>> classcap;
    std::vector<std::vector<DropStrategy>> droprule;

    /**
     * `sn.rt` and `sn.rtnodes`: the class-expanded routing.
     *
     * rt is (nstateful*nclasses) square over the STATEFUL nodes -- the
     * stochastic complement that removes Fork, ClassSwitch and Router nodes --
     * and rtnodes is (nnodes*nclasses) square over every node. Row (i-1)*K + r
     * is node i in class r, which is MATLAB's ordering.
     */
    Matrix<T> rt, rtnodes;

    /**
     * `sn.nvars`, (nnodes x 3R+1): the LOCAL VARIABLE columns each node appends
     * to its state, beyond the buffer and the servers. Layout, verbatim from
     * refreshLocalVars:
     *
     *   1 .. R        modulating phase, one per class with a MAP/MMPP2 process
     *   R+1 .. 2R     routing variable, one per class routed round-robin
     *   2R+1          the SHARED node block: cache width, or the BAS blocked
     *                 marker, or the breakdown status, or the polling
     *                 controller. They share one column and are therefore
     *                 mutually exclusive -- refreshLocalVars rejects the
     *                 combinations rather than widening the state.
     *   2R+1+r        the REPLY blocked-server counter for calling class r
     *
     * The trailing sum(nvars(ind,:)) columns are what every state slicer must
     * take clear of before reading the server block.
     */
    std::vector<std::vector<std::size_t>> nvars;

    /**
     * `sn.isbasblocking`, per NODE: true where the node is the BLOCKING
     * (upstream) side of a true-BAS relation and therefore carries the blocked
     * marker in `nvars` column 2R+1.
     *
     * TRUE BAS is not a drop rule, it is a HELD JOB. Under blocking after
     * service the job that finished at this station cannot leave because its
     * destination is full; it stays at the server, occupying it, until room
     * frees. That needs a state bit -- "the front job here is completed and
     * waiting" -- which no queue length can express, and it is what separates
     * BAS from DROP (where the arrival is lost) and from WAITQ (where the
     * arrival waits at the region gate instead).
     *
     * IT IS KEYED ON THE NODE, NOT ON THE DROP RULE, and the distinction is
     * load-bearing. BAS may be DECLARED on either side: upstream, on the station
     * that will hold the job, or on the full destination (the JMT/LDES
     * convention). Both resolve to the same blocking station, because the held
     * job sits upstream either way, so keying enumeration on the station's own
     * `droprule` misses every destination-declared model.
     */
    std::vector<bool> isbasblocking;
    /**
     * `sn.isbasdestination`, (nstations x nclasses): true where a refusal at
     * this station must BLOCK an upstream BAS station rather than drop the job.
     *
     * `arrival_is_lost` needs it because it sees only the station where the
     * refusal happens. Under the upstream declaration form that station carries
     * no BAS rule of its own, so without this mask an open class refused there
     * would be declared lost, the become-blocked edge would never fire, and the
     * blocking station would behave as if its destination were unbounded.
     */
    std::vector<std::vector<bool>> isbasdestination;

    /**
     * The G-network signal declaration, per CLASS. `issignal` is the gate: the
     * remaining vectors are only read where it is true.
     *
     * `signaltarget` is the 1-based class a NEGATIVE signal may remove, or 0
     * for the classic untargeted Gelenbe customer, which is eligible against
     * every non-signal class. `signalremdist` is the batch-size pmf indexed by
     * batch size 0,1,2,...; an empty entry means "remove exactly one".
     */
    /**
     * The pass-and-swap / order-independent parameters of a PAS station, keyed
     * by station index (Dorsman and Gardner 2024, Sect. 2).
     *
     * `svc_rate_fun` is the total service rate mu(c) of an ordered list of
     * class indices; the rate of the token at position p is the INCREMENT
     * mu(c1..cp) - mu(c1..c_{p-1}), which is what makes the station
     * order-independent. `swap_graph` is the (R x R) adjacency saying which
     * class may take another's place. An OI station is the special case of an
     * empty swap graph.
     */
    struct PasParam {
        std::function<T(const std::vector<std::size_t>&)> svc_rate_fun;
        std::vector<std::vector<bool>> swap_graph;
    };
    std::map<std::size_t, PasParam> pasparam;

    /**
     * `sn.replyblock` (nnodes x nclasses) and `sn.syncreply` (nclasses).
     *
     * A SYNCHRONOUS call: a job of a calling class leaves this station for the
     * callee but KEEPS its server, released only when the matching REPLY class
     * arrives back. `replyblock` marks the (node, calling class) pairs that
     * hold such a server, and `syncreply[r]` is the 1-based reply class the
     * calling class r expects (0 = none). A CTMC has no job identity to key on
     * as LDES does, so the state carries COUNTS of held servers.
     */
    std::vector<std::vector<bool>> replyblock;
    std::vector<std::size_t> syncreply;

    /**
     * The polling controller of a POLLING station, keyed by station index.
     *
     * `switchover[r]` is the walk into buffer r, taken from the buffer the
     * server LEAVES to get there. An Immediate switchover is NOT represented as
     * a state: taking its ~1e8 rate literally would make the generator stiff
     * and add a spurious state per buffer, so `pollingNext` folds it into the
     * enclosing transition instead.
     */
    struct PollingParam {
        lang::PollingType ptype = lang::PollingType::EXHAUSTIVE;
        std::size_t pk = 1;  ///< the K of K-LIMITED
        std::vector<lang::Distrib<T>> switchover;
    };
    std::map<std::size_t, PollingParam> pollingparam;

    /**
     * The polling controller of station `ist`, from whichever API declared it.
     *
     * TWO APIS WRITE THE SAME CONTROLLER. `set_polling` fills `pollingparam`;
     * `Queue.setPollingType` / `Queue.setSwitchover` -- the MATLAB-faithful pair
     * the JSON reader also uses -- fill `Station::polling_type`, `switchover`
     * and `polling_par`. Reading only the first left a station built the second
     * way with NO controller at all, which is not a degraded model but an
     * unrepresentable one: the state handlers then index an empty `polled`.
     *
     * A POLLING station that declares neither still HAS a controller, exactly as
     * `State.pollingInfo` builds one: EXHAUSTIVE service, every switchover
     * immediate. That is a discipline, not an absence.
     */
    PollingParam effective_polling(std::size_t ist) const {
        PollingParam pp;
        if (ist == 0 || ist > stations.size()) return pp;
        const typename std::map<std::size_t, PollingParam>::const_iterator it =
            pollingparam.find(ist);
        if (it != pollingparam.end()) return it->second;
        const Station<T>& st = stations[ist - 1];
        if (!st.polling_type.empty()) {
            pp.ptype = st.polling_type[0];
            if (pp.ptype == lang::PollingType::KLIMITED && st.polling_par >= 1)
                pp.pk = static_cast<std::size_t>(st.polling_par);
        }
        pp.switchover = st.switchover;
        return pp;
    }

    std::vector<bool> issignal;
    std::vector<lang::SignalType> signaltype;
    std::vector<std::size_t> signaltarget;
    std::vector<lang::RemovalPolicy> signalrempolicy;
    std::vector<std::vector<T>> signalremdist;

    /** Total local-variable width of node `ind` (1-based). */
    std::size_t nvars_of(std::size_t ind) const {
        // A node index is 1-BASED, so 0 means "no node" and must return 0. The
        // old guard `nvars.size() < ind` is vacuously false at ind == 0 because
        // both sides are unsigned, and the loop below then indexed nvars[-1]:
        // heap corruption, surfacing as a SIGABRT far from here.
        if (ind == 0 || ind > nvars.size()) return 0;
        std::size_t w = 0;
        for (std::size_t j = 0; j < nvars[ind - 1].size(); ++j) w += nvars[ind - 1][j];
        return w;
    }

    std::size_t nof_stations() const { return stations.size(); }
    std::size_t nof_classes() const { return classes.size(); }
    std::size_t nof_nodes() const { return nodes.size(); }
    std::size_t nof_stateful() const { return stateful_nodes.size(); }

    /** `sn.njobs`: the population of each class, infinite for an open one. */
    std::vector<double> njobs() const {
        std::vector<double> v;
        v.reserve(classes.size());
        for (const JobClass& c : classes) v.push_back(c.population);
        return v;
    }

    /** `sn.nclosedjobs`: the total population of the closed classes. */
    double nclosedjobs() const { return total_jobs(); }

    /** `sn.procid(i,r)`: the process type of a (station, class) pair. */
    ProcessType procid(std::size_t ist, std::size_t r) const {
        return service[ist - 1][r - 1].type;
    }

    /**
     * `sn.phases(i,r)`: the order of the process representation.
     *
     * A PLACE HOLDS TOKENS IN ONE PHASE, and that is not cosmetic. MATLAB's
     * `refreshProcessRepresentations` decides the count from what `ph{ist}{r}`
     * turned out to be, and it has TWO ways of having no service:
     *   - `isempty(ph{ist}{r})` -> `phases = 1`, its "fluid fails otherwise" arm
     *   - a process that is all NaN -> `phases = 0`
     * A non-queueing Place carries a `ServiceTunnel`, so its entry is EMPTY and
     * it gets 1; a Join is handed an explicit `Coxian(NaN,NaN)` and gets 0.
     * C++ spells both "no service" as a disabled `Distrib`, whose `phases()` is
     * 0, so the Place silently took the Join's answer.
     *
     * The cost was total: `from_marginal_node_first` can build no row for a node
     * with 0 phases, so `default_init_state` failed on EVERY Place and every
     * closed SPN was refused with "the model's initial marking admits no state"
     * -- spn_inhibiting, spn_basic_closed and spn_fourmodes alike, under a
     * message about Place populations that were in fact correct.
     *
     * A QUEUEING Place is untouched: it has a real server and a real law, so it
     * is not disabled and answers from its representation as before.
     */
    std::size_t phases_of(std::size_t ist, std::size_t r) const {
        const Distrib<T>& d = service[ist - 1][r - 1];
        if (d.disabled && stations[ist - 1].nodetype == NodeType::Place) return 1;
        return d.phases();
    }
    /**
     * `sn.phasessz(i,r) = max(sn.phases(i,r),1)`: THE WIDTH of class r's phase
     * block in a state row, as opposed to `phases_of`, which is its CONTENT.
     *
     * The two differ exactly where a class is disabled at a station, and the
     * reference keeps the column anyway: a disabled process is stored as a
     * `1 x 1 NaN`, never as an empty, so `length(sn.proc{ist}{r}{1})` is 1 and
     * `State.fromMarginal` emits one always-zero column for it. A Source
     * serving one of six classes therefore writes `[Inf 1 0 0 0 0 0]`, and this
     * port used to write `[Inf 1]`.
     *
     * The narrow row is internally consistent, so it is invisible until it
     * CROSSES A BRIDGE: a `model.json` exported by MATLAB carries the wide row,
     * and decoding it against a narrow layout put the arrival one-hot in
     * another class's column. `solver_ssa_serial` then walked into a state with
     * no enabled transition and reported a deadlock rather than a refusal.
     *
     * Use this wherever a WIDTH or an OFFSET is computed -- `row_layout`,
     * `to_marginal`, every `from_marginal*` builder -- and `phases_of` wherever
     * the question is whether the class has a process at all, or how many real
     * phases it has to iterate over.
     */
    std::size_t phasessz_of(std::size_t ist, std::size_t r) const {
        const std::size_t p = phases_of(ist, r);
        return p > 0 ? p : 1;
    }
    bool has_fork() const {
        for (const NodeDef& n : nodes)
            if (n.nodetype == NodeType::Fork) return true;
        return false;
    }
    /**
     * MATLAB's `any(sn.isstatedep(:,3))` NARROWED TO THE ONE STRATEGY THIS PORT
     * EVALUATES PER STATE: Krzesinski's SDR.
     *
     * RROBIN, WRROBIN, JSQ and SQ are state dependent too and `sn_has_sd_routing`
     * reports them, but their per-state tables need auxiliary state this port
     * does not carry (a round-robin pointer) or are refused outright, so a
     * generator asking "must I re-evaluate the routing at every state?" must ask
     * this, not that.
     */
    bool has_sdr_routing() const {
        for (const NodeDef& n : nodes)
            for (RoutingStrategy rs : n.routing)
                if (rs == RoutingStrategy::SDR) return true;
        return false;
    }
    /**
     * ROUND-ROBIN DISPATCH, the state that makes it deterministic.
     *
     * A round-robin dispatcher is not a coin: which link the next job takes is
     * decided by a POINTER the node carries between departures, and without it
     * `refresh_routing`'s uniform expansion answers a random-routing model under
     * a dispatcher's name. The pointer lives in the node's local-variable block,
     * one column per class that routes RROBIN or WRROBIN
     * (`refresh_local_vars` allocates it as `nvars[ind-1][R+r-1]`).
     *
     * RROBIN stores the DESTINATION NODE INDEX in its slot; WRROBIN stores a
     * POSITION in the weighted cycle, because a repeated outlink must advance
     * once per repetition and a destination value could not tell the copies
     * apart. Both are the reference's encodings (`refreshLocalVars.m:318-355`,
     * `afterEventStation.m:632-658`, `afterEventRouter.m:22-52`).
     */
    std::vector<std::size_t> rr_outlinks(std::size_t ind, std::size_t r) const {
        std::vector<std::size_t> out;
        const T zero = num_traits<T>::from_int(0);
        const std::size_t K = classes.size(), I = nodes.size();
        for (std::size_t j = 1; j <= I; ++j) {
            bool linked = false;
            for (std::size_t s = 1; s <= K && !linked; ++s)
                if (get_route(r, s, ind, j) > zero) linked = true;
            if (linked) out.push_back(j);
        }
        return out;
    }

    /** The WRROBIN cycle: each outlink repeated by its weight, weight 0 once. */
    std::vector<std::size_t> rr_weighted_outlinks(std::size_t ind, std::size_t r) const {
        const std::vector<std::size_t> ol = rr_outlinks(ind, r);
        const std::map<std::size_t, double>* w =
            (ind <= nodes.size() && nodes[ind - 1].routing_weights.size() >= r)
                ? &nodes[ind - 1].routing_weights[r - 1]
                : NULL;
        std::vector<std::size_t> cycle;
        for (std::size_t d = 0; d < ol.size(); ++d) {
            long reps = 1;
            if (w != NULL) {
                const std::map<std::size_t, double>::const_iterator it = w->find(ol[d]);
                if (it != w->end() && it->second > 0.0)
                    reps = std::lround(it->second) < 1 ? 1 : std::lround(it->second);
            }
            for (long q = 0; q < reps; ++q) cycle.push_back(ol[d]);
        }
        return cycle;
    }

    /**
     * 1-BASED index of the pointer of (ind, r) INSIDE the node's local-variable
     * block, or 0 when that pair does not dispatch round-robin.
     *
     * `sum(nvars(ind, 1:(R+class)))`, which counts over the phase columns first:
     * the reference indexes `space_var` with exactly that, and reading the
     * pointer as "the r-th trailing column" instead is wrong the moment a class
     * carries a modulating phase or only some classes dispatch.
     */
    std::size_t rr_var_slot(std::size_t ind, std::size_t r) const {
        const std::size_t R = classes.size();
        if (ind == 0 || ind > nodes.size() || r == 0 || r > R) return 0;
        const std::vector<RoutingStrategy>& rt_i = nodes[ind - 1].routing;
        if (rt_i.size() < r) return 0;
        if (rt_i[r - 1] != RoutingStrategy::RROBIN && rt_i[r - 1] != RoutingStrategy::WRROBIN)
            return 0;
        if (nvars.size() < ind) return 0;
        std::size_t slot = 0;
        for (std::size_t j = 0; j < R + r && j < nvars[ind - 1].size(); ++j)
            slot += nvars[ind - 1][j];
        return slot;
    }

    /**
     * The destination node the pointer in VARROW names, or 0 when (ind, r) does
     * not dispatch round-robin or the pointer is out of range.
     *
     * @param varrow the node's local-variable columns, as `after_event` slices them
     */
    std::size_t rr_dest(std::size_t ind, std::size_t r, const std::vector<T>& varrow) const {
        const std::size_t slot = rr_var_slot(ind, r);
        if (slot == 0 || slot > varrow.size()) return 0;
        const long v = std::lround(num_traits<T>::to_double(varrow[slot - 1]));
        if (nodes[ind - 1].routing[r - 1] == RoutingStrategy::RROBIN)
            return v > 0 ? static_cast<std::size_t>(v) : 0;
        const std::vector<std::size_t> cycle = rr_weighted_outlinks(ind, r);
        if (v < 1 || static_cast<std::size_t>(v) > cycle.size()) return 0;
        return cycle[static_cast<std::size_t>(v) - 1];
    }

    /**
     * Advance the pointer of (ind, r) in VARROW by one position, cyclically.
     *
     * A no-op where the pair does not dispatch round-robin, so a caller can run
     * it unconditionally on every departure.
     */
    void rr_advance(std::size_t ind, std::size_t r, std::vector<T>& varrow) const {
        const std::size_t slot = rr_var_slot(ind, r);
        if (slot == 0 || slot > varrow.size()) return;
        if (nodes[ind - 1].routing[r - 1] == RoutingStrategy::RROBIN) {
            const std::vector<std::size_t> ol = rr_outlinks(ind, r);
            if (ol.empty()) return;
            const long cur = std::lround(num_traits<T>::to_double(varrow[slot - 1]));
            std::size_t idx = ol.size();  // "not found" -> restart at the first
            for (std::size_t d = 0; d < ol.size(); ++d)
                if (static_cast<long>(ol[d]) == cur) { idx = d; break; }
            const std::size_t nxt = (idx + 1 < ol.size()) ? idx + 1 : 0;
            varrow[slot - 1] = num_traits<T>::from_int(static_cast<long>(ol[nxt]));
            return;
        }
        const std::vector<std::size_t> cycle = rr_weighted_outlinks(ind, r);
        if (cycle.empty()) return;
        const long pos = std::lround(num_traits<T>::to_double(varrow[slot - 1]));
        const long nxt = (pos < 1 || static_cast<std::size_t>(pos) >= cycle.size()) ? 1 : pos + 1;
        varrow[slot - 1] = num_traits<T>::from_int(nxt);
    }

    /** Whether ANY (node, class) pair dispatches round-robin. */
    bool has_rr_routing() const {
        for (const NodeDef& n : nodes)
            for (RoutingStrategy rs : n.routing)
                if (rs == RoutingStrategy::RROBIN || rs == RoutingStrategy::WRROBIN) return true;
        return false;
    }

    /** `any(sn.immfeed(:))`: whether any (station, class) pair feeds back. */
    bool has_immediate_feedback() const {
        for (const std::vector<bool>& row : immfeed)
            for (bool b : row)
                if (b) return true;
        return false;
    }
    /** 1-based node index of a station, and the reverse; 0 when absent. */
    /** 1-based node index of station `st`, 0 when the map has no entry. */
    std::size_t node_of_station(std::size_t st) const {
        // A Layer built by the LQN path can carry fewer map entries than
        // stations; returning 0 says "no node" rather than reading past the end.
        if (st == 0 || st > station_to_node.size()) return 0;
        return station_to_node[st - 1];
    }
    /** 1-based stateful index of node `ind`, 0 when the node is not stateful. */
    std::size_t stateful_index(std::size_t ind) const {
        for (std::size_t k = 0; k < stateful_nodes.size(); ++k)
            if (stateful_nodes[k] == ind) return k + 1;
        return 0;
    }

    std::size_t stateful_of_station(std::size_t st) const {
        const std::size_t nd = station_to_node[st - 1];
        for (std::size_t k = 0; k < stateful_nodes.size(); ++k)
            if (stateful_nodes[k] == nd) return k + 1;
        throw InputError("network: a station is not a stateful node");
    }

    /** Add a station, which is also a node, and grow the service table. */
    std::size_t add_station(const Station<T>& st) {
        stations.push_back(st);
        service.emplace_back(classes.size(), Distrib<T>::disabled_dist());
        NodeDef nd;
        nd.name = st.name;
        nd.nodetype = st.nodetype;
        nd.stateful = true;
        nd.station = stations.size();
        nodes.push_back(nd);
        station_to_node.push_back(nodes.size());
        stateful_nodes.push_back(nodes.size());
        grow_routing();
        return stations.size();
    }

    /** Add a non-station node (a Fork, a Router). Returns its 1-based index. */
    std::size_t add_node(const std::string& nm, NodeType ty, bool stateful) {
        NodeDef nd;
        nd.name = nm;
        nd.nodetype = ty;
        nd.stateful = stateful;
        nd.station = 0;
        nodes.push_back(nd);
        if (stateful) stateful_nodes.push_back(nodes.size());
        grow_routing();
        return nodes.size();
    }

    /** Add a class and grow the service table. */
    std::size_t add_class(const JobClass& cl) {
        classes.push_back(cl);
        for (auto& row : service) row.emplace_back(Distrib<T>::disabled_dist());
        return classes.size();
    }

    void set_service(std::size_t station, std::size_t cls, const Distrib<T>& d) {
        service[station - 1][cls - 1] = d;
    }

    /** P{r,s}(i,j) = p, with 1-based NODE and class indices. */
    void set_route(std::size_t r, std::size_t s, std::size_t i, std::size_t j, const T& p) {
        auto key = std::make_pair(r, s);
        auto it = P.find(key);
        if (it == P.end())
            it = P.emplace(key, Matrix<T>(nodes.size(), nodes.size(), num_traits<T>::from_int(0)))
                     .first;
        if (it->second.rows() != nodes.size()) grow_block(it->second);
        it->second(i - 1, j - 1) = p;
    }

    /**
     * Write into the routing the consumers actually read.
     *
     * The Sink -> Source closure is derived by the refresh, not given by the
     * user, so it belongs in the effective routing; writing it into P alone
     * would make it invisible on any model whose routing was expanded.
     */
    void set_route_effective(std::size_t r, std::size_t s, std::size_t i, std::size_t j,
                             const T& p) {
        if (Peff.empty()) {
            set_route(r, s, i, j, p);
            return;
        }
        auto key = std::make_pair(r, s);
        auto it = Peff.find(key);
        if (it == Peff.end())
            it = Peff.emplace(key, Matrix<T>(nodes.size(), nodes.size(), num_traits<T>::from_int(0)))
                     .first;
        if (it->second.rows() != nodes.size()) {
            Matrix<T> g(nodes.size(), nodes.size(), num_traits<T>::from_int(0));
            for (std::size_t a = 0; a < it->second.rows(); ++a)
                for (std::size_t b = 0; b < it->second.cols(); ++b) g(a, b) = it->second(a, b);
            it->second = g;
        }
        it->second(i - 1, j - 1) = p;
    }

    /**
     * P{r,s}(i,j), AS THE USER SET IT.
     *
     * Every consumer of the routing reads route_eff() instead, which is this
     * matrix once refresh_routing() has expanded the strategies that are not
     * literal probabilities. The two agree on a model whose routing is entirely
     * PROB, which is every model the LN layer builder produces.
     */
    T get_route(std::size_t r, std::size_t s, std::size_t i, std::size_t j) const {
        auto it = P.find(std::make_pair(r, s));
        if (it == P.end() || it->second.rows() < nodes.size()) return num_traits<T>::from_int(0);
        return it->second(i - 1, j - 1);
    }

    /** The routing actually in force: the expansion when there is one, else P. */
    T route_eff(std::size_t r, std::size_t s, std::size_t i, std::size_t j) const {
        if (Peff.empty()) return get_route(r, s, i, j);
        auto it = Peff.find(std::make_pair(r, s));
        if (it == Peff.end() || it->second.rows() < nodes.size()) return num_traits<T>::from_int(0);
        return it->second(i - 1, j - 1);
    }

    // -----------------------------------------------------------------------
    // Refresh
    // -----------------------------------------------------------------------

    /**
     * The whole chain, in MATLAB's refreshStruct order.
     *
     * The order is load bearing: the routing expansion must precede the chains
     * (they are read off the class-switch structure of the routing), the chains
     * must precede the capacities (a station's buffer is bounded by the
     * population of the chains that reach it), and the visits must precede
     * nothing but must follow the sink closure, which needs the chains.
     */
    void refresh_struct() {
        line::util::LineConsole::compiling(this->name);
        line::util::LineConsole::compile_detail("refreshing service and arrival processes");
        refresh_rates();
        refresh_sched_param();
        // BEFORE the routing table: a round-robin Router must be STATEFUL, or
        // the stochastic complement that builds `rt` erases it and the pointer
        // has nowhere to live.
        refresh_router_stateful();
        line::util::LineConsole::compile_detail("computing the routing table");
        refresh_routing();
        line::util::LineConsole::compile_detail("computing the chains and the visit ratios");
        refresh_chains();
        line::util::LineConsole::compile_detail(
            "found %s over %s",
            line::util::LineConsole::plural(static_cast<long>(this->nchains), "chain", "chains").c_str(),
            line::util::LineConsole::plural(static_cast<long>(this->classes.size()), "class", "classes").c_str());
        refresh_capacity();
        refresh_rt();
        line::util::LineConsole::compile_detail("refreshing node parameters and state-dependent routing");
        refresh_local_vars();
        refresh_immfeed();
        // Needs the visit ratios, so it runs after refresh_chains.
        check_service_reachable();
    }

    /**
     * True when node `ind` (1-based) holds a server across a synchronous call
     * whose reply class is `r` (1-based), i.e. `r` returns here to release a
     * server rather than to be served by one.
     */
    bool holds_reply_for(std::size_t ind, std::size_t r) const {
        if (ind == 0 || ind > replyblock.size()) return false;
        for (std::size_t k = 0; k < syncreply.size(); ++k)
            if (syncreply[k] == r && k < replyblock[ind - 1].size() && replyblock[ind - 1][k])
                return true;
        return false;
    }

    /**
     * Whether a job of class `r` (0-based) can LEAVE node `ind` (1-based) again.
     *
     * True for anything that is not a service station, and for a station that
     * serves `r`, declares heterogeneous server types (its per-class process is
     * disabled by construction there) or holds a server across a synchronous
     * call whose reply class is `r`. False only for the flow sink itself, which
     * is what stops the walk in `reached_node_classes` -- the same three
     * exemptions the guard applies, kept in one place so the walk and the
     * verdict cannot drift apart.
     */
    bool serves_class(std::size_t ind, std::size_t r) const {
        if (ind == 0 || ind > nodes.size()) return true;
        const NodeType nt = nodes[ind - 1].nodetype;
        if (nt != NodeType::Queue && nt != NodeType::Delay) return true;
        const std::size_t sti = nodes[ind - 1].station;
        if (sti == 0 || sti > stations.size()) return true;
        if (!stations[sti - 1].server_types.empty()) return true;
        if (procid(sti, r + 1) != ProcessType::DISABLED) return true;
        return holds_reply_for(ind, r + 1);
    }

    /**
     * The (node, class) pairs a job can actually ARRIVE at, 0-based on both axes.
     *
     * A forward walk of `rtnodes` from the feed points. See
     * `check_service_reachable` for why the evidence is the UNMASKED routing
     * kernel rather than `nodevisits`, and for the three rules -- absorbing
     * Sink, Source expanded only as a seed, unservable pair reached but not
     * expanded -- that keep the walk from over-approximating.
     *
     * Returns an empty vector when `rtnodes` is not the expected (N*K) square,
     * which leaves the caller checking nothing, exactly as before.
     */
    std::vector<std::vector<bool>> reached_node_classes() const {
        const std::size_t N = nodes.size(), K = nclasses;
        std::vector<std::vector<bool>> reached;
        if (N == 0 || K == 0) return reached;
        if (rtnodes.rows() < N * K || rtnodes.cols() < N * K) return reached;
        reached.assign(N, std::vector<bool>(K, false));
        std::vector<std::vector<bool>> seed(N, std::vector<bool>(K, false));
        for (std::size_t i = 0; i < stations.size(); ++i) {
            if (stations[i].nodetype != NodeType::Source) continue;
            const std::size_t ind = node_of_station(i + 1);
            if (ind == 0 || ind > N) continue;
            for (std::size_t r = 0; r < K; ++r)
                if (procid(i + 1, r + 1) != ProcessType::DISABLED) seed[ind - 1][r] = true;
        }
        for (std::size_t r = 0; r < K && r < classes.size(); ++r) {
            const double pop = classes[r].population;
            if (!std::isfinite(pop) || pop <= 0.0) continue;
            const std::size_t ind = node_of_station(classes[r].refstat);
            if (ind == 0 || ind > N) continue;
            seed[ind - 1][r] = true;
        }
        std::vector<std::pair<std::size_t, std::size_t>> stack;
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t r = 0; r < K; ++r)
                if (seed[i][r]) {
                    reached[i][r] = true;
                    stack.push_back(std::make_pair(i, r));
                }
        while (!stack.empty()) {
            const std::size_t i = stack.back().first;
            const std::size_t r = stack.back().second;
            stack.pop_back();
            const NodeType nt = nodes[i].nodetype;
            if (nt == NodeType::Sink) continue;
            if (nt == NodeType::Source && !seed[i][r]) continue;
            if (!serves_class(i + 1, r)) continue;
            const std::size_t row = i * K + r;
            for (std::size_t col = 0; col < N * K; ++col) {
                if (num_traits<T>::to_double(rtnodes(row, col)) <=
                    lang::GlobalConstants::Zero)
                    continue;
                const std::size_t j = col / K, sIdx = col % K;
                if (!reached[j][sIdx]) {
                    reached[j][sIdx] = true;
                    stack.push_back(std::make_pair(j, sIdx));
                }
            }
        }
        return reached;
    }

    /**
     * Refuses a class that is ROUTED TO a station which cannot serve it.
     *
     * The reference's `sanitize` disables the OUTGOING routing of a class a
     * station cannot serve, which is what keeps it out of that station's visit
     * ratios -- but nothing stopped the class being routed IN, and a class that
     * arrives where it cannot be served is a flow sink: it enters and never
     * leaves. The station-level guard next to it cannot see this, because it
     * asks whether the station serves ANY class, not whether it serves the
     * classes that reach it.
     *
     * One such model gave three different wrong answers, none flagged, on a
     * closed cycle D <-> Q whose class C2 has no service at Q: MVA reported
     * Q/C2 with ArvR 1 against Tput 0, CTMC dropped C2 entirely, and SSA
     * returned D/C2 QLen 2e-06 with the Q rows absent.
     *
     * READS `rtnodes`, WALKED FORWARD FROM THE FEED POINTS -- not `nodevisits`,
     * which this guard read until 2026-09-02 and which 94d5570f3 had made blind
     * to the very case it exists for. That commit extended the `served` mask
     * from the station chain to the NODE chain, and it had to: on a materialised
     * LQN replica the unserved states close into a spurious cycle. But the mask
     * zeroes exactly the (station, class) cell a flow sink shows up in. On
     * Source -> Q -> Sink with class B unservable at Q, B's chain went from
     * Q = 1 to Q = 0 and the guard fell silent, while the Sink still read 1 --
     * flow arriving downstream of a node it never visited. A MASKED VISIT VECTOR
     * CANNOT ANSWER THIS QUESTION, because the mask IS the answer being looked
     * for. Do not route this guard back through `nodevisits` or `visits`; both
     * carry that mask. See _kb/07-cross-language-parity.md.
     *
     * `rtnodes` on its own over-approximates -- it says where a class WOULD go
     * if one existed -- and the WALK is what removes the slack. It starts only
     * at (Source, class) pairs whose arrival process is not DISABLED, and at the
     * reference station of each closed class with a positive population, so the
     * disabled-arrival row a class-switching Source carries is never entered.
     * Three rules keep it honest: a Sink is ABSORBING (`rtnodes` wraps it back
     * to the Source to close the kernel, and following that wrap re-enters every
     * Source row, including the disabled ones the seeding just excluded); a
     * Source is expanded ONLY AS A SEED, for the same reason; and an unservable
     * (station, class) is REACHED BUT NOT EXPANDED, since nothing leaves it --
     * that is the whole complaint -- so nothing downstream of it is evidence of
     * anything.
     *
     * This SUBSUMES the fed-chain precondition the guard used to carry
     * separately: a chain no job can enter has no seed, so its rows are never
     * walked at all. That is strictly finer than the per-chain test it replaces,
     * which admitted every class of a chain any one of whose classes was fed.
     */
    void check_service_reachable() const {
        if (nclasses == 0) return;
        // Same exemption as the reference's sanitize checks: a station of a
        // cache, Petri-net or fork-join model legitimately carries no per-class
        // service.
        for (std::size_t ind = 0; ind < nodes.size(); ++ind) {
            const NodeType nt = nodes[ind].nodetype;
            if (nt == NodeType::Cache || nt == NodeType::Place ||
                nt == NodeType::Transition || nt == NodeType::Fork ||
                nt == NodeType::Join)
                return;
        }
        const std::size_t K = nclasses;
        const std::vector<std::vector<bool>> reached = reached_node_classes();
        if (reached.empty()) return;
        for (std::size_t i = 0; i < nstations; ++i) {
            const NodeType nt = stations[i].nodetype;
            if (nt != NodeType::Queue && nt != NodeType::Delay) continue;
            // A heterogeneous pool carries its service on the server types, so
            // the per-class process is legitimately disabled there. Skipped in
            // all four codebases.
            if (!stations[i].server_types.empty()) continue;
            const std::size_t ind = node_of_station(i + 1);
            if (ind == 0) continue;
            for (std::size_t r = 0; r < K; ++r) {
                if (procid(i + 1, r + 1) != ProcessType::DISABLED) continue;
                // A SYNCHRONOUS REPLY is not served by the station it returns
                // to: it releases the server that station held across the call,
                // which is the whole content of set_sync_reply. Its Disabled
                // service there is the marker of the feature, not a flow sink.
                if (holds_reply_for(ind, r + 1)) continue;
                if (ind - 1 >= reached.size() || r >= reached[ind - 1].size()) continue;
                if (!reached[ind - 1][r]) continue;
                const std::string kind = (nt == NodeType::Delay) ? "Delay" : "Queue";
                throw InputError(kind + " '" + stations[i].name +
                                 "' has no service configured for job class '" +
                                 classes[r].name +
                                 "', but the class is routed to it. Jobs would arrive and "
                                 "never leave. Configure a service for that class, or route "
                                 "it elsewhere.");
            }
        }
    }

    /**
     * Port of the `sn.immfeed` block of `@@MNetwork/refreshStruct.m`.
     *
     * The station's own per-class setting OR the class-level one, which is the
     * reference's `stationHas || classHas`. Sized here rather than in
     * `refresh_rates` because it is a property of the model's topology and not
     * of the service processes, so it must not be cleared when only the rates
     * are re-derived (SolverLN calls `refresh_rates` per layer, per iteration).
     */
    void refresh_immfeed() {
        const std::size_t M = stations.size(), R = classes.size();
        immfeed.assign(M, std::vector<bool>(R, false));
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                const bool station_has =
                    r < stations[i].immfeed.size() && stations[i].immfeed[r];
                immfeed[i][r] = station_has || classes[r].immfeed;
            }
    }

    /**
     * Port of MNetwork.refreshLocalVars: the per-node local-variable widths.
     *
     * A column is reserved only where the model actually needs one, because the
     * width is part of the state encoding -- a spurious column widens every
     * state row and makes it unmatchable against the enumerated space.
     */
    /**
     * A Router that DISPATCHES ROUND-ROBIN holds state, so it must survive the
     * stochastic complement that removes the stateless nodes.
     *
     * The reference makes EVERY Router stateful (`refreshStruct.m:155`); this
     * port promotes only the dispatching ones, because a PROB or RAND Router
     * genuinely holds nothing and complementing it away is both correct and
     * cheaper -- it keeps the state space of every existing model unchanged.
     * The promotion happens at refresh rather than at `add_router`, since
     * `set_routing` is called after the node exists.
     */
    void refresh_router_stateful() {
        for (std::size_t ind = 1; ind <= nodes.size(); ++ind) {
            NodeDef& nd = nodes[ind - 1];
            if (nd.stateful || nd.nodetype != NodeType::Router) continue;
            bool rr = false;
            for (std::size_t r = 0; r < nd.routing.size(); ++r)
                if (nd.routing[r] == RoutingStrategy::RROBIN ||
                    nd.routing[r] == RoutingStrategy::WRROBIN)
                    rr = true;
            if (!rr) continue;
            nd.stateful = true;
            stateful_nodes.insert(
                std::lower_bound(stateful_nodes.begin(), stateful_nodes.end(), ind), ind);
        }
    }

    void refresh_local_vars() {
        const std::size_t R = classes.size();
        refresh_replyblock();
        nvars.assign(nodes.size(), std::vector<std::size_t>(3 * R + 1, 0));
        for (std::size_t ind = 1; ind <= nodes.size(); ++ind) {
            const std::size_t ist = nodes[ind - 1].station;
            // A Markov-modulated process restarts from the phase it was left
            // in, so that phase has to survive in the state between services.
            if (ist != 0)
                for (std::size_t r = 1; r <= R; ++r) {
                    const ProcessType pt = procid(ist, r);
                    if (pt == ProcessType::MAP || pt == ProcessType::MMPP2)
                        nvars[ind - 1][r - 1] += 1;
                }
            // Round-robin routing is stateful: the pointer to the next outgoing
            // link is what makes it deterministic rather than random.
            const std::vector<RoutingStrategy>& rt_i = nodes[ind - 1].routing;
            for (std::size_t r = 1; r <= R && r <= rt_i.size(); ++r)
                if (rt_i[r - 1] == RoutingStrategy::RROBIN ||
                    rt_i[r - 1] == RoutingStrategy::WRROBIN)
                    nvars[ind - 1][R + r - 1] += 1;
            // The shared node block. A Cache stores its list contents here;
            // the BAS marker, the breakdown status and the polling controller
            // claim the same column, which is why the reference rejects those
            // combinations instead of widening the state.
            // The polling controller shares the node block with the cache,
            // the BAS marker and the breakdown status, which is why those
            // combinations are rejected rather than the state widened.
            // The SAME resolution polling_info reads, or its offset misses nvars.
            if (ist != 0 && stations[ist - 1].sched == SchedStrategy::POLLING) {
                const PollingParam pp = effective_polling(ist);
                bool anysw = false;
                for (std::size_t r = 0; r < pp.switchover.size(); ++r) {
                    const lang::Distrib<T>& d = pp.switchover[r];
                    if (!d.disabled && d.D0.rows() > 0 && d.type != ProcessType::IMMEDIATE)
                        anysw = true;
                }
                std::size_t w = anysw ? 2 : 0;   // pos and swk
                if (pp.ptype != lang::PollingType::EXHAUSTIVE) w += 1;  // ctr
                nvars[ind - 1][2 * R] = w;
            }
            // The REPLY block: one counter column per (node, calling class).
            for (std::size_t r = 1; r <= R; ++r)
                if (replyblock.size() >= ind && replyblock[ind - 1].size() >= r &&
                    replyblock[ind - 1][r - 1])
                    nvars[ind - 1][2 * R + r] = 1;
            const typename std::map<std::size_t, CacheParam<T>>::const_iterator ci =
                nodeparam.find(ind);
            if (ci != nodeparam.end()) {
                // Contents, plus block A (a per-item occupancy bitmap) and block
                // B (a per-retrieval-class count of merged secondary requests)
                // when a delayed-hit retrieval system is attached.
                std::size_t w = 0;
                for (std::size_t u = 0; u < ci->second.itemcap.size(); ++u)
                    if (ci->second.itemcap[u] > 0)
                        w += static_cast<std::size_t>(ci->second.itemcap[u]);
                if (ci->second.retrieval_capacity > 0) {
                    w += ci->second.nitems;
                    std::vector<std::size_t> rcl, rci, rco;
                    cache_retrieval_class_map(ci->second, rcl, rci, rco);
                    w += rcl.size();
                }
                nvars[ind - 1][2 * R] = w;
            }
        }
        refresh_bas_blocking();
    }

    /**
     * `sn.replyblock`, DERIVED from `sn.syncreply` and the routing.
     *
     * Port of `refreshLocalVars.m:340-386`. A server is held wherever the REPLY
     * class can arrive: every node that is a station, is not the Source, and is
     * not an infinite server -- an INF station has a server per job, so holding
     * one is immaterial and needs no state. Every such station must be FCFS,
     * because a held server is encoded as a per-class COUNTER, which is exact
     * only where the servers are interchangeable.
     *
     * Derived rather than declared because the model layer carries only the
     * class-to-class binding (`JobClass.setReplySignalClass`), and a caller
     * naming the nodes itself would silently disagree with the reference on any
     * model whose routing sends the reply somewhere it did not think of.
     */
    void refresh_replyblock() {
        const std::size_t R = classes.size(), N = nodes.size();
        bool any = false;
        for (std::size_t r = 0; r < R && r < syncreply.size(); ++r)
            if (syncreply[r] >= 1 && syncreply[r] <= R) any = true;
        if (!any) return;
        replyblock.assign(N, std::vector<bool>(R, false));
        for (std::size_t r = 1; r <= R; ++r) {
            if (r > syncreply.size()) break;
            const std::size_t s = syncreply[r - 1];
            if (s < 1 || s > R) continue;
            for (std::size_t ind = 1; ind <= N; ++ind) {
                const std::size_t ist = nodes[ind - 1].station;
                if (ist == 0 || nodes[ind - 1].nodetype == lang::NodeType::Source) continue;
                if (stations[ist - 1].sched == SchedStrategy::INF) continue;
                bool arrives_here = false;
                for (std::size_t i = 1; i <= N && !arrives_here; ++i)
                    for (std::size_t q = 1; q <= R; ++q)
                        if (rtnodes.rows() >= N * R &&
                            num_traits<T>::to_double(
                                rtnodes((i - 1) * R + q - 1, (ind - 1) * R + s - 1)) > 0) {
                            arrives_here = true;
                            break;
                        }
                if (!arrives_here) continue;
                if (stations[ist - 1].sched != SchedStrategy::FCFS)
                    throw InputError(
                        "network '" + name +
                        "': synchronous calls (REPLY signals) are supported only at FCFS "
                        "stations, but '" +
                        nodes[ind - 1].name +
                        "' uses another discipline. A held server is encoded as a per-class "
                        "counter, which is exact only where servers are interchangeable (FCFS) "
                        "or unlimited (INF).");
                replyblock[ind - 1][r - 1] = true;
            }
        }
    }

    /**
     * The nodes directly downstream of `ind`, walking THROUGH stateless nodes
     * and stopping at the first station on each path.
     *
     * `downstreamStations` in `refreshLocalVars.m`. A blocked job is held for its
     * IMMEDIATE destination, so the walk stops at the first station: a Router or
     * ClassSwitch in between is a routing decision, not a place to wait.
     *
     * The reference reads `sn.connmatrix`; this port has no such field and reads
     * `rtnodes` instead, which is the same graph after the refresh has resolved
     * the routing strategies. The two differ only on a link the model declares
     * and then routes zero mass over, and such a link cannot fill a destination,
     * so it cannot block anything either.
     */
    std::vector<std::size_t> downstream_stations(std::size_t ind) const {
        const std::size_t R = classes.size();
        const std::size_t N = nodes.size();
        std::vector<std::size_t> out;
        if (rtnodes.rows() < N * R) return out;
        std::vector<bool> seen(N + 1, false);
        std::vector<std::size_t> frontier(1, ind);
        seen[ind] = true;
        while (!frontier.empty()) {
            const std::size_t i = frontier.back();
            frontier.pop_back();
            for (std::size_t j = 1; j <= N; ++j) {
                if (seen[j]) continue;
                bool linked = false;
                for (std::size_t r = 0; r < R && !linked; ++r)
                    for (std::size_t s = 0; s < R && !linked; ++s)
                        if (num_traits<T>::to_double(
                                rtnodes((i - 1) * R + r, (j - 1) * R + s)) > 0)
                            linked = true;
                if (!linked) continue;
                seen[j] = true;
                if (nodes[j - 1].station != 0) {
                    out.push_back(j);  // a station terminates this path
                } else {
                    frontier.push_back(j);  // walk through the stateless node
                }
            }
        }
        return out;
    }

    /**
     * Port of `refreshLocalVars`' true-BAS block and its `declaresBlockedMarker`
     * helper: which nodes carry the blocked marker, and where a refusal blocks.
     */
    void refresh_bas_blocking() {
        const std::size_t R = classes.size();
        const std::size_t N = nodes.size();
        isbasblocking.assign(N, false);
        isbasdestination.assign(stations.size(), std::vector<bool>(R, false));
        if (droprule.empty()) return;
        for (std::size_t ind = 1; ind <= N; ++ind) {
            const std::size_t ist = nodes[ind - 1].station;
            if (ist == 0) continue;
            const NodeType nt = nodes[ind - 1].nodetype;
            if (nt == NodeType::Source || nt == NodeType::Cache) continue;
            const std::vector<std::size_t> dests = downstream_stations(ind);
            bool declares = false;
            for (std::size_t r = 1; r <= R; ++r) {
                const bool here_bas = droprule.size() >= ist && droprule[ist - 1].size() >= r &&
                                      droprule[ist - 1][r - 1] == DropStrategy::BAS;
                for (std::size_t d = 0; d < dests.size(); ++d) {
                    const std::size_t jst = nodes[dests[d] - 1].station;
                    if (jst == 0 || nodes[dests[d] - 1].nodetype == NodeType::Source) continue;
                    // The destination must be able to FILL. An unbounded queue
                    // never refuses, so nothing upstream of it can ever block,
                    // and reserving the marker would widen the state for nothing.
                    //
                    // THE DECLARED CAPACITY, not the refreshed `cap`. The refresh
                    // clamps an unbounded station to the total closed population,
                    // which is finite, so reading `cap` here would find every
                    // destination in a closed model "able to fill" and reserve a
                    // marker on every upstream station. The reference reads
                    // `dnode.cap`, the value the user set.
                    const double dcap = stations[jst - 1].cap;
                    if (!std::isfinite(dcap) || !(dcap > 0)) continue;
                    const bool there_bas =
                        droprule.size() >= jst && droprule[jst - 1].size() >= r &&
                        droprule[jst - 1][r - 1] == DropStrategy::BAS;
                    if (!here_bas && !there_bas) continue;
                    // Do NOT stop at the first hit: every such destination has to
                    // be recorded, since a refusal at ANY of them must block here
                    // rather than drop.
                    declares = true;
                    isbasdestination[jst - 1][r - 1] = true;
                }
            }
            if (!declares) continue;
            // The marker shares nvars column 2R+1 with the cache contents, the
            // breakdown status and the polling controller, so the reference
            // rejects the combinations rather than widening the state.
            if (breakdownparam.find(ist) != breakdownparam.end())
                throw UnsupportedError(
                    "station '" + nodes[ind - 1].name +
                    "' combines server breakdowns with true-BAS blocking: the breakdown status "
                    "and the BAS blocked marker share one local-state column. Remove the BAS drop "
                    "rule or the breakdown");
            if (stations[ist - 1].sched == SchedStrategy::POLLING)
                throw UnsupportedError(
                    "true BAS blocking is not supported at the polling station '" +
                    nodes[ind - 1].name +
                    "': the polling controller and the BAS blocked marker share one local-state "
                    "column. Use a non-polling discipline at the blocking station, or remove the "
                    "BAS drop rule");
            for (std::size_t r = 1; r <= R; ++r)
                if (replyblock.size() >= ind && replyblock[ind - 1].size() >= r &&
                    replyblock[ind - 1][r - 1])
                    throw UnsupportedError(
                        "true BAS blocking is not supported at station '" + nodes[ind - 1].name +
                        "', which also holds servers for a synchronous reply: fromMarginal appends "
                        "the reply counters AFTER the blocked marker, so the marker would no "
                        "longer be the trailing column the departure handler reads");
            nvars[ind - 1][2 * R] = 1;
            isbasblocking[ind - 1] = true;
        }
    }

    /**
     * Port of MNetwork.refreshScheduling's schedparam half.
     *
     * DPS and GPS take a per-class weight, defaulting to 1; SEPT and LEPT take
     * the rank of the class's mean service time among the distinct means, which
     * is what the reference computes here rather than at the solver.
     *
     * EVERY OTHER DISCIPLINE DEFAULTS TO 1, NOT 0. Measured against MATLAB
     * R2025a: FCFS, PS, LPS, SIRO, LCFS, HOL and SJF all report schedparam 1
     * per class on a two-class model, and only SEPT/LEPT differ. Defaulting to
     * zero here made every consumer that normalizes by the weight total divide
     * 0/0: a load-dependent PS station reported Util exactly 0 against the
     * reference's 0.400822578299582, and `ctmc_signal_busy` shares the pattern.
     */
    void refresh_sched_param() {
        const T one = num_traits<T>::from_int(1);
        for (std::size_t i = 0; i < stations.size(); ++i) {
            Station<T>& st = stations[i];
            if (st.sched == SchedStrategy::DPS || st.sched == SchedStrategy::GPS) {
                if (st.schedparam.size() != classes.size())
                    st.schedparam.assign(classes.size(), one);
                continue;
            }
            if (st.sched == SchedStrategy::SEPT || st.sched == SchedStrategy::LEPT) {
                if (st.schedparam.size() == classes.size()) continue;
                std::vector<double> means;
                for (std::size_t r = 0; r < classes.size(); ++r)
                    means.push_back(num_traits<T>::to_double(service[i][r].mean));
                std::vector<double> sorted = means;
                std::sort(sorted.begin(), sorted.end());
                sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
                if (st.sched == SchedStrategy::LEPT)
                    std::reverse(sorted.begin(), sorted.end());
                st.schedparam.assign(classes.size(), num_traits<T>::from_int(0));
                for (std::size_t r = 0; r < classes.size(); ++r)
                    for (std::size_t k = 0; k < sorted.size(); ++k)
                        if (sorted[k] == means[r])
                            st.schedparam[r] = num_traits<T>::from_int(static_cast<long>(k) + 1);
                continue;
            }
            if (st.schedparam.empty()) st.schedparam.assign(classes.size(), one);
        }
    }

    /**
     * Port of the part of MNetwork.refreshRoutingMatrix this port reaches: the
     * expansion of a routing STRATEGY into the probabilities `rt` is built from.
     *
     * PROB is already probabilities and is copied through. RAND and RROBIN
     * spread the mass uniformly over the nodes the user connected this one to,
     * which is what MATLAB's `RoutingStrategy.RAND` means once `link` has
     * recorded the connections, and what `getRoutingMatrix.m:117` does for both
     * in the same branch: a round-robin pointer visits every outgoing link
     * equally often, so the ROUTING PROBABILITIES it induces are uniform and
     * only the higher moments of the split are deterministic. Recovering that
     * determinism is the consumer's job -- `npfqn_traffic_split_rr` gives QNA
     * and MNA the split degree, and a solver that needs the pointer itself must
     * carry it in the state. A solver that does neither must not declare
     * `RoutingStrategy_RROBIN` in its feature set, or it answers a random-
     * routing model under a round-robin name.
     * A ClassSwitch node's outgoing mass is multiplied by its
     * class-switch matrix, so a job leaving it in class r continues as class s
     * with probability C(r,s) -- the node itself is not stateful, and the
     * stochastic complement folds it into the edges around it.
     *
     * Every other strategy is state dependent (WRROBIN, JSQ, SQ, FIRING) and is
     * REFUSED by name: silently treating one as PROB returns a product-form
     * answer for a model that does not have one. WRROBIN is not RROBIN with
     * weights for this purpose -- its uniform expansion would be wrong even in
     * the first moment.
     */
    void refresh_routing() {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t K = classes.size(), I = nodes.size();
        Peff.clear();
        bool trivial = true;
        for (const NodeDef& nd : nodes)
            for (RoutingStrategy rs : nd.routing)
                if (rs != RoutingStrategy::PROB && rs != RoutingStrategy::DISABLED) trivial = false;
        // A (node, class) PAIR THE CALLER NEVER ROUTED IS NOT TRIVIAL EITHER,
        // because the reference FILLS it -- see the unrouted branch below. This
        // port has no DISABLED strategy to key on (every pair defaults to PROB),
        // so the condition is the one MATLAB's DISABLED case IS: no outgoing
        // mass for this class at this node, at a node that has connections.

        bool has_cs = false;
        for (const NodeDef& nd : nodes)
            if (nd.nodetype == NodeType::ClassSwitch) has_cs = true;
        // A Cache node is an implicit class switch: its read class routes on to
        // the hit and miss classes, so the routing is never trivial.
        const bool has_cache = !nodeparam.empty();
        if (trivial && !has_cs && !has_cache) return;  // P is already the effective routing

        for (const auto& kv : P) Peff[kv.first] = kv.second;
        auto eff_at = [&](std::size_t r, std::size_t s) -> Matrix<T>& {
            auto key = std::make_pair(r, s);
            auto it = Peff.find(key);
            if (it == Peff.end())
                it = Peff.emplace(key, Matrix<T>(I, I, zero)).first;
            if (it->second.rows() != I) {
                Matrix<T> g(I, I, zero);
                for (std::size_t a = 0; a < it->second.rows(); ++a)
                    for (std::size_t b = 0; b < it->second.cols(); ++b) g(a, b) = it->second(a, b);
                it->second = g;
            }
            return it->second;
        };

        for (std::size_t i = 1; i <= I; ++i) {
            const NodeDef& nd = nodes[i - 1];
            for (std::size_t r = 1; r <= K; ++r) {
                const RoutingStrategy rs =
                    nd.routing.size() >= r ? nd.routing[r - 1] : RoutingStrategy::PROB;
                if (rs == RoutingStrategy::PROB || rs == RoutingStrategy::DISABLED) continue;
                if (rs != RoutingStrategy::RAND && rs != RoutingStrategy::RROBIN &&
                    rs != RoutingStrategy::JSQ && rs != RoutingStrategy::SQ &&
                    rs != RoutingStrategy::WRROBIN && rs != RoutingStrategy::SDR)
                    throw UnsupportedError(std::string("network: routing strategy '") +
                                           routing_to_text(rs) + "' at node '" + nd.name +
                                           "' has no routing-matrix expansion in this port");
                // THE DESTINATIONS ARE THE NODE'S CONNECTIONS, NOT THIS CLASS'S
                // OWN ARCS. `getRoutingMatrix.m:117-135` enumerates `sn.connmatrix`
                // -- a class-INDEPENDENT topology -- and routes same-class over
                // it, which is what makes a class that declared no arc at a node
                // it nevertheless reaches get a routing at all. Reading the
                // class's own entries instead left such a row EMPTY, and an empty
                // row is not a visit: on `gallery_erlerl1`, where Class1 declares
                // Source->Queue and Class2 declares Queue->Sink, the Queue came
                // out unvisited by Class1 and every metric of the model was zero.
                // The `served` mask in `refresh_visits` is what keeps the fill
                // honest -- it drops the (station, class) pairs with no service.
                //
                // A CLOSED CLASS IS NOT ROUTED INTO A SINK, and not out of a
                // Source at all, which is the same branch's rule: a sink would
                // absorb a job the population must conserve.
                const bool open_class = !std::isfinite(
                    num_traits<T>::to_double(classes[r - 1].population));
                const bool from_source = nd.nodetype == NodeType::Source;
                const bool from_sink = nd.nodetype == NodeType::Sink;
                if (!open_class && (from_source || from_sink)) continue;
                std::vector<std::pair<std::size_t, std::size_t>> dest;  // (node, class)
                for (std::size_t j = 1; j <= I; ++j) {
                    if (!open_class && nodes[j - 1].nodetype == NodeType::Sink) continue;
                    bool connected = false;
                    for (std::size_t a = 1; a <= K && !connected; ++a)
                        for (std::size_t b = 1; b <= K && !connected; ++b)
                            if (get_route(a, b, i, j) > zero) connected = true;
                    if (connected) dest.emplace_back(j, r);
                }
                if (dest.empty()) continue;
                // WRROBIN SPREADS BY ITS WEIGHTS, everything else uniformly.
                // That is `getRoutingMatrix.m`'s split and it keeps the FIRST
                // MOMENT right for a weighted dispatcher, which a uniform
                // expansion would not.
                std::vector<T> share(dest.size(),
                                     T(num_traits<T>::from_int(1) /
                                       num_traits<T>::from_int(static_cast<long>(dest.size()))));
                if (rs == RoutingStrategy::WRROBIN) {
                    const std::map<std::size_t, double>* w =
                        nd.routing_weights.size() >= r ? &nd.routing_weights[r - 1] : NULL;
                    double total = 0.0;
                    std::vector<double> raw(dest.size(), 0.0);
                    if (w != NULL)
                        for (std::size_t d = 0; d < dest.size(); ++d) {
                            const std::map<std::size_t, double>::const_iterator it =
                                w->find(dest[d].first);
                            raw[d] = (it == w->end()) ? 0.0 : it->second;
                            total += raw[d];
                        }
                    if (total > 0.0)
                        for (std::size_t d = 0; d < dest.size(); ++d)
                            share[d] = num_traits<T>::from_double(raw[d] / total);
                }
                for (std::size_t s = 1; s <= K; ++s) {
                    Matrix<T>& B = eff_at(r, s);
                    for (std::size_t j = 1; j <= I; ++j) B(i - 1, j - 1) = zero;
                }
                for (std::size_t d = 0; d < dest.size(); ++d)
                    eff_at(r, dest[d].second)(i - 1, dest[d].first - 1) = share[d];
            }
        }

        // ClassSwitch: split the outgoing mass across the class-switch matrix.
        for (const auto& kv : csmatrix) {
            const std::size_t cs = kv.first;
            const Matrix<T>& C = kv.second;
            if (C.rows() != K || C.cols() != K)
                throw InputError("network: the class-switch matrix of node '" +
                                 nodes[cs - 1].name + "' is not (nclasses x nclasses)");
            // WHERE A SWITCHED JOB GOES IS THE ARRIVAL CLASS'S ROUTING, NOT THE
            // DEPARTURE CLASS'S. `getRoutingMatrix.m`'s StatelessClassSwitcher
            // block sets rtnodes(r, (j-1)*K+s) = Pcs(r,s) * Pij(s,s): the
            // destination distribution is read off the DIAGONAL, i.e. from the
            // row class s routes on with, and only the mass is Pcs(r,s).
            // Reading class r's own outgoing row instead is the same matrix
            // only when every class of the chain leaves the switch the same
            // way. Where they do not -- `Delay -> CS1 -> Queue -> CS2 -> Delay`
            // with each class Disabled at the station it does not visit, so
            // each class is routed on one arc of the cycle only -- class r has
            // NO outgoing arc at the switch it enters, the whole block comes
            // out zero, and the two classes then fall into separate chains with
            // empty visits. The state space collapses to its initial state and
            // the CTMC, NC and MVA answers are an empty table, not an error.
            std::vector<std::vector<T>> diag(K, std::vector<T>(I, zero));
            for (std::size_t s = 1; s <= K; ++s)
                for (std::size_t j = 1; j <= I; ++j)
                    diag[s - 1][j - 1] =
                        Peff.empty() ? get_route(s, s, cs, j) : eff_at(s, s)(cs - 1, j - 1);
            for (std::size_t r = 1; r <= K; ++r)
                for (std::size_t s = 1; s <= K; ++s) {
                    Matrix<T>& B = eff_at(r, s);
                    for (std::size_t j = 1; j <= I; ++j) B(cs - 1, j - 1) = zero;
                }
            for (std::size_t r = 1; r <= K; ++r)
                for (std::size_t s = 1; s <= K; ++s) {
                    if (!(C(r - 1, s - 1) > zero)) continue;
                    Matrix<T>& B = eff_at(r, s);
                    for (std::size_t j = 1; j <= I; ++j)
                        if (diag[s - 1][j - 1] > zero)
                            B(cs - 1, j - 1) = T(diag[s - 1][j - 1] * C(r - 1, s - 1));
                }
        }

        // Cache: the read (input) class self-switches at the cache node to the
        // hit and the miss class, which then follow their own routing. The
        // reference leaves the split unresolved (NaN) in the base struct and the
        // solver decides it; for the visit equations that back the offered
        // arrival rate and residence time it resolves to a uniform 1/2 - 1/2,
        // which is what makes ArvR the OFFERED rate (Tput is the carried one the
        // cacheqn decomposition produces). A half of the read mass reaching each
        // of hit/miss is a same-node class switch, hence a self-loop edge.
        const T half = T(num_traits<T>::from_int(1) / num_traits<T>::from_int(2));
        for (const auto& kv : nodeparam) {
            const std::size_t ci = kv.first;  // 1-based cache node
            const CacheParam<T>& cp = kv.second;
            for (std::size_t r = 0; r < cp.hitclass.size() && r < K; ++r) {
                if (cp.hitclass[r] == 0) continue;
                // clear the read class's own outgoing routing at the cache
                for (std::size_t s = 1; s <= K; ++s) {
                    Matrix<T>& B = eff_at(r + 1, s);
                    for (std::size_t j = 1; j <= I; ++j) B(ci - 1, j - 1) = zero;
                }
                eff_at(r + 1, cp.hitclass[r])(ci - 1, ci - 1) = half;
                if (r < cp.missclass.size() && cp.missclass[r] != 0)
                    eff_at(r + 1, cp.missclass[r])(ci - 1, ci - 1) = half;
            }
        }
    }

    /**
     * Port of MNetwork.refreshCapacity.
     *
     * `classcap(i,r)` is the population of r's CHAIN, cut down by any explicit
     * per-class or per-station buffer, and 0 where the class does not visit the
     * station; `cap(i)` is the explicit station buffer when there is one, and
     * otherwise the smaller of the chain and class sums.
     *
     * `chaincap` IS K COLUMNS WIDE, NOT nchains, exactly as the reference sizes
     * it (`chaincap = Inf*ones(M,K)`). The last class of a chain wins the write,
     * so a chain with a class disabled at station i can leave chaincap(i,c) at
     * 0; under class switching, where nchains < K, the untouched columns stay
     * Inf and carry the sum, which is what stops that 0 from capping a station
     * that holds the whole chain population at nothing.
     *
     * The derived drop rule follows the reference exactly, including that WAITQ
     * means "never consulted" at an unbounded station and "unsettled" for a
     * closed class at a bounded one: only an OPEN class at a real finite buffer
     * gets DROP.
     *
     * A PLACE IS EXEMPT FROM THE ZEROING. `disabled` is this port's marker for
     * the reference's `isnan(sn.rates(i,r))`, and a Place holds a marking
     * rather than serving, so it has no service process and every class reads
     * as disabled there. Zeroing it leaves a token container that cannot hold a
     * token: `cap` and `classcap` both come out 0, JMT is handed a Storage
     * section of capacity 0, and the net is dead on arrival. `refreshCapacity.m`
     * carries the same `~= NodeType.Place` guard on its `isnan` test.
     */
    /**
     * The number of siblings the fork-join pair ending at Join node `joinNode`
     * (1-based) emits per parent job: the matched Fork's out-degree times its
     * tasksPerLink, with the Join's in-degree as the fallback when no Fork is
     * matched. Port of `matlab/src/api/fj/sn_join_siblings.m`.
     *
     * Read off the DECLARED routing `P` rather than off `rtnodes`, because
     * `refresh_capacity` asks this question and runs BEFORE `refresh_rt`. The
     * two differ only on a declared link that carries no mass, which forks no
     * sibling either.
     */
    std::size_t join_siblings(std::size_t joinNode, std::size_t r = 0) const {
        const std::size_t I = nodes.size();
        if (joinNode == 0 || joinNode > I) return 0;
        std::size_t forkNode = 0;
        for (std::size_t a = 0; a < fj.size(); ++a)
            if (fj[a].second == joinNode) {
                forkNode = fj[a].first;
                break;
            }
        const T zero = num_traits<T>::from_int(0);
        std::vector<bool> seen(I + 1, false);
        std::size_t deg = 0;
        for (typename std::map<std::pair<std::size_t, std::size_t>, Matrix<T> >::const_iterator
                 it = P.begin();
             it != P.end(); ++it) {
            if (it->second.rows() < I || it->second.cols() < I) continue;
            for (std::size_t b = 1; b <= I; ++b) {
                if (seen[b]) continue;
                const T v = (forkNode != 0) ? it->second(forkNode - 1, b - 1)
                                            : it->second(b - 1, joinNode - 1);
                if (v > zero) {
                    seen[b] = true;
                    ++deg;
                }
            }
        }
        if (forkNode == 0 || forkNode > I) return deg;

        // THE COUNT IS PER LINK when the fork carries a VARIABLE FORKING LEVEL.
        // `fan_out_link(d,r)` is the expected tasks towards destination d for
        // class r (the DISTRIBUTION case stores its mean there) and
        // `fan_out_prob(d,r)` whether the branch fires at all, so the expected
        // sibling count is sum_d fan_out_prob(d,r)*fan_out_link(d,r). A fork
        // with no override has no `forkparam` entry, which is how the classic
        // out-degree times `tasks_per_link` case is told apart.
        const ForkParam<T>* fp = fork_param_of(forkNode);
        if (fp != 0 && fp->fan_out_link.rows() > 0) {
            const Matrix<T>& fol = fp->fan_out_link;
            const bool haveProb = fp->fan_out_prob.rows() == fol.rows() &&
                                  fp->fan_out_prob.cols() == fol.cols();
            const std::size_t lo = (r != 0 && r <= fol.cols()) ? r - 1 : 0;
            const std::size_t hi = (r != 0 && r <= fol.cols()) ? r - 1 : fol.cols() - 1;
            double best = 0.0;
            for (std::size_t c = lo; c <= hi && c < fol.cols(); ++c) {
                double acc = 0.0;
                for (std::size_t d = 0; d < fol.rows(); ++d) {
                    const double link = num_traits<T>::to_double(fol(d, c));
                    if (!(link > 0.0)) continue;
                    acc += link * (haveProb ? num_traits<T>::to_double(fp->fan_out_prob(d, c))
                                            : 1.0);
                }
                if (acc > best) best = acc;
            }
            if (best > 0.0) return static_cast<std::size_t>(best + 0.5);
        }

        double w = nodes[forkNode - 1].tasks_per_link;
        if (!(w >= 1.0)) w = 1.0;
        return deg * static_cast<std::size_t>(w + 0.5);
    }

    /**
     * The 1-based Join nodes that fire on a STRICT quorum, i.e. on FEWER
     * siblings than are forked. A model holding one is not
     * population-conserving at the sibling level: the join releases the parent
     * at the k-th of n siblings and the n-k stragglers stay in their branches,
     * so the parent forks again while they are still in flight.
     */
    std::vector<std::size_t> quorum_joins() const {
        std::vector<std::size_t> out;
        for (typename std::map<std::size_t, JoinDecl>::const_iterator it = joindecl.begin();
             it != joindecl.end(); ++it) {
            if (it->second.strategy != lang::JoinStrategy::PARTIAL) continue;
            if (it->second.quorum <= 0.0) continue;
            const std::size_t n = join_siblings(it->first);
            const std::size_t k = static_cast<std::size_t>(it->second.quorum + 0.5);
            if (n == 0 || k < n) out.push_back(it->first);
        }
        return out;
    }

    void refresh_capacity() {
        const std::size_t M = stations.size(), K = classes.size();
        const double inf = std::numeric_limits<double>::infinity();
        cap.assign(M, 0.0);
        classcap.assign(M, std::vector<double>(K, inf));
        droprule.assign(M, std::vector<DropStrategy>(K, DropStrategy::WAITQ));
        std::vector<std::vector<double>> chaincap(M, std::vector<double>(std::max(nchains, K), inf));

        // A class routed into a QUORUM Join is not population-conserving: the
        // stragglers of an already-fired parent are still in their branches
        // when it forks again, and nothing bounds that backlog, so a branch
        // station holds no more than the class population only under a STANDARD
        // join. Capping it at the chain population makes the engine drop a
        // closed job. see _kb/05-solvers-overview.md
        std::vector<char> quorum_class(K + 1, 0);
        {
            const std::vector<std::size_t> qj = quorum_joins();
            const T zero = num_traits<T>::from_int(0);
            for (std::size_t x = 0; x < qj.size(); ++x)
                for (typename std::map<std::pair<std::size_t, std::size_t>,
                                       Matrix<T> >::const_iterator it = P.begin();
                     it != P.end(); ++it) {
                    if (it->second.rows() < nodes.size() || qj[x] > nodes.size()) continue;
                    for (std::size_t a = 1; a <= nodes.size(); ++a)
                        if (it->second(a - 1, qj[x] - 1) > zero) {
                            if (it->first.first <= K) quorum_class[it->first.first] = 1;
                            if (it->first.second <= K) quorum_class[it->first.second] = 1;
                        }
                }
        }

        // A fork with tasksPerLink = w > 1 puts w tasks of the SAME parent on one
        // link, so a branch station can hold w jobs per circulating parent and the
        // chain population is no longer its bound. The multiplier is the PRODUCT
        // over the forks, because a fork nested in another's branch multiplies
        // again; that is an upper bound for forks in series, where a cap that never
        // binds costs nothing, and exact for the single-fork case. Without it this
        // engine drops a closed job at a branch station.
        double fork_task_factor = 1.0;
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            if (nodes[i].nodetype != NodeType::Fork) continue;
            const double w = nodes[i].tasks_per_link;
            if (w >= 1.0 && std::isfinite(w)) fork_task_factor *= std::floor(w + 0.5);
        }

        for (std::size_t c = 0; c < nchains; ++c) {
            double chain_cap = 0.0;
            bool open = false;
            bool quorum = false;
            for (std::size_t r : inchain[c]) {
                if (std::isinf(classes[r - 1].population)) open = true;
                else chain_cap += classes[r - 1].population;
                if (r <= K && quorum_class[r]) quorum = true;
            }
            chain_cap *= fork_task_factor;
            if (open || quorum) chain_cap = inf;
            for (std::size_t r : inchain[c])
                for (std::size_t i = 0; i < M; ++i) {
                    const Station<T>& st = stations[i];
                    const bool user_rule = st.droprule.size() >= r && st.droprule[r - 1] != 0;
                    if (st.nodetype != NodeType::Source) {
                        const bool cap_finite = st.cap >= 0.0 && !std::isinf(st.cap);
                        const bool classcap_finite = st.classcap.size() >= r &&
                                                     st.classcap[r - 1] > 0.0 &&
                                                     !std::isinf(st.classcap[r - 1]);
                        if (user_rule &&
                            st.droprule[r - 1] == static_cast<int>(DropStrategy::WAITQ) &&
                            std::isinf(classes[r - 1].population) &&
                            (cap_finite || classcap_finite))
                            throw UnsupportedError(
                                "station '" + st.name + "' declares setDropRule(WAITQ) for the "
                                "open class '" + classes[r - 1].name + "' at a finite capacity: "
                                "LINE does not implement waiting-room blocking for an open "
                                "arrival at a plain finite buffer. Use DropStrategy.DROP for a "
                                "loss station, or BAS / BBS / RSRD for blocking between "
                                "stations");
                        if (user_rule) {
                            droprule[i][r - 1] = static_cast<DropStrategy>(st.droprule[r - 1]);
                        } else if (std::isinf(st.cap)) {
                            droprule[i][r - 1] = DropStrategy::WAITQ;
                        } else if (!std::isinf(classes[r - 1].population)) {
                            droprule[i][r - 1] = DropStrategy::WAITQ;
                        } else {
                            droprule[i][r - 1] = DropStrategy::DROP;
                        }
                    }
                    // A Place has no service process, so `disabled` is not absence.
                    if (disabled[i][r - 1] && st.nodetype != NodeType::Place) {
                        classcap[i][r - 1] = 0.0;
                        chaincap[i][c] = 0.0;
                        continue;
                    }
                    chaincap[i][c] = chain_cap;
                    classcap[i][r - 1] = chain_cap;
                    if (st.classcap.size() >= r && st.classcap[r - 1] >= 0.0)
                        classcap[i][r - 1] = std::min(classcap[i][r - 1], st.classcap[r - 1]);
                    if (st.cap >= 0.0) classcap[i][r - 1] = std::min(classcap[i][r - 1], st.cap);
                }
        }
        for (std::size_t i = 0; i < M; ++i) {
            if (stations[i].cap >= 0.0 && !std::isinf(stations[i].cap)) {
                cap[i] = stations[i].cap;
                continue;
            }
            double sc = 0.0, scl = 0.0;
            for (std::size_t c = 0; c < chaincap[i].size(); ++c) sc += chaincap[i][c];
            for (std::size_t r = 0; r < K; ++r) scl += classcap[i][r];
            cap[i] = std::min(sc, scl);
        }
    }

    /**
     * `sn.rt` and `sn.rtnodes`: the class-expanded routing matrices.
     *
     * rtnodes is the routing over every node; rt is its stochastic complement
     * over the stateful ones, which is the same construction the per-chain
     * visits use -- but over ALL classes at once rather than one chain at a
     * time, because that is what `sn.rt` means and what solver_qna reads.
     */
    void refresh_rt() {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t K = classes.size(), I = nodes.size();
        std::vector<std::size_t> all(K);
        for (std::size_t r = 0; r < K; ++r) all[r] = r + 1;
        rtnodes = Matrix<T>(I * K, I * K, zero);
        for (std::size_t a = 0; a < I; ++a)
            for (std::size_t r = 0; r < K; ++r)
                for (std::size_t b = 0; b < I; ++b)
                    for (std::size_t s = 0; s < K; ++s)
                        rtnodes(a * K + r, b * K + s) = route_eff(r + 1, s + 1, a + 1, b + 1);
        rt = station_routing(all);
    }

    /**
     * Re-resolve the cache read self-switch from the offered 1/2-1/2 to the
     * ACTUAL hit/miss probabilities the cacheqn decomposition converged on, then
     * recompute `rt` and the visits. The runner reads ArvR and ResidT off the
     * result, matching MATLAB whose `sn.visits` carry the actual
     * (setResultHitProb) split, not the offered one. Served metrics are
     * unaffected: they come from the analyzer's own over-routed inner solve.
     * `hitprob`/`missprob` are (ncaches x nclasses), indexed by cache order and
     * input class, exactly as `da_cacheqn` returns them.
     */
    void refresh_cacheqn_actual_visits(const Matrix<T>& hitprob, const Matrix<T>& missprob) {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t K = classes.size();
        std::size_t cidx = 0;
        for (const auto& kv : nodeparam) {
            const std::size_t ci = kv.first;  // 1-based cache node
            const CacheParam<T>& cp = kv.second;
            for (std::size_t r = 0; r < cp.hitclass.size() && r < K; ++r) {
                if (cp.hitclass[r] == 0) continue;
                set_route_effective(r + 1, cp.hitclass[r], ci, ci, hitprob(cidx, r));
                if (r < cp.missclass.size() && cp.missclass[r] != 0)
                    set_route_effective(r + 1, cp.missclass[r], ci, ci, missprob(cidx, r));
            }
            ++cidx;
        }
        refresh_rt();
        const bool fork = has_fork();
        visits.assign(nchains, Matrix<T>(stateful_nodes.size(), nclasses, zero));
        nodevisits.assign(nchains, Matrix<T>(nodes.size(), nclasses, zero));
        for (std::size_t c = 0; c < nchains; ++c) {
            visits[c] = chain_visits(c, stateful_nodes, true, fork);
            nodevisits[c] = chain_visits(c, all_nodes(), false, fork);
        }
    }

    /**
     * Port of MNetwork.refreshRates: lower each service process onto a rate and
     * an SCV. A disabled (station, class) pair becomes NaN, which is the marker
     * every downstream consumer keys on to mean "this class never visits here".
     */
    void refresh_rates() {
        nstations = stations.size();
        nclasses = classes.size();
        const T zero = num_traits<T>::from_int(0);
        rates = Matrix<T>(nstations, nclasses, zero);
        scv = Matrix<T>(nstations, nclasses, zero);
        disabled.assign(nstations, std::vector<bool>(nclasses, true));
        for (std::size_t i = 0; i < nstations; ++i)
            for (std::size_t r = 0; r < nclasses; ++r) {
                // Join infinite-rate rationale: see _kb/04-networkstruct.md (cpp port notes)
                if (stations[i].nodetype == NodeType::Join) {
                    disabled[i][r] = false;
                    rates(i, r) = num_traits<T>::from_double(
                        std::numeric_limits<double>::infinity());
                    scv(i, r) = zero;
                    continue;
                }
                const Distrib<T>& d = service[i][r];
                if (d.disabled) continue;
                disabled[i][r] = false;
                rates(i, r) = d.rate();
                scv(i, r) = d.scv;
            }
    }

    /**
     * Port of MNetwork.refreshChains followed by sn_refresh_visits.
     *
     * The class-switch mask is (r == s) or "some link carries r into s", which
     * is what link() records in csMatrix and what refreshChains then re-derives
     * from `rt`; the chains are the connected components of that mask read as
     * an undirected graph. MATLAB orders the chains by `sortrows(...,'descend')`
     * on the indicator rows, which for disjoint components is the same as
     * ordering by their smallest class index, and that is what is done here.
     */
    void refresh_chains() {
        refresh_rates();
        const T zero = num_traits<T>::from_int(0);
        const std::size_t K = nclasses;

        // class-switch mask rationale: see _kb/04-networkstruct.md (cpp port notes)
        std::vector<std::vector<bool>> cs(K, std::vector<bool>(K, false));
        for (std::size_t r = 0; r < K; ++r) cs[r][r] = true;
        for (const auto& kv : (Peff.empty() ? P : Peff)) {
            bool any = false;
            for (std::size_t a = 0; a < kv.second.rows() && !any; ++a)
                for (std::size_t b = 0; b < kv.second.cols() && !any; ++b)
                    if (kv.second(a, b) > zero) any = true;
            if (any) cs[kv.first.first - 1][kv.first.second - 1] = true;
        }
        // A ClassSwitch NODE couples classes too, and once `link()` synthesizes
        // one the coupling lives ONLY here: the rewrite folds `P{r,s}(i,j)` into
        // a same-class pair of legs, so the r != s block that used to carry it
        // is gone from P. `link.m` accumulates the same thing --
        // `csMatrix = csMatrix | nodes{ind}.server.csMatrix > 0` -- and without
        // it every switched class falls into its own chain.
        for (const auto& kv : csmatrix) {
            const Matrix<T>& C = kv.second;
            for (std::size_t r = 0; r < K && r < C.rows(); ++r)
                for (std::size_t s = 0; s < K && s < C.cols(); ++s)
                    if (C(r, s) > zero) cs[r][s] = true;
        }
        // A Cache node couples each input class with the classes it switches to
        // on a hit and a miss; that coupling lives in nodeparam, not in P, so it
        // must be added here or the read/hit/miss classes fall into separate
        // chains and the switched (open) classes get no arrivals.
        for (const auto& kv : nodeparam) {
            const CacheParam<T>& cp = kv.second;
            for (std::size_t r = 0; r < cp.hitclass.size() && r < K; ++r) {
                if (cp.hitclass[r] != 0) cs[r][cp.hitclass[r] - 1] = true;
                if (r < cp.missclass.size() && cp.missclass[r] != 0)
                    cs[r][cp.missclass[r] - 1] = true;
            }
            // A RETRIEVAL class is also a class the cache switches the read
            // class into on a miss -- it just lives in `retrieval_classes`
            // rather than in `missclass`, because there is one per ITEM. The
            // coupling above missed them, so on a closed delayed-hit model each
            // retrieval class formed its own singleton chain: 4 chains where
            // MATLAB has 1, with zero visits at the fetch station and the whole
            // population parked at the delay. The rule in the comment above
            // applies to them unchanged.
            for (const std::vector<std::size_t>& row : cp.retrieval_classes)
                for (std::size_t r = 0; r < row.size() && r < K; ++r)
                    if (row[r] != 0) cs[r][row[r] - 1] = true;
        }
        std::vector<std::size_t> comp(K, K);
        std::size_t ncomp = 0;
        for (std::size_t r = 0; r < K; ++r) {
            if (comp[r] != K) continue;
            std::vector<std::size_t> stack{r};
            comp[r] = ncomp;
            while (!stack.empty()) {
                const std::size_t v = stack.back();
                stack.pop_back();
                for (std::size_t w = 0; w < K; ++w)
                    if (comp[w] == K && (cs[v][w] || cs[w][v])) {
                        comp[w] = ncomp;
                        stack.push_back(w);
                    }
            }
            ++ncomp;
        }
        // components are already discovered in order of their smallest member,
        // which is the order sortrows(...,'descend') produces
        nchains = ncomp;
        chains.assign(nchains, std::vector<bool>(K, false));
        inchain.assign(nchains, {});
        for (std::size_t r = 0; r < K; ++r) {
            chains[comp[r]][r] = true;
            inchain[comp[r]].push_back(r + 1);
        }

        // ---- reference class per chain -------------------------------------
        refclass.assign(nchains, 0);
        for (std::size_t c = 0; c < nchains; ++c)
            for (std::size_t k : inchain[c])
                if (classes[k - 1].is_ref_class) refclass[c] = k;

        for (std::size_t c = 0; c < nchains; ++c) {
            const std::size_t rs = classes[inchain[c][0] - 1].refstat;
            for (std::size_t k : inchain[c])
                if (classes[k - 1].refstat != rs)
                    throw InputError("network '" + name + "': classes within a chain have different "
                                     "reference stations");
        }

        apply_sink_closure();

        // ---- visits, at stateful-node and at node level ---------------------
        const bool fork = has_fork();
        visits.assign(nchains, Matrix<T>(stateful_nodes.size(), nclasses, zero));
        nodevisits.assign(nchains, Matrix<T>(nodes.size(), nclasses, zero));
        for (std::size_t c = 0; c < nchains; ++c) {
            visits[c] = chain_visits(c, stateful_nodes, true, fork);
            nodevisits[c] = chain_visits(c, all_nodes(), false, fork);
        }
    }

    /**
     * Route every open chain from the Sink back into the Source, as the tail of
     * MATLAB's getRoutingMatrix does before it takes the stochastic complement.
     *
     * WHY IT IS NOT OPTIONAL. `visits` is the stationary vector of the chain's
     * routing DTMC. An open chain that ends at the Sink has no such vector --
     * the Sink is absorbing, so the solve returns the point mass there and
     * every station comes out with zero visits. The closing arc turns the open
     * chain into a recurrent one whose stationary vector, renormalised by the
     * Source, is the visit count PER ARRIVAL, which is what the chain demands
     * are built from.
     *
     * The destination class is drawn in proportion to the ARRIVAL RATES of the
     * chain's classes, so a job leaving the Sink re-enters as the class the
     * Source would have generated. A chain whose arrivals are all disabled has
     * no such proportion, and the reference falls back to the uniform choice
     * rather than dividing by zero -- those chains carry no traffic and only
     * need to stay well posed.
     *
     * The arcs are rebuilt on every refresh, and the previous ones cleared
     * first: the rates move between passes of the fork-join fixed point, and a
     * stale arc would leave two closures with different weights in place.
     */
    void apply_sink_closure() {
        if (sourceIdx == 0 || sinkNode == 0) return;
        const T zero = num_traits<T>::from_int(0);
        const T one = num_traits<T>::from_int(1);
        const std::size_t src = station_to_node[sourceIdx - 1];
        // closure-clearing rationale: see _kb/04-networkstruct.md (cpp port notes)
        for (auto& kv : P) {
            if (kv.second.rows() < nodes.size()) continue;
            kv.second(sinkNode - 1, src - 1) = zero;
        }
        for (auto& kv : Peff) {
            if (kv.second.rows() < nodes.size()) continue;
            kv.second(sinkNode - 1, src - 1) = zero;
        }
        for (std::size_t c = 0; c < nchains; ++c) {
            bool open = false;
            for (std::size_t k : inchain[c])
                if (std::isinf(classes[k - 1].population)) open = true;
            if (!open) continue;
            T tot = zero;
            for (std::size_t k : inchain[c])
                if (!disabled[sourceIdx - 1][k - 1]) tot += rates(sourceIdx - 1, k - 1);
            const T uniform =
                T(one / num_traits<T>::from_int(static_cast<long>(inchain[c].size())));
            for (std::size_t s : inchain[c]) {
                T p = uniform;
                if (tot > zero) {
                    const T ar =
                        disabled[sourceIdx - 1][s - 1] ? zero : rates(sourceIdx - 1, s - 1);
                    p = T(ar / tot);
                }
                if (!(p > zero)) continue;
                for (std::size_t r : inchain[c]) set_route_effective(r, s, sinkNode, src, p);
            }
        }
    }

    /** Every node index, 1-based, for the node-level visit computation. */
    std::vector<std::size_t> all_nodes() const {
        std::vector<std::size_t> v(nodes.size());
        for (std::size_t i = 0; i < nodes.size(); ++i) v[i] = i + 1;
        return v;
    }

    /**
     * Port of the per-chain body of sn_refresh_visits, over an arbitrary node
     * subset -- the stateful nodes for `visits`, all nodes for `nodevisits`.
     *
     * The FORK CORRECTIONS are the reference's, and they are blunt on purpose.
     * A Fork row sums to its fan-out rather than to one, so the rows are
     * renormalised before the DTMC solve; afterwards, rather than trusting the
     * resulting stationary vector, a population-preserving SPN argument sets
     * EVERY visited entry to 1 (and a Join entry to its in-degree, one per
     * incoming branch). Reproduced exactly: the corrected visits are what the
     * chain demands are built from, so an "improved" version would disagree
     * with every other codebase.
     */
    Matrix<T> chain_visits(std::size_t c, const std::vector<std::size_t>& sel, bool complement,
                           bool fork) const {
        const T zero = num_traits<T>::from_int(0);
        const std::vector<std::size_t>& ic = inchain[c];
        const std::size_t nIC = ic.size();
        const std::size_t dim = sel.size() * nIC;
        Matrix<T> Pc(dim, dim, zero);
        if (complement) {
            const Matrix<T> rt = station_routing(ic);
            Pc = rt;
        } else {
            for (std::size_t a = 0; a < sel.size(); ++a)
                for (std::size_t x = 0; x < nIC; ++x)
                    for (std::size_t b = 0; b < sel.size(); ++b)
                        for (std::size_t y = 0; y < nIC; ++y)
                            Pc(a * nIC + x, b * nIC + y) = route_eff(ic[x], ic[y], sel[a], sel[b]);
        }

        // THE `served` MASK of sn_refresh_visits.m, ON BOTH CHAINS.
        // `getRoutingMatrix` leaves a JMT-oriented UNIFORM FILL on (station,
        // class) pairs the class has no service at, and those states are not
        // reachable: a class with no service law cannot be served there. Left
        // in, they leak mass between what are otherwise separate recurrent
        // branches -- on cs_transient_class that moved Queue1 from 0.1875 to
        // 0.22917 and Queue2 from 0.3125 to 0.27083, against an absorption of
        // 0.375/0.625 halved over each two-station cycle.
        //
        // It matters MORE on the node chain, which the reference used to leave
        // untouched: at station level an unserved (station,class) is a dead end,
        // while the node kernel keeps the class-switch nodes between the
        // stations, so the disabled states close into a whole spurious CYCLE. A
        // materialised LQN replica is exactly that -- replica 2's stations still
        // carry replica 1's classes in rtnodes -- and the reducible solve then
        // splits the mass between the real chain and the phantom one, giving
        // every node of replica 2 a visit in replica 1's classes.
        //
        // A Place, a Transition and a station declaring heterogeneous server
        // types are EXEMPT, exactly as the reference exempts them: all three
        // carry NaN station rates by construction, so the test would read as
        // "not served" for a station that plainly is.
        std::vector<bool> served(dim, true);
        {
            for (std::size_t a = 0; a < sel.size(); ++a) {
                const std::size_t sti = nodes[sel[a] - 1].station;
                if (sti == 0 || sti > stations.size()) continue;
                const NodeType nt = stations[sti - 1].nodetype;
                if (nt == NodeType::Place || nt == NodeType::Transition) continue;
                if (!stations[sti - 1].server_types.empty()) continue;
                if (sti > disabled.size()) continue;
                // The reference tests `isnan(sn.rates(...))`; THIS PORT'S MARKER
                // IS `disabled`. refresh_rates leaves a disabled pair's rate at
                // 0 rather than NaN (its comment says otherwise and is stale),
                // so an isnan test here would never fire and the mask would be
                // silently inert. `disabled` is what every other consumer keys
                // on and it carries exactly MATLAB's meaning.
                for (std::size_t x = 0; x < nIC; ++x)
                    if (ic[x] <= disabled[sti - 1].size() && disabled[sti - 1][ic[x] - 1])
                        served[a * nIC + x] = false;
            }
            for (std::size_t row = 0; row < dim; ++row)
                if (!served[row])
                    for (std::size_t col = 0; col < dim; ++col) Pc(row, col) = zero;
            for (std::size_t col = 0; col < dim; ++col)
                if (!served[col])
                    for (std::size_t row = 0; row < dim; ++row) Pc(row, col) = zero;
        }

        std::vector<std::size_t> visited;
        std::vector<T> rowsum(dim, zero);
        for (std::size_t row = 0; row < dim; ++row) {
            T s = zero;
            for (std::size_t col = 0; col < dim; ++col) s += Pc(row, col);
            rowsum[row] = s;
            if (s > zero) visited.push_back(row);
        }
        bool oversum = false;
        if (fork) {
            for (std::size_t row = 0; row < dim; ++row) {
                if (num_traits<T>::to_double(rowsum[row]) > 1.0 + GlobalConstants::FineTol)
                    oversum = true;
                if (num_traits<T>::to_double(rowsum[row]) > GlobalConstants::FineTol)
                    for (std::size_t col = 0; col < dim; ++col)
                        Pc(row, col) = T(Pc(row, col) / rowsum[row]);
            }
        }

        // Detect a genuinely reducible chain (more than one recurrent class)
        // up front, so its DISABLED-routing fillers and the reducible solver
        // apply ONLY here -- every irreducible chain, i.e. the whole existing
        // test surface, stays on the exact dtmc_solve path below. In MATLAB the
        // singular normalization solve returns NaN here, tripping the same
        // fallback; this reproduces that without changing the shared solver.
        bool reducible = false;
        if (!fork && !visited.empty()) {
            Matrix<T> Pv0(visited.size(), visited.size(), zero);
            for (std::size_t a = 0; a < visited.size(); ++a)
                for (std::size_t b = 0; b < visited.size(); ++b) Pv0(a, b) = Pc(visited[a], visited[b]);
            const mc::SccResult scc = mc::stronglyconncomp(Pv0);
            std::size_t nrec = 0;
            for (bool rc : scc.recurrent)
                if (rc) ++nrec;
            reducible = (nrec > 1);
        }
        if (reducible) {
            // DISABLED-routing filler (getRoutingMatrix.m DISABLED case): a
            // (node,class) with no outgoing routing at a physically connected node
            // routes SAME-CLASS to each neighbour at 1/nconn. These 0-visit
            // transient states give the reducible solver the multiple transient
            // SCCs MATLAB has, so it weights the recurrent classes correctly
            // (11:13) instead of the single-transient absorption split (3:5).
            std::vector<std::size_t> nconn(sel.size(), 0);
            std::vector<std::vector<bool>> conn(sel.size(), std::vector<bool>(sel.size(), false));
            for (std::size_t a = 0; a < sel.size(); ++a)
                for (std::size_t b = 0; b < sel.size(); ++b) {
                    bool cbit = false;
                    for (std::size_t x = 0; x < nIC && !cbit; ++x)
                        for (std::size_t y = 0; y < nIC && !cbit; ++y)
                            if (Pc(a * nIC + x, b * nIC + y) > zero) cbit = true;
                    conn[a][b] = cbit;
                    if (cbit) ++nconn[a];
                }
            for (std::size_t a = 0; a < sel.size(); ++a) {
                if (nconn[a] == 0) continue;
                for (std::size_t x = 0; x < nIC; ++x) {
                    T outsum = zero;
                    for (std::size_t col = 0; col < dim; ++col) outsum += Pc(a * nIC + x, col);
                    if (outsum > zero) continue;  // class already routes onward here
                    const T p = num_traits<T>::from_int(1) /
                                num_traits<T>::from_int(static_cast<long>(nconn[a]));
                    for (std::size_t b = 0; b < sel.size(); ++b)
                        if (conn[a][b]) Pc(a * nIC + x, b * nIC + x) = p;
                }
            }
            // AND THE FILLER IS MASKED AGAIN, which is the whole point of the
            // mask. The reference's order is: getRoutingMatrix lays the fill
            // down, THEN `served` zeroes those rows and columns -- so an
            // unserved state never carries fill in the matrix that is solved.
            // Synthesising the fill here and stopping would restore exactly what
            // the mask was ported to remove: on cs_transient_class the two
            // recurrent branches came out 11:13 (0.22917 / 0.27083) instead of
            // the absorption split 3:5 (0.1875 / 0.3125) that MATLAB, Java and
            // Python all report and that the routing implies.
            for (std::size_t row = 0; row < dim; ++row)
                if (!served[row])
                    for (std::size_t col = 0; col < dim; ++col) Pc(row, col) = zero;
            for (std::size_t col = 0; col < dim; ++col)
                if (!served[col])
                    for (std::size_t row = 0; row < dim; ++row) Pc(row, col) = zero;
            visited.clear();
            for (std::size_t row = 0; row < dim; ++row) {
                T s = zero;
                for (std::size_t col = 0; col < dim; ++col) s += Pc(row, col);
                if (s > zero) visited.push_back(row);
            }
        }

        std::vector<T> alpha(dim, zero);
        if (!visited.empty()) {
            Matrix<T> Pv(visited.size(), visited.size(), zero);
            for (std::size_t a = 0; a < visited.size(); ++a)
                for (std::size_t b = 0; b < visited.size(); ++b)
                    Pv(a, b) = Pc(visited[a], visited[b]);
            std::vector<T> av;
            bool ok = false;
            if (!reducible) {
                try {
                    av = mc::dtmc_solve(Pv);
                    ok = true;
                    bool allzero = true, hasnan = false;
                    for (const T& x : av) {
                        if (x != zero) allzero = false;
                        if (num_traits<T>::to_double(x) != num_traits<T>::to_double(x)) hasnan = true;
                    }
                    if (allzero || hasnan) ok = false;
                } catch (const Error&) {
                    ok = false;
                }
            }
            if (!ok) {
                // DTMC solver order: see _kb/04-networkstruct.md (cpp port notes) and _kb/11-conventions-and-gotchas.md
                av = mc::dtmc_solve_reducible(Pv, GlobalConstants::FineTol).pi;
            }
            for (std::size_t a = 0; a < visited.size(); ++a) alpha[visited[a]] = av[a];
        }

        if (fork && oversum) {
            for (std::size_t idx = 0; idx < dim; ++idx) {
                if (!(num_traits<T>::to_double(alpha[idx]) > GlobalConstants::FineTol)) continue;
                const std::size_t nd = sel[idx / nIC];
                if (!complement && nodes[nd - 1].nodetype == NodeType::Join) {
                    // a Join is entered once per incoming branch. The reference
                    // counts POSITIVE ROWS of rtnodes in column (nd,r), i.e.
                    // (source node, source class) PAIRS over every class -- two
                    // classes entering the Join as r from the same predecessor
                    // are two branches, so neither the break nor the restriction
                    // to the chain's own classes belongs here.
                    const std::size_t r = ic[idx % nIC];
                    std::size_t nsrc = 0;
                    for (std::size_t src = 1; src <= nodes.size(); ++src)
                        for (std::size_t q = 1; q <= nclasses; ++q)
                            if (num_traits<T>::to_double(route_eff(q, r, src, nd)) >
                                GlobalConstants::FineTol)
                                ++nsrc;
                    alpha[idx] = num_traits<T>::from_int(static_cast<long>(nsrc));
                } else {
                    alpha[idx] = num_traits<T>::from_int(1);
                }
            }
        }

        Matrix<T> out(sel.size(), nclasses, zero);
        for (std::size_t a = 0; a < sel.size(); ++a)
            for (std::size_t x = 0; x < nIC; ++x) out(a, ic[x] - 1) = alpha[a * nIC + x];

        // normalise by the total visits of the chain at its reference node
        const std::size_t rstat = classes[ic[0] - 1].refstat;
        const std::size_t refnode = station_to_node[rstat - 1];
        std::size_t refrow = sel.size();
        for (std::size_t a = 0; a < sel.size(); ++a)
            if (sel[a] == refnode) refrow = a;
        if (refrow < sel.size()) {
            T norm = zero;
            for (std::size_t x = 0; x < nIC; ++x) norm += out(refrow, ic[x] - 1);
            if (num_traits<T>::to_double(norm) > GlobalConstants::FineTol)
                for (std::size_t a = 0; a < sel.size(); ++a)
                    for (std::size_t x = 0; x < nIC; ++x)
                        out(a, ic[x] - 1) = T(out(a, ic[x] - 1) / norm);
        }
        for (std::size_t a = 0; a < sel.size(); ++a)
            for (std::size_t x = 0; x < nIC; ++x)
                if (out(a, ic[x] - 1) < zero) out(a, ic[x] - 1) = T(-out(a, ic[x] - 1));
        return out;
    }

    /**
     * The visit ratios of one chain from an already-formed chain routing block
     * `Pc` (dim = sel.size() * nIC), the no-fork body of sn_refresh_visits: solve
     * the embedded DTMC over the visited states and normalise by the total visits
     * at the reference node. NaN entries (a Cache class switch can leave one) are
     * given the row's residual probability spread uniformly, as the reference
     * does. Fork models keep the route_eff path (`chain_visits`); this is used by
     * the cacheqn driver, which rewrites `rtnodes` directly and has no fork.
     */
    Matrix<T> visits_from_block(const Matrix<T>& Pc_in, const std::vector<std::size_t>& sel,
                                const std::vector<std::size_t>& ic, std::size_t refnode) const {
        const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
        const std::size_t nIC = ic.size();
        const std::size_t dim = Pc_in.rows();
        Matrix<T> Pc = Pc_in;
        // NaN handling: distribute the row's residual mass over the NaN columns.
        for (std::size_t row = 0; row < dim; ++row) {
            T nonnan = zero;
            std::size_t nnan = 0;
            for (std::size_t col = 0; col < dim; ++col) {
                const double v = num_traits<T>::to_double(Pc(row, col));
                if (v != v) ++nnan;
                else nonnan = T(nonnan + Pc(row, col));
            }
            if (nnan == 0) continue;
            const double rem = 1.0 - num_traits<T>::to_double(nonnan);
            const T fill = (rem > 0.0) ? num_traits<T>::from_double(rem / double(nnan)) : zero;
            for (std::size_t col = 0; col < dim; ++col) {
                const double v = num_traits<T>::to_double(Pc(row, col));
                if (v != v) Pc(row, col) = fill;
            }
        }
        std::vector<std::size_t> visited;
        for (std::size_t row = 0; row < dim; ++row) {
            T s = zero;
            for (std::size_t col = 0; col < dim; ++col) s += Pc(row, col);
            if (s > zero) visited.push_back(row);
        }
        std::vector<T> alpha(dim, zero);
        if (!visited.empty()) {
            Matrix<T> Pv(visited.size(), visited.size(), zero);
            for (std::size_t a = 0; a < visited.size(); ++a)
                for (std::size_t b = 0; b < visited.size(); ++b) Pv(a, b) = Pc(visited[a], visited[b]);
            std::vector<T> av;
            bool ok = false;
            try {
                av = mc::dtmc_solve(Pv);
                ok = true;
                bool allzero = true, hasnan = false;
                for (const T& x : av) {
                    if (x != zero) allzero = false;
                    if (num_traits<T>::to_double(x) != num_traits<T>::to_double(x)) hasnan = true;
                }
                if (allzero || hasnan) ok = false;
            } catch (const Error&) {
                ok = false;
            }
            if (!ok) av = mc::dtmc_solve_reducible(Pv, GlobalConstants::FineTol).pi;
            for (std::size_t a = 0; a < visited.size(); ++a) alpha[visited[a]] = av[a];
        }
        Matrix<T> out(sel.size(), nclasses, zero);
        for (std::size_t a = 0; a < sel.size(); ++a)
            for (std::size_t x = 0; x < nIC; ++x) out(a, ic[x] - 1) = alpha[a * nIC + x];
        std::size_t refrow = sel.size();
        for (std::size_t a = 0; a < sel.size(); ++a)
            if (sel[a] == refnode) refrow = a;
        if (refrow < sel.size()) {
            T norm = zero;
            for (std::size_t x = 0; x < nIC; ++x) norm += out(refrow, ic[x] - 1);
            if (num_traits<T>::to_double(norm) > GlobalConstants::FineTol)
                for (std::size_t a = 0; a < sel.size(); ++a)
                    for (std::size_t x = 0; x < nIC; ++x)
                        out(a, ic[x] - 1) = T(out(a, ic[x] - 1) / norm);
        }
        for (std::size_t a = 0; a < sel.size(); ++a)
            for (std::size_t x = 0; x < nIC; ++x)
                if (out(a, ic[x] - 1) < zero) out(a, ic[x] - 1) = T(-out(a, ic[x] - 1));
        return out;
    }

    /**
     * Recompute `rt`, `visits` and `nodevisits` after the caller has rewritten
     * `rtnodes` in place -- the cacheqn driver's per-sweep refresh. `rt` is the
     * stochastic complement of `rtnodes` over the stateful (node, class) rows;
     * the per-chain visits then read directly off `rt` (stateful) and `rtnodes`
     * (all nodes), with the chains held FIXED (the reference does not recompute
     * them inside the sweep). No-fork only; a fork keeps the refresh path.
     */
    void da_recompute_visits_from_rtnodes() {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t K = nclasses, I = nodes.size(), M = stateful_nodes.size();
        for (const NodeDef& nd : nodes)
            if (nd.nodetype == NodeType::Fork)
                throw UnsupportedError("da_recompute_visits_from_rtnodes: fork models keep the "
                                       "route_eff visit path");
        std::vector<std::size_t> keep;
        keep.reserve(M * K);
        for (std::size_t p = 0; p < M; ++p) {
            const std::size_t base = (stateful_nodes[p] - 1) * K;
            for (std::size_t r = 0; r < K; ++r) keep.push_back(base + r);
        }
        rt = mc::dtmc_stochcomp(rtnodes, keep);

        std::vector<std::size_t> allnodes(I);
        for (std::size_t a = 0; a < I; ++a) allnodes[a] = a + 1;
        visits.assign(nchains, Matrix<T>(M, K, zero));
        nodevisits.assign(nchains, Matrix<T>(I, K, zero));
        for (std::size_t c = 0; c < nchains; ++c) {
            const std::vector<std::size_t>& ic = inchain[c];
            const std::size_t nIC = ic.size();
            // stateful visits from rt (ordered by stateful position)
            Matrix<T> Ps(M * nIC, M * nIC, zero);
            for (std::size_t a = 0; a < M; ++a)
                for (std::size_t x = 0; x < nIC; ++x)
                    for (std::size_t b = 0; b < M; ++b)
                        for (std::size_t y = 0; y < nIC; ++y)
                            Ps(a * nIC + x, b * nIC + y) = rt(a * K + (ic[x] - 1), b * K + (ic[y] - 1));
            const std::size_t rstat = classes[ic[0] - 1].refstat;
            const std::size_t refnode = station_to_node[rstat - 1];
            // sel for the stateful block is the stateful node list; its "refnode"
            // position must be matched by node index, so pass stateful_nodes.
            visits[c] = visits_from_block(Ps, stateful_nodes, ic, refnode);
            // node visits from rtnodes (ordered by node)
            Matrix<T> Pn(I * nIC, I * nIC, zero);
            for (std::size_t a = 0; a < I; ++a)
                for (std::size_t x = 0; x < nIC; ++x)
                    for (std::size_t b = 0; b < I; ++b)
                        for (std::size_t y = 0; y < nIC; ++y)
                            Pn(a * nIC + x, b * nIC + y) =
                                rtnodes(a * K + (ic[x] - 1), b * K + (ic[y] - 1));
            nodevisits[c] = visits_from_block(Pn, allnodes, ic, refnode);
        }
    }

    /**
     * The chain-restricted routing over the STATEFUL nodes: the node-level
     * routing with the non-stateful nodes eliminated by a stochastic
     * complement, S = P11 + P12 (I - P22)^-1 P21. That is MATLAB's
     * dtmc_stochcomp, and it is what removes a Fork (and, there, the auto-added
     * ClassSwitch nodes) from the visit equations.
     */
    Matrix<T> station_routing(const std::vector<std::size_t>& ic) const {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t nIC = ic.size();
        const std::size_t I = nodes.size();
        Matrix<T> full(I * nIC, I * nIC, zero);
        for (std::size_t a = 0; a < I; ++a)
            for (std::size_t x = 0; x < nIC; ++x)
                for (std::size_t b = 0; b < I; ++b)
                    for (std::size_t y = 0; y < nIC; ++y)
                        full(a * nIC + x, b * nIC + y) = route_eff(ic[x], ic[y], a + 1, b + 1);
        return stoch_comp_stateful(full, nIC);
    }

    /**
     * The stochastic complement of a NODE-level routing block over the stateful
     * nodes, S = P11 + P12 (I - P22)^-1 P21.
     *
     * Split out of `station_routing` because the STATE-DEPENDENT routing table
     * (`rt_state`, state.h) is the same complement of a block whose SDR rows
     * have been re-evaluated at one state: the two must eliminate the stateless
     * nodes identically, or the per-state table and `rt` would disagree on a
     * model that merely has a Router in it.
     *
     * @param full node-major (nnodes * nIC) square block
     * @param nIC  number of classes carried per node in that block
     */
    Matrix<T> stoch_comp_stateful(const Matrix<T>& full, std::size_t nIC) const {
        const T zero = num_traits<T>::from_int(0);
        const T one = num_traits<T>::from_int(1);
        const std::size_t I = nodes.size();
        std::vector<std::size_t> keep, drop;
        for (std::size_t a = 0; a < I; ++a) {
            const bool st = nodes[a].stateful;
            for (std::size_t x = 0; x < nIC; ++x) (st ? keep : drop).push_back(a * nIC + x);
        }
        if (drop.empty()) return full;

        const std::size_t nk = keep.size(), nd = drop.size();
        Matrix<T> P11(nk, nk, zero), P12(nk, nd, zero), P21(nd, nk, zero), A(nd, nd, zero);
        for (std::size_t a = 0; a < nk; ++a) {
            for (std::size_t b = 0; b < nk; ++b) P11(a, b) = full(keep[a], keep[b]);
            for (std::size_t b = 0; b < nd; ++b) P12(a, b) = full(keep[a], drop[b]);
        }
        for (std::size_t a = 0; a < nd; ++a) {
            for (std::size_t b = 0; b < nk; ++b) P21(a, b) = full(drop[a], keep[b]);
            for (std::size_t b = 0; b < nd; ++b)
                A(a, b) = T((a == b ? one : zero) - full(drop[a], drop[b]));
        }
        // X = A^-1 P21 by Gaussian elimination with partial pivoting
        Matrix<T> X = P21;
        std::vector<std::size_t> piv(nd);
        for (std::size_t i = 0; i < nd; ++i) piv[i] = i;
        for (std::size_t col = 0; col < nd; ++col) {
            std::size_t best = col;
            double bv = std::fabs(num_traits<T>::to_double(A(col, col)));
            for (std::size_t r2 = col + 1; r2 < nd; ++r2) {
                const double v = std::fabs(num_traits<T>::to_double(A(r2, col)));
                if (v > bv) { bv = v; best = r2; }
            }
            if (best != col) {
                for (std::size_t b = 0; b < nd; ++b) std::swap(A(col, b), A(best, b));
                for (std::size_t b = 0; b < nk; ++b) std::swap(X(col, b), X(best, b));
            }
            if (A(col, col) == zero)
                throw NumericError("network: the stochastic complement is singular");
            for (std::size_t r2 = 0; r2 < nd; ++r2) {
                if (r2 == col) continue;
                const T f = T(A(r2, col) / A(col, col));
                if (f == zero) continue;
                for (std::size_t b = 0; b < nd; ++b) A(r2, b) = T(A(r2, b) - f * A(col, b));
                for (std::size_t b = 0; b < nk; ++b) X(r2, b) = T(X(r2, b) - f * X(col, b));
            }
        }
        for (std::size_t r2 = 0; r2 < nd; ++r2)
            for (std::size_t b = 0; b < nk; ++b) X(r2, b) = T(X(r2, b) / A(r2, r2));

        Matrix<T> S = P11;
        for (std::size_t a = 0; a < nk; ++a)
            for (std::size_t b = 0; b < nk; ++b) {
                T acc = zero;
                for (std::size_t m = 0; m < nd; ++m) acc += P12(a, m) * X(m, b);
                S(a, b) = T(S(a, b) + acc);
            }
        return S;
    }

    // ---- predicates, ports of the sn_has_* family -------------------------

    bool has_open_classes() const {
        for (const JobClass& c : classes)
            if (std::isinf(c.population)) return true;
        return false;
    }

    /**
     * `sn_is_open_model`: EVERY class is open, which is not `has_open_classes`.
     * A mixed model passes that predicate and fails this one, and an analyzer
     * that confuses the two hands a closed chain to a solver with no level for
     * its population. An empty class list is not an open model either.
     */
    bool is_open_model() const {
        for (const JobClass& c : classes)
            if (!std::isinf(c.population)) return false;
        return !classes.empty();
    }

    bool has_multi_server() const {
        for (const Station<T>& s : stations)
            if (std::isfinite(s.nservers) && s.nservers > 1.0) return true;
        return false;
    }

    bool has_fractional_populations() const {
        for (const JobClass& c : classes)
            if (std::isfinite(c.population) && c.population != std::floor(c.population + 0.5))
                return true;
        return false;
    }

    bool has_class_switching() const { return nclasses != nchains; }

    bool has_priorities() const {
        for (const JobClass& c : classes)
            if (c.prio > 0) return true;
        return false;
    }

    /**
     * Whether the classes carry more than one priority level.
     *
     * NOT `has_priorities()`, which asks whether any priority is nonzero: a
     * model whose classes all sit at level 1 has priorities by that test and
     * nothing to distinguish, and the reference's warning below keys on the
     * distinction rather than on the magnitude.
     */
    bool has_distinct_priorities() const {
        if (classes.empty()) return false;
        for (const JobClass& c : classes)
            if (c.prio != classes.front().prio) return true;
        return false;
    }

    /**
     * Whether some station runs a discipline that READS the class priorities.
     *
     * The list is exactly the `*PRIO` family. Priority-awareness is a property
     * of the DECLARED policy and is never inferred from the data: a base policy
     * is not upgraded because the classes it was handed carry unequal
     * priorities. `afterEventStation` did infer it once, and plain LCFS with
     * distinct priorities then behaved as none of the three policies involved.
     */
    bool sched_has_priority_aware() const {
        for (const Station<T>& s : stations) {
            switch (s.sched) {
                case SchedStrategy::HOL:
                case SchedStrategy::PSPRIO:
                case SchedStrategy::DPSPRIO:
                case SchedStrategy::GPSPRIO:
                case SchedStrategy::LCFSPRIO:
                case SchedStrategy::LCFSPRPRIO:
                case SchedStrategy::LCFSPIPRIO:
                case SchedStrategy::FCFSPRPRIO:
                case SchedStrategy::FCFSPIPRIO:
                case SchedStrategy::SRPTPRIO:
                    return true;
                default:
                    break;
            }
        }
        return false;
    }

    /** Priorities were declared and no station will read them. */
    bool priorities_ignored() const {
        return has_distinct_priorities() && !sched_has_priority_aware();
    }

    /**
     * Port of sn_has_homogeneous_scheduling.
     *
     * The MATLAB function is `length(findstring(sn.sched, strategy)) ==
     * sn.nstations`, and findstring matches STRINGS: on the numeric sched
     * vector its strcmp is false, so it returns the sentinel -1, whose length
     * is 1. The predicate therefore reduces to nstations == 1 whatever the
     * disciplines are, which is what the reference actually computes and what
     * the AMVA dispatch in solver_amva actually sees. Reproduced rather than
     * corrected: fixing it here would send homogeneous-delay layers down a
     * different branch than every other codebase takes.
     */
    bool has_homogeneous_scheduling(SchedStrategy) const { return nstations == 1; }

    bool has_multi_class_heter_fcfs() const {
        bool bad = false;
        for (std::size_t i = 0; i < nstations; ++i) {
            if (stations[i].sched != SchedStrategy::FCFS) continue;
            bool any = false;
            T lo = num_traits<T>::from_int(0), hi = num_traits<T>::from_int(0);
            for (std::size_t r = 0; r < nclasses; ++r) {
                if (disabled[i][r]) continue;  // MATLAB drops the NaN entries here
                if (!any) {
                    lo = hi = rates(i, r);
                    any = true;
                } else {
                    if (rates(i, r) < lo) lo = rates(i, r);
                    if (rates(i, r) > hi) hi = rates(i, r);
                }
            }
            if (any && hi > lo) bad = true;
        }
        return bad;
    }

    bool sched_is_product_form() const {
        for (const Station<T>& s : stations)
            if (!(s.sched == SchedStrategy::INF || s.sched == SchedStrategy::PS ||
                  s.sched == SchedStrategy::FCFS || s.sched == SchedStrategy::LCFSPR ||
                  s.sched == SchedStrategy::LCFS || s.sched == SchedStrategy::EXT))
                return false;
        return true;
    }

    /**
     * BCMP type 1 asks the FCFS service to be exponential. `has_multi_class_heter_fcfs`
     * compares the class MEANS only, so a class-homogeneous Erlang, hyper-exponential or
     * deterministic FCFS station used to pass this gate and be dispatched to exact MVA,
     * which reads the means alone and returns the exponential answer with no warning.
     */
    bool has_exponential_fcfs() const {
        for (std::size_t i = 0; i < nstations; ++i) {
            if (stations[i].sched != SchedStrategy::FCFS) continue;
            for (std::size_t r = 0; r < nclasses; ++r) {
                const double v = num_traits<T>::to_double(scv(i, r));
                if (std::isinf(v) || !(v > 0.0)) continue;
                if (!(v > 1.0 - GlobalConstants::FineTol && v < 1.0 + GlobalConstants::FineTol))
                    return false;
            }
        }
        return true;
    }

    /**
     * Kendall's K of station IST (1-based), +inf when unbounded. Implementation of
     * api/sn/sn_get_buffer_size.h, which delegates here so the solvers' member calls and
     * the api free function cannot drift; the traps are documented on that header.
     */
    double buffer_size(std::size_t ist) const {
        double k = std::numeric_limits<double>::infinity();
        if (cap[ist - 1] >= 0.0) k = std::min(k, cap[ist - 1]);
        double ccap = 0.0;
        bool anyServed = false;
        for (std::size_t r = 0; r < nclasses; ++r)
            if (classcap[ist - 1][r] > 0.0) {
                ccap += classcap[ist - 1][r];
                anyServed = true;
            }
        if (anyServed) k = std::min(k, ccap);
        double reachable = 0.0;
        for (std::size_t r = 0; r < nclasses; ++r)
            if (!anyServed || classcap[ist - 1][r] > 0.0) reachable += classes[r].population;
        return (k >= reachable) ? std::numeric_limits<double>::infinity() : k;
    }

    /** Implementation of api::sn_is_mm1k_loss; that free function delegates here. */
    bool is_mm1k_loss() const {
        if (nclasses != 1 || nodes.size() != 3) return false;
        for (std::size_t k = 0; k < classes.size(); ++k)
            if (std::isfinite(classes[k].population)) return false;  // nclosedjobs ~= 0
        std::size_t nq = 0, nsrc = 0, nsink = 0, qnode = 0, snode = 0;
        for (std::size_t a = 0; a < nodes.size(); ++a) {
            switch (nodes[a].nodetype) {
                case NodeType::Queue: ++nq; qnode = a + 1; break;
                case NodeType::Source: ++nsrc; snode = a + 1; break;
                case NodeType::Sink: ++nsink; break;
                default: return false;
            }
        }
        if (nq != 1 || nsrc != 1 || nsink != 1) return false;
        const std::size_t qist = nodes[qnode - 1].station;
        const std::size_t sist = nodes[snode - 1].station;
        if (qist == 0 || sist == 0) return false;
        if (stations[qist - 1].nservers != 1.0) return false;
        if (droprule.size() < qist || droprule[qist - 1].empty() ||
            droprule[qist - 1][0] != DropStrategy::DROP)
            return false;
        if (!(qist <= cap.size()) || !std::isfinite(cap[qist - 1]) || !(cap[qist - 1] > 0.0))
            return false;
        if (std::fabs(num_traits<T>::to_double(scv(sist - 1, 0)) - 1.0) > 1e-6) return false;
        if (std::fabs(num_traits<T>::to_double(scv(qist - 1, 0)) - 1.0) > 1e-6) return false;
        return true;
    }

    /**
     * Some station can REFUSE a job: its own buffer BINDS, or a finite capacity region
     * caps a set of stations jointly. Implementation of api::sn_has_blocking, which
     * delegates here; the rule and its two exemptions are documented on that function.
     */
    bool has_blocking() const {
        if (!regions.empty()) return true;
        for (std::size_t a = 0; a < nodes.size(); ++a)
            if (nodes[a].nodetype == NodeType::Cache) return false;
        if (is_mm1k_loss()) return false;
        for (std::size_t ist = 1; ist <= stations.size(); ++ist)
            if (std::isfinite(buffer_size(ist))) return true;
        return false;
    }

    bool has_product_form() const {
        // BCMP asks for infinite buffers: without this conjunct a BAS-blocked station or
        // any binding finite buffer read as product form, though its truncation couples
        // the station occupancies.
        return sched_is_product_form() && !has_multi_class_heter_fcfs() && !has_priorities() &&
               !has_blocking() && has_exponential_fcfs();
    }

    /**
     * Port of sn_has_product_form_not_het_fcfs: LCFS is excluded, and at FCFS the service
     * must be exponential AND class-independent, which is what BCMP type 1 asks for. With
     * unequal per-class means the product-form solve returns a wait proportional to each
     * class's own demand where FCFS makes every class wait behind the same queue. The mean
     * comparison is between CHAIN service times (visit-weighted over the classes that
     * actually visit the station): a class that never visits cannot break product form,
     * and within-chain heterogeneity is invisible to both the product-form and the qd
     * branch, which deaggregate a chain result proportionally to each class's own demand,
     * so only between-chain heterogeneity warrants the divert. LN layers carry seeded
     * rates for classes with zero visits, which a raw per-class comparison mistakes for
     * heterogeneity.
     *
     * CHECK_MEANS drops the mean test; pass false only for an algorithm
     * that models class-dependent FCFS itself (ab, schmidt, schmidt-ext).
     */
    bool has_product_form_not_het_fcfs(bool check_means = true) const {
        for (const Station<T>& s : stations)
            if (!(s.sched == SchedStrategy::INF || s.sched == SchedStrategy::PS ||
                  s.sched == SchedStrategy::FCFS || s.sched == SchedStrategy::LCFSPR ||
                  s.sched == SchedStrategy::EXT))
                return false;
        if (has_priorities()) return false;
        for (std::size_t i = 0; i < nstations; ++i) {
            if (stations[i].sched != SchedStrategy::FCFS) continue;
            for (std::size_t r = 0; r < nclasses; ++r) {
                if (disabled[i][r]) continue;
                const double v = num_traits<T>::to_double(scv(i, r));
                if (!std::isinf(v) && v > 0.0 &&
                    !(v > 1.0 - GlobalConstants::FineTol && v < 1.0 + GlobalConstants::FineTol))
                    return false;
            }
            if (!check_means || visits.empty()) continue;
            const std::size_t isf = stateful_of_station(i + 1) - 1;
            double stmin = 0.0, stmax = 0.0;
            bool anyserved = false;
            for (std::size_t c = 0; c < nchains && c < visits.size(); ++c) {
                double num = 0.0, den = 0.0;
                for (std::size_t r = 0; r < nclasses; ++r) {
                    if (!chains.empty() && !chains[c][r]) continue;
                    if (disabled[i][r]) continue;
                    const double w = num_traits<T>::to_double(visits[c](isf, r));
                    const double rate = num_traits<T>::to_double(rates(i, r));
                    if (w > GlobalConstants::Zero && std::isfinite(rate) && rate > 0.0) {
                        num += w / rate;
                        den += w;
                    }
                }
                if (den > 0.0) {
                    const double st = num / den;
                    if (!anyserved) {
                        stmin = stmax = st;
                        anyserved = true;
                    } else {
                        stmin = std::min(stmin, st);
                        stmax = std::max(stmax, st);
                    }
                }
            }
            if (anyserved && stmax - stmin > GlobalConstants::CoarseTol * stmax) return false;
        }
        return true;
    }

    /** Grow every routing block to the current node count. */
    void grow_routing() {
        for (auto& kv : P) grow_block(kv.second);
    }
    void grow_block(Matrix<T>& B) const {
        if (B.rows() == nodes.size()) return;
        Matrix<T> g(nodes.size(), nodes.size(), num_traits<T>::from_int(0));
        for (std::size_t a = 0; a < B.rows(); ++a)
            for (std::size_t b = 0; b < B.cols(); ++b) g(a, b) = B(a, b);
        B = g;
    }

    /** Total population, as MATLAB's getNumberOfJobs summed. */
    double total_jobs() const {
        double s = 0.0;
        for (const JobClass& c : classes)
            if (std::isfinite(c.population)) s += c.population;
        return s;
    }
};

/**
 * The swap graph of a PAS / OI station, with the defaults `refreshLocalVars.m`
 * installs applied.
 *
 * MATLAB fills `sn.nodeparam{ind}.swapGraph` at refresh time and NEVER leaves it
 * empty: an OI station always gets `zeros(R,R)`, a PAS station with no explicit
 * graph gets the complete compatibility graph `ones(R,R) - eye(R)`, and an
 * explicit graph is taken as given. The builder here stores the raw user graph,
 * so the defaulting has to happen on read; doing it here rather than in the
 * builder keeps it independent of whether the classes were added before or
 * after `set_service_rate_function`.
 *
 * @param sn  the struct
 * @param ist 1-based station index
 */
template <class T>
Matrix<T> station_swap_graph(const NetworkStruct<T>& sn, std::size_t ist) {
    const Station<T>& st = sn.stations[ist - 1];
    const std::size_t R = sn.nclasses;
    const T zero = num_traits<T>::from_int(0);
    if (st.sched == SchedStrategy::OI) return Matrix<T>(R, R, zero);
    if (!st.swap_graph.empty()) return st.swap_graph;
    if (st.sched != SchedStrategy::PAS) return Matrix<T>(R, R, zero);
    Matrix<T> g(R, R, num_traits<T>::from_int(1));
    for (std::size_t r = 0; r < R; ++r) g(r, r) = zero;
    return g;
}

/** True when the station's materialized swap graph is entirely zero. */
template <class T>
bool station_swap_graph_is_zero(const NetworkStruct<T>& sn, std::size_t ist) {
    const Matrix<T> g = station_swap_graph(sn, ist);
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t a = 0; a < g.rows(); ++a)
        for (std::size_t b = 0; b < g.cols(); ++b)
            if (g(a, b) != zero) return false;
    return true;
}

}  // namespace qn
}  // namespace line

#endif  // LINE_LANG_QN_NETWORK_STRUCT_H
