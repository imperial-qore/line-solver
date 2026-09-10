/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_LQN_LQN_BUILDER_H
#define LINE_LANG_LQN_LQN_BUILDER_H

/**
 * Build a layered queueing network in code, as the MATLAB constructors do.
 *
 * Port of the Processor / Task / Entry / Activity / ActivityPrecedence API in
 * matlab/src/lang/layered/, feeding the same `lqn_finalize` (getStruct) that
 * the .lqnx reader feeds.
 *
 * WHY THIS EXISTS rather than always going through a file. The .lqnx
 * interchange is LOSSY for models a script can express. The clearest case is a
 * think time on a non-reference task: `writeXML` omits it because lqns rejects
 * the attribute there (see _kb, "lqnx cannot carry non-ref think time"), so
 * exporting such a model and reading it back silently drops the think time and
 * changes the answer. matlab/examples/basic/layeredModel/lqn_workflows.m is
 * exactly that shape -- its T3 has both an infinite-server discipline and a
 * think time of 10 -- so it cannot be reached through the file format at all.
 *
 * ORDER MATTERS, as in the reader: element indices are assigned in declaration
 * order (all processors, then all tasks, then all entries, then all
 * activities), and within a task the activities are ordered by declaration.
 * Declare in the same order as the reference script and the indices agree.
 *
 * The builder validates references by name and throws on an unknown one; it
 * does not silently create elements, because a typo that creates a second
 * disconnected task is a model that still solves.
 */

#include <string>
#include <vector>

#include "line/lang/lqn/lqn_reader.h"
#include "line/util/error.h"

namespace line {
namespace lqn {

template <class T>
class LqnBuilder {
public:
    /** Add a processor. `mult` may be infinite; INF scheduling forces it so. */
    std::size_t processor(const std::string& name, double mult, SchedStrategy sched,
                          double repl = 1.0) {
        detail::RawProc p;
        p.name = name;
        p.sched = sched;
        p.repl = repl;
        p.mult = sched == SchedStrategy::INF ? std::numeric_limits<double>::infinity() : mult;
        m_.procs.push_back(p);
        return m_.procs.size() - 1;
    }

    /** Add a task on a processor. */
    std::size_t task(const std::string& name, double mult, SchedStrategy sched,
                     const std::string& on_processor, double repl = 1.0) {
        detail::RawTask<T> t;
        t.name = name;
        t.sched = sched;
        t.repl = repl;
        t.mult = sched == SchedStrategy::INF ? std::numeric_limits<double>::infinity() : mult;
        t.thinktime = Distrib<T>::immediate();
        t.proc_slot = find_proc(on_processor);
        m_.tasks.push_back(t);
        return m_.tasks.size() - 1;
    }

    /**
     * Set a task's think time.
     *
     * Accepted on ANY task, not only a reference one. That is what the MATLAB
     * API allows and what the .lqnx writer cannot express; refusing it here to
     * match the file format would make the builder strictly weaker than the
     * reference for no gain.
     */
    void think_time(const std::string& task_name, const Distrib<T>& d) {
        m_.tasks[find_task(task_name)].thinktime = d;
    }

    /**
     * A CacheTask: a task whose entries are looked up in a cache of `nitems`.
     *
     * `itemcap` is the capacity of each cache list, so a plain single-level
     * cache passes one value. The task itself is an ordinary server in its own
     * layer; the Cache NODE appears in its HOST's layer, which is where
     * buildLayersRecursive puts it (`iscachelayer` is a host-layer test).
     */
    std::size_t cache_task(const std::string& name, double mult, SchedStrategy sched,
                           const std::string& on_processor, std::size_t nitems,
                           const std::vector<int>& itemcap, ReplacementStrategy replacestrat,
                           double repl = 1.0) {
        if (nitems == 0)
            throw InputError("LqnBuilder::cache_task: '" + name + "' caches no item");
        if (itemcap.empty())
            throw InputError("LqnBuilder::cache_task: '" + name + "' has no cache list");
        int total = 0;
        for (std::size_t i = 0; i < itemcap.size(); ++i) total += itemcap[i];
        if (total <= 0 || static_cast<std::size_t>(total) >= nitems)
            throw InputError("LqnBuilder::cache_task: '" + name +
                             "' has a total capacity that is not between 1 and nitems-1; "
                             "a cache that holds every item never misses");
        const std::size_t t = task(name, mult, sched, on_processor, repl);
        m_.tasks[t].nitems = nitems;
        m_.tasks[t].itemcap = itemcap;
        m_.tasks[t].replacestrat = replacestrat;
        return t;
    }

    /**
     * A SetupTask: a server that powers down when idle and pays to restart.
     *
     * `setup` is charged to the first arrival that finds the server off;
     * `delayoff` is the idle timer that has to expire before it goes off, so a
     * long delay-off makes the setup rare. Both are ordinary Task properties in
     * the reference (`Task.setSetupTime`/`setDelayOffTime`), so this is a plain
     * task with the two times attached rather than a distinct kind.
     */
    void setup_time(const std::string& task_name, const Distrib<T>& setup,
                    const Distrib<T>& delayoff) {
        if (setup.disabled)
            throw InputError("LqnBuilder::setup_time: '" + task_name + "' has no setup time");
        if (delayoff.disabled)
            throw InputError(
                "LqnBuilder::setup_time: '" + task_name +
                "' has a setup time but no delay-off time; a server that never shuts down pays "
                "the setup at most once and the reference declines to model it");
        detail::RawTask<T>& t = m_.tasks[find_task(task_name)];
        t.setuptime = setup;
        t.delayofftime = delayoff;
    }

    /**
     * An ItemEntry: the entry a cache read enters, over `cardinality` items.
     *
     * `popularity` is the pmf of the read over those items. The reference takes
     * a discrete Distribution (a Zipf, typically) and reads only its pmf; this
     * port has no discrete-distribution type, so the pmf is given directly.
     */
    std::size_t item_entry(const std::string& name, const std::string& on_task,
                           std::size_t cardinality, const std::vector<T>& popularity) {
        if (cardinality == 0)
            throw InputError("LqnBuilder::item_entry: '" + name + "' indexes no item");
        if (popularity.size() != cardinality)
            throw InputError("LqnBuilder::item_entry: '" + name +
                             "' needs one popularity per item");
        const std::size_t e = entry(name, on_task);
        m_.entries[e].cardinality = cardinality;
        m_.entries[e].popularity = popularity;
        return e;
    }

    /** Add an entry on a task. */
    std::size_t entry(const std::string& name, const std::string& on_task) {
        detail::RawEntry<T> e;
        e.name = name;
        e.task_slot = find_task(on_task);
        m_.entries.push_back(e);
        return m_.entries.size() - 1;
    }

    /** An open arrival stream at an entry. */
    void open_arrival(const std::string& entry_name, const Distrib<T>& d) {
        detail::RawEntry<T>& e = m_.entries[find_entry(entry_name)];
        e.has_arrival = true;
        e.arrival = d;
    }

    /** Add an activity on a task, with its host demand. */
    std::size_t activity(const std::string& name, const Distrib<T>& hostdem,
                         const std::string& on_task) {
        detail::RawActivity<T> a;
        a.name = name;
        a.hostdem = hostdem;
        a.thinktime = Distrib<T>::immediate();
        a.task_slot = find_task(on_task);
        a.phase = 1;
        m_.acts.push_back(a);
        return m_.acts.size() - 1;
    }

    /** Bind an activity to an entry: it is the entry's first activity. */
    void bound_to(const std::string& act, const std::string& entry_name) {
        m_.acts[find_act(act)].bound_to_entry = entry_name;
    }

    /** A synchronous call from an activity to an entry of another task. */
    void sync_call(const std::string& act, const std::string& dest_entry, const T& mean) {
        m_.acts[find_act(act)].sync_calls.push_back({dest_entry, mean});
    }

    /** An asynchronous call from an activity to an entry of another task. */
    void async_call(const std::string& act, const std::string& dest_entry, const T& mean) {
        m_.acts[find_act(act)].async_calls.push_back({dest_entry, mean});
    }

    /**
     * ONE synchronous call per invocation, its destination CYCLING over the
     * targets in the order given (`synchCallRoundRobin`).
     *
     * The members are ordinary sync calls of mean `mean/n`, so the aggregate
     * call rate is `mean` either way; what round robin removes is the variance
     * of the branching, which is what smooths the target queues. Needs at least
     * two targets -- a group of one is not a dispatch decision.
     */
    void sync_call_round_robin(const std::string& act,
                               const std::vector<std::string>& dest_entries, const T& mean) {
        add_call_group(act, dest_entries, mean, lang::RoutingStrategy::RROBIN,
                       "sync_call_round_robin");
    }

    /** As above, with the least loaded target taking the call (`synchCallJSQ`). */
    void sync_call_jsq(const std::string& act, const std::vector<std::string>& dest_entries,
                       const T& mean) {
        add_call_group(act, dest_entries, mean, lang::RoutingStrategy::JSQ, "sync_call_jsq");
    }

    /**
     * Forwarding: whenever `src_entry` is invoked, with probability `prob` the
     * request is handed onward to `dest_entry` instead of `src_entry` replying
     * -- lqn_finalize (shared with the .lqnx reader, see lqn_reader.h) turns
     * this into a CallType::FWD call and SolverLN's lqn_fwd_rendezvous rewrite
     * (lqn_helpers.h) flattens it into a caller-side pseudo rendezvous.
     */
    void forward(const std::string& src_entry, const std::string& dest_entry, const T& prob) {
        detail::RawEntry<T>& e = m_.entries[find_entry(src_entry)];
        e.fwd_dest.push_back(dest_entry);
        e.fwd_prob.push_back(prob);
    }

    /** Mark an activity as the one that replies to an entry. */
    void replies_to(const std::string& act, const std::string& entry_name) {
        m_.entries[find_entry(entry_name)].reply_activities.push_back(act);
    }

    /** An activity think time, in series with the host demand. */
    void act_think_time(const std::string& act, const Distrib<T>& d) {
        m_.acts[find_act(act)].thinktime = d;
    }

    // ---- precedences ------------------------------------------------------

    /** pre -> post, a plain sequence. */
    void serial(const std::string& pre, const std::string& post) {
        detail::RawPrecedence<T> p;
        p.pretype = PrecedenceType::PRE_SEQ;
        p.posttype = PrecedenceType::POST_SEQ;
        p.preacts.push_back(pre);
        p.postacts.push_back(post);
        add_prec(pre, p);
    }

    /** pre -> every post, concurrently. */
    void and_fork(const std::string& pre, const std::vector<std::string>& posts) {
        detail::RawPrecedence<T> p;
        p.pretype = PrecedenceType::PRE_SEQ;
        p.posttype = PrecedenceType::POST_AND;
        p.preacts.push_back(pre);
        p.postacts = posts;
        add_prec(pre, p);
    }

    /** all pres (or `quorum` of them) -> post. */
    void and_join(const std::vector<std::string>& pres, const std::string& post,
                  std::size_t quorum = 0) {
        detail::RawPrecedence<T> p;
        p.pretype = PrecedenceType::PRE_AND;
        p.posttype = PrecedenceType::POST_SEQ;
        p.preacts = pres;
        p.postacts.push_back(post);
        if (quorum > 0) {
            p.has_quorum = true;
            p.quorum = quorum;
        }
        add_prec(pres.at(0), p);
    }

    /** pre -> one of the posts, with the given branch probabilities. */
    void or_fork(const std::string& pre, const std::vector<std::string>& posts,
                 const std::vector<T>& probs) {
        if (posts.size() != probs.size())
            throw InputError("LqnBuilder::or_fork: one probability per branch is required");
        detail::RawPrecedence<T> p;
        p.pretype = PrecedenceType::PRE_SEQ;
        p.posttype = PrecedenceType::POST_OR;
        p.preacts.push_back(pre);
        p.postacts = posts;
        p.postparams = probs;
        add_prec(pre, p);
    }

    /**
     * `ActivityPrecedence.CacheAccess(pre, {hit, miss})`.
     *
     * `pre` reads the cache; the job leaves it on the HIT branch or the MISS
     * branch. POST_CACHE lands on the two successors, not on `pre` -- that is
     * how getStruct.m records a post type, and what SolverLN keys on.
     */
    void cache_access(const std::string& pre, const std::string& hit, const std::string& miss) {
        detail::RawPrecedence<T> p;
        p.pretype = PrecedenceType::PRE_SEQ;
        p.posttype = PrecedenceType::POST_CACHE;
        p.preacts.push_back(pre);
        p.postacts.push_back(hit);
        p.postacts.push_back(miss);
        add_prec(pre, p);
    }

    /** any of the pres -> post. */
    void or_join(const std::vector<std::string>& pres, const std::string& post) {
        detail::RawPrecedence<T> p;
        p.pretype = PrecedenceType::PRE_OR;
        p.posttype = PrecedenceType::POST_SEQ;
        p.preacts = pres;
        p.postacts.push_back(post);
        // PRE_OR carries a probability per branch in MATLAB; a plain or-join
        // leaves them unset, which getStruct reads as absent
        add_prec(pres.at(0), p);
    }

    /**
     * pre -> body, repeated `count` times in expectation, then -> end.
     *
     * The body list is the loop body in order; `end` is the activity the loop
     * exits to. This is MATLAB's ActivityPrecedence.Loop(pre, body, count),
     * whose postacts vector is body followed by end.
     */
    void loop(const std::string& pre, const std::vector<std::string>& body,
              const std::string& end, const T& count) {
        detail::RawPrecedence<T> p;
        p.pretype = PrecedenceType::PRE_SEQ;
        p.posttype = PrecedenceType::POST_LOOP;
        p.preacts.push_back(pre);
        p.postacts = body;
        p.postacts.push_back(end);
        p.postparams.assign(body.size(), count);
        add_prec(pre, p);
    }

    // ---- admission constraints --------------------------------------------

    /**
     * `elem.addConstraint(operands, coeffs, cap)`: sum(coeffs .* n(operands)) <= cap.
     *
     * `elem` is a task (whose operands are its entries) or a host (whose
     * operands are its tasks). Port of LayeredNetworkElement.addConstraint;
     * the operands stay NAMED until build(), because their column order is the
     * task's entriesof / the host's tasksof and neither exists yet.
     */
    void add_constraint(const std::string& elem, const std::vector<std::string>& operands,
                        const std::vector<T>& coeffs, const T& cap) {
        if (operands.size() != coeffs.size())
            throw InputError("LqnBuilder::add_constraint: one coefficient per operand is required");
        if (operands.empty())
            throw InputError("LqnBuilder::add_constraint: the constraint names no operand");
        for (std::size_t i = 0; i < operands.size(); ++i)
            for (std::size_t j = i + 1; j < operands.size(); ++j)
                if (operands[i] == operands[j])
                    throw InputError("LqnBuilder::add_constraint: operand '" + operands[i] +
                                     "' appears twice");
        detail::RawLinConRow<T> row;
        row.names = operands;
        row.coeffs = coeffs;
        row.cap = cap;
        linconrows_of(elem).push_back(row);
    }

    /**
     * The positional form, `elem.setConstraint(A, b)`.
     *
     * Only the column COUNT is checked, at build() time, as the reference does:
     * the columns are the element's entries or tasks in declaration order.
     */
    void set_constraint(const std::string& elem, const Matrix<T>& A, const std::vector<T>& b) {
        if (A.rows() != b.size())
            throw InputError("LqnBuilder::set_constraint: A and b disagree on the number of rows");
        for (std::size_t i = 0; i < m_.tasks.size(); ++i)
            if (m_.tasks[i].name == elem) {
                m_.tasks[i].lincon_A = A;
                m_.tasks[i].lincon_b = b;
                return;
            }
        for (std::size_t i = 0; i < m_.procs.size(); ++i)
            if (m_.procs[i].name == elem) {
                m_.proc_lincon[i] = std::make_pair(A, b);
                return;
            }
        throw InputError("LqnBuilder::set_constraint: unknown task or processor '" + elem + "'");
    }

    /**
     * `elem.setLoadDependence(alpha)`: alpha(n) scales the rate of the layer
     * station of ELEM when it holds n jobs in total, on top of the multiplicity.
     */
    void set_load_dependence(const std::string& elem, const std::vector<T>& alpha) {
        assert_rate_dependent(elem, "Load");
        if (alpha.empty())
            throw InputError("LqnBuilder::set_load_dependence: alpha is empty");
        if (is_task(elem)) {
            m_.tasks[find_task(elem)].lldscaling = alpha;
        } else {
            m_.proc_lldscaling[find_proc(elem)] = alpha;
        }
    }

    /**
     * `elem.setClassDependence(beta, peak)`: the product-form handle, whose
     * argument is the per-OPERAND population of ELEM -- task j of a processor,
     * entry j of a task, in declaration order.
     */
    void set_class_dependence(const std::string& elem, const CdScaling<T>& beta,
                              const std::vector<T>& peak) {
        assert_rate_dependent(elem, "Class");
        assert_dependence_handle(beta, peak, "Class");
        if (is_task(elem)) {
            const std::size_t t = find_task(elem);
            m_.tasks[t].cdscaling = beta;
            m_.tasks[t].cdscalingpeak = peak;
        } else {
            const std::size_t p = find_proc(elem);
            m_.proc_cdscaling[p] = beta;
            m_.proc_cdscalingpeak[p] = peak;
        }
    }

    /**
     * `elem.setJointDependence(eta, peak)`: the non-product-form handle, read at
     * the whole per-operand vector, so solvers treat it as an approximation.
     */
    void set_joint_dependence(const std::string& elem, const CdScaling<T>& eta,
                              const std::vector<T>& peak) {
        assert_rate_dependent(elem, "Joint");
        assert_dependence_handle(eta, peak, "Joint");
        assert_no_pools(elem, "a joint dependence");
        if (is_task(elem)) {
            const std::size_t t = find_task(elem);
            m_.tasks[t].jdscaling = eta;
            m_.tasks[t].jdscalingpeak = peak;
        } else {
            const std::size_t p = find_proc(elem);
            m_.proc_jdscaling[p] = eta;
            m_.proc_jdscalingpeak[p] = peak;
        }
    }

    /**
     * `elem.addServerType(ServerType(pool, count, compatible))`: one pool of
     * COUNT identical servers, each running at RATE, eligible for the operands
     * named in COMPATIBLE.
     *
     * The operands are resolved by name against the element's own operand list
     * at build() time, so a pool may name a task or an entry that is declared
     * later. Pools accumulate; SolverLN lowers the whole declaration to the
     * activated-server rate.
     */
    void add_server_type(const std::string& elem, const std::string& pool, double count,
                         const std::vector<std::string>& compatible, const T& rate) {
        assert_rate_dependent(elem, "Compatibility");
        assert_no_jd(elem);
        if (count < 1)
            throw InputError("LqnBuilder::add_server_type: pool '" + pool +
                             "' must hold at least one server");
        if (compatible.empty())
            throw InputError("LqnBuilder::add_server_type: pool '" + pool +
                             "' is compatible with no operand, so it can never serve");
        detail::RawServerPool<T> sp;
        sp.name = pool;
        sp.count = count;
        sp.rate = rate;
        sp.compatible = compatible;
        raw_pools_of(elem).push_back(sp);
    }

    /** Flatten into the struct SolverLN consumes. */
    LqnStruct<T> build() const { return lqn_finalize(m_); }

    const LqnModel<T>& model() const { return m_; }

private:
    LqnModel<T> m_;

    std::size_t find_proc(const std::string& n) const {
        for (std::size_t i = 0; i < m_.procs.size(); ++i)
            if (m_.procs[i].name == n) return i;
        throw InputError("LqnBuilder: unknown processor '" + n + "'");
    }
    std::size_t find_task(const std::string& n) const {
        for (std::size_t i = 0; i < m_.tasks.size(); ++i)
            if (m_.tasks[i].name == n) return i;
        throw InputError("LqnBuilder: unknown task '" + n + "'");
    }
    std::size_t find_entry(const std::string& n) const {
        for (std::size_t i = 0; i < m_.entries.size(); ++i)
            if (m_.entries[i].name == n) return i;
        throw InputError("LqnBuilder: unknown entry '" + n + "'");
    }
    std::size_t find_act(const std::string& n) const {
        for (std::size_t i = 0; i < m_.acts.size(); ++i)
            if (m_.acts[i].name == n) return i;
        throw InputError("LqnBuilder: unknown activity '" + n + "'");
    }
    /** The constraint-row list of a task or a host, by name. */
    std::vector<detail::RawLinConRow<T>>& linconrows_of(const std::string& n) {
        for (std::size_t i = 0; i < m_.tasks.size(); ++i)
            if (m_.tasks[i].name == n) return m_.tasks[i].linconrows;
        for (std::size_t i = 0; i < m_.procs.size(); ++i)
            if (m_.procs[i].name == n) return m_.proc_linconrows[i];
        throw InputError("LqnBuilder: unknown task or processor '" + n + "'");
    }

    bool is_task(const std::string& n) const {
        for (std::size_t i = 0; i < m_.tasks.size(); ++i)
            if (m_.tasks[i].name == n) return true;
        return false;
    }

    /** The declared pool list of a task or a host, by name. */
    std::vector<detail::RawServerPool<T>>& raw_pools_of(const std::string& n) {
        for (std::size_t i = 0; i < m_.tasks.size(); ++i)
            if (m_.tasks[i].name == n) return m_.tasks[i].pools;
        for (std::size_t i = 0; i < m_.procs.size(); ++i)
            if (m_.procs[i].name == n) return m_.proc_pools[i];
        throw InputError("LqnBuilder: unknown task or processor '" + n + "'");
    }

    /**
     * Only a Task or a Host becomes a layer STATION, and only a PS or FCFS one
     * admits a rate scaling. Twin of LayeredNetworkElement.assertRateDependent.
     */
    void assert_rate_dependent(const std::string& elem, const char* what) const {
        SchedStrategy sched = SchedStrategy::NONE;
        bool found = false;
        for (std::size_t i = 0; i < m_.tasks.size() && !found; ++i)
            if (m_.tasks[i].name == elem) {
                sched = m_.tasks[i].sched;
                found = true;
            }
        for (std::size_t i = 0; i < m_.procs.size() && !found; ++i)
            if (m_.procs[i].name == elem) {
                sched = m_.procs[i].sched;
                found = true;
            }
        if (!found)
            throw InputError(std::string(what) +
                             "-dependence can only be set on a Task or a Host, which are the only "
                             "elements that become server stations in a layer; '" +
                             elem + "' is neither");
        if (sched != SchedStrategy::PS && sched != SchedStrategy::FCFS)
            throw InputError(std::string(what) +
                             "-dependence supported only for processor sharing (PS) and "
                             "first-come first-serve (FCFS) servers, but '" +
                             elem + "' is scheduled otherwise");
    }

    /** The peak rate is a model input; without it utilization has no normalizer. */
    void assert_dependence_handle(const CdScaling<T>& f, const std::vector<T>& peak,
                                  const char* what) const {
        if (!f)
            throw InputError(std::string(what) + "-dependence needs a handle");
        if (peak.empty())
            throw InputError(std::string(what) +
                             "-dependence needs a peak rate per operand, which normalizes "
                             "utilization as U = T*S/peak");
    }

    /** One rate law per server: pools and an explicit handle would both claim it. */
    void assert_no_pools(const std::string& elem, const char* what) {
        if (!raw_pools_of(elem).empty())
            throw InputError("LqnBuilder: '" + elem + "' already declares server pools, which are "
                             "themselves a rate law, so it cannot also take " + what);
    }

    void assert_no_jd(const std::string& elem) {
        bool has = false;
        for (std::size_t i = 0; i < m_.tasks.size(); ++i)
            if (m_.tasks[i].name == elem && m_.tasks[i].jdscaling) has = true;
        for (std::size_t i = 0; i < m_.procs.size(); ++i)
            if (m_.procs[i].name == elem && m_.proc_jdscaling.count(i)) has = true;
        if (has)
            throw InputError("LqnBuilder: '" + elem + "' already declares a joint dependence, so "
                             "it cannot also declare server pools, which are a rate law of their "
                             "own");
    }

    /** A precedence belongs to the task owning its activities. */
    void add_prec(const std::string& anchor_act, const detail::RawPrecedence<T>& p) {
        m_.tasks[m_.acts[find_act(anchor_act)].task_slot].precedences.push_back(p);
    }

    /**
     * A routed group: n ordinary sync calls of mean/n, plus the record that
     * they are one dispatch. Splitting the mean is the reference's own
     * `addCallGroup` (Activity.m), and it is what keeps the aggregate call rate
     * equal to the probabilistic twin's.
     */
    void add_call_group(const std::string& act, const std::vector<std::string>& dest_entries,
                        const T& mean, lang::RoutingStrategy rs, const char* who) {
        if (dest_entries.size() < 2)
            throw InputError(std::string("LqnBuilder::") + who +
                             " needs at least two target entries: a group of one is not a "
                             "dispatch decision");
        const T share = T(mean / num_traits<T>::from_int(int(dest_entries.size())));
        for (const std::string& d : dest_entries) sync_call(act, d, share);
        m_.acts[find_act(act)].call_groups.push_back(std::make_pair(rs, dest_entries));
    }
};

}  // namespace lqn
}  // namespace line

#endif  // LINE_LANG_LQN_LQN_BUILDER_H
