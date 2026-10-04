/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_LQN2QN_H
#define LINE_IO_LQN2QN_H

/**
 * @file
 * @ingroup line_io
 * Convert a layered queueing network into a single queueing network in which
 * synchronous call blocking is carried by REPLY signals.
 *
 * Port of matlab/src/io/LQN2QN.m (the reference), mirrored by
 * jar/src/main/java/jline/io/LQN2QN.java and python/line_solver/io/__init__.py.
 * The construction is the reference's, step for step, and so is the order in
 * which nodes and classes are created, so that a model converted here and one
 * converted by MATLAB or Python have the same node and class tables:
 *
 * - One station per host processor replica (a Delay when the processor is an
 *   infinite server), one Delay per reference task replica for its think time.
 * - One class per STEP of the expanded activity graph: an activity, one call
 *   stage of an activity, or a merge/trigger step. Reference-task steps are
 *   closed classes of population 0, except the think class, which carries the
 *   task multiplicity; open-arrival steps are open classes.
 * - A synchronous call site that holds its caller's server gets a REPLY signal
 *   bound to its class (sn.syncreply); the callee's replying step switches into
 *   that signal, which returns to the caller station and releases the server.
 * - A call mean m is unrolled into floor(m) mandatory stages plus one stage
 *   taken with probability m - floor(m), capped at 20 stages.
 * - OR branches and loops follow the graph weights; AND forks become a Fork plus
 *   one Router per branch and AND joins a Join, with a PARTIAL strategy when the
 *   declared quorum is below the branch count.
 * - A CacheTask becomes a Cache node whose read step switches into the hit and
 *   the miss class; delayed-hit retrieval adds a PS fetch station per replica.
 * - Asynchronous calls are non-blocking visits, forwarding splits the reply
 *   exits of the forwarding entry, and the multiplicity of a non-reference task
 *   is a thread pool enforced by one finite capacity region with one linear
 *   admission row per task replica.
 * - Phase-2 activities run after the reply, spawned by the replying step's
 *   completions (sn.classspawn) and destroyed at the chain end by a NEGATIVE
 *   signal (closed chains) or at the Sink (open chains).
 * - Activity think times are steps on a shared ActivityThink delay; setup tasks
 *   carry their setup/delay-off pair onto their host station.
 * - Replication is materialised (one station and one step-graph copy per
 *   replica, calls reaching the fan-out block {(i*f+k) mod r}) or pooled (one
 *   station of r times the servers), selected by 'auto', 'materialize' or
 *   'pool' exactly as in the reference.
 *
 * REPLIES. The reference reads which activity replies to which entry from
 * `lsn.replygraph`, which holds the EXPLICIT replies plus the implicit ones
 * getStruct.m infers for leaf activities. The C++ LqnStruct keeps no replygraph,
 * so the LqnModel overload recovers it with lqn::detail::reply_activities, the
 * port of that same pass. The LqnStruct overload can only infer the implicit
 * replies: an entry whose reply is declared on a non-leaf activity, which is
 * what a phase-2 activity is, needs the LqnModel overload.
 *
 * Construction choices that are the reference's own and not yet represented
 * there either (retrieval on a cache read with phase-2 successors, the thread
 * pool of a task with an internal AND-fork, and the others listed in LQN2QN.m)
 * are reported as warnings, exactly where the reference calls line_warning.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "line/lang/lqn/lqn_reader.h"
#include "line/lang/lqn/lqn_struct.h"
#include "line/lang/lqn/lqn_writer.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace io {

namespace detail {

/** The LQN2QN walk: the MATLAB nested functions as members, their shared state as fields. */
template <class T>
class Lqn2Qn {
  public:
    typedef std::pair<std::size_t, std::size_t> Key;  ///< (element index, 0-based replica)
    typedef lang::Distrib<T> Dist;

    Lqn2Qn(const lqn::LqnStruct<T>& lsn, const std::set<Key>& replies,
           const std::string& replication, const std::string& name,
           std::vector<std::string>* warnings)
        : l_(lsn), replies_(replies), replication_(replication), warnings_(warnings),
          net_(name + "-QN") {}

    qn::Network<T> run() {
        build();
        return net_;
    }

  private:
    static const std::size_t NONE = static_cast<std::size_t>(-1);
    static const std::size_t MAXCALLSTAGES = 20;     // guard against unrolling a huge call mean
    static const std::size_t MAXREPLINSTANCES = 128;  // above this 'auto' pools instead

    /** The port a step is left through; `sig` means through its REPLY signal. */
    struct Port {
        std::size_t step;
        bool sig;
        Port(std::size_t s = NONE, bool g = false) : step(s), sig(g) {}
        bool operator==(const Port& o) const { return step == o.step && sig == o.sig; }
    };
    struct Exit {
        std::size_t step;
        bool sig;
        T prob;
    };
    struct Flow {
        std::size_t from, to;
        T p;
        bool fromSig, inTarget;
    };
    struct Reply {
        std::size_t exit, owner;
        bool viaSig;
        T p;
    };
    struct Ph2Exit {
        std::size_t step;
        bool sig;
        T p;
        Key ref;
    };
    struct CacheWire {
        std::size_t node, readStep, hitStep, missStep, eidx, fetch;
        Dist fetchSvc;
    };
    struct Stage {
        std::size_t target;
        T prob;
        bool isasync;
    };
    struct Step {
        std::size_t aidx;
        Key host;
        Dist svc;  ///< disabled = no service declared by the step
        std::string name;
        bool blocks, isthink;
        Key ref;           ///< (reference task, replica); (0,0) on an open chain
        std::size_t node;  ///< Fork/Join/Router/Cache/ActivityThink node, 0 for stations
        std::size_t owner;
        std::vector<Key> tasks;  ///< thread-pool task replicas holding a thread here
    };
    struct Expansion {
        std::size_t first = NONE;
        std::vector<Exit> replies, terms;
    };
    /** State of one expandActivities call, the closure of the MATLAB walk. */
    struct WalkCtx {
        std::size_t eidx, tidx, trep;
        Key ref;
        std::vector<std::size_t> localActs;
        std::map<std::size_t, std::pair<std::size_t, Port> > visited;
        std::map<std::size_t, std::size_t> joinOf;
        std::map<std::size_t, std::size_t> joinFork;  ///< join activity -> fork step it closes
        std::vector<std::size_t> forkOwnerStack;
        bool sawReply = false;
        std::vector<Exit> replyExits, terminals;
    };

    const lqn::LqnStruct<T>& l_;
    const std::set<Key>& replies_;
    std::string replication_;
    std::vector<std::string>* warnings_;
    qn::Network<T> net_;

    std::vector<double> replRaw_;
    bool materialize_ = false;
    std::vector<bool> fcrTask_;
    std::map<std::size_t, std::size_t> hostNTasks_;  // host element -> number of tasks it runs
    std::vector<std::size_t> refTasks_, openEntries_;
    std::map<Key, std::size_t> hostStation_, thinkNode_, cacheNodeOf_, fetchNodeOf_;
    std::map<Key, bool> hostIsDelay_;
    std::vector<Step> steps_;
    std::vector<Flow> flow_;
    std::vector<Reply> reply_;
    std::vector<std::pair<std::size_t, std::size_t> > spawnPairs_, joinQuorum_;
    std::vector<Ph2Exit> ph2Exits_;
    std::vector<Key> entryStack_, threadStack_;
    std::vector<CacheWire> cacheWiring_;
    std::map<std::string, int> usedNames_;
    std::size_t actThinkNode_ = 0;
    std::size_t srcNode_ = 0, snkNode_ = 0;

    // ------------------------------------------------------------ small helpers

    static T tnum(double v) { return num_traits<T>::from_double(v); }
    static double dbl(const T& v) { return num_traits<T>::to_double(v); }

    void warn(const std::string& msg) {
        if (warnings_ != nullptr)
            warnings_->push_back("LQN2QN: " + msg);
        else
            std::cerr << "[LINE] Warning: LQN2QN: " << msg << std::endl;
    }

    /** A distribution that is declared, not Immediate and above tolerance. */
    static bool timed(const Dist& d) {
        return !d.disabled && !d.is_immediate() &&
               dbl(d.mean) > lang::GlobalConstants::FineTol;
    }

    std::size_t parent(std::size_t idx) const { return l_.parent[idx]; }
    const std::string& nm(std::size_t idx) const { return l_.names[idx]; }

    bool isAndJoinPre(std::size_t aidx) const {
        return aidx < l_.actpretype.size() && l_.actpretype[aidx] == lang::PrecedenceType::PRE_AND;
    }
    bool isPostAnd(std::size_t aidx) const {
        return aidx < l_.actposttype.size() &&
               l_.actposttype[aidx] == lang::PrecedenceType::POST_AND;
    }
    bool isAndFork(const std::vector<std::size_t>& succ) const {
        if (succ.size() < 2) return false;
        for (std::size_t s : succ)
            if (!isPostAnd(s)) return false;
        return true;
    }
    bool repliesTo(std::size_t aidx, std::size_t eidx) const {
        return replies_.count(Key(aidx, eidx)) > 0;
    }
    bool isRouter(std::size_t node) const {
        return node != 0 && net_raw().nodes[node - 1].nodetype == lang::NodeType::Router;
    }
    const qn::NetworkStruct<T>& net_raw() const {
        return const_cast<qn::Network<T>&>(net_).raw_struct();
    }

    std::size_t nrep(std::size_t idx) const {
        return materialize_ ? static_cast<std::size_t>(replRaw_[idx]) : 1;
    }
    double poolFactor(std::size_t idx) const { return materialize_ ? 1.0 : replRaw_[idx]; }

    /** Callee replicas one caller replica reaches; an unset fan-out solves repl(a)*f = repl(b). */
    std::size_t fanOutOf(std::size_t a, std::size_t b, std::size_t rb) const {
        if (rb <= 1) return 1;
        double f = l_.fanout_at(a, b);
        if (f <= 0.0) {
            const double ra = replRaw_[a];
            f = double(rb) > ra ? std::max(1.0, std::floor(double(rb) / ra)) : 1.0;
        }
        f = std::min(std::max(1.0, std::round(f)), double(rb));
        return static_cast<std::size_t>(f);
    }

    std::vector<std::size_t> targetReplicas(std::size_t a, std::size_t arep, std::size_t b) const {
        const std::size_t rb = nrep(b);
        if (rb <= 1) return std::vector<std::size_t>(1, 0);
        const std::size_t f = fanOutOf(a, b, rb);
        std::vector<std::size_t> out;
        for (std::size_t k = 0; k < f; ++k) out.push_back((arep * f + k) % rb);
        return out;
    }

    Key hostKey(std::size_t tidx, std::size_t trep) const {
        const std::size_t h = parent(tidx);
        return Key(h, trep % nrep(h));
    }

    /** Replica 1 keeps the plain name; the suffix stays in [A-Za-z0-9_] (JSON object keys). */
    static std::string suffixed(const std::string& name, std::size_t rep) {
        return rep == 0 ? name : name + "_r" + std::to_string(rep + 1);
    }

    /** A subgraph copied per call site or replica repeats names; a repeat gets _d<n>. */
    std::string uniqueName(const std::string& base) {
        std::map<std::string, int>::iterator it = usedNames_.find(base);
        if (it == usedNames_.end()) {
            usedNames_[base] = 1;
            return base;
        }
        it->second += 1;
        return base + "_d" + std::to_string(it->second);
    }

    /** Local activity successors of aidx within an entry, in ascending index order. */
    std::vector<std::size_t> localSucc(std::size_t aidx, const std::vector<std::size_t>& local) const {
        std::vector<std::size_t> out;
        for (std::size_t s : l_.graph.succ(aidx))
            if (std::find(local.begin(), local.end(), s) != local.end()) out.push_back(s);
        return out;
    }

    /** Rows (target entry, probability) of the forwarding calls of an entry. */
    std::vector<std::pair<std::size_t, T> > forwardingOf(std::size_t eidx) const {
        std::vector<std::pair<std::size_t, T> > fwd;
        const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
        for (std::size_t c = 1; c <= l_.ncalls; ++c) {
            if (l_.calltype[c] != lang::CallType::FWD || l_.callpair_src[c] != eidx) continue;
            T p = l_.callproc_mean[c];
            if (p < zero) p = zero;
            if (p > one) p = one;
            if (dbl(p) > lang::GlobalConstants::FineTol) fwd.push_back(std::make_pair(l_.callpair_dst[c], p));
        }
        return fwd;
    }

    bool taskHasAndFork(std::size_t tidx) const {
        for (std::size_t a = l_.ashift + 1; a <= l_.ashift + l_.nacts; ++a)
            if (parent(a) == tidx && isPostAnd(a)) return true;
        return false;
    }

    // -------------------------------------------------------------- replication

    /** Copies of the task step graphs a materialised expansion would create. */
    double replInstantiations() const {
        std::map<std::size_t, double> memo;
        std::set<std::size_t> seeds(refTasks_.begin(), refTasks_.end());
        for (std::size_t e : openEntries_) seeds.insert(parent(e));
        double n = 0.0;
        for (std::size_t t : seeds) n += replRaw_[t] * taskCost(t, std::vector<std::size_t>(), memo);
        return n;
    }

    double taskCost(std::size_t tidx, const std::vector<std::size_t>& stack,
                    std::map<std::size_t, double>& memo) const {
        if (std::find(stack.begin(), stack.end(), tidx) != stack.end()) return 1.0;
        std::map<std::size_t, double>::const_iterator it = memo.find(tidx);
        if (it != memo.end()) return it->second;
        std::vector<std::size_t> deeper = stack;
        deeper.push_back(tidx);
        double c = 1.0;
        for (std::size_t eidx : l_.entriesof[tidx]) {
            std::vector<std::size_t> targets;
            for (std::size_t a : l_.actsof[eidx])
                for (std::size_t cidx : l_.callsof[a])
                    if (l_.calltype[cidx] == lang::CallType::SYNC ||
                        l_.calltype[cidx] == lang::CallType::ASYNC)
                        targets.push_back(l_.callpair_dst[cidx]);
            // A forwarded entry is expanded per replica just as a call is
            const std::vector<std::pair<std::size_t, T> > fwd = forwardingOf(eidx);
            for (std::size_t f = 0; f < fwd.size(); ++f) targets.push_back(fwd[f].first);
            for (std::size_t te : targets) {
                const std::size_t b = parent(te);
                c += double(fanOutOf(tidx, b, static_cast<std::size_t>(replRaw_[b]))) *
                     taskCost(b, deeper, memo);
            }
        }
        memo[tidx] = c;
        return c;
    }

    // -------------------------------------------------------------------- steps

    std::size_t addStep(std::size_t aidx, Key host, const Dist& svc, const std::string& name,
                        bool blocks, bool isthink, Key ref) {
        Step s;
        s.aidx = aidx;
        s.host = host;
        s.svc = svc;
        s.name = name;
        s.blocks = blocks;
        s.isthink = isthink;
        s.ref = ref;
        s.node = 0;
        s.owner = steps_.size();
        // Every thread-pool task replica on the stack holds a thread here; a sync caller releases it only on reply
        std::set<Key> held;
        for (const Key& k : threadStack_)
            if (fcrTask_[k.first]) held.insert(k);
        s.tasks.assign(held.begin(), held.end());
        steps_.push_back(s);
        return steps_.size() - 1;
    }

    /** A step on a Fork, Join or Router node: no station, no service, no class of its own. */
    std::size_t addAuxStep(std::size_t node, std::size_t ownerStep, const std::string& name, Key ref) {
        const std::size_t id = addStep(0, Key(0, 0), Dist::disabled_dist(), name, false, false, ref);
        steps_[id].node = node;
        steps_[id].owner = steps_[ownerStep].owner;
        return id;
    }

    void addRoute(const Port& from, std::size_t to, const T& p) {
        Flow f;
        f.from = from.step;
        f.to = to;
        f.p = p;
        f.fromSig = from.sig;
        f.inTarget = false;
        flow_.push_back(f);
    }

    /** Leaving a Cache node: the node switches into hit or miss, so the route is in the target class. */
    void addCacheRoute(std::size_t from, std::size_t to) {
        Flow f;
        f.from = from;
        f.to = to;
        f.p = num_traits<T>::from_int(1);
        f.fromSig = false;
        f.inTarget = true;
        flow_.push_back(f);
    }

    std::size_t cacheNodeOf(std::size_t tidx, std::size_t trep) {
        const Key k(tidx, trep);
        std::map<Key, std::size_t>::const_iterator it = cacheNodeOf_.find(k);
        if (it != cacheNodeOf_.end()) return it->second;
        qn::CacheParam<T> cp;
        cp.nitems = l_.nitems[tidx];
        cp.itemcap = l_.itemcap[tidx];
        cp.replacestrat = l_.replacestrat[tidx];
        const std::size_t nd = net_.add_cache(suffixed(nm(tidx) + "_Cache", trep), cp);
        cacheNodeOf_[k] = nd;
        return nd;
    }

    bool hasRetrieval(std::size_t tidx) const {
        return tidx < l_.hasretrieval.size() && l_.hasretrieval[tidx];
    }

    /** One PS fetch station per cache replica, as in the LN cache sublayer. */
    std::size_t fetchNodeOf(std::size_t tidx, std::size_t trep) {
        const Key k(tidx, trep);
        std::map<Key, std::size_t>::const_iterator it = fetchNodeOf_.find(k);
        if (it != fetchNodeOf_.end()) return it->second;
        const std::size_t nd =
            net_.add_queue(suffixed(nm(tidx) + "_Cache_Fetch", trep), lang::SchedStrategy::PS);
        fetchNodeOf_[k] = nd;
        return nd;
    }

    std::size_t actThinkStation() {
        if (actThinkNode_ == 0) actThinkNode_ = net_.add_delay("ActivityThink");
        return actThinkNode_;
    }

    std::vector<Stage> synchCallStages(std::size_t aidx) {
        std::vector<Stage> stages;
        const T one = num_traits<T>::from_int(1);
        for (std::size_t cidx : l_.callsof[aidx]) {
            const lang::CallType ct = l_.calltype[cidx];
            if (ct != lang::CallType::SYNC && ct != lang::CallType::ASYNC) continue;
            const bool isasync = ct == lang::CallType::ASYNC;
            const T m = l_.callproc_mean[cidx];
            std::size_t nfull =
                static_cast<std::size_t>(std::floor(dbl(m) + lang::GlobalConstants::FineTol));
            T frac = T(m - num_traits<T>::from_int(static_cast<long>(nfull)));
            if (nfull > MAXCALLSTAGES) {
                warn("Call multiplicity " + lqn::detail::fmt_num(dbl(m)) + " on " + l_.callnames[cidx] +
                     " truncated to " + std::to_string(MAXCALLSTAGES) + " stages.");
                nfull = MAXCALLSTAGES;
                frac = num_traits<T>::from_int(0);
            }
            for (std::size_t k = 0; k < nfull; ++k) stages.push_back(Stage{l_.callpair_dst[cidx], one, isasync});
            if (dbl(frac) > lang::GlobalConstants::FineTol)
                stages.push_back(Stage{l_.callpair_dst[cidx], frac, isasync});
        }
        return stages;
    }

    // ---------------------------------------------------------------- expansion

    Expansion expandEntry(std::size_t eidx, Key ref, std::size_t trep) {
        for (const Key& k : entryStack_)
            if (k.first == eidx && k.second == trep) {
                warn("Recursive call cycle at entry " + nm(eidx) + " truncated.");
                return Expansion();
            }
        entryStack_.push_back(Key(eidx, trep));
        threadStack_.push_back(Key(parent(eidx), trep));
        Expansion out = expandEntryBody(eidx, ref, trep);
        entryStack_.pop_back();
        threadStack_.pop_back();
        return out;
    }

    Expansion expandEntryBody(std::size_t eidx, Key ref, std::size_t trep) {
        Expansion out;
        if (eidx >= l_.actsof.size() || l_.actsof[eidx].empty()) return out;
        // The bound activity is the activity successor of the entry
        std::size_t bound = NONE;
        for (std::size_t s : l_.graph.succ(eidx)) {
            if (l_.type[s] != lang::LqnElement::ACTIVITY) continue;
            if (std::find(l_.actsof[eidx].begin(), l_.actsof[eidx].end(), s) == l_.actsof[eidx].end())
                continue;
            bound = s;
            break;
        }
        if (bound == NONE) return out;
        out = expandActivities(bound, eidx, ref, trep);

        // Forwarding: with prob p the entry hands off to a target that replies to the original caller
        const std::vector<std::pair<std::size_t, T> > fwd = forwardingOf(eidx);
        if (!fwd.empty() && !out.replies.empty()) {
            const std::vector<Exit> ownPorts = out.replies;
            std::vector<Exit> fwdExits;
            T pforw = num_traits<T>::from_int(0);
            // The forwarder's thread is released at the handoff
            const Key fwdThread = threadStack_.back();
            threadStack_.pop_back();
            const std::size_t fwdTidx = parent(eidx);
            for (std::size_t f = 0; f < fwd.size(); ++f) {
                const std::vector<std::size_t> reps = targetReplicas(fwdTidx, trep, parent(fwd[f].first));
                const T p = fwd[f].second;
                const T nreps = num_traits<T>::from_int(static_cast<long>(reps.size()));
                std::size_t reached = 0;
                for (std::size_t mrep : reps) {
                    Expansion fe = expandEntry(fwd[f].first, ref, mrep);
                    if (fe.first == NONE) continue;
                    ++reached;
                    for (const Exit& op : ownPorts) addRoute(Port(op.step, op.sig), fe.first, T(op.prob * p / nreps));
                    fwdExits.insert(fwdExits.end(), fe.replies.begin(), fe.replies.end());
                    fwdExits.insert(fwdExits.end(), fe.terms.begin(), fe.terms.end());
                }
                pforw = T(pforw + p * num_traits<T>::from_int(static_cast<long>(reached)) / nreps);
            }
            threadStack_.push_back(fwdThread);
            // What is left of each of this entry's own ports still replies
            T resid = T(num_traits<T>::from_int(1) - pforw);
            if (resid < num_traits<T>::from_int(0)) resid = num_traits<T>::from_int(0);
            for (Exit& ex : out.replies) ex.prob = T(ex.prob * resid);
            out.replies.insert(out.replies.end(), fwdExits.begin(), fwdExits.end());
        }
        return out;
    }

    Expansion expandActivities(std::size_t a0, std::size_t eidx, Key ref, std::size_t trep) {
        WalkCtx c;
        c.eidx = eidx;
        c.tidx = parent(a0);
        c.trep = trep;
        c.ref = ref;
        c.localActs = l_.actsof[eidx];
        Expansion out;
        out.first = walk(c, a0).first;
        // Every chain end returns the token to the caller: as a reply if the entry replies anywhere
        if (c.sawReply) {
            c.replyExits.insert(c.replyExits.end(), c.terminals.begin(), c.terminals.end());
            c.terminals.clear();
        }
        out.replies = c.replyExits;
        out.terms = c.terminals;
        return out;
    }

    std::string sname(const WalkCtx& c, const std::string& name) {
        return uniqueName(suffixed(name, c.trep));
    }

    /** Move the terminals appended since n0 into the phase-2 chain ends. */
    void takePh2Terminals(WalkCtx& c, std::size_t n0) {
        for (std::size_t r = n0; r < c.terminals.size(); ++r) {
            Ph2Exit pe;
            pe.step = c.terminals[r].step;
            pe.sig = c.terminals[r].sig;
            pe.p = c.terminals[r].prob;
            pe.ref = c.ref;
            ph2Exits_.push_back(pe);
        }
        c.terminals.resize(n0);
    }

    std::pair<std::size_t, Port> walk(WalkCtx& c, std::size_t aidx) {
        {
            typename std::map<std::size_t, std::pair<std::size_t, Port> >::const_iterator it =
                c.visited.find(aidx);
            if (it != c.visited.end()) return it->second;
        }
        const T one = num_traits<T>::from_int(1);
        const std::size_t tidx = c.tidx;
        // The activity bound to an ItemEntry of a CacheTask is the read step, on the Cache node
        std::size_t cacheNode = 0;
        if (tidx < l_.iscache.size() && l_.iscache[tidx] && dbl(l_.graph.get(c.eidx, aidx)) > 0.0)
            cacheNode = cacheNodeOf(tidx, c.trep);
        std::pair<std::size_t, Port> se = makeActivitySteps(c, aidx, cacheNode);
        const std::size_t entryStep = se.first;
        Port exitPort = se.second;
        c.visited[aidx] = se;

        // The reply is deferred to the chain ends, not emitted here
        const bool repliesHere = repliesTo(aidx, c.eidx);
        if (repliesHere) c.sawReply = true;

        const std::vector<std::size_t> succ = localSucc(aidx, c.localActs);
        if (succ.empty()) {
            c.terminals.push_back(Exit{exitPort.step, exitPort.sig, one});
            return se;
        }

        // Phase 2 runs after the reply via sn.classspawn and needs a station departure in the step's own class
        if (repliesHere) {
            const Key hk = hostKey(tidx, c.trep);
            if (cacheNode != 0 && succ.size() >= 2) {
                // Phase 2 at a cache read: each hit/miss outcome routes through its own immediate trigger step
                const std::size_t trigH =
                    addStep(aidx, hk, Dist::disabled_dist(), sname(c, nm(aidx) + "_ph2h"), false, false, c.ref);
                const std::size_t trigM =
                    addStep(aidx, hk, Dist::disabled_dist(), sname(c, nm(aidx) + "_ph2m"), false, false, c.ref);
                addCacheRoute(entryStep, trigH);
                addCacheRoute(entryStep, trigM);
                if (hasRetrieval(tidx))
                    warn("Delayed-hit retrieval of " + nm(tidx) +
                         " is not represented on a cache read with phase-2 successors.");
                cacheWiring_.push_back(CacheWire{cacheNode, entryStep, trigH, trigM, c.eidx, 0,
                                                 Dist::disabled_dist()});
                c.replyExits.push_back(Exit{trigH, false, one});
                c.replyExits.push_back(Exit{trigM, false, one});
                const std::vector<Key> savedStack = threadStack_;
                threadStack_.assign(1, Key(tidx, c.trep));
                const std::size_t nT0 = c.terminals.size();
                for (std::size_t hm = 0; hm < 2; ++hm) {
                    std::size_t sEntry = walk(c, succ[hm]).first;
                    if (steps_[sEntry].node != 0) {
                        const std::size_t head = addStep(aidx, hk, Dist::disabled_dist(),
                                                         sname(c, nm(aidx) + "_ph2b" + std::to_string(hm + 1)),
                                                         false, false, c.ref);
                        addRoute(Port(head, false), sEntry, one);
                        sEntry = head;
                    }
                    spawnPairs_.push_back(std::make_pair(hm == 0 ? trigH : trigM, sEntry));
                }
                takePh2Terminals(c, nT0);
                threadStack_ = savedStack;
                return se;
            }
            // Inside an AND-fork branch the lift applies only at the branch tail
            const bool okCtx = cacheNode == 0 && (c.forkOwnerStack.empty() || isAndJoinPre(aidx));
            // A merge step at the host station normalises a call-site exit to a station departure
            if (okCtx && (exitPort.sig || isRouter(steps_[exitPort.step].node))) {
                const std::size_t trig =
                    addStep(aidx, hk, Dist::disabled_dist(), sname(c, nm(aidx) + "_ph2t"), false, false, c.ref);
                addRoute(exitPort, trig, one);
                exitPort = Port(trig, false);
            }
            const bool ph2Spawn = okCtx && !exitPort.sig && steps_[exitPort.step].node == 0;
            if (!ph2Spawn) {
                warn("Phase-2 activities of " + nm(aidx) +
                     " run before the reply: the boundary is not a station departure, a degenerate "
                     "cache read, or mid-branch inside an AND-fork.");
            } else {
                c.replyExits.push_back(Exit{exitPort.step, exitPort.sig, one});
                const std::vector<Key> savedStack = threadStack_;
                threadStack_.assign(1, Key(tidx, c.trep));
                const std::size_t nT0 = c.terminals.size();
                std::vector<std::size_t> posSucc;
                for (std::size_t s : succ)
                    if (dbl(l_.graph.get(aidx, s)) > 0.0) posSucc.push_back(s);
                std::size_t target = NONE;
                if (isAndFork(succ)) {
                    // Phase 2 opens with an AND-fork: spawn into an immediate head and fork from there
                    const std::size_t head =
                        addStep(aidx, hk, Dist::disabled_dist(), sname(c, nm(aidx) + "_ph2"), false, false, c.ref);
                    wireAndFork(c, Port(head, false), succ, aidx);
                    target = head;
                } else if (isAndJoinPre(aidx)) {
                    // At an AND-join branch tail: the spawned head stands in for this branch at the Join
                    const std::size_t head =
                        addStep(aidx, hk, Dist::disabled_dist(), sname(c, nm(aidx) + "_ph2"), false, false, c.ref);
                    wireAndJoin(c, Port(head, false), succ[0]);
                    target = head;
                } else if (posSucc.size() == 1) {
                    const std::size_t sEntry = walk(c, posSucc[0]).first;
                    if (steps_[sEntry].node == 0) target = sEntry;
                }
                if (target == NONE) {
                    // Branching phase 2, or a head on a non-station node: an immediate head carries the branches
                    const std::size_t head =
                        addStep(aidx, hk, Dist::disabled_dist(), sname(c, nm(aidx) + "_ph2"), false, false, c.ref);
                    for (std::size_t s2 : posSucc) {
                        const std::size_t sEntry = walk(c, s2).first;
                        addRoute(Port(head, false), sEntry, l_.graph.get(aidx, s2));
                    }
                    target = head;
                }
                spawnPairs_.push_back(std::make_pair(exitPort.step, target));
                takePh2Terminals(c, nT0);
                threadStack_ = savedStack;
                return se;
            }
        }

        if (cacheNode != 0) {
            // CacheAccess: successors are the hit then the miss branch, and the Cache node decides
            if (succ.size() < 2) {
                warn("Cache read " + nm(aidx) + " has no hit/miss pair; treated as an ordinary activity.");
            } else {
                const std::size_t hEntry = walk(c, succ[0]).first;
                const std::size_t mEntry = walk(c, succ[1]).first;
                addCacheRoute(entryStep, hEntry);
                addCacheRoute(entryStep, mEntry);
                std::size_t fetch = 0;
                Dist fetchSvc = Dist::disabled_dist();
                if (hasRetrieval(tidx)) {
                    // The fetch is what the miss branch does, so its demand moves to the fetch station
                    fetch = fetchNodeOf(tidx, c.trep);
                    fetchSvc = steps_[mEntry].svc;
                    steps_[mEntry].svc = Dist::disabled_dist();
                    if (!l_.callsof[succ[1]].empty())
                        warn("Calls of miss activity " + nm(succ[1]) +
                             " stay outside the retrieval system, so they are not coalesced across "
                             "concurrent misses.");
                }
                cacheWiring_.push_back(CacheWire{cacheNode, entryStep, hEntry, mEntry, c.eidx, fetch, fetchSvc});
                return se;
            }
        }

        if (isAndFork(succ)) {
            wireAndFork(c, exitPort, succ, aidx);
            return se;
        }
        if (isAndJoinPre(aidx)) {
            wireAndJoin(c, exitPort, succ[0]);
            return se;
        }
        for (std::size_t s : succ) {
            const T p = l_.graph.get(aidx, s);
            if (!(dbl(p) > 0.0)) continue;
            const std::size_t sEntry = walk(c, s).first;
            addRoute(exitPort, sEntry, p);
        }
        return se;
    }

    /** AND-fork: a Fork replicates the job, one Router per branch since a Fork cannot switch class per link. */
    void wireAndFork(WalkCtx& c, const Port& from, std::vector<std::size_t> fsucc, std::size_t aidx) {
        const T one = num_traits<T>::from_int(1);
        const std::string forkName = sname(c, "Fork_" + nm(aidx));
        const std::size_t forkNode = net_.add_fork(forkName);
        const std::size_t forkStep = addAuxStep(forkNode, from.step, forkName, c.ref);
        addRoute(from, forkStep, one);
        c.forkOwnerStack.push_back(forkStep);
        // Walk a replying branch first so the Join and post-join subgraph are created in its phase-2 context
        std::vector<std::size_t> ordered;
        for (std::size_t s : fsucc)
            if (branchReplies(c, s)) ordered.push_back(s);
        for (std::size_t s : fsucc)
            if (!branchReplies(c, s)) ordered.push_back(s);
        for (std::size_t b = 0; b < ordered.size(); ++b) {
            const std::string routerName = sname(c, "Fork_" + nm(aidx) + "_" + std::to_string(b + 1));
            const std::size_t routerNode = net_.add_router(routerName);
            const std::size_t routerStep = addAuxStep(routerNode, forkStep, routerName, c.ref);
            addRoute(Port(forkStep, false), routerStep, one);
            const std::size_t sEntry = walk(c, ordered[b]).first;
            addRoute(Port(routerStep, false), sEntry, one);
        }
        c.forkOwnerStack.pop_back();
    }

    /** Route a branch tail into the AND-join, created on first arrival, in the class that entered the fork. */
    void wireAndJoin(WalkCtx& c, const Port& from, std::size_t joinAidx) {
        const T one = num_traits<T>::from_int(1);
        std::map<std::size_t, std::size_t>::const_iterator it = c.joinOf.find(joinAidx);
        if (it != c.joinOf.end()) {
            // A Join closes one fork, so every branch tail must arrive with that fork innermost
            if (c.forkOwnerStack.empty() || c.forkOwnerStack.back() != c.joinFork.at(joinAidx))
                throw InputError("LQN2QN: AND-join at " + nm(joinAidx) +
                                 " merges branches of different AND-forks: the fork-join structure is not"
                                 " nested, so it has no Fork/Join representation.");
            addRoute(from, it->second, one);
            return;
        }
        if (c.forkOwnerStack.empty()) {
            warn("AND-join at " + nm(joinAidx) + " has no enclosing AND-fork; branches are serialised.");
            const std::size_t sEntry = walk(c, joinAidx).first;
            addRoute(from, sEntry, one);
            return;
        }
        const std::string joinName = sname(c, "Join_" + nm(joinAidx));
        const std::size_t forkOwner = c.forkOwnerStack.back();
        const std::size_t joinNode = net_.add_join(joinName, steps_[forkOwner].node);
        const std::size_t joinStep = addAuxStep(joinNode, forkOwner, joinName, c.ref);
        addRoute(from, joinStep, one);
        // Registered before the walk, which runs outside the fork the Join closes
        c.joinOf[joinAidx] = joinStep;
        c.joinFork[joinAidx] = forkOwner;
        const std::vector<std::size_t> savedForks = c.forkOwnerStack;
        c.forkOwnerStack.pop_back();
        const std::size_t sEntry = walk(c, joinAidx).first;
        c.forkOwnerStack = savedForks;
        addRoute(Port(joinStep, false), sEntry, one);
        joinQuorum_.push_back(std::make_pair(joinStep, joinAidx));
    }

    /** True if the branch rooted at a0 replies to the current entry, searching up to the closing AND-join. */
    bool branchReplies(const WalkCtx& c, std::size_t a0) const {
        std::vector<std::size_t> stack(1, a0), seen;
        while (!stack.empty()) {
            const std::size_t a = stack.back();
            stack.pop_back();
            if (std::find(seen.begin(), seen.end(), a) != seen.end()) continue;
            seen.push_back(a);
            if (repliesTo(a, c.eidx)) return true;
            if (isAndJoinPre(a)) continue;
            const std::vector<std::size_t> nxt = localSucc(a, c.localActs);
            stack.insert(stack.end(), nxt.begin(), nxt.end());
        }
        return false;
    }

    /** One step for the host demand, plus one per unrolled call stage; returns (entry step, exit port). */
    std::pair<std::size_t, Port> makeActivitySteps(WalkCtx& c, std::size_t aidx, std::size_t cacheNode) {
        const T one = num_traits<T>::from_int(1);
        const std::size_t tidx = c.tidx, trep = c.trep;
        const Key hidx = hostKey(tidx, trep);
        if (cacheNode != 0) {
            // A read step holds no demand and issues no call: the work is on the hit or miss branch
            const std::size_t entryStep = addStep(aidx, hidx, Dist::disabled_dist(),
                                                  uniqueName(suffixed(nm(aidx), trep)), false, false, c.ref);
            steps_[entryStep].node = cacheNode;
            if (!l_.callsof[aidx].empty())
                warn("Calls issued by cache read activity " + nm(aidx) + " are ignored.");
            if (timed(l_.hostdem[aidx]))
                warn("Host demand of cache read activity " + nm(aidx) + " is ignored.");
            return std::make_pair(entryStep, Port(entryStep, false));
        }
        const Dist svc = timed(l_.hostdem[aidx]) ? l_.hostdem[aidx] : Dist::disabled_dist();
        const std::vector<Stage> callStages = synchCallStages(aidx);
        const std::size_t entryStep =
            addStep(aidx, hidx, svc, uniqueName(suffixed(nm(aidx), trep)), false, false, c.ref);
        Port cur(entryStep, false);

        // Activity think time: a delay in series with the host demand on a shared INF station
        if (aidx < l_.actthink.size() && timed(l_.actthink[aidx])) {
            const std::size_t thinkStep = addStep(aidx, hidx, l_.actthink[aidx],
                                                  uniqueName(suffixed(nm(aidx) + "_think", trep)), false,
                                                  false, c.ref);
            steps_[thinkStep].node = actThinkStation();
            addRoute(cur, thinkStep, one);
            cur = Port(thinkStep, false);
        }

        // A call holds the caller's server only at a finite, non-pool, closed-chain task alone on its host
        const bool hostBlocks = !hostIsDelay_[hidx] && !fcrTask_[tidx] && c.ref.first != 0 &&
                                std::isfinite(l_.mult[tidx]) && l_.sched[tidx] != lang::SchedStrategy::INF &&
                                hostNTasks_.at(l_.parent[tidx]) == 1;

        for (std::size_t k = 0; k < callStages.size(); ++k) {
            const Stage& stage = callStages[k];
            // An async call does not hold the caller's server; it stays serialised, the approximation
            const bool blocks = hostBlocks && !stage.isasync;
            // An async send also releases the caller's thread for the callee expansion
            Key asyncThread(0, 0);
            if (stage.isasync) {
                asyncThread = threadStack_.back();
                threadStack_.pop_back();
            }
            // The stage reaches the callee replicas this caller replica addresses, sharing the mean uniformly
            const std::vector<std::size_t> reps = targetReplicas(tidx, trep, parent(stage.target));
            std::vector<std::size_t> calleeFirsts;
            std::vector<Exit> calleeReplies;
            for (std::size_t mrep : reps) {
                Expansion ce = expandEntry(stage.target, c.ref, mrep);
                if (ce.first == NONE) continue;
                calleeFirsts.push_back(ce.first);
                // A callee path that neither replies nor continues still returns the token
                calleeReplies.insert(calleeReplies.end(), ce.replies.begin(), ce.replies.end());
                calleeReplies.insert(calleeReplies.end(), ce.terms.begin(), ce.terms.end());
            }
            if (stage.isasync) threadStack_.push_back(asyncThread);
            if (calleeFirsts.empty()) continue;  // callee not expandable: drop the call, never block
            const T share = T(one / num_traits<T>::from_int(static_cast<long>(calleeFirsts.size())));

            // A merge step is needed when the call may be skipped, another stage follows, the call does
            // not block, or an AND-join branch tail must reach the Join in an ordinary class
            const bool needsMerge =
                stage.prob < one || k + 1 < callStages.size() || !blocks || isAndJoinPre(aidx);
            std::size_t nxt = NONE;
            if (needsMerge) {
                const std::string retName =
                    uniqueName(suffixed(nm(aidx) + "_c" + std::to_string(k + 1) + "_ret", trep));
                nxt = addStep(aidx, hidx, Dist::disabled_dist(), retName, false, false, c.ref);
                // A non-blocking return carries no signal, so its merge point can sit on a Router
                if (!blocks) steps_[nxt].node = net_.add_router(retName);
            }

            if (blocks) {
                // The first mandatory call binds to the service class itself: a class switch would release the server
                std::size_t blk;
                if (k == 0 && !(stage.prob < one) && cur == Port(entryStep, false)) {
                    blk = entryStep;
                    steps_[blk].blocks = true;
                } else {
                    blk = addStep(aidx, hidx, Dist::disabled_dist(),
                                  uniqueName(suffixed(nm(aidx) + "_c" + std::to_string(k + 1), trep)), true,
                                  false, c.ref);
                    addRoute(cur, blk, stage.prob);
                    if (stage.prob < one) addRoute(cur, nxt, T(one - stage.prob));
                }
                // Every reached replica replies into the same signal, so the call site blocks once
                for (std::size_t cf : calleeFirsts) addRoute(Port(blk, false), cf, share);
                for (const Exit& r : calleeReplies) reply_.push_back(Reply{r.step, blk, r.sig, r.prob});
                if (needsMerge) {
                    addRoute(Port(blk, true), nxt, one);
                    cur = Port(nxt, false);
                } else {
                    cur = Port(blk, true);
                }
            } else {
                for (std::size_t cf : calleeFirsts) addRoute(cur, cf, T(stage.prob * share));
                if (stage.prob < one) addRoute(cur, nxt, T(one - stage.prob));
                for (const Exit& r : calleeReplies) addRoute(Port(r.step, r.sig), nxt, r.prob);
                cur = Port(nxt, false);
            }
        }
        return std::make_pair(entryStep, cur);
    }

    std::size_t stationOf(std::size_t i) {
        const Step& s = steps_[i];
        if (s.node != 0) return s.node;
        if (s.isthink) return thinkNode_.at(s.ref);
        return hostStation_.at(s.host);
    }

    // -------------------------------------------------------------------- build

    void build() {
        const lqn::LqnStruct<T>& l = l_;
        const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
        const double FineTol = lang::GlobalConstants::FineTol;
        if (replication_ != "auto" && replication_ != "materialize" && replication_ != "pool")
            throw InputError("LQN2QN: replication must be 'auto', 'materialize' or 'pool'.");
        const std::size_t NT = l.nhosts + l.ntasks;

        for (std::size_t t = l.tshift + 1; t <= l.tshift + l.ntasks; ++t)
            if (l.isref[t]) refTasks_.push_back(t);

        // Entries with an open arrival process: Source -> open classes -> Sink
        for (std::size_t e = l.eshift + 1; e <= l.eshift + l.nentries; ++e) {
            if (e >= l.arrival.size() || !l.has_arrival[e] || l.arrival[e].disabled) continue;
            const double m = dbl(l.arrival[e].mean);
            if (std::isfinite(m) && m > FineTol) openEntries_.push_back(e);
        }
        if (refTasks_.empty() && openEntries_.empty())
            throw InputError("LQN2QN: LQN must have at least one reference task or open arrival.");

        for (std::size_t c = 1; c <= l.ncalls; ++c)
            if (l.calltype[c] == lang::CallType::ASYNC) {
                warn("Asynchronous calls are represented by LQN2QN as non-blocking visits: the caller "
                     "releases its server but remains serialised behind the callee.");
                break;
            }

        // ---- replication
        replRaw_.assign(NT + 1, 1.0);
        std::vector<std::size_t> replicated;
        for (std::size_t i = 1; i <= NT && i < l.repl.size(); ++i) {
            replRaw_[i] = std::max(1.0, std::round(l.repl[i]));
            if (replRaw_[i] > 1.0) replicated.push_back(i);
        }
        materialize_ = !replicated.empty() &&
                       (replication_ == "materialize" ||
                        (replication_ == "auto" && replInstantiations() <= double(MAXREPLINSTANCES)));
        if (!replicated.empty() && !materialize_)
            warn("Replication of " + nm(replicated[0]) +
                 " is pooled: its replicas become one station of r times the servers, one admission row "
                 "of r times the bound and one reference class of r times the population, which is exact "
                 "at an infinite-server host and optimistic elsewhere. Pass 'materialize' for one station "
                 "and one step-graph copy per replica.");

        // A processor hosting any task besides the caller is not held across a synchronous call: LQN releases the
        // processor while a thread waits, and holding it blocks the other tasks and deadlocks a call back onto it
        hostNTasks_.clear();
        for (std::size_t t = l.tshift + 1; t <= l.tshift + l.ntasks; ++t) ++hostNTasks_[l.parent[t]];

        // ---- tasks whose multiplicity is a thread pool
        fcrTask_.assign(NT + 1, false);
        for (std::size_t t = l.tshift + 1; t <= l.tshift + l.ntasks; ++t) {
            if (l.isref[t] || (t < l.iscache.size() && l.iscache[t])) continue;
            if (!std::isfinite(l.mult[t]) || l.sched[t] == lang::SchedStrategy::INF) continue;
            if (taskHasAndFork(t)) {
                warn("Multiplicity of task " + nm(t) +
                     " is not enforced: an AND-fork inside a task cannot be capped by a finite capacity "
                     "region, whose job count would double-count the forked siblings.");
                continue;
            }
            fcrTask_[t] = true;
        }
        // A task called transitively from inside an AND-fork branch is excluded too
        std::vector<std::size_t> branchActs;
        for (std::size_t a = l.ashift + 1; a <= l.ashift + l.nacts; ++a) {
            if (!isPostAnd(a)) continue;
            std::vector<std::size_t> frontier(1, a);
            while (!frontier.empty()) {
                const std::size_t cur = frontier.front();
                frontier.erase(frontier.begin());
                if (std::find(branchActs.begin(), branchActs.end(), cur) != branchActs.end()) continue;
                branchActs.push_back(cur);
                if (isAndJoinPre(cur)) continue;  // branch tail: do not traverse past the join
                for (std::size_t s : l.graph.succ(cur))
                    if (s > l.ashift && parent(s) == parent(cur)) frontier.push_back(s);
            }
        }
        if (!branchActs.empty() && std::find(fcrTask_.begin(), fcrTask_.end(), true) != fcrTask_.end()) {
            std::vector<std::size_t> front;
            for (std::size_t a : branchActs)
                for (std::size_t cidx : l.callsof[a]) front.push_back(parent(l.callpair_dst[cidx]));
            std::vector<bool> shadow(NT + 1, false);
            while (!front.empty()) {
                const std::size_t t = front.front();
                front.erase(front.begin());
                if (shadow[t]) continue;
                shadow[t] = true;
                for (std::size_t a = l.ashift + 1; a <= l.ashift + l.nacts; ++a)
                    if (parent(a) == t)
                        for (std::size_t cidx : l.callsof[a]) front.push_back(parent(l.callpair_dst[cidx]));
            }
            for (std::size_t t = 1; t <= NT; ++t)
                if (shadow[t] && fcrTask_[t]) {
                    fcrTask_[t] = false;
                    warn("Multiplicity of task " + nm(t) +
                         " is not enforced: it is called from inside an AND-fork branch, whose flows the "
                         "fork-join transformation retags outside the admission constraint.");
                }
        }

        // ---- stations: one per host processor replica
        for (std::size_t h = 1; h <= l.nhosts; ++h) {
            const double nservers = l.mult[h];
            for (std::size_t m = 0; m < nrep(h); ++m) {
                const Key k(h, m);
                if (std::isinf(nservers) || l.sched[h] == lang::SchedStrategy::INF) {
                    hostStation_[k] = net_.add_delay(suffixed(nm(h), m));
                    hostIsDelay_[k] = true;
                } else {
                    const std::size_t q = net_.add_queue(suffixed(nm(h), m), l.sched[h]);
                    // A pooled processor carries the servers of all its replicas
                    net_.set_number_of_servers(q, nservers * poolFactor(h));
                    hostStation_[k] = q;
                    hostIsDelay_[k] = false;
                }
            }
        }
        // ---- think delays, one per reference task replica
        for (std::size_t rt : refTasks_)
            for (std::size_t m = 0; m < nrep(rt); ++m)
                thinkNode_[Key(rt, m)] = net_.add_delay(suffixed(nm(rt) + "_Think", m));

        // ---- pass 1: one closed chain per reference task replica
        for (std::size_t rt : refTasks_)
            for (std::size_t rep = 0; rep < nrep(rt); ++rep) {
                const Key refKey(rt, rep);
                const std::size_t thinkStep =
                    addStep(0, Key(0, 0), Dist::disabled_dist(),
                            uniqueName(suffixed(nm(rt) + "_Think", rep)), false, true, refKey);
                for (std::size_t eidx : l.entriesof[rt]) {
                    const Expansion ex = expandEntry(eidx, refKey, rep);
                    if (ex.first == NONE) continue;
                    addRoute(Port(thinkStep, false), ex.first, one);
                    // A reference task has no caller: replies and dead ends both close at the think delay
                    for (const Exit& s : ex.replies) addRoute(Port(s.step, s.sig), thinkStep, s.prob);
                    for (const Exit& s : ex.terms) addRoute(Port(s.step, s.sig), thinkStep, s.prob);
                }
            }

        // ---- pass 1b: open arrival chains
        struct OpenWire {
            std::size_t eidx, first;
            std::vector<Exit> exits;
        };
        std::vector<OpenWire> openWiring;
        if (!openEntries_.empty()) {
            srcNode_ = net_.add_source("Source");
            snkNode_ = net_.add_sink("Sink");
        }
        for (std::size_t eidx : openEntries_) {
            // Each replica of the entry's task receives its own arrival stream
            for (std::size_t rep = 0; rep < nrep(parent(eidx)); ++rep) {
                Expansion ex = expandEntry(eidx, Key(0, 0), rep);
                if (ex.first == NONE) {
                    warn("Open arrival entry " + nm(eidx) + " has no bound activity; ignored.");
                    continue;
                }
                OpenWire ow;
                ow.eidx = eidx;
                ow.first = ex.first;
                ow.exits = ex.replies;
                ow.exits.insert(ow.exits.end(), ex.terms.begin(), ex.terms.end());
                openWiring.push_back(ow);
            }
        }

        // ---- pass 2: classes and reply signals
        const std::size_t nsteps = steps_.size();
        std::vector<std::size_t> stepClass(nsteps, 0), stepSignal(nsteps, 0);
        for (std::size_t i = 0; i < nsteps; ++i) {
            const Step& s = steps_[i];
            if (s.owner != i) continue;  // Fork/Join/Router steps travel in the class that entered the fork
            if (s.ref.first == 0) {
                stepClass[i] = net_.add_open_class(s.name);
            } else if (s.isthink) {
                // A pooled reference task holds the population of all its replicas
                stepClass[i] = net_.add_closed_class(s.name, l.mult[s.ref.first] * poolFactor(s.ref.first),
                                                     thinkNode_.at(s.ref));
            } else {
                stepClass[i] = net_.add_closed_class(s.name, 0.0, thinkNode_.at(s.ref));
            }
        }
        for (std::size_t i = 0; i < nsteps; ++i) stepClass[i] = stepClass[steps_[i].owner];

        // AND-join quorum, in the class the siblings are matched in
        for (const std::pair<std::size_t, std::size_t>& jq : joinQuorum_) {
            const std::size_t joinAidx = jq.second;
            if (joinAidx >= l.actquorum.size()) continue;
            const std::size_t quorum = l.actquorum[joinAidx];
            std::size_t nb = 0;
            for (std::size_t p : l.graph.pred(joinAidx))
                if (p != joinAidx && isAndJoinPre(p)) ++nb;
            // A quorum equal to the branch count is the default wait-for-all
            if (quorum < 1 || nb < 1 || quorum >= nb) continue;
            net_.set_join_strategy(steps_[jq.first].node, lang::JoinStrategy::PARTIAL, double(quorum));
        }

        for (std::size_t i = 0; i < nsteps; ++i) {
            if (!steps_[i].blocks) continue;
            const std::size_t sig = net_.add_closed_class(steps_[i].name + "_Reply", 0.0, thinkNode_.at(steps_[i].ref));
            net_.set_reply_signal_class(stepClass[i], sig);
            stepSignal[i] = sig;
        }

        // Spawn bindings for phase-2 continuations
        for (const std::pair<std::size_t, std::size_t>& sp : spawnPairs_)
            net_.set_class_spawn(stepClass[sp.first], stepClass[sp.second]);

        // Phase-2 token destructor on closed chains: a NEGATIVE signal at a station nothing visits
        std::size_t ph2Dump = 0;
        std::map<Key, std::size_t> ph2DestructorOf;
        {
            std::set<Key> ph2Ref;
            for (const Ph2Exit& pe : ph2Exits_)
                if (pe.ref.first > 0) ph2Ref.insert(pe.ref);
            if (!ph2Ref.empty()) {
                ph2Dump = net_.add_queue("Ph2Sink", lang::SchedStrategy::FCFS);
                for (const Key& k : ph2Ref) {
                    const std::size_t sig =
                        net_.add_closed_class(suffixed("Ph2End_" + nm(k.first), k.second), 0.0, thinkNode_.at(k));
                    net_.set_signal(sig, lang::SignalType::NEGATIVE);
                    net_.set_service(ph2Dump, sig, Dist::immediate());
                    ph2DestructorOf[k] = sig;
                }
            }
        }

        // ---- pass 3: service times
        for (std::size_t i = 0; i < nsteps; ++i) {
            const Step& s = steps_[i];
            if (s.node != 0) {
                if (s.owner == i && actThinkNode_ != 0 && s.node == actThinkNode_) {
                    net_.set_service(actThinkNode_, stepClass[i], s.svc);
                    continue;
                }
                // A Router-hosted merge step owns a class, declared Immediate where the class is referenced
                if (s.owner == i && isRouter(s.node)) {
                    if (s.ref.first == 0)
                        net_.set_service(hostStation_.at(s.host), stepClass[i], Dist::immediate());
                    else
                        net_.set_service(thinkNode_.at(s.ref), stepClass[i], Dist::immediate());
                }
                continue;
            }
            if (s.isthink) {
                const Dist& th = l.think[s.ref.first];
                net_.set_service(thinkNode_.at(s.ref), stepClass[i], timed(th) ? th : Dist::immediate());
            } else {
                net_.set_service(hostStation_.at(s.host), stepClass[i], s.svc.disabled ? Dist::immediate() : s.svc);
            }
        }

        // A SetupTask is the Queue setup/delay-off pair at its host station
        std::set<std::size_t> warnedSetup;
        for (std::size_t i = 0; i < nsteps; ++i) {
            const Step& s = steps_[i];
            if (s.node != 0 || s.isthink || s.aidx == 0) continue;
            const std::size_t t = parent(s.aidx);
            if (t >= l.hassetup.size() || !l.hassetup[t] || !timed(l.setuptime[t])) continue;
            if (l.delayofftime[t].disabled) {
                if (warnedSetup.insert(t).second)
                    warn("Setup of setup task " + nm(t) +
                         " is not represented: it has no delay-off time, so its server never shuts down and "
                         "never sets up again.");
                continue;
            }
            if (hostIsDelay_[s.host]) {
                if (warnedSetup.insert(t).second)
                    warn("Setup of setup task " + nm(t) +
                         " is not represented: its processor is an infinite server, which never shuts down.");
                continue;
            }
            net_.set_setup_delayoff(hostStation_.at(s.host), stepClass[i], l.setuptime[t], l.delayofftime[t]);
        }

        // A reply signal is consumed at the caller's station, declared everywhere it may fall back
        for (std::size_t i = 0; i < nsteps; ++i) {
            if (stepSignal[i] == 0) continue;
            for (const std::pair<const Key, std::size_t>& hs : hostStation_)
                net_.set_service(hs.second, stepSignal[i], Dist::immediate());
            for (const std::pair<const Key, std::size_t>& tn : thinkNode_)
                net_.set_service(tn.second, stepSignal[i], Dist::immediate());
        }

        // ---- cache read/hit/miss wiring, now that the classes exist
        for (const CacheWire& cw : cacheWiring_) {
            const std::size_t rc = stepClass[cw.readStep];
            {
                qn::CacheParam<T>& cp = net_.raw_struct().nodeparam.at(cw.node);
                resizeCacheParam(cp);
                // The item pmf of the entry, over the cache's items: an item beyond the entry's cardinality is never read
                std::vector<T> pmf(cp.nitems, zero);
                const std::vector<T>& ip = l.itemproc[cw.eidx];
                for (std::size_t k = 0; k < cp.nitems && k < ip.size(); ++k) pmf[k] = ip[k];
                cp.pread[rc - 1] = pmf;
                // setReadItemEntry copies the entry's DiscreteSampler, the law the JMT export needs beside the pmf
                cp.preadkind[rc - 1].type = lang::ProcessType::DISCRETESAMPLER;
                cp.preadkind[rc - 1].n = cp.nitems;
                cp.hitclass[rc - 1] = stepClass[cw.hitStep];
                cp.missclass[rc - 1] = stepClass[cw.missStep];
            }
            if (cw.fetch != 0) {
                // Service and routing of the retrieval system are read off the read class
                net_.set_service(cw.fetch, rc, cw.fetchSvc.disabled ? Dist::immediate() : cw.fetchSvc);
                net_.set_retrieval_system(cw.node, rc, stepClass[cw.missStep], std::vector<std::size_t>(1, cw.fetch));
            }
        }
        for (std::pair<const std::size_t, qn::CacheParam<T> >& np : net_.raw_struct().nodeparam)
            resizeCacheParam(np.second);

        // ---- pass 4: routing
        qn::RoutingMatrix<T> P = net_.init_routing_matrix();
        for (const Flow& f : flow_) {
            const std::size_t si = stationOf(f.from), sj = stationOf(f.to);
            if (f.inTarget)
                P.set(stepClass[f.to], stepClass[f.to], si, sj, f.p);
            else if (f.fromSig)
                P.set(stepSignal[f.from], stepClass[f.to], si, sj, f.p);
            else
                P.set(stepClass[f.from], stepClass[f.to], si, sj, f.p);
        }
        for (const Reply& r : reply_) {
            // A nested call returns through its own reply signal, switching into the outer call site's signal
            const std::size_t src = r.viaSig ? stepSignal[r.exit] : stepClass[r.exit];
            P.set(src, stepSignal[r.owner], stationOf(r.exit), stationOf(r.owner), r.p);
        }
        // Retrieval systems: the read class circulates cache -> fetch -> cache
        for (const CacheWire& cw : cacheWiring_)
            if (cw.fetch != 0) {
                const std::size_t rc = stepClass[cw.readStep];
                P.set(rc, rc, cw.node, cw.fetch, one);
                P.set(rc, rc, cw.fetch, cw.node, one);
            }
        // Phase-2 chain ends: destroy the spawned token
        for (const Ph2Exit& pe : ph2Exits_) {
            if (pe.ref.first == 0) {
                const std::size_t ec = stepClass[pe.step];
                P.set(ec, ec, stationOf(pe.step), snkNode_, pe.p);
            } else {
                const std::size_t src = pe.sig ? stepSignal[pe.step] : stepClass[pe.step];
                P.set(src, ph2DestructorOf.at(pe.ref), stationOf(pe.step), ph2Dump, pe.p);
            }
        }
        // Open arrival wiring: Source into the first step, exits into the Sink
        for (const OpenWire& ow : openWiring) {
            const std::size_t fc = stepClass[ow.first];
            net_.set_arrival(srcNode_, fc, l.arrival[ow.eidx]);
            P.set(fc, fc, srcNode_, stationOf(ow.first), one);
            // Open chains carry no signals, so every exit is an ordinary class
            for (const Exit& ex : ow.exits) {
                const std::size_t ec = stepClass[ex.step];
                P.set(ec, ec, stationOf(ex.step), snkNode_, ex.prob);
            }
        }
        net_.link(P);

        // ---- thread pools: one finite capacity region, one admission row per task replica
        std::set<Key> fcrSet;
        for (const Step& s : steps_) fcrSet.insert(s.tasks.begin(), s.tasks.end());
        if (!fcrSet.empty()) {
            const std::vector<Key> fcrList(fcrSet.begin(), fcrSet.end());
            const std::size_t K = net_.raw_struct().classes.size();
            Matrix<T> A(fcrList.size(), K, zero);
            std::vector<T> b(fcrList.size(), zero);
            std::vector<std::size_t> regionNodes;
            bool anyA = false;
            for (std::size_t ts = 0; ts < fcrList.size(); ++ts) {
                for (std::size_t i = 0; i < nsteps; ++i) {
                    const std::vector<Key>& tk = steps_[i].tasks;
                    if (std::find(tk.begin(), tk.end(), fcrList[ts]) == tk.end()) continue;
                    A(ts, stepClass[i] - 1) = one;
                    anyA = true;
                    const std::size_t nd = stationOf(i);
                    if (isStation(nd) && std::find(regionNodes.begin(), regionNodes.end(), nd) == regionNodes.end())
                        regionNodes.push_back(nd);
                    // The caller holds its thread for the whole fetch
                    for (const CacheWire& cw : cacheWiring_)
                        if (cw.readStep == i && cw.fetch != 0 &&
                            std::find(regionNodes.begin(), regionNodes.end(), cw.fetch) == regionNodes.end())
                            regionNodes.push_back(cw.fetch);
                }
                // A pooled task keeps one row whose bound covers all its replicas
                const std::size_t t = fcrList[ts].first;
                b[ts] = tnum(l.mult[t] * poolFactor(t));
            }
            if (anyA && !regionNodes.empty()) {
                const std::size_t rg = net_.add_region(regionNodes, std::vector<double>());
                net_.set_region_constraint(rg, A, b);
            }
        }
    }

    bool isStation(std::size_t node) const {
        const qn::NetworkStruct<T>& sn = net_raw();
        return sn.nodes[node - 1].station != 0 && sn.nodes[node - 1].nodetype != lang::NodeType::Source;
    }

    void resizeCacheParam(qn::CacheParam<T>& cp) {
        const std::size_t K = net_.raw_struct().classes.size();
        if (cp.pread.size() < K) cp.pread.resize(K);
        if (cp.preadkind.size() < K) cp.preadkind.resize(K);
        if (cp.hitclass.size() < K) cp.hitclass.resize(K, 0);
        if (cp.missclass.size() < K) cp.missclass.resize(K, 0);
        if (cp.classitem.size() < K) cp.classitem.resize(K, 0);
    }
};

/** Entry/activity index pairs that reply, from the explicit and the inferred replies of each entry. */
template <class T>
std::set<std::pair<std::size_t, std::size_t> > lqn2qn_replies(const lqn::LqnModel<T>& m,
                                                              const lqn::LqnStruct<T>& l) {
    std::map<std::string, std::size_t> ent, act;
    for (std::size_t i = 1; i <= l.nidx; ++i) {
        if (l.type[i] == lang::LqnElement::ENTRY) ent[l.names[i]] = i;
        if (l.type[i] == lang::LqnElement::ACTIVITY) act[l.names[i]] = i;
    }
    std::set<std::pair<std::size_t, std::size_t> > out;
    const std::map<std::string, std::vector<std::string> > rep = lqn::detail::reply_activities(m, l);
    for (const std::pair<const std::string, std::vector<std::string> >& kv : rep) {
        std::map<std::string, std::size_t>::const_iterator ei = ent.find(kv.first);
        if (ei == ent.end()) continue;
        for (const std::string& a : kv.second) {
            std::map<std::string, std::size_t>::const_iterator ai = act.find(a);
            if (ai != act.end()) out.insert(std::make_pair(ai->second, ei->second));
        }
    }
    return out;
}

}  // namespace detail

/**
 * Port of MATLAB `LQN2QN(lqn, replication)`: flatten a layered network into a
 * queueing network whose synchronous calls block through REPLY signals.
 *
 * @param lqn          the layered model, from lqn::LqnBuilder::model() or lqn::read_lqnx_model
 * @param replication  'auto' (default), 'materialize' or 'pool', as in the reference
 * @param name         the layered model's name; the network is called `<name>-QN`
 * @param warnings     when non-null, the reference's line_warning messages are appended here
 *                     instead of being printed to stderr
 * @return             the queueing network, routed and ready for get_struct()
 */
template <class T>
qn::Network<T> lqn2qn(const lqn::LqnModel<T>& lqn, const std::string& replication = "auto",
                      const std::string& name = "model", std::vector<std::string>* warnings = nullptr) {
    const lqn::LqnStruct<T> l = lqn::lqn_finalize(lqn);
    const std::set<std::pair<std::size_t, std::size_t> > replies = detail::lqn2qn_replies(lqn, l);
    detail::Lqn2Qn<T> conv(l, replies, replication, name, warnings);
    return conv.run();
}

/**
 * LQN2QN over an already flattened LqnStruct.
 *
 * The struct carries no replygraph, so only the IMPLICIT replies (leaf activities) are known here;
 * a model with an explicit reply on a non-leaf activity (phase 2) must use the LqnModel overload.
 */
template <class T>
qn::Network<T> lqn2qn(const lqn::LqnStruct<T>& lsn, const std::string& replication = "auto",
                      const std::string& name = "model", std::vector<std::string>* warnings = nullptr) {
    const std::set<std::pair<std::size_t, std::size_t> > replies =
        detail::lqn2qn_replies(lqn::LqnModel<T>(), lsn);
    detail::Lqn2Qn<T> conv(lsn, replies, replication, name, warnings);
    return conv.run();
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_LQN2QN_H
