/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_LQN_LQN_BALANCE_EQUATIONS_H
#define LINE_API_LQN_LQN_BALANCE_EQUATIONS_H

/**
 * Conservation laws of a layered queueing network, enumerated from its structure.
 *
 * Port of `matlab/src/api/lqn/lqn_balance_equations.m`, twin of the JAR
 * `jline.api.lqn.LqnBalanceEquations` and of the Python
 * `line_solver.api.lqn.balance_equations`.
 *
 * A layered model is not free to report any tuple of throughputs, think times and
 * utilizations: five families of relations tie them together, and every one of them
 * is fixed by the STRUCTURE of the model alone. This walks an `LqnStruct` and emits
 * them, one record per relation, with the index sets to aggregate over, the constant
 * coefficients and a printable form. NOTHING IS SOLVED HERE.
 *
 * The families, with `kind` as emitted:
 *
 * `little`    Little's law on a task's THREAD POOL. The threads of task t form a
 *             closed cycle of one delay stage (the surrogate think time SolverLN
 *             imputes to the task, plus the declared think time of a reference task)
 *             and one service stage (holding a request from above). With B(t,k) the
 *             mean number of threads of t busy serving caller class k -- the
 *             per-class utilization in JOB units --
 *
 *                 X(t)*(Z(t) + z(t)) + sum_k B(t,k) = N(t)
 *
 *             which is the update `updateThinkTimes` iterates on. In the
 *             [0,1]-normalized utilization LINE reports for a queueing station,
 *             B(t,k) = N(t)*U(t,k), giving X(t)*(Z(t)+z(t)) = N(t)*(1 - sum_k
 *             U(t,k)); at an infinite server the utilization is already a job count,
 *             so B = U. The caller classes k are the CALLS targeting an entry of t --
 *             the in-edges of t in the call graph -- plus each entry of t carrying an
 *             OPEN ARRIVAL, a stream that holds a thread exactly as a call does and
 *             that a task can have alongside its callers. Both are structural
 *             neighbours of t, which makes the relation node-local; `termisentry`
 *             says which of the two a term is.
 *
 * `callflow`  X(c) = X(src(c))*y(c), with y(c) the mean number of calls and src(c)
 *             the dispatching activity (the dispatching ENTRY for a forwarding call).
 *
 * `entryflow` X(e) = sum_c X(c) + lambda(e): the requests an entry serves are the
 *             calls reaching it plus its open-arrival stream.
 *
 * `actflow`   X(a) = X(e)*v(a): an activity executes v(a) times per invocation of its
 *             entry. An AND-JOIN is the one place where flow does not add up -- its
 *             target executes once per fork, not once per branch -- so the arcs into
 *             a join are scaled by 1/(number of joined branches).
 *
 * `hostutil`  sum_a X(a)*D(a) = m(h)*U(h), the utilization law at a processor;
 *             m(h)*U(h) is a job count and the factor m(h) drops at an infinite
 *             server.
 *
 * Together these close the system: `little` alone is one equation per task and admits
 * the all-zero solution, so a physics-informed loss built on it should carry the flow
 * and utilization families as well.
 *
 * CONVENTIONS. Rates and populations in a `little` record are PER REPLICA, matching
 * `updateThinkTimes`: X is tput/repl and N is the multiplicity of one copy. Elsewhere
 * throughputs and utilizations are as the solver reports them, totalled over
 * replicas, which is why the server count in `hostutil` is `mult` and NOT
 * `mult*repl`. N(t) is `lqn.mult`; SolverLN iterates on `njobs`, which carries the
 * interlocking corrections and may be `maxmult` under replication, so both are
 * returned per record. S(k) is the entry SERVICE time (phase 1 plus phase 2), the
 * time a thread is held, not the residence time the caller waits for; the difference
 * is the phase-2 tail, flagged by `phase2`.
 *
 * SATURATION. SolverLN clamps the think time at zero, so the `little` equality is an
 * INEQUALITY at a saturated task: when sum_k U(t,k) -> 1 the right-hand side reaches
 * zero and Z can no longer absorb the imbalance. `clamped` marks those once
 * instantiated, and a loss built on these relations should use a one-sided (hinge)
 * form there.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/lqn/lqn_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace lqn {

/** One conservation law, as an aggregation over a node's structural neighbours. */
template <class T>
struct LqnRelation {
    std::string kind;    ///< little | callflow | entryflow | actflow | hostutil
    std::string branch;  ///< the case that produced it: ref/inf/queueing/fwd/arrival, or a call type
    std::size_t target = 0;  ///< absolute index the relation is anchored on
    std::string targetname;
    std::vector<std::size_t> terms;  ///< absolute indices to aggregate over
    /// per term, true for an ENTRY class (open arrival, or a reference task's own
    /// cycle) and false for a CALL class; the two read their rate from different vectors
    std::vector<bool> termisentry;
    std::vector<T> coeff;                        ///< constant coefficient of each term
    T rhsconst = num_traits<T>::from_int(0);     ///< right-hand-side constant
    double mult = 0.0, maxmult = 0.0, repl = 1.0;
    bool scaled = false;      ///< the per-class utilization needs *mult to reach job units
    bool phase2 = false, setup = false;
    bool degenerate = false;  ///< unusable as a residual: infinite mult, no terms, zero call count
    bool clamped = false;     ///< instantiated and saturated: the equality is unattainable
    std::string text;
    /// NaN until instantiated. `has_residual` is false when no solution was supplied,
    /// or when this relation needs a datum the solution does not carry.
    bool has_residual = false;
    T lhs = num_traits<T>::from_int(0), rhs = num_traits<T>::from_int(0);
    T residual = num_traits<T>::from_int(0), relresidual = num_traits<T>::from_int(0);
    std::vector<T> perclassutil;
};

/**
 * The iterates of a solved layered model, indexed by ABSOLUTE element index.
 *
 * The five vectors mirror the same-named members of the layered solver. `un` is the
 * REPORTED utilization, which is not the same quantity as `util`: `util` holds a
 * task's utilization as a SERVER in its own task layer -- the U that closes the
 * thread-pool cycle -- and is left at zero on a host. The `hostutil` family needs the
 * reported one, so it is supplied separately rather than read off the solver:
 * `getEnsembleAvg` re-enters the fixed point in every codebase, and a diagnostic must
 * not re-run one as a side effect of being asked a question. A `hostutil` record
 * whose host has no `un` entry is emitted symbolically with `has_residual` false.
 */
template <class T>
struct LqnSolution {
    std::vector<T> tput, util, thinkt, servt, residt, un;
};

/** The relation set of one layered model. */
template <class T>
struct LqnBalanceEquations {
    std::vector<LqnRelation<T>> eqs;
    std::vector<T> visits;  ///< (nidx+1) executions of each activity per invocation of its entry
    std::vector<std::string> convention;
    Matrix<T> A_little;  ///< (ntasks+1, ncalls+1) 1 where call c is a caller class of task t
    Matrix<T> A_flow;    ///< (ncalls+1, nidx+1) call-to-source incidence weighted by y(c)
    Matrix<T> A_host;    ///< (nhosts+1, nidx+1) host-to-activity incidence weighted by D(a)
    bool has_maxresidual = false;
    T maxresidual = num_traits<T>::from_int(0);
    std::vector<std::string> text;

    std::string str() const {
        std::ostringstream os;
        for (std::size_t i = 0; i < text.size(); ++i) {
            if (i) os << '\n';
            os << text[i];
        }
        return os.str();
    }
};

namespace detail {

/** A number as the reference prints it: an integer bare, Inf as "Inf". */
template <class T>
inline std::string lqn_num(const T& x) {
    const double d = num_traits<T>::to_double(x);
    if (std::isinf(d)) return d > 0 ? "Inf" : "-Inf";
    if (std::isnan(d)) return "NaN";
    if (d == std::floor(d) && std::fabs(d) < 1e15) {
        std::ostringstream os;
        os << static_cast<long long>(d);
        return os.str();
    }
    std::ostringstream os;
    os.precision(6);
    os << d;
    return os.str();
}

inline std::string lqn_num_d(double d) {
    if (std::isinf(d)) return d > 0 ? "Inf" : "-Inf";
    if (std::isnan(d)) return "NaN";
    if (d == std::floor(d) && std::fabs(d) < 1e15) {
        std::ostringstream os;
        os << static_cast<long long>(d);
        return os.str();
    }
    std::ostringstream os;
    os.precision(6);
    os << d;
    return os.str();
}

template <class T>
inline std::string lqn_elem_name(const LqnStruct<T>& lqn, std::size_t idx) {
    if (idx < lqn.hashnames.size() && !lqn.hashnames[idx].empty()) return lqn.hashnames[idx];
    std::ostringstream os;
    os << '#' << idx;
    return os.str();
}

template <class T>
inline std::string lqn_call_name(const LqnStruct<T>& lqn, std::size_t cidx) {
    if (cidx < lqn.callhashnames.size() && !lqn.callhashnames[cidx].empty())
        return lqn.callhashnames[cidx];
    std::ostringstream os;
    os << 'C' << cidx;
    return os.str();
}

template <class T>
inline std::string lqn_term_name(const LqnStruct<T>& lqn, bool isentry, std::size_t idx) {
    return isentry ? lqn_elem_name(lqn, idx) : lqn_call_name(lqn, idx);
}

/** Element of a solution vector, zero where the vector does not reach. */
template <class T>
inline T lqn_sol_at(const std::vector<T>& v, std::size_t i) {
    if (i < v.size()) {
        const double d = num_traits<T>::to_double(v[i]);
        if (!std::isnan(d)) return v[i];
    }
    return num_traits<T>::from_int(0);
}

/** Declared multiplicity, or an out-of-range index as a NaN-like sentinel. */
template <class T>
inline double lqn_mult_of(const std::vector<double>& v, std::size_t i) {
    return i < v.size() ? v[i] : std::numeric_limits<double>::quiet_NaN();
}

/**
 * The precedence graph with the arcs into every AND-join target divided by the number
 * of branches the join waits for.
 *
 * An AND-JOIN is the one place where flow does not add up: its target executes ONCE
 * per fork, not once per branch, so summing the inbound arcs would count it as many
 * times as there are branches. Dividing recovers the rate of one branch exactly when
 * the branches carry equal rate, the case for a well-formed fork/join block.
 * `actpretype` marks the joined PREDECESSORS, so a join target is any successor of one.
 */
template <class T>
inline SparseGraph<T> lqn_join_scaled_graph(const LqnStruct<T>& lqn) {
    SparseGraph<T> G = lqn.graph;
    std::vector<std::size_t> andpre;
    for (std::size_t i = 1; i < lqn.actpretype.size() && i <= G.n; ++i)
        if (lqn.actpretype[i] == lang::PrecedenceType::PRE_AND) andpre.push_back(i);
    if (andpre.empty()) return G;
    const T zero = num_traits<T>::from_int(0);
    std::map<std::size_t, std::vector<std::size_t>> joined;
    for (std::size_t i : andpre)
        for (const std::pair<std::size_t, T>& e : G.row[i])
            if (e.second != zero) joined[e.first].push_back(i);
    for (std::map<std::size_t, std::vector<std::size_t>>::const_iterator it = joined.begin();
         it != joined.end(); ++it) {
        if (it->second.size() < 2) continue;
        const T k = num_traits<T>::from_int(static_cast<int>(it->second.size()));
        for (std::size_t i : it->second) G.set(i, it->first, T(G.get(i, it->first) / k));
    }
    return G;
}

/**
 * Expected executions of every activity per invocation of its entry.
 *
 * The activity precedence arcs of a task are a transient Markov chain whose absorbing
 * state is the reply, so the visit counts of the block solve v = e0*(I-P)^-1: a loop
 * back-edge of weight 1-1/count returns count, and an AND-fork row summing above one
 * returns the branching expectation. Call arcs leave the block and drop out.
 */
template <class T>
inline std::vector<T> lqn_act_visits(const LqnStruct<T>& lqn) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<T> v(lqn.nidx + 1, zero);
    const SparseGraph<T> G = lqn_join_scaled_graph(lqn);
    for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
        const std::size_t tidx = lqn.tshift + t;
        if (tidx >= lqn.actsof.size()) continue;
        std::vector<std::size_t> A = lqn.actsof[tidx];
        std::sort(A.begin(), A.end());
        A.erase(std::unique(A.begin(), A.end()), A.end());
        if (A.empty()) continue;
        std::map<std::size_t, std::size_t> pos;
        for (std::size_t i = 0; i < A.size(); ++i) pos[A[i]] = i;
        // (I-P) transposed, so that solving yields the row vector v as a column
        Matrix<T> M(A.size(), A.size(), zero);
        for (std::size_t i = 0; i < A.size(); ++i)
            for (std::size_t j = 0; j < A.size(); ++j)
                M(j, i) = T((i == j ? one : zero) - (A[i] <= G.n ? G.get(A[i], A[j]) : zero));
        for (std::size_t eidx : lqn.entriesof[tidx]) {
            std::vector<T> e0(A.size(), zero);
            bool any = false;
            if (eidx <= G.n) {
                for (const std::pair<std::size_t, T>& e : G.row[eidx]) {
                    std::map<std::size_t, std::size_t>::const_iterator it = pos.find(e.first);
                    if (it != pos.end() && e.second != zero) {
                        e0[it->second] = one;
                        any = true;
                    }
                }
            }
            if (!any) continue;
            const std::vector<T> x = line::solve(M, e0);
            for (std::size_t i = 0; i < A.size(); ++i) v[A[i]] = T(v[A[i]] + x[i]);
        }
    }
    return v;
}

template <class T>
inline T lqn_arrival_rate(const LqnStruct<T>& lqn, std::size_t eidx) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (eidx >= lqn.has_arrival.size() || !lqn.has_arrival[eidx]) return zero;
    if (eidx >= lqn.arrival.size() || lqn.arrival[eidx].disabled) return zero;
    const T m = lqn.arrival[eidx].mean;
    const double d = num_traits<T>::to_double(m);
    if (std::isfinite(d) && d > 1e-12) return T(one / m);
    return zero;
}

/** Calls targeting each task, i.e. the caller classes of its thread pool. */
template <class T>
inline std::map<std::size_t, std::vector<std::size_t>> lqn_calls_into(const LqnStruct<T>& lqn) {
    std::map<std::size_t, std::vector<std::size_t>> out;
    for (std::size_t cidx = 1; cidx <= lqn.ncalls; ++cidx) {
        const std::size_t dst = lqn.callpair_dst[cidx];
        if (dst < 1 || dst > lqn.nidx) continue;
        out[lqn.parent[dst]].push_back(cidx);
    }
    return out;
}

template <class T>
inline std::vector<std::size_t> lqn_incoming_calls(const LqnStruct<T>& lqn, std::size_t eidx) {
    std::vector<std::size_t> inc;
    for (std::size_t cidx = 1; cidx <= lqn.ncalls; ++cidx)
        if (lqn.callpair_dst[cidx] == eidx) inc.push_back(cidx);
    return inc;
}

/**
 * Entry an activity belongs to: the one whose activity block contains it. An activity
 * shared by two entries is attributed to the first, matching how lqn_act_visits
 * accumulates its visit counts.
 */
template <class T>
inline std::size_t lqn_entry_of_activity(const LqnStruct<T>& lqn, std::size_t aidx) {
    const std::size_t tidx = lqn.parent[aidx];
    if (tidx >= lqn.entriesof.size()) return 0;
    for (std::size_t eidx : lqn.entriesof[tidx]) {
        if (eidx >= lqn.actsof.size()) continue;
        const std::vector<std::size_t>& as = lqn.actsof[eidx];
        if (std::find(as.begin(), as.end(), aidx) != as.end()) return eidx;
    }
    return 0;
}

template <class T>
inline bool lqn_has_phase2(const LqnStruct<T>& lqn, std::size_t tidx) {
    if (tidx >= lqn.entriesof.size()) return false;
    for (std::size_t eidx : lqn.entriesof[tidx]) {
        if (eidx >= lqn.actsof.size()) continue;
        for (std::size_t aidx : lqn.actsof[eidx]) {
            const std::size_t a = aidx - lqn.ashift;
            if (a >= 1 && a < lqn.actphase.size() && lqn.actphase[a] > 1) return true;
        }
    }
    return false;
}

/**
 * Declared think time of a task as it enters the thread cycle: the value for a
 * REFERENCE task, zero for any other. Twin of MATLAB lqn_ref_thinktime.
 */
template <class T>
inline T lqn_ref_thinktime(const LqnStruct<T>& lqn, std::size_t tidx) {
    const T zero = num_traits<T>::from_int(0);
    if (tidx >= lqn.isref.size() || !lqn.isref[tidx]) return zero;
    if (tidx >= lqn.think.size() || lqn.think[tidx].disabled) return zero;
    const double d = num_traits<T>::to_double(lqn.think[tidx].mean);
    if (!std::isfinite(d) || d < 0) return zero;
    return lqn.think[tidx].mean;
}

/** A thread pool's case label and its caller classes. */
struct LqnPool {
    std::string branch;
    std::vector<std::size_t> terms;
    std::vector<bool> isentry;
    bool present = false;
};

/**
 * Which thread-pool case task TIDX falls in, and the caller classes of its pool.
 *
 * The term set is uniform across the cases: every request stream that can hold a
 * thread of the task contributes one class. That is each CALL targeting one of its
 * entries, plus each entry carrying an OPEN ARRIVAL, plus, for a reference task, its
 * own entries, since the cycle of a reference task closes on itself and no layer above
 * drives it. The branch label follows the case analysis of `updateThinkTimes`, which
 * is what decides whether the utilization is a job count (infinite server) or is
 * normalized to [0,1] (every other discipline).
 */
template <class T>
inline LqnPool lqn_thread_pool_branch(const LqnStruct<T>& lqn, std::size_t tidx,
                                      const std::vector<std::size_t>& cin) {
    const T zero = num_traits<T>::from_int(0);
    LqnPool p;
    std::vector<std::size_t> ents;
    if (tidx < lqn.entriesof.size()) ents = lqn.entriesof[tidx];
    std::vector<std::size_t> arv;
    for (std::size_t eidx : ents)
        if (lqn_arrival_rate(lqn, eidx) != zero) arv.push_back(eidx);
    for (std::size_t c : cin) {
        p.terms.push_back(c);
        p.isentry.push_back(false);
    }
    for (std::size_t e : arv) {
        p.terms.push_back(e);
        p.isentry.push_back(true);
    }
    if (tidx < lqn.isref.size() && lqn.isref[tidx]) {
        p.branch = "ref";
        for (std::size_t e : ents)
            if (std::find(arv.begin(), arv.end(), e) == arv.end()) {
                p.terms.push_back(e);
                p.isentry.push_back(true);
            }
        p.present = true;
        return p;
    }
    bool blocking = false;
    for (std::size_t c : cin)
        if (lqn.calltype[c] != lang::CallType::FWD) {
            blocking = true;
            break;
        }
    if (blocking) {
        p.branch = (lqn.sched[tidx] == lang::SchedStrategy::INF) ? "inf" : "queueing";
    } else if (!cin.empty()) {
        p.branch = "fwd";  // a forwarded request holds a thread too
    } else if (!arv.empty()) {
        p.branch = "arrival";
    } else {
        return p;  // no caller, no arrival: no cycle to close
    }
    p.present = true;
    return p;
}

}  // namespace detail

/**
 * Enumerate the conservation laws of the layered model `lqn`.
 *
 * When `sol` is non-null every relation is also instantiated and its residual
 * reported. On a converged fixed point they vanish to the solver's own tolerance; a
 * residual that does not is either a documented convention difference (see `mult`
 * against `maxmult`) or a defect in the solution.
 */
template <class T>
LqnBalanceEquations<T> lqn_balance_equations(const LqnStruct<T>& lqn,
                                            const LqnSolution<T>* sol = nullptr) {
    using namespace detail;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t nidx = lqn.nidx;

    LqnBalanceEquations<T> out;
    out.visits = lqn_act_visits(lqn);
    const std::map<std::size_t, std::vector<std::size_t>> callsinto = lqn_calls_into(lqn);

    // A call throughput is not an iterate of the layered solver: a call inherits the
    // rate of its dispatching element scaled by the mean call count.
    std::vector<T> calltput(lqn.ncalls + 1, zero);
    if (sol) {
        for (std::size_t cidx = 1; cidx <= lqn.ncalls; ++cidx) {
            const std::size_t src = lqn.callpair_src[cidx];
            if (src >= 1 && src <= nidx)
                calltput[cidx] = T(lqn_sol_at(sol->tput, src) * lqn.callproc_mean[cidx]);
        }
    }

    // ---------------- Little's law on each thread pool -----------------------
    for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
        const std::size_t tidx = lqn.tshift + t;
        std::vector<std::size_t> cin;
        std::map<std::size_t, std::vector<std::size_t>>::const_iterator ci = callsinto.find(tidx);
        if (ci != callsinto.end()) cin = ci->second;
        const LqnPool pool = lqn_thread_pool_branch(lqn, tidx, cin);
        if (!pool.present) continue;

        LqnRelation<T> r;
        r.kind = "little";
        r.branch = pool.branch;
        r.target = tidx;
        r.targetname = lqn_elem_name(lqn, tidx);
        r.terms = pool.terms;
        r.termisentry = pool.isentry;
        r.coeff.assign(pool.terms.size(), one);
        r.mult = lqn_mult_of<T>(lqn.mult, tidx);
        r.maxmult = lqn_mult_of<T>(lqn.maxmult, tidx);
        r.repl = (tidx < lqn.repl.size() && lqn.repl[tidx] > 1.0) ? lqn.repl[tidx] : 1.0;
        r.rhsconst = num_traits<T>::from_double(r.mult);
        r.scaled = pool.branch != "inf";
        r.phase2 = lqn_has_phase2(lqn, tidx);
        r.setup = tidx < lqn.hassetup.size() && lqn.hassetup[tidx];
        r.degenerate = !std::isfinite(r.mult) || pool.terms.empty();

        const T z = lqn_ref_thinktime(lqn, tidx);
        {
            std::ostringstream lhs;
            lhs << "X(" << r.targetname << ")*(Z(" << r.targetname << ")";
            if (num_traits<T>::to_double(z) > 0) lhs << " + " << lqn_num(z);
            lhs << ")";
            for (std::size_t i = 0; i < r.terms.size(); ++i)
                lhs << " + B(" << lqn_term_name(lqn, r.termisentry[i], r.terms[i]) << ")";
            std::ostringstream txt;
            txt << lhs.str() << " = " << lqn_num_d(r.mult);
            if (r.scaled && !r.terms.empty()) {
                txt << "\n           equivalently  X(" << r.targetname << ")*Z(" << r.targetname
                    << ") = " << lqn_num_d(r.mult) << "*(1";
                for (std::size_t i = 0; i < r.terms.size(); ++i)
                    txt << " - U(" << r.targetname << ","
                        << lqn_term_name(lqn, r.termisentry[i], r.terms[i]) << ")";
                txt << ")";
            }
            r.text = txt.str();
        }

        if (sol) {
            const T X = T(lqn_sol_at(sol->tput, tidx) / num_traits<T>::from_double(r.repl));
            T sumB = zero;
            r.perclassutil.clear();
            const bool jobunits = lqn.sched[tidx] == lang::SchedStrategy::INF ||
                                  !std::isfinite(r.mult) || r.mult <= 0;
            for (std::size_t i = 0; i < r.terms.size(); ++i) {
                T b;
                if (r.termisentry[i]) {
                    b = T(lqn_sol_at(sol->tput, r.terms[i]) * lqn_sol_at(sol->servt, r.terms[i]));
                } else {
                    b = T(calltput[r.terms[i]] *
                          lqn_sol_at(sol->servt, lqn.callpair_dst[r.terms[i]]));
                }
                sumB = T(sumB + b);
                r.perclassutil.push_back(jobunits ? b : T(b / num_traits<T>::from_double(r.mult)));
            }
            const T zt = lqn_sol_at(sol->thinkt, tidx);
            r.lhs = T(X * T(zt + z) + sumB);
            r.rhs = r.rhsconst;
            r.residual = T(r.lhs - r.rhs);
            const double den = std::max(std::fabs(num_traits<T>::to_double(r.rhs)), 1e-12);
            r.relresidual = T(r.residual / num_traits<T>::from_double(den));
            r.clamped = std::isfinite(r.mult) &&
                        num_traits<T>::to_double(T(r.rhsconst - sumB - T(X * z))) < 0;
            r.has_residual = true;
        }
        out.eqs.push_back(r);
    }

    // ---------------- flow conservation --------------------------------------
    for (std::size_t cidx = 1; cidx <= lqn.ncalls; ++cidx) {
        LqnRelation<T> r;
        r.kind = "callflow";
        const lang::CallType ct = lqn.calltype[cidx];
        r.branch = ct == lang::CallType::SYNC    ? "sync"
                   : ct == lang::CallType::ASYNC ? "async"
                   : ct == lang::CallType::FWD   ? "fwd"
                                                 : "";
        r.target = lqn.cshift + cidx;
        r.targetname = lqn_call_name(lqn, cidx);
        const std::size_t src = lqn.callpair_src[cidx];
        const T y = lqn.callproc_mean[cidx];
        r.terms.push_back(src);
        r.coeff.push_back(y);
        r.degenerate = (y == zero);
        std::ostringstream txt;
        txt << "X(" << r.targetname << ") = X(" << lqn_elem_name(lqn, src) << ") * " << lqn_num(y);
        r.text = txt.str();
        if (sol) {
            r.lhs = calltput[cidx];
            r.rhs = T(lqn_sol_at(sol->tput, src) * y);
            r.residual = T(r.lhs - r.rhs);
            const double den = std::max(std::fabs(num_traits<T>::to_double(r.rhs)), 1e-12);
            r.relresidual = T(r.residual / num_traits<T>::from_double(den));
            r.has_residual = true;
        }
        out.eqs.push_back(r);
    }

    for (std::size_t e = 1; e <= lqn.nentries; ++e) {
        const std::size_t eidx = lqn.eshift + e;
        const std::vector<std::size_t> inc = lqn_incoming_calls(lqn, eidx);
        const T lam = lqn_arrival_rate(lqn, eidx);
        if (inc.empty() && lam == zero) continue;  // a reference entry is driven by its own cycle
        LqnRelation<T> r;
        r.kind = "entryflow";
        r.target = eidx;
        r.targetname = lqn_elem_name(lqn, eidx);
        r.terms = inc;
        r.coeff.assign(inc.size(), one);
        r.rhsconst = lam;
        std::ostringstream txt;
        txt << "X(" << r.targetname << ") =";
        for (std::size_t i = 0; i < inc.size(); ++i)
            txt << (i == 0 ? " X(" : " + X(") << lqn_call_name(lqn, inc[i]) << ")";
        if (lam != zero)
            txt << (inc.empty() ? " " : " + ") << lqn_num(lam) << "   (open arrival)";
        r.text = txt.str();
        if (sol) {
            r.lhs = lqn_sol_at(sol->tput, eidx);
            T rhs = lam;
            for (std::size_t c : inc) rhs = T(rhs + calltput[c]);
            r.rhs = rhs;
            r.residual = T(r.lhs - r.rhs);
            const double den = std::max(std::fabs(num_traits<T>::to_double(r.rhs)), 1e-12);
            r.relresidual = T(r.residual / num_traits<T>::from_double(den));
            r.has_residual = true;
        }
        out.eqs.push_back(r);
    }

    for (std::size_t a = 1; a <= lqn.nacts; ++a) {
        const std::size_t aidx = lqn.ashift + a;
        const std::size_t eidx = lqn_entry_of_activity(lqn, aidx);
        if (eidx == 0) continue;
        LqnRelation<T> r;
        r.kind = "actflow";
        r.target = aidx;
        r.targetname = lqn_elem_name(lqn, aidx);
        r.terms.push_back(eidx);
        r.coeff.push_back(out.visits[aidx]);
        r.degenerate = !std::isfinite(num_traits<T>::to_double(out.visits[aidx]));
        std::ostringstream txt;
        txt << "X(" << r.targetname << ") = X(" << lqn_elem_name(lqn, eidx) << ") * "
            << lqn_num(out.visits[aidx]);
        r.text = txt.str();
        if (sol) {
            r.lhs = lqn_sol_at(sol->tput, aidx);
            r.rhs = T(lqn_sol_at(sol->tput, eidx) * out.visits[aidx]);
            r.residual = T(r.lhs - r.rhs);
            const double den = std::max(std::fabs(num_traits<T>::to_double(r.rhs)), 1e-12);
            r.relresidual = T(r.residual / num_traits<T>::from_double(den));
            r.has_residual = true;
        }
        out.eqs.push_back(r);
    }

    // ---------------- utilization law at each processor ----------------------
    for (std::size_t h = 1; h <= lqn.nhosts; ++h) {
        const std::size_t hidx = lqn.hshift + h;
        LqnRelation<T> r;
        r.kind = "hostutil";
        r.target = hidx;
        r.targetname = lqn_elem_name(lqn, hidx);
        r.mult = lqn_mult_of<T>(lqn.mult, hidx);
        r.repl = (hidx < lqn.repl.size() && lqn.repl[hidx] > 1.0) ? lqn.repl[hidx] : 1.0;
        r.scaled = lqn.sched[hidx] != lang::SchedStrategy::INF;
        if (hidx < lqn.tasksof.size()) {
            for (std::size_t tidx : lqn.tasksof[hidx]) {
                if (tidx >= lqn.actsof.size()) continue;
                for (std::size_t aidx : lqn.actsof[tidx]) {
                    if (aidx >= lqn.hostdem.size() || lqn.hostdem[aidx].disabled) continue;
                    const T d = lqn.hostdem[aidx].mean;
                    if (d == zero) continue;
                    r.terms.push_back(aidx);
                    r.coeff.push_back(d);
                }
            }
        }
        // The server count is the declared multiplicity ALONE, not mult*repl: a
        // replicated host reports its throughputs and its utilization as TOTALS over
        // the copies, so the extra factor would double-count the replication.
        double m = r.mult;
        if (!r.scaled || !std::isfinite(m)) m = 1.0;
        r.rhsconst = num_traits<T>::from_double(m);
        r.branch = r.scaled ? "queueing" : "inf";
        r.degenerate = r.terms.empty();
        std::ostringstream txt;
        for (std::size_t i = 0; i < r.terms.size(); ++i) {
            if (i) txt << " + ";
            txt << "X(" << lqn_elem_name(lqn, r.terms[i]) << ")*" << lqn_num(r.coeff[i]);
        }
        std::string lhs = txt.str();
        if (lhs.empty()) lhs = "0";
        std::ostringstream full;
        full << lhs << " = " << lqn_num_d(m) << "*U(" << r.targetname << ")";
        r.text = full.str();
        const bool haveun =
            sol && hidx < sol->un.size() && !std::isnan(num_traits<T>::to_double(sol->un[hidx]));
        if (haveun) {
            T v = zero;
            for (std::size_t i = 0; i < r.terms.size(); ++i)
                v = T(v + lqn_sol_at(sol->tput, r.terms[i]) * r.coeff[i]);
            r.lhs = v;
            r.rhs = T(r.rhsconst * sol->un[hidx]);
            r.residual = T(r.lhs - r.rhs);
            const double den = std::max(std::fabs(num_traits<T>::to_double(r.rhs)), 1e-12);
            r.relresidual = T(r.residual / num_traits<T>::from_double(den));
            r.has_residual = true;
        }
        out.eqs.push_back(r);
    }

    // ---------------- aggregation incidences ---------------------------------
    out.A_little = Matrix<T>(lqn.ntasks + 1, lqn.ncalls + 1, zero);
    out.A_flow = Matrix<T>(lqn.ncalls + 1, nidx + 1, zero);
    out.A_host = Matrix<T>(lqn.nhosts + 1, nidx + 1, zero);
    for (const LqnRelation<T>& r : out.eqs) {
        if (r.kind == "little") {
            // the call classes only; an entry class (open arrival, or the self-driven
            // cycle of a reference task) is not a call
            for (std::size_t i = 0; i < r.terms.size(); ++i)
                if (!r.termisentry[i]) out.A_little(r.target - lqn.tshift, r.terms[i]) = one;
        } else if (r.kind == "callflow") {
            out.A_flow(r.target - lqn.cshift, r.terms[0]) = r.coeff[0];
        } else if (r.kind == "hostutil") {
            for (std::size_t i = 0; i < r.terms.size(); ++i)
                out.A_host(r.target - lqn.hshift, r.terms[i]) = r.coeff[i];
        }
    }

    for (const LqnRelation<T>& r : out.eqs) {
        if (r.degenerate || !r.has_residual) continue;
        const T v = num_traits<T>::from_double(std::fabs(num_traits<T>::to_double(r.residual)));
        if (!out.has_maxresidual || v > out.maxresidual) {
            out.maxresidual = v;
            out.has_maxresidual = true;
        }
    }

    out.convention.push_back("Conventions:");
    out.convention.push_back(
        "  little   : rates and populations PER REPLICA (X = tput/repl, N = mult of one copy).");
    out.convention.push_back(
        "             B(t,k) is the per-class utilization in JOB units; for a queueing task");
    out.convention.push_back(
        "             B = mult*U with U in [0,1], at an infinite server B = U directly.");
    out.convention.push_back(
        "             S(k) is the entry SERVICE time (phase 1 + phase 2), the thread hold time.");
    out.convention.push_back(
        "  other    : throughputs and utilizations as the solver reports them, totalled over "
        "replicas.");
    out.convention.push_back(
        "  N(t)     : lqn.mult. SolverLN iterates on njobs (interlocking corrections, maxmult "
        "under");
    out.convention.push_back("             replication), reported per record as mult/maxmult.");

    // ---------------- report --------------------------------------------------
    {
        std::ostringstream hdr;
        hdr << "LQN balance equations: " << lqn.nhosts << " hosts, " << lqn.ntasks << " tasks, "
            << lqn.nentries << " entries, " << lqn.nacts << " activities, " << lqn.ncalls
            << " calls";
        out.text.push_back(hdr.str());
    }
    for (const std::string& c : out.convention) out.text.push_back(c);
    const char* kinds[5] = {"little", "callflow", "entryflow", "actflow", "hostutil"};
    const char* titles[5] = {"thread-pool Little's law", "call-flow balance", "entry-flow balance",
                             "activity-flow balance", "host utilization law"};
    for (int k = 0; k < 5; ++k) {
        std::vector<std::size_t> sel;
        for (std::size_t i = 0; i < out.eqs.size(); ++i)
            if (out.eqs[i].kind == kinds[k]) sel.push_back(i);
        if (sel.empty()) continue;
        out.text.push_back("");
        {
            std::ostringstream os;
            os << "--- " << titles[k] << "  (kind='" << kinds[k] << "', " << sel.size()
               << " relations) ---";
            out.text.push_back(os.str());
        }
        for (std::size_t i : sel) {
            const LqnRelation<T>& r = out.eqs[i];
            std::ostringstream head;
            head << "[" << (i + 1) << "] " << r.targetname;
            if (r.kind == "little" || r.kind == "hostutil")
                head << " " << r.branch << " mult=" << lqn_num_d(r.mult)
                     << " repl=" << lqn_num_d(r.repl);
            else if (!r.branch.empty())
                head << " " << r.branch;
            out.text.push_back(head.str());
            std::string body = r.text;
            std::size_t start = 0;
            while (start <= body.size()) {
                const std::size_t nl = body.find('\n', start);
                const std::string line =
                    body.substr(start, nl == std::string::npos ? std::string::npos : nl - start);
                out.text.push_back("       " + line);
                if (nl == std::string::npos) break;
                start = nl + 1;
            }
            std::vector<std::string> notes;
            if (r.phase2) notes.push_back("phase-2 tail on Z");
            if (r.setup) notes.push_back("setup charge on Z");
            if (r.degenerate) notes.push_back("DEGENERATE (not usable as a residual)");
            if (r.clamped) notes.push_back("SATURATED (equality unattainable, use a hinge)");
            if (sol && !r.degenerate && !r.has_residual)
                notes.push_back("not instantiated (pass un for the host utilization law)");
            if (!notes.empty()) {
                std::ostringstream os;
                os << "       note: ";
                for (std::size_t j = 0; j < notes.size(); ++j) {
                    if (j) os << ", ";
                    os << notes[j];
                }
                out.text.push_back(os.str());
            }
            if (r.has_residual && !r.degenerate) {
                std::ostringstream os;
                os.precision(6);
                os << "       lhs=" << num_traits<T>::to_double(r.lhs)
                   << "  rhs=" << num_traits<T>::to_double(r.rhs)
                   << "  residual=" << num_traits<T>::to_double(r.residual)
                   << "  rel=" << num_traits<T>::to_double(r.relresidual);
                out.text.push_back(os.str());
            }
        }
    }
    if (sol) {
        std::size_t nd = 0;
        for (const LqnRelation<T>& r : out.eqs)
            if (!r.degenerate) ++nd;
        out.text.push_back("");
        std::ostringstream os;
        os.precision(3);
        os << "max |residual| over " << nd << " non-degenerate relations: "
           << (out.has_maxresidual ? num_traits<T>::to_double(out.maxresidual)
                                   : std::numeric_limits<double>::quiet_NaN());
        out.text.push_back(os.str());
    }

    return out;
}

}  // namespace lqn
}  // namespace line

#endif  // LINE_API_LQN_LQN_BALANCE_EQUATIONS_H
