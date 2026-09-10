/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_WORKFLOW_WORKFLOW_H
#define LINE_LANG_WORKFLOW_WORKFLOW_H

/**
 * An activity workflow reduced to one phase-type law.
 *
 * Port of matlab/src/lang/workflow/Workflow.m and WorkflowActivity.m, and of
 * their twins jline.lang.workflow.Workflow (JAR) and
 * line_solver.lang.workflow.Workflow (Python). A workflow is a precedence
 * graph over activities, each carrying a host demand; `to_ph` composes those
 * laws into a single phase-type distribution.
 *
 * WHAT THIS IS NOT. `include/line/api/wf/` holds the workflow pattern
 * DETECTORS of `jline.api.wf`, which read a link matrix and report sequences,
 * branches, parallel blocks and loops. This header is the algebra those
 * patterns feed: the composition that turns a precedence graph into an
 * (alpha, T) pair. An LQN fork-join is also a different object -- its branches
 * contend for a host, so it yields a response time under contention rather
 * than the order statistic of independent activity times returned here.
 *
 * A precedence graph that is SERIES-PARALLEL is reduced exactly, by recursive
 * composition of the series-parallel tree, which handles arbitrary nesting (a
 * fork inside a loop, a branch that is itself a fork-join). Any other graph
 * falls back to the block composition, which is a heuristic and is documented
 * as one.
 *
 * A LOOP repeats its body a GEOMETRIC number of times of mean COUNT, the
 * POST_LOOP semantics of an activity graph (`getStruct.m:528-567`): a count of
 * at least one runs the body once and takes the back edge with probability
 * 1-1/COUNT, and a fractional count runs the body at most once, with
 * probability COUNT. The COUNT-fold convolution is a different law with the
 * same mean; it stays available as `compose_repeat` for a caller that wants
 * exactly that.
 *
 * An EXTERNAL CALL has no separate representation: an activity whose host
 * demand is the law of the call response time composes exactly like a local
 * computation, and an asynchronous call blocks the caller for no time and is
 * simply left out of the workflow.
 *
 * A QUORUM (partial) AND-join is REFUSED by name rather than served as a full
 * join, because the first k of n branches to finish is not their maximum.
 */

#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "line/lang/dist_fitters.h"
#include "line/lang/dist_scale_rate.h"
#include "line/lang/distribution.h"
#include "line/lang/lang_types.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace workflow {

using lang::Distrib;
using lang::GlobalConstants;
// The precedence kinds of `ActivityPrecedenceType.m`, shared with the LQN layer
using lang::PrecedenceType;
using lang::ProcessType;

/**
 * One precedence of the activity graph.
 *
 * `pre_params` carries the AND-join quorum when there is one; `post_params`
 * carries the OR-fork probabilities or the loop count.
 */
template <class T>
struct Precedence {
    std::vector<std::string> pre_acts;
    std::vector<std::string> post_acts;
    PrecedenceType pre_type = PrecedenceType::PRE_SEQ;
    PrecedenceType post_type = PrecedenceType::POST_SEQ;
    std::vector<T> pre_params;
    std::vector<T> post_params;
};

/** A phase-type law as the composition rules pass it around. */
template <class T>
struct PhLaw {
    std::vector<T> alpha;
    Matrix<T> S;
};

/**
 * A computational activity.
 *
 * Carries no call list: an external call is an activity whose host demand is
 * the law of the call response time.
 */
template <class T>
class WorkflowActivity {
public:
    WorkflowActivity() = default;
    WorkflowActivity(const std::string& name, const Distrib<T>& host_demand)
        : name_(name), host_demand_(host_demand) {}

    const std::string& name() const { return name_; }

    const Distrib<T>& host_demand() const { return host_demand_; }
    void set_host_demand(const Distrib<T>& d) { host_demand_ = d; }

    T host_demand_mean() const { return host_demand_.mean; }
    T host_demand_scv() const { return host_demand_.scv; }

    /**
     * The (alpha, T) pair of this activity.
     *
     * A zero-time activity is one immediate phase, the representation used
     * throughout this class; a Markovian law hands back its own pair; anything
     * else is fitted to an acyclic phase-type on its first two moments, as the
     * MATLAB and JAR twins do.
     */
    PhLaw<T> ph_representation() const {
        const T zero = num_traits<T>::from_int(0);
        PhLaw<T> out;

        if (host_demand_.type == ProcessType::IMMEDIATE ||
            !(host_demand_.mean > num_traits<T>::from_double(GlobalConstants::FineTol))) {
            out.alpha.assign(1, num_traits<T>::from_int(1));
            out.S = Matrix<T>(1, 1, zero);
            out.S(0, 0) = T(-num_traits<T>::from_double(GlobalConstants::Immediate));
            return out;
        }

        if (host_demand_.has_map()) {
            out.S = host_demand_.D0;
            out.alpha = init_prob_of(host_demand_);
            return out;
        }

        T scv = host_demand_.scv;
        if (!(scv > num_traits<T>::from_double(GlobalConstants::FineTol)))
            scv = num_traits<T>::from_int(1);
        const Distrib<T> aph = lang::aph_fit_mean_scv(host_demand_.mean, scv);
        out.S = aph.D0;
        out.alpha = init_prob_of(aph);
        return out;
    }

    std::size_t num_phases() const { return ph_representation().S.rows(); }

private:
    /**
     * The initial vector of a Markovian law.
     *
     * `Distrib::phase_type` stores alpha in `params`, so a PH/APH/ME hands it
     * back directly. For every other Markovian family D1 = (-D0 e) alpha, so
     * ANY row with a positive exit rate recovers alpha -- and it cannot be row
     * 0, since a canonical bidiagonal APH never completes from phase 1.
     */
    static std::vector<T> init_prob_of(const Distrib<T>& d) {
        const std::size_t n = d.D0.rows();
        const T zero = num_traits<T>::from_int(0);
        if ((d.type == ProcessType::PH || d.type == ProcessType::APH ||
             d.type == ProcessType::ME) &&
            d.params.size() == n) {
            return std::vector<T>(d.params.begin(), d.params.begin() + static_cast<long>(n));
        }
        for (std::size_t i = 0; i < n; ++i) {
            T out = zero;
            for (std::size_t j = 0; j < n; ++j) out += d.D1(i, j);
            if (!(out > zero)) continue;
            std::vector<T> alpha(n);
            for (std::size_t j = 0; j < n; ++j) alpha[j] = T(d.D1(i, j) / out);
            return alpha;
        }
        // No phase completes: start in phase 0, which is what a degenerate
        // representation leaves as the only defensible reading
        std::vector<T> alpha(n, zero);
        if (n > 0) alpha[0] = num_traits<T>::from_int(1);
        return alpha;
    }

    std::string name_;
    Distrib<T> host_demand_;
};

/** A node of the series-parallel tree. */
enum class SPNodeType { LEAF, SERIAL, PAR, OR, LOOP };

template <class T>
struct SPNode {
    SPNodeType type = SPNodeType::LEAF;
    /** Activity index for a LEAF, npos otherwise. */
    std::size_t act = static_cast<std::size_t>(-1);
    std::vector<std::size_t> kids;
    std::size_t parent = static_cast<std::size_t>(-1);
    /** Branch probabilities for an OR node. */
    std::vector<T> probs;
    /** Loop count for a LOOP node. */
    T count = num_traits<T>::from_int(0);
    PhLaw<T> law;
    bool valid = false;
};

/**
 * The flat series-parallel tree.
 *
 * `execs` carries the expected number of executions of each node per workflow
 * execution, which is the weight by which an LQN metric reconstruction splits
 * a layer result back over entries, activities and calls.
 */
template <class T>
struct SPTree {
    std::vector<SPNode<T>> nodes;
    std::size_t root = static_cast<std::size_t>(-1);
    /** Node index of each activity's leaf, npos when the activity has none. */
    std::vector<std::size_t> leaf_of;
    std::vector<T> execs;
};

template <class T>
class Workflow {
public:
    static constexpr std::size_t npos = static_cast<std::size_t>(-1);

    explicit Workflow(const std::string& name) : name_(name) {}

    const std::string& name() const { return name_; }

    /** Add an activity; the name must be unique. */
    std::size_t add_activity(const std::string& name, const Distrib<T>& host_demand) {
        if (activity_map_.find(name) != activity_map_.end())
            throw InputError("Workflow: activity '" + name + "' is already declared");
        activities_.push_back(WorkflowActivity<T>(name, host_demand));
        const std::size_t idx = activities_.size() - 1;
        activity_map_[name] = idx;
        invalidate_topology();
        return idx;
    }

    /** Add an activity with an exponential host demand of the given mean. */
    std::size_t add_activity(const std::string& name, const T& mean) {
        return add_activity(name, Distrib<T>::exp_mean(mean));
    }

    void add_precedence(const Precedence<T>& prec) {
        precedences_.push_back(prec);
        invalidate_topology();
    }

    std::size_t num_activities() const { return activities_.size(); }
    const std::vector<WorkflowActivity<T>>& activities() const { return activities_; }
    const std::vector<Precedence<T>>& precedences() const { return precedences_; }

    std::size_t activity_index(const std::string& name) const {
        auto it = activity_map_.find(name);
        return it == activity_map_.end() ? npos : it->second;
    }

    WorkflowActivity<T>& activity(const std::string& name) {
        const std::size_t i = activity_index(name);
        if (i == npos) throw InputError("Workflow: activity '" + name + "' not found");
        return activities_[i];
    }

    /** Activity by index, the form a series-parallel LEAF node names it in. */
    const WorkflowActivity<T>& activity_at(std::size_t i) const {
        if (i >= activities_.size()) throw InputError("Workflow: activity index out of range");
        return activities_[i];
    }

    /**
     * Validate the workflow, throwing on the first defect.
     *
     * A quorum AND-join is refused here rather than composed as a full join.
     */
    void validate() const {
        if (activities_.empty())
            throw InputError("Workflow must have at least one activity.");

        for (const Precedence<T>& prec : precedences_) {
            for (const std::string& nm : prec.pre_acts)
                if (activity_index(nm) == npos)
                    throw InputError("Activity '" + nm +
                                     "' referenced in precedence not found in workflow.");
            for (const std::string& nm : prec.post_acts)
                if (activity_index(nm) == npos)
                    throw InputError("Activity '" + nm +
                                     "' referenced in precedence not found in workflow.");
        }

        for (const Precedence<T>& prec : precedences_) {
            if (prec.post_type == PrecedenceType::POST_OR) {
                if (prec.post_params.empty())
                    throw InputError("OR-fork must have probabilities specified.");
                T total = num_traits<T>::from_int(0);
                for (const T& p : prec.post_params) total += p;
                const double gap =
                    std::abs(num_traits<T>::to_double(total) - 1.0);
                if (gap > GlobalConstants::FineTol)
                    throw InputError("OR-fork probabilities must sum to 1.");
            }
            if (prec.post_type == PrecedenceType::POST_LOOP) {
                if (prec.post_params.size() != 1)
                    throw InputError("Loop count must be a single positive number.");
                if (!(prec.post_params[0] > num_traits<T>::from_int(0)))
                    throw InputError("Loop count must be a positive number.");
            }
            if (prec.pre_type == PrecedenceType::PRE_AND && !prec.pre_params.empty()) {
                const double quorum = num_traits<T>::to_double(prec.pre_params[0]);
                const double nb = static_cast<double>(prec.pre_acts.size());
                if (quorum > 0.0 && quorum < nb)
                    throw UnsupportedError(
                        "AND-join with quorum " + std::to_string(static_cast<long>(quorum)) +
                        " of " + std::to_string(static_cast<long>(nb)) +
                        " is not supported by Workflow: a partial join is not the maximum of the "
                        "branches. Use a full join, or SolverLN with method='default', which "
                        "routes the join explicitly.");
            }
        }
    }

    /**
     * The composed law of the workflow.
     *
     * The generator is acyclic unless a geometric loop closes a cycle over a
     * multi-phase body, which `is_acyclic_generator` reports; the returned
     * Distrib is typed APH or PH accordingly.
     */
    Distrib<T> to_ph() {
        if (cached_valid_) return cached_ph_;
        validate();

        PhLaw<T> law;
        if (!compose_series_parallel(law)) {
            // Not series-parallel: fall back to the block composition
            law = build_ctmc();
        }

        cached_ph_ = Distrib<T>::phase_type(law.alpha, law.S, is_acyclic_generator(law.S));
        cached_valid_ = true;
        return cached_ph_;
    }

    /**
     * Recompose the law after a demand change.
     *
     * Only the series-parallel nodes on the path from a dirty leaf to the root
     * recompose; nodes whose subtree is unchanged keep their cached law. The
     * topology is not rebuilt.
     */
    Distrib<T> refresh_ph() { return to_ph(); }

    /**
     * Change the host demand of one activity.
     *
     * Marks only that leaf dirty, which is the entry point an iterative solver
     * uses when it updates call-response laws at each iteration.
     */
    void set_activity_demand(const std::string& name, const Distrib<T>& host_demand) {
        const std::size_t i = activity_index(name);
        if (i == npos) throw InputError("Workflow: activity '" + name + "' not found");
        activities_[i].set_host_demand(host_demand);
        invalidate_activity(i);
    }

    /**
     * Change only the mean of one activity, preserving its shape.
     *
     * The law is scaled in time rather than refitted, so its SCV, its skewness
     * and its order are preserved and the cached tree keeps its shape.
     */
    void set_activity_demand_mean(const std::string& name, const T& mean_value) {
        if (!(mean_value > num_traits<T>::from_int(0)))
            throw InputError("The activity mean must be a positive finite scalar.");
        const std::size_t i = activity_index(name);
        if (i == npos) throw InputError("Workflow: activity '" + name + "' not found");

        const Distrib<T>& d = activities_[i].host_demand();
        const T old_mean = d.mean;
        if (d.type == ProcessType::IMMEDIATE ||
            !(old_mean > num_traits<T>::from_double(GlobalConstants::FineTol))) {
            activities_[i].set_host_demand(Distrib<T>::exp_mean(mean_value));
            invalidate_activity(i);
            return;
        }

        const T factor = T(old_mean / mean_value);
        activities_[i].set_host_demand(lang::dist_scale_rate(d, factor));
        // The SCV is invariant under a time scaling
        rescale_activity_leaf(i, factor);
    }

    /** Discard the cached law and the decomposition. */
    void invalidate_topology() {
        cached_valid_ = false;
        sp_tree_.nodes.clear();
        sp_tree_.root = npos;
        sp_tree_.leaf_of.clear();
        sp_tree_.execs.clear();
        sp_built_ = false;
        sp_failed_ = false;
    }

    /** Mark one activity law dirty, keeping every other cached block. */
    void invalidate_activity(std::size_t act_idx) {
        cached_valid_ = false;
        if (!sp_built_) return;
        if (act_idx >= sp_tree_.leaf_of.size() || sp_tree_.leaf_of[act_idx] == npos) {
            // Activity outside the decomposition: rebuild it entirely
            invalidate_topology();
            return;
        }
        invalidate_branch(sp_tree_, sp_tree_.leaf_of[act_idx]);
    }

    /**
     * Time-scale a cached leaf in place: S -> S*FACTOR with alpha fixed.
     *
     * Ancestors still recompose, because they mix phases of several leaves,
     * but the tree keeps its shape and no acyclic phase-type is refitted.
     */
    void rescale_activity_leaf(std::size_t act_idx, const T& factor) {
        cached_valid_ = false;
        if (!sp_built_) return;
        if (act_idx >= sp_tree_.leaf_of.size() || sp_tree_.leaf_of[act_idx] == npos) {
            invalidate_topology();
            return;
        }
        const std::size_t k = sp_tree_.leaf_of[act_idx];
        invalidate_branch(sp_tree_, k);
        SPNode<T>& node = sp_tree_.nodes[k];
        if (node.law.S.rows() > 0) {
            for (std::size_t i = 0; i < node.law.S.rows(); ++i)
                for (std::size_t j = 0; j < node.law.S.cols(); ++j)
                    node.law.S(i, j) = T(node.law.S(i, j) * factor);
            node.valid = true;
        }
    }

    /**
     * The cached series-parallel decomposition, or null when the precedence
     * graph is not series-parallel.
     */
    const SPTree<T>* sp_tree() {
        if (!sp_built_ && !sp_failed_) build_sp_tree();
        return sp_built_ ? &sp_tree_ : nullptr;
    }

    // -----------------------------------------------------------------
    // Composition rules
    // -----------------------------------------------------------------

    /**
     * Serial composition: the second law starts when the first absorbs.
     *
     *   S = [S1, (-S1 e) alpha2; 0, S2]
     *
     * A defective alpha1 carries an atom at zero, which starts the second law
     * immediately; this is `aph_simplify` pattern 1.
     */
    static PhLaw<T> compose_serial(const PhLaw<T>& a, const PhLaw<T>& b) {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t n1 = a.S.rows(), n2 = b.S.rows();

        PhLaw<T> out;
        out.S = Matrix<T>(n1 + n2, n1 + n2, zero);
        for (std::size_t i = 0; i < n1; ++i)
            for (std::size_t j = 0; j < n1; ++j) out.S(i, j) = a.S(i, j);
        for (std::size_t i = 0; i < n2; ++i)
            for (std::size_t j = 0; j < n2; ++j) out.S(n1 + i, n1 + j) = b.S(i, j);

        for (std::size_t i = 0; i < n1; ++i) {
            T rate = zero;
            for (std::size_t j = 0; j < n1; ++j) rate += a.S(i, j);
            rate = T(-rate);
            for (std::size_t j = 0; j < n2; ++j) out.S(i, n1 + j) = T(rate * b.alpha[j]);
        }

        T defect = num_traits<T>::from_int(1);
        for (const T& v : a.alpha) defect -= v;
        out.alpha.assign(n1 + n2, zero);
        for (std::size_t i = 0; i < n1; ++i) out.alpha[i] = a.alpha[i];
        for (std::size_t j = 0; j < n2; ++j) out.alpha[n1 + j] = T(defect * b.alpha[j]);
        return out;
    }

    /**
     * Parallel (AND-fork/join) composition: the time until BOTH complete.
     *
     * States (i,j) with both active carry the Kronecker sum S1 (+) S2; when one
     * branch absorbs the chain moves into that branch's own survivor block, so
     * the law is the maximum of the two.
     */
    static PhLaw<T> compose_parallel(const PhLaw<T>& a, const PhLaw<T>& b) {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t n1 = a.S.rows(), n2 = b.S.rows();
        const std::size_t nboth = n1 * n2;
        const std::size_t ntot = nboth + n1 + n2;

        std::vector<T> abs1(n1, zero), abs2(n2, zero);
        for (std::size_t i = 0; i < n1; ++i) {
            T r = zero;
            for (std::size_t j = 0; j < n1; ++j) r += a.S(i, j);
            abs1[i] = T(-r);
        }
        for (std::size_t i = 0; i < n2; ++i) {
            T r = zero;
            for (std::size_t j = 0; j < n2; ++j) r += b.S(i, j);
            abs2[i] = T(-r);
        }

        PhLaw<T> out;
        out.S = Matrix<T>(ntot, ntot, zero);

        // Kronecker sum on the block where both are still running
        for (std::size_t i = 0; i < n1; ++i)
            for (std::size_t j = 0; j < n2; ++j) {
                const std::size_t r = i * n2 + j;
                for (std::size_t ii = 0; ii < n1; ++ii)
                    out.S(r, ii * n2 + j) += a.S(i, ii);
                for (std::size_t jj = 0; jj < n2; ++jj)
                    out.S(r, i * n2 + jj) += b.S(j, jj);
                // one branch absorbs, the other keeps running
                out.S(r, nboth + i) += abs2[j];
                out.S(r, nboth + n1 + j) += abs1[i];
            }

        for (std::size_t i = 0; i < n1; ++i)
            for (std::size_t j = 0; j < n1; ++j) out.S(nboth + i, nboth + j) = a.S(i, j);
        for (std::size_t i = 0; i < n2; ++i)
            for (std::size_t j = 0; j < n2; ++j)
                out.S(nboth + n1 + i, nboth + n1 + j) = b.S(i, j);

        out.alpha.assign(ntot, zero);
        for (std::size_t i = 0; i < n1; ++i)
            for (std::size_t j = 0; j < n2; ++j)
                out.alpha[i * n2 + j] = T(a.alpha[i] * b.alpha[j]);
        return out;
    }

    /**
     * Probabilistic mixture: a block-diagonal generator whose initial vector
     * picks branch i with probability PROBS[i]. `aph_simplify` pattern 3,
     * generalised to any number of branches.
     */
    static PhLaw<T> compose_mixture(const std::vector<PhLaw<T>>& laws,
                                    const std::vector<T>& probs) {
        const T zero = num_traits<T>::from_int(0);
        std::size_t total = 0;
        for (const PhLaw<T>& l : laws) total += l.S.rows();

        PhLaw<T> out;
        out.S = Matrix<T>(total, total, zero);
        out.alpha.assign(total, zero);

        std::size_t off = 0;
        for (std::size_t b = 0; b < laws.size(); ++b) {
            const std::size_t nb = laws[b].S.rows();
            for (std::size_t i = 0; i < nb; ++i) {
                for (std::size_t j = 0; j < nb; ++j) out.S(off + i, off + j) = laws[b].S(i, j);
                out.alpha[off + i] = T(probs[b] * laws[b].alpha[i]);
            }
            off += nb;
        }
        return out;
    }

    /**
     * Geometric repetition of a phase-type law, the POST_LOOP semantics.
     *
     * For COUNT >= 1 the body runs at least once and repeats on absorption
     * with probability P = 1-1/COUNT, so
     *
     *   S_out = S + P/D (-S e) alpha,   alpha_out = alpha / D
     *
     * with D = 1 - P(1 - alpha e) the correction for an atom at zero in alpha.
     * The ORDER is that of the body, unlike the COUNT-fold convolution, and the
     * mean is COUNT times the body mean in both cases.
     *
     * For COUNT < 1 the body runs at most once, with probability COUNT; the
     * skipped branch is an immediate phase.
     */
    static PhLaw<T> compose_loop_geometric(const PhLaw<T>& body, const T& count) {
        const T zero = num_traits<T>::from_int(0);
        const T one = num_traits<T>::from_int(1);
        const std::size_t n = body.S.rows();
        const double c = num_traits<T>::to_double(count);

        PhLaw<T> out;
        if (!(c > 0.0)) {
            out.alpha.assign(1, one);
            out.S = Matrix<T>(1, 1, zero);
            out.S(0, 0) = T(-num_traits<T>::from_double(GlobalConstants::Immediate));
            return out;  // zero-time branch
        }

        if (std::abs(c - 1.0) <= GlobalConstants::FineTol) return body;

        if (c < 1.0) {
            // Executed with probability COUNT, skipped otherwise
            out.alpha.assign(n + 1, zero);
            for (std::size_t i = 0; i < n; ++i) out.alpha[i] = T(count * body.alpha[i]);
            out.alpha[n] = T(one - count);
            out.S = Matrix<T>(n + 1, n + 1, zero);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j) out.S(i, j) = body.S(i, j);
            out.S(n, n) = T(-num_traits<T>::from_double(GlobalConstants::Immediate));
            return out;  // zero-time skip branch
        }

        const T p = T(one - one / count);
        T defect = one;
        for (const T& v : body.alpha) defect -= v;
        const T denom = T(one - p * defect);

        out.alpha.assign(n, zero);
        for (std::size_t i = 0; i < n; ++i) out.alpha[i] = T(body.alpha[i] / denom);

        out.S = body.S;
        for (std::size_t i = 0; i < n; ++i) {
            T rate = zero;
            for (std::size_t j = 0; j < n; ++j) rate += body.S(i, j);
            rate = T(-rate * p / denom);
            for (std::size_t j = 0; j < n; ++j) out.S(i, j) += T(rate * body.alpha[j]);
        }
        return out;
    }

    /**
     * COUNT-fold convolution of a phase-type law.
     *
     * The DETERMINISTIC repetition, kept for a caller that genuinely wants an
     * exact number of executions. POST_LOOP is geometric and uses
     * `compose_loop_geometric` instead.
     */
    static PhLaw<T> compose_repeat(const PhLaw<T>& body, long count) {
        if (count <= 0) {
            PhLaw<T> out;
            out.alpha.assign(1, num_traits<T>::from_int(1));
            out.S = Matrix<T>(1, 1, num_traits<T>::from_int(0));
            out.S(0, 0) = T(-num_traits<T>::from_double(GlobalConstants::Immediate));
            return out;
        }
        PhLaw<T> out = body;
        for (long i = 1; i < count; ++i) out = compose_serial(out, body);
        return out;
    }

    /**
     * True when the phase graph of S has no cycle.
     *
     * A geometric loop over a body of two or more phases closes a cycle, so the
     * composed law is a PH and not an APH.
     */
    static bool is_acyclic_generator(const Matrix<T>& S) {
        const std::size_t n = S.rows();
        std::vector<std::vector<bool>> A(n, std::vector<bool>(n, false));
        std::vector<std::size_t> in_deg(n, 0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                if (i == j) continue;
                if (std::abs(num_traits<T>::to_double(S(i, j))) > GlobalConstants::ArcTol) {
                    A[i][j] = true;
                    ++in_deg[j];
                }
            }

        std::vector<std::size_t> queue;
        for (std::size_t i = 0; i < n; ++i)
            if (in_deg[i] == 0) queue.push_back(i);
        std::size_t visited = 0, head = 0;
        while (head < queue.size()) {
            const std::size_t cur = queue[head++];
            ++visited;
            for (std::size_t j = 0; j < n; ++j) {
                if (!A[cur][j]) continue;
                if (--in_deg[j] == 0) queue.push_back(j);
            }
        }
        return visited == n;
    }

    // -----------------------------------------------------------------
    // Precedence factories, named as in MATLAB and the JAR
    // -----------------------------------------------------------------

    static Precedence<T> Serial(const std::string& pre, const std::string& post) {
        Precedence<T> p;
        p.pre_acts.push_back(pre);
        p.post_acts.push_back(post);
        p.pre_type = PrecedenceType::PRE_SEQ;
        p.post_type = PrecedenceType::POST_SEQ;
        return p;
    }

    static std::vector<Precedence<T>> SerialSequence(const std::vector<std::string>& acts) {
        std::vector<Precedence<T>> out;
        for (std::size_t i = 0; i + 1 < acts.size(); ++i)
            out.push_back(Serial(acts[i], acts[i + 1]));
        return out;
    }

    static Precedence<T> AndFork(const std::string& pre,
                                 const std::vector<std::string>& posts) {
        Precedence<T> p;
        p.pre_acts.push_back(pre);
        p.post_acts = posts;
        p.pre_type = PrecedenceType::PRE_SEQ;
        p.post_type = PrecedenceType::POST_AND;
        return p;
    }

    /** An empty QUORUM means a full join; a partial one is refused by validate. */
    static Precedence<T> AndJoin(const std::vector<std::string>& pres,
                                 const std::string& post,
                                 const std::vector<T>& quorum = std::vector<T>()) {
        Precedence<T> p;
        p.pre_acts = pres;
        p.post_acts.push_back(post);
        p.pre_type = PrecedenceType::PRE_AND;
        p.post_type = PrecedenceType::POST_SEQ;
        p.pre_params = quorum;
        return p;
    }

    static Precedence<T> OrFork(const std::string& pre, const std::vector<std::string>& posts,
                                const std::vector<T>& probs) {
        Precedence<T> p;
        p.pre_acts.push_back(pre);
        p.post_acts = posts;
        p.pre_type = PrecedenceType::PRE_SEQ;
        p.post_type = PrecedenceType::POST_OR;
        p.post_params = probs;
        return p;
    }

    static Precedence<T> OrJoin(const std::vector<std::string>& pres,
                                const std::string& post) {
        Precedence<T> p;
        p.pre_acts = pres;
        p.post_acts.push_back(post);
        p.pre_type = PrecedenceType::PRE_OR;
        p.post_type = PrecedenceType::POST_SEQ;
        return p;
    }

    /**
     * A loop: PRE runs once, then POSTS[0..n-2] repeat a geometric number of
     * times of mean COUNT and POSTS[n-1] continues after the loop.
     */
    static Precedence<T> Loop(const std::string& pre, const std::vector<std::string>& posts,
                              const T& count) {
        Precedence<T> p;
        p.pre_acts.push_back(pre);
        p.post_acts = posts;
        p.pre_type = PrecedenceType::PRE_SEQ;
        p.post_type = PrecedenceType::POST_LOOP;
        p.post_params.push_back(count);
        return p;
    }

private:
    // -----------------------------------------------------------------
    // Series-parallel decomposition
    // -----------------------------------------------------------------

    /** How a parsed sequence terminated. */
    enum class ParseStatus { END, STOP, JOIN, FAIL };

    struct ParseState {
        /** Index of the precedence each activity heads / is reached by. */
        std::vector<std::size_t> out_p, in_p;
        std::vector<bool> consumed;
    };

    bool compose_series_parallel(PhLaw<T>& out) {
        if (!sp_built_) {
            if (sp_failed_) return false;
            build_sp_tree();
            if (!sp_built_) return false;
        }
        out = compose_node(sp_tree_.root);
        return true;
    }

    /**
     * Decompose the precedence graph into a series-parallel tree.
     *
     * The decomposition is attempted once per topology: on failure `sp_failed_`
     * is set and the block path takes over.
     */
    bool build_sp_tree() {
        sp_tree_.nodes.clear();
        sp_tree_.root = npos;
        sp_tree_.leaf_of.clear();
        sp_tree_.execs.clear();
        sp_built_ = false;
        sp_failed_ = true;

        const std::size_t n = activities_.size();
        if (n == 0) return false;

        ParseState S;
        S.out_p.assign(n, npos);
        S.in_p.assign(n, npos);
        S.consumed.assign(n, false);

        // An activity may head at most one precedence and be reached by at most
        // one precedence; otherwise the graph is not series-parallel
        for (std::size_t p = 0; p < precedences_.size(); ++p) {
            for (const std::string& nm : precedences_[p].pre_acts) {
                const std::size_t i = activity_index(nm);
                if (i == npos || S.out_p[i] != npos) return false;
                S.out_p[i] = p;
            }
            for (const std::string& nm : precedences_[p].post_acts) {
                const std::size_t j = activity_index(nm);
                if (j == npos || S.in_p[j] != npos) return false;
                S.in_p[j] = p;
            }
        }

        std::vector<std::size_t> starts;
        for (std::size_t i = 0; i < n; ++i)
            if (S.in_p[i] == npos) starts.push_back(i);
        if (starts.size() != 1) return false;

        std::vector<std::size_t> kids;
        std::size_t stop_at = npos;
        const ParseStatus st = sp_parse_seq(S, starts[0], std::vector<std::size_t>(), kids, stop_at);
        if (st != ParseStatus::END) return false;
        for (std::size_t i = 0; i < n; ++i)
            if (!S.consumed[i]) return false;

        const std::size_t root = sp_serial_node(kids);
        if (root == npos) return false;

        sp_tree_.root = root;
        sp_tree_.leaf_of.assign(n, npos);
        for (std::size_t k = 0; k < sp_tree_.nodes.size(); ++k)
            if (sp_tree_.nodes[k].type == SPNodeType::LEAF)
                sp_tree_.leaf_of[sp_tree_.nodes[k].act] = k;
        sp_tree_.execs = sp_execution_counts(sp_tree_);

        sp_built_ = true;
        sp_failed_ = false;
        return true;
    }

    /** Parse a maximal sequence of blocks starting at CUR. */
    ParseStatus sp_parse_seq(ParseState& S, std::size_t cur,
                             const std::vector<std::size_t>& stop_set,
                             std::vector<std::size_t>& kids, std::size_t& stop_at) {
        stop_at = npos;
        while (true) {
            if (cur == npos) return ParseStatus::END;
            if (std::find(stop_set.begin(), stop_set.end(), cur) != stop_set.end()) {
                stop_at = cur;
                return ParseStatus::STOP;
            }
            if (S.consumed[cur]) return ParseStatus::FAIL;
            S.consumed[cur] = true;
            kids.push_back(sp_add_node(SPNodeType::LEAF, cur, std::vector<std::size_t>(),
                                       std::vector<T>(), num_traits<T>::from_int(0)));

            const std::size_t p = S.out_p[cur];
            if (p == npos) return ParseStatus::END;
            const Precedence<T>& prec = precedences_[p];
            if (prec.pre_acts.size() > 1) {
                // CUR is the tail of a branch: the caller composes the join
                stop_at = p;
                return ParseStatus::JOIN;
            }

            std::vector<std::size_t> post_inds = sp_indices_of(prec.post_acts);
            for (std::size_t i : post_inds)
                if (i == npos) return ParseStatus::FAIL;

            std::size_t knode = npos, next_act = npos;
            switch (prec.post_type) {
                case PrecedenceType::POST_AND:
                    if (!sp_parse_fork(S, post_inds, std::vector<T>(), stop_set, true, knode,
                                       next_act))
                        return ParseStatus::FAIL;
                    kids.push_back(knode);
                    cur = next_act;
                    break;
                case PrecedenceType::POST_OR:
                    if (prec.post_params.size() != post_inds.size()) return ParseStatus::FAIL;
                    if (!sp_parse_fork(S, post_inds, prec.post_params, stop_set, false, knode,
                                       next_act))
                        return ParseStatus::FAIL;
                    kids.push_back(knode);
                    cur = next_act;
                    break;
                case PrecedenceType::POST_LOOP:
                    if (!sp_parse_loop(S, post_inds, prec.post_params, stop_set, knode, next_act))
                        return ParseStatus::FAIL;
                    kids.push_back(knode);
                    cur = next_act;
                    break;
                case PrecedenceType::POST_SEQ:
                    if (post_inds.size() != 1) return ParseStatus::FAIL;
                    cur = post_inds[0];
                    break;
                default:
                    // POST_CACHE and any other pattern is not a workflow
                    // composition rule
                    return ParseStatus::FAIL;
            }
        }
    }

    /** Parse the branches of a fork and their join. */
    bool sp_parse_fork(ParseState& S, const std::vector<std::size_t>& branch_heads,
                       const std::vector<T>& probs, const std::vector<std::size_t>& stop_set,
                       bool is_and, std::size_t& knode, std::size_t& next_act) {
        const std::size_t nb = branch_heads.size();
        std::vector<std::size_t> branch_nodes(nb, npos), bstop(nb, npos);
        std::vector<ParseStatus> bstatus(nb, ParseStatus::FAIL);

        for (std::size_t b = 0; b < nb; ++b) {
            std::vector<std::size_t> bkids;
            std::size_t sa = npos;
            const ParseStatus st = sp_parse_seq(S, branch_heads[b], stop_set, bkids, sa);
            if (st == ParseStatus::FAIL) return false;
            const std::size_t bn = sp_serial_node(bkids);
            if (bn == npos) return false;
            branch_nodes[b] = bn;
            bstatus[b] = st;
            bstop[b] = sa;
        }

        const bool all_join =
            std::all_of(bstatus.begin(), bstatus.end(),
                        [](ParseStatus s) { return s == ParseStatus::JOIN; });
        const bool all_end = std::all_of(bstatus.begin(), bstatus.end(),
                                         [](ParseStatus s) { return s == ParseStatus::END; });
        const bool all_stop = std::all_of(bstatus.begin(), bstatus.end(),
                                          [](ParseStatus s) { return s == ParseStatus::STOP; });
        const bool same_stop =
            std::all_of(bstop.begin(), bstop.end(),
                        [&bstop](std::size_t x) { return x == bstop[0]; });

        if (all_join) {
            if (!same_stop) return false;
            const Precedence<T>& join_prec = precedences_[bstop[0]];
            if (join_prec.pre_acts.size() != nb) return false;
            if (is_and) {
                if (join_prec.pre_type != PrecedenceType::PRE_AND) return false;
            } else {
                if (join_prec.pre_type != PrecedenceType::PRE_OR) return false;
            }
            const std::vector<std::size_t> post_inds = sp_indices_of(join_prec.post_acts);
            if (post_inds.size() != 1 || post_inds[0] == npos) return false;
            next_act = post_inds[0];
        } else if (all_end) {
            // Branches terminate the workflow. An AND-fork with no join still
            // synchronises at the end of the workflow
            next_act = npos;
        } else if (!is_and && all_stop && same_stop) {
            next_act = bstop[0];
        } else {
            return false;
        }

        knode = sp_add_node(is_and ? SPNodeType::PAR : SPNodeType::OR, npos, branch_nodes, probs,
                            num_traits<T>::from_int(0));
        return true;
    }

    /**
     * Parse a loop block: the last post activity continues after the loop, the
     * others form the body.
     */
    bool sp_parse_loop(ParseState& S, const std::vector<std::size_t>& post_inds,
                       const std::vector<T>& counts, const std::vector<std::size_t>& stop_set,
                       std::size_t& knode, std::size_t& next_act) {
        if (counts.size() != 1) return false;
        const T count = counts[0];

        std::vector<std::size_t> body_acts;
        std::size_t end_act = npos;
        if (post_inds.size() >= 2) {
            body_acts.assign(post_inds.begin(), post_inds.end() - 1);
            end_act = post_inds.back();
        } else {
            body_acts.push_back(post_inds[0]);
        }

        std::vector<std::size_t> loop_stop = stop_set;
        loop_stop.insert(loop_stop.end(), body_acts.begin(), body_acts.end());
        if (end_act != npos) loop_stop.push_back(end_act);

        std::vector<std::size_t> body_kids;
        std::size_t j = 0;
        while (j < body_acts.size()) {
            const std::size_t a = body_acts[j];
            if (S.consumed[a]) {
                ++j;
                continue;
            }
            std::vector<std::size_t> this_stop;
            for (std::size_t x : loop_stop)
                if (x != a) this_stop.push_back(x);

            std::vector<std::size_t> kk;
            std::size_t sa = npos;
            const ParseStatus st = sp_parse_seq(S, a, this_stop, kk, sa);
            if (st == ParseStatus::FAIL) return false;
            body_kids.insert(body_kids.end(), kk.begin(), kk.end());

            if (st == ParseStatus::END) {
                ++j;
            } else if (st == ParseStatus::STOP) {
                const auto it = std::find(body_acts.begin(), body_acts.end(), sa);
                if (it != body_acts.end()) {
                    j = static_cast<std::size_t>(it - body_acts.begin());
                } else if (end_act != npos && sa == end_act) {
                    j = body_acts.size();
                } else {
                    return false;
                }
            } else {
                // A join reached from inside the body crosses the loop
                // boundary, so the graph is not series-parallel
                return false;
            }
        }

        const std::size_t body_node = sp_serial_node(body_kids);
        if (body_node == npos) return false;

        knode = sp_add_node(SPNodeType::LOOP, npos, std::vector<std::size_t>(1, body_node),
                            std::vector<T>(), count);
        next_act = end_act;
        return true;
    }

    /** Wrap a list of nodes in a serial node; a single node is returned as is. */
    std::size_t sp_serial_node(const std::vector<std::size_t>& kids) {
        if (kids.empty()) return npos;
        if (kids.size() == 1) return kids[0];
        return sp_add_node(SPNodeType::SERIAL, npos, kids, std::vector<T>(),
                           num_traits<T>::from_int(0));
    }

    std::size_t sp_add_node(SPNodeType type, std::size_t act,
                            const std::vector<std::size_t>& kids, const std::vector<T>& probs,
                            const T& count) {
        SPNode<T> node;
        node.type = type;
        node.act = act;
        node.kids = kids;
        node.probs = probs;
        node.count = count;
        sp_tree_.nodes.push_back(node);
        const std::size_t k = sp_tree_.nodes.size() - 1;
        for (std::size_t c : kids) sp_tree_.nodes[c].parent = k;
        return k;
    }

    std::vector<std::size_t> sp_indices_of(const std::vector<std::string>& names) const {
        std::vector<std::size_t> out(names.size(), npos);
        for (std::size_t i = 0; i < names.size(); ++i) out[i] = activity_index(names[i]);
        return out;
    }

    /**
     * Compose one node. A cached node is returned untouched, so a demand change
     * only recomposes the path from the dirty leaf to the root.
     */
    PhLaw<T> compose_node(std::size_t k) {
        SPNode<T>& node = sp_tree_.nodes[k];
        if (node.valid) return node.law;

        PhLaw<T> law;
        switch (node.type) {
            case SPNodeType::LEAF:
                law = activities_[node.act].ph_representation();
                break;
            case SPNodeType::SERIAL: {
                const std::vector<std::size_t> kids = node.kids;
                law = compose_node(kids[0]);
                for (std::size_t i = 1; i < kids.size(); ++i)
                    law = compose_serial(law, compose_node(kids[i]));
                break;
            }
            case SPNodeType::PAR: {
                const std::vector<std::size_t> kids = node.kids;
                law = compose_node(kids[0]);
                for (std::size_t i = 1; i < kids.size(); ++i)
                    law = compose_parallel(law, compose_node(kids[i]));
                break;
            }
            case SPNodeType::OR: {
                const std::vector<std::size_t> kids = node.kids;
                std::vector<PhLaw<T>> laws;
                for (std::size_t c : kids) laws.push_back(compose_node(c));
                law = compose_mixture(laws, sp_tree_.nodes[k].probs);
                break;
            }
            case SPNodeType::LOOP: {
                const std::size_t child = node.kids[0];
                law = compose_loop_geometric(compose_node(child), sp_tree_.nodes[k].count);
                break;
            }
        }

        sp_tree_.nodes[k].law = law;
        sp_tree_.nodes[k].valid = true;
        return law;
    }

    static void invalidate_branch(SPTree<T>& tree, std::size_t k) {
        while (k != npos) {
            tree.nodes[k].valid = false;
            k = tree.nodes[k].parent;
        }
    }

    /**
     * Expected executions of each node per workflow run.
     *
     * Serial and parallel children inherit the count of their parent, an OR
     * branch is weighted by its probability, and a loop body by the loop count.
     */
    static std::vector<T> sp_execution_counts(const SPTree<T>& tree) {
        std::vector<T> execs(tree.nodes.size(), num_traits<T>::from_int(0));
        if (tree.root == npos) return execs;
        execs[tree.root] = num_traits<T>::from_int(1);
        std::vector<std::size_t> stack(1, tree.root);
        while (!stack.empty()) {
            const std::size_t k = stack.back();
            stack.pop_back();
            const SPNode<T>& node = tree.nodes[k];
            for (std::size_t i = 0; i < node.kids.size(); ++i) {
                const std::size_t c = node.kids[i];
                if (node.type == SPNodeType::OR) {
                    execs[c] = T(execs[k] * node.probs[i]);
                } else if (node.type == SPNodeType::LOOP) {
                    execs[c] = T(execs[k] * node.count);
                } else {
                    execs[c] = execs[k];
                }
                stack.push_back(c);
            }
        }
        return execs;
    }

    // -----------------------------------------------------------------
    // Block composition, the fallback for a graph that is not series-parallel
    // -----------------------------------------------------------------

    struct ForkInfo {
        bool is_and = true;
        std::size_t pre_act = npos;
        std::vector<std::size_t> post_acts;
        std::vector<T> probs;
    };

    struct JoinInfo {
        bool is_and = true;
        std::vector<std::size_t> pre_acts;
        std::size_t post_act = npos;
    };

    struct LoopInfo {
        std::size_t pre_act = npos;
        std::vector<std::size_t> body_acts;
        std::size_t end_act = npos;
        T count = num_traits<T>::from_int(1);
    };

    struct Structure {
        std::vector<std::vector<std::size_t>> adj;
        std::vector<std::size_t> in_deg, out_deg;
        std::vector<ForkInfo> forks;
        std::vector<JoinInfo> joins;
        std::vector<LoopInfo> loops;
    };

    Structure analyze_structure() const {
        const std::size_t n = activities_.size();
        Structure st;
        st.adj.assign(n, std::vector<std::size_t>());
        st.in_deg.assign(n, 0);
        st.out_deg.assign(n, 0);

        for (const Precedence<T>& prec : precedences_) {
            const std::vector<std::size_t> pre = sp_indices_of(prec.pre_acts);
            const std::vector<std::size_t> post = sp_indices_of(prec.post_acts);
            for (std::size_t i : pre)
                for (std::size_t j : post) {
                    st.adj[i].push_back(j);
                    ++st.out_deg[i];
                    ++st.in_deg[j];
                }

            if (prec.post_type == PrecedenceType::POST_AND) {
                ForkInfo f;
                f.is_and = true;
                f.pre_act = pre[0];
                f.post_acts = post;
                st.forks.push_back(f);
            } else if (prec.post_type == PrecedenceType::POST_OR) {
                ForkInfo f;
                f.is_and = false;
                f.pre_act = pre[0];
                f.post_acts = post;
                f.probs = prec.post_params;
                st.forks.push_back(f);
            } else if (prec.post_type == PrecedenceType::POST_LOOP) {
                LoopInfo l;
                l.pre_act = pre[0];
                if (post.size() > 1) {
                    l.body_acts.assign(post.begin(), post.end() - 1);
                    l.end_act = post.back();
                } else {
                    l.body_acts = post;
                }
                l.count = prec.post_params.empty() ? num_traits<T>::from_int(1)
                                                   : prec.post_params[0];
                st.loops.push_back(l);
            }

            if (prec.pre_type == PrecedenceType::PRE_AND ||
                prec.pre_type == PrecedenceType::PRE_OR) {
                JoinInfo j;
                j.is_and = prec.pre_type == PrecedenceType::PRE_AND;
                j.pre_acts = pre;
                j.post_act = post[0];
                st.joins.push_back(j);
            }
        }
        return st;
    }

    std::vector<std::size_t> topological_sort(
        const std::vector<std::vector<std::size_t>>& adj) const {
        const std::size_t n = activities_.size();
        std::vector<std::size_t> in_deg(n, 0);
        for (const std::vector<std::size_t>& nbrs : adj)
            for (std::size_t j : nbrs) ++in_deg[j];

        std::vector<std::size_t> queue, order;
        for (std::size_t i = 0; i < n; ++i)
            if (in_deg[i] == 0) queue.push_back(i);
        std::size_t head = 0;
        while (head < queue.size()) {
            const std::size_t cur = queue[head++];
            order.push_back(cur);
            for (std::size_t nx : adj[cur])
                if (--in_deg[nx] == 0) queue.push_back(nx);
        }
        std::vector<bool> seen(n, false);
        for (std::size_t i : order) seen[i] = true;
        for (std::size_t i = 0; i < n; ++i)
            if (!seen[i]) order.push_back(i);
        return order;
    }

    /**
     * The block composition.
     *
     * A HEURISTIC, and the reason the series-parallel reduction exists: it
     * folds each fork/join and loop block independently and then chains what is
     * left in topological order, so a block nested inside another is not
     * reduced exactly.
     */
    PhLaw<T> build_ctmc() {
        const std::size_t n = activities_.size();
        if (n == 1) return activities_[0].ph_representation();

        const Structure st = analyze_structure();

        std::vector<PhLaw<T>> block(n);
        std::vector<bool> absorbed(n, false);
        for (std::size_t i = 0; i < n; ++i) block[i] = activities_[i].ph_representation();

        for (const LoopInfo& loop : st.loops) {
            PhLaw<T> body = activities_[loop.body_acts[0]].ph_representation();
            for (std::size_t j = 1; j < loop.body_acts.size(); ++j)
                body = compose_serial(body, activities_[loop.body_acts[j]].ph_representation());

            PhLaw<T> res = compose_serial(block[loop.pre_act],
                                          compose_loop_geometric(body, loop.count));
            if (loop.end_act != npos) {
                res = compose_serial(res, activities_[loop.end_act].ph_representation());
                absorbed[loop.end_act] = true;
            }
            block[loop.pre_act] = res;
            for (std::size_t idx : loop.body_acts) absorbed[idx] = true;
        }

        for (const ForkInfo& fork : st.forks) {
            const JoinInfo* join = find_matching_join(fork.post_acts, st.joins, fork.is_and);
            if (fork.is_and && join == nullptr) continue;

            PhLaw<T> inner;
            if (fork.is_and) {
                inner = block[fork.post_acts[0]];
                for (std::size_t i = 1; i < fork.post_acts.size(); ++i)
                    inner = compose_parallel(inner, block[fork.post_acts[i]]);
            } else {
                std::vector<PhLaw<T>> laws;
                for (std::size_t idx : fork.post_acts) laws.push_back(block[idx]);
                inner = compose_mixture(laws, fork.probs);
            }

            PhLaw<T> res = compose_serial(block[fork.pre_act], inner);
            if (join != nullptr && !absorbed[join->post_act]) {
                res = compose_serial(res, block[join->post_act]);
                absorbed[join->post_act] = true;
            }
            block[fork.pre_act] = res;
            for (std::size_t idx : fork.post_acts) absorbed[idx] = true;
        }

        const std::vector<std::size_t> order = topological_sort(st.adj);
        bool started = false;
        PhLaw<T> out;
        for (std::size_t idx : order) {
            if (absorbed[idx]) continue;
            if (!started) {
                out = block[idx];
                started = true;
            } else {
                out = compose_serial(out, block[idx]);
            }
        }
        if (!started) out = activities_[0].ph_representation();
        return out;
    }

    static const JoinInfo* find_matching_join(const std::vector<std::size_t>& post_acts,
                                              const std::vector<JoinInfo>& joins, bool is_and) {
        std::vector<std::size_t> want = post_acts;
        std::sort(want.begin(), want.end());
        for (const JoinInfo& j : joins) {
            if (j.is_and != is_and) continue;
            std::vector<std::size_t> have = j.pre_acts;
            std::sort(have.begin(), have.end());
            if (have == want) return &j;
        }
        return nullptr;
    }

    std::string name_;
    std::vector<WorkflowActivity<T>> activities_;
    std::map<std::string, std::size_t> activity_map_;
    std::vector<Precedence<T>> precedences_;

    Distrib<T> cached_ph_;
    bool cached_valid_ = false;
    SPTree<T> sp_tree_;
    bool sp_built_ = false;
    bool sp_failed_ = false;
};

}  // namespace workflow
}  // namespace line

#endif  // LINE_LANG_WORKFLOW_WORKFLOW_H
