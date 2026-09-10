/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_QN_ENVIRONMENT_H
#define LINE_LANG_QN_ENVIRONMENT_H

/**
 * A random environment: a port of `matlab/src/lang/Environment.m`, restricted
 * to what `SolverENV` reads out of it.
 *
 * WHAT IT MODELS. A queueing network whose PARAMETERS change over time because
 * the world around it does: a server that breaks down and is repaired, a
 * workload that has a day phase and a night phase, a system that degrades and
 * is reset. Each of those is a STAGE, holding its own complete network, and the
 * environment is a semi-Markov process over the stages. The point is that the
 * network never restarts empty at a switch -- the jobs in it at the moment of
 * the switch are carried into the next stage, which is what couples the stages
 * and what makes this more than solving each stage separately.
 *
 * THE HOLDING TIME IS A COMPETING RISK, and this is the part worth reading the
 * reference for. Every enabled transition e -> h has its own distribution, and
 * they all run at once: the stage ends when the FIRST of them fires, and which
 * one fired decides the next stage. `init()` therefore superposes the outgoing
 * transitions of e by a Kronecker sum into a single marked process
 * `hold_time[e]`, whose per-mark rates give the embedded jump chain `Pemb`.
 * The collapse that follows -- each block's row sums moved into its FIRST
 * column -- is the reference's, and it makes the superposed process restart in
 * phase one after every jump, i.e. renewal at each stage entry.
 *
 * WHAT init() PRODUCES, and what the solver consumes:
 *   proc[e][h]   the e -> h transition process, whose CDF weights the stage
 *                transient when computing the metrics AT a switch to h
 *   hold_time[e] the superposed holding time of stage e, whose CDF weights the
 *                transient when computing the metrics OVER the whole stage
 *   prob_env[e]  the stationary probability of being in stage e
 *   prob_orig    prob_orig(h, e) = P(the previous stage was h | now entering e)
 *
 * A DISABLED transition is the 1 x 1 zero pair, exactly as in the reference:
 * krons(A, 0) leaves A unchanged, so a disabled arc costs nothing and needs no
 * special case anywhere below.
 *
 * NODE BREAKDOWNS are the one stage pattern the reference gives a name to
 * (`addNodeBreakdown` / `addNodeRepair` / `addNodeFailureRepair`): the UP stage
 * holds the base network, the DOWN_<node> stage holds the same network with one
 * node's service replaced by its degraded distribution, and the two arcs
 * between them are the time to failure and the time to repair. It is a macro
 * over `set_stage` and `add_transition` and nothing more, EXCEPT for the pair of
 * queue-length reset policies it attaches, which are the only part of a
 * breakdown that the expanded stages cannot express -- hence `NodeFailure`,
 * which records them beside the stages so an environment read back from
 * `model.json` is the same model it was written from.
 */

#include <cctype>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mc/ctmc_courtois.h"
#include "line/lang/distribution.h"
#include "line/lang/lang_types.h"
#include "line/lang/lqn/lqn_struct.h"
#include "line/lang/prior.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace env {

/**
 * The reset policy of a transition, `resetFun` in the reference.
 *
 * It maps the mean queue lengths at the moment of the switch onto the mean
 * queue lengths the next stage starts from. Identity means the jobs are simply
 * carried over; zero means the buffer is flushed on the switch.
 */
using ResetMarginal = std::function<Matrix<double>(const Matrix<double>&)>;

/**
 * The two NAMED reset policies of `Environment.resolveResetPolicy`, which are
 * the only ones the JSON interchange can carry (a function handle is written as
 * `custom` and warned about by both writers, never reloaded).
 *
 * `keep` resolves to the EMPTY function rather than to an explicit identity,
 * because empty is how this port spells identity everywhere a reset is read:
 * `SolverEnv::post` skips the call, and the compression's macro-arc fold
 * compares resets only by whether one is PRESENT, so an explicit identity on
 * one arc and nothing on another would be refused as a disagreement although
 * the two mean the same thing.
 */
inline ResetMarginal env_reset_policy(const std::string& name) {
    std::string low;
    low.reserve(name.size());
    for (std::size_t i = 0; i < name.size(); ++i)
        low.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(name[i]))));
    if (low == "keep") return ResetMarginal();
    if (low == "clear")
        return [](const Matrix<double>& q) { return Matrix<double>(q.rows(), q.cols(), 0.0); };
    throw InputError("Environment: unknown reset policy '" + name +
                     "'. The serializable policies are 'keep' (carry the queue lengths across "
                     "the switch) and 'clear' (empty the queues); a custom policy is a function "
                     "and is installed through add_transition, not by name");
}

/**
 * `Environment.nodeFailures{k}`: the declarative record of one node breakdown.
 *
 * The expanded stages and transitions already carry the STRUCTURE of a
 * breakdown; what only this record carries is the pair of reset policies, which
 * decide what the next stage starts from and therefore change the numbers. It
 * is kept on the environment for the same reason the reference keeps it: so an
 * environment read from `model.json` in its expanded form can be given back its
 * policies, and so it serializes out again as the same model.
 */
template <class T>
struct NodeFailure {
    std::string node;            ///< the node that breaks down
    lang::Distrib<T> breakdown;  ///< time to failure, the UP -> DOWN transition
    lang::Distrib<T> down_service;  ///< the node's service while it is down
    lang::Distrib<T> repair;     ///< time to repair, the DOWN -> UP transition
    bool has_repair = false;     ///< false for a breakdown with no repair arc
    std::string breakdown_reset = "keep";
    std::string repair_reset;    ///< empty when there is no repair
};

/**
 * One stage: a name, a category, and the model in force while it lasts.
 *
 * THE MODEL IS ONE OF TWO KINDS, and which one is the stage's own property
 * rather than the environment's: a flat `qn::NetworkStruct`, or a
 * `lqn::LqnStruct` -- a LayeredNetwork, whose stations and classes are those of
 * the LAYERS SolverLN builds out of it. The reference draws no distinction at
 * this level either (`Environment.addStage` stores whatever model it is handed
 * and `SolverENV` branches on its class), and the two fields are kept apart
 * rather than unified because nothing about them is shared below the name: a
 * flat stage is integrated directly, a layered one through a whole fixed point
 * over its layers.
 *
 * `has_model` and `has_lqn` are mutually exclusive; `set_stage` and
 * `set_lqn_stage` each clear the other, so a stage re-declared as the other
 * kind cannot leave a stale twin behind for a later consumer to read.
 */
template <class T>
struct EnvStage {
    std::string name;
    std::string type;  ///< the stage category, informational only
    qn::NetworkStruct<T> model;
    bool has_model = false;
    /** The layered model, when this stage holds a LayeredNetwork. */
    lqn::LqnStruct<T> lqn_model;
    bool has_lqn = false;
    /** True when the stage carries a model of either kind. */
    bool declared() const { return has_model || has_lqn; }
};

/**
 * `resetEnvRatesFun` in the reference: the state-dependent environment rate.
 *
 * It is given the arc's CURRENT transition distribution and the exit metrics of
 * the stage the arc leaves -- mean queue lengths, utilizations and throughputs,
 * averaged over when that arc fires -- and returns the distribution the arc
 * should carry next. It is what makes the environment process depend on the
 * network it modulates, and `method = "statedep"` is what applies it.
 */
template <class T>
using ResetEnvRates =
    std::function<lang::Distrib<T>(const lang::Distrib<T>&, const Matrix<double>&,
                                   const Matrix<double>&, const Matrix<double>&)>;

/** One arc of the environment process. */
template <class T>
struct EnvArc {
    bool enabled = false;
    lang::Distrib<T> dist;   ///< the e -> h transition time
    ResetMarginal reset;     ///< empty means the identity
    ResetEnvRates<T> reset_rates;  ///< empty means the rate does not depend on the state
};

/**
 * The DOWN stage network of `addNodeBreakdown`: the base network with ONE
 * node's service replaced by its degraded distribution, for EVERY class.
 *
 * Every class, and not only the ones that were enabled there, is what the
 * reference does (`for c = 1:length(classes), nodes{nodeIdx}.setService(...)`),
 * so a class that was disabled at the node while it was up is served at the
 * degraded rate while it is down. The whole refresh chain is rerun afterwards
 * because `rates`, `scv` and the chain-derived tables are all read off the
 * service table; editing the table alone would leave the struct describing the
 * UP stage and the solver reading the DOWN one.
 */
template <class T>
qn::NetworkStruct<T> env_degraded_model(const qn::NetworkStruct<T>& base,
                                        const std::string& node_name,
                                        const lang::Distrib<T>& down_service) {
    qn::NetworkStruct<T> sn = base;
    std::size_t node = 0;
    for (std::size_t i = 0; i < sn.nodes.size(); ++i)
        if (sn.nodes[i].name == node_name) {
            node = i + 1;
            break;
        }
    if (node == 0)
        throw InputError("Environment: node '" + node_name +
                         "' is not in the base model, so it has no service to degrade");
    const std::size_t ist = sn.nodes[node - 1].station;
    if (ist == 0)
        throw InputError("Environment: node '" + node_name +
                         "' is not a station and carries no service distribution, so it cannot "
                         "break down into a degraded service");
    lang::Distrib<T> d = down_service;
    if (d.is_prior())
        lang::prior_refresh_moments(d);
    else
        lang::dist_refresh_moments(d);
    for (std::size_t r = 1; r <= sn.classes.size(); ++r) sn.set_service(ist, r, d);
    sn.refresh_struct();
    return sn;
}

template <class T>
class Environment {
public:
    Environment(const std::string& nm, std::size_t nstages)
        : name_(nm), stages_(nstages), arcs_(nstages, std::vector<EnvArc<T>>(nstages)) {
        if (nstages == 0) throw InputError("Environment: a random environment needs a stage");
    }

    std::size_t nstages() const { return stages_.size(); }
    const std::string& name() const { return name_; }

    /** `addStage`: name the stage and give it its network. */
    void set_stage(std::size_t e, const std::string& nm, const std::string& type,
                   const qn::NetworkStruct<T>& model) {
        check(e);
        stages_[e].name = nm;
        stages_[e].type = type;
        stages_[e].model = model;
        stages_[e].has_model = true;
        stages_[e].has_lqn = false;
        stages_[e].lqn_model = lqn::LqnStruct<T>();
    }

    /**
     * `addStage` with a LayeredNetwork: the stage holds a LAYERED model.
     *
     * The environment itself does nothing with the difference -- probEnv,
     * probOrig and the holding times are read off the ARCS and never off a
     * stage model -- so the whole content of this overload is that the stage
     * records which kind of model it carries and `SolverEnv` runs the matching
     * stage solver. What the two kinds must still agree on is the (station,
     * class) SHAPE of the metrics being blended, and for a layered stage that
     * shape is the block-diagonal union of its layers; `SolverEnv::init` is
     * where the shapes are compared, because only there is a SolverLN built and
     * the layer blocks known.
     */
    void set_lqn_stage(std::size_t e, const std::string& nm, const std::string& type,
                       const lqn::LqnStruct<T>& model) {
        check(e);
        stages_[e].name = nm;
        stages_[e].type = type;
        stages_[e].lqn_model = model;
        stages_[e].has_lqn = true;
        stages_[e].has_model = false;
        stages_[e].model = qn::NetworkStruct<T>();
    }

    /** True when stage `e` holds a LayeredNetwork rather than a flat network. */
    bool is_lqn(std::size_t e) const {
        check(e);
        return stages_[e].has_lqn;
    }

    /** True when ANY stage holds a LayeredNetwork. */
    bool has_lqn_stages() const {
        for (std::size_t e = 0; e < stages_.size(); ++e)
            if (stages_[e].has_lqn) return true;
        return false;
    }

    /**
     * Refuse an environment carrying a LayeredNetwork stage, by name.
     *
     * Shared by every consumer that reads `stage(e).model` directly -- the
     * state-vector coupling, the closed-form limits, the compression -- so that
     * each says the same thing about the same gap rather than reading an EMPTY
     * NetworkStruct and reporting a confident answer about a model that is not
     * there. `who` names the caller and `why` says what about it needs a flat
     * stage.
     */
    void reject_lqn_stages(const std::string& who, const std::string& why) const {
        for (std::size_t e = 0; e < stages_.size(); ++e)
            if (stages_[e].has_lqn)
                throw UnsupportedError(
                    who + ": stage " + std::to_string(e + 1) + " ('" + stages_[e].name +
                    "') holds a LayeredNetwork, and " + why +
                    ". The mean-field coupling (method 'default' / 'meanfield') is the one that "
                    "solves a layered stage, through SolverLN");
    }

    /** `addTransition`: enable e -> h with a distribution and a reset policy. */
    void add_transition(std::size_t e, std::size_t h, const lang::Distrib<T>& d,
                        const ResetMarginal& reset = ResetMarginal()) {
        check(e);
        check(h);
        arcs_[e][h].enabled = true;
        arcs_[e][h].dist = d;
        arcs_[e][h].reset = reset;
    }

    const EnvStage<T>& stage(std::size_t e) const {
        check(e);
        return stages_[e];
    }
    const EnvArc<T>& arc(std::size_t e, std::size_t h) const {
        check(e);
        check(h);
        return arcs_[e][h];
    }

    /**
     * `resetEnvRatesFun{e,h}`: make the e -> h transition depend on the state
     * the stage is left in. Only `method = "statedep"` reads it.
     */
    void set_env_rate_reset(std::size_t e, std::size_t h, const ResetEnvRates<T>& f) {
        check(e);
        check(h);
        if (!arcs_[e][h].enabled)
            throw InputError("Environment: no transition " + stages_[e].name + " -> " +
                             stages_[h].name +
                             " is declared, so its rate cannot be made state dependent");
        arcs_[e][h].reset_rates = f;
    }

    /**
     * Replace the distribution of an arc that is already declared.
     *
     * This is what the state-dependent method writes back each iteration; every
     * other caller declares the arc once through `add_transition`. `init()` must
     * be rerun afterwards, because the superposed holding times and the stage
     * probabilities are all derived from these distributions.
     */
    void set_transition_dist(std::size_t e, std::size_t h, const lang::Distrib<T>& d) {
        check(e);
        check(h);
        if (!arcs_[e][h].enabled)
            throw InputError("Environment: no transition " + stages_[e].name + " -> " +
                             stages_[h].name + " is declared, so its distribution cannot be set");
        arcs_[e][h].dist = d;
    }

    /** The index of the stage called `nm`, or `nstages()` when there is none. */
    std::size_t find_stage(const std::string& nm) const {
        for (std::size_t e = 0; e < stages_.size(); ++e)
            if (stages_[e].has_model && stages_[e].name == nm) return e;
        return stages_.size();
    }

    /**
     * Install a reset policy on an arc that is already declared.
     *
     * The arc must exist: a reset on a disabled arc is a policy for a switch
     * the environment cannot make, and it would sit there reporting nothing.
     * The reference reaches this through `setBreakdownResetPolicy`, which
     * likewise errors when the stages it names are absent.
     */
    void set_reset(std::size_t e, std::size_t h, const ResetMarginal& reset) {
        check(e);
        check(h);
        if (!arcs_[e][h].enabled)
            throw InputError("Environment: no transition " + stages_[e].name + " -> " +
                             stages_[h].name +
                             " is declared, so it cannot be given a reset policy");
        arcs_[e][h].reset = reset;
    }

    // ---- node breakdown and repair ---------------------------------------

    /** `nodeFailures`, the declarative record of the breakdowns declared here. */
    const std::vector<NodeFailure<T>>& node_failures() const { return node_failures_; }

    /** `findNodeFailure`: the descriptor for `nm`, or `node_failures().size()`. */
    std::size_t find_node_failure(const std::string& nm) const {
        for (std::size_t i = 0; i < node_failures_.size(); ++i)
            if (node_failures_[i].node == nm) return i;
        return node_failures_.size();
    }

    /** The name `addNodeBreakdown` gives the stage in which `nm` is down. */
    static std::string down_stage_name(const std::string& nm) { return "DOWN_" + nm; }

    /**
     * Port of `addNodeBreakdown`, on a FIXED stage count.
     *
     * The reference grows its stage graph as breakdowns are declared; this
     * environment is sized at construction, exactly as `set_stage` is, so the
     * caller says which slot is UP and which is the DOWN stage of this node.
     * Everything else is the reference's: the UP stage holds the base model and
     * is named `UP`, the DOWN stage holds the degraded copy and is named
     * `DOWN_<node>`, and the UP -> DOWN arc carries the breakdown time and the
     * breakdown reset policy.
     */
    void add_node_breakdown(std::size_t up, std::size_t down, const qn::NetworkStruct<T>& base,
                            const std::string& node_name, const lang::Distrib<T>& breakdown,
                            const lang::Distrib<T>& down_service,
                            const std::string& reset_policy = "keep") {
        check(up);
        check(down);
        if (up == down)
            throw InputError("Environment: the UP and DOWN stages of node '" + node_name +
                             "' must be different stages");
        if (!stages_[up].has_model) set_stage(up, "UP", "operational", base);
        set_stage(down, down_stage_name(node_name), "failed",
                  env_degraded_model(base, node_name, down_service));
        add_transition(up, down, breakdown, env_reset_policy(reset_policy));

        NodeFailure<T> nf;
        nf.node = node_name;
        nf.breakdown = breakdown;
        nf.down_service = down_service;
        nf.breakdown_reset = reset_policy;
        record_node_failure(nf);
    }

    /** Port of `addNodeRepair`: the DOWN_<node> -> UP arc and its policy. */
    void add_node_repair(const std::string& node_name, const lang::Distrib<T>& repair,
                         const std::string& reset_policy = "keep") {
        const std::size_t down = find_stage(down_stage_name(node_name));
        const std::size_t up = find_stage("UP");
        if (down == stages_.size())
            throw InputError("Environment: stage '" + down_stage_name(node_name) +
                             "' is not defined; declare the breakdown of node '" + node_name +
                             "' before its repair");
        if (up == stages_.size())
            throw InputError("Environment: no UP stage is defined, so node '" + node_name +
                             "' has nothing to be repaired into");
        add_transition(down, up, repair, env_reset_policy(reset_policy));

        const std::size_t idx = find_node_failure(node_name);
        if (idx == node_failures_.size())
            throw InputError("Environment: no breakdown is recorded for node '" + node_name +
                             "', so its repair would describe a failure that was never declared");
        node_failures_[idx].repair = repair;
        node_failures_[idx].has_repair = true;
        node_failures_[idx].repair_reset = reset_policy;
    }

    /** `addNodeFailureRepair`: the two calls above, in order. */
    void add_node_failure_repair(std::size_t up, std::size_t down,
                                 const qn::NetworkStruct<T>& base, const std::string& node_name,
                                 const lang::Distrib<T>& breakdown, const lang::Distrib<T>& repair,
                                 const lang::Distrib<T>& down_service,
                                 const std::string& breakdown_reset = "keep",
                                 const std::string& repair_reset = "keep") {
        add_node_breakdown(up, down, base, node_name, breakdown, down_service, breakdown_reset);
        add_node_repair(node_name, repair, repair_reset);
    }

    /**
     * Port of `registerNodeFailure`: attach a breakdown descriptor, and its
     * reset policies, to stages that ALREADY exist.
     *
     * This is the read path of an environment saved in its expanded form: the
     * `UP` and `DOWN_<node>` stages and their two arcs came off the wire, and
     * the only thing the wire could not carry is the pair of policies, which is
     * what this installs.
     */
    void register_node_failure(const std::string& node_name, const lang::Distrib<T>& breakdown,
                               const lang::Distrib<T>& down_service, bool has_repair,
                               const lang::Distrib<T>& repair,
                               const std::string& breakdown_reset,
                               const std::string& repair_reset) {
        const std::size_t up = find_stage("UP");
        const std::size_t down = find_stage(down_stage_name(node_name));
        if (up == stages_.size())
            throw InputError("Environment: cannot register a node failure on '" + node_name +
                             "': no UP stage is defined in this environment");
        if (down == stages_.size())
            throw InputError("Environment: cannot register a node failure on '" + node_name +
                             "': no '" + down_stage_name(node_name) +
                             "' stage is defined in this environment");
        set_reset(up, down, env_reset_policy(breakdown_reset));

        NodeFailure<T> nf;
        nf.node = node_name;
        nf.breakdown = breakdown;
        nf.down_service = down_service;
        nf.breakdown_reset = breakdown_reset;
        if (has_repair) {
            set_reset(down, up, env_reset_policy(repair_reset));
            nf.repair = repair;
            nf.has_repair = true;
            nf.repair_reset = repair_reset;
        }
        record_node_failure(nf);
    }

    /**
     * Port of `Environment.init()`.
     *
     * Superpose the outgoing transitions of each stage, read the embedded jump
     * chain off the per-destination rates, and solve the resulting semi-Markov
     * process for its stationary stage probabilities.
     */
    void init() {
        const std::size_t E = stages_.size();
        for (std::size_t e = 0; e < E; ++e)
            if (!stages_[e].declared())
                throw InputError("Environment: stage " + std::to_string(e + 1) + " has no model");

        // proc[e][h], the transition process, as an MMAP marked by destination.
        proc.assign(E, std::vector<mam::Mmap<double>>(E));
        for (std::size_t e = 0; e < E; ++e)
            for (std::size_t h = 0; h < E; ++h) proc[e][h] = arc_mmap(e, h, E);

        hold_time.assign(E, mam::Mmap<double>());
        Matrix<double> Pemb(E, E, 0.0);
        std::vector<double> lambda(E, 0.0);
        for (std::size_t e = 0; e < E; ++e) {
            // The reference seeds the superposition with the SELF transition
            // and folds in every other destination.
            mam::Mmap<double> ht = proc[e][e];
            for (std::size_t h = 0; h < E; ++h) {
                if (h == e) continue;
                ht = superpose_collapsed(ht, proc[e][h]);
            }
            hold_time[e] = ht;
            const std::vector<double> cl = mam::mmap_count_lambda(ht);
            double tot = 0.0;
            for (double v : cl) tot += v;
            if (tot > 0.0)
                for (std::size_t h = 0; h < E; ++h) Pemb(e, h) = cl[h] / tot;
            const double m = mam::map_mean(ht.map());
            lambda[e] = (m > 0.0) ? 1.0 / m : 0.0;
        }

        bool all_positive = true;
        for (double v : lambda)
            if (!(v > 0.0)) all_positive = false;
        if (!all_positive)
            throw UnsupportedError(
                "Environment: a stage has no finite holding time, so the environment is "
                "absorbing; the reference leaves that case unimplemented (Environment.init has "
                "no branch for it)");

        // A = -lambda_e (I - Pemb): the generator of the semi-Markov process
        // observed at its jumps, whose stationary law is the stage probability.
        Matrix<double> A(E, E, 0.0);
        for (std::size_t e = 0; e < E; ++e)
            for (std::size_t h = 0; h < E; ++h)
                A(e, h) = -lambda[e] * ((e == h ? 1.0 : 0.0) - Pemb(e, h));
        prob_env = mc::ctmc_solve_reducible(A).pi;

        prob_orig = Matrix<double>(E, E, 0.0);
        for (std::size_t e = 0; e < E; ++e) {
            for (std::size_t h = 0; h < E; ++h)
                prob_orig(h, e) = prob_env[h] * lambda[h] * Pemb(h, e);
            if (prob_env[e] > 0.0) {
                double s = 0.0;
                for (std::size_t h = 0; h < E; ++h) s += prob_orig(h, e);
                if (s > 0.0)
                    for (std::size_t h = 0; h < E; ++h) prob_orig(h, e) /= s;
            }
        }
        pemb = Pemb;
        rate = lambda;
    }

    /**
     * `getReliabilityTable`: MTTF, MTTR, MTBF and availability of an
     * environment built out of node breakdowns.
     *
     * It reads the UP and DOWN_* stages BY NAME, as the reference does, so it
     * applies to a breakdown environment and to nothing else; a stage set that
     * carries no such names is refused rather than answered about.
     *
     * The four are not independent readings of the same thing. MTTF is the
     * COMPETING-RISK mean of the arcs out of UP, i.e. one over the SUM of the
     * breakdown rates, so a second failing node shortens it. MTTR is the mean
     * repair time weighted by the CONDITIONAL probability of being in each DOWN
     * stage given that the system is down. Availability is read off probEnv and
     * not from MTTF/(MTTF+MTTR): the two agree for a Markovian environment and
     * the stationary law is the one that stays right when the arcs are not.
     */
    struct Reliability {
        double mttf = 0.0;          ///< mean time to failure, UP -> any DOWN
        double mttr = 0.0;          ///< mean time to repair, DOWN -> UP
        double mtbf = 0.0;          ///< MTTF + MTTR
        double availability = 0.0;  ///< stationary probability of being UP
    };

    Reliability reliability() const {
        if (prob_env.empty())
            throw InputError(
                "Environment: reliability metrics are read from probEnv; call init() first");
        const std::size_t E = stages_.size();
        const std::size_t up = find_stage("UP");
        if (up == E)
            throw InputError(
                "Environment: no UP stage, so there is nothing for a breakdown to leave; "
                "reliability metrics apply to an environment built with add_node_breakdown");
        std::vector<std::size_t> down;
        for (std::size_t e = 0; e < E; ++e)
            if (stages_[e].name.compare(0, 5, "DOWN_") == 0) down.push_back(e);
        if (down.empty())
            throw InputError(
                "Environment: no DOWN_<node> stage, so no node breaks down; reliability metrics "
                "apply to an environment built with add_node_breakdown");

        double lambda_total = 0.0;
        for (std::size_t h : down)
            if (arcs_[up][h].enabled) lambda_total += arc_rate(up, h);
        if (!(lambda_total > 0.0))
            throw InputError(
                "Environment: no UP -> DOWN transition, so the system never fails and its mean "
                "time to failure is not a number");

        std::vector<double> mu, p;
        for (std::size_t e : down)
            if (arcs_[e][up].enabled) {
                mu.push_back(arc_rate(e, up));
                p.push_back(prob_env[e]);
            }
        if (mu.empty())
            throw InputError(
                "Environment: no DOWN -> UP transition, so the system is never repaired and its "
                "mean time to repair is not a number");

        Reliability r;
        r.mttf = 1.0 / lambda_total;
        double ptot = 0.0;
        for (double v : p) ptot += v;
        r.mttr = 0.0;
        for (std::size_t i = 0; i < mu.size(); ++i)
            r.mttr += (ptot > 0.0 ? p[i] / ptot : 1.0 / static_cast<double>(mu.size())) / mu[i];
        r.mtbf = r.mttf + r.mttr;
        double pup = prob_env[up], pdown = 0.0;
        for (std::size_t e : down) pdown += prob_env[e];
        r.availability = (pup + pdown > 0.0) ? pup / (pup + pdown) : 0.0;
        return r;
    }

    // ---- what init() produced --------------------------------------------
    std::vector<std::vector<mam::Mmap<double>>> proc;  ///< proc[e][h]
    std::vector<mam::Mmap<double>> hold_time;          ///< holdTime[e]
    std::vector<double> prob_env;                      ///< probEnv
    Matrix<double> prob_orig;                          ///< probOrig(h, e)
    Matrix<double> pemb;                               ///< the embedded jump chain
    std::vector<double> rate;                          ///< 1/E[holding time]

private:
    void check(std::size_t e) const {
        if (e >= stages_.size()) throw InputError("Environment: stage index out of range");
    }

    /** `1 / env{e,h}.getMean()`: the arc rate the reliability reading uses. */
    double arc_rate(std::size_t e, std::size_t h) const {
        const double m = num_traits<T>::to_double(arcs_[e][h].dist.mean);
        if (!(m > 0.0))
            throw InputError("Environment: the transition " + stages_[e].name + " -> " +
                             stages_[h].name +
                             " has no positive mean, so it carries no rate to report");
        return 1.0 / m;
    }

    /** Record a descriptor, replacing the one this node already had. */
    void record_node_failure(const NodeFailure<T>& nf) {
        const std::size_t idx = find_node_failure(nf.node);
        if (idx < node_failures_.size())
            node_failures_[idx] = nf;
        else
            node_failures_.push_back(nf);
    }

    /**
     * `emmap{e}{h}`: the arc's process marked by destination -- its D1 in slot
     * h and zero in every other slot. A disabled arc is the 1 x 1 zero pair.
     */
    mam::Mmap<double> arc_mmap(std::size_t e, std::size_t h, std::size_t E) const {
        mam::Mmap<double> m;
        if (!arcs_[e][h].enabled) {
            m.D0 = Matrix<double>(1, 1, 0.0);
            m.D1 = Matrix<double>(1, 1, 0.0);
            m.Dc.assign(E, Matrix<double>(1, 1, 0.0));
            return m;
        }
        const lang::Distrib<T>& d = arcs_[e][h].dist;
        const std::size_t n = d.D0.rows();
        if (n == 0)
            throw InputError("Environment: the transition distribution has no representation");
        m.D0 = Matrix<double>(n, n, 0.0);
        m.D1 = Matrix<double>(n, n, 0.0);
        for (std::size_t a = 0; a < n; ++a)
            for (std::size_t b = 0; b < n; ++b) {
                m.D0(a, b) = num_traits<T>::to_double(d.D0(a, b));
                m.D1(a, b) = num_traits<T>::to_double(d.D1(a, b));
            }
        m.Dc.assign(E, Matrix<double>(n, n, 0.0));
        m.Dc[h] = m.D1;
        return m;
    }

    /**
     * One fold of the reference's superposition loop: Kronecker-sum the blocks,
     * then move each block's row sums into its FIRST column.
     *
     * The collapse is what makes the superposed holding time RENEW at every
     * jump: without it the phase surviving from the losing risks would carry
     * into the next stage, which is not the semi-Markov model the solver then
     * solves.
     */
    static mam::Mmap<double> superpose_collapsed(const mam::Mmap<double>& a,
                                                 const mam::Mmap<double>& b) {
        mam::Mmap<double> s;
        s.D0 = mam::krons(a.D0, b.D0);
        s.D1 = collapse(mam::krons(a.D1, b.D1));
        s.Dc.reserve(a.Dc.size());
        for (std::size_t c = 0; c < a.Dc.size(); ++c)
            s.Dc.push_back(collapse(mam::krons(a.Dc[c], b.Dc[c])));
        return mam::mmap_normalize(s);
    }

    /** Row sums into column one, the rest zeroed. */
    static Matrix<double> collapse(const Matrix<double>& X) {
        Matrix<double> Y(X.rows(), X.cols(), 0.0);
        for (std::size_t i = 0; i < X.rows(); ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < X.cols(); ++j) s += X(i, j);
            Y(i, 0) = s;
        }
        return Y;
    }

    std::string name_;
    std::vector<EnvStage<T>> stages_;
    std::vector<std::vector<EnvArc<T>>> arcs_;
    std::vector<NodeFailure<T>> node_failures_;
};

}  // namespace env
}  // namespace line

#endif  // LINE_LANG_QN_ENVIRONMENT_H
