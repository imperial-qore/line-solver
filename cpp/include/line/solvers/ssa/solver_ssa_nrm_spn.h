/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_SSA_SOLVER_SSA_NRM_SPN_H
#define LINE_SOLVERS_SSA_SOLVER_SSA_NRM_SPN_H

/**
 * @file
 * @ingroup line_solvers
 * SolverSSA on a STOCHASTIC PETRI NET: the port of `solver_ssa_nrm_spn`, the
 * sub-engine `solver_ssa_nrm.m:174` hands a net carrying Transition nodes to.
 *
 * WHY IT IS A SEPARATE ENGINE. The queueing reaction builder of
 * `solver_ssa_nrm.h` reads a state vector of JOBS PER (node, class, phase) and
 * a departure is the absorption of a service process. A Petri net has neither:
 * a Place holds TOKENS, a Transition mode is a reaction whose stoichiometry is
 * the arc incidence, and nothing in the net has a service phase. The reference
 * therefore branches on `any(sn.nodetype == NodeType.Transition)` before it
 * builds a single reaction, and so does this port.
 *
 * THE MAPPING IS EXACT, not an approximation. A Place holds a per-class token
 * count (one slot of the state vector), and a timed Transition mode is a
 * reaction: input (enabling) arcs consume, output (firing) arcs produce.
 * Enabling is a propensity gate -- every input place at or above its arc
 * weight, every inhibitor place strictly below its threshold -- and a
 * single-server mode fires at its exponential rate while an infinite or
 * k-server mode fires at that rate times its enabling degree. One firing
 * applies the stoichiometry once, which is the atomic GSPN firing the exact
 * CTMC, JMT and the standard GSPN tools all take.
 *
 * IMMEDIATE MODES ARE NOT REACTIONS. They fire in zero time, so they are
 * resolved by VANISHING-MARKING ELIMINATION: after every timed firing (and once
 * on the initial marking) every enabled immediate mode is fired -- highest
 * firing priority first and, among equal priority, drawn in proportion to
 * firing weight -- until the marking is tangible. The timed race therefore only
 * ever samples from tangible markings and the immediate modes consume no
 * simulated time.
 *
 * A SOURCE IS NOT A TRANSITION and needs a reaction of its own, or the Place it
 * feeds stays empty and the net deadlocks on the first draw. Splitting a
 * Poisson stream by independent routing probabilities yields independent
 * Poisson streams, so the edge of probability p carries rate lambda*p exactly;
 * the reaction has an EMPTY enabling set (a constant propensity) and deposits
 * one token into the routed Place slot.
 *
 * WHAT IT REFUSES, and each refusal is the reference's: a non-exponential timed
 * firing (representing an in-flight firing's phase would need per-mode phase
 * state the reaction network does not carry), a non-exponential Source arrival
 * into a Place, a Source that reaches no Place, an infinite initial marking, and
 * a marking-dependent firing rate -- the last under EVERY SSA method, since no
 * SSA engine applies the g(marking) multiplier and answering with the nominal
 * rate would be silently wrong. `spn_nrm_supported` is the predicate; the
 * `default` dispatch reads it to decide whether to prefer this engine, and the
 * explicit `nrm` arm reads it to refuse by name.
 *
 * THE SLOT LAYOUT IS FLAT, one slot per (node, class), where the queueing
 * engine's is per (node, class, PHASE). Nothing here has a phase: a Place holds
 * tokens rather than jobs in service, and the reference's own SPN path only
 * ever addresses `phOff(place, class) + 1`, the first phase slot of the pair.
 * The two layouts therefore agree on every slot this engine touches.
 *
 * NOT A SAMPLE-PATH TWIN OF THE REFERENCE. The random stream is this port's
 * MT19937 (`SsaRng`) and the draw order is the algorithm's, so a seeded run
 * matches the reference STATISTICALLY and never bit for bit, exactly as the
 * queueing NRM does.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/solvers/ssa/ssa_types.h"
#include "line/util/error.h"
#include "line/util/line_console.h"
#include "line/util/matrix.h"

namespace line {
namespace ssa {

namespace detail {

/**
 * One reaction of the net: a timed Transition mode, an immediate one, or a
 * Source arrival. The reference's `spnEmptyRx` record, field for field.
 */
struct SpnReaction {
    std::size_t node = 0;  ///< 1-based node the reaction belongs to
    std::size_t mode = 0;  ///< 1-based mode; 0 marks a Source arrival
    /** Stoichiometry over the (node, class) slots: consume negative, produce positive. */
    std::vector<double> S;
    std::vector<std::size_t> en_slot;   ///< input arcs: the slots read
    std::vector<double> en_w;           ///< and the tokens each one needs
    std::vector<std::size_t> inh_slot;  ///< inhibitor arcs: the slots read
    std::vector<double> inh_thr;        ///< and the count at which each blocks
    double base_rate = 0.0;             ///< exponential firing rate, sum(D1)
    double nservers = 1.0;              ///< concurrent firings the mode may run
    double weight = 1.0;                ///< immediate: weight among equal priorities
    double prio = 1.0;                  ///< immediate: firing priority
    std::vector<std::size_t> dep_slots;  ///< slots this reaction DEPOSITS into
};

/**
 * Enabling degree: how many concurrent firings the marking supports, the
 * minimum over the input arcs of floor(tokens / weight), zeroed by any active
 * inhibitor. A mode with no input arc is single-degree -- which is what makes a
 * Source arrival a constant-propensity reaction.
 */
inline double spn_en_degree(const std::vector<double>& n, const SpnReaction& rx) {
    for (std::size_t i = 0; i < rx.inh_slot.size(); ++i)
        if (n[rx.inh_slot[i]] >= rx.inh_thr[i]) return 0.0;
    if (rx.en_slot.empty()) return 1.0;
    double d = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < rx.en_slot.size(); ++i)
        d = std::min(d, std::floor(n[rx.en_slot[i]] / rx.en_w[i]));
    return d;
}

/**
 * Propensity of a timed mode: the exponential rate times the effective server
 * count, min(enabling degree, mode servers). A single-server mode therefore
 * fires at its rate whenever enabled and an infinite-server one at the rate
 * scaled by the enabling degree.
 */
inline double spn_propensity(const std::vector<double>& n, const SpnReaction& rx) {
    const double eff = std::min(spn_en_degree(n, rx), rx.nservers);
    return eff <= 0.0 ? 0.0 : rx.base_rate * eff;
}

/**
 * `ssa_firingdep_refusal.m`: the NRM cannot answer a net whose firing rates
 * depend on the marking, and says so by name.
 *
 * The NRM builds one CONSTANT-propensity reaction per timed mode, so the
 * g(marking) multiplier has no place in it and answering with the nominal rate
 * would be a wrong number under the caller's own method name.
 *
 * NARROWER THAN THE REFERENCE'S, DELIBERATELY. MATLAB raises this from
 * `solver_ssa_analyzer` before an engine is chosen, so it refuses under every
 * SSA method; this port refuses only the NRM, because its SERIAL engine reaches
 * the multiplier for real -- `state_events.h:3649` scales the mode's firing rate
 * by `firingdep[m](marking)` on the state it is leaving, which is exactly the
 * semantics SolverCTMC gives it. Raising here as well would withdraw an answer
 * this port computes correctly. The gate therefore reads it with `raise=false`,
 * so `default` FALLS BACK to the serial engine, and with `raise=true` only when
 * the caller named `nrm`.
 */
template <class T>
bool ssa_check_firingdep(const qn::NetworkStruct<T>& sn, bool raise = true) {
    for (typename std::map<std::size_t, qn::TransitionParam<T>>::const_iterator it =
             sn.transparam.begin();
         it != sn.transparam.end(); ++it) {
        const qn::TransitionParam<T>& tp = it->second;
        for (std::size_t m = 0; m < tp.nmodes && m < tp.firingdep.size(); ++m) {
            if (!tp.firingdep[m]) continue;
            if (!raise) return false;
            const std::size_t ind = it->first;
            throw UnsupportedError(
                "SolverSSA: transition '" +
                (ind >= 1 && ind <= sn.nodes.size() ? sn.nodes[ind - 1].name : std::string("?")) +
                "' mode " + std::to_string(m + 1) +
                " uses a marking-dependent firing rate (setFiringRateDependence), which the "
                "NRM does not support: it builds one constant-propensity reaction per timed "
                "mode and never evaluates the handle. Ask for method='serial', which applies "
                "it, or use SolverCTMC or SolverLDES");
        }
    }
    return true;
}

/**
 * `ssa_nrm_guards.m`'s `g.spn`: the shapes the SPN reaction builder can read.
 *
 * FOUR CONDITIONS, and all four are the reference's. Each one the reference
 * RAISES from inside the builder, and each is a condition the serial engine
 * serves, so `default` falls back rather than refusing -- which is why this is
 * asked with `raise = false` by the eligibility test and with `raise = true`
 * only when the caller named `nrm` itself.
 */
template <class T>
bool spn_nrm_supported(const qn::NetworkStruct<T>& sn, bool raise = true) {
    const std::size_t I = sn.nodes.size();
    const std::size_t K = sn.nclasses;
    bool any_transition = false;
    for (std::size_t i = 0; i < I; ++i)
        if (sn.nodes[i].nodetype == qn::NodeType::Transition) any_transition = true;
    if (!any_transition) return true;
    if (!ssa_check_firingdep(sn, raise)) return false;

    for (std::size_t i = 0; i < I; ++i) {
        if (sn.nodes[i].nodetype == qn::NodeType::Transition) {
            const typename std::map<std::size_t, qn::TransitionParam<T>>::const_iterator it =
                sn.transparam.find(i + 1);
            if (it == sn.transparam.end()) continue;
            const qn::TransitionParam<T>& tp = it->second;
            for (std::size_t m = 0; m < tp.nmodes; ++m) {
                if (m < tp.timing.size() && tp.timing[m] == lang::TimingStrategy::IMMEDIATE)
                    continue;
                const bool one_phase = m < tp.firingphases.size() && tp.firingphases[m] == 1;
                const bool has_proc = m < tp.firingproc.size() && tp.firingproc[m].D1.rows() > 0;
                if (one_phase && has_proc) continue;
                if (!raise) return false;
                throw UnsupportedError(
                    "SolverSSA(method='nrm'): transition '" + sn.nodes[i].name + "' mode " +
                    std::to_string(m + 1) +
                    " has a non-exponential timed firing. The reaction network carries one "
                    "constant-rate reaction per mode and no per-mode phase, so an in-flight "
                    "firing has nowhere to keep its phase; use method='serial'");
            }
        } else if (sn.nodes[i].nodetype == qn::NodeType::Source) {
            const std::size_t ist = sn.nodes[i].station;
            if (ist == 0) continue;
            for (std::size_t r = 0; r < K; ++r) {
                const double lambda = num_traits<T>::to_double(sn.rates(ist - 1, r));
                if (std::isnan(lambda) || lambda <= 0.0) continue;
                if (sn.service[ist - 1][r].type != lang::ProcessType::EXP) {
                    if (!raise) return false;
                    throw UnsupportedError(
                        "SolverSSA(method='nrm'): source '" + sn.nodes[i].name + "' class '" +
                        sn.classes[r].name +
                        "' has a non-exponential arrival. An arrival into a Place is one "
                        "constant-propensity reaction here, which only a Poisson stream is; "
                        "use method='serial'");
                }
                bool feeds = false;
                for (std::size_t j = 0; j < I && !feeds; ++j) {
                    if (sn.nodes[j].nodetype != qn::NodeType::Place) continue;
                    for (std::size_t s = 0; s < K && !feeds; ++s)
                        if (num_traits<T>::to_double(sn.rtnodes(i * K + r, j * K + s)) > 0.0)
                            feeds = true;
                }
                if (!feeds) {
                    if (!raise) return false;
                    throw UnsupportedError(
                        "SolverSSA(method='nrm'): source '" + sn.nodes[i].name + "' class '" +
                        sn.classes[r].name +
                        "' routes to no Place. The SPN path deposits an arrival into a Place "
                        "slot and has nowhere else to put one; use method='serial'");
                }
            }
        } else if (sn.nodes[i].nodetype == qn::NodeType::Place) {
            const typename std::map<std::size_t, std::vector<T>>::const_iterator im =
                sn.initmarking.find(i + 1);
            if (im == sn.initmarking.end()) continue;
            for (std::size_t r = 0; r < im->second.size(); ++r) {
                if (std::isfinite(num_traits<T>::to_double(im->second[r]))) continue;
                if (!raise) return false;
                throw UnsupportedError("SolverSSA(method='nrm'): place '" + sn.nodes[i].name +
                                       "' declares an infinite initial marking, which the "
                                       "reaction network cannot count; use method='serial'");
            }
        }
    }
    return true;
}

}  // namespace detail

/**
 * The SPN Next-Reaction-Method engine, `solver_ssa_nrm_spn` of the reference.
 *
 * DOUBLE ONLY, for the reason `NrmEngine` is: the sample path is generated from
 * exponential clocks, which are logarithms of uniform draws, so there is no
 * exact value a wider arithmetic could carry. The analyzer refuses a non-double
 * backend by name before this class is instantiated.
 */
template <class T>
class NrmSpnEngine {
public:
    NrmSpnEngine(const qn::NetworkStruct<T>& sn, const SsaOptions& opt)
        : sn_(sn), opt_(opt), rng_(opt.seed) {
        build();
    }

    /** Run `opt.samples` firings and return the time-averaged metrics. */
    SsaSolution run();

    /** The reaction count, for the tests that assert the builder. */
    std::size_t nreactions() const { return rx_.size(); }
    /** The immediate-mode count, likewise. */
    std::size_t nimmediate() const { return imm_.size(); }

private:
    using Rx = detail::SpnReaction;

    const qn::NetworkStruct<T>& sn_;
    SsaOptions opt_;
    SsaRng rng_;

    std::size_t I_ = 0, K_ = 0, M_ = 0, NS_ = 0;

    std::vector<Rx> rx_;   ///< the timed reactions, the ones that race
    std::vector<Rx> imm_;  ///< the immediate modes, resolved between races
    /** Timed reactions consuming from each (0-based node, class): the Place throughput. */
    std::vector<std::vector<std::vector<std::size_t>>> consumers_;
    /** Arrival reactions injecting at each (0-based Source node, class). */
    std::vector<std::vector<std::vector<std::size_t>>> producers_;

    std::vector<double> nvec0_;  ///< the initial marking, before the collapse
    std::vector<double> pcap_slot_;  ///< per-(place, class) capacity, infinite when none
    /** Per-place total capacities: the bound and the slots it is taken over. */
    std::vector<std::pair<double, std::vector<std::size_t>>> place_total_caps_;
    bool has_caps_ = false;

    /** The livelock guard of the vanishing-marking collapse, the reference's 1e5. */
    static const std::size_t kMaxImmSteps = 100000;

    std::size_t slot(std::size_t node0, std::size_t cls) const { return node0 * K_ + cls; }

    void build();
    Rx build_mode(std::size_t ind0, std::size_t m) const;
    void apply_caps(std::vector<double>& n, const std::vector<std::size_t>& deposited) const;
    void collapse(std::vector<double>& n);
};

/**
 * The reaction list, the initial marking and the capacity tables.
 *
 * ONE PASS PER NODE KIND, in the reference's order: the Transition modes become
 * reactions (timed) or immediate records, then the Source arrivals become
 * constant-propensity reactions, then the Places contribute the marking and
 * their capacities.
 */
template <class T>
void NrmSpnEngine<T>::build() {
    I_ = sn_.nodes.size();
    K_ = sn_.nclasses;
    M_ = sn_.nstations;
    NS_ = I_ * K_;
    consumers_.assign(I_, std::vector<std::vector<std::size_t>>(K_));
    producers_.assign(I_, std::vector<std::vector<std::size_t>>(K_));

    // Refused here as well as in the gate, as the reference raises it from
    // inside its own builder too (`solver_ssa_nrm.m:3426`): that is the path a
    // caller reaches with the checks disabled, and this engine must not run a
    // net whose rates it would silently ignore however it was entered.
    detail::ssa_check_firingdep(sn_, true);

    for (std::size_t ind = 0; ind < I_; ++ind) {
        if (sn_.nodes[ind].nodetype != qn::NodeType::Transition) continue;
        const typename std::map<std::size_t, qn::TransitionParam<T>>::const_iterator it =
            sn_.transparam.find(ind + 1);
        if (it == sn_.transparam.end()) continue;
        const qn::TransitionParam<T>& tp = it->second;
        for (std::size_t m = 0; m < tp.nmodes; ++m) {
            const Rx rec = build_mode(ind, m);
            if (m < tp.timing.size() && tp.timing[m] == lang::TimingStrategy::IMMEDIATE) {
                imm_.push_back(rec);
                continue;
            }
            rx_.push_back(rec);
            const std::size_t ridx = rx_.size() - 1;
            for (std::size_t a = 0; a < rec.en_slot.size(); ++a)
                consumers_[rec.en_slot[a] / K_][rec.en_slot[a] % K_].push_back(ridx);
        }
    }

    // Source arrivals. Thinning a Poisson stream by the independent routing
    // probabilities leaves independent Poisson streams, so one reaction per
    // routed (Source, class) -> (Place, class) edge carries rate lambda*p
    // exactly. `producers_` is what makes the Source station report that rate as
    // its throughput, which is the open class's reference-station throughput.
    for (std::size_t ind = 0; ind < I_; ++ind) {
        if (sn_.nodes[ind].nodetype != qn::NodeType::Source) continue;
        const std::size_t ist = sn_.nodes[ind].station;
        if (ist == 0) continue;
        for (std::size_t r = 0; r < K_; ++r) {
            const double lambda = num_traits<T>::to_double(sn_.rates(ist - 1, r));
            if (std::isnan(lambda) || lambda <= 0.0) continue;
            if (sn_.service[ist - 1][r].type != lang::ProcessType::EXP)
                throw UnsupportedError(
                    "solver_ssa_nrm_spn: source '" + sn_.nodes[ind].name + "' class '" +
                    sn_.classes[r].name +
                    "' has a non-exponential arrival, which the SPN reaction network cannot "
                    "express; use method='serial' or SolverJMT");
            bool found_place = false;
            for (std::size_t jnd = 0; jnd < I_; ++jnd) {
                if (sn_.nodes[jnd].nodetype != qn::NodeType::Place) continue;
                for (std::size_t s = 0; s < K_; ++s) {
                    const double p =
                        num_traits<T>::to_double(sn_.rtnodes(ind * K_ + r, jnd * K_ + s));
                    if (p <= 0.0) continue;
                    found_place = true;
                    Rx rec;
                    rec.node = ind + 1;
                    rec.mode = 0;
                    rec.S.assign(NS_, 0.0);
                    rec.S[slot(jnd, s)] += 1.0;
                    rec.base_rate = lambda * p;
                    rec.nservers = 1.0;
                    rec.dep_slots.push_back(slot(jnd, s));
                    rx_.push_back(rec);
                    producers_[ind][r].push_back(rx_.size() - 1);
                }
            }
            if (!found_place)
                throw UnsupportedError("solver_ssa_nrm_spn: source '" + sn_.nodes[ind].name +
                                       "' class '" + sn_.classes[r].name +
                                       "' does not route to any Place; the SPN path needs a "
                                       "Source->Place arc");
        }
    }

    if (rx_.empty())
        throw InputError(
            "solver_ssa_nrm_spn: the stochastic Petri net has no timed reaction; there is "
            "nothing to simulate");

    // THE INITIAL MARKING, in the reference's two steps. MATLAB reads it off
    // `sn.state`, which `initDefault` fills by putting each closed class's whole
    // population at its reference station and which `Place.setState` then
    // overrides. This port has no State package -- the same position
    // `NrmEngine::build_initial_state` takes, and for the same reason -- so the
    // two steps are taken directly: a closed class whose reference station IS a
    // Place seeds that Place with its population, and a DECLARED `initmarking`
    // then replaces that Place's whole row, because an explicit marking is a
    // statement about every class of the place and not an addition to one.
    //
    // WITHOUT THE FIRST STEP a closed Petri net starts empty and deadlocks on the
    // first draw, since `initmarking` is only written by an explicit
    // `set_initial_marking` and a net that declares its tokens through a
    // ClosedClass population carries none.
    nvec0_.assign(NS_, 0.0);
    for (std::size_t r = 0; r < K_; ++r) {
        const double pop = sn_.classes[r].population;
        if (!std::isfinite(pop) || pop <= 0.0) continue;
        const std::size_t rs = sn_.classes[r].refstat;
        if (rs < 1 || rs > M_) continue;
        const std::size_t ind = sn_.station_to_node[rs - 1] - 1;
        if (sn_.nodes[ind].nodetype != qn::NodeType::Place) continue;
        nvec0_[slot(ind, r)] = pop;
    }
    for (std::size_t ind = 0; ind < I_; ++ind) {
        if (sn_.nodes[ind].nodetype != qn::NodeType::Place) continue;
        const typename std::map<std::size_t, std::vector<T>>::const_iterator im =
            sn_.initmarking.find(ind + 1);
        if (im == sn_.initmarking.end()) continue;
        for (std::size_t r = 0; r < K_; ++r) {
            const double v =
                r < im->second.size() ? num_traits<T>::to_double(im->second[r]) : 0.0;
            if (!std::isfinite(v))
                throw UnsupportedError("solver_ssa_nrm_spn: place '" + sn_.nodes[ind].name +
                                       "' declares an infinite initial marking, which the "
                                       "reaction network cannot count");
            nvec0_[slot(ind, r)] = v > 0.0 ? v : 0.0;
        }
    }

    // FINITE-CAPACITY PLACES DROP. A Place with a finite per-class (`classcap`)
    // or total (`cap`) capacity loses any arriving token that would exceed it,
    // which is the loss semantics JMT and the exact CTMC both give it: an
    // M/M/1/1 Place at rho = 0.5 holds a mean 1/3, not the unbounded 1. Without
    // the clamp the deposit accumulates past the bound.
    const double inf = std::numeric_limits<double>::infinity();
    pcap_slot_.assign(NS_, inf);
    for (std::size_t ind = 0; ind < I_; ++ind) {
        if (sn_.nodes[ind].nodetype != qn::NodeType::Place) continue;
        const std::size_t ist = sn_.nodes[ind].station;
        if (ist == 0) continue;
        std::vector<std::size_t> slots_here;
        for (std::size_t r = 0; r < K_; ++r) {
            slots_here.push_back(slot(ind, r));
            if (ist - 1 < sn_.classcap.size() && r < sn_.classcap[ist - 1].size()) {
                const double cc = sn_.classcap[ist - 1][r];
                if (std::isfinite(cc)) pcap_slot_[slot(ind, r)] = cc;
            }
        }
        if (ist - 1 < sn_.cap.size() && std::isfinite(sn_.cap[ist - 1]))
            place_total_caps_.push_back(std::make_pair(sn_.cap[ist - 1], slots_here));
    }
    has_caps_ = !place_total_caps_.empty();
    for (std::size_t j = 0; j < NS_ && !has_caps_; ++j)
        if (std::isfinite(pcap_slot_[j])) has_caps_ = true;
}

/**
 * `spnBuildMode`: the reaction record of transition `ind0` mode `m`.
 *
 * The enabling, firing and inhibiting matrices are (nnodes x nclasses), so a
 * nonzero entry (p, c) is an arc on class c of place p and its slot is the
 * (p, c) pair. An inhibiting entry is an arc only when it is FINITE: infinity
 * is how "this place never blocks the mode" is written, here as in the
 * reference.
 */
template <class T>
typename NrmSpnEngine<T>::Rx NrmSpnEngine<T>::build_mode(std::size_t ind0, std::size_t m) const {
    const qn::TransitionParam<T>& tp = sn_.transparam.at(ind0 + 1);
    Rx rec;
    rec.node = ind0 + 1;
    rec.mode = m + 1;
    rec.S.assign(NS_, 0.0);

    if (m < tp.enabling.size()) {
        const Matrix<T>& en = tp.enabling[m];
        for (std::size_t p = 0; p < en.rows() && p < I_; ++p)
            for (std::size_t c = 0; c < en.cols() && c < K_; ++c) {
                const double w = num_traits<T>::to_double(en(p, c));
                if (w == 0.0) continue;
                rec.en_slot.push_back(slot(p, c));
                rec.en_w.push_back(w);
                rec.S[slot(p, c)] -= w;
            }
    }
    if (m < tp.firing.size()) {
        const Matrix<T>& fir = tp.firing[m];
        for (std::size_t p = 0; p < fir.rows() && p < I_; ++p)
            for (std::size_t c = 0; c < fir.cols() && c < K_; ++c) {
                const double w = num_traits<T>::to_double(fir(p, c));
                if (w == 0.0) continue;
                rec.S[slot(p, c)] += w;
            }
    }
    if (m < tp.inhibiting.size()) {
        const Matrix<T>& inh = tp.inhibiting[m];
        for (std::size_t p = 0; p < inh.rows() && p < I_; ++p)
            for (std::size_t c = 0; c < inh.cols() && c < K_; ++c) {
                const double thr = num_traits<T>::to_double(inh(p, c));
                if (!std::isfinite(thr)) continue;
                rec.inh_slot.push_back(slot(p, c));
                rec.inh_thr.push_back(thr);
            }
    }
    for (std::size_t j = 0; j < NS_; ++j)
        if (rec.S[j] > 0.0) rec.dep_slots.push_back(j);

    // The exponential firing rate is the single-phase completion rate sum(D1).
    // A non-exponential firing is refused by `spn_nrm_supported` upstream and
    // again here, so no route reaches the run loop with a rate it invented.
    if (!(m < tp.timing.size() && tp.timing[m] == lang::TimingStrategy::IMMEDIATE)) {
        const bool one_phase = m < tp.firingphases.size() && tp.firingphases[m] == 1;
        if (!one_phase || m >= tp.firingproc.size() || tp.firingproc[m].D1.rows() == 0)
            throw UnsupportedError("solver_ssa_nrm_spn: transition '" + sn_.nodes[ind0].name +
                                   "' mode " + std::to_string(m + 1) +
                                   " has a non-exponential firing, which the SPN reaction "
                                   "network cannot express; use method='serial'");
        const Matrix<T>& d1 = tp.firingproc[m].D1;
        double s = 0.0;
        for (std::size_t a = 0; a < d1.rows(); ++a)
            for (std::size_t b = 0; b < d1.cols(); ++b) s += num_traits<T>::to_double(d1(a, b));
        rec.base_rate = s;
    }
    rec.nservers = m < tp.nmodeservers.size() ? tp.nmodeservers[m] : 1.0;
    if (std::isinf(rec.nservers)) rec.nservers = lang::GlobalConstants::MaxInt;
    rec.weight = m < tp.fireweight.size() ? num_traits<T>::to_double(tp.fireweight[m]) : 1.0;
    rec.prio = m < tp.firingprio.size() ? tp.firingprio[m] : 1.0;
    return rec;
}

/**
 * `applyPlaceCaps`: drop the tokens a firing pushed above a Place's bound.
 *
 * Only the just-deposited slots can overflow, so the clamp is local to them --
 * which also keeps it from disturbing a marking that was already over a bound
 * when the engine started.
 */
template <class T>
void NrmSpnEngine<T>::apply_caps(std::vector<double>& n,
                                 const std::vector<std::size_t>& deposited) const {
    for (std::size_t a = 0; a < deposited.size(); ++a) {
        const std::size_t j = deposited[a];
        if (n[j] > pcap_slot_[j]) n[j] = pcap_slot_[j];
    }
    for (std::size_t p = 0; p < place_total_caps_.size(); ++p) {
        const double tcap = place_total_caps_[p].first;
        const std::vector<std::size_t>& slots = place_total_caps_[p].second;
        double total = 0.0;
        for (std::size_t a = 0; a < slots.size(); ++a) total += n[slots[a]];
        double excess = total - tcap;
        for (std::size_t a = 0; a < deposited.size() && excess > 0.0; ++a) {
            const std::size_t j = deposited[a];
            if (std::find(slots.begin(), slots.end(), j) == slots.end()) continue;
            if (n[j] <= 0.0) continue;
            const double d = std::min(excess, n[j]);
            n[j] -= d;
            excess -= d;
        }
    }
}

/**
 * `spnCollapse`: vanishing-marking elimination.
 *
 * Fire enabled immediate modes until the marking is tangible, highest firing
 * priority first and ties drawn in proportion to firing weight. These firings
 * take zero time and advance no clock, so the timed race only ever resumes from
 * a tangible marking.
 */
template <class T>
void NrmSpnEngine<T>::collapse(std::vector<double>& n) {
    if (imm_.empty()) return;
    std::size_t steps = 0;
    while (true) {
        std::vector<std::size_t> enabled;
        for (std::size_t m = 0; m < imm_.size(); ++m)
            if (detail::spn_en_degree(n, imm_[m]) >= 1.0) enabled.push_back(m);
        if (enabled.empty()) return;
        double top_prio = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < enabled.size(); ++i)
            top_prio = std::max(top_prio, imm_[enabled[i]].prio);
        std::vector<std::size_t> top;
        std::vector<double> w;
        for (std::size_t i = 0; i < enabled.size(); ++i)
            if (imm_[enabled[i]].prio == top_prio) {
                top.push_back(enabled[i]);
                w.push_back(imm_[enabled[i]].weight);
            }
        const std::size_t pick = top.size() == 1 ? top[0] : top[rng_.draw(w)];
        for (std::size_t j = 0; j < NS_; ++j) n[j] += imm_[pick].S[j];
        if (++steps > kMaxImmSteps)
            throw NumericError(
                "solver_ssa_nrm_spn: immediate-transition livelock -- the vanishing-marking "
                "collapse did not reach a tangible marking");
    }
}

/**
 * The Next-Reaction-Method run loop over the tangible markings.
 *
 * EVERY PROPENSITY IS REFRESHED after a firing rather than a dependency subset:
 * a firing plus its immediate cascade can change any place, and the reaction
 * count of a Petri net is small, so refreshing all of them removes any
 * dependency-graph blind spot at no cost worth measuring.
 *
 * A PLACE IS AN INF STATION, so its utilization IS its mean token count -- the
 * SPN convention the CTMC analyzer also reports -- and its throughput is the
 * summed firing rate of the modes consuming from it, counted once per firing.
 */
template <class T>
SsaSolution NrmSpnEngine<T>::run() {
    const std::size_t nrx = rx_.size();
    SsaSolution out;
    out.QN = Matrix<double>(M_, K_, 0.0);
    out.UN = Matrix<double>(M_, K_, 0.0);
    out.RN = Matrix<double>(M_, K_, 0.0);
    out.TN = Matrix<double>(M_, K_, 0.0);
    out.CN.assign(K_, 0.0);
    out.XN.assign(K_, 0.0);
    out.StartN = Matrix<double>(M_, K_, 0.0);
    out.PreemptN = Matrix<double>(M_, K_, 0.0);
    out.method = "nrm";

    std::vector<double> nvec = nvec0_;
    collapse(nvec);
    if (has_caps_) {
        std::vector<std::size_t> all(NS_);
        for (std::size_t j = 0; j < NS_; ++j) all[j] = j;
        apply_caps(nvec, all);
    }

    std::vector<double> Ak(nrx, 0.0), Pk(nrx, 0.0), Tk(nrx, 0.0), tau(nrx, 0.0);
    const double inf = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < nrx; ++k) {
        Ak[k] = detail::spn_propensity(nvec, rx_[k]);
        Pk[k] = -std::log(rng_.uniform());
        tau[k] = Ak[k] > 0.0 ? (Pk[k] - Tk[k]) / Ak[k] : inf;
    }

    double total_time = 0.0;
    std::size_t n = 0;
    line::util::LineConsole::loop("drawing the sample path: %zu samples requested",
                                  static_cast<std::size_t>(opt_.samples));
    const std::size_t console_every = std::max<std::size_t>(1, opt_.samples / 20);
    for (; n < opt_.samples; ++n) {
        if ((n + 1) % console_every == 0)
            line::util::LineConsole::iter(
                static_cast<long>((n + 1) / console_every),
                "simulated %zu of %zu samples (%.0f%%), simulated time %.4g", n + 1,
                static_cast<std::size_t>(opt_.samples),
                100.0 * static_cast<double>(n + 1) / static_cast<double>(opt_.samples),
                total_time);
        std::size_t kfire = 0;
        double dt = inf;
        for (std::size_t k = 0; k < nrx; ++k)
            if (tau[k] < dt) {
                dt = tau[k];
                kfire = k;
            }
        if (std::isinf(dt))
            throw NumericError(
                "solver_ssa_nrm_spn: deadlock -- no transition is enabled, so the sample path "
                "cannot advance");
        total_time += dt;

        for (std::size_t ist = 0; ist < M_; ++ist) {
            const std::size_t ind = sn_.station_to_node[ist] - 1;
            for (std::size_t c = 0; c < K_; ++c) {
                const double tokens = nvec[slot(ind, c)];
                out.QN(ist, c) += tokens * dt;
                out.UN(ist, c) += tokens * dt;
                double depr = 0.0;
                const std::vector<std::size_t>& cons = consumers_[ind][c];
                for (std::size_t a = 0; a < cons.size(); ++a) depr += Ak[cons[a]];
                // A Source has no consuming transition: its throughput is the
                // aggregate arrival rate it injects, which is what makes the
                // reference station report the open class's arrival rate.
                const std::vector<std::size_t>& prod = producers_[ind][c];
                for (std::size_t a = 0; a < prod.size(); ++a) depr += Ak[prod[a]];
                out.TN(ist, c) += depr * dt;
            }
        }

        // One atomic firing, then the immediate cascade the new marking enabled.
        // A capacity-bound Place loses the tokens pushed above its bound BEFORE
        // the cascade sees the marking.
        for (std::size_t j = 0; j < NS_; ++j) nvec[j] += rx_[kfire].S[j];
        if (has_caps_) apply_caps(nvec, rx_[kfire].dep_slots);
        collapse(nvec);

        // Advance the Gibson and Bruck clocks with the PRE-firing propensities,
        // then refresh every propensity from the new marking.
        for (std::size_t k = 0; k < nrx; ++k) Tk[k] += Ak[k] * dt;
        for (std::size_t k = 0; k < nrx; ++k) Ak[k] = detail::spn_propensity(nvec, rx_[k]);
        Pk[kfire] -= std::log(rng_.uniform());
        for (std::size_t k = 0; k < nrx; ++k)
            tau[k] = Ak[k] > 0.0 ? (Pk[k] - Tk[k]) / Ak[k] : inf;
    }

    if (total_time > 0.0)
        for (std::size_t ist = 0; ist < M_; ++ist)
            for (std::size_t c = 0; c < K_; ++c) {
                out.QN(ist, c) /= total_time;
                out.UN(ist, c) /= total_time;
                out.TN(ist, c) /= total_time;
            }
    for (std::size_t c = 0; c < K_; ++c) {
        out.XN[c] = out.TN(sn_.classes[c].refstat - 1, c);
        for (std::size_t ist = 0; ist < M_; ++ist)
            out.RN(ist, c) = out.TN(ist, c) > 0.0 ? out.QN(ist, c) / out.TN(ist, c) : 0.0;
        if (out.XN[c] > 0.0) out.CN[c] = sn_.classes[c].population / out.XN[c];
    }
    out.simulated_time = total_time;
    out.samples = n;
    return out;
}

}  // namespace ssa
}  // namespace line

#endif  // LINE_SOLVERS_SSA_SOLVER_SSA_NRM_SPN_H
