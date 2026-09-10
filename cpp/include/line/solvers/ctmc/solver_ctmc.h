/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Port of `solver_ctmc.m`: the infinitesimal generator of a queueing network,
 * assembled from the enumerated state space and the synchronization list.
 *
 * THE ASSEMBLY. Every transition of the chain is one SYNCHRONIZATION: an active
 * event that sets the rate, and a passive event that says where the job lands.
 * For each (synchronization, state) pair the active handler is applied at its
 * node, the passive handler at its node, and the two local successors are
 * spliced back into a full network state whose index gives the column. The
 * generator entry is rate * probability, accumulated -- one (s, ns) pair can be
 * reached by several synchronizations, and each contributes.
 *
 * WHY THE DIAGONAL COMES LAST. Self-loops are generated deliberately (a lost
 * arrival is one), and they must cancel: `ctmc_makeinfgen` drops the diagonal
 * and then sets it to minus the row sum, so a self-loop contributes nothing to
 * the balance equations while the event still fired for rate-counting purposes.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_H

#include <cmath>
#include <cstddef>
#include <map>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/state.h"
#include "line/api/mam/map_moment.h"
#include "line/lang/qn/state_events.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/util/error.h"
#include "line/util/line_console.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace ctmc {

using qn::EventOutcome;
using qn::NetState;
using qn::NetworkStruct;
using qn::Sync;
using lang::EventType;
using lang::GlobalConstants;
using lang::NodeType;
using lang::SchedStrategy;

/** The generator, the state space it is indexed by, and the event rates. */
template <class T>
struct CtmcResult {
    Matrix<T> Q;                        ///< (n x n) infinitesimal generator
    std::vector<NetState<T>> space;     ///< row i of Q is space[i]
    /**
     * `arvRates` / `depRates`, indexed [state][stateful-1][class-1]: the total
     * rate of arrivals into, and departures out of, each stateful node in each
     * class, from each state. They are accumulated per synchronization rather
     * than read off Q, because Q has already summed the contributions of every
     * synchronization into one entry and they cannot be separated afterwards.
     */
    std::vector<std::vector<std::vector<T>>> arv_rates, dep_rates;
    /**
     * `Dfilt`, MATLAB's EVENT FILTRATION: `filt[a]` holds only the rates that
     * synchronization `a` contributed, so `sum_a filt[a]` is the off-diagonal
     * part of Q.
     *
     * IT CANNOT BE RECOVERED FROM Q, which is why it is carried rather than
     * recomputed: Q has already summed every synchronization's contribution
     * into one entry. The response-time CDF is built by splitting the generator
     * on ONE event -- the tagged job's arrival, then its departure -- into a
     * MAP (Q - D1, D1), and that split needs the per-event matrix.
     *
     * Empty unless `solver_ctmc` was asked for it: it costs one n x n matrix per
     * synchronization, which on a model with many routing pairs dwarfs Q itself.
     */
    std::vector<Matrix<T>> filt;
    /**
     * The DERIVED START and PREEMPT filtrations, indexed [station-1][class-1]:
     * the rate at which a transition starts a class-r service at station i, and
     * the rate at which it pushes a class-r job in service back into the buffer.
     *
     * They are NOT part of `filt`, which pairs one-to-one with the
     * synchronization list and whose sum is the off-diagonal of Q: a START rides
     * on the SAME arc as the ARV or DEP that causes it, so adding it there would
     * double-count the generator. Filled whenever `filt` is.
     */
    std::vector<std::vector<Matrix<T>>> start_filt, preempt_filt;
    /**
     * `Qimm`, the IMMEDIATE-ONLY part of Q: the arcs contributed by a Router or
     * Fork pass-through, by a Join firing on an FJ-augmented struct, and by an
     * SPN ENABLE or an IMMEDIATE-mode firing. Empty when the model has no such
     * source.
     *
     * It exists for the VANISHING-ROW PURGE. A zero-sojourn state is left the
     * instant it is entered, so a timed arc out of one describes an event that
     * cannot occur; the purge replaces such a row by its immediate part. Q has
     * already summed the two together, so the split cannot be recovered from it
     * afterwards -- the same reason `filt` is carried rather than recomputed.
     */
    Matrix<T> Qimm;
    /**
     * The parts of `arv_rates` / `dep_rates` contributed by those same immediate
     * sources. Kept alongside the totals for two reasons, both from
     * `solver_ctmc.m`: the purge has to restate a vanishing row's rates as its
     * immediate part, and the RATE COMPLEMENT applies to the immediate part
     * alone -- an event that fires only from vanishing states would otherwise be
     * lost when those rows are eliminated. Empty when `Qimm` is.
     */
    std::vector<std::vector<std::vector<T>>> arv_rates_imm, dep_rates_imm;
    /**
     * The rows the purge restated, i.e. the vanishing states that `Qimm` gave an
     * immediate exit. `ctmc_eliminate_vanishing` complements exactly these out;
     * carrying them avoids re-running the predicate, which costs one event
     * evaluation per (state, source).
     */
    std::vector<std::size_t> vanishing;
};

namespace ctmc_detail {

/** The flattened key of a network state, for exact index lookup. */
template <class T>
std::vector<double> state_key(const NetState<T>& ns) {
    std::vector<double> key;
    for (std::size_t i = 0; i < ns.local.size(); ++i) {
        // The separator keeps two different splits of the same concatenation
        // from colliding, which a plain flatten would allow.
        key.push_back(-2.0);
        for (std::size_t j = 0; j < ns.local[i].size(); ++j)
            key.push_back(num_traits<T>::to_double(ns.local[i][j]));
    }
    return key;
}

/**
 * The largest per-class count any single stateful node holds in a state.
 *
 * The peak and not the sum, because the bound it feeds is per node.
 *
 * A TRANSITION IS SKIPPED, not counted. Its row is per MODE -- idle servers,
 * firing phases, fired counts -- so the leading columns `to_marginal_aggr`
 * would read as class counts are mode counts, and an infinite-server mode
 * carries MaxInt there. It holds no jobs; the tokens are in the places.
 */
template <class T>
std::vector<double> state_peak_occupancy(const NetworkStruct<T>& sn, const NetState<T>& st) {
    std::vector<double> pk(sn.nclasses, 0.0);
    const std::vector<std::size_t>& sfn = sn.stateful_nodes;
    for (std::size_t f = 0; f < sfn.size() && f < st.local.size(); ++f) {
        if (sn.nodes[sfn[f] - 1].nodetype == NodeType::Transition) continue;
        const std::pair<T, std::vector<T>> mg = qn::to_marginal_aggr(sn, sfn[f], st.local[f]);
        for (std::size_t r = 0; r < sn.nclasses && r < mg.second.size(); ++r) {
            const double v = num_traits<T>::to_double(mg.second[r]);
            if (v > pk[r]) pk[r] = v;
        }
    }
    return pk;
}

/**
 * True when no stateful node holds more of an OPEN class than `lim` allows.
 *
 * PER NODE AND NOT IN TOTAL, which is the reference's `capacityc(ind,r)` of
 * `State.spaceGeneratorNodes`: an open class is capped at the cutoff AT EACH
 * node. A total bound is wrong for an SPN whose firings do not conserve tokens
 * -- `spn_open_sevenplaces` has a mode consuming one token and producing two --
 * because the sum then crosses the bound on a firing that no place overflows,
 * the arc vanishes, and the truncated chain absorbs at the boundary and reports
 * Tput 0. A per-node bound censors only where a place itself overflows.
 */
template <class T>
bool within_cutoff(const NetworkStruct<T>& sn, const NetState<T>& st,
                   const std::vector<double>& njobs, const std::vector<std::size_t>& lim,
                   const std::vector<std::vector<std::size_t>>& lim_mat =
                       std::vector<std::vector<std::size_t>>()) {
    const std::vector<std::size_t>& sfn = sn.stateful_nodes;
    for (std::size_t f = 0; f < sfn.size() && f < st.local.size(); ++f) {
        if (sn.nodes[sfn[f] - 1].nodetype == NodeType::Transition) continue;
        // The reference's cutoff may be a (station x class) MATRIX, in which
        // case the bound at THIS node is its own row rather than the per-class
        // maximum; a node with no station (Cache, Router) keeps the vector.
        const std::size_t ist = sn.nodes[sfn[f] - 1].station;
        const std::pair<T, std::vector<T>> mg = qn::to_marginal_aggr(sn, sfn[f], st.local[f]);
        for (std::size_t r = 0; r < sn.nclasses && r < njobs.size(); ++r) {
            if (std::isfinite(njobs[r])) continue;
            std::size_t bound = r < lim.size() ? lim[r] : 0;
            if (!lim_mat.empty() && ist != 0 && ist - 1 < lim_mat.size() &&
                r < lim_mat[ist - 1].size())
                bound = lim_mat[ist - 1][r];
            if (bound == 0) continue;
            if (r < mg.second.size() &&
                num_traits<T>::to_double(mg.second[r]) > static_cast<double>(bound))
                return false;
        }
    }
    return true;
}

/**
 * Accumulate the START/PREEMPT annotation of one successor row into the derived
 * filtrations of the station behind NODE. W is the same weight the caller added
 * to Q and to `filt`, so the filtration integrates rate * count and pi*F*e is a
 * rate of starts (or of preemptions) per unit time.
 */
template <class T, class R>
inline void add_aux_filt(R& res, const NetworkStruct<T>& sn, std::size_t node, std::size_t s,
                         std::size_t ns, const T& w, const qn::EventOutcome<T>& oc,
                         std::size_t row) {
    if (num_traits<T>::to_double(w) == 0) return;
    if (node == 0 || node > sn.nodes.size()) return;
    const std::size_t ist = sn.nodes[node - 1].station;
    if (ist == 0 || ist > res.start_filt.size()) return;
    if (row < oc.start.size())
        for (std::size_t j = 0; j < oc.start[row].size(); ++j) {
            const std::size_t cls = oc.start[row][j];
            if (cls >= 1 && cls <= res.start_filt[ist - 1].size())
                res.start_filt[ist - 1][cls - 1](s, ns) += w;
        }
    if (row < oc.preempt.size())
        for (std::size_t j = 0; j < oc.preempt[row].size(); ++j) {
            const std::size_t cls = oc.preempt[row][j];
            if (cls >= 1 && cls <= res.preempt_filt[ist - 1].size())
                res.preempt_filt[ist - 1][cls - 1](s, ns) += w;
        }
}

/**
 * Solve (-Q22) X = B, one elimination shared by every column.
 *
 * The vanishing machinery needs the SAME solve three times over -- inside the
 * complement, once per event filtration and once per rate vector -- so the
 * censored block is factorized once and back-substituted per right-hand side.
 * `ctmc_stochcomp` performs its own factorization for S; this is the one every
 * OTHER right-hand side rides on.
 */
template <class T>
Matrix<T> censored_solve(const Matrix<T>& Q22, const Matrix<T>& B) {
    const std::size_t nd = Q22.rows();
    if (B.rows() != nd) throw InputError("censored_solve: the right-hand side is misshapen");
    Matrix<T> A(nd, nd, num_traits<T>::from_int(0));
    for (std::size_t a = 0; a < nd; ++a)
        for (std::size_t b = 0; b < nd; ++b) A(a, b) = T(-Q22(a, b));
    const std::vector<std::size_t> piv = lu_factor(A);
    Matrix<T> X = B;
    std::vector<T> rhs(nd);
    for (std::size_t c = 0; c < B.cols(); ++c) {
        for (std::size_t d = 0; d < nd; ++d) rhs[d] = B(d, c);
        lu_solve(A, piv, rhs);
        for (std::size_t d = 0; d < nd; ++d) X(d, c) = rhs[d];
    }
    return X;
}

/**
 * Complement ONE filtration onto the tangible states:
 *   Dnew = D(nonimm, nonimm) + Q12 (-Q22)^-1 D(imm, nonimm)
 *
 * All right-hand sides of one matrix go through a single elimination, since the
 * censored block is the same for every filtration in the model.
 */
template <class T>
Matrix<T> complement_one_filt(const Matrix<T>& D, const std::vector<std::size_t>& nonimm,
                              const std::vector<std::size_t>& imm,
                              const mc::StochCompResult<T>& sc) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t nk = nonimm.size(), nd = imm.size();
    Matrix<T> out(nk, nk, zero);
    for (std::size_t a = 0; a < nk; ++a)
        for (std::size_t b = 0; b < nk; ++b) out(a, b) = D(nonimm[a], nonimm[b]);
    if (nd == 0) return out;
    Matrix<T> B(nd, nk, zero);
    bool any = false;
    for (std::size_t d = 0; d < nd; ++d)
        for (std::size_t b = 0; b < nk; ++b) {
            B(d, b) = D(imm[d], nonimm[b]);
            if (num_traits<T>::to_double(B(d, b)) != 0) any = true;
        }
    if (!any) return out;
    const Matrix<T> X = censored_solve(sc.Q22, B);
    for (std::size_t a = 0; a < nk; ++a)
        for (std::size_t b = 0; b < nk; ++b) {
            T acc = zero;
            for (std::size_t d = 0; d < nd; ++d) acc = T(acc + sc.Q12(a, d) * X(d, b));
            out(a, b) = T(out(a, b) + acc);
        }
    return out;
}

/** Apply `complement_one_filt` to the event, START and PREEMPT filtrations. */
template <class T, class R>
void complement_filtrations(R& res, const std::vector<std::size_t>& nonimm,
                            const std::vector<std::size_t>& imm,
                            const mc::StochCompResult<T>& sc) {
    for (std::size_t a = 0; a < res.filt.size(); ++a)
        res.filt[a] = complement_one_filt(res.filt[a], nonimm, imm, sc);
    for (std::size_t i = 0; i < res.start_filt.size(); ++i)
        for (std::size_t r = 0; r < res.start_filt[i].size(); ++r)
            res.start_filt[i][r] = complement_one_filt(res.start_filt[i][r], nonimm, imm, sc);
    for (std::size_t i = 0; i < res.preempt_filt.size(); ++i)
        for (std::size_t r = 0; r < res.preempt_filt[i].size(); ++r)
            res.preempt_filt[i][r] = complement_one_filt(res.preempt_filt[i][r], nonimm, imm, sc);
}

}  // namespace ctmc_detail

/**
 * Port of `ctmc_makeinfgen`: turn an off-diagonal rate matrix into a generator.
 *
 * The diagonal is discarded first and then set to minus the row sum, so any
 * self-loop that was accumulated cancels exactly.
 */
template <class T>
void make_infgen(Matrix<T>& Q) {
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < Q.rows(); ++i) Q(i, i) = zero;
    for (std::size_t i = 0; i < Q.rows(); ++i) {
        T s = zero;
        for (std::size_t j = 0; j < Q.cols(); ++j) s += Q(i, j);
        Q(i, i) = T(-s);
    }
}

template <class T>
Matrix<T> ctmc_state_space_aggr(const NetworkStruct<T>& sn, const std::vector<NetState<T>>& space);

/**
 * Tabulates the globally state-dependent rate scaling phi(n) declared through
 * `set_global_dependence`, ONE evaluation per state.
 *
 * Returns an (nstates x nstations*nclasses) matrix, row-major in (station,
 * class), of the scaling applying at each state; empty when the model declares
 * no global dependence. Evaluating once per state rather than per transition is
 * the whole point: phi may be expensive (a bandwidth-sharing allocation solves a
 * convex program per call), and within a state it is a CONSTANT multiplying every
 * rate there, which is why it factors out of the generator assembly below.
 */
template <class T>
Matrix<T> ctmc_gd_factor(const NetworkStruct<T>& sn, const std::vector<NetState<T>>& space) {
    const T zero = num_traits<T>::from_int(0);
    if (!static_cast<bool>(sn.gdscaling)) return Matrix<T>(0, 0, zero);
    if (!sn.regions.empty())
        throw InputError(
            "setGlobalDependence cannot be combined with finite capacity regions: the region "
            "generator builds its own transitions and would ignore the scaling");
    const std::size_t M = sn.stations.size(), K = sn.nclasses, n = space.size();
    const std::size_t max_entries = 30000000u;
    if (n * M * K > max_entries)
        throw InputError(
            "the global dependence table would exceed the state budget; lower the cutoff");
    const Matrix<T> aggr = ctmc_state_space_aggr(sn, space);
    Matrix<T> out(n, M * K, num_traits<T>::from_int(1));
    std::vector<T> npop(M * K, zero);
    for (std::size_t s = 0; s < n; ++s) {
        for (std::size_t i = 0; i < M * K; ++i) npop[i] = aggr(s, i);
        const std::vector<T> v = sn.gdscaling(npop);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) {
                const T f = v.size() == 1 ? v[0] : (v.size() == M ? v[i] : v[i * K + r]);
                if (!(num_traits<T>::to_double(f) >= 0))
                    throw InputError(
                        "the global dependence handle returned a non-finite or negative scaling");
                out(s, i * K + r) = f;
            }
    }
    return out;
}

/**
 * Port of `ctmc_find_vanishing_states` (`solver_ctmc.m:928`): the indices of the
 * VANISHING (zero-sojourn) global states.
 *
 * A state is vanishing when the model leaves it at the `GlobalConstants`
 * Immediate scale rather than at a modelled rate, so its sojourn is an artefact
 * of realising "instantaneous" as a very fast exponential. Four sources, which
 * are the whole list:
 *
 *  - a Router or Fork holding a job. Neither performs service: the job is in
 *    transit and leaves on the next event.
 *  - a Join whose sibling set is COMPLETE for some original class, so the
 *    rendezvous can fire. Only on an FJ-augmented struct -- without the tag
 *    classes a Join buffers nothing and its departures are timed elsewhere.
 *  - an SPN marking from which an ENABLE moves the Transition's own row.
 *  - an SPN marking enabling a `TimingStrategy::IMMEDIATE` firing mode.
 *
 * The same predicate drives BOTH the vanishing-row purge and the stochastic
 * complementation, which is why it is one function: a row purged as vanishing
 * and then left in the chain would carry only its immediate arcs and dominate
 * the stationary vector with a 1e-8 sojourn, and a row complemented out without
 * being purged would push its timed arcs into the tangible states.
 *
 * @return 0-based row indices, ascending and unique
 */
template <class T>
std::vector<std::size_t> ctmc_find_vanishing_states(
    const NetworkStruct<T>& sn, const std::vector<NetState<T>>& space,
    const std::vector<qn::GlobalSync<T>>& gsync, bool isfjaug) {
    const std::size_t n = space.size();
    const std::size_t R = sn.nclasses;
    std::vector<bool> mark(n, false);

    // Router / Fork pass-through occupancy.
    for (std::size_t ind = 1; ind <= sn.nodes.size(); ++ind) {
        const NodeType nt = sn.nodes[ind - 1].nodetype;
        if (nt != NodeType::Router && nt != NodeType::Fork) continue;
        if (sn.nodes[ind - 1].station != 0) continue;  // stateful and NOT a station
        const std::size_t isf = sn.stateful_index(ind);
        if (isf == 0) continue;
        for (std::size_t s = 0; s < n; ++s) {
            if (mark[s]) continue;
            const std::pair<T, std::vector<T>> mg =
                qn::to_marginal_aggr(sn, ind, space[s].local[isf - 1]);
            for (std::size_t r = 0; r < R && r < mg.second.size(); ++r) {
                const double v = num_traits<T>::to_double(mg.second[r]);
                if (std::isfinite(v) && v > 0.0) {
                    mark[s] = true;
                    break;
                }
            }
        }
    }

    // Join rendezvous ready to fire.
    if (isfjaug) {
        for (std::size_t ind = 1; ind <= sn.nodes.size(); ++ind) {
            if (sn.nodes[ind - 1].nodetype != NodeType::Join) continue;
            const std::size_t isf = sn.stateful_index(ind);
            if (isf == 0) continue;
            const typename std::map<std::size_t, qn::FjJoinParam>::const_iterator jit =
                sn.fjjoinparam.find(ind);
            if (jit == sn.fjjoinparam.end()) continue;
            const std::vector<std::size_t>& origcl = jit->second.origclasses;
            for (std::size_t s = 0; s < n; ++s) {
                if (mark[s]) continue;
                for (std::size_t x = 0; x < origcl.size(); ++x) {
                    const qn::EventOutcome<T> oj = qn::after_event_join(
                        sn, ind, space[s].local[isf - 1], EventType::DEP, origcl[x]);
                    if (!oj.space.empty()) {
                        mark[s] = true;
                        break;
                    }
                }
            }
        }
    }

    // SPN: an ENABLE that moves the Transition's row, and an IMMEDIATE firing.
    for (std::size_t g = 0; g < gsync.size(); ++g) {
        const qn::ModeEvent<T>& ae = gsync[g].active;
        const std::size_t isf_t = sn.stateful_index(ae.node);
        if (isf_t == 0) continue;
        bool is_enable = ae.event == EventType::ENABLE;
        bool is_imm_fire = false;
        if (!is_enable && ae.event == EventType::FIRE) {
            const typename std::map<std::size_t, qn::TransitionParam<T>>::const_iterator it =
                sn.transparam.find(ae.node);
            is_imm_fire = it != sn.transparam.end() && ae.mode >= 1 &&
                          ae.mode <= it->second.timing.size() &&
                          it->second.timing[ae.mode - 1] == lang::TimingStrategy::IMMEDIATE;
        }
        if (!is_enable && !is_imm_fire) continue;
        for (std::size_t s = 0; s < n; ++s) {
            if (mark[s]) continue;
            const qn::GlobalOutcome<T> go = qn::after_global_event(sn, space[s], gsync[g]);
            for (std::size_t io = 0; io < go.space.size(); ++io) {
                if (num_traits<T>::to_double(go.rate[io]) <= 0) continue;
                // AN ENABLE THAT LEAVES THE ROW ALONE IS NOT A MOVE. The
                // reference tests the Transition's OWN row rather than the whole
                // marking (`solver_ctmc.m:983`): an enabling that finds the mode
                // already enabled re-emits the same state and takes no time to
                // do nothing, which is not a vanishing state.
                if (is_imm_fire || go.space[io].local[isf_t - 1] != space[s].local[isf_t - 1]) {
                    mark[s] = true;
                    break;
                }
            }
        }
    }

    std::vector<std::size_t> imm;
    for (std::size_t s = 0; s < n; ++s)
        if (mark[s]) imm.push_back(s);
    return imm;
}

/**
 * Port of the generator assembly of `solver_ctmc.m`.
 *
 * @param sn     the network struct
 * @param space  the enumerated state space, from `space_generator`
 * @param sync   the synchronization list, from `refresh_sync`
 * @param gsync  the SPN global synchronizations, from `refresh_gsync`
 * @param fjsync the fork firing list, from `fj_tag`
 * @param want_filtration also return the per-synchronisation rate matrices (the filtration), which sampling and reward paths need
 */
template <class T>
CtmcResult<T> solver_ctmc(const NetworkStruct<T>& sn, const std::vector<NetState<T>>& space,
                          const std::vector<Sync<T>>& sync,
                          const std::vector<qn::GlobalSync<T>>& gsync =
                              std::vector<qn::GlobalSync<T>>(),
                          bool want_filtration = false,
                          const std::vector<qn::FjSync<T>>& fjsync =
                              std::vector<qn::FjSync<T>>()) {
    const std::size_t n = space.size();
    const std::size_t local = sn.nodes.size() + 1;  // the dummy passive node
    const T zero = num_traits<T>::from_int(0);
    CtmcResult<T> res;
    res.space = space;
    res.Q = Matrix<T>(n, n, zero);

    const std::size_t NF = sn.stateful_nodes.size();
    const std::size_t R = sn.nclasses;
    res.arv_rates.assign(n, std::vector<std::vector<T>>(NF, std::vector<T>(R, zero)));
    res.dep_rates.assign(n, std::vector<std::vector<T>>(NF, std::vector<T>(R, zero)));

    if (want_filtration) {
        res.filt.assign(sync.size(), Matrix<T>(n, n, zero));
        res.start_filt.assign(sn.nstations, std::vector<Matrix<T>>(R, Matrix<T>(n, n, zero)));
        res.preempt_filt.assign(sn.nstations, std::vector<Matrix<T>>(R, Matrix<T>(n, n, zero)));
    }

    // `immAction` / `immGsync` of `solver_ctmc.m:430-467`: which sources emit at
    // the GlobalConstants::Immediate scale, and therefore land in `Qimm`. A Join
    // counts only on an FJ-AUGMENTED struct -- without the tag classes a Join is
    // an ordinary pass-through whose departures are timed by the siblings.
    const bool isfjaug = !fjsync.empty();
    std::vector<bool> imm_action(sync.size(), false);
    bool has_imm = !fjsync.empty();
    for (std::size_t a = 0; a < sync.size(); ++a) {
        const std::size_t na = sync[a].active.node;
        if (na == 0 || na > sn.nodes.size()) continue;
        const NodeType nt = sn.nodes[na - 1].nodetype;
        imm_action[a] = nt == NodeType::Router || nt == NodeType::Fork ||
                        (isfjaug && nt == NodeType::Join);
        if (imm_action[a]) has_imm = true;
    }
    std::vector<bool> imm_gsync(gsync.size(), false);
    for (std::size_t g = 0; g < gsync.size(); ++g) {
        const qn::ModeEvent<T>& ae = gsync[g].active;
        if (ae.event == EventType::ENABLE) {
            imm_gsync[g] = true;
        } else if (ae.event == EventType::FIRE) {
            const typename std::map<std::size_t, qn::TransitionParam<T>>::const_iterator it =
                sn.transparam.find(ae.node);
            imm_gsync[g] = it != sn.transparam.end() && ae.mode >= 1 &&
                           ae.mode <= it->second.timing.size() &&
                           it->second.timing[ae.mode - 1] == lang::TimingStrategy::IMMEDIATE;
        }
        if (imm_gsync[g]) has_imm = true;
    }
    if (has_imm) {
        res.Qimm = Matrix<T>(n, n, zero);
        res.arv_rates_imm.assign(n, std::vector<std::vector<T>>(NF, std::vector<T>(R, zero)));
        res.dep_rates_imm.assign(n, std::vector<std::vector<T>>(NF, std::vector<T>(R, zero)));
    }

    // The true-BAS become-blocked edges, kept out of Q until the end because they
    // are not departures: they must reach the generator but not `dep_rates`.
    bool any_bas = false;
    for (std::size_t i = 0; i < sn.isbasblocking.size(); ++i)
        if (sn.isbasblocking[i]) any_bas = true;
    Matrix<T> bas_block(any_bas ? n : 0, any_bas ? n : 0, zero);

    std::map<std::vector<double>, std::size_t> index;
    for (std::size_t s = 0; s < n; ++s) index[ctmc_detail::state_key(space[s])] = s;

    // phi(n) is constant within a state, so it factors out of every rate there
    const Matrix<T> gd = ctmc_gd_factor(sn, space);
    const bool has_gd = gd.rows() != 0;

    // The state-dependent routing table, one per state. It is tabulated here for
    // the same reason phi(n) is: the loop below runs synchronization-outer and
    // state-inner, so evaluating eq. (10) inside it would redo one stochastic
    // complement per (sync, state) pair rather than one per state.
    std::vector<Matrix<T>> rt_by_state;
    if (sn.has_sdr_routing()) {
        rt_by_state.reserve(n);
        for (std::size_t s = 0; s < n; ++s) rt_by_state.push_back(qn::rt_state(sn, space[s].local));
    }

    for (std::size_t a = 0; a < sync.size(); ++a) {
        const Sync<T>& sy = sync[a];
        const std::size_t node_a = sy.active.node;
        const std::size_t isf_a = sn.stateful_index(node_a);
        if (isf_a == 0) continue;  // a stateless node schedules nothing
        const std::size_t node_p = sy.passive.node;
        const std::size_t isf_p = node_p == local ? 0 : sn.stateful_index(node_p);
        if (node_p != local && isf_p == 0) continue;

        // PHASE is scaled too, or phase-type service would advance unscaled
        const bool gd_here = has_gd && sn.nodes[node_a - 1].station != 0 &&
                             (sy.active.event == EventType::DEP ||
                              sy.active.event == EventType::PHASE);
        const std::size_t gd_col = gd_here ? (sn.nodes[node_a - 1].station - 1) * sn.nclasses +
                                                 (sy.active.cls - 1)
                                           : 0;
        // A round-robin dispatcher decides the destination from its own pointer;
        // see the proute branch below.
        const bool rr_here = sy.active.event == EventType::DEP &&
                             sn.rr_var_slot(node_a, sy.active.cls) != 0;
        for (std::size_t s = 0; s < n; ++s) {
            const NetState<T>& st = space[s];
            T fired = zero;  // the rate this synchronization contributes here
            const EventOutcome<T> oa =
                qn::after_event(sn, node_a, st.local[isf_a - 1], sy.active.event, sy.active.cls);
            for (std::size_t ia = 0; ia < oa.space.size(); ++ia) {
                const T rate = gd_here ? T(oa.rate[ia] * gd(s, gd_col)) : oa.rate[ia];
                // A zero-rate successor is a state the reference still emits so
                // that the event exists; it contributes nothing to the balance
                // equations, so skip it here rather than adding a zero.
                if (num_traits<T>::to_double(rate) == 0) continue;

                if (node_p == local) {
                    // A local action moves no job elsewhere: only the active
                    // node's block changes.
                    NetState<T> nsx = st;
                    nsx.local[isf_a - 1] = oa.space[ia];
                    const typename std::map<std::vector<double>, std::size_t>::const_iterator it =
                        index.find(ctmc_detail::state_key(nsx));
                    if (it == index.end()) continue;
                    res.Q(s, it->second) += T(rate * oa.prob[ia]);
                    if (imm_action[a]) res.Qimm(s, it->second) += T(rate * oa.prob[ia]);
                    if (want_filtration) {
                        res.filt[a](s, it->second) += T(rate * oa.prob[ia]);
                        // local action: only the active node can tag
                        ctmc_detail::add_aux_filt(res, sn, node_a, s, it->second,
                                                  T(rate * oa.prob[ia]), oa, ia);
                    }
                    fired += T(rate * oa.prob[ia]);
                    continue;
                }

                // A self-loop synchronization reads the passive node's state
                // AFTER the active half has been applied, since they are the
                // same node; otherwise the two halves see independent blocks.
                const std::vector<T>& src =
                    node_p == node_a ? oa.space[ia] : st.local[isf_p - 1];
                const EventOutcome<T> op =
                    qn::after_event(sn, node_p, src, sy.passive.event, sy.passive.cls);
                // The routing probability, read at the state the job LEAVES from,
                // which is what `sub_sdr` reads: the branch it is admitted to is
                // decided by the populations the departing customer sees.
                //
                // ROUND-ROBIN IS THE ONE THAT READS THE STATE AFTER. Its
                // destination is the pointer the ACTIVE node carries once its own
                // departure has advanced it, so the probability is the 0/1
                // indicator of that pointer and not the uniform mask
                // `refresh_routing` wrote into `rt`. Reading `rt` here instead
                // would spread the job over every outlink, which is the random
                // routing the dispatcher exists not to be. This is the
                // reference's `sub_rr`, which likewise takes `state_after`.
                T proute = sy.passive.statedep
                               ? rt_by_state[s](sy.passive.rt_row, sy.passive.rt_col)
                               : sy.passive.prob;
                if (rr_here) {
                    const std::size_t w = sn.nvars_of(node_a);
                    const std::vector<T>& arow = oa.space[ia];
                    std::size_t dest = 0;
                    if (arow.size() >= w) {
                        const std::vector<T> var(arow.end() - w, arow.end());
                        dest = sn.rr_dest(node_a, sy.active.cls, var);
                    }
                    proute = (dest == node_p && sy.passive.cls == sy.active.cls)
                                 ? num_traits<T>::from_int(1)
                                 : zero;
                }
                bool placed = false;
                for (std::size_t ip = 0; ip < op.space.size(); ++ip) {
                    NetState<T> nsx = st;
                    nsx.local[isf_a - 1] = oa.space[ia];
                    nsx.local[isf_p - 1] = op.space[ip];
                    const typename std::map<std::vector<double>, std::size_t>::const_iterator it =
                        index.find(ctmc_detail::state_key(nsx));
                    if (it == index.end()) continue;
                    placed = true;
                    // THE ACTIVE HALF'S OWN PROBABILITY COUNTS TOO. `oa.prob[ia]`
                    // is the share of the completion that leads to THIS successor
                    // -- SIRO's random pick of the next job to promote is the only
                    // branch that returns it below 1 -- and the LOCAL branch above
                    // already multiplies by it. Omitting it here gave every
                    // promotion candidate the FULL service rate, so a SIRO station
                    // with two waiting classes left the state at 2*mu: on
                    // prio_hol_open that is Source Tput 0.052334 against the
                    // reference's 0.052281. The reference folds the same share
                    // into its `outrate` instead (afterEventStation.m, `pick_prob`)
                    // and keeps `outprob` at 1, which is the same product.
                    const T w = T(rate * oa.prob[ia] * proute * op.prob[ip]);
                    res.Q(s, it->second) += w;
                    if (imm_action[a]) res.Qimm(s, it->second) += w;
                    if (want_filtration) {
                        res.filt[a](s, it->second) += w;
                        // Both halves of the synchronization are tagged: a DEP
                        // promotes at the sender while the paired ARV starts or
                        // preempts at the receiver.
                        ctmc_detail::add_aux_filt(res, sn, node_a, s, it->second, w, oa, ia);
                        ctmc_detail::add_aux_filt(res, sn, node_p, s, it->second, w, op, ip);
                    }
                    fired += w;
                }
                // TRUE BAS, the become-blocked half. The passive arrival was
                // refused at every outcome, so the completing job cannot leave.
                // Only the generator can emit this edge: the event layer sees one
                // node at a time and cannot know the destination is full.
                //
                // The successor is the CURRENT state with the marker set, not the
                // post-departure one: the job stays in the server it completed in,
                // which is the whole content of blocking after service. It is
                // accumulated separately from Q and folded in below because it is
                // NOT a departure -- counting it in `dep_rates` would inflate
                // throughput by the blocked transitions.
                if (!placed && sy.active.event == EventType::DEP &&
                    node_a <= sn.isbasblocking.size() && sn.isbasblocking[node_a - 1] &&
                    !st.local[isf_a - 1].empty() &&
                    num_traits<T>::to_double(st.local[isf_a - 1].back()) == 0) {
                    NetState<T> nsb = st;
                    nsb.local[isf_a - 1].back() = num_traits<T>::from_int(1);
                    const typename std::map<std::vector<double>, std::size_t>::const_iterator ib =
                        index.find(ctmc_detail::state_key(nsb));
                    // THE ROUTING PROBABILITY IS NOT APPLIED, verbatim from
                    // `solver_ctmc.m:361`, which adds `rate_a(ia)` alone. It makes
                    // no difference where the blocked destination is the only one,
                    // and `solver_ctmc_avg_from_pi` declines to shift queue lengths
                    // at a station with several destinations anyway.
                    if (ib != index.end()) bas_block(s, ib->second) += rate;
                }
            }
            // A DEP synchronization is one job LEAVING the active node and
            // ENTERING the passive one, so the same accumulated rate is both a
            // departure there and an arrival here. The passive half of a LOCAL
            // action is the dummy node, which is nobody's arrival.
            if (sy.active.event == EventType::DEP && num_traits<T>::to_double(fired) != 0) {
                res.dep_rates[s][isf_a - 1][sy.active.cls - 1] += fired;
                if (isf_p != 0) res.arv_rates[s][isf_p - 1][sy.passive.cls - 1] += fired;
                // ONLY AN FJ-AUGMENTED JOIN takes the rate complement, verbatim
                // from `solver_ctmc.m:841`: a Router or Fork DEP is restricted
                // to the tangible rows like any timed action, and only a Join
                // firing -- which exists nowhere else -- is complemented back.
                if (isfjaug && sn.nodes[node_a - 1].nodetype == NodeType::Join) {
                    res.dep_rates_imm[s][isf_a - 1][sy.active.cls - 1] += fired;
                    if (isf_p != 0)
                        res.arv_rates_imm[s][isf_p - 1][sy.passive.cls - 1] += fired;
                }
            }
        }
    }

    // SPN global synchronizations. A firing is ATOMIC across all its arcs, so
    // unlike an ordinary sync it rewrites several nodes in one transition and
    // cannot be decomposed into per-node halves.
    //
    // A PLACE'S FLOW IS COUNTED HERE OR NOWHERE. `refresh_sync` emits no DEP
    // sync touching a Place -- a token crosses an arc of a firing, not a routing
    // edge -- so leaving this loop to write only Q left `dep_rates` and
    // `arv_rates` identically zero at every Place, and the AvgTable reported
    // Tput 0 (hence RespT 0) for a Place whose QLen and Util were right. The
    // reference accumulates the same two from `Dfilt_gsync_comp`
    // (solver_ctmc.m:624-649): the PRE passives of a FIRE are that place's
    // departures and the POST passives its arrivals.
    //
    // ONLY A COMPLETION COUNTS, which is what `GlobalOutcome::completion`
    // records and what the reference's `is_comp` gates `Dfilt_gsync_comp` on. A
    // FIRE outcome that merely starts a firing phase moves no token, and an
    // ENABLE outcome never does, so neither is a flow.
    for (std::size_t g = 0; g < gsync.size(); ++g) {
        const bool is_fire = gsync[g].active.event == EventType::FIRE;
        for (std::size_t s = 0; s < n; ++s) {
            const qn::GlobalOutcome<T> go = qn::after_global_event(sn, space[s], gsync[g]);
            T completed = zero;
            for (std::size_t io = 0; io < go.space.size(); ++io) {
                const T contrib = T(go.rate[io] * go.prob[io]);
                if (num_traits<T>::to_double(contrib) == 0) continue;
                const typename std::map<std::vector<double>, std::size_t>::const_iterator it =
                    index.find(ctmc_detail::state_key(go.space[io]));
                if (it == index.end()) continue;
                res.Q(s, it->second) += contrib;
                if (imm_gsync[g]) res.Qimm(s, it->second) += contrib;
                if (is_fire && io < go.completion.size() && go.completion[io]) completed += contrib;
            }
            if (!is_fire || num_traits<T>::to_double(completed) == 0) continue;
            for (std::size_t j = 0; j < gsync[g].passive.size(); ++j) {
                const qn::ModeEvent<T>& pev = gsync[g].passive[j];
                if (pev.node == 0 || pev.node > sn.nodes.size()) continue;
                const std::size_t isf_v = sn.stateful_index(pev.node);
                if (isf_v == 0 || pev.cls == 0 || pev.cls > sn.nclasses) continue;
                // THE ARC MULTIPLICITY IS PART OF THE FLOW. A place loses (or
                // gains) `weight` tokens per firing, not one, so a rate counted
                // per firing is a firing rate and not a job rate. Measured
                // against the reference: `spn_closed_fourplaces`, whose cycle moves 2
                // tokens a firing and the unweighted count reported exactly half
                // its throughput, while `spn_twomodes`, whose two arcs have
                // multiplicities 4 and 2, was off by exactly those two factors
                // at its two places. Multiplicity 1 is the common case and
                // leaves `spn_basic_closed`/`_open` unchanged.
                const T flow = T(completed * pev.weight);
                // EVERY FIRE completion takes the rate complement, not only an
                // immediate one: `solver_ctmc.m:833-855` complements
                // `Dfilt_gsync_comp{g}` for all of them, and for a timed mode the
                // vanishing rows contribute nothing so the two agree.
                if (pev.event == EventType::PRE) {
                    res.dep_rates[s][isf_v - 1][pev.cls - 1] += flow;
                    if (has_imm) res.dep_rates_imm[s][isf_v - 1][pev.cls - 1] += flow;
                } else if (pev.event == EventType::POST) {
                    res.arv_rates[s][isf_v - 1][pev.cls - 1] += flow;
                    if (has_imm) res.arv_rates_imm[s][isf_v - 1][pev.cls - 1] += flow;
                }
            }
        }
    }

    // FORK FIRINGS. Like an SPN firing this is atomic across several nodes, so it
    // takes the whole network state and cannot be decomposed into sync halves.
    //
    // The rate statistics are accumulated by hand rather than through the DEP
    // path: a firing is one DEPARTURE of the parent class at the fork and one
    // ARRIVAL of the tag's auxiliary class at each branch head, which is B
    // arrivals for one departure and is exactly what an ordinary sync cannot
    // express. `refresh_sync` therefore emits no DEP sync for a Fork.
    for (std::size_t k = 0; k < fjsync.size(); ++k) {
        const qn::FjSync<T>& e = fjsync[k];
        const std::size_t isf_f = sn.stateful_index(e.fork);
        if (isf_f == 0) continue;
        for (std::size_t s = 0; s < n; ++s) {
            const qn::GlobalOutcome<T> fo = qn::after_fj_event(sn, e, space[s]);
            T fired = zero;
            for (std::size_t io = 0; io < fo.space.size(); ++io) {
                const T contrib = T(fo.rate[io] * fo.prob[io]);
                if (num_traits<T>::to_double(contrib) <= 0) continue;
                const typename std::map<std::vector<double>, std::size_t>::const_iterator it =
                    index.find(ctmc_detail::state_key(fo.space[io]));
                if (it == index.end()) continue;
                res.Q(s, it->second) += contrib;
                res.Qimm(s, it->second) += contrib;
                fired += contrib;
            }
            if (num_traits<T>::to_double(fired) == 0) continue;
            res.dep_rates[s][isf_f - 1][e.cls - 1] += fired;
            res.dep_rates_imm[s][isf_f - 1][e.cls - 1] += fired;
            for (std::size_t b = 0; b < e.branchheads.size(); ++b) {
                const std::size_t isf_b = sn.stateful_index(e.branchheads[b]);
                if (isf_b == 0) continue;
                res.arv_rates[s][isf_b - 1][e.auxclasses[b] - 1] += fired;
                res.arv_rates_imm[s][isf_b - 1][e.auxclasses[b] - 1] += fired;
            }
        }
    }

    if (any_bas)
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) res.Q(i, j) += bas_block(i, j);

    // THE VANISHING-ROW PURGE, `solver_ctmc.m:592-620`. A zero-sojourn state is
    // left the instant it is entered, so a TIMED arc out of one describes an
    // event that cannot occur -- the state is gone before the clock advances.
    // Restating such a row as its immediate part is what makes the stochastic
    // complement below the elimination of the immediate transitions rather than
    // a censoring of a chain that still contains them.
    if (has_imm) {
        res.vanishing = ctmc_find_vanishing_states(sn, space, gsync, isfjaug);
        // A vanishing state with NO immediate exit is a disagreement between the
        // predicate and the arc tagging, not a model: purging it would leave an
        // absorbing row and complementing it out would make the censored block
        // singular. The reference purges neither but complements it anyway and
        // lands on a NaN; drop it from both and say so.
        std::vector<std::size_t> keep_van;
        std::size_t gap = 0;
        for (std::size_t x = 0; x < res.vanishing.size(); ++x) {
            const std::size_t s = res.vanishing[x];
            T out_imm = zero;
            for (std::size_t j = 0; j < n; ++j)
                if (j != s) out_imm += res.Qimm(s, j);
            if (num_traits<T>::to_double(out_imm) > 0)
                keep_van.push_back(s);
            else
                ++gap;
        }
        if (gap > 0)
            line::util::LineConsole::step(
                "CTMC: %zu vanishing state(s) have no immediate outgoing arc; the vanishing "
                "predicate and the immediate-arc tagging disagree, so those rows keep their "
                "timed arcs",
                gap);
        res.vanishing.swap(keep_van);
        for (std::size_t x = 0; x < res.vanishing.size(); ++x) {
            const std::size_t s = res.vanishing[x];
            for (std::size_t j = 0; j < n; ++j) res.Q(s, j) = res.Qimm(s, j);
            res.arv_rates[s] = res.arv_rates_imm[s];
            res.dep_rates[s] = res.dep_rates_imm[s];
            for (std::size_t a = 0; a < res.filt.size(); ++a)
                if (!imm_action[a])
                    for (std::size_t j = 0; j < n; ++j) res.filt[a](s, j) = zero;
        }
    }

    make_infgen(res.Q);
    return res;
}

/**
 * Port of the "now remove immediate transitions" block of `solver_ctmc.m`
 * (:812-870): eliminate the vanishing states by stochastic complementation.
 *
 * WHAT IT CHANGES AND WHAT IT MUST NOT. The chain shrinks to its TANGIBLE states
 * and every metric read from it is unchanged to the digits printed, because the
 * vanishing states carry ~1e-8 of the mass -- which is exactly why the omission
 * was invisible until a caller asked for the chain itself. `-a states` and
 * `-a gen` are the callers that see it: on `fj_tiny_closed` this takes the six
 * enumerated states to the four the other three codebases return.
 *
 * THE RATE COMPLEMENT IS NOT OPTIONAL. An action that fires ONLY from vanishing
 * states -- a fork firing, a join rendezvous, an immediate SPN mode -- has no
 * tangible row to be read off, so restricting `arv_rates` / `dep_rates` to the
 * tangible states alone would silently zero its flow. The rate observed from a
 * tangible state is its own plus the expected number of firings along the
 * vanishing excursion entered from it,
 *
 *   r = total(nonimm) + Q12 (-Q22)^-1 imm_part(imm),
 *
 * which is `solver_ctmc_ratecomplement` applied to the immediate part and added
 * to the plain restriction of the total. The two agree term by term with the
 * reference's per-filtration form, since the totals are already the row sums.
 *
 * A no-op when the model declared no immediate source, which is the common case.
 */
template <class T>
void ctmc_eliminate_vanishing(CtmcResult<T>& res) {
    const T zero = num_traits<T>::from_int(0);
    if (res.vanishing.empty()) {
        res.Qimm = Matrix<T>();
        res.arv_rates_imm.clear();
        res.dep_rates_imm.clear();
        return;
    }
    const std::size_t n = res.Q.rows();
    std::vector<bool> is_van(n, false);
    for (std::size_t x = 0; x < res.vanishing.size(); ++x) is_van[res.vanishing[x]] = true;
    std::vector<std::size_t> nonimm;
    nonimm.reserve(n - res.vanishing.size());
    for (std::size_t s = 0; s < n; ++s)
        if (!is_van[s]) nonimm.push_back(s);
    if (nonimm.empty())
        throw NumericError(
            "SolverCTMC: every state is vanishing; the chain has no tangible state to observe");

    const mc::StochCompResult<T> sc = mc::ctmc_stochcomp(res.Q, nonimm);
    const std::size_t nk = nonimm.size(), nd = res.vanishing.size();

    // The rate complement, batched: one elimination of (-Q22) serves every
    // (stateful, class) pair of both directions, which is what keeps the cost at
    // one factorization rather than 2*NF*R of them.
    const std::size_t NF = res.arv_rates.empty() ? 0 : res.arv_rates[0].size();
    const std::size_t R = NF == 0 ? 0 : res.arv_rates[0][0].size();
    const std::size_t ncol = 2 * NF * R;
    Matrix<T> corr;
    if (ncol > 0) {
        Matrix<T> B(nd, ncol, zero);
        for (std::size_t d = 0; d < nd; ++d) {
            const std::size_t s = res.vanishing[d];
            for (std::size_t f = 0; f < NF; ++f)
                for (std::size_t r = 0; r < R; ++r) {
                    B(d, f * R + r) = res.arv_rates_imm[s][f][r];
                    B(d, NF * R + f * R + r) = res.dep_rates_imm[s][f][r];
                }
        }
        const Matrix<T> X = ctmc_detail::censored_solve(sc.Q22, B);
        corr = Matrix<T>(nk, ncol, zero);
        for (std::size_t a = 0; a < nk; ++a)
            for (std::size_t c = 0; c < ncol; ++c) {
                T acc = zero;
                for (std::size_t d = 0; d < nd; ++d) acc = T(acc + sc.Q12(a, d) * X(d, c));
                corr(a, c) = acc;
            }
    }

    std::vector<NetState<T>> space;
    std::vector<std::vector<std::vector<T>>> arv, dep;
    space.reserve(nk);
    arv.reserve(nk);
    dep.reserve(nk);
    for (std::size_t a = 0; a < nk; ++a) {
        space.push_back(res.space[nonimm[a]]);
        arv.push_back(res.arv_rates[nonimm[a]]);
        dep.push_back(res.dep_rates[nonimm[a]]);
        for (std::size_t f = 0; f < NF; ++f)
            for (std::size_t r = 0; r < R; ++r) {
                arv[a][f][r] = T(arv[a][f][r] + corr(a, f * R + r));
                dep[a][f][r] = T(dep[a][f][r] + corr(a, NF * R + f * R + r));
            }
    }

    // Each filtration is complemented exactly as Q is: an event that reaches a
    // tangible state THROUGH a vanishing excursion belongs on the arc it induces,
    // and dropping the excursion rows would lose it from every CDF built on the
    // split. The derived START and PREEMPT filtrations take the same treatment --
    // a service start that lands on a vanishing state would otherwise undercount.
    ctmc_detail::complement_filtrations(res, nonimm, res.vanishing, sc);

    res.Q = sc.S;
    res.space.swap(space);
    res.arv_rates.swap(arv);
    res.dep_rates.swap(dep);
    res.Qimm = Matrix<T>();
    res.arv_rates_imm.clear();
    res.dep_rates_imm.clear();
    res.vanishing.clear();
}


/**
 * Port of `State.reachableSpaceGenerator`: the states reachable from `init`.
 *
 * `space_generator` enumerates every state the ENCODING admits; this walks the
 * ones the DYNAMICS can actually occupy. The two differ whenever the encoding
 * is wider than the model -- a retrial station's idle-server states are
 * reachable, whereas an ordinary queue's are not, and enumerating the latter
 * leaves a generator with absorbing junk that perturbs the stationary vector
 * after normalization.
 *
 * The walk applies exactly the same handlers the generator does, so a state is
 * included precisely when some synchronization produces it at a positive rate.
 *
 * IT IS THE ONLY GENERATOR AN SPN HAS. `from_marginal_node` emits a single row
 * for a Transition -- every mode's servers free, nothing firing -- because a
 * transition's state is per-MODE and no population marginal determines it. The
 * lattice enumeration therefore never produces a state in which a mode is
 * firing, and a generator built over that space has every ENABLE landing
 * outside it. Walking `gsync` from the idle state is what materializes them,
 * which is why the reference forces `state_space_gen='reachable'` for any model
 * whose firings break per-chain population conservation.
 *
 * `cutoff` TRUNCATES AN OPEN CLASS, and without it this walk does not terminate.
 * The lattice generator bounds an open class's total population by the cutoff;
 * this walk had no such bound, so on ANY open SPN -- `spn_basic_open` at cutoff
 * 1, `spn_pareto_service`, `spn_open_sevenplaces` -- the Source kept producing
 * tokens and the walk ran to the `maxst` cap instead of answering. Passing the
 * same cutoff makes the two paths mean the same thing by "cutoff": a candidate
 * whose open-class population would exceed it is not a state of the truncated
 * chain, so the arc to it simply does not exist and `make_infgen` re-closes the
 * row, which is exactly what the lattice path leaves behind. Empty means
 * unbounded, which is right for a closed model and is what every existing
 * caller passes.
 *
 * THE INITIAL MARKING RAISES THE BOUND WHERE IT EXCEEDS IT. A Place may start
 * with more tokens than the cutoff -- `spn_open_sevenplaces` puts 2 in P1
 * against a default cutoff of 2 -- and a bound below the state the walk starts
 * from censors every successor of it, leaving the initial state alone in a
 * chain that is not the model's. A state space that cannot contain its own
 * initial state is empty by construction, so the floor is the initial marking.
 */
template <class T>
std::vector<NetState<T>> reachable_space_generator(
    const NetworkStruct<T>& sn, const NetState<T>& init, const std::vector<Sync<T>>& sync,
    const std::vector<qn::GlobalSync<T>>& gsync = std::vector<qn::GlobalSync<T>>(),
    std::size_t maxst = 3000000,
    const std::vector<qn::FjSync<T>>& fjsync = std::vector<qn::FjSync<T>>(),
    const std::vector<std::size_t>& cutoff = std::vector<std::size_t>(),
    const std::vector<std::vector<std::size_t>>& cutoff_mat =
        std::vector<std::vector<std::size_t>>()) {
    const std::size_t local = sn.nodes.size() + 1;
    std::vector<NetState<T>> out;
    std::map<std::vector<double>, std::size_t> seen;
    std::vector<std::size_t> stack;
    const std::vector<double> njobs = sn.njobs();
    // Only an OPEN class can leave the bound, so a closed model pays nothing.
    bool bound = false;
    for (std::size_t r = 0; r < sn.nclasses && r < njobs.size(); ++r)
        if (!std::isfinite(njobs[r]) && r < cutoff.size() && cutoff[r] > 0) bound = true;
    const bool bounded = bound;
    std::vector<std::size_t> lim = cutoff;
    if (bounded) {
        const std::vector<double> n0 = ctmc_detail::state_peak_occupancy(sn, init);
        for (std::size_t r = 0; r < lim.size() && r < n0.size(); ++r) {
            const std::size_t p0 =
                n0[r] > 0 ? static_cast<std::size_t>(std::floor(n0[r] + 0.5)) : 0;
            if (p0 > lim[r]) lim[r] = p0;
        }
    }

    seen[ctmc_detail::state_key(init)] = 0;
    out.push_back(init);
    stack.push_back(0);

    while (!stack.empty()) {
        const std::size_t si = stack.back();
        stack.pop_back();
        const NetState<T> st = out[si];  // by value: `out` grows inside the loop

        for (std::size_t a = 0; a < sync.size(); ++a) {
            const Sync<T>& sy = sync[a];
            const std::size_t isf_a = sn.stateful_index(sy.active.node);
            if (isf_a == 0) continue;
            const std::size_t isf_p =
                sy.passive.node == local ? 0 : sn.stateful_index(sy.passive.node);
            if (sy.passive.node != local && isf_p == 0) continue;

            const EventOutcome<T> oa = qn::after_event(sn, sy.active.node, st.local[isf_a - 1],
                                                       sy.active.event, sy.active.cls);
            for (std::size_t ia = 0; ia < oa.space.size(); ++ia) {
                if (num_traits<T>::to_double(oa.rate[ia]) <= 0) continue;

                std::vector<NetState<T>> cand;
                if (sy.passive.node == local) {
                    NetState<T> nsx = st;
                    nsx.local[isf_a - 1] = oa.space[ia];
                    cand.push_back(nsx);
                } else {
                    const std::vector<T>& src =
                        sy.passive.node == sy.active.node ? oa.space[ia] : st.local[isf_p - 1];
                    const EventOutcome<T> op = qn::after_event(sn, sy.passive.node, src,
                                                               sy.passive.event, sy.passive.cls);
                    for (std::size_t ip = 0; ip < op.space.size(); ++ip) {
                        if (num_traits<T>::to_double(op.prob[ip]) <= 0) continue;
                        NetState<T> nsx = st;
                        nsx.local[isf_a - 1] = oa.space[ia];
                        nsx.local[isf_p - 1] = op.space[ip];
                        cand.push_back(nsx);
                    }
                }
                for (std::size_t c = 0; c < cand.size(); ++c) {
                    if (bounded && !ctmc_detail::within_cutoff(sn, cand[c], njobs, lim, cutoff_mat))
                        continue;
                    const std::vector<double> key = ctmc_detail::state_key(cand[c]);
                    if (seen.find(key) != seen.end()) continue;
                    if (out.size() >= maxst)
                        throw UnsupportedError(
                            "reachable_space_generator: the reachable state space exceeds the "
                            "cap of " + std::to_string(maxst) + " states");
                    seen[key] = out.size();
                    out.push_back(cand[c]);
                    stack.push_back(out.size() - 1);
                }
            }
        }

        // The SPN half of the walk. A global synchronization already returns
        // WHOLE network states, since a firing is atomic across every arc it
        // touches and cannot be decomposed into an active and a passive half.
        for (std::size_t g = 0; g < gsync.size(); ++g) {
            const qn::GlobalOutcome<T> go = qn::after_global_event(sn, st, gsync[g]);
            for (std::size_t io = 0; io < go.space.size(); ++io) {
                if (num_traits<T>::to_double(go.rate[io]) <= 0) continue;
                if (num_traits<T>::to_double(go.prob[io]) <= 0) continue;
                if (bounded && !ctmc_detail::within_cutoff(sn, go.space[io], njobs, lim, cutoff_mat))
                    continue;
                const std::vector<double> key = ctmc_detail::state_key(go.space[io]);
                if (seen.find(key) != seen.end()) continue;
                if (out.size() >= maxst)
                    throw UnsupportedError(
                        "reachable_space_generator: the reachable state space exceeds the cap of " +
                        std::to_string(maxst) + " states");
                seen[key] = out.size();
                out.push_back(go.space[io]);
                stack.push_back(out.size() - 1);
            }
        }

        // The fork-join half. A firing is atomic in the same sense and returns
        // whole network states too; it is walked here rather than folded into the
        // sync loop because it has no active/passive decomposition at all.
        for (std::size_t k = 0; k < fjsync.size(); ++k) {
            const qn::GlobalOutcome<T> fo = qn::after_fj_event(sn, fjsync[k], st);
            for (std::size_t io = 0; io < fo.space.size(); ++io) {
                if (num_traits<T>::to_double(fo.rate[io]) <= 0) continue;
                if (num_traits<T>::to_double(fo.prob[io]) <= 0) continue;
                if (bounded && !ctmc_detail::within_cutoff(sn, fo.space[io], njobs, lim, cutoff_mat))
                    continue;
                const std::vector<double> key = ctmc_detail::state_key(fo.space[io]);
                if (seen.find(key) != seen.end()) continue;
                if (out.size() >= maxst)
                    throw UnsupportedError(
                        "reachable_space_generator: the reachable state space exceeds the cap of " +
                        std::to_string(maxst) + " states");
                seen[key] = out.size();
                out.push_back(fo.space[io]);
                stack.push_back(out.size() - 1);
            }
        }
    }
    return out;
}

/**
 * Port of `StateSpaceAggr`: the per-(station, class) job counts of every state,
 * as an (nstates x nstations*nclasses) matrix in column block order
 * `(ist-1)*K + k`.
 *
 * It is what `@@SolverCTMC/getStateSpaceAggr` returns and what the transient
 * analyzer, the reward analyzer and the BAS shift all index; building it once
 * keeps the three from re-deriving the same marginal decode with three chances
 * to disagree about the buffer encoding.
 *
 * A SOURCE ROW IS ZERO, not Inf. `to_marginal` reports an infinite reservoir for
 * an EXT station, which describes the encoding rather than a queue length, and
 * an Inf here would propagate into every aggregate that sums this matrix.
 */
template <class T>
Matrix<T> ctmc_state_space_aggr(const NetworkStruct<T>& sn,
                                const std::vector<NetState<T>>& space) {
    const std::size_t M = sn.stations.size(), K = sn.nclasses;
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> A(space.size(), M * K, zero);
    for (std::size_t ist = 1; ist <= M; ++ist) {
        const std::size_t isf = sn.stateful_of_station(ist);
        const std::size_t ind = sn.node_of_station(ist);
        if (isf == 0) continue;
        if (sn.stations[ist - 1].nodetype == NodeType::Source) continue;
        std::vector<std::size_t> ph(K, 1), shift(K, 0);
        std::size_t w = 0;
        for (std::size_t k = 0; k < K; ++k) {
            ph[k] = sn.phasessz_of(ist, k + 1);
            shift[k] = w;
            w += ph[k];
        }
        const std::size_t nvar = sn.nvars_of(ind);
        for (std::size_t s = 0; s < space.size(); ++s) {
            const qn::Marginal<T> m =
                qn::to_marginal(sn, ist, space[s].local[isf - 1], ph, shift, nvar);
            for (std::size_t k = 0; k < K; ++k) A(s, (ist - 1) * K + k) = m.nir[k];
        }
    }
    return A;
}

/**
 * Port of `ctmc_signal_lossy`: classes a G-network signal can annihilate here.
 *
 * Such a job leaves the station WITHOUT a service completion, so the
 * arrival-based (offered-load) utilization estimator is invalid for it and only
 * the departure-based carried load is meaningful -- the same reasoning as the
 * finite-capacity `canDropClass`, reached by a different route.
 *
 * A signal class is active at this node when its stationary arrival rate there
 * is positive. A TARGETED signal removes only its target class; an untargeted
 * one is class-agnostic and removes any non-signal class, matching
 * `after_event_station_signal`, MAM and LDES.
 */
template <class T>
std::vector<bool> ctmc_signal_lossy(const NetworkStruct<T>& sn, const CtmcResult<T>& r,
                                    const std::vector<T>& p, std::size_t isf) {
    const std::size_t K = sn.nclasses;
    std::vector<bool> lossy(K, false);
    bool any = false;
    for (std::size_t k = 0; k < K && k < sn.issignal.size(); ++k) any = any || sn.issignal[k];
    if (!any) return lossy;

    for (std::size_t r2 = 1; r2 <= K; ++r2) {
        if (sn.issignal.size() < r2 || !sn.issignal[r2 - 1]) continue;
        T arv = num_traits<T>::from_int(0);
        for (std::size_t s = 0; s < r.space.size(); ++s)
            arv += T(p[s] * r.arv_rates[s][isf - 1][r2 - 1]);
        if (num_traits<T>::to_double(arv) <= 0) continue;
        const std::size_t tgt =
            sn.signaltarget.size() >= r2 ? sn.signaltarget[r2 - 1] : 0;
        if (tgt >= 1 && tgt <= K) {
            lossy[tgt - 1] = true;
        } else {
            for (std::size_t k = 0; k < K; ++k)
                if (k >= sn.issignal.size() || !sn.issignal[k]) lossy[k] = true;
        }
    }
    return lossy;
}

/**
 * Port of `ctmc_signal_busy`: the exact per-class busy-server fraction, read
 * off the enumerated state space.
 *
 * WHY THE DEPARTURE ESTIMATOR IS NOT ENOUGH. `T*E[S]/c` is exact for a lossy
 * class only under EXPONENTIAL service: a job destroyed mid-service leaves busy
 * time behind with no completion to account for it, so with phase-type service
 * the carried-rate estimator under-counts. Measured on an M/Er2/1 with
 * lambda+ = 0.5 and lambda- = 0.4 it gives 0.34941 against a true 0.37696. The
 * in-service occupancy below is exact for any service process.
 *
 * A PS-like discipline shares the servers among every resident job, so class k
 * takes the weighted share n_k w_k / sum_j n_j w_j of the busy servers; every
 * other discipline exposes the in-service indicator directly through
 * `to_marginal`.
 */
template <class T>
std::vector<T> ctmc_signal_busy(const NetworkStruct<T>& sn, std::size_t ist,
                                const CtmcResult<T>& r, const std::vector<T>& p,
                                std::size_t isf) {
    const std::size_t K = sn.nclasses;
    const T zero = num_traits<T>::from_int(0);
    std::vector<T> unb(K, zero);
    const SchedStrategy sched = sn.stations[ist - 1].sched;
    const bool is_ps = sched == SchedStrategy::PS || sched == SchedStrategy::DPS ||
                       sched == SchedStrategy::GPS || sched == SchedStrategy::LPS;
    const double S = sn.stations[ist - 1].nservers;
    const std::size_t ind = sn.node_of_station(ist);

    std::vector<std::size_t> ph(K, 1), shift(K, 0);
    std::size_t w = 0;
    for (std::size_t k = 0; k < K; ++k) {
        ph[k] = sn.phasessz_of(ist, k + 1);
        shift[k] = w;
        w += ph[k];
    }
    const std::size_t nvar = sn.nvars_of(ind);

    for (std::size_t s = 0; s < r.space.size(); ++s) {
        if (num_traits<T>::to_double(p[s]) == 0) continue;
        const qn::Marginal<T> m =
            qn::to_marginal(sn, ist, r.space[s].local[isf - 1], ph, shift, nvar);
        double ni = 0;
        for (std::size_t k = 0; k < K; ++k) ni += num_traits<T>::to_double(m.nir[k]);
        if (ni <= 0) continue;
        if (is_ps) {
            T wtot = zero;
            for (std::size_t k = 0; k < K; ++k)
                wtot += T(m.nir[k] * sn.stations[ist - 1].schedparam[k]);
            if (num_traits<T>::to_double(wtot) <= 0) continue;
            const double busy = std::min(ni, S) / S;
            for (std::size_t k = 0; k < K; ++k)
                unb[k] += T(p[s] * (m.nir[k] * sn.stations[ist - 1].schedparam[k] / wtot) *
                            num_traits<T>::from_double(busy));
        } else {
            for (std::size_t k = 0; k < K; ++k)
                unb[k] += T(p[s] * m.sir[k] / num_traits<T>::from_double(S));
        }
    }
    return unb;
}

/**
 * MATLAB's `all(sn.isph(:))`, read off the matrices instead of off a flag.
 *
 * False as soon as one station-class process is a matrix exponential: its D0
 * carries negative off-diagonal entries, or its D1 negative entries, neither of
 * which a phase-type has. A Source keeps its arrival process in `sn.service`,
 * so this one scan covers arrivals and services alike.
 */
template <class T>
bool ctmc_all_phasetype(const NetworkStruct<T>& sn) {
    for (std::size_t i = 0; i < sn.service.size(); ++i)
        for (std::size_t r = 0; r < sn.service[i].size(); ++r) {
            const lang::Distrib<T>& d = sn.service[i][r];
            if (d.disabled || d.D0.rows() == 0) continue;
            for (std::size_t a = 0; a < d.D0.rows(); ++a)
                for (std::size_t b = 0; b < d.D0.cols(); ++b)
                    if (a != b && num_traits<T>::to_double(d.D0(a, b)) < 0) return false;
            for (std::size_t a = 0; a < d.D1.rows(); ++a)
                for (std::size_t b = 0; b < d.D1.cols(); ++b)
                    if (num_traits<T>::to_double(d.D1(a, b)) < 0) return false;
        }
    return true;
}

/** The mean performance metrics a stationary vector maps to. */
template <class T>
struct CtmcAvg {
    Matrix<T> QN, UN, RN, TN;  ///< (nstations x nclasses)
    std::vector<T> XN, CN;     ///< (nclasses) system throughput and response time
    /**
     * The DERIVED rates, (nstations x nclasses): how often per unit time a
     * class-r service STARTS at station i, and how often a class-r job in
     * service is PUSHED BACK into the buffer there. All zero unless the result
     * carried the filtrations they reduce (`solver_ctmc` was asked for them).
     *
     * At a lossless station with no in-service abandonment
     *     StartN == TN + PreemptN,
     * because every job starts service once per entry into a server and every
     * preemption is followed by exactly one later resume or restart.
     */
    Matrix<T> StartN, PreemptN;
};

/**
 * Port of `solver_ctmc_avg_from_pi`: map a state distribution to mean metrics.
 *
 * Factored from the analyzer exactly as the reference factors it, so a caller
 * holding its own distribution -- a time-averaged transient one, say -- reuses
 * the same discipline-aware reduction instead of re-deriving it.
 */
template <class T>
CtmcAvg<T> solver_ctmc_avg_from_pi(const NetworkStruct<T>& sn, const CtmcResult<T>& r,
                                   const std::vector<T>& pivec) {
    const std::size_t M = sn.stations.size(), R = sn.nclasses, n = r.space.size();
    const T zero = num_traits<T>::from_int(0);
    CtmcAvg<T> a;
    a.QN = Matrix<T>(M, R, zero);
    a.UN = Matrix<T>(M, R, zero);
    a.RN = Matrix<T>(M, R, zero);
    a.TN = Matrix<T>(M, R, zero);
    a.XN.assign(R, zero);
    a.CN.assign(R, zero);
    a.StartN = Matrix<T>(M, R, zero);
    a.PreemptN = Matrix<T>(M, R, zero);

    // Renormalize, clamping the numerical dust an eigenvector solve leaves
    // below the zero threshold. WITH A MATRIX EXPONENTIAL THE CLAMP IS SKIPPED:
    // the stationary vector is then a genuinely SIGNED measure, only its
    // aggregates over each phase block are probabilities, and deleting the
    // negative entries deletes real mass -- an M/CME/1 at rho = 0.3 came out
    // with QLen 0.4469 against the Pollaczek-Khinchine 0.3772. Every metric
    // below is linear in the vector and stays exact without the clamp.
    const bool signed_measure = !ctmc_all_phasetype(sn);
    std::vector<T> p = pivec;
    T tot = zero;
    for (std::size_t s = 0; s < n; ++s) {
        if (!signed_measure && num_traits<T>::to_double(p[s]) < GlobalConstants::Zero) p[s] = zero;
        tot += p[s];
    }
    if (num_traits<T>::to_double(tot) > 0)
        for (std::size_t s = 0; s < n; ++s) p[s] = T(p[s] / tot);

    // System throughput is the arrival rate seen at each class's REFERENCE
    // station, which is what makes X a per-class quantity rather than a sum.
    for (std::size_t k = 1; k <= R; ++k) {
        const std::size_t refsf = sn.stateful_of_station(sn.classes[k - 1].refstat);
        for (std::size_t s = 0; s < n; ++s) a.XN[k - 1] += T(p[s] * r.arv_rates[s][refsf - 1][k - 1]);
    }

    // The derived rates: pi * F * e over each filtration, the same reduction
    // the departure rates use, so the three are directly comparable.
    for (std::size_t ist = 1; ist <= M && ist <= r.start_filt.size(); ++ist)
        for (std::size_t k = 1; k <= R && k <= r.start_filt[ist - 1].size(); ++k) {
            T accs = zero, accp = zero;
            for (std::size_t s = 0; s < n; ++s) {
                T rows = zero, rowp = zero;
                for (std::size_t ns = 0; ns < n; ++ns) {
                    rows += r.start_filt[ist - 1][k - 1](s, ns);
                    rowp += r.preempt_filt[ist - 1][k - 1](s, ns);
                }
                accs += T(p[s] * rows);
                accp += T(p[s] * rowp);
            }
            a.StartN(ist - 1, k - 1) = accs;
            a.PreemptN(ist - 1, k - 1) = accp;
        }

    for (std::size_t ist = 1; ist <= M; ++ist) {
        const std::size_t isf = sn.stateful_of_station(ist);
        const std::size_t ind = sn.node_of_station(ist);
        const bool is_source = sn.stations[ist - 1].nodetype == NodeType::Source;
        const double S = sn.stations[ist - 1].nservers;
        std::vector<std::size_t> ph(R, 1), shift(R, 0);
        std::size_t w = 0;
        for (std::size_t r2 = 0; r2 < R; ++r2) {
            ph[r2] = sn.phasessz_of(ist, r2 + 1);
            shift[r2] = w;
            w += ph[r2];
        }
        const std::size_t nvar = sn.nvars_of(ind);

        for (std::size_t k = 1; k <= R; ++k)
            for (std::size_t s = 0; s < n; ++s)
                a.TN(ist - 1, k - 1) += T(p[s] * r.dep_rates[s][isf - 1][k - 1]);

        if (is_source) {
            // `to_marginal` encodes a Source as nir = Inf, an infinite
            // reservoir. That is a statement about the ENCODING, not a queue
            // length; reading it as one gives Q = Inf and then R = Q/T = Inf.
            continue;
        }

        for (std::size_t s = 0; s < n; ++s) {
            if (num_traits<T>::to_double(p[s]) == 0) continue;
            const qn::Marginal<T> m =
                qn::to_marginal(sn, ist, r.space[s].local[isf - 1], ph, shift, nvar);
            for (std::size_t k = 1; k <= R; ++k) a.QN(ist - 1, k - 1) += T(p[s] * m.nir[k - 1]);
        }

        const SchedStrategy sched = sn.stations[ist - 1].sched;
        // PAS / order-independent: utilization is the IN-SERVICE occupancy, not
        // the offered load. The two coincide only when a job engages a single
        // server, which is exactly what a pass-and-swap station does not do.
        if (sched == SchedStrategy::PAS) {
            for (std::size_t s = 0; s < n; ++s) {
                if (num_traits<T>::to_double(p[s]) == 0) continue;
                const qn::Marginal<T> m =
                    qn::to_marginal(sn, ist, r.space[s].local[isf - 1], ph, shift, nvar);
                for (std::size_t k = 0; k < R; ++k)
                    a.UN(ist - 1, k) += T(p[s] * m.sir[k] / num_traits<T>::from_double(S));
            }
            continue;
        }
        // A load-dependent station has no single service rate, so `T*E[S]/c` is
        // not its utilization: the reference accumulates the per-state capacity
        // share weighted by the lld factor and divides by the EFFECTIVE server
        // count max(c, max lld), which is what the scaling can deliver.
        // A class- or joint-dependent station takes the SAME branch: none of the
        // three has a single service rate, so `T*E[S]/c` is not its utilization.
        // The reference gates all three on one `isempty(lld) && isempty(cd) &&
        // isempty(jd)` test, and the cd/jd cases are then OVERWRITTEN by the
        // declared-peak normalization at the end of this function -- this branch
        // is what a jd-only station keeps, and what a cd station holds until the
        // peak pass replaces it.
        const std::vector<T>& lld = sn.stations[ist - 1].lldscaling;
        if (!lld.empty() || sn.stations[ist - 1].cdscaling || sn.stations[ist - 1].jdscaling) {
            double ceff = S;
            for (std::size_t j = 0; j < lld.size(); ++j)
                ceff = std::max(ceff, num_traits<T>::to_double(lld[j]));
            const bool share = sched == SchedStrategy::PS || sched == SchedStrategy::DPS ||
                               sched == SchedStrategy::GPS || sched == SchedStrategy::LPS;
            for (std::size_t s = 0; s < n; ++s) {
                if (num_traits<T>::to_double(p[s]) == 0) continue;
                const qn::Marginal<T> m =
                    qn::to_marginal(sn, ist, r.space[s].local[isf - 1], ph, shift, nvar);
                double ni = 0;
                for (std::size_t k = 0; k < R; ++k) ni += num_traits<T>::to_double(m.nir[k]);
                if (ni <= 0) continue;
                double lldnow = 1.0;
                if (!lld.empty()) {
                    const std::size_t li = std::min<std::size_t>(
                        lld.size(), std::max<std::size_t>(1, static_cast<std::size_t>(ni)));
                    lldnow = num_traits<T>::to_double(lld[li - 1]);
                }
                if (share) {
                    T wtot = zero;
                    for (std::size_t k = 0; k < R; ++k)
                        wtot += T(m.nir[k] * sn.stations[ist - 1].schedparam[k]);
                    if (num_traits<T>::to_double(wtot) <= 0) continue;
                    for (std::size_t k = 0; k < R; ++k)
                        a.UN(ist - 1, k) +=
                            T(p[s] * (m.nir[k] * sn.stations[ist - 1].schedparam[k] / wtot) *
                              num_traits<T>::from_double(lldnow / ceff));
                } else {
                    double sirtot = 0;
                    for (std::size_t k = 0; k < R; ++k)
                        sirtot += num_traits<T>::to_double(m.sir[k]);
                    if (sirtot <= 0) continue;
                    for (std::size_t k = 0; k < R; ++k)
                        a.UN(ist - 1, k) += T(p[s] * num_traits<T>::from_double(
                                                         num_traits<T>::to_double(m.sir[k]) /
                                                         sirtot * lldnow / ceff));
                }
            }
            continue;
        }
        if (sched == SchedStrategy::INF) {
            // An infinite server is "utilized" by every job it holds: there is
            // no queueing, so utilization and queue length coincide.
            for (std::size_t k = 1; k <= R; ++k) a.UN(ist - 1, k - 1) = a.QN(ist - 1, k - 1);
        } else {
            // A class that can be DROPPED here -- an open class at a station
            // with a finite capacity -- must be measured on the CARRIED rate
            // alone. The offered rate counts arrivals that never entered
            // service, so `max` would report the offered load as utilization:
            // for an M/M/1/4 with lambda = 0.6 that is 0.6 against the true
            // 1 - p0 = 0.566. Where nothing can be dropped the two estimates
            // agree in steady state and the max only guards numerical noise.
            // A class a G-network SIGNAL can annihilate is droppable for the
            // same reason a capacity-limited one is: the job leaves without a
            // completion, so the offered rate is not what the server did.
            const std::vector<bool> lossy = ctmc_signal_lossy(sn, r, p, isf);
            // A station inside a DROP region loses arrivals it cannot admit,
            // exactly as a finite per-station capacity does, so the same
            // carried-rate rule applies there.
            bool in_drop = false;
            for (std::size_t f = 0; f < sn.regions.size() && !in_drop; ++f) {
                bool has_drop = false;
                for (std::size_t rr = 0; rr < sn.regions[f].rule.size(); ++rr)
                    if (sn.regions[f].rule[rr] == lang::DropStrategy::DROP) has_drop = true;
                if (has_drop && ist - 1 < sn.regions[f].members.size() &&
                    sn.regions[f].members[ist - 1])
                    in_drop = true;
            }
            for (std::size_t k = 1; k <= R; ++k) {
                const bool can_drop =
                    (!std::isfinite(sn.njobs()[k - 1]) &&
                     (std::isfinite(sn.cap[ist - 1]) ||
                      std::isfinite(sn.classcap[ist - 1][k - 1]) || in_drop)) ||
                    lossy[k - 1];
                const lang::Distrib<T>& d = sn.service[ist - 1][k - 1];
                if (d.disabled || d.D0.rows() == 0) continue;
                mam::Map<T> mp;
                mp.D0 = d.D0;
                mp.D1 = d.D1;
                const T mean = mam::map_mean(mp);
                const T u_dep = T(a.TN(ist - 1, k - 1) * mean / num_traits<T>::from_double(S));
                if (can_drop) {
                    a.UN(ist - 1, k - 1) = u_dep;
                    continue;
                }
                T arv = zero;
                for (std::size_t s = 0; s < n; ++s)
                    arv += T(p[s] * r.arv_rates[s][isf - 1][k - 1]);
                const T u_arv = T(arv * mean / num_traits<T>::from_double(S));
                a.UN(ist - 1, k - 1) =
                    num_traits<T>::to_double(u_arv) > num_traits<T>::to_double(u_dep) ? u_arv
                                                                                     : u_dep;
            }
            // For a lossy class the carried rate is still not exact unless the
            // service is exponential, so the in-service occupancy read off the
            // state space REPLACES it -- see `ctmc_signal_busy`.
            bool anylossy = false;
            for (std::size_t k = 0; k < R; ++k) anylossy = anylossy || lossy[k];
            if (anylossy) {
                const std::vector<T> unb = ctmc_signal_busy(sn, ist, r, p, isf);
                for (std::size_t k = 0; k < R; ++k)
                    if (lossy[k]) a.UN(ist - 1, k) = unb[k];
            }
        }
    }

    // TRUE BAS: the held job is counted at its DESTINATION, not where it sits.
    //
    // The state has it at the blocking station -- that is what the marker means --
    // but it has FINISHED service there and is on its way out, so reporting it in
    // the upstream queue length would double-count the time it spends waiting for
    // room. The reference moves it, and only when the destination is unambiguous:
    // with several downstream stations there is no single place to move it to, so
    // it is left where it sits rather than assigned arbitrarily.
    for (std::size_t ist = 1; ist <= M && !sn.isbasblocking.empty(); ++ist) {
        const std::size_t ind = sn.node_of_station(ist);
        const std::size_t isf = sn.stateful_of_station(ist);
        if (ind == 0 || isf == 0 || ind > sn.isbasblocking.size()) continue;
        if (!sn.isbasblocking[ind - 1]) continue;
        const std::vector<std::size_t> dests = sn.downstream_stations(ind);
        if (dests.size() != 1) continue;
        const std::size_t jst = sn.nodes[dests[0] - 1].station;
        if (jst == 0) continue;
        std::vector<std::size_t> ph2(R, 1), sh2(R, 0);
        std::size_t w2 = 0;
        for (std::size_t k = 0; k < R; ++k) {
            ph2[k] = sn.phasessz_of(ist, k + 1);
            sh2[k] = w2;
            w2 += ph2[k];
        }
        const std::size_t nv2 = sn.nvars_of(ind);
        for (std::size_t s = 0; s < n; ++s) {
            if (num_traits<T>::to_double(p[s]) == 0) continue;
            const std::vector<T>& row = r.space[s].local[isf - 1];
            if (row.empty() || num_traits<T>::to_double(row.back()) != 1) continue;
            const qn::Marginal<T> m = qn::to_marginal(sn, ist, row, ph2, sh2, nv2);
            for (std::size_t k = 0; k < R; ++k) {
                // A blocked state holds EXACTLY ONE completed job, so the count
                // is capped at 1: the rest of the queue has not finished and
                // stays where it is.
                const double nk = num_traits<T>::to_double(m.nir[k]);
                if (!(nk > 0)) continue;
                const T shift = T(p[s] * num_traits<T>::from_double(std::min(nk, 1.0)));
                a.QN(ist - 1, k) -= shift;
                a.QN(jst - 1, k) += shift;
            }
        }
    }

    // DECLARED-PEAK UTILIZATION at a class- or joint-dependent station. The
    // per-state capacity share accumulated above is a busy-server probability,
    // and at a dependent station that is not what utilization means: a beta_r(n)
    // emulating extra servers would report more than one. The reference REPLACES
    // it by T/mu/peak, restoring the T*S/c convention against the DECLARED peak.
    //
    // The `all njobs finite` guard is the reference's and is not a convenience:
    // with an open class the state space is a TRUNCATION, so the throughput this
    // divides is the truncated chain's and the ratio is not a utilization of the
    // model. The busy-server estimate above at least stays a probability, so it
    // is what an open dependent station keeps.
    //
    // THE TWO PEAKS MULTIPLY. beta_r(n) and eta_i(n) scale the SAME nominal rate
    // and the event layer folds them multiplicatively (`cd_factor`), so the peak
    // attainable rate is `rates * cdpeak * jdpeak`. The reference used to run two
    // independent blocks each ASSIGNING UN, which let the jd peak overwrite the cd
    // one and dropped a factor `cdpeak` at a station carrying both; that was fixed
    // in `solver_ctmc_analyzer.m` / `solver_ctmc_avg_from_pi.m` together with this
    // port, and it is what SolverSSA already did
    // (`solver_ssa_analyzer_serial.m:87-90`).
    bool all_closed = true;
    for (std::size_t k = 0; k < R; ++k)
        if (!std::isfinite(sn.njobs()[k])) all_closed = false;
    if (all_closed) {
        for (std::size_t ist = 1; ist <= M; ++ist) {
            const qn::Station<T>& st = sn.stations[ist - 1];
            const bool has_cd = static_cast<bool>(st.cdscaling);
            const bool has_jd = static_cast<bool>(st.jdscaling);
            if (!has_cd && !has_jd) continue;
            for (std::size_t k = 0; k < R; ++k) {
                double bmax = 1.0;
                if (has_cd)
                    bmax *= k < st.cdscalingpeak.size()
                                ? num_traits<T>::to_double(st.cdscalingpeak[k])
                                : 0.0;
                if (has_jd)
                    bmax *= k < st.jdscalingpeak.size()
                                ? num_traits<T>::to_double(st.jdscalingpeak[k])
                                : 0.0;
                const double mu = num_traits<T>::to_double(sn.rates(ist - 1, k));
                a.UN(ist - 1, k) = (std::isfinite(mu) && mu > 0 && bmax > 0)
                                       ? T(a.TN(ist - 1, k) /
                                           num_traits<T>::from_double(mu * bmax))
                                       : zero;
            }
        }
    }

    // SYNCHRONOUS CALLS. A caller blocked waiting for its reply still HOLDS its
    // server, so those jobs belong in QLen and Util even though they are not in
    // service here. They are NOT in the response time: response time is time
    // spent AT the station, and the blocked job is at its callee.
    Matrix<T> qn_blocked(M, R, zero);
    if (!sn.replyblock.empty()) {
        for (std::size_t ist = 1; ist <= M; ++ist) {
            const std::size_t ind = sn.node_of_station(ist);
            const std::size_t isf = sn.stateful_of_station(ist);
            if (isf == 0 || sn.replyblock.size() < ind) continue;
            bool any = false;
            for (std::size_t k = 0; k < sn.replyblock[ind - 1].size(); ++k)
                any = any || sn.replyblock[ind - 1][k];
            if (!any) continue;
            const qn::ReplyBlockInfo info = qn::reply_block_info(sn, ind);
            if (info.width == 0) continue;
            const double S2 = sn.stations[ist - 1].nservers;
            // The counters are the LAST columns of the local row, one per
            // calling class, in the order `reply_block_info` lists them.
            for (std::size_t s = 0; s < n; ++s) {
                if (num_traits<T>::to_double(p[s]) == 0) continue;
                const std::vector<T>& row = r.space[s].local[isf - 1];
                for (std::size_t pos = 0; pos < info.classes.size(); ++pos) {
                    const std::size_t col = row.size() - info.width + pos;
                    const std::size_t cls = info.classes[pos];
                    qn_blocked(ist - 1, cls - 1) += T(p[s] * row[col]);
                }
            }
            for (std::size_t k = 0; k < R; ++k) {
                a.QN(ist - 1, k) += qn_blocked(ist - 1, k);
                a.UN(ist - 1, k) += T(qn_blocked(ist - 1, k) / num_traits<T>::from_double(S2));
            }
        }
    }

    // Little's law per station, then per class at the reference station.
    for (std::size_t k = 1; k <= R; ++k) {
        for (std::size_t ist = 1; ist <= M; ++ist)
            a.RN(ist - 1, k - 1) =
                num_traits<T>::to_double(a.TN(ist - 1, k - 1)) > 0
                    ? T((a.QN(ist - 1, k - 1) - qn_blocked(ist - 1, k - 1)) /
                        a.TN(ist - 1, k - 1))
                    : zero;
        const double nk = sn.njobs()[k - 1];
        a.CN[k - 1] = std::isfinite(nk) && num_traits<T>::to_double(a.XN[k - 1]) > 0
                          ? T(num_traits<T>::from_double(nk) / a.XN[k - 1])
                          : zero;
    }
    return a;
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_H
