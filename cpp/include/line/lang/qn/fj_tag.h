/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_QN_FJ_TAG_H
#define LINE_LANG_QN_FJ_TAG_H

/**
 * Port of `matlab/src/api/fj/sn_fj_validate.m` and
 * `matlab/src/io/@@ModelAdapter/fjtag.m`: the model-to-model transform that makes
 * a closed fork-join network EXACTLY solvable by a Markov chain.
 *
 * WHY A TRANSFORM AND NOT A JOIN HANDLER. A join must release one parent job
 * once ALL siblings of THAT parent have arrived. A per-class count vector at the
 * join cannot express "of that parent": with two parents outstanding and two
 * branches, the counts (1,1) are consistent both with one sibling from each
 * parent (nothing may join) and with both siblings of one parent (a join must
 * fire). Any solver that reads only counts has to guess, and guessing is what
 * makes a join an approximation.
 *
 * WHAT THE TRANSFORM DOES. For each (fork f, class r) with matched join j, each
 * branch b = 1..B and each TAG t = 1..T, it mints an auxiliary closed class
 * A(f,r,b,t) with population 0. The tag names the parent job. A fork firing
 * consumes one parent held at the fork and emits one sibling per branch in the
 * classes of ONE tag -- the lowest free one -- and the join fires only when
 * every branch of some tag is present. Identity matching is then exact even when
 * siblings overtake one another, and the price is B*T extra classes per (fork,
 * class): the state space grows accordingly, which is why this is the exact
 * solver's route and `fj_mmt` remains the mean-value one.
 *
 * WHY T IS THE CHAIN POPULATION AND NOT THE CLASS POPULATION. Class switching
 * OUTSIDE the fork-join section can concentrate the whole chain in class r -- a
 * class switch on the edge into the fork does exactly that -- so the number of
 * concurrently outstanding forked jobs is bounded by the chain, not by the
 * declared population of r.
 *
 * WHY THE LOWEST FREE TAG. Tags are interchangeable, so without a canonical
 * choice each firing would produce T! equivalent successors and the chain would
 * carry a permutation group's worth of duplicate states. Allocating the lowest
 * free tag makes exactly one `fjsync` entry enabled per (fork, class) in any
 * state.
 *
 * HOW THIS DIFFERS FROM THE REFERENCE'S ROUTE. `fjtag.m` copies the Network
 * OBJECT, adds `ClosedClass` objects, re-links it and calls `getStruct`. Here
 * the transform is applied to the refreshed `NetworkStruct` directly, as
 * `tag_chain` and `fj_mmt` already are: the auxiliary routing is written into P
 * and one `refresh_struct()` re-derives the chains, capacities and rt. The
 * post-edits that follow are the reference's own -- it overrides the visits and
 * the auxiliary capacities after `getStruct` for the same reason, because the
 * engines read them only as a zero-versus-nonzero support gate.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"

namespace line {
namespace qn {

/** The augmented struct and everything needed to read its results back. */
template <class T>
struct FjTagged {
    NetworkStruct<T> V;
    /** fjclassmap[a-1] is the ORIGINAL class of auxiliary class a, 0 for originals. */
    std::vector<std::size_t> fjclassmap;
    std::vector<std::size_t> fjforkmap, fjjoinmap, fjbranchmap, fjtagmap;
    std::vector<FjSync<T>> fjsync;
    std::map<std::size_t, FjJoinParam> joinparam;  ///< keyed by 1-based Join node
    std::size_t korig = 0;
};

namespace fj_detail {

/** Total per-class node visits, `cellsum(sn.nodevisits)`. */
template <class T>
std::vector<std::vector<double>> node_visit_sum(const NetworkStruct<T>& sn) {
    const std::size_t I = sn.nodes.size(), K = sn.nclasses;
    std::vector<std::vector<double>> V(I, std::vector<double>(K, 0.0));
    for (std::size_t c = 0; c < sn.nodevisits.size(); ++c)
        for (std::size_t i = 0; i < I; ++i)
            for (std::size_t r = 0; r < K; ++r)
                V[i][r] += num_traits<T>::to_double(sn.nodevisits[c](i, r));
    return V;
}

}  // namespace fj_detail

/**
 * Port of `sn_fj_validate`: is this fork-join model inside the exact solver's
 * reach?
 *
 * The checks here are the STRUCTURAL ones that can be asked of the original
 * struct; the branch-level ones (class switching on a branch, nesting, sibling
 * traps) can only be asked during branch discovery and are raised there.
 */
template <class T>
void sn_fj_validate(const NetworkStruct<T>& sn) {
    const std::vector<std::vector<double>> Vn = fj_detail::node_visit_sum(sn);
    for (std::size_t f = 1; f <= sn.nodes.size(); ++f) {
        if (sn.nodes[f - 1].nodetype != NodeType::Fork) continue;
        std::size_t nj = 0;
        for (std::size_t p = 0; p < sn.fj.size(); ++p)
            if (sn.fj[p].first == f) ++nj;
        if (nj == 0)
            throw UnsupportedError(
                "sn_fj_validate: the Fork node '" + sn.nodes[f - 1].name +
                "' has no matched Join; the exact fork-join implementation needs the pair, since "
                "the tag pool is sized by the join that consumes it");
        if (nj > 1)
            throw UnsupportedError("sn_fj_validate: the Fork node '" + sn.nodes[f - 1].name +
                                   "' has several matched Joins, which the exact fork-join "
                                   "implementation does not support");
        // TASKS PER LINK. The tag-augmented construction carries an integer
        // weight per branch, so a link whose count is an integer is served
        // exactly, whether or not it differs from its siblings'. What it cannot
        // carry is a count that is not fixed at build time: a DISTRIBUTION has
        // to be drawn per firing, and the auxiliary class capacity would have to
        // be its support's maximum with a different tag occupancy per draw.
        const qn::NodeDef& fk = sn.nodes[f - 1];
        const qn::ForkParam<T>* fp = sn.fork_param_of(f);
        const double w = fk.tasks_per_link;
        if (fp == 0) {
            if (w != std::floor(w))
                throw UnsupportedError(
                    "sn_fj_validate: the Fork node '" + fk.name +
                    "' has a non-integer tasksPerLink; a sibling is a job and cannot be emitted "
                    "in fractions");
        } else {
            bool uncertain = false;
            for (std::size_t k = 0; k < fp->fan_out_link.rows(); ++k)
                for (std::size_t r = 0; r < fp->fan_out_link.cols(); ++r) {
                    const double p = num_traits<T>::to_double(fp->fan_out_prob(k, r));
                    if (p == 0.0) continue;  // link not taken
                    if (!fp->fan_out_dist[k][r].disabled)
                        throw UnsupportedError(
                            "sn_fj_validate: the Fork node '" + fk.name +
                            "' draws its tasks per link from a distribution, which the exact "
                            "fork-join implementation cannot carry: the tag occupancy would "
                            "differ per firing. Use SolverJMT or SolverLDES, which draw the "
                            "degree at the fork epoch");
                    const double lv = num_traits<T>::to_double(fp->fan_out_link(k, r));
                    if (lv != std::floor(lv))
                        throw UnsupportedError(
                            "sn_fj_validate: the Fork node '" + fk.name +
                            "' has a non-integer tasksPerLink; a sibling is a job and cannot be "
                            "emitted in fractions");
                    // A PER-LINK COUNT is carried: `ForkInfo::wlink` holds one
                    // integer per branch and it reaches the auxiliary capacity,
                    // the join's required count and the emission list. The only
                    // thing that has to be fixed at build time is that it is an
                    // integer and does not change per firing, which is what the
                    // two tests above ask.
                    if (p != 1.0) uncertain = true;
                }
            // BRANCH PROBABILITIES. A branch that may decline makes the SET of
            // siblings random, so a firing has 2^B outcomes and the join's
            // required count is a function of which subset fired. The tag
            // construction records no such per-firing state -- `required[r]` is
            // fixed when the state space is built -- so this path would emit
            // every branch anyway and answer with the CERTAIN-fork number. It is
            // refused by name instead: silently returning the answer to a
            // different model is the one outcome worth avoiding. MATLAB and
            // Python refuse it here for the same reason.
            if (uncertain)
                throw UnsupportedError(
                    "sn_fj_validate: the Fork node '" + fk.name +
                    "' has a branch activation probability below one, so the SET of siblings is "
                    "random and the tag construction fixes it when the state space is built. Use "
                    "SolverJMT or SolverLDES, which draw the activation at the fork epoch, or "
                    "SolverMVA, whose MMT transform sees the expected degree");
        }
        for (std::size_t r = 1; r <= sn.nclasses; ++r) {
            if (!(Vn[f - 1][r - 1] > 0)) continue;
            // The tag pool is the chain population, so the WHOLE chain must be
            // closed -- an open class anywhere in it makes the pool unbounded.
            std::size_t c = 0;
            for (std::size_t cc = 0; cc < sn.chains.size(); ++cc)
                if (sn.chains[cc][r - 1]) { c = cc + 1; break; }
            bool open = !std::isfinite(sn.njobs()[r - 1]);
            if (c != 0)
                for (std::size_t s = 0; s < sn.nclasses; ++s)
                    if (sn.chains[c - 1][s] && !std::isfinite(sn.njobs()[s])) open = true;
            if (open)
                throw UnsupportedError(
                    "sn_fj_validate: class '" + sn.classes[r - 1].name +
                    "' is routed through the Fork node '" + sn.nodes[f - 1].name +
                    "' and belongs to an OPEN chain. The tag pool is sized by the chain "
                    "population, which is unbounded here; use SolverMVA (fj_mmt) or SolverSSA");
        }
    }
    for (std::size_t j = 1; j <= sn.nodes.size(); ++j) {
        if (sn.nodes[j - 1].nodetype != NodeType::Join) continue;
        bool matched = false;
        std::size_t fsrc = 0;
        for (std::size_t p = 0; p < sn.fj.size(); ++p)
            if (sn.fj[p].second == j) { matched = true; fsrc = sn.fj[p].first; }
        if (!matched)
            throw UnsupportedError("sn_fj_validate: the Join node '" + sn.nodes[j - 1].name +
                                   "' has no matched Fork");
        // JOIN STRATEGY. PARTIAL is served: `FjJoinParam::required[r][b]` is a
        // per-branch count, so a quorum is that count LOWERED and not a
        // different mechanism. What it must be is REACHABLE -- a quorum above
        // what the fork is certain to emit would never fire, and a join that
        // never fires turns the chain into an absorbing set rather than an
        // error.
        typename std::map<std::size_t, typename NetworkStruct<T>::JoinDecl>::const_iterator jd =
            sn.joindecl.find(j);
        if (jd != sn.joindecl.end()) {
            if (jd->second.strategy != lang::JoinStrategy::STD &&
                jd->second.strategy != lang::JoinStrategy::PARTIAL)
                throw UnsupportedError(
                    "sn_fj_validate: only JoinStrategy STD and PARTIAL are supported by the exact "
                    "fork-join implementation; the Join node '" + sn.nodes[j - 1].name +
                    "' declares neither");
            if (jd->second.strategy == lang::JoinStrategy::PARTIAL && jd->second.quorum > 0.0 &&
                fsrc != 0) {
                const qn::NodeDef& fk2 = sn.nodes[fsrc - 1];
                const qn::ForkParam<T>* fp2 = sn.fork_param_of(fsrc);
                double emitted = 0.0;
                if (fp2 == 0) {
                    for (std::size_t nd = 1; nd <= sn.nodes.size(); ++nd)
                        for (std::size_t r = 1; r <= sn.nclasses; ++r)
                            if (num_traits<T>::to_double(sn.rtnodes((fsrc - 1) * sn.nclasses + r - 1,
                                                                    (nd - 1) * sn.nclasses + r - 1)) >
                                0) {
                                emitted += fk2.tasks_per_link;
                                break;
                            }
                } else {
                    for (std::size_t k = 0; k < fp2->fan_out_link.rows(); ++k)
                        for (std::size_t r = 0; r < fp2->fan_out_link.cols(); ++r)
                            if (num_traits<T>::to_double(fp2->fan_out_prob(k, r)) == 1.0)
                                emitted += num_traits<T>::to_double(fp2->fan_out_link(k, r));
                }
                if (jd->second.quorum > emitted)
                    throw UnsupportedError(
                        "sn_fj_validate: the Join node '" + sn.nodes[j - 1].name +
                        "' asks for more siblings than the Fork node '" + sn.nodes[fsrc - 1].name +
                        "' is certain to emit, so it could never fire");
            }
        }
    }
}

/**
 * Can the exact fork-join construction be asked for this model?
 *
 * The fork-join model class `sn_fj_validate` admits, asked as a predicate
 * rather than thrown, so that a REPORT can reach it: `solver_ctmc_analyzer` and
 * the SSA runner reach the SAME rules through `fj_tag`, which calls the
 * validator on its first line. The message the validator throws is about "the
 * exact fork-join implementation", which both share, so this predicate belongs
 * beside it rather than inside either solver.
 *
 * IT WRAPS THE VALIDATOR RATHER THAN RESTATING IT, and that is the point: the
 * rules are eight and they move (pairing, join strategy, tasks-per-link, branch
 * probability, open classes through a fork), so a second copy would be a second
 * thing to keep in step. There is exactly one body of rules and two ways in --
 * one that throws, for the run, and this one, which answers.
 *
 * WHAT IT REFUSES AND WHY THE ANALYZER IS RIGHT TO. The fork-join PAIRING is a
 * declaration carried by the Join (`Join(model, name, fork)` in all four
 * codebases), not a derivation from the routing: a nested model such as
 * fj_basic_nesting has two forks and two joins whose pairing the routing alone
 * does not determine. So a Join built without naming its fork leaves `sn.fj`
 * empty, and "has no matched Join" is the honest answer to a model that
 * declares none -- not a topology test that failed to see one.
 *
 * @param sn the refreshed struct of the model
 * @return an empty string when the fork-join construction may run
 */
template <class T>
std::string sn_fj_supports(const NetworkStruct<T>& sn) {
    bool any_fj = sn.has_fork();
    for (std::size_t i = 0; i < sn.nodes.size() && !any_fj; ++i)
        if (sn.nodes[i].nodetype == NodeType::Join) any_fj = true;
    if (!any_fj) return std::string();
    try {
        sn_fj_validate(sn);
    } catch (const line::Error& e) {
        // the validator's own refusal, turned into an answer
        return std::string(e.what());
    }
    return std::string();
}

/**
 * Port of `ModelAdapter.fjtag`.
 *
 * @param sn a refreshed struct whose Fork nodes all have a matched Join
 */
template <class T>
FjTagged<T> fj_tag(const NetworkStruct<T>& sn) {
    sn_fj_validate(sn);

    const std::size_t K = sn.nclasses;
    const std::size_t I = sn.nodes.size();
    const std::vector<std::vector<double>> Vn = fj_detail::node_visit_sum(sn);
    const T zero = num_traits<T>::from_int(0);
    const double inf = std::numeric_limits<double>::infinity();

    FjTagged<T> out;
    out.korig = K;
    out.V = sn;
    NetworkStruct<T>& V = out.V;
    V.isfjaugmented = true;

    // A Fork holds the parent job for the instant between its arrival and the
    // firing, so in the augmented model it is STATEFUL. `stateful_nodes` must
    // stay ASCENDING, because `stateful_index` is a position in it and every
    // network state row is indexed by that position.
    for (std::size_t i = 0; i < I; ++i)
        if (V.nodes[i].nodetype == NodeType::Fork) V.nodes[i].stateful = true;
    V.stateful_nodes.clear();
    for (std::size_t i = 0; i < I; ++i)
        if (V.nodes[i].stateful) V.stateful_nodes.push_back(i + 1);

    // One row per (fork, class): the branch structure the post-edits and the
    // firing list are both built from.
    struct ForkInfo {
        std::size_t f = 0, j = 0, r = 0, w = 1;
        /** Per-branch tasksPerLink, empty when every branch carries `w`. */
        std::vector<std::size_t> wlink;
        std::vector<std::size_t> branchheads;
        std::vector<std::vector<std::size_t>> branchsets;
        std::vector<std::vector<std::size_t>> auxmatrix;  // B x T
    };
    std::vector<ForkInfo> info;

    for (std::size_t f = 1; f <= I; ++f) {
        if (sn.nodes[f - 1].nodetype != NodeType::Fork) continue;
        std::size_t j = 0;
        for (std::size_t p = 0; p < sn.fj.size(); ++p)
            if (sn.fj[p].first == f) j = sn.fj[p].second;
        const std::size_t w = static_cast<std::size_t>(sn.nodes[f - 1].tasks_per_link);

        for (std::size_t r = 1; r <= K; ++r) {
            if (!(Vn[f - 1][r - 1] > 0)) continue;
            ForkInfo fi;
            fi.f = f;
            fi.j = j;
            fi.r = r;
            fi.w = w == 0 ? 1 : w;

            // Branch heads: the nodes the fork routes class r to directly.
            for (std::size_t nd = 1; nd <= I; ++nd)
                if (num_traits<T>::to_double(sn.rtnodes((f - 1) * K + (r - 1),
                                                        (nd - 1) * K + (r - 1))) > 0)
                    fi.branchheads.push_back(nd);
            const std::size_t B = fi.branchheads.size();
            if (B < 2)
                throw UnsupportedError(
                    "fj_tag: the Fork node '" + sn.nodes[f - 1].name +
                    "' has a single output link for class '" + sn.classes[r - 1].name +
                    "'. A degenerate fork is not a fork -- remove it, or give it a second branch");

            // PER-BRANCH TASKS PER LINK. `w` is the node-wide count and stays
            // the mean the scalar slot carries; `wlink` is the count of the link
            // that actually feeds branch b, and it is filled only when the
            // branches disagree, so a plain fork keeps exactly the shape it had.
            {
                const qn::ForkParam<T>* fkn = sn.fork_param_of(f);
                if (fkn != 0) {
                    std::vector<std::size_t> wv(B, 0);
                    bool uniform = true;
                    for (std::size_t b = 0; b < B; ++b) {
                        const double lv = num_traits<T>::to_double(
                            fkn->fan_out_link(fi.branchheads[b] - 1, r - 1));
                        wv[b] = static_cast<std::size_t>(lv + 0.5);
                        if (wv[b] == 0) wv[b] = 1;
                        if (wv[b] != wv[0]) uniform = false;
                    }
                    if (uniform) {
                        fi.w = wv[0];
                    } else {
                        fi.wlink = wv;
                    }
                }
            }

            // Branch discovery: the class-r reachable closure from each head, up
            // to but excluding the join.
            fi.branchsets.assign(B, std::vector<std::size_t>());
            for (std::size_t b = 0; b < B; ++b) {
                std::vector<std::size_t> visitset(1, fi.branchheads[b]);
                std::vector<std::size_t> frontier(1, fi.branchheads[b]);
                while (!frontier.empty()) {
                    const std::size_t cn = frontier.front();
                    frontier.erase(frontier.begin());
                    if (sn.nodes[cn - 1].nodetype == NodeType::Fork)
                        throw UnsupportedError(
                            "fj_tag: nested fork-join (the Fork node '" + sn.nodes[cn - 1].name +
                            "' sits on a branch of '" + sn.nodes[f - 1].name +
                            "') is not supported: a sibling would need a tag from each enclosing "
                            "fork and this transform mints one tag dimension");
                    if (sn.nodes[cn - 1].nodetype == NodeType::Join && cn != j)
                        throw UnsupportedError(
                            "fj_tag: overlapping fork-join pairs (the Join node '" +
                            sn.nodes[cn - 1].name + "' sits on a branch of '" +
                            sn.nodes[f - 1].name + "', which closes at '" + sn.nodes[j - 1].name +
                            "') are not supported");
                    for (std::size_t nd = 1; nd <= I; ++nd)
                        for (std::size_t s = 1; s <= K; ++s) {
                            if (!(num_traits<T>::to_double(
                                      sn.rtnodes((cn - 1) * K + (r - 1), (nd - 1) * K + (s - 1))) >
                                  0))
                                continue;
                            if (s != r)
                                throw UnsupportedError(
                                    "fj_tag: class switching between the Fork node '" +
                                    sn.nodes[f - 1].name + "' and its Join is not supported: the "
                                    "sibling classes are minted per ORIGINAL class, so a sibling "
                                    "that switches has no auxiliary twin to switch into");
                            if (nd == j) continue;
                            if (std::find(visitset.begin(), visitset.end(), nd) == visitset.end()) {
                                visitset.push_back(nd);
                                frontier.push_back(nd);
                            }
                        }
                }
                // Every branch node must REACH the join. A sibling that can be
                // trapped away from it never joins, so the parent never
                // completes and the chain has an absorbing set the reference
                // refuses rather than solves.
                std::vector<std::size_t> canreach(1, j);
                bool changed = true;
                while (changed) {
                    changed = false;
                    for (std::size_t a = 0; a < visitset.size(); ++a) {
                        const std::size_t cn = visitset[a];
                        if (std::find(canreach.begin(), canreach.end(), cn) != canreach.end())
                            continue;
                        for (std::size_t b2 = 0; b2 < canreach.size(); ++b2)
                            if (num_traits<T>::to_double(
                                    sn.rtnodes((cn - 1) * K + (r - 1),
                                               (canreach[b2] - 1) * K + (r - 1))) > 0) {
                                canreach.push_back(cn);
                                changed = true;
                                break;
                            }
                    }
                }
                for (std::size_t a = 0; a < visitset.size(); ++a)
                    if (std::find(canreach.begin(), canreach.end(), visitset[a]) == canreach.end())
                        throw UnsupportedError(
                            "fj_tag: the Join node '" + sn.nodes[j - 1].name +
                            "' is unreachable from the branch node '" +
                            sn.nodes[visitset[a] - 1].name +
                            "'; a sibling trapped there never joins and its parent never completes");
                fi.branchsets[b] = visitset;
            }

            // The tag pool: the population of the CHAIN, not of the class.
            std::size_t c = 0;
            for (std::size_t cc = 0; cc < sn.chains.size(); ++cc)
                if (sn.chains[cc][r - 1]) { c = cc + 1; break; }
            double tot = 0;
            if (c == 0) {
                tot = sn.njobs()[r - 1];
            } else {
                for (std::size_t s = 0; s < K; ++s)
                    if (sn.chains[c - 1][s]) tot += sn.njobs()[s];
            }
            const std::size_t Tt = static_cast<std::size_t>(tot + 0.5);
            if (Tt == 0)
                throw InputError("fj_tag: the chain of class '" + sn.classes[r - 1].name +
                                 "' routed through '" + sn.nodes[f - 1].name +
                                 "' carries no jobs, so no fork firing can ever occur");

            fi.auxmatrix.assign(B, std::vector<std::size_t>(Tt, 0));
            for (std::size_t t = 1; t <= Tt; ++t)
                for (std::size_t b = 0; b < B; ++b) {
                    JobClass ac = sn.classes[r - 1];
                    ac.name = sn.classes[r - 1].name + "_f" + std::to_string(f) + "_b" +
                              std::to_string(b + 1) + "_t" + std::to_string(t);
                    ac.type = JobClassType::CLOSED;
                    ac.population = 0.0;
                    ac.completes = false;
                    ac.is_ref_class = false;
                    const std::size_t a = V.add_class(ac);
                    fi.auxmatrix[b][t - 1] = a;

                    // The sibling is served exactly as the parent would have
                    // been, at every station of its own branch.
                    for (std::size_t x = 0; x < fi.branchsets[b].size(); ++x) {
                        const std::size_t cn = fi.branchsets[b][x];
                        const std::size_t ist = sn.nodes[cn - 1].station;
                        if (ist == 0 || sn.nodes[cn - 1].nodetype == NodeType::Join) continue;
                        V.set_service(ist, a, sn.service[ist - 1][r - 1]);
                    }
                    // The sibling routing is the parent's, restricted to the
                    // branch. The auxiliary class TERMINATES at the join: it has
                    // no outgoing row there, because the join consumes it and
                    // releases a parent instead.
                    for (std::size_t x = 0; x < fi.branchsets[b].size(); ++x) {
                        const std::size_t cn = fi.branchsets[b][x];
                        for (std::size_t nd = 1; nd <= I; ++nd) {
                            const T p = sn.get_route(r, r, cn, nd);
                            if (num_traits<T>::to_double(p) != 0) V.set_route(a, a, cn, nd, p);
                        }
                    }
                }
            info.push_back(fi);
        }
    }

    const std::size_t Kaug = V.classes.size();
    out.fjclassmap.assign(Kaug, 0);
    out.fjforkmap.assign(Kaug, 0);
    out.fjjoinmap.assign(Kaug, 0);
    out.fjbranchmap.assign(Kaug, 0);
    out.fjtagmap.assign(Kaug, 0);
    for (std::size_t row = 0; row < info.size(); ++row) {
        const ForkInfo& fi = info[row];
        for (std::size_t b = 0; b < fi.auxmatrix.size(); ++b)
            for (std::size_t t = 0; t < fi.auxmatrix[b].size(); ++t) {
                const std::size_t a = fi.auxmatrix[b][t];
                out.fjclassmap[a - 1] = fi.r;
                out.fjforkmap[a - 1] = fi.f;
                out.fjjoinmap[a - 1] = fi.j;
                out.fjbranchmap[a - 1] = b + 1;
                out.fjtagmap[a - 1] = t + 1;
            }
    }

    // Pad every per-class table to the widened class set BEFORE the refresh, so
    // that an untouched auxiliary class looks exactly as the refresh would have
    // found it: unbounded capacity, derived drop rule, no routing strategy of its
    // own. `add_class` grows only the service table.
    for (std::size_t i = 0; i < V.stations.size(); ++i) {
        Station<T>& st = V.stations[i];
        if (st.classcap.size() < Kaug) st.classcap.resize(Kaug, inf);
        if (st.droprule.size() < Kaug) st.droprule.resize(Kaug, 0);
        if (!st.schedparam.empty() && st.schedparam.size() < Kaug)
            st.schedparam.resize(Kaug, zero);
        if (!st.cdscalingpeak.empty() && st.cdscalingpeak.size() < Kaug)
            st.cdscalingpeak.resize(Kaug, zero);
        if (!st.jdscalingpeak.empty() && st.jdscalingpeak.size() < Kaug)
            st.jdscalingpeak.resize(Kaug, zero);
    }
    for (std::size_t row = 0; row < info.size(); ++row) {
        const ForkInfo& fi = info[row];
        for (std::size_t b = 0; b < fi.auxmatrix.size(); ++b)
            for (std::size_t t = 0; t < fi.auxmatrix[b].size(); ++t) {
                const std::size_t a = fi.auxmatrix[b][t];
                for (std::size_t x = 0; x < fi.branchsets[b].size(); ++x) {
                    const std::size_t cn = fi.branchsets[b][x];
                    if (sn.nodes[cn - 1].station == 0) continue;
                    Station<T>& st = V.stations[sn.nodes[cn - 1].station - 1];
                    if (!st.schedparam.empty())
                        st.schedparam[a - 1] = st.schedparam[fi.r - 1];
                }
            }
    }
    for (std::size_t i = 0; i < I; ++i) {
        NodeDef& nd = V.nodes[i];
        if (nd.routing.empty()) continue;
        if (nd.routing.size() < Kaug) nd.routing.resize(Kaug, RoutingStrategy::PROB);
        for (std::size_t row = 0; row < info.size(); ++row) {
            const ForkInfo& fi = info[row];
            for (std::size_t b = 0; b < fi.auxmatrix.size(); ++b)
                for (std::size_t t = 0; t < fi.auxmatrix[b].size(); ++t)
                    nd.routing[fi.auxmatrix[b][t] - 1] =
                        sn.nodes[i].routing.size() >= fi.r ? sn.nodes[i].routing[fi.r - 1]
                                                           : RoutingStrategy::PROB;
        }
    }
    // An explicit ClassSwitch matrix is (nclasses x nclasses); widen it with the
    // identity on the auxiliary block. A sibling never switches class -- the
    // transform refuses a branch that switches -- so the identity is not a
    // default, it is the only admissible row.
    for (typename std::map<std::size_t, Matrix<T>>::iterator it = V.csmatrix.begin();
         it != V.csmatrix.end(); ++it) {
        const Matrix<T> old = it->second;
        Matrix<T> C(Kaug, Kaug, zero);
        for (std::size_t x = 0; x < old.rows() && x < Kaug; ++x)
            for (std::size_t y = 0; y < old.cols() && y < Kaug; ++y) C(x, y) = old(x, y);
        for (std::size_t a = K; a < Kaug; ++a) C(a, a) = num_traits<T>::from_int(1);
        it->second = C;
    }

    V.nclasses = Kaug;
    V.refresh_struct();

    // ---- post-edits, exactly as `fjtag.m` applies them after getStruct ----
    //
    // WHY THE VISITS ARE OVERWRITTEN. An auxiliary class has population 0 and a
    // routing block that terminates at the join, so it is a TRANSIENT class whose
    // visit ratios are not defined by any traffic equation -- the refresh either
    // divides by a zero visit at the reference station or spreads mass over a
    // sub-stochastic block. The engines read auxiliary visits only as a
    // zero-versus-nonzero support gate ("may this class appear at this node"), so
    // that is what is written: 1 on the branch support, 0 elsewhere. The ORIGINAL
    // classes' visits are restored from the pre-augmentation struct for the same
    // reason: making the fork stateful changes the stateful index space and the
    // fork's own row, which the original model's visits never had.
    for (std::size_t r = 1; r <= K; ++r) {
        std::size_t corig = 0, cnew = 0;
        for (std::size_t cc = 0; cc < sn.chains.size(); ++cc)
            if (sn.chains[cc][r - 1]) { corig = cc + 1; break; }
        for (std::size_t cc = 0; cc < V.chains.size(); ++cc)
            if (V.chains[cc][r - 1]) { cnew = cc + 1; break; }
        if (corig == 0 || cnew == 0) continue;
        for (std::size_t isf = 1; isf <= V.stateful_nodes.size(); ++isf) {
            const std::size_t ind = V.stateful_nodes[isf - 1];
            const std::size_t isf_old = sn.stateful_index(ind);
            if (isf_old != 0) {
                V.visits[cnew - 1](isf - 1, r - 1) = sn.visits[corig - 1](isf_old - 1, r - 1);
            } else {
                // The stateful Fork had no row before. Its visit is read only as
                // a capacity gate, so a nonzero marker is all it needs.
                V.visits[cnew - 1](isf - 1, r - 1) =
                    num_traits<T>::from_double(Vn[ind - 1][r - 1] > 0 ? 1.0 : 0.0);
            }
        }
        for (std::size_t i = 0; i < I; ++i)
            V.nodevisits[cnew - 1](i, r - 1) = sn.nodevisits[corig - 1](i, r - 1);
    }
    for (std::size_t row = 0; row < info.size(); ++row) {
        const ForkInfo& fi = info[row];
        std::size_t cnew = 0;
        for (std::size_t cc = 0; cc < V.chains.size(); ++cc)
            if (V.chains[cc][fi.r - 1]) { cnew = cc + 1; break; }
        if (cnew == 0) continue;
        for (std::size_t b = 0; b < fi.auxmatrix.size(); ++b) {
            std::vector<std::size_t> support = fi.branchsets[b];
            support.push_back(fi.j);
            // THIS BRANCH's count, not the node-wide one. Capping the auxiliary
            // class of a branch that was sent 3 tasks at the node-wide mean of 2
            // makes the third sibling unplaceable, and the chain deadlocks with
            // no error to say why.
            const std::size_t wb = fi.wlink.empty() ? fi.w : fi.wlink[b];
            for (std::size_t t = 0; t < fi.auxmatrix[b].size(); ++t) {
                const std::size_t a = fi.auxmatrix[b][t];
                for (std::size_t c2 = 0; c2 < V.nchains; ++c2) {
                    for (std::size_t x = 0; x < V.visits[c2].rows(); ++x)
                        V.visits[c2](x, a - 1) = zero;
                    for (std::size_t x = 0; x < V.nodevisits[c2].rows(); ++x)
                        V.nodevisits[c2](x, a - 1) = zero;
                }
                for (std::size_t i = 0; i < V.stations.size(); ++i)
                    V.classcap[i][a - 1] = 0.0;
                for (std::size_t x = 0; x < support.size(); ++x) {
                    const std::size_t cn = support[x];
                    V.nodevisits[cnew - 1](cn - 1, a - 1) = num_traits<T>::from_int(1);
                    const std::size_t isf = V.stateful_index(cn);
                    if (isf != 0)
                        V.visits[cnew - 1](isf - 1, a - 1) = num_traits<T>::from_int(1);
                    // Each auxiliary class holds at most `tasksPerLink` siblings
                    // network-wide: one tag is outstanding at a time under STD.
                    if (V.nodes[cn - 1].station != 0)
                        V.classcap[V.nodes[cn - 1].station - 1][a - 1] =
                            static_cast<double>(wb);
                }
            }
        }
        FjJoinParam& jp = out.joinparam[fi.j];
        jp.fork = fi.f;
        if (std::find(jp.origclasses.begin(), jp.origclasses.end(), fi.r) == jp.origclasses.end())
            jp.origclasses.push_back(fi.r);
        jp.auxmatrix[fi.r] = fi.auxmatrix;
        // Siblings of branch b a firing consumes. Under STD it is the tasks that
        // branch was sent -- the per-destination count when the fork declares
        // one, the node-wide count otherwise. Under PARTIAL the Join's own
        // quorum is what it waits for, so the quorum LOWERS this count rather
        // than replacing the mechanism.
        std::vector<std::size_t> reqb(fi.auxmatrix.size(), fi.w);
        if (!fi.wlink.empty())
            for (std::size_t b = 0; b < reqb.size() && b < fi.wlink.size(); ++b)
                reqb[b] = fi.wlink[b];
        typename std::map<std::size_t, typename NetworkStruct<T>::JoinDecl>::const_iterator jdq =
            sn.joindecl.find(fi.j);
        if (jdq != sn.joindecl.end() && jdq->second.strategy == lang::JoinStrategy::PARTIAL &&
            jdq->second.quorum > 0.0) {
            const std::size_t q = static_cast<std::size_t>(jdq->second.quorum + 0.5);
            for (std::size_t b = 0; b < reqb.size(); ++b)
                if (q < reqb[b]) reqb[b] = q;
        }
        jp.required[fi.r] = reqb;
    }

    // The firing list: one entry per (fork, class, tag).
    for (std::size_t row = 0; row < info.size(); ++row) {
        const ForkInfo& fi = info[row];
        const std::size_t Tt = fi.auxmatrix.empty() ? 0 : fi.auxmatrix[0].size();
        for (std::size_t t = 1; t <= Tt; ++t) {
            FjSync<T> e;
            e.fork = fi.f;
            e.join = fi.j;
            e.cls = fi.r;
            e.tag = t;
            e.branchheads = fi.branchheads;
            for (std::size_t b = 0; b < fi.auxmatrix.size(); ++b)
                e.auxclasses.push_back(fi.auxmatrix[b][t - 1]);
            e.auxall = fi.auxmatrix;
            e.weight = fi.w;
            e.weightlink = fi.wlink;  // empty whenever every branch agrees
            out.fjsync.push_back(e);
        }
    }

    V.fjclassmap = out.fjclassmap;
    V.fjjoinparam = out.joinparam;
    return out;
}

}  // namespace qn
}  // namespace line

#endif  // LINE_LANG_QN_FJ_TAG_H
