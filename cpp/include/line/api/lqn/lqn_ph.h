/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_LQN_LQN_PH_H
#define LINE_API_LQN_LQN_PH_H

/**
 * Phase-type composition of an LQN activity graph, the machinery behind
 * SolverLN method 'srvn.ph'.
 *
 * An entry becomes a Workflow whose leaves are its activities and, when
 * requested, its synchronous calls; the series-parallel reduction of that
 * workflow is then the exact law of the entry service time. Port of the MATLAB
 * lqn_entry_workflow.m, lqn_ph_serial_law.m and lqn_ph_moments.m, of the JAR
 * jline.api.lqn.LqnPh and of the Python line_solver.api.lqn.lqn_ph.
 */

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/lqn/lqn_struct.h"
#include "line/lang/workflow/workflow.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/lu.h"

namespace line {
namespace api {
namespace lqn {

using line::lang::Distrib;
using line::lang::GlobalConstants;
using line::lang::PrecedenceType;
using line::workflow::PhLaw;
using line::workflow::Precedence;
using line::workflow::SPNodeType;
using line::workflow::SPTree;
using line::workflow::Workflow;

/**
 * The (alpha, S) pair of a phase-type Distrib, the form the composition rules
 * take. For every Markovian family D1 = (-D0 e) alpha, so any row with a
 * positive exit rate recovers alpha when the parameters do not carry it.
 */
template <class T>
PhLaw<T> ph_law_of(const Distrib<T>& d) {
    const T zero = num_traits<T>::from_int(0);
    PhLaw<T> out;
    out.S = d.D0;
    const std::size_t n = d.D0.rows();
    if ((d.type == line::lang::ProcessType::PH || d.type == line::lang::ProcessType::APH ||
         d.type == line::lang::ProcessType::ME) &&
        d.params.size() == n) {
        out.alpha.assign(d.params.begin(), d.params.begin() + static_cast<long>(n));
        return out;
    }
    for (std::size_t i = 0; i < n; ++i) {
        T tot = zero;
        for (std::size_t j = 0; j < n; ++j) tot += d.D1(i, j);
        if (!(tot > zero)) continue;
        out.alpha.assign(n, zero);
        for (std::size_t j = 0; j < n; ++j) out.alpha[j] = T(d.D1(i, j) / tot);
        return out;
    }
    out.alpha.assign(n, zero);
    if (n > 0) out.alpha[0] = num_traits<T>::from_int(1);
    return out;
}

/** Activity graph of one entry, as a workflow plus its execution counts. */
template <class T>
struct EntryWorkflow {
    Workflow<T> wf{"entry"};
    /** Workflow activity index of each LQN activity, npos when absent. */
    std::unordered_map<std::size_t, std::size_t> act_idx_of;
    /** Workflow activity index of each call, npos when absent. */
    std::unordered_map<std::size_t, std::size_t> call_idx_of;
    /** Expected executions of each LQN activity per entry invocation. */
    std::unordered_map<std::size_t, T> execs;
    /** Expected executions of each call per entry invocation. */
    std::unordered_map<std::size_t, T> callexecs;
};

/**
 * Activity graph of LQN entry EIDX as a Workflow.
 *
 * The precedences are read from LqnStruct::precedences rather than
 * reconstructed from `graph`, whose loop back-edges carry probabilities and not
 * counts.
 *
 * WITH_CALLS true expands every synchronous call of an activity into a leaf of
 * its own, placed in series after the activity, so that the call response law
 * and the host demand law stay separable across iterations. False keeps only
 * the host demands, which is the processor-demand law of the entry: the host is
 * released while a call is outstanding.
 */
template <class T>
EntryWorkflow<T> entry_workflow(const ::line::lqn::LqnStruct<T>& lqn, std::size_t eidx, bool with_calls) {
    const std::size_t tidx = lqn.parent[eidx];
    const std::vector<std::size_t>& acts = lqn.actsof[eidx];
    if (acts.empty())
        throw InputError("Entry " + lqn.hashnames[eidx] + " binds no activity.");

    EntryWorkflow<T> out;
    out.wf = Workflow<T>(lqn.names[eidx] + ".Workflow");

    std::unordered_map<std::size_t, std::string> head_name, tail_name;
    std::unordered_map<std::size_t, std::size_t> act_of_name_idx;  // by activity index
    std::unordered_map<std::string, std::size_t> act_of_name;

    for (std::size_t aidx : acts) {
        const std::string& nm = lqn.names[aidx];
        out.act_idx_of[aidx] = out.wf.add_activity(nm, lqn.hostdem[aidx]);
        head_name[aidx] = nm;
        tail_name[aidx] = nm;
        act_of_name[nm] = aidx;
        act_of_name_idx[aidx] = aidx;
    }

    if (with_calls) {
        for (std::size_t aidx : acts) {
            std::vector<std::string> chain;
            chain.push_back(head_name[aidx]);
            for (std::size_t cidx : lqn.callsof[aidx]) {
                if (lqn.calltype[cidx] != line::lang::CallType::SYNC) {
                    continue;  // an asynchronous call blocks the caller for no time
                }
                const std::string& cnm = lqn.callhashnames[cidx];
                out.call_idx_of[cidx] = out.wf.add_activity(cnm, Distrib<T>::immediate());
                chain.push_back(cnm);
            }
            if (chain.size() > 1) {
                for (std::size_t k = 1; k < chain.size(); ++k) {
                    Precedence<T> p;
                    p.pre_acts.push_back(chain[k - 1]);
                    p.post_acts.push_back(chain[k]);
                    p.pre_type = PrecedenceType::PRE_SEQ;
                    p.post_type = PrecedenceType::POST_SEQ;
                    out.wf.add_precedence(p);
                }
                tail_name[aidx] = chain.back();
            }
        }
    }

    // Precedences of the task, restricted to the activities of this entry and
    // rewritten so that a predecessor is entered at its head and left at its tail
    if (tidx < lqn.precedences.size()) {
        for (const auto& prec : lqn.precedences[tidx]) {
            bool mine = true;
            for (std::size_t a : prec.preacts)
                if (act_of_name_idx.find(a) == act_of_name_idx.end()) mine = false;
            for (std::size_t a : prec.postacts)
                if (act_of_name_idx.find(a) == act_of_name_idx.end()) mine = false;
            if (!mine) continue;  // the precedence belongs to another entry of the same task
            Precedence<T> p;
            p.pre_type = prec.pretype;
            p.post_type = prec.posttype;
            p.pre_params = prec.preparams;
            p.post_params = prec.postparams;
            // A loop carries ONE count. LqnBuilder::loop repeats it once per body
            // activity, which the arc expansion reads position 0 of and ignores the
            // rest; the series-parallel parser instead requires exactly one.
            if (prec.posttype == PrecedenceType::POST_LOOP && p.post_params.size() > 1)
                p.post_params.resize(1);
            for (std::size_t a : prec.preacts) p.pre_acts.push_back(tail_name[a]);
            for (std::size_t a : prec.postacts) p.post_acts.push_back(head_name[a]);
            out.wf.add_precedence(p);
        }
    }

    // sp_tree() does not validate -- only to_ph() does -- and a quorum AND-join
    // is exactly what validation refuses, so ask for it explicitly.
    out.wf.validate();
    const SPTree<T>* tree = out.wf.sp_tree();
    if (tree == nullptr)
        throw UnsupportedError(
            "Entry " + lqn.hashnames[eidx] +
            " has a precedence graph that is not series-parallel, so its activity graph has no "
            "exact phase-type reduction. Use method='default'.");

    for (std::size_t aidx : acts)
        out.execs[aidx] = tree->execs[tree->leaf_of[out.act_idx_of[aidx]]];
    for (const auto& kv : out.call_idx_of)
        out.callexecs[kv.first] = tree->execs[tree->leaf_of[kv.second]];
    return out;
}

/** Recursive body of serial_law, declared first so serial_law can call it. */
template <class T>
PhLaw<T> detail_compose_serialized(Workflow<T>& wf, const SPTree<T>& tree,
                                   std::size_t k) {
    const auto& node = tree.nodes[k];
    switch (node.type) {
        case SPNodeType::LEAF:
            return wf.activity_at(node.act).ph_representation();
        case SPNodeType::SERIAL:
        case SPNodeType::PAR: {
            PhLaw<T> acc = detail_compose_serialized(wf, tree, node.kids[0]);
            for (std::size_t i = 1; i < node.kids.size(); ++i)
                acc = Workflow<T>::compose_serial(acc,
                                                  detail_compose_serialized(wf, tree, node.kids[i]));
            return acc;
        }
        case SPNodeType::OR: {
            std::vector<PhLaw<T>> laws;
            for (std::size_t kid : node.kids)
                laws.push_back(detail_compose_serialized(wf, tree, kid));
            return Workflow<T>::compose_mixture(laws, node.probs);
        }
        case SPNodeType::LOOP:
            return Workflow<T>::compose_loop_geometric(
                detail_compose_serialized(wf, tree, node.kids[0]), node.count);
    }
    throw InputError("Unknown series-parallel node type.");
}

/**
 * Composed law of a workflow in which the branches of an AND fork are SERIAL
 * rather than concurrent, that is, the total work the branches request rather
 * than the elapsed time until the last of them finishes.
 *
 * This is the law of the PROCESSOR demand of an LQN entry. Two branches of an
 * AND fork are two activity threads of the same task instance: they overlap in
 * time, so the entry response time is the maximum of the branches, but they run
 * on ONE processor, so the demand they place on it is the sum. Composing the
 * host law with Workflow::to_ph would charge the processor the maximum and let
 * the layer report a utilization below the true one, which no amount of
 * iterating recovers.
 *
 * Every other node keeps its own composition rule: an OR fork is a mixture, a
 * loop is a geometric compound, so the correlation within a branch survives.
 */
template <class T>
PhLaw<T> serial_law(Workflow<T>& wf) {
    const SPTree<T>* tree = wf.sp_tree();
    if (tree == nullptr)
        throw UnsupportedError("Workflow " + wf.name() +
                               " is not series-parallel, so it has no exact phase-type reduction.");
    return detail_compose_serialized(wf, *tree, tree->root);
}

/**
 * First two moments of a phase-type law without building a Distrib, which is
 * what the layered fixed point needs at every iteration for every composed
 * entry law. A defective ALPHA carries an atom at zero and contributes nothing
 * to either moment.
 *
 * Returns {mean, squared coefficient of variation}.
 */
template <class T>
std::pair<T, T> ph_moments(const std::vector<T>& alpha, const Matrix<T>& S) {
    const std::size_t n = S.rows();
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> LU = S;
    std::vector<std::size_t> piv;
    std::pair<T, T> bad(num_traits<T>::from_double(GlobalConstants::FineTol),
                        num_traits<T>::from_int(1));
    try {
        piv = lu_factor(LU);
    } catch (const std::exception&) {
        return bad;
    }
    std::vector<T> x1(n, num_traits<T>::from_int(1));
    lu_solve(LU, piv, x1);
    for (std::size_t i = 0; i < n; ++i) x1[i] = T(-x1[i]);
    T m1 = zero;
    for (std::size_t i = 0; i < n && i < alpha.size(); ++i) m1 += alpha[i] * x1[i];

    std::vector<T> x2 = x1;
    lu_solve(LU, piv, x2);
    T m2 = zero;
    for (std::size_t i = 0; i < n && i < alpha.size(); ++i) m2 += alpha[i] * x2[i];
    m2 = T(num_traits<T>::from_int(-2) * m2);

    const double m1d = num_traits<T>::to_double(m1);
    if (!std::isfinite(m1d) || m1d <= GlobalConstants::FineTol) return bad;
    T scv = T(m2 / (m1 * m1) - num_traits<T>::from_int(1));
    const double scvd = num_traits<T>::to_double(scv);
    if (!std::isfinite(scvd) || scvd <= GlobalConstants::FineTol)
        scv = num_traits<T>::from_double(GlobalConstants::FineTol);
    return std::make_pair(m1, scv);
}

}  // namespace lqn
}  // namespace api
}  // namespace line

#endif  // LINE_API_LQN_LQN_PH_H
