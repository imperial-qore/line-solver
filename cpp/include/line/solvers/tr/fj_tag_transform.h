// Copyright (c) 2012-2026, QORE Lab, Imperial College London
// All rights reserved.
#ifndef LINE_SOLVERS_TR_FJ_TAG_TRANSFORM_H
#define LINE_SOLVERS_TR_FJ_TAG_TRANSFORM_H

/**
 * @file
 * Fork-join TAG AUGMENTATION: the fold-back half of the transform/lift pair
 * that CTMC and SSA share.
 *
 * `qn::fj_tag` is the OTHER fork-join route. Where `mmt` and `ht` drive an
 * outer fixed point for MVA, NC and Fluid (`solvers/mva/fj_driver.h`), the tag
 * augmentation is EXACT and single pass: it rewrites the struct so each sibling
 * branch carries its own auxiliary class, the engine runs unchanged on that
 * struct, and the auxiliary columns are folded back at the end.
 *
 * This is deliberately NOT built on the `fj_driver.h` shape. That driver owns a
 * loop and takes the inner solve as a template parameter because MMT re-solves
 * a transformed model repeatedly. `fj_tag` substitutes the struct and then the
 * caller's own engine runs on it to completion: there is no callback seam and
 * no second pass, so the reusable unit is the fold-back, not a driver.
 *
 * THE ONE C++-SPECIFIC POINT. CTMC and SSA carry DIFFERENT result containers
 * (`ctmc::CtmcAvg<T>` over `Matrix<T>`, `ssa::SsaSolution` over
 * `Matrix<double>`), which is why this is templated on the container rather
 * than taking one. It is the same obstacle python meets, where SolverMVA keeps
 * a dict and SolverNC a dataclass.
 *
 * Mirrors MATLAB `matlab/src/solvers/TR/solver_tr_fjtag_analyzer.m`, python
 * `line_solver/solvers/fjtag_transform.py` and the JAR
 * `jline.solvers.tr.FJTagTransform`.
 */

#include <cstddef>
#include <type_traits>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/util/matrix.h"

namespace line {
namespace tr {

/**
 * Whether the model needs the tag augmentation at all.
 *
 * A Fork or a Join sends the model down the `qn::fj_tag` route; neither means
 * it stays on the ordinary one. Four analyzers carried their own copy of this
 * loop (CTMC steady-state and transient, SSA serial, the chain tables), which
 * is the C++ share of the bookkeeping that MATLAB, the JAR and python absorbed
 * into an `expand` phase.
 *
 * There is deliberately NO `expand` here to match those three. `qn::fj_tag`
 * ALREADY returns the context they had to build by hand: `FjTagged` carries the
 * augmented struct, the sync list, `fjclassmap` and `korig` together, so
 * wrapping it would add a layer without removing a duplication.
 */
template <class T>
bool has_fork_join(const qn::NetworkStruct<T>& sn) {
    if (sn.has_fork()) return true;
    for (std::size_t i = 0; i < sn.nodes.size(); ++i)
        if (sn.nodes[i].nodetype == lang::NodeType::Join) return true;
    return false;
}

/**
 * Reduce the augmented metrics onto the original classes.
 *
 * A sibling class is a PART of the class it was forked from, so its queue
 * length, utilization and throughput are exact aggregates and simply add.
 * Response time is NOT additive and is recomputed by Little's law afterwards.
 * The SYSTEM metrics are truncated rather than summed: `XN` is the departure
 * rate at the parent class's reference station, which no sibling visits, and
 * adding a branch's throughput to it would count each forked task once per
 * branch.
 *
 * @tparam T   the struct's numeric type
 * @tparam Avg the caller's result container, with QN, UN, TN, RN, CN, XN
 */
template <class T, class Avg>
void fj_foldback(const qn::NetworkStruct<T>& sn, Avg& a,
                 const std::vector<std::size_t>& fjclassmap, std::size_t korig) {
    typedef typename std::remove_cv<
        typename std::remove_reference<decltype(a.QN(0, 0))>::type>::type S;
    const S zero = num_traits<S>::from_int(0);
    const std::size_t M = a.QN.rows();
    for (std::size_t x = 0; x < fjclassmap.size(); ++x) {
        const std::size_t r = fjclassmap[x];
        if (r == 0) continue;
        for (std::size_t i = 0; i < M; ++i) {
            a.QN(i, r - 1) += a.QN(i, x);
            a.UN(i, r - 1) += a.UN(i, x);
            a.TN(i, r - 1) += a.TN(i, x);
        }
    }
    Matrix<S> Q(M, korig, zero), U(M, korig, zero), TT(M, korig, zero), Rr(M, korig, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < korig; ++k) {
            Q(i, k) = a.QN(i, k);
            U(i, k) = a.UN(i, k);
            TT(i, k) = a.TN(i, k);
            Rr(i, k) = num_traits<S>::to_double(TT(i, k)) > 0 ? S(Q(i, k) / TT(i, k)) : zero;
        }
    a.QN = Q;
    a.UN = U;
    a.TN = TT;
    a.RN = Rr;
    // A Join's response time is QLen over the SIBLING arrival rate, not over
    // its own firing rate: a Join sees one arrival per sibling for every job it
    // releases. Applied here, on the FOLDED table, so both engines and every
    // in-process caller read the same convention.
    mva::sn_apply_join_respt(sn, a.QN, mva::sn_get_arvr_from_tput(sn, a.TN), a.RN);
    a.CN.resize(korig);
    a.XN.resize(korig);
}

}  // namespace tr
}  // namespace line

#endif  // LINE_SOLVERS_TR_FJ_TAG_TRANSFORM_H
