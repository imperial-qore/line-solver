/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_FJ_DRIVER_H
#define LINE_SOLVERS_MVA_FJ_DRIVER_H

/**
 * The fork-join fixed point that drives one inner MVA solve.
 *
 * Port of @@NetworkSolver/fjFixedPoint.m. The transform (fj_mmt.h) turns a model
 * with a Fork into a plain mixed queueing network carrying auxiliary open
 * classes; this driver solves that network repeatedly, each pass re-setting the
 * auxiliary arrival rates from the current forkLambda, recomputing the
 * synchronisation delay at the join from the branch response times the solve
 * reports, moving forkLambda halfway towards the join throughput, and merging
 * the auxiliary classes back into the ones they stand for.
 *
 * The inner solve is a callback so the same driver serves both entry points that
 * need it: SolverLN, which solves a layer with solver_mva_analyzer, and
 * SolverMVA over a general Network, which solves with mva_dispatch. The callback
 * takes the transformed model by reference (the driver mutates its services and
 * chains between passes) and returns an MvaSolution over the auxiliary-expanded
 * class set; everything else here -- the merge-back, the auxiliary-column drop
 * and the non-finite guard -- is codebase-independent and lives here once.
 */

#include "line/util/line_console.h"
#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/fj_mmt.h"
#include "line/solvers/mva/solver_mva.h"
#include "line/util/error.h"

namespace line {
namespace mva {

using lang::Distrib;
using lang::GlobalConstants;

/**
 * Port of ModelAdapter.findPathsCS' per-node charge: the queue length of the
 * merge set at this station divided by its throughput, which is the time a job
 * of that set spends at the node. See the reference comment carried in
 * solver_ln.h: a non-station node, and a station the merge set never completes
 * at, contribute nothing.
 */
template <class T>
T fj_node_time(const qn::NetworkStruct<T>& V, std::size_t nd,
               const std::vector<std::size_t>& merge, const Matrix<T>& QN,
               const Matrix<T>& TN) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t st = V.nodes[nd - 1].station;
    if (st == 0) return zero;
    T q = zero, t = zero;
    for (std::size_t k : merge) {
        q += QN(st - 1, k - 1);
        t += TN(st - 1, k - 1);
    }
    if (!(t > zero)) return zero;
    return T(q / t);
}

/**
 * Port of ModelAdapter.findPathsCS: the response time along every path from a
 * fork to ITS join, in the given class. `merge`'s first entry is the class being
 * traversed and its remaining entries the auxiliary classes of the forks whose
 * branches are being walked; the end node's own time is subtracted back off,
 * since the synchronisation delay is what is being computed and must not be
 * counted into the branch it delays.
 *
 * NESTED FORKS ARE COLLAPSED IN PLACE. Meeting a fork part-way along a branch,
 * the walk recurses to THAT fork's join, forms its `E[max]` there and then
 * continues from its join with `t0 + E[max]` -- so an inner fork contributes the
 * time its own branches take in parallel, not the time they would take in series.
 * It also WRITES the inner join's synchronisation delay for the merge set as a
 * side effect, which is why `V` is taken by non-const reference here: the
 * reference does the same, and it is the only place an inner join's delay is set
 * (the driver's own loop sets it only for outer forks).
 *
 * The `visited` set of (node, class) pairs is the reference's own cycle guard: a
 * path returning to a pair it already holds is a routing loop, such as the .Aux
 * self-loop of a call whose mean exceeds one, not a new branch. It replaces the
 * depth cap this port used while it handled a single fork, which could not tell a
 * legitimate deep nest from a cycle.
 */
template <class T>
void fj_find_paths(FjMmt<T>& tr, std::size_t curNode, std::size_t endNode, std::size_t curClass,
                   const std::vector<std::size_t>& merge, const Matrix<T>& QN, const Matrix<T>& TN,
                   const T& t0, std::vector<T>& out,
                   std::vector<std::pair<std::size_t, std::size_t>> visited) {
    qn::NetworkStruct<T>& V = tr.V;
    const T zero = num_traits<T>::from_int(0);
    if (curNode == endNode) {
        out.push_back(T(t0 - fj_node_time(V, curNode, merge, QN, TN)));
        return;
    }
    const std::pair<std::size_t, std::size_t> here(curNode, curClass);
    if (std::find(visited.begin(), visited.end(), here) != visited.end()) return;
    visited.push_back(here);

    for (std::size_t s = 1; s <= V.classes.size(); ++s)
        for (std::size_t nd = 1; nd <= V.nodes.size(); ++nd) {
            // CS-aware: a synthesized ClassSwitch on a branch carries its switch
            // in `csmatrix`, not in `P`, so a walk over `P` alone would end the
            // path there and lose the rest of the branch's response time.
            if (!(detail::fj_route_cs(V, curClass, s, curNode, nd) > zero)) continue;
            std::vector<std::size_t> m2 = merge;
            m2[0] = s;

            // Is `nd` another fork of this transform? Its node type is Router by
            // now, so the fork records are the only way to tell.
            std::size_t inner = tr.forks.size();
            for (std::size_t b = 0; b < tr.forks.size(); ++b)
                if (tr.forks[b].node == nd && tr.forks[b].joinNode != 0) inner = b;

            if (inner < tr.forks.size()) {
                const std::size_t ijoin = tr.forks[inner].joinNode;
                const std::size_t istat = tr.forks[inner].joinStation;
                // the inner fork's own auxiliary class for the class arriving at it
                std::size_t iaux = 0;
                for (std::size_t x : tr.auxclasses)
                    if (tr.fjforkmap[x] == inner && tr.fjclassmap[x] == s) iaux = x;
                std::vector<std::size_t> m3 = m2;
                if (iaux != 0) m3.push_back(iaux);

                std::vector<T> paths;
                fj_find_paths(tr, nd, ijoin, s, m3, QN, TN, zero, paths, visited);
                T d0 = zero;
                if (!paths.empty()) {
                    T mean = zero;
                    for (const T& x : paths) mean += x;
                    mean = T(mean / num_traits<T>::from_int(static_cast<long>(paths.size())));
                    // The join fires on the k-th branch completion, k = the branch
                    // count on a standard join and the declared quorum on a PARTIAL one.
                    d0 = fj_expected_ordstat(paths, fj_join_quorum(V, ijoin, paths.size()));
                    // The inner join's delay, for every class of the merge set.
                    // Note the reference charges E[max] - mean here WITHOUT the
                    // fanOut factor it applies at an outer fork. Under a quorum the
                    // k-th completion can precede a branch's own, and then the join
                    // adds no delay: the parent left on a sibling.
                    const Distrib<T> d =
                        fj_exp_fit_mean(d0 > mean ? T(d0 - mean) : num_traits<T>::from_int(0));
                    if (istat != 0)
                        for (std::size_t cls : m3) V.set_service(istat, cls, d);
                }
                fj_find_paths(tr, ijoin, endNode, s, m2, QN, TN, T(t0 + d0), out, visited);
            } else {
                fj_find_paths(tr, nd, endNode, s, m2, QN, TN,
                              T(t0 + fj_node_time(V, nd, m2, QN, TN)), out, visited);
            }
        }
}

/** The leading `n` columns of a class-indexed metric. */
template <class T>
Matrix<T> fj_leading_cols(const Matrix<T>& A, std::size_t n) {
    Matrix<T> out(A.rows(), n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < n && j < A.cols(); ++j) out(i, j) = A(i, j);
    return out;
}

/** The mixed absolute/relative stopping test of the fork-join loop. */
template <class T>
bool fj_converged(const Matrix<T>& A, const Matrix<T>& B, double iter_tol) {
    if (A.rows() != B.rows() || A.cols() != B.cols()) return false;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) {
            const double d = std::fabs(num_traits<T>::to_double(A(i, j)) -
                                       num_traits<T>::to_double(B(i, j)));
            if (d > GlobalConstants::Zero +
                        iter_tol * std::fabs(num_traits<T>::to_double(B(i, j))))
                return false;
        }
    return true;
}

/**
 * Port of `fjFixedPoint.m:130-136`: the firing throughput of fork `fa`, per
 * class of the BASE model.
 *
 * A fork is not a station and has no throughput of its own, so the reference
 * derives one from the visit ratios: the fork's node visits, divided by the
 * class's total visits to its reference station, times the throughput observed
 * there. The divisor is the CHAIN's visit total at the reference station and
 * the multiplier the CHAIN's throughput there, because visits are normalised
 * per chain and a class-only ratio would not be scale-free under class
 * switching.
 *
 * TWO INDEX SPACES MEET HERE, exactly as in the reference. `visits` is indexed
 * by STATEFUL node and `refstat` is a STATION, so the reference station must be
 * mapped through `stateful_of_station` before it indexes `visits`; `nodevisits`
 * is indexed by NODE, which is what `parent` names. The throughput `TN` comes
 * from the solve of the TRANSFORMED model and is indexed by its station rows,
 * which the transform leaves aligned with the base model's -- the same
 * assumption `fjFixedPoint.m` makes when it writes `TN(nodeToStation(joinIdx))`.
 *
 * `parent` is this fork's own index when it is not nested: an inner fork fires
 * as often as the outer one it sits behind, so the visits that measure it are
 * the OUTER fork's, which is what `sortForks` returns in `parent_forks`.
 */
template <class T>
std::vector<T> fj_fork_tput(const qn::NetworkStruct<T>& L, const FjMmt<T>& tr, const Matrix<T>& TN,
                            std::size_t fa) {
    const T zero = num_traits<T>::from_int(0);
    std::vector<T> tnfork(L.nclasses, zero);
    const std::size_t pnode = tr.forks[tr.forks[fa].parent].node;
    for (std::size_t c = 0; c < L.nchains; ++c) {
        const std::vector<std::size_t>& ic = L.inchain[c];
        for (std::size_t r : ic) {
            const std::size_t rs = L.classes[r - 1].refstat;
            const std::size_t isf = L.stateful_of_station(rs);
            T den = zero, tsum = zero;
            for (std::size_t k : ic) {
                den += L.visits[c](isf - 1, k - 1);
                tsum += TN(rs - 1, k - 1);
            }
            if (!(den > zero)) continue;
            tnfork[r - 1] = T(T(L.nodevisits[c](pnode - 1, r - 1) / den) * tsum);
        }
    }
    return tnfork;
}

/** A metric read out of a solve, with a non-finite entry read as zero. */
template <class T>
T fj_finite_or_zero(const T& x) {
    if (!std::isfinite(num_traits<T>::to_double(x))) return num_traits<T>::from_int(0);
    return x;
}

/**
 * Port of the `heidelberger-trivedi` arm of `fjFixedPoint.m:212-255`: the
 * synchronisation delays of one pass, written onto the transformed model.
 *
 * The branch response time of an auxiliary class is its TOTAL response time over
 * the model, less what it spends at the join (the residual synchronisation delay
 * this function is about to overwrite) and at the auxiliary delay (the original
 * class's time outside the span). That leaves exactly the time the branch takes.
 * From those,
 *
 *   d0 = E[X_(k)]  fires the join. k is the branch count on a standard join and
 *                  the declared quorum on a PARTIAL one.
 *   di = d0 - ri   is what branch i still waits at the join once it has finished.
 *   r0             is the original class's cycle response time less its own time
 *                  at the join, i.e. everything outside the span, which is what
 *                  the auxiliary delay must hold so that the auxiliary token
 *                  cycles at the original's rate.
 *
 * The join charges `d0 * fanOut` to the ORIGINAL class, which is the whole span:
 * `fj_ht` routed it straight past the branches. `fanOut` is 1 on every model that
 * reaches here -- `fj_ht` refuses tasksPerLink > 1 by name -- and is kept as a
 * named factor because the reference multiplies by it.
 *
 * UNDER A QUORUM `d0` can precede a branch's own completion, and then that branch
 * waits no further; `di` floors at zero, which is a boundary of the transform and
 * not a choice. The reference raises a line_warning there; this port has no
 * warning channel out of the fixed point (see mva_dispatch.h).
 */
template <class T>
void fj_ht_sync_delays(const qn::NetworkStruct<T>& L, FjMmt<T>& tr, const MvaSolution<T>& out) {
    const T zero = num_traits<T>::from_int(0);
    qn::NetworkStruct<T>& V = tr.V;
    for (std::size_t fa = 0; fa < tr.forks.size(); ++fa) {
        const std::size_t fnode = tr.forks[fa].node;
        const std::size_t jnode = tr.forks[fa].joinNode;
        if (jnode == 0) continue;
        const std::size_t jstat = tr.forks[fa].joinStation;
        const std::size_t adstat = tr.auxDelayStation[jnode];
        for (std::size_t c = 0; c < L.nchains; ++c) {
            const std::vector<std::size_t>& ic = L.inchain[c];
            for (std::size_t xi = 0; xi < ic.size(); ++xi) {
                const std::size_t r = ic[xi];
                if (num_traits<T>::to_double(L.nodevisits[c](fnode - 1, r - 1)) == 0.0) continue;
                // THIS FORK's auxiliary classes for r. The reference selects them
                // by `fjclassmap == r` alone, which on a model with a SECOND fork
                // over the same class would mix the two forks' branches into one
                // order statistic; the fork is added here because H-T charges each
                // fork at its own join.
                std::vector<std::size_t> br;
                for (std::size_t s : tr.auxclasses)
                    if (tr.fjforkmap[s] == fa && tr.fjclassmap[s] == r) br.push_back(s);
                if (br.empty()) continue;

                std::vector<T> ri(br.size(), zero);
                for (std::size_t b = 0; b < br.size(); ++b) {
                    T acc = zero;
                    for (std::size_t i = 0; i < out.R.rows(); ++i)
                        acc += fj_finite_or_zero(out.R(i, br[b] - 1));
                    // the two subtracted terms are read RAW, as in the reference:
                    // its NaN sweep applies to the summed matrix only
                    acc = T(acc - out.R(adstat - 1, br[b] - 1));
                    ri[b] = T(acc - out.R(jstat - 1, br[b] - 1));
                }
                const T d0 = fj_expected_ordstat(ri, fj_join_quorum(L, jnode, br.size()));

                T r0 = zero;
                for (std::size_t i = 0; i < out.R.rows(); ++i) {
                    T row = zero;
                    for (std::size_t yi = 0; yi < ic.size(); ++yi) row += out.R(i, ic[yi] - 1);
                    r0 += fj_finite_or_zero(row);
                }
                r0 = T(r0 - out.R(jstat - 1, r - 1));

                V.set_service(jstat, r,
                              fj_exp_fit_mean(T(d0 * num_traits<T>::from_double(
                                                        tr.forks[fa].fanOut))));
                const Distrib<T> outside = fj_exp_fit_mean(r0);
                for (std::size_t b = 0; b < br.size(); ++b) {
                    T di = T(d0 - ri[b]);
                    if (di < zero) di = zero;
                    V.set_service(jstat, br[b], fj_exp_fit_mean(di));
                    V.set_service(adstat, br[b], outside);
                }
            }
        }
    }
}

/**
 * Port of the `heidelberger-trivedi` arm of `fjFixedPoint.m:263-291`: fold the
 * auxiliary columns back into the classes they stand for.
 *
 * EVERY JOIN LOSES ITS ORIGINAL-CLASS METRICS FIRST. The original class was
 * routed straight past the branches and charged the whole span at the join, so
 * what the solve reports for it there is the transform's own bookkeeping and not
 * a queue the model has; the branches' figures, which the auxiliary classes
 * carry, are what belongs there. Throughput is the exception: the join's
 * original-class rate is how often the fork-join completes, so it is saved and
 * put back after the merge, where the reference puts it back too.
 *
 * The auxiliary DELAY rows are left in place rather than deleted. The reference
 * deletes them; they are appended after every base station, so the leading block
 * the caller reads is the same either way, and deleting them would only move the
 * join rows the restore step above indexes by their pre-deletion position.
 */
template <class T>
void fj_ht_merge(const qn::NetworkStruct<T>& L, const FjMmt<T>& tr, MvaSolution<T>& out) {
    const T zero = num_traits<T>::from_int(0);
    std::vector<bool> orig(tr.V.classes.size() + 1, false);
    for (std::size_t s : tr.auxclasses) orig[tr.fjclassmap[s]] = true;

    std::vector<std::vector<T>> tnJoin(tr.joinStations.size(),
                                       std::vector<T>(L.nclasses + 1, zero));
    for (std::size_t a = 0; a < tr.joinStations.size(); ++a) {
        const std::size_t js = tr.joinStations[a];
        for (std::size_t r = 1; r <= L.nclasses; ++r) {
            if (!orig[r]) continue;
            tnJoin[a][r] = out.Tp(js - 1, r - 1);
            out.Q(js - 1, r - 1) = zero;
            out.R(js - 1, r - 1) = zero;
            out.Tp(js - 1, r - 1) = zero;
            out.U(js - 1, r - 1) = zero;
        }
    }
    const std::size_t nrows = out.Tp.rows();
    for (std::size_t s : tr.auxclasses) {
        const std::size_t r = tr.fjclassmap[s];
        for (std::size_t i = 0; i < nrows; ++i) {
            out.Q(i, r - 1) = T(out.Q(i, r - 1) + out.Q(i, s - 1));
            out.U(i, r - 1) = T(out.U(i, r - 1) + out.U(i, s - 1));
            out.Tp(i, r - 1) = T(out.Tp(i, r - 1) + out.Tp(i, s - 1));
            out.R(i, r - 1) =
                out.Tp(i, r - 1) != zero ? T(out.Q(i, r - 1) / out.Tp(i, r - 1)) : zero;
        }
    }
    for (std::size_t a = 0; a < tr.joinStations.size(); ++a) {
        const std::size_t js = tr.joinStations[a];
        for (std::size_t r = 1; r <= L.nclasses; ++r)
            if (orig[r]) out.Tp(js - 1, r - 1) = tnJoin[a][r];
    }
    out.Q = fj_leading_cols(out.Q, L.nclasses);
    out.U = fj_leading_cols(out.U, L.nclasses);
    out.R = fj_leading_cols(out.R, L.nclasses);
    out.Tp = fj_leading_cols(out.Tp, L.nclasses);
}

/**
 * Drive the fork-join fixed point of a transformed model to convergence.
 *
 * @param L    the base model, whose services are copied back into the transform
 *             at the start of every solve (ModelAdapter.refreshServicesFromBase)
 * @param tr   the transform, mutated in place: its join and Source services carry
 *             the current pass's synchronisation delays and auxiliary arrivals
 * @param lam  the auxiliary arrival rates, `self.fjForkLambda`; warm-started by
 *             the caller across outer iterations and updated here in place
 * @param opt  MVA options (iter_max, iter_tol)
 * @param inner the inner solve over the auxiliary-expanded model
 * @return the merged MvaSolution over the base class set (auxiliary columns
 *         dropped, join and Source throughputs kept at their original-class value)
 */
template <class T, class InnerSolve>
MvaSolution<T> fj_fixed_point(const qn::NetworkStruct<T>& L, FjMmt<T>& tr,
                              std::vector<T>& lam, const MvaOptions& opt,
                              InnerSolve inner) {
    const T zero = num_traits<T>::from_int(0);
    const T two = num_traits<T>::from_int(2);
    qn::NetworkStruct<T>& V = tr.V;

    // ---- ModelAdapter.refreshServicesFromBase --------------------------------
    // MMT ONLY. `fjFixedPoint.m` reuses the MMT transform across outer
    // iterations and has to undo the previous one's converged services; the H-T
    // arm rebuilds its transform on every call (`ModelAdapter.ht` is invoked
    // unconditionally at forkIter == 1) and has no auxiliary open stream, so
    // there is nothing here for it to reset.
    if (!tr.heidelberger_trivedi) {
    // Base-derived slots are re-read from the base model. The
    // transformation-owned ones -- the join's synchronisation delays and the
    // auxiliary arrivals -- are RESET to what a cold transform would hold, not
    // left alone: they still carry the previous outer iteration's converged
    // values, and keeping them would silently warm-start the fork loop from a
    // different point than the reference does.
    for (std::size_t i = 1; i <= L.stations.size(); ++i)
        for (std::size_t k = 1; k <= L.classes.size(); ++k)
            V.set_service(i, k, L.service[i - 1][k - 1]);
    for (std::size_t s : tr.auxclasses) {
        const std::size_t r = tr.fjclassmap[s];
        for (std::size_t i = 1; i <= L.stations.size(); ++i) {
            if (tr.is_join_station(i)) continue;
            V.set_service(i, s, L.service[i - 1][r - 1]);
        }
    }
    for (std::size_t js : tr.joinStations)
        for (std::size_t k = 1; k <= V.classes.size(); ++k)
            V.set_service(js, k, Distrib<T>::immediate());
    for (std::size_t s : tr.auxclasses)
        V.set_service(tr.sourceStation, s,
                      tr.auxdisabled[s]
                          ? Distrib<T>::disabled_dist()
                          : Distrib<T>::exp_rate(
                                num_traits<T>::from_double(GlobalConstants::FineTol)));
    }

    // QN starts at the Immediate sentinel and QN_1 at zero, so the first two
    // passes can never be mistaken for a converged pair.
    Matrix<T> QN(1, L.nclasses, num_traits<T>::from_double(GlobalConstants::Immediate));
    Matrix<T> QN_1(1, L.nclasses, zero);
    MvaSolution<T> out;
    bool forkLoop = true;
    int forkIter = 0;
    while (forkLoop && forkIter < opt.iter_max) {
        ++forkIter;
        // `fjFixedPoint.m:77`: the auxiliary arrival update is guarded on the
        // method NOT being H-T, whose auxiliary classes are closed and carry no
        // arrival rate to re-set.
        if (forkIter > 1 && !tr.heidelberger_trivedi) {
            // the auxiliary stream carries the fanout-1 branches the circulating
            // job did not take
            for (std::size_t s : tr.auxclasses) {
                if (tr.auxdisabled[s] || !(tr.fanout[s] > 0.0)) continue;
                V.set_service(
                    tr.sourceStation, s,
                    Distrib<T>::exp_rate(
                        T(num_traits<T>::from_double(tr.fanout[s] - 1.0) * lam[s])));
            }
        }
        V.refresh_chains();

        if (line::util::LineConsole::owns_log()) {
            if (QN_1.rows() == QN.rows() && QN_1.cols() == QN.cols()) {
                double moved = 0.0;
                for (std::size_t i = 0; i < QN.rows(); ++i)
                    for (std::size_t j = 0; j < QN.cols(); ++j)
                        moved = std::max(moved, std::fabs(num_traits<T>::to_double(QN_1(i, j)) -
                                                          num_traits<T>::to_double(QN(i, j))));
                line::util::LineConsole::step(
                    "fork-join iteration %d: queue lengths moved by at most %.3e", forkIter, moved);
            } else {
                line::util::LineConsole::step(
                    "fork-join iteration %d: transformed model rebuilt", forkIter);
            }
        }

        if (fj_converged(QN_1, QN, opt.iter_tol) && forkIter > 2)
            forkLoop = false;
        else
            QN_1 = QN;

        out = inner(V);

        if (tr.heidelberger_trivedi) {
            fj_ht_sync_delays(L, tr, out);
            fj_ht_merge(L, tr, out);
            QN = out.Q;
            continue;
        }

        // ---- synchronisation delays and the forkLambda update ----------------
        // One pass per fork, as in the reference: each fork drives only ITS OWN
        // auxiliary classes, off its own join.
        for (std::size_t fa = 0; fa < tr.forks.size(); ++fa) {
            const std::size_t fnode = tr.forks[fa].node;
            const std::size_t jstat = tr.forks[fa].joinStation;

            // A FORK WITH NO JOIN. `fjFixedPoint.m:142-152` drives forkLambda
            // from the fork's OWN firing rate rather than a join's throughput,
            // and charges no synchronisation delay, there being no station to
            // charge it to. `fj_sort_forks` has already left such a fork outer
            // with itself as parent, which is what `sortForks` returns for it.
            if (jstat == 0) {
                const std::vector<T> tnfork = fj_fork_tput(L, tr, out.Tp, fa);
                for (std::size_t s : tr.auxclasses) {
                    if (tr.fjforkmap[s] != fa) continue;
                    lam[s] = T((lam[s] + tnfork[tr.fjclassmap[s] - 1]) / two);
                }
                continue;
            }

            for (std::size_t s : tr.auxclasses) {
                if (tr.fjforkmap[s] != fa) continue;
                const std::size_t r = tr.fjclassmap[s];
                T acc = zero;
                for (std::size_t s2 : tr.auxclasses)
                    if (tr.fjclassmap[s2] == r) acc += out.Tp(jstat - 1, s2 - 1);
                out.Tp(jstat - 1, r - 1) =
                    T(out.Tp(jstat - 1, r - 1) + acc - out.Tp(jstat - 1, s - 1));
                lam[s] = T((lam[s] + out.Tp(jstat - 1, r - 1)) / two);

                // An INNER fork's delay is charged by fj_find_paths while the
                // enclosing branch is walked, so computing it again here would
                // count it twice; only an outer fork's own delay is set here.
                if (!tr.forks[fa].outer[r]) continue;

                std::vector<T> ri;
                std::vector<std::size_t> merge{r, s};
                fj_find_paths(tr, fnode, tr.forks[fa].joinNode, r, merge, out.Q, out.Tp, zero, ri,
                              {});
                T sync = zero;
                if (!ri.empty()) {
                    // tasksPerLink = w sends w IDENTICAL tasks down each link, so the
                    // join synchronises on w*B siblings and not on B: the sibling set
                    // is each branch's completion time REPLICATED w times, and the
                    // order statistic is taken over that multiset. Scaling E[X_(k)] by
                    // w instead (what this did before) is w*H_B/mu where the answer is
                    // H_(w*B)/mu, which OVER-states the delay by more the larger w is.
                    // w = 1 replicates to itself, so nothing moves on an ordinary fork.
                    double wd = tr.forks[fa].fanOut;
                    std::size_t w = (wd >= 1.0) ? static_cast<std::size_t>(wd + 0.5) : 1;
                    if (w > 1) {
                        const std::vector<T> base = ri;
                        for (std::size_t k = 1; k < w; ++k)
                            ri.insert(ri.end(), base.begin(), base.end());
                    }
                    T mean = zero;
                    for (const T& x : ri) mean += x;
                    mean = T(mean / num_traits<T>::from_int(static_cast<long>(ri.size())));
                    // The join fires on the k-th sibling completion, k = the sibling
                    // count on a standard join and the declared quorum on a PARTIAL
                    // one. The quorum is declared against the SIBLING count w*B, which
                    // is the replicated length.
                    const T d0 = fj_expected_ordstat(
                        ri, fj_join_quorum(L, tr.forks[fa].joinNode, ri.size()));
                    const T raw = T(d0 - mean);
                    // The quorum can be met BEFORE the branch this transform's own
                    // token walks, and then the parent ought to leave ahead of it.
                    // The MMT cannot express that: its token is a job of the closed
                    // chain and must finish its branch, and that closed token is what
                    // keeps the branch stable, so it cannot be made open either. The
                    // delay floors at zero, which OVER-states the cycle time.
                    // The reference raises a line_warning here; this port has no
                    // warning channel out of the fixed point (see mva_dispatch.h).
                    sync = raw > zero ? raw : zero;
                }
                const Distrib<T> d = fj_exp_fit_mean(sync);
                V.set_service(jstat, s, d);
                V.set_service(jstat, r, d);
            }
        }

        // ---- merge the auxiliary classes back --------------------------------
        // EVERY join and the Source keep their ORIGINAL-class throughputs: the
        // auxiliary tokens are extra fork branches, not extra completions of the
        // class, so adding them there would double the rate at which the fork is
        // seen to fire.
        const std::size_t nrows = out.Tp.rows();
        std::vector<std::vector<T>> tnJoin(tr.joinStations.size(),
                                           std::vector<T>(L.nclasses, zero));
        std::vector<T> tnSource(L.nclasses, zero);
        for (std::size_t a = 0; a < tr.joinStations.size(); ++a)
            for (std::size_t k = 0; k < L.nclasses; ++k)
                tnJoin[a][k] = out.Tp(tr.joinStations[a] - 1, k);
        for (std::size_t k = 0; k < L.nclasses; ++k)
            tnSource[k] = out.Tp(tr.sourceStation - 1, k);
        for (std::size_t s : tr.auxclasses) {
            const std::size_t r = tr.fjclassmap[s];
            for (std::size_t i = 0; i < nrows; ++i) {
                out.Q(i, r - 1) = T(out.Q(i, r - 1) + out.Q(i, s - 1));
                out.U(i, r - 1) = T(out.U(i, r - 1) + out.U(i, s - 1));
                out.Tp(i, r - 1) = T(out.Tp(i, r - 1) + out.Tp(i, s - 1));
                out.R(i, r - 1) = out.Tp(i, r - 1) != zero
                                      ? T(out.Q(i, r - 1) / out.Tp(i, r - 1))
                                      : zero;
            }
        }
        for (std::size_t a = 0; a < tr.joinStations.size(); ++a)
            for (std::size_t k = 0; k < L.nclasses; ++k)
                out.Tp(tr.joinStations[a] - 1, k) = tnJoin[a][k];
        for (std::size_t k = 0; k < L.nclasses; ++k)
            out.Tp(tr.sourceStation - 1, k) = tnSource[k];

        // drop the auxiliary columns; they are all appended after the original
        // ones, so the original block is the leading submatrix
        out.Q = fj_leading_cols(out.Q, L.nclasses);
        out.U = fj_leading_cols(out.U, L.nclasses);
        out.R = fj_leading_cols(out.R, L.nclasses);
        out.Tp = fj_leading_cols(out.Tp, L.nclasses);
        QN = out.Q;
    }
    // A fork model that comes back non-finite has not been solved. It happens
    // when the auxiliary open stream saturates a branch station: the mixed MVA
    // then divides by a non-positive slack and the NaN propagates, where it
    // reads as a solved model with a few blank entries. Refuse by name instead.
    for (std::size_t i = 0; i < out.Q.rows(); ++i)
        for (std::size_t k = 0; k < out.Q.cols(); ++k)
            if (!std::isfinite(num_traits<T>::to_double(out.Q(i, k))) ||
                !std::isfinite(num_traits<T>::to_double(out.U(i, k))) ||
                !std::isfinite(num_traits<T>::to_double(out.R(i, k))) ||
                !std::isfinite(num_traits<T>::to_double(out.Tp(i, k))))
                throw NumericError(
                    "SolverMVA: the fork-join fixed point of model '" + L.name +
                    "' returned a non-finite metric; the auxiliary open classes the transform "
                    "adds have saturated one of the branch stations");
    return out;
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_FJ_DRIVER_H
